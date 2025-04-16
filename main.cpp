#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <future>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <exception>
#include <csignal>
#include <cstdlib>
#include <system_error>

// Include GDAL/OGR headers.
#include "gdal_priv.h"
#include "ogrsf_frmts.h"

namespace fs = std::filesystem;

// --------------------
// Helper Functions
// --------------------
std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\n\r\"");
    size_t end = s.find_last_not_of(" \t\n\r\"");
    if (start == std::string::npos || end == std::string::npos)
        return "";
    return s.substr(start, end - start + 1);
}

// --------------------
// Data Structures
// --------------------
struct Record {
    std::string gbifID;
    std::string species;
    std::string countryCode;
    double lon;          // decimalLongitude
    double lat;          // decimalLatitude
    std::string eventDate;
    int row_id;
};

struct RasterCell {
    int id;
    double centroid_lon;
    double centroid_lat;
    int count;
};

// --------------------
// Global Cleaning Summary
// --------------------
struct CleaningSummary {
    std::string functionName;
    size_t beforeCount;
    size_t afterCount;
    bool skipped;
    double durationSeconds; // in seconds
};

std::vector<CleaningSummary> cleaningSummaries;

// --------------------
// GDAL Geometry Helpers
// --------------------
std::vector<OGRGeometry*> loadPolygons(const std::string &shapefilePath) {
    std::vector<OGRGeometry*> polygons;
    GDALDataset *poDS = (GDALDataset*)GDALOpenEx(shapefilePath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
    if (poDS == nullptr) {
        std::cerr << "Failed to open shapefile: " << shapefilePath << std::endl;
        return polygons;
    }
    OGRLayer *poLayer = poDS->GetLayer(0);
    poLayer->ResetReading();
    OGRFeature *poFeature;
    while ((poFeature = poLayer->GetNextFeature()) != nullptr) {
        OGRGeometry *poGeometry = poFeature->GetGeometryRef();
        if (poGeometry != nullptr) {
            OGRGeometry *poClone = poGeometry->clone();
            polygons.push_back(poClone);
        }
        OGRFeature::DestroyFeature(poFeature);
    }
    GDALClose(poDS);
    return polygons;
}

void freePolygons(std::vector<OGRGeometry*> &polygons) {
    for (auto geom : polygons) {
        if (geom) {
            OGRGeometryFactory::destroyGeometry(geom);
        }
    }
    polygons.clear();
}

// --------------------
// Timing Helper Template
// --------------------
template<typename F, typename... Args>
auto measureDuration(F func, Args&&... args)
    -> std::pair<decltype(func(std::forward<Args>(args)...)), double> {
    auto start = std::chrono::steady_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return {result, elapsed.count()};
}

// --------------------
// Spatial and Utility Functions
// --------------------
double haversine(double lon1, double lat1, double lon2, double lat2) {
    const double R = 6371.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(std::sqrt(a), std::sqrt(1 - a));
    return R * c;
}

std::vector<Record> reassign_row_ids(const std::vector<Record> &records) {
    std::vector<Record> newRecords = records;
    for (size_t i = 0; i < newRecords.size(); i++) {
        newRecords[i].row_id = static_cast<int>(i);
    }
    return newRecords;
}

// --------------------
// Full Cleaning Functions
// --------------------

// cc_val: Remove records with invalid coordinates.
std::vector<Record> cc_val(const std::vector<Record> &records,
                           const std::string &lonCol = "decimalLongitude",
                           const std::string &latCol = "decimalLatitude",
                           const std::string &value = "clean",
                           bool verbose = true) {
    std::vector<Record> valid;
    for (const auto &rec : records) {
        bool invalid = false;
        if (std::isnan(rec.lon) || std::isnan(rec.lat))
            invalid = true;
        if (rec.lon < -180.0 || rec.lon > 180.0)
            invalid = true;
        if (rec.lat < -90.0 || rec.lat > 90.0)
            invalid = true;
        if (!invalid)
            valid.push_back(rec);
    }
    size_t removed = records.size() - valid.size();
    if (verbose) {
        if (value == "clean")
            std::cout << "cc_val: Removed " << removed << " records." << std::endl;
        else
            std::cout << "cc_val: Flagged " << removed << " records." << std::endl;
    }
    return valid;
}

// cc_outl: Remove outlier records for each species based on distance.
// This version prints progress for each species group processed.
std::vector<Record> cc_outl(const std::vector<Record> &records,
                            const std::string &method = "quantile",
                            double mltpl = 5.0,
                            double tdi = 1000.0,
                            int min_occs = 7,
                            bool thinning = false,
                            double thinning_res = 0.5) {
    std::unordered_map<std::string, std::vector<Record>> speciesMap;
    for (const auto &rec : records) {
        speciesMap[rec.species].push_back(rec);
    }
    std::vector<bool> keep(records.size(), true);
    int groupCounter = 0;
    for (const auto &kv : speciesMap) {
        groupCounter++;
        const auto &group = kv.second;
        if (group.size() < static_cast<size_t>(min_occs)) {
            std::cerr << "warning: species " << kv.first << " has fewer than " << min_occs << " records. Skipping outlier test for this species." << std::endl;
            continue;
        }
        size_t n = group.size();
        std::vector<double> values(n, 0.0);
        bool raster_flag = (group.size() >= 10000) || thinning;
        if (raster_flag) {
            std::vector<double> minDistances(n, std::numeric_limits<double>::max());
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    if (i == j) continue;
                    double d = haversine(group[i].lon, group[i].lat,
                                         group[j].lon, group[j].lat);
                    if (d < minDistances[i])
                        minDistances[i] = d;
                }
            }
            values = minDistances;
        } else {
            if (method == "distance") {
                std::vector<double> minDistances(n, std::numeric_limits<double>::max());
                for (size_t i = 0; i < n; i++) {
                    for (size_t j = 0; j < n; j++) {
                        if (i == j) continue;
                        double d = haversine(group[i].lon, group[i].lat,
                                             group[j].lon, group[j].lat);
                        if (d < minDistances[i])
                            minDistances[i] = d;
                    }
                }
                values = minDistances;
            } else {
                std::vector<double> meanDistances(n, 0.0);
                for (size_t i = 0; i < n; i++) {
                    double sum = 0.0;
                    int count = 0;
                    for (size_t j = 0; j < n; j++) {
                        if (i == j) continue;
                        double d = haversine(group[i].lon, group[i].lat,
                                             group[j].lon, group[j].lat);
                        sum += d;
                        count++;
                    }
                    meanDistances[i] = (count > 0 ? sum / count : 0.0);
                }
                values = meanDistances;
            }
        }
        if (method == "distance") {
            for (size_t i = 0; i < n; i++) {
                if (values[i] > tdi)
                    keep[group[i].row_id] = false;
            }
        } else if (method == "quantile") {
            std::vector<double> sorted = values;
            std::sort(sorted.begin(), sorted.end());
            double q1 = sorted[sorted.size() / 4];
            double q3 = sorted[(3 * sorted.size()) / 4];
            double iqr = q3 - q1;
            double threshold = q3 + mltpl * iqr;
            for (size_t i = 0; i < n; i++) {
                if (values[i] > threshold)
                    keep[group[i].row_id] = false;
            }
        } else if (method == "mad") {
            std::vector<double> sorted = values;
            std::sort(sorted.begin(), sorted.end());
            double median = sorted[sorted.size() / 2];
            std::vector<double> absDev;
            for (double v : values) {
                absDev.push_back(std::abs(v - median));
            }
            std::sort(absDev.begin(), absDev.end());
            double mad = absDev[absDev.size() / 2];
            double threshold = median + mltpl * mad;
            for (size_t i = 0; i < n; i++) {
                if (values[i] > threshold)
                    keep[group[i].row_id] = false;
            }
        }
        std::cout << "Processed species group \"" << kv.first << "\" (" << group.size() << " records)" << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    std::cout << "cc_outl: removed " << (records.size() - result.size()) << " records." << std::endl;
    return result;
}

/*
 * cc_cen: Removes records too close to specified centroids.
 */
std::vector<Record> cc_cen(const std::vector<Record> &records,
                           const std::vector<std::pair<double, double>> &centroids,
                           double buffer = 1000.0) {
    double buffer_km = buffer / 1000.0;
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        for (const auto &cen : centroids) {
            double d = haversine(records[i].lon, records[i].lat, cen.first, cen.second);
            if (d <= buffer_km) {
                keep[i] = false;
                break;
            }
        }
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    std::cout << "cc_cen: removed " << (records.size() - result.size()) << " records." << std::endl;
    return result;
}

/*
 * cc_gbif: Removes records near the GBIF headquarters (12.58, 55.67).
 */
std::vector<Record> cc_gbif(const std::vector<Record> &records,
                            double buffer = 1000.0,
                            bool geod = true) {
    double gbif_lon = 12.58;
    double gbif_lat = 55.67;
    double buffer_km = buffer / 1000.0;
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        double d = haversine(records[i].lon, records[i].lat, gbif_lon, gbif_lat);
        if (d <= buffer_km)
            keep[i] = false;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    std::cout << "cc_gbif: removed " << (records.size() - result.size()) << " records." << std::endl;
    return result;
}

/*
 * cc_cap: Removes records that are within a buffer (in meters) of any given capital coordinates.
 */
std::vector<Record> cc_cap(const std::vector<Record> &records,
                           const std::vector<std::pair<double, double>> &capitals,
                           double buffer = 10000,
                           bool geod = true,
                           bool verify = false,
                           const std::string &value = "clean",
                           bool verbose = true) {
    double buffer_km = buffer / 1000.0;
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        for (const auto &cap : capitals) {
            double d = haversine(records[i].lon, records[i].lat, cap.first, cap.second);
            if (d <= buffer_km) {
                keep[i] = false;
                break;
            }
        }
    }
    if (verbose) {
        std::cout << "cc_cap: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    return result;
}

/*
 * cc_coun: Removes records where the reported country code does not match the country determined by coordinates.
 * (Heuristic: if latitude > 0 assume "USA", otherwise "RUS".)
 */
std::vector<Record> cc_coun(const std::vector<Record> &records,
                            const std::string &value = "clean",
                            bool verbose = true) {
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        std::string ref_code = (records[i].lat > 0 ? "USA" : "RUS");
        if (ref_code != records[i].countryCode)
            keep[i] = false;
    }
    if (verbose) {
        std::cout << "cc_coun: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    return result;
}

/*
 * cc_dupl: Removes duplicate records based on species and coordinates.
 */
std::vector<Record> cc_dupl(const std::vector<Record> &records,
                            const std::string &value = "clean",
                            bool verbose = true) {
    std::unordered_map<std::string, bool> seen;
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        std::ostringstream key;
        key << records[i].species << "_" << records[i].lon << "_" << records[i].lat;
        std::string keyStr = key.str();
        if (seen.find(keyStr) != seen.end())
            keep[i] = false;
        else
            seen[keyStr] = true;
    }
    if (verbose) {
        std::cout << "cc_dupl: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " duplicate records." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    return result;
}

/*
 * cc_equ: Removes records with identical latitude and longitude.
 * "absolute" compares absolute values; "identical" compares the actual values.
 */
std::vector<Record> cc_equ(const std::vector<Record> &records,
                           const std::string &test = "absolute",
                           const std::string &value = "clean",
                           bool verbose = true) {
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        if (test == "absolute") {
            if (std::abs(records[i].lon) == std::abs(records[i].lat))
                keep[i] = false;
        } else if (test == "identical") {
            if (records[i].lon == records[i].lat)
                keep[i] = false;
        }
    }
    if (verbose) {
        std::cout << "cc_equ: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records with equal lat/lon." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    return result;
}

/*
 * cc_inst: Removes records that fall within a buffer of institutional coordinates.
 */
std::vector<Record> cc_inst(const std::vector<Record> &records,
                            const std::vector<std::pair<double, double>> &instCoords,
                            double buffer = 1000.0,
                            bool verbose = true) {
    double buffer_km = buffer / 1000.0;
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        for (const auto &inst : instCoords) {
            double d = haversine(records[i].lon, records[i].lat, inst.first, inst.second);
            if (d <= buffer_km) {
                keep[i] = false;
                break;
            }
        }
    }
    if (verbose) {
        std::cout << "cc_inst: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records (institutional coordinates)." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    return result;
}

/*
 * cc_zero: Removes records where both longitude and latitude are exactly zero.
 */
std::vector<Record> cc_zero(const std::vector<Record> &records,
                           bool verbose = true) {
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        if (records[i].lon == 0.0 && records[i].lat == 0.0)
            keep[i] = false;
    }
    if (verbose) {
        std::cout << "cc_zero: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records (zero coordinates)." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    return result;
}

// --------------------
// New Functions using GDAL/OGR for Spatial Tests
// --------------------

/*
 * cc_sea_gdal: Removes records that are non-terrestrial (not on land)
 * by checking points against a land shapefile.
 */
std::vector<Record> cc_sea_gdal(const std::vector<Record> &records,
                                const std::string &shapefilePath,
                                bool verbose = true,
                                bool *skipped = nullptr) {
    if (verbose) {
        std::cout << "Testing sea coordinates using GDAL (land reference)" << std::endl;
    }
    std::vector<OGRGeometry*> landPolygons = loadPolygons(shapefilePath);
    if (landPolygons.empty()) {
        std::cerr << "No land polygons loaded in cc_sea_gdal. Skipping sea test." << std::endl;
        if (skipped) *skipped = true;
        return records;
    }
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        OGRPoint pt(records[i].lon, records[i].lat);
        bool onLand = false;
        for (auto poly : landPolygons) {
            if (pt.Within(poly)) {
                onLand = true;
                break;
            }
        }
        if (!onLand)
            keep[i] = false;
    }
    if (verbose) {
        std::cout << "cc_sea_gdal: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records (non-terrestrial)." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    freePolygons(landPolygons);
    return result;
}

/*
 * cc_sea_buffland: Alternative approach using a buffered coastline reference.
 * Loads a coastline shapefile and buffers each geometry by 1 degree,
 * then removes records that fall outside any buffered polygon.
 */
std::vector<Record> cc_sea_buffland(const std::vector<Record> &records,
                                    const std::string &coastlineShapefilePath,
                                    bool verbose = true,
                                    bool *skipped = nullptr) {
    if (verbose) {
        std::cout << "Testing sea coordinates using GDAL (buffered coastline)" << std::endl;
    }
    std::vector<OGRGeometry*> originalPolygons = loadPolygons(coastlineShapefilePath);
    if (originalPolygons.empty()) {
        std::cerr << "No coastline polygons loaded in cc_sea_buffland. Skipping buffered sea test." << std::endl;
        if (skipped) *skipped = true;
        return records;
    }
    // Buffer each polygon by 1 degree.
    std::vector<OGRGeometry*> bufferedPolygons;
    for (auto geom : originalPolygons) {
        OGRGeometry *buff = geom->Buffer(1.0);
        if (buff)
            bufferedPolygons.push_back(buff);
    }
    freePolygons(originalPolygons);
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        OGRPoint pt(records[i].lon, records[i].lat);
        bool onLand = false;
        for (auto poly : bufferedPolygons) {
            if (pt.Within(poly)) {
                onLand = true;
                break;
            }
        }
        if (!onLand)
            keep[i] = false;
    }
    if (verbose) {
        std::cout << "cc_sea_buffland: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records (outside buffered coastline)." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    freePolygons(bufferedPolygons);
    return result;
}

/*
 * cc_urb_gdal: Removes records that fall inside urban areas by checking points
 * against an urban areas shapefile.
 */
std::vector<Record> cc_urb_gdal(const std::vector<Record> &records,
                                const std::string &shapefilePath,
                                bool verbose = true,
                                bool *skipped = nullptr) {
    if (verbose) {
        std::cout << "Testing urban areas using GDAL" << std::endl;
    }
    std::vector<OGRGeometry*> urbanPolygons = loadPolygons(shapefilePath);
    if (urbanPolygons.empty()) {
        std::cerr << "No urban polygons loaded in cc_urb_gdal. Skipping urban test." << std::endl;
        if (skipped) *skipped = true;
        return records;
    }
    std::vector<bool> keep(records.size(), true);
    for (size_t i = 0; i < records.size(); i++) {
        OGRPoint pt(records[i].lon, records[i].lat);
        for (auto poly : urbanPolygons) {
            if (pt.Within(poly)) {
                keep[i] = false;
                break;
            }
        }
    }
    if (verbose) {
        std::cout << "cc_urb_gdal: removed " 
                  << (records.size() - std::count(keep.begin(), keep.end(), true))
                  << " records (urban areas)." << std::endl;
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    freePolygons(urbanPolygons);
    return result;
}

// --------------------
// CSV Loading and Utility Functions
// --------------------
std::vector<std::string> split_line(const std::string &line, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<Record> loadCSV(const std::string &filename) {
    std::vector<Record> records;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "error: cannot open csv file " << filename << std::endl;
        return records;
    }
    std::string header;
    if (!std::getline(file, header)) {
        std::cerr << "error: csv file " << filename << " is empty." << std::endl;
        return records;
    }
    char delim = (header.find('\t') != std::string::npos) ? '\t' : ',';
    std::vector<std::string> headerTokens = split_line(header, delim);
    std::vector<std::string> lowerHeader;
    for (auto &t : headerTokens) {
        std::string trimmed = trim(t);
        std::transform(trimmed.begin(), trimmed.end(), trimmed.begin(), ::tolower);
        lowerHeader.push_back(trimmed);
    }
    int idx_gbifID = -1, idx_species = -1, idx_countryCode = -1;
    int idx_decimalLatitude = -1, idx_decimalLongitude = -1, idx_eventDate = -1;
    for (size_t i = 0; i < lowerHeader.size(); i++) {
        if (lowerHeader[i] == "gbifid")
            idx_gbifID = static_cast<int>(i);
        else if (lowerHeader[i] == "species")
            idx_species = static_cast<int>(i);
        else if (lowerHeader[i] == "countrycode")
            idx_countryCode = static_cast<int>(i);
        else if (lowerHeader[i] == "decimallatitude")
            idx_decimalLatitude = static_cast<int>(i);
        else if (lowerHeader[i] == "decimallongitude")
            idx_decimalLongitude = static_cast<int>(i);
        else if (lowerHeader[i] == "eventdate")
            idx_eventDate = static_cast<int>(i);
    }
    if (idx_gbifID == -1 || idx_species == -1 || idx_countryCode == -1 ||
        idx_decimalLatitude == -1 || idx_decimalLongitude == -1 || idx_eventDate == -1) {
        std::cerr << "error: one or more required columns are missing in the CSV header." << std::endl;
        return records;
    }
    int row_id = 0;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> tokens = split_line(line, delim);
        if (tokens.size() < lowerHeader.size())
            continue;
        Record rec;
        rec.gbifID = trim(tokens[idx_gbifID]);
        rec.species = trim(tokens[idx_species]);
        rec.countryCode = trim(tokens[idx_countryCode]);
        try {
            rec.lat = std::stod(tokens[idx_decimalLatitude]);
        } catch (...) {
            rec.lat = std::numeric_limits<double>::quiet_NaN();
        }
        try {
            rec.lon = std::stod(tokens[idx_decimalLongitude]);
        } catch (...) {
            rec.lon = std::numeric_limits<double>::quiet_NaN();
        }
        rec.eventDate = trim(tokens[idx_eventDate]);
        rec.row_id = row_id++;
        records.push_back(rec);
    }
    std::cout << "loaded " << records.size() << " records from csv." << std::endl;
    return records;
}

std::vector<std::pair<double, double>> loadCentroidsCSV(const std::string &filename) {
    std::vector<std::pair<double, double>> centroids;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "error: cannot open centroids csv file " << filename << std::endl;
        return centroids;
    }
    std::string headerLine;
    if (!std::getline(file, headerLine)) {
        std::cerr << "error: empty centroids file or missing header." << std::endl;
        return centroids;
    }
    std::vector<std::string> headerTokens;
    std::stringstream headerStream(headerLine);
    std::string token;
    while (std::getline(headerStream, token, ',')) {
        headerTokens.push_back(trim(token));
    }
    int lonIndex = -1, latIndex = -1;
    for (size_t i = 0; i < headerTokens.size(); i++) {
        std::string lowerToken = headerTokens[i];
        std::transform(lowerToken.begin(), lowerToken.end(), lowerToken.begin(), ::tolower);
        if (lowerToken == "centroid.lon")
            lonIndex = static_cast<int>(i);
        else if (lowerToken == "centroid.lat")
            latIndex = static_cast<int>(i);
    }
    if (lonIndex == -1 || latIndex == -1) {
        std::cerr << "error: header does not contain required columns 'centroid.lon' and 'centroid.lat'" << std::endl;
        return centroids;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> tokens;
        while (std::getline(ss, token, ',')) {
            tokens.push_back(trim(token));
        }
        if (tokens.size() <= static_cast<size_t>(std::max(lonIndex, latIndex)))
            continue;
        double lonVal = std::numeric_limits<double>::quiet_NaN();
        double latVal = std::numeric_limits<double>::quiet_NaN();
        try {
            lonVal = std::stod(tokens[lonIndex]);
        } catch (...) {}
        try {
            latVal = std::stod(tokens[latIndex]);
        } catch (...) {}
        if (!std::isnan(lonVal) && !std::isnan(latVal))
            centroids.push_back(std::make_pair(lonVal, latVal));
    }
    std::cout << "loaded " << centroids.size() << " centroids from " << filename << std::endl;
    return centroids;
}

void segfault_handler(int signum) {
    std::cerr << "segmentation fault (signal " << signum << ") occurred." << std::endl;
    std::exit(signum);
}

// --------------------
// Main Function
// --------------------
int main() {
    std::signal(SIGSEGV, segfault_handler);

    GDALAllRegister();

    auto overallStart = std::chrono::steady_clock::now();

    // File paths for CSV data and centroids.
    std::string csvFile = "downloads/0023500-241107131044228.csv";
    std::string centroidsCSV = "countryref.csv";

    // Shapefile paths for spatial tests.
    std::string landShapefile = "Downloads/ne_10m_land/ne_10m_land.shp";
    std::string coastlineShapefile = "Downloads/ne_10m_coastline/ne_10m_coastline.shp";
    std::string urbanShapefile = "Downloads/ne_10m_urban_areas/ne_10m_urban_areas.shp";

    size_t initialCount = 0;
    std::vector<Record> records = loadCSV(csvFile);
    initialCount = records.size();

    // --- Cleaning Steps with Duration Timing ---
    {
        auto [res, dur] = measureDuration(cc_val, records, "decimalLongitude", "decimalLatitude", "clean", true);
        records = res;
        cleaningSummaries.push_back({"cc_val", initialCount, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        auto [res, dur] = measureDuration(cc_zero, records, true);
        records = res;
        cleaningSummaries.push_back({"cc_zero", before, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        auto [res, dur] = measureDuration(cc_dupl, records, "clean", true);
        records = res;
        cleaningSummaries.push_back({"cc_dupl", before, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        auto [res, dur] = measureDuration(cc_inst, records, std::vector<std::pair<double, double>>{{-77.0365, 38.8977}, {-0.1278, 51.5074}}, 1000.0, true);
        records = res;
        cleaningSummaries.push_back({"cc_inst", before, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        auto [res, dur] = measureDuration(cc_coun, records, "clean", true);
        records = res;
        cleaningSummaries.push_back({"cc_coun", before, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        // Run simpler functions first, then outlier (cc_outl) as a heavy step.
        auto [res, dur] = measureDuration(cc_outl, records, "quantile", 5.0, 1000.0, 7, false, 0.5);
        records = res;
        cleaningSummaries.push_back({"cc_outl", before, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        std::vector<std::pair<double, double>> capitals = { {-77.0369, 38.9072}, {2.3522, 48.8566} };
        auto [res, dur] = measureDuration(cc_cap, records, capitals, 10000, true, false, "clean", true);
        records = res;
        cleaningSummaries.push_back({"cc_cap", before, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        auto [res, dur] = measureDuration(cc_gbif, records, 1000.0, true);
        records = res;
        cleaningSummaries.push_back({"cc_gbif", before, records.size(), false, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        // cc_sea_gdal using land shapefile.
        bool skippedSea = false;
        auto [res, dur] = measureDuration(cc_sea_gdal, records, landShapefile, true, &skippedSea);
        records = res;
        cleaningSummaries.push_back({"cc_sea_gdal", before, records.size(), skippedSea, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        // cc_sea_buffland using coastline shapefile.
        bool skippedBuff = false;
        auto [res, dur] = measureDuration(cc_sea_buffland, records, coastlineShapefile, true, &skippedBuff);
        records = res;
        cleaningSummaries.push_back({"cc_sea_buffland", before, records.size(), skippedBuff, dur});
        records = reassign_row_ids(records);
    }
    {
        size_t before = records.size();
        // cc_urb_gdal using urban areas shapefile.
        bool skippedUrb = false;
        auto [res, dur] = measureDuration(cc_urb_gdal, records, urbanShapefile, true, &skippedUrb);
        records = res;
        cleaningSummaries.push_back({"cc_urb_gdal", before, records.size(), skippedUrb, dur});
        records = reassign_row_ids(records);
    }

    auto overallEnd = std::chrono::steady_clock::now();
    std::chrono::duration<double> overallDur = overallEnd - overallStart;

    std::cout << "\n--- Cleaning Summary ---\n";
    for (const auto &summary : cleaningSummaries) {
        std::cout << summary.functionName << ": "
                  << summary.beforeCount << " records before, "
                  << summary.afterCount << " records after, "
                  << (summary.beforeCount - summary.afterCount) << " removed, took "
                  << summary.durationSeconds << " seconds";
        if (summary.skipped)
            std::cout << " (SKIPPED)";
        std::cout << std::endl;
    }
    std::cout << "Overall: " << initialCount << " records at start, " << records.size()
              << " records after cleaning.\n";
    std::cout << "Total cleaning time: " << overallDur.count() << " seconds" << std::endl;

    std::ofstream outfile("cleaned_data1.4.csv");
    if (!outfile.is_open()) {
        std::cerr << "error: cannot open output file for writing." << std::endl;
        return 1;
    }
    outfile << "gbifID,species,countryCode,decimalLongitude,decimalLatitude,eventDate\n";
    for (const auto &rec : records) {
        outfile << rec.gbifID << "," << rec.species << "," << rec.countryCode << ","
                << rec.lon << "," << rec.lat << "," << rec.eventDate << "\n";
    }
    outfile.close();

    std::cout << "\ncleaned data written to cleaned_data1.4.csv" << std::endl;
    return 0;
}
