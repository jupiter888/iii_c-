#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

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

namespace fs = std::filesystem;

// structure for a record (row)
struct Record {
    std::string species;
    double lon;
    double lat;
    int row_id;
};

// structure for a raster cell (used in raster approximation)
struct RasterCell {
    int id;                   // unique cell id (derived from grid coordinates)
    double centroid_lon;      // center longitude of the cell
    double centroid_lat;      // center latitude of the cell
    int count;                // number of points in this cell
};

// haversine formula (returns distance in kilometers)
double haversine(double lon1, double lat1, double lon2, double lat2) {
    const double R = 6371.0; // earth's radius in km
    double dlon = (lon2 - lon1) * M_PI / 180.0;
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

/*
 * raster approximation function:
 *
 * for a given species group (vector of record) and a thinning resolution (in degrees),
 * this function builds a grid covering the group's extent.
 * each record is assigned to a grid cell.
 * for each unique cell, it computes the centroid (center of the cell) and counts
 * the number of records in that cell.
 * then it computes the distance matrix among these unique cells using the haversine formula.
 * if compute_min is true, the minimum distance from the cell to any other cell is computed;
 * otherwise, the mean distance is computed (weighted by cell count if weights is true).
 * each record is then assigned the computed value for its cell.
 */
std::vector<double> ras_dist(const std::vector<Record>& group,
                             double thinning_res,
                             bool weights,
                             bool compute_min) {
    double min_lon = std::numeric_limits<double>::max();
    double max_lon = std::numeric_limits<double>::lowest();
    double min_lat = std::numeric_limits<double>::max();
    double max_lat = std::numeric_limits<double>::lowest();
    for (const auto& rec : group) {
        if (rec.lon < min_lon) min_lon = rec.lon;
        if (rec.lon > max_lon) max_lon = rec.lon;
        if (rec.lat < min_lat) min_lat = rec.lat;
        if (rec.lat > max_lat) max_lat = rec.lat;
    }

    int ncols = static_cast<int>(std::ceil((max_lon - min_lon) / thinning_res));
    int nrows = static_cast<int>(std::ceil((max_lat - min_lat) / thinning_res));

    std::unordered_map<int, RasterCell> cellMap;
    std::vector<int> recordCellId(group.size(), -1);

    for (size_t i = 0; i < group.size(); i++) {
        int cell_x = static_cast<int>(std::floor((group[i].lon - min_lon) / thinning_res));
        int cell_y = static_cast<int>(std::floor((group[i].lat - min_lat) / thinning_res));
        int cell_id = cell_y * ncols + cell_x;
        recordCellId[i] = cell_id;
        if (cellMap.find(cell_id) == cellMap.end()) {
            double centroid_lon = min_lon + (cell_x + 0.5) * thinning_res;
            double centroid_lat = min_lat + (cell_y + 0.5) * thinning_res;
            cellMap[cell_id] = RasterCell{cell_id, centroid_lon, centroid_lat, 1};
        } else {
            cellMap[cell_id].count += 1;
        }
    }

    std::vector<RasterCell> cells;
    std::unordered_map<int, int> cellIdToIndex;
    for (const auto& kv : cellMap) {
        cellIdToIndex[kv.first] = static_cast<int>(cells.size());
        cells.push_back(kv.second);
    }

    int U = cells.size();
    std::vector<double> cellValues(U, 0.0);
    if (U < 2) {
        std::vector<double> result(group.size(), 0.0);
        return result;
    }

    std::vector<std::vector<double>> dist(U, std::vector<double>(U, 0.0));
    for (int i = 0; i < U; i++) {
        for (int j = 0; j < U; j++) {
            if (i == j) {
                dist[i][j] = 0.0;
            } else {
                dist[i][j] = haversine(cells[i].centroid_lon, cells[i].centroid_lat,
                                         cells[j].centroid_lon, cells[j].centroid_lat);
            }
        }
    }

    for (int i = 0; i < U; i++) {
        if (compute_min) {
            double minD = std::numeric_limits<double>::max();
            for (int j = 0; j < U; j++) {
                if (i == j) continue;
                if (dist[i][j] < minD)
                    minD = dist[i][j];
            }
            cellValues[i] = minD;
        } else {
            double sum = 0.0;
            double totalWeight = 0.0;
            int count = 0;
            for (int j = 0; j < U; j++) {
                if (i == j) continue;
                if (weights) {
                    sum += dist[i][j] * cells[j].count;
                    totalWeight += cells[j].count;
                } else {
                    sum += dist[i][j];
                    count++;
                }
            }
            cellValues[i] = (weights ? (totalWeight > 0 ? sum / totalWeight : 0.0)
                                     : (count > 0 ? sum / count : 0.0));
        }
    }

    std::vector<double> result(group.size(), 0.0);
    for (size_t i = 0; i < group.size(); i++) {
        int cell_id = recordCellId[i];
        int idx = cellIdToIndex[cell_id];
        result[i] = cellValues[idx];
    }
    return result;
}

// function to load country centroids from a csv file
std::vector<std::pair<double,double>> loadCentroidsCSV(const std::string& filename) {
    std::vector<std::pair<double,double>> centroids;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "error: cannot open centroids csv file " << filename << std::endl;
        return centroids;
    }

    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "error: empty centroids file or missing header." << std::endl;
        return centroids;
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        double lonVal = std::numeric_limits<double>::quiet_NaN();
        double latVal = std::numeric_limits<double>::quiet_NaN();

        // assuming the csv columns are:
        // type,iso3,area_sqkm,centroid.lon,centroid.lat,...
        if (!std::getline(ss, token, ',')) continue; // type
        if (!std::getline(ss, token, ',')) continue; // iso3
        if (!std::getline(ss, token, ',')) continue; // area_sqkm
        if (!std::getline(ss, token, ',')) continue; // centroid.lon
        try {
            lonVal = std::stod(token);
        } catch(...) {}
        if (!std::getline(ss, token, ',')) continue; // centroid.lat
        try {
            latVal = std::stod(token);
        } catch(...) {}

        if (!std::isnan(lonVal) && !std::isnan(latVal)) {
            centroids.push_back(std::make_pair(lonVal, latVal));
        }
    }
    std::cout << "loaded " << centroids.size() << " centroids from " << filename << std::endl;
    return centroids;
}

// cc_val: validate coordinates
std::vector<Record> cc_val(const std::vector<Record>& records) {
    std::vector<Record> valid;
    for (const auto& rec : records) {
        if (!std::isnan(rec.lon) && !std::isnan(rec.lat) &&
            rec.lon >= -180.0 && rec.lon <= 180.0 &&
            rec.lat >= -90.0 && rec.lat <= 90.0) {
            valid.push_back(rec);
        }
    }
    std::cout << "cc_val: removed " << (records.size() - valid.size()) << " records." << std::endl;
    return valid;
}

// cc_outl: identify geographic outliers with support for raster estimation
std::vector<Record> cc_outl(const std::vector<Record>& records,
                            const std::string &method = "quantile",
                            double mltpl = 5.0,
                            double tdi = 1000.0,
                            int min_occs = 7,
                            bool thinning = false,
                            double thinning_res = 0.5) {
    std::unordered_map<std::string, std::vector<Record>> speciesMap;
    for (const auto& rec : records) {
        speciesMap[rec.species].push_back(rec);
    }
    std::vector<bool> keep(records.size(), true);

    for (const auto& kv : speciesMap) {
        const auto &group = kv.second;
        if (group.size() < static_cast<size_t>(min_occs)) {
            std::cerr << "warning: species " << kv.first << " has fewer than min_occs records. skipping outlier test for this species." << std::endl;
            continue;
        }
        size_t n = group.size();
        std::vector<double> values(n, 0.0);

        bool raster_flag = (group.size() >= 10000) || thinning;
        if (raster_flag) {
            bool useWeights = (!thinning);
            if (method == "distance") {
                values = ras_dist(group, thinning_res, useWeights, true);
            } else {
                values = ras_dist(group, thinning_res, useWeights, false);
            }
        } else {
            if (method == "distance") {
                std::vector<double> minDistances(n, std::numeric_limits<double>::max());
                for (size_t i = 0; i < n; i++) {
                    for (size_t j = 0; j < n; j++) {
                        if (i == j) continue;
                        double d = haversine(group[i].lon, group[i].lat, group[j].lon, group[j].lat);
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
                        double d = haversine(group[i].lon, group[i].lat, group[j].lon, group[j].lat);
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
                if (values[i] > tdi) {
                    keep[group[i].row_id] = false;
                }
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
    }
    std::vector<Record> result;
    for (size_t i = 0; i < records.size(); i++) {
        if (keep[i])
            result.push_back(records[i]);
    }
    std::cout << "cc_outl: removed " << (records.size() - result.size()) << " records." << std::endl;
    return result;
}

// cc_cen: remove records near country/province centroids
std::vector<Record> cc_cen(const std::vector<Record>& records,
                           const std::vector<std::pair<double,double>>& centroids,
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

// cc_gbif: remove records near gbif headquarters (fixed coordinate)
std::vector<Record> cc_gbif(const std::vector<Record>& records,
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

// load data from csv (expects header: species,decimalLongitude,decimalLatitude)
std::vector<Record> loadCSV(const std::string &filename) {
    std::vector<Record> records;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "error: cannot open csv file " << filename << std::endl;
        return records;
    }
    std::string line;
    std::getline(file, line);
    int row_id = 0;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        Record rec;
        if (!std::getline(ss, token, ','))
            continue;
        rec.species = token;
        if (!std::getline(ss, token, ','))
            continue;
        try {
            rec.lon = std::stod(token);
        } catch (...) {
            rec.lon = std::numeric_limits<double>::quiet_NaN();
        }
        if (!std::getline(ss, token, ','))
            continue;
        try {
            rec.lat = std::stod(token);
        } catch (...) {
            rec.lat = std::numeric_limits<double>::quiet_NaN();
        }
        rec.row_id = row_id++;
        records.push_back(rec);
    }
    std::cout << "loaded " << records.size() << " records from csv." << std::endl;
    return records;
}

// load data from a directory of parquet files
std::vector<Record> loadParquetDirectory(const std::string &directory) {
    std::vector<Record> records;
    int row_id = 0;
    for (const auto &entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".parquet") {
            std::string file_path = entry.path().string();
            std::shared_ptr<arrow::io::ReadableFile> infile;
            auto status = arrow::io::ReadableFile::Open(file_path, &infile);
            if (!status.ok()) {
                std::cerr << "error opening parquet file: " << file_path << std::endl;
                continue;
            }
            std::unique_ptr<parquet::arrow::FileReader> reader;
            status = parquet::arrow::OpenFile(infile, arrow::default_memory_pool(), &reader);
            if (!status.ok()) {
                std::cerr << "error creating parquet reader for file: " << file_path << std::endl;
                continue;
            }
            std::shared_ptr<arrow::Table> table;
            status = reader->ReadTable(&table);
            if (!status.ok()) {
                std::cerr << "error reading parquet table from file: " << file_path << std::endl;
                continue;
            }
            auto speciesArray = table->GetColumnByName("species");
            auto lonArray = table->GetColumnByName("decimalLongitude");
            auto latArray = table->GetColumnByName("decimalLatitude");
            if (!speciesArray || !lonArray || !latArray) {
                std::cerr << "error: missing required columns in parquet file: " << file_path << std::endl;
                continue;
            }
            for (int c = 0; c < speciesArray->num_chunks(); c++) {
                auto speciesChunk = speciesArray->chunk(c);
                auto lonChunk = lonArray->chunk(c);
                auto latChunk = latArray->chunk(c);
                int64_t num_rows = speciesChunk->length();
                for (int64_t i = 0; i < num_rows; i++) {
                    Record rec;
                    auto speciesScalarResult = speciesChunk->GetScalar(i);
                    if (!speciesScalarResult.ok())
                        continue;
                    auto speciesScalar = std::static_pointer_cast<arrow::StringScalar>(speciesScalarResult.ValueOrDie());
                    rec.species = speciesScalar->value->ToString();
                    
                    auto lonScalarResult = lonChunk->GetScalar(i);
                    if (!lonScalarResult.ok())
                        continue;
                    auto lonScalar = std::static_pointer_cast<arrow::DoubleScalar>(lonScalarResult.ValueOrDie());
                    rec.lon = lonScalar->value;
                    
                    auto latScalarResult = latChunk->GetScalar(i);
                    if (!latScalarResult.ok())
                        continue;
                    auto latScalar = std::static_pointer_cast<arrow::DoubleScalar>(latScalarResult.ValueOrDie());
                    rec.lat = latScalar->value;
                    
                    rec.row_id = row_id++;
                    records.push_back(rec);
                }
            }
            std::cout << "loaded records from file: " << file_path << std::endl;
        }
    }
    std::cout << "total records loaded from parquet: " << records.size() << std::endl;
    return records;
}

int main() {
    // input paths
    std::string parquetDir = "/users/njord888/desktop/myco_pkgs/bench_this/data/parquet_output";
    std::string csvFile = "/users/njord888/downloads/0023500-241107131044228.csv";
    std::string centroidsCSV = "/users/njord888/countryref.csv";
    
    std::vector<Record> records;
    if (fs::exists(parquetDir) && fs::is_directory(parquetDir)) {
        records = loadParquetDirectory(parquetDir);
        if (records.empty()) {
            std::cout << "no records loaded from parquet. falling back to csv." << std::endl;
            records = loadCSV(csvFile);
        }
    } else {
        std::cout << "parquet directory not found. loading csv." << std::endl;
        records = loadCSV(csvFile);
    }
    
    // step 1: validate coordinates
    records = cc_val(records);
    
    // step 2: remove geographic outliers
    records = cc_outl(records, "quantile", 5.0, 1000.0, 7, false, 0.5);
    
    // step 3: load real centroids from csv and remove records near centroids
    std::vector<std::pair<double,double>> centroids = loadCentroidsCSV(centroidsCSV);
    records = cc_cen(records, centroids, 1000.0);
    
    // step 4: remove records near gbif headquarters
    records = cc_gbif(records, 1000.0, true);
    
    // write final cleaned data to csv
    std::ofstream outfile("cleaned_data.csv");
    if (!outfile.is_open()) {
        std::cerr << "error: cannot open output file for writing." << std::endl;
        return 1;
    }
    outfile << "species,decimalLongitude,decimalLatitude\n";
    for (const auto &rec : records) {
        outfile << rec.species << "," << rec.lon << "," << rec.lat << "\n";
    }
    outfile.close();
    
    std::cout << "cleaned data written to cleaned_data.csv" << std::endl;
    return 0;
}
