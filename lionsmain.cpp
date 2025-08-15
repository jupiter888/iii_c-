// brew install gdal gcc@14
// COMPILE
// g++-14 -std=c++17 -O3 -march=native -fopenmp \ main.cpp \ $(gdal-config --cflags) $(gdal-config --libs) \ -Wl,-rpath,"$(gdal-config --prefix)/lib" \ -o data_cleaner
// export OMP_NUM_THREADS=$(sysctl -n hw.physicalcpu)
// export OMP_PROC_BIND=close
// export OMP_PLACES=cores
// caffeinate -i ./data_cleaner
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <csignal>
#include <cstdlib>
#include <stdexcept>

// OpenMP
#include <omp.h>

// GDAL/OGR headers
#include "gdal_priv.h"
#include "gdal_alg.h"
#include "ogrsf_frmts.h"
#include "cpl_conv.h" // CPLFree

namespace fs = std::filesystem;

// --------------------
// Helper Functions
// --------------------
std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\n\r\"");
    size_t end   = s.find_last_not_of(" \t\n\r\"");
    if (start == std::string::npos || end == std::string::npos)
        return "";
    return s.substr(start, end - start + 1);
}

std::string upper(const std::string &s) {
    std::string u = s;
    std::transform(u.begin(), u.end(), u.begin(), ::toupper);
    return u;
}

// --------------------
// Data Structure
// --------------------
struct Record {
    std::string gbifID;
    std::string species;
    std::string countryCode; // always ISO3 after conversion
    double      lon = std::numeric_limits<double>::quiet_NaN();
    double      lat = std::numeric_limits<double>::quiet_NaN();
    std::string eventDate;
    int         row_id = -1;
};

// --------------------
// Cleaning Summary
// --------------------
struct CleaningSummary {
    std::string functionName;
    size_t      beforeCount;
    size_t      afterCount;
    bool        skipped;
    double      durationSeconds;
};
std::vector<CleaningSummary> cleaningSummaries;

// --------------------
// GDAL Geometry Helpers
// --------------------
std::vector<OGRGeometry*> loadPolygons(const std::string &shapefilePath) {
    std::vector<OGRGeometry*> polygons;
    GDALDataset *ds = (GDALDataset*)GDALOpenEx(
        shapefilePath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr
    );
    if (!ds) {
        std::cerr << "Failed to open shapefile: " << shapefilePath << "\n";
        return polygons;
    }
    OGRLayer *layer = ds->GetLayer(0);
    if (!layer) {
        std::cerr << "No layer in shapefile: " << shapefilePath << "\n";
        GDALClose(ds);
        return polygons;
    }
    layer->ResetReading();
    OGRFeature *feat;
    while ((feat = layer->GetNextFeature())) {
        OGRGeometry *geom = feat->GetGeometryRef();
        if (geom) polygons.push_back(geom->clone());
        OGRFeature::DestroyFeature(feat);
    }
    GDALClose(ds);
    return polygons;
}

void freePolygons(std::vector<OGRGeometry*> &polygons) {
    for (auto *g : polygons) {
        OGRGeometryFactory::destroyGeometry(g);
    }
    polygons.clear();
}

// --------------------
// Load Country Features
// --------------------
std::vector<std::pair<OGRGeometry*, std::pair<std::string,std::string>>>
loadCountryFeaturesWithCodes(const std::string &shapefilePath,
                             const std::string &iso2Field = "ISO_A2",
                             const std::string &iso3Field = "ISO_A3") {
    std::vector<std::pair<OGRGeometry*, std::pair<std::string,std::string>>> feats;
    GDALDataset *ds = (GDALDataset*)GDALOpenEx(
        shapefilePath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr
    );
    if (!ds) {
        std::cerr << "Failed to open country shapefile: " << shapefilePath << "\n";
        return feats;
    }
    OGRLayer *layer = ds->GetLayer(0);
    if (!layer) {
        std::cerr << "No layer in shapefile: " << shapefilePath << "\n";
        GDALClose(ds);
        return feats;
    }
    layer->ResetReading();
    OGRFeature *feat;
    while ((feat = layer->GetNextFeature())) {
        OGRGeometry *geom = feat->GetGeometryRef();
        if (!geom) { OGRFeature::DestroyFeature(feat); continue; }
        OGRGeometry *clone = geom->clone();
        const char *a2 = feat->GetFieldAsString(iso2Field.c_str());
        const char *a3 = feat->GetFieldAsString(iso3Field.c_str());
        feats.emplace_back(
            clone,
            std::make_pair(
                a2 ? upper(a2) : std::string{},
                a3 ? upper(a3) : std::string{}
            )
        );
        OGRFeature::DestroyFeature(feat);
    }
    GDALClose(ds);
    return feats;
}

// Build a map from ISO2→ISO3
std::unordered_map<std::string,std::string>
buildIso2To3Map(const std::string &shapefilePath) {
    auto feats = loadCountryFeaturesWithCodes(shapefilePath);
    std::unordered_map<std::string,std::string> m;
    for (auto &p : feats) {
        auto &codes = p.second;
        if (!codes.first.empty() && !codes.second.empty()) {
            m[codes.first] = codes.second;
        }
        OGRGeometryFactory::destroyGeometry(p.first);
    }
    return m;
}

// Convert any two-letter codes in records to three-letter
void convertIso2To3(std::vector<Record> &records,
                    const std::unordered_map<std::string,std::string> &m) {
    for (auto &r : records) {
        if (r.countryCode.size() == 2) {
            auto it = m.find(upper(r.countryCode));
            if (it != m.end()) {
                r.countryCode = it->second;
            }
        }
    }
}

// --------------------
// Count ISO3 Codes
// --------------------
std::unordered_map<std::string,int>
countIso3(const std::vector<Record> &records) {
    std::unordered_map<std::string,int> counts;
    for (auto &r : records) {
        if (r.countryCode.size() == 3) {
            counts[r.countryCode]++;
        }
    }
    return counts;
}

// --------------------
// Timing Helper
// --------------------
template<typename F, typename... Args>
auto measureDuration(F func, Args&&... args)
    -> std::pair<decltype(func(std::forward<Args>(args)...)), double> {
    auto start  = std::chrono::steady_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end    = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return {result, elapsed.count()};
}

// --------------------
// Spatial Utility
// --------------------
double haversine(double lon1, double lat1, double lon2, double lat2) {
    const double R = 6371.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double a = sin(dlat/2)*sin(dlat/2)
             + cos(lat1*M_PI/180.0)*cos(lat2*M_PI/180.0)
               * sin(dlon/2)*sin(dlon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

std::vector<Record> reassign_row_ids(const std::vector<Record> &records) {
    std::vector<Record> out = records;
    for (size_t i = 0; i < out.size(); ++i) {
        out[i].row_id = static_cast<int>(i);
    }
    return out;
}

// --------------------
// Cleaning Functions
// --------------------

// cc_val
std::vector<Record> cc_val(const std::vector<Record> &records,
                           bool verbose = true) {
    std::vector<Record> valid;
    valid.reserve(records.size());
    for (auto &r : records) {
        if (std::isnan(r.lon) || std::isnan(r.lat)) continue;
        if (r.lon < -180 || r.lon > 180)    continue;
        if (r.lat < -90  || r.lat > 90)     continue;
        valid.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_val: Removed "<<(records.size()-valid.size())
                 <<" records.\n";
    return valid;
}

// cc_zero
std::vector<Record> cc_zero(const std::vector<Record> &records,
                            bool verbose = true) {
    std::vector<Record> out;
    out.reserve(records.size());
    for (auto &r : records) {
        if (!(r.lon == 0.0 && r.lat == 0.0))
            out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_zero: Removed "<<(records.size()-out.size())
                 <<" records (zero coords).\n";
    return out;
}

// cc_dupl with 10 m jitter + exact eventDate (kept single-threaded for determinism)
std::vector<Record> cc_dupl(const std::vector<Record> &records,
                            bool verbose = true) {
    const double jitter_deg = 10.0 / 111320.0;
    std::unordered_set<std::string> seen;
    seen.reserve(records.size() * 2);
    std::vector<Record> out;
    out.reserve(records.size());
    for (auto &r : records) {
        int cellX = static_cast<int>(std::floor(r.lon  / jitter_deg));
        int cellY = static_cast<int>(std::floor(r.lat  / jitter_deg));
        std::ostringstream key;
        key << upper(r.species) << "_" 
            << cellX << "_" << cellY << "_"
            << trim(r.eventDate);
        std::string ks = key.str();
        if (seen.insert(ks).second) {
            out.push_back(r);
        }
    }
    if (verbose)
        std::cout<<"cc_dupl: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    return out;
}

// cc_equ — equal lat/lon removal (absolute or identical)
std::vector<Record> cc_equ(const std::vector<Record> &records,
                           const std::string &test = "absolute",
                           bool verbose = true) {
    bool absolute;
    if (test == "absolute") {
        absolute = true;
    } else if (test == "identical") {
        absolute = false;
    } else {
        std::cerr << "cc_equ: unknown test \"" << test
                  << "\" (use \"absolute\" or \"identical\"). Using \"absolute\".\n";
        absolute = true;
    }

    std::vector<Record> out;
    out.reserve(records.size());
    size_t removed = 0;

    for (const auto &r : records) {
        if (std::isnan(r.lon) || std::isnan(r.lat)) {
            out.push_back(r); // leave NaNs to cc_val
            continue;
        }
        bool equal = absolute
            ? (std::fabs(r.lon) == std::fabs(r.lat))
            : (r.lon == r.lat);
        if (equal) ++removed;
        else out.push_back(r);
    }

    if (verbose) {
        std::cout << "cc_equ (" << (absolute ? "absolute" : "identical")
                  << "): Removed " << removed << " records.\n";
    }
    return out;
}

// cc_inst (parallel)
std::vector<Record> cc_inst(const std::vector<Record> &records,
                            const std::vector<std::pair<double,double>> &inst,
                            double buffer = 1000.0,
                            bool verbose = true) {
    double buf_km = buffer / 1000.0;
    std::vector<Record> out;
    out.reserve(records.size());

    #pragma omp parallel
    {
        std::vector<Record> local;
        local.reserve(1024);

        #pragma omp for nowait
        for (int i = 0; i < static_cast<int>(records.size()); ++i) {
            const auto &r = records[i];
            bool keep = true;
            for (auto &iPt : inst) {
                if (haversine(r.lon,r.lat,iPt.first,iPt.second) <= buf_km) {
                    keep = false; break;
                }
            }
            if (keep) local.push_back(r);
        }

        #pragma omp critical
        out.insert(out.end(), local.begin(), local.end());
    }

    if (verbose)
        std::cout<<"cc_inst: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    return out;
}

// cc_coun with debug (kept single-threaded due to OGR geometry checks)
std::vector<Record> cc_coun(const std::vector<Record> &records,
                            const std::string &shapefilePath,
                            const std::string &isoField = "ISO_A3",
                            bool verbose = true) {
    auto feats = loadCountryFeaturesWithCodes(shapefilePath, "ISO_A2", isoField);
    std::vector<Record> out;
    out.reserve(records.size());
    size_t removed = 0, debugCount = 0;
    for (auto &r : records) {
        OGRPoint pt(r.lon, r.lat);
        bool match = false;
        for (auto &p : feats) {
            if (pt.Within(p.first) && p.second.second == upper(r.countryCode)) {
                match = true;
                break;
            }
        }
        if (match) {
            out.push_back(r);
        } else {
            if (debugCount < 10) {
                std::cerr<<"Flagged: "<<r.gbifID
                         <<" @("<<r.lon<<","<<r.lat<<") "
                         <<"countryCode="<<r.countryCode<<"\n";
                ++debugCount;
            }
            ++removed;
        }
    }
    if (verbose)
        std::cout<<"cc_coun: Removed "<<removed
                 <<" records.\n";
    for (auto &p : feats)
        OGRGeometryFactory::destroyGeometry(p.first);
    return out;
}

// cc_outl per-species (parallel across species; thread-safe writes)
std::vector<Record> cc_outl(const std::vector<Record> &records,
                            const std::string &method = "quantile",
                            double mltpl = 5.0,
                            double tdi = 1000.0,
                            int min_occs = 7,
                            bool use_tdi = false,
                            double tdi_weight = 0.0) {
    (void)method; (void)tdi; (void)use_tdi; (void)tdi_weight; // behavior unchanged

    // group by species (serial)
    std::unordered_map<std::string,std::vector<int>> groups;
    groups.reserve(records.size()/2+1);
    for (int i=0;i<(int)records.size();++i) {
        groups[upper(records[i].species)].push_back(i);
    }

    // char mask is safer than vector<bool> for parallel writes
    std::vector<char> keep(records.size(), 1);

    std::vector<std::string> speciesKeys;
    speciesKeys.reserve(groups.size());
    for (auto &kv : groups) speciesKeys.push_back(kv.first);

    // parallel across species
    #pragma omp parallel for schedule(dynamic,1)
    for (int si = 0; si < (int)speciesKeys.size(); ++si) {
        const auto &key = speciesKeys[si];
        const auto &idxs = groups[key];
        if (idxs.size() < (size_t)min_occs) {
            continue;
        }
        size_t n = idxs.size();
        std::vector<double> vals(n, 0.0);
        for (size_t i=0;i<n;i++) {
            const auto &ri = records[idxs[i]];
            double minD = std::numeric_limits<double>::infinity();
            for (size_t j=0;j<n;j++) {
                if (i==j) continue;
                const auto &rj = records[idxs[j]];
                double d = haversine(ri.lon,ri.lat,rj.lon,rj.lat);
                if (d<minD) minD = d;
            }
            vals[i] = minD;
        }
        std::vector<double> sorted = vals;
        std::sort(sorted.begin(), sorted.end());
        double q1 = sorted[n/4], q3 = sorted[(3*n)/4], iqr = q3 - q1;
        double thr = q3 + mltpl*iqr;
        for (size_t i=0;i<n;i++) {
            if (vals[i] > thr) {
                keep[records[idxs[i]].row_id] = 0;
            }
        }
    }

    std::vector<Record> out;
    out.reserve(records.size());
    for (size_t i=0;i<records.size();i++)
        if (keep[i]) out.push_back(records[i]);

    std::cout<<"cc_outl: Removed "<<(records.size()-out.size())
             <<" records.\n";
    return out;
}

// cc_cap (parallel)
std::vector<Record> cc_cap(const std::vector<Record> &records,
                           const std::vector<std::pair<double,double>> &caps,
                           double buffer = 10000.0,
                           bool verbose = true) {
    double buf_km = buffer / 1000.0;
    std::vector<Record> out;
    out.reserve(records.size());

    #pragma omp parallel
    {
        std::vector<Record> local;
        local.reserve(1024);

        #pragma omp for nowait
        for (int i = 0; i < static_cast<int>(records.size()); ++i) {
            const auto &r = records[i];
            bool keep = true;
            for (auto &c : caps) {
                if (haversine(r.lon,r.lat,c.first,c.second) <= buf_km) {
                    keep = false; break;
                }
            }
            if (keep) local.push_back(r);
        }

        #pragma omp critical
        out.insert(out.end(), local.begin(), local.end());
    }

    if (verbose)
        std::cout<<"cc_cap: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    return out;
}

// cc_gbif (parallel)
std::vector<Record> cc_gbif(const std::vector<Record> &records,
                            double buffer = 1000.0,
                            bool verbose = true) {
    double buf_km = buffer / 1000.0;
    double glon = 12.58, glat = 55.67;
    std::vector<Record> out;
    out.reserve(records.size());

    #pragma omp parallel
    {
        std::vector<Record> local;
        local.reserve(1024);

        #pragma omp for nowait
        for (int i = 0; i < static_cast<int>(records.size()); ++i) {
            const auto &r = records[i];
            if (haversine(r.lon,r.lat,glon,glat) > buf_km)
                local.push_back(r);
        }

        #pragma omp critical
        out.insert(out.end(), local.begin(), local.end());
    }

    if (verbose)
        std::cout<<"cc_gbif: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    return out;
}

// --------------------
// Raster Mask Helpers (for cc_sea_*)
// --------------------
struct RasterMask {
    std::vector<uint8_t> data;
    int xSize = 0, ySize = 0;
    double west = -180.0, east = 180.0, south = -90.0, north = 90.0;
    double res = 0.1; // degrees/pixel
};

static RasterMask rasterizeGeometriesToMask(
        const std::vector<OGRGeometry*>& geoms,
        double res = 0.1,
        double west = -180.0, double east = 180.0,
        double south = -90.0, double north = 90.0)
{
    RasterMask rm;
    rm.res = res; rm.west = west; rm.east = east; rm.south = south; rm.north = north;
    rm.xSize = static_cast<int>((east - west) / res);
    rm.ySize = static_cast<int>((north - south) / res);

    GDALDriver *memDrv = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!memDrv) throw std::runtime_error("MEM driver not available.");

    GDALDataset *maskDs = memDrv->Create("", rm.xSize, rm.ySize, 1, GDT_Byte, nullptr);
    if (!maskDs) throw std::runtime_error("Failed to create in-memory mask.");

    double gt[6] = {west, res, 0, north, 0, -res};
    maskDs->SetGeoTransform(gt);

    // Set WGS84 projection
    {
        OGRSpatialReference geoSRS;
        geoSRS.SetWellKnownGeogCS("WGS84");
        char *wkt = nullptr;
        geoSRS.exportToWkt(&wkt);
        maskDs->SetProjection(wkt);
        CPLFree(wkt);
    }

    // Prepare geometry handles
    std::vector<OGRGeometryH> hGeoms;
    hGeoms.reserve(geoms.size());
    for (auto *g : geoms) hGeoms.push_back((OGRGeometryH)g);

    // Burn value 1 into all pixels touched by geometries
    char **opts = nullptr;
    opts = CSLSetNameValue(opts, "ALL_TOUCHED", "TRUE");
    int bandList[1] = {1};
    double burnValues[1] = {1.0};

    CPLErr err = GDALRasterizeGeometries(
        maskDs,
        1, bandList,
        static_cast<int>(hGeoms.size()),
        hGeoms.data(),
        nullptr, nullptr,
        burnValues,
        nullptr, nullptr, opts);

    CSLDestroy(opts);
    if (err != CE_None) {
        GDALClose(maskDs);
        throw std::runtime_error("Rasterization failed.");
    }

    // Read mask into RAM
    rm.data.assign(rm.xSize * (size_t)rm.ySize, 0);
    CPLErr rc = maskDs->GetRasterBand(1)->RasterIO(
        GF_Read, 0, 0, rm.xSize, rm.ySize,
        rm.data.data(), rm.xSize, rm.ySize,
        GDT_Byte, 0, 0, nullptr);

    GDALClose(maskDs);
    if (rc != CE_None) throw std::runtime_error("Error reading mask data.");

    return rm;
}

static std::vector<Record> filterRecordsByMask(
        const std::vector<Record>& records,
        const RasterMask& rm)
{
    std::vector<Record> kept;
    kept.reserve(records.size());
    const int xSize = rm.xSize, ySize = rm.ySize;
    const double west = rm.west, north = rm.north, res = rm.res;

    #pragma omp parallel
    {
        std::vector<Record> local;
        local.reserve(1024);

        #pragma omp for nowait
        for (int i = 0; i < static_cast<int>(records.size()); ++i) {
            const auto &r = records[i];
            if (std::isnan(r.lon) || std::isnan(r.lat)) continue;
            int px = static_cast<int>((r.lon - west)  / res);
            int py = static_cast<int>((north - r.lat) / res);
            if (px < 0 || px >= xSize || py < 0 || py >= ySize) continue;
            if (rm.data[py * (size_t)xSize + px]) local.push_back(r);
        }

        #pragma omp critical
        kept.insert(kept.end(), local.begin(), local.end());
    }

    return kept;
}

// --------------------
// cc_sea_* (raster-mask versions; parallel lookup)
// --------------------

// cc_sea_gdal (mask)
std::vector<Record> cc_sea_gdal(const std::vector<Record> &records,
                                const std::string &landShp,
                                bool verbose = true)
{
    if (verbose) std::cout<<"cc_sea_gdal (raster mask): land test\n";
    std::vector<OGRGeometry*> landGeoms;
    {
        GDALDataset *ds = (GDALDataset*)GDALOpenEx(landShp.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
        if (!ds) { std::cerr << "Unable to open shapefile: " << landShp << "\n"; return records; }
        OGRLayer *layer = ds->GetLayer(0);
        if (!layer) { std::cerr << "No layer in shapefile: " << landShp << "\n"; GDALClose(ds); return records; }
        layer->ResetReading();
        OGRFeature *feat;
        while ((feat = layer->GetNextFeature())) {
            if (auto *g = feat->GetGeometryRef()) landGeoms.push_back(g->clone());
            OGRFeature::DestroyFeature(feat);
        }
        GDALClose(ds);
    }
    if (landGeoms.empty()) { std::cerr << "No land geometries loaded.\n"; return records; }

    std::vector<Record> out;
    try {
        auto mask = rasterizeGeometriesToMask(landGeoms, /*res*/0.1, -180, 180, -90, 90);
        out = filterRecordsByMask(records, mask);
    } catch (const std::exception &e) {
        std::cerr << "cc_sea_gdal error: " << e.what() << "\n";
        out = records; // fallback: no change
    }

    for (auto *g : landGeoms) OGRGeometryFactory::destroyGeometry(g);

    if (verbose)
        std::cout<<"cc_sea_gdal: Removed "<<(records.size()-out.size())<<" records.\n";
    return out;
}

// cc_sea_buffland (mask of buffered coastline)
std::vector<Record> cc_sea_buffland(const std::vector<Record> &records,
                                    const std::string &coastShp,
                                    bool verbose = true,
                                    double bufferDistDeg = 1.0)
{
    if (verbose) std::cout<<"cc_sea_buffland (raster mask): buffered coastline test\n";
    std::vector<OGRGeometry*> bufGeoms;
    {
        GDALDataset *ds = (GDALDataset*)GDALOpenEx(coastShp.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
        if (!ds) { std::cerr << "Unable to open shapefile: " << coastShp << "\n"; return records; }
        OGRLayer *layer = ds->GetLayer(0);
        if (!layer) { std::cerr << "No layer in shapefile: " << coastShp << "\n"; GDALClose(ds); return records; }
        layer->ResetReading();
        OGRFeature *feat;
        while ((feat = layer->GetNextFeature())) {
            if (auto *g = feat->GetGeometryRef()) {
                OGRGeometry* cloned = g->clone();
                OGRGeometry* buffered = cloned->Buffer(bufferDistDeg);
                OGRGeometryFactory::destroyGeometry(cloned);
                if (buffered) bufGeoms.push_back(buffered);
            }
            OGRFeature::DestroyFeature(feat);
        }
        GDALClose(ds);
    }
    if (bufGeoms.empty()) { std::cerr << "No buffered coastline geometries.\n"; return records; }

    std::vector<Record> out;
    try {
        auto mask = rasterizeGeometriesToMask(bufGeoms, /*res*/0.1, -180, 180, -90, 90);
        out = filterRecordsByMask(records, mask);
    } catch (const std::exception &e) {
        std::cerr << "cc_sea_buffland error: " << e.what() << "\n";
        out = records; // fallback
    }

    for (auto *g : bufGeoms) OGRGeometryFactory::destroyGeometry(g);

    if (verbose)
        std::cout<<"cc_sea_buffland: Removed "<<(records.size()-out.size())<<" records.\n";
    return out;
}

// --------------------
// CSV Loading
// --------------------
std::vector<std::string> split_line(const std::string &line, char d) {
    std::vector<std::string> tok; std::stringstream ss(line); std::string t;
    while (std::getline(ss,t,d)) tok.push_back(t);
    return tok;
}
std::vector<Record> loadCSV(const std::string &fn) {
    std::ifstream f(fn);
    if(!f.is_open()){ std::cerr<<"Cannot open "<<fn<<"\n"; return{}; }
    std::string hdr; std::getline(f,hdr);
    char d = hdr.find('\t')!=std::string::npos?'\t':',';
    auto cols = split_line(hdr,d);
    std::vector<std::string> lh;
    for(auto &c:cols){
        auto t=trim(c);
        std::transform(t.begin(),t.end(),t.begin(),::tolower);
        lh.push_back(t);
    }
    int ig=-1,is=-1,ic=-1,ilat=-1,ilon=-1,ie=-1;
    for(int i=0;i<(int)lh.size();++i){
        if(lh[i]=="gbifid")         ig=i;
        else if(lh[i]=="species")   is=i;
        else if(lh[i]=="countrycode")ic=i;
        else if(lh[i]=="decimallatitude") ilat=i;
        else if(lh[i]=="decimallongitude")ilon=i;
        else if(lh[i]=="eventdate") ie=i;
    }
    if(ig<0||is<0||ic<0||ilat<0||ilon<0||ie<0){
        std::cerr<<"Missing column\n";return{};
    }
    std::vector<Record> recs; int rid=0;
    std::string ln;
    while(std::getline(f,ln)){
        auto t=split_line(ln,d);
        if((int)t.size() < (int)lh.size()) continue;
        Record r;
        r.gbifID      = trim(t[ig]);
        r.species     = trim(t[is]);
        r.countryCode = upper(trim(t[ic]));
        try{ r.lat = std::stod(t[ilat]); } catch(...){ r.lat = NAN; }
        try{ r.lon = std::stod(t[ilon]); } catch(...){ r.lon = NAN; }
        r.eventDate = trim(t[ie]);
        r.row_id    = rid++;
        recs.push_back(r);
    }
    std::cout<<"loaded "<<recs.size()<<" records from csv.\n";
    return recs;
}

// --------------------
// Segfault Handler
// --------------------
void segfault_handler(int sig) {
    std::cerr<<"Segfault (signal "<<sig<<")\n";
    std::exit(sig);
}

// --------------------
// main()
// --------------------
int main() {
    std::signal(SIGSEGV, segfault_handler);
    GDALAllRegister();

    auto t0 = std::chrono::steady_clock::now();

    std::string csvFile    = "/Users/njord888/downloads/0023500-241107131044228.csv";
    std::string landShp    = "/Users/njord888/Downloads/ne_10m_land/ne_10m_land.shp";
    std::string coastShp   = "/Users/njord888/Downloads/ne_10m_coastline/ne_10m_coastline.shp";
    std::string urbShp     = "/Users/njord888/Downloads/ne_10m_urban_areas/ne_10m_urban_areas.shp";
    std::string countryShp = "/Users/njord888/Downloads/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp";

    auto records = loadCSV(csvFile);
    size_t initial = records.size();

    // Convert any ISO2 to ISO3 (cleaning process necessary for data integrity)
    auto iso2to3 = buildIso2To3Map(countryShp);
    convertIso2To3(records, iso2to3);

    // Count ISO3 codes before cleaning
    auto iso3Counts = countIso3(records);
    std::cout<<"ISO Alpha-3 countryCode counts:\n";
    for (auto &kv : iso3Counts) {
        std::cout<<"  "<<kv.first<<": "<<kv.second<<"\n";
    }

    {
        auto [r,d] = measureDuration(cc_val, records, true);
        cleaningSummaries.push_back({"cc_val", initial, r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_zero, records, true);
        cleaningSummaries.push_back({"cc_zero", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_dupl, records, true);
        cleaningSummaries.push_back({"cc_dupl", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_equ, records, "absolute", true);
        cleaningSummaries.push_back({"cc_equ", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    { // cc_inst (parallel)
        std::vector<std::pair<double,double>> inst = {{-77.0365,38.8977},{-0.1278,51.5074}};
        auto [r,d] = measureDuration(cc_inst, records, inst, 1000.0, true);
        cleaningSummaries.push_back({"cc_inst", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_coun, records, countryShp, "ISO_A3", true);
        cleaningSummaries.push_back({"cc_coun", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    { // cc_outl (parallel across species)
        auto [r,d] = measureDuration(cc_outl, records, "quantile",5.0,1000.0,7,false,0.0);
        cleaningSummaries.push_back({"cc_outl", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        std::vector<std::pair<double,double>> caps = {{-77.0369,38.9072},{2.3522,48.8566}};
        auto [r,d] = measureDuration(cc_cap, records, caps, 10000.0, true);
        cleaningSummaries.push_back({"cc_cap", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_gbif, records, 1000.0, true);
        cleaningSummaries.push_back({"cc_gbif", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_sea_gdal, records, landShp, true);
        cleaningSummaries.push_back({"cc_sea_gdal", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_sea_buffland, records, coastShp, true);
        cleaningSummaries.push_back({"cc_sea_buffland", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r,d] = measureDuration(cc_urb_gdal, records, urbShp, true);
        cleaningSummaries.push_back({"cc_urb_gdal", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> total = t1 - t0;

    std::cout<<"\n--- Cleaning Summary ---\n";
    for (auto &s : cleaningSummaries) {
        std::cout<<s.functionName<<": "
                 <<s.beforeCount<<"→"<<s.afterCount
                 <<" ("<<(s.beforeCount-s.afterCount)<<" removed)"
                 <<" in "<<s.durationSeconds<<"s\n";
    }
    std::cout<<"Overall: "<<initial<<"→"<<records.size()
             <<" in "<<total.count()<<"s\n";

    std::ofstream ofs("cleaned_data1.6.csv");
    ofs<<"gbifID,species,countryCode,decimalLongitude,decimalLatitude,eventDate\n";
    for (auto &r : records) {
        ofs<<r.gbifID<<","<<r.species<<","<<r.countryCode<<","
           <<r.lon<<","<<r.lat<<","<<r.eventDate<<"\n";
    }
    std::cout<<"Wrote cleaned_data1.6.csv\n";
    return 0;
}
