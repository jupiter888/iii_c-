#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
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

#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "gdal_alg.h" // For RTree API

namespace fs = std::filesystem;

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

struct Record {
    std::string gbifID;
    std::string species;
    std::string countryCode;
    double      lon;
    double      lat;
    std::string eventDate;
    int         row_id;
};

struct CleaningSummary {
    std::string functionName;
    size_t      beforeCount;
    size_t      afterCount;
    bool        skipped;
    double      durationSeconds;
};
std::vector<CleaningSummary> cleaningSummaries;

std::vector<OGRGeometry*> loadPolygons(const std::string &shapefilePath) {
    std::vector<OGRGeometry*> polygons;
    GDALDataset *ds = (GDALDataset*)GDALOpenEx(shapefilePath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
    if (!ds) {
        std::cerr << "Failed to open shapefile: " << shapefilePath << "\n";
        return polygons;
    }
    OGRLayer *layer = ds->GetLayer(0);
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

template<typename F, typename... Args>
auto measureDuration(F func, Args&&... args)
    -> std::pair<decltype(func(std::forward<Args>(args)...)), double> {
    auto start  = std::chrono::steady_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end    = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return {result, elapsed.count()};
}

std::vector<Record> reassign_row_ids(const std::vector<Record> &records) {
    std::vector<Record> out = records;
    for (size_t i = 0; i < out.size(); ++i) {
        out[i].row_id = static_cast<int>(i);
    }
    return out;
}

std::vector<Record> cc_sea_gdal(const std::vector<Record> &records,
                                const std::string &shp,
                                bool verbose = true) {
    if (verbose) std::cout << "cc_sea_gdal: land test\n";
    auto land = loadPolygons(shp);
    if (land.empty()) {
        std::cerr << "no land loaded\n";
        return records;
    }
    std::vector<Record> out;
    for (auto &r : records) {
        OGRPoint pt(r.lon, r.lat);
        bool onLand = false;
        for (auto *g : land) {
            if (pt.Within(g)) {
                onLand = true;
                break;
            }
        }
        if (onLand) out.push_back(r);
    }
    if (verbose)
        std::cout << "cc_sea_gdal: Removed " << (records.size() - out.size()) << " records.\n";
    for (auto *g : land) OGRGeometryFactory::destroyGeometry(g);
    return out;
}

std::vector<Record> cc_sea_rtree(const std::vector<Record> &records,
                                 const std::string &shp,
                                 bool verbose = true) {
    if (verbose) std::cout << "cc_sea_rtree: high-performance coastline water removal\n";

    auto coast = loadPolygons(shp);
    if (coast.empty()) {
        std::cerr << "no coast loaded\n";
        return records;
    }

    const double simplifyTolerance = 0.01;   // Reduce geometry detail
    const double bufferDistance = 0.1;       // ~11 km buffer

    std::vector<OGRGeometry*> bufferedGeoms;
    std::vector<OGRPreparedGeometry*> preparedGeoms;

    void* hRTree = CPLRTreeCreate(8);

    for (size_t i = 0; i < coast.size(); ++i) {
        OGRGeometry* simplified = coast[i]->Simplify(simplifyTolerance);
        OGRGeometry* buffered = simplified->Buffer(bufferDistance);

        if (buffered) {
            OGREnvelope env;
            buffered->getEnvelope(&env);

            double minmax[4] = { env.MinX, env.MinY, env.MaxX, env.MaxY };
            CPLRTreeInsert(hRTree, minmax, (int)i);

            bufferedGeoms.push_back(buffered);
            preparedGeoms.push_back(OGRCreatePreparedGeometry(buffered));
        }

        OGRGeometryFactory::destroyGeometry(simplified);
        OGRGeometryFactory::destroyGeometry(coast[i]);
    }
    coast.clear();

    std::vector<Record> out;
    out.reserve(records.size());

    int hits[256];

    for (const auto &r : records) {
        double ptEnv[4] = { r.lon, r.lat, r.lon, r.lat };
        int numHits = CPLRTreeSearch(hRTree, ptEnv, hits, 256);

        bool keep = false;
        if (numHits > 0) {
            OGRPoint pt(r.lon, r.lat);
            for (int j = 0; j < numHits; ++j) {
                if (OGRPreparedGeometryIntersects(preparedGeoms[hits[j]], &pt)) {
                    keep = true;
                    break;
                }
            }
        }

        if (keep) out.push_back(r);
    }

    CPLRTreeDestroy(hRTree);
    for (auto *g : bufferedGeoms) OGRGeometryFactory::destroyGeometry(g);
    for (auto *pg : preparedGeoms) OGRDestroyPreparedGeometry(pg);

    if (verbose)
        std::cout << "cc_sea_rtree: Removed " << (records.size() - out.size()) << " records.\n";

    return out;
}

std::vector<std::string> split_line(const std::string &line, char d) {
    std::vector<std::string> tok;
    std::stringstream ss(line);
    std::string t;
    while (std::getline(ss, t, d)) tok.push_back(t);
    return tok;
}

std::vector<Record> loadCSV(const std::string &fn) {
    std::ifstream f(fn);
    if (!f.is_open()) {
        std::cerr << "Cannot open " << fn << "\n";
        return {};
    }
    std::string hdr;
    std::getline(f, hdr);
    char d = hdr.find('\t') != std::string::npos ? '\t' : ',';
    auto cols = split_line(hdr, d);
    std::vector<std::string> lh;
    for (auto &c : cols) {
        auto t = trim(c);
        std::transform(t.begin(), t.end(), t.begin(), ::tolower);
        lh.push_back(t);
    }
    int ig = -1, is = -1, ic = -1, ilat = -1, ilon = -1, ie = -1;
    for (int i = 0; i < (int)lh.size(); ++i) {
        if (lh[i] == "gbifid") ig = i;
        else if (lh[i] == "species") is = i;
        else if (lh[i] == "countrycode") ic = i;
        else if (lh[i] == "decimallatitude") ilat = i;
        else if (lh[i] == "decimallongitude") ilon = i;
        else if (lh[i] == "eventdate") ie = i;
    }
    if (ig < 0 || is < 0 || ic < 0 || ilat < 0 || ilon < 0 || ie < 0) {
        std::cerr << "Missing column\n";
        return {};
    }
    std::vector<Record> recs;
    int rid = 0;
    std::string ln;
    while (std::getline(f, ln)) {
        auto t = split_line(ln, d);
        if ((int)t.size() < (int)lh.size()) continue;
        Record r;
        r.gbifID      = trim(t[ig]);
        r.species     = trim(t[is]);
        r.countryCode = upper(trim(t[ic]));
        try { r.lat = std::stod(t[ilat]); } catch (...) { r.lat = NAN; }
        try { r.lon = std::stod(t[ilon]); } catch (...) { r.lon = NAN; }
        r.eventDate = trim(t[ie]);
        r.row_id    = rid++;
        recs.push_back(r);
    }
    std::cout << "loaded " << recs.size() << " records from csv.\n";
    return recs;
}

void segfault_handler(int sig) {
    std::cerr << "Segfault (signal " << sig << ")\n";
    std::exit(sig);
}

int main() {
    std::signal(SIGSEGV, segfault_handler);
    GDALAllRegister();

    auto t0 = std::chrono::steady_clock::now();

    std::string csvFile    = "/Users/njord888/downloads/10k_germany.csv";
    std::string landShp    = "/Users/njord888/Downloads/ne_10m_land/ne_10m_land.shp";
    std::string coastShp   = "/Users/njord888/Downloads/ne_10m_coastline/ne_10m_coastline.shp";

    auto records = loadCSV(csvFile);
    size_t initial = records.size();

    {
        auto [r, d] = measureDuration(cc_sea_gdal, records, landShp, true);
        cleaningSummaries.push_back({"cc_sea_gdal", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }
    {
        auto [r, d] = measureDuration(cc_sea_rtree, records, coastShp, true);
        cleaningSummaries.push_back({"cc_sea_rtree", records.size(), r.size(), false, d});
        records = reassign_row_ids(r);
    }

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> total = t1 - t0;

    std::cout << "\n--- Cleaning Summary ---\n";
    for (auto &s : cleaningSummaries) {
        std::cout << s.functionName << ": "
                  << s.beforeCount << "→" << s.afterCount
                  << " (" << (s.beforeCount - s.afterCount) << " removed)"
                  << " in " << s.durationSeconds << "s\n";
    }
    std::cout << "Overall: " << initial << "→" << records.size()
              << " in " << total.count() << "s\n";

    std::ofstream ofs("rtree_seas_cleaned.csv");
    ofs << "gbifID,species,countryCode,decimalLongitude,decimalLatitude,eventDate\n";
    for (auto &r : records) {
        ofs << r.gbifID << "," << r.species << "," << r.countryCode << ","
            << r.lon << "," << r.lat << "," << r.eventDate << "\n";
    }
    std::cout << "Wrote rtree_seas_cleaned.csv\n";
    return 0;
}
