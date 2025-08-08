// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>______SEAS ONLY__(Improved4 with Rtrees from GEOS(!GDAL) )____<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Compile: g++ -std=c++17 main.cpp -o data_cleaner $(gdal-config --cflags) $(gdal-config --libs)
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

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "gdal_priv.h"
#include "ogrsf_frmts.h"

namespace fs = std::filesystem;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef bg::model::box<BoostPoint> BoostBox;
typedef std::pair<BoostBox, size_t> RTreeValue;

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
    if (verbose) std::cout << "cc_sea_rtree: fast land test using Boost R-tree\n";
    auto land = loadPolygons(shp);
    if (land.empty()) {
        std::cerr << "no land polygons loaded\n";
        return records;
    }

    bgi::rtree<RTreeValue, bgi::quadratic<16>> rtree;
    std::vector<OGRPreparedGeometry*> preparedGeoms;

    for (size_t i = 0; i < land.size(); ++i) {
        OGREnvelope env;
        land[i]->getEnvelope(&env);
        BoostBox bbox(BoostPoint(env.MinX, env.MinY), BoostPoint(env.MaxX, env.MaxY));
        rtree.insert(std::make_pair(bbox, i));
        preparedGeoms.push_back(OGRCreatePreparedGeometry(land[i]));
    }

    std::vector<Record> out;
    for (const auto &r : records) {
        BoostPoint pt(r.lon, r.lat);
        std::vector<RTreeValue> results;
        BoostBox queryBox(pt, pt);
        rtree.query(bgi::intersects(queryBox), std::back_inserter(results));

        bool onLand = false;
        OGRPoint ogrPt(r.lon, r.lat);
        for (const auto &val : results) {
            if (OGRPreparedGeometryIntersects(preparedGeoms[val.second], &ogrPt)) {
                onLand = true;
                break;
            }
        }
        if (onLand) out.push_back(r);  // Keep only if on land
    }

    for (auto *g : preparedGeoms) OGRDestroyPreparedGeometry(g);
    for (auto *g : land) OGRGeometryFactory::destroyGeometry(g);

    if (verbose)
        std::cout << "cc_sea_rtree: Removed " << (records.size() - out.size()) << " records.\n";
    return out;
}
