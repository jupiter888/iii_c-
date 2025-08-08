// ____Using Homebrew’s GCC 14 (native OpenMP support)____
// brew list libomp                                       
// brew list gcc                               
// brew info  libomp                                  
//  brew info  gcc

// g++-14 -std=c++17 main.cpp -o data_cleaner \\n    $(gdal-config --cflags) $(gdal-config --libs) \\n    -fopenmp\n
// caffeinate -i ./data_cleaner \\n  ~/Downloads/0023500-241107131044228.csv \\n  ~/Downloads/ne_10m_land/ne_10m_land.shp\n
// main.cpp
#include <gdal_priv.h>
#include <gdal_alg.h>
#include <ogrsf_frmts.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <algorithm>

struct Record {
    std::string gbifID, species, countryCode, eventDate;
    double lon = NAN, lat = NAN;
};

static std::string trim(const std::string &s) {
    size_t b = s.find_first_not_of(" \t\n\r\"");
    size_t e = s.find_last_not_of(" \t\n\r\"");
    return (b == std::string::npos || e == std::string::npos)
               ? ""
               : s.substr(b, e - b + 1);
}

static std::string upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

std::vector<Record> loadCSV(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open CSV file: " + path);

    std::string hdr;
    std::getline(in, hdr);
    char delim = (hdr.find('\t') != std::string::npos ? '\t' : ',');

    std::vector<std::string> cols;
    {
        std::stringstream ss(hdr);
        std::string c;
        while (std::getline(ss, c, delim))
            cols.push_back(trim(c));
    }

    int ig=-1, is=-1, ic=-1, ilat=-1, ilon=-1, ie=-1;
    for (int i = 0; i < (int)cols.size(); ++i) {
        std::string name = cols[i];
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        if (name == "gbifid")         ig = i;
        else if (name == "species")   is = i;
        else if (name == "countrycode") ic = i;
        else if (name == "decimallatitude")  ilat = i;
        else if (name == "decimallongitude") ilon = i;
        else if (name == "eventdate") ie = i;
    }
    if (ig<0||is<0||ic<0||ilat<0||ilon<0||ie<0)
        throw std::runtime_error("CSV missing required column.");

    std::vector<Record> recs;
    recs.reserve(10000);
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> f;
        std::stringstream ls(line);
        std::string tok;
        while (std::getline(ls, tok, delim))
            f.push_back(tok);
        if (f.size() < cols.size()) continue;
        Record r;
        r.gbifID      = trim(f[ig]);
        r.species     = trim(f[is]);
        r.countryCode = upper(trim(f[ic]));
        try { r.lat = std::stod(f[ilat]); } catch (...) {}
        try { r.lon = std::stod(f[ilon]); } catch (...) {}
        r.eventDate   = trim(f[ie]);
        recs.push_back(r);
    }
    std::cout << "Loaded " << recs.size() << " records.\n";
    return recs;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <points.csv> <land_shapefile.shp>\n";
        return 1;
    }
    const char *csvPath = argv[1];
    const char *shpPath = argv[2];

    GDALAllRegister();
    OGRRegisterAll();

    auto t0 = std::chrono::steady_clock::now();
    std::vector<Record> records;
    try {
        records = loadCSV(csvPath);
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    std::cout << "Opening shapefile: " << shpPath << "\n";
    GDALDataset *shpDs = static_cast<GDALDataset*>(
        GDALOpenEx(shpPath, GDAL_OF_VECTOR, nullptr, nullptr, nullptr)
    );
    if (!shpDs) {
        std::cerr << "Unable to open shapefile: " << shpPath << "\n";
        return 1;
    }
    std::cout << "Shapefile opened.\n";

    // C API layer handle
    OGRLayerH layerH = OGR_DS_GetLayer(shpDs, 0);
    if (!layerH) {
        std::cerr << "Failed to fetch layer from shapefile.\n";
        GDALClose(shpDs);
        return 1;
    }

    // Grid parameters
    const double res = 0.1;
    const double west=-180.0, east=180.0, south=-90.0, north=90.0;
    const int xSize = int((east - west)/res);
    const int ySize = int((north - south)/res);

    std::cout << "Creating mask raster (" << xSize << "x" << ySize << ")...\n";
    GDALDriver *memDrv = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!memDrv) {
        std::cerr << "MEM driver not available.\n";
        GDALClose(shpDs);
        return 1;
    }
    GDALDataset *maskDs = memDrv->Create("", xSize, ySize, 1, GDT_Byte, nullptr);
    if (!maskDs) {
        std::cerr << "Failed to create in-memory mask.\n";
        GDALClose(shpDs);
        return 1;
    }
    double gt[6] = {west, res, 0, north, 0, -res};
    maskDs->SetGeoTransform(gt);

    // Force WGS84 even if shapefile .prj is missing
    {
        OGRSpatialReference geoSRS;
        geoSRS.SetWellKnownGeogCS("WGS84");
        char *projWKT = nullptr;
        geoSRS.exportToWkt(&projWKT);
        maskDs->SetProjection(projWKT);
        CPLFree(projWKT);
    }

    std::cout << "Rasterizing land layer into mask...\n";
    char **opts = nullptr;
    opts = CSLSetNameValue(opts, "ALL_TOUCHED", "TRUE");
    int bandList[1]      = {1};
    OGRLayerH layers[1]  = {layerH};
    double burnValues[1] = {1.0};

    CPLErr err = GDALRasterizeLayers(
        maskDs,
        1,          // one band
        bandList,
        1,          // one layer
        layers,
        nullptr,    // no geo-transformer
        nullptr,    // no transformer args
        burnValues, // burn 1.0 where polygons exist
        opts,
        nullptr,    // no progress callback
        nullptr
    );
    CSLDestroy(opts);
    GDALClose(shpDs);
    if (err != CE_None) {
        std::cerr << "Rasterization failed.\n";
        GDALClose(maskDs);
        return 1;
    }
    std::cout << "Rasterization complete.\n";

    std::cout << "Reading mask into memory...\n";
    std::vector<uint8_t> mask(xSize * (size_t)ySize);
    {
        CPLErr rc = maskDs->GetRasterBand(1)->RasterIO(
            GF_Read, 0, 0, xSize, ySize,
            mask.data(), xSize, ySize,
            GDT_Byte, 0, 0, nullptr
        );
        if (rc != CE_None) {
            std::cerr << "Error reading mask data.\n";
            GDALClose(maskDs);
            return 1;
        }
    }
    GDALClose(maskDs);
    std::cout << "Mask loaded.\n";

    std::cout << "Filtering points...\n";
    std::vector<Record> kept;
    kept.reserve(records.size());
    #pragma omp parallel
    {
        std::vector<Record> local;
        #pragma omp for nowait
        for (int i = 0; i < (int)records.size(); ++i) {
            const auto &r = records[i];
            if (std::isnan(r.lon) || std::isnan(r.lat)) continue;
            int px = int((r.lon - west)  / res);
            int py = int((north - r.lat) / res);
            if (px < 0 || px >= xSize || py < 0 || py >= ySize) continue;
            if (mask[py * (size_t)xSize + px])
                local.push_back(r);
        }
        #pragma omp critical
        kept.insert(kept.end(), local.begin(), local.end());
    }
    std::cout << "Filtering complete.\n";

    auto t1 = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Kept " << kept.size() << " of "
              << records.size() << " in "
              << dt << " seconds.\n";

    std::ofstream out("filtered_land.csv");
    out << "gbifID,species,countryCode,decimalLongitude,decimalLatitude,eventDate\n";
    for (auto &r : kept) {
        out << r.gbifID << ',' << r.species << ','
            << r.countryCode << ',' << r.lon << ','
            << r.lat << ',' << r.eventDate << '\n';
    }
    std::cout << "Wrote filtered_land88.1313.csv\n";
    return 0;
}
