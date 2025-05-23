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

// GDAL/OGR headers
#include "gdal_priv.h"
#include "ogrsf_frmts.h"

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
    double      lon;
    double      lat;
    std::string eventDate;
    int         row_id;
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
    for (auto &r : records) {
        if (!(r.lon == 0.0 && r.lat == 0.0))
            out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_zero: Removed "<<(records.size()-out.size())
                 <<" records (zero coords).\n";
    return out;
}

// cc_dupl with 10 m jitter + exact eventDate
std::vector<Record> cc_dupl(const std::vector<Record> &records,
                            bool verbose = true) {
    const double jitter_deg = 10.0 / 111320.0;
    std::unordered_set<std::string> seen;
    std::vector<Record> out;
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

// cc_inst
std::vector<Record> cc_inst(const std::vector<Record> &records,
                            const std::vector<std::pair<double,double>> &inst,
                            double buffer = 1000.0,
                            bool verbose = true) {
    double buf_km = buffer / 1000.0;
    std::vector<Record> out;
    for (auto &r : records) {
        bool keep = true;
        for (auto &i : inst) {
            if (haversine(r.lon,r.lat,i.first,i.second) <= buf_km) {
                keep = false; break;
            }
        }
        if (keep) out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_inst: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    return out;
}

// cc_coun with debug
std::vector<Record> cc_coun(const std::vector<Record> &records,
                            const std::string &shapefilePath,
                            const std::string &isoField = "ISO_A3",
                            bool verbose = true) {
    auto feats = loadCountryFeaturesWithCodes(shapefilePath, "ISO_A2", isoField);
    std::vector<Record> out;
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

// cc_outl per-species
std::vector<Record> cc_outl(const std::vector<Record> &records,
                            const std::string &method = "quantile",
                            double mltpl = 5.0,
                            double tdi = 1000.0,
                            int min_occs = 7,
                            bool        = false,
                            double      = 0.0) {
    std::unordered_map<std::string,std::vector<Record>> speciesMap;
    for (auto &r : records) speciesMap[upper(r.species)].push_back(r);
    std::vector<bool> keep(records.size(), true);
    size_t idx = 0;
    for (auto &kv : speciesMap) {
        ++idx;
        auto &grp = kv.second;
        if (grp.size() < (size_t)min_occs) {
            std::cerr<<"warning: species "<<kv.first
                     <<" <"<<min_occs<<"; skipping\n";
            continue;
        }
        size_t n = grp.size();
        std::vector<double> vals(n, 0.0);
        for (size_t i=0;i<n;i++) {
            double minD = std::numeric_limits<double>::infinity();
            for (size_t j=0;j<n;j++) {
                if (i==j) continue;
                double d = haversine(grp[i].lon,grp[i].lat,
                                     grp[j].lon,grp[j].lat);
                if (d<minD) minD = d;
            }
            vals[i] = minD;
        }
        std::vector<double> sorted = vals;
        std::sort(sorted.begin(), sorted.end());
        double q1 = sorted[n/4], q3 = sorted[(3*n)/4], iqr = q3 - q1;
        double thr = q3 + mltpl*iqr;
        for (size_t i=0;i<n;i++) {
            if (vals[i] > thr)
                keep[grp[i].row_id] = false;
        }
        std::cout<<"Processed outlier "<<idx<<"/"<<speciesMap.size()
                 <<"\r"<<std::flush;
    }
    std::vector<Record> out;
    for (size_t i=0;i<records.size();i++)
        if (keep[i]) out.push_back(records[i]);
    std::cout<<"\ncc_outl: Removed "<<(records.size()-out.size())
             <<" records.\n";
    return out;
}

// cc_cap
std::vector<Record> cc_cap(const std::vector<Record> &records,
                           const std::vector<std::pair<double,double>> &caps,
                           double buffer = 10000.0,
                           bool verbose = true) {
    double buf_km = buffer / 1000.0;
    std::vector<Record> out;
    for (auto &r : records) {
        bool keep = true;
        for (auto &c : caps) {
            if (haversine(r.lon,r.lat,c.first,c.second) <= buf_km) {
                keep = false; break;
            }
        }
        if (keep) out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_cap: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    return out;
}

// cc_gbif
std::vector<Record> cc_gbif(const std::vector<Record> &records,
                            double buffer = 1000.0,
                            bool verbose = true) {
    double buf_km = buffer / 1000.0;
    double glon = 12.58, glat = 55.67;
    std::vector<Record> out;
    for (auto &r : records) {
        if (haversine(r.lon,r.lat,glon,glat) > buf_km)
            out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_gbif: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    return out;
}

// cc_sea_gdal
std::vector<Record> cc_sea_gdal(const std::vector<Record> &records,
                                const std::string &shp,
                                bool verbose = true) {
    if (verbose) std::cout<<"cc_sea_gdal: land test\n";
    auto land = loadPolygons(shp);
    if (land.empty()) { std::cerr<<"no land loaded\n"; return records; }
    std::vector<Record> out;
    for (auto &r : records) {
        OGRPoint pt(r.lon,r.lat);
        bool onLand = false;
        for (auto *g : land) {
            if (pt.Within(g)) { onLand = true; break; }
        }
        if (onLand) out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_sea_gdal: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    for (auto *g : land) OGRGeometryFactory::destroyGeometry(g);
    return out;
}

// cc_sea_buffland
std::vector<Record> cc_sea_buffland(const std::vector<Record> &records,
                                    const std::string &shp,
                                    bool verbose = true) {
    if (verbose) std::cout<<"cc_sea_buffland: buffered test\n";
    auto coast = loadPolygons(shp);
    if (coast.empty()) { std::cerr<<"no coast loaded\n"; return records; }
    std::vector<OGRGeometry*> buf;
    for (auto *g : coast) buf.push_back(g->Buffer(1.0));
    for (auto *g : coast) OGRGeometryFactory::destroyGeometry(g);
    std::vector<Record> out;
    for (auto &r : records) {
        OGRPoint pt(r.lon,r.lat);
        bool onLand = false;
        for (auto *g : buf) {
            if (pt.Within(g)) { onLand = true; break; }
        }
        if (onLand) out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_sea_buffland: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    for (auto *g : buf) OGRGeometryFactory::destroyGeometry(g);
    return out;
}

// cc_urb_gdal
std::vector<Record> cc_urb_gdal(const std::vector<Record> &records,
                                const std::string &shp,
                                bool verbose = true) {
    if (verbose) std::cout<<"cc_urb_gdal: urban test\n";
    auto urb = loadPolygons(shp);
    if (urb.empty()) { std::cerr<<"no urban loaded\n"; return records; }
    std::vector<Record> out;
    for (auto &r : records) {
        OGRPoint pt(r.lon,r.lat);
        bool inU = false;
        for (auto *g : urb) {
            if (pt.Within(g)) { inU = true; break; }
        }
        if (!inU) out.push_back(r);
    }
    if (verbose)
        std::cout<<"cc_urb_gdal: Removed "<<(records.size()-out.size())
                 <<" records.\n";
    for (auto *g : urb) OGRGeometryFactory::destroyGeometry(g);
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
    std::string countryShp = "/Users/njord888/Downloads/ne_110m_admin_0_countries/"
                             "ne_110m_admin_0_countries.shp";

    auto records = loadCSV(csvFile);
    size_t initial = records.size();

    // Convert any ISO2 to ISO3
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
    {
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
