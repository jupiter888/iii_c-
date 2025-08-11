#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct Record {
    std::string gbifID, species, countryCode, eventDate;
    double lon, lat;
};

struct RasterMask {
    double west = -180.0, east = 180.0, south = -90.0, north = 90.0;
    double res  = 0.1; // degrees
    int xSize = 0, ySize = 0;
    std::vector<uint8_t> data; // 1 = land, 0 = sea

    inline bool is_land(double lon, double lat) const {
        if (!(lon == lon) || !(lat == lat)) return false; // NaN check
        if (lon < west || lon >= east || lat < south || lat >= north) return false;
        const int px = int((lon - west) / res);
        const int py = int((north - lat) / res);
        if (px < 0 || px >= xSize || py < 0 || py >= ySize) return false;
        return data[py * (size_t)xSize + px] != 0;
    }
};

// Build once per run (expects land polygons in WGS84).
RasterMask build_land_mask_from_shapefile(const std::string& shpPath, double res_deg = 0.1);

// Keep only records on land (i.e., removes sea/ocean points).
size_t cc_seas_rtrees(const std::vector<Record>& in,
                      std::vector<Record>& out,
                      const RasterMask& mask);
