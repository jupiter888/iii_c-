#include "cc_seas_rtrees.h"
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <stdexcept>
#include <iostream>
#include <chrono>
#ifdef _OPENMP
  #include <omp.h>
#endif

static void force_wgs84(GDALDataset* ds) {
    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    ds->SetProjection(wkt);
    CPLFree(wkt);
}

RasterMask build_land_mask_from_shapefile(const std::string& shpPath, double res_deg) {
    if (res_deg <= 0.0) throw std::invalid_argument("Resolution must be > 0.");
    GDALAllRegister();
    OGRRegisterAll();

    GDALDataset* shpDs = static_cast<GDALDataset*>(
        GDALOpenEx(shpPath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!shpDs) throw std::runtime_error("Unable to open land shapefile: " + shpPath);

    OGRLayerH layerH = OGR_DS_GetLayer(shpDs, 0);
    if (!layerH) { GDALClose(shpDs); throw std::runtime_error("Failed to get layer 0."); }

    RasterMask mask;
    mask.res = res_deg;
    mask.xSize = static_cast<int>((mask.east  - mask.west)/mask.res);
    mask.ySize = static_cast<int>((mask.north - mask.south)/mask.res);

    GDALDriver* memDrv = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!memDrv) { GDALClose(shpDs); throw std::runtime_error("GDAL MEM driver not available."); }

    GDALDataset* maskDs = memDrv->Create("", mask.xSize, mask.ySize, 1, GDT_Byte, nullptr);
    if (!maskDs) { GDALClose(shpDs); throw std::runtime_error("Failed to create in-memory raster."); }

    double gt[6] = { mask.west, mask.res, 0.0, mask.north, 0.0, -mask.res };
    if (maskDs->SetGeoTransform(gt) != CE_None) {
        GDALClose(maskDs); GDALClose(shpDs);
        throw std::runtime_error("Failed to set geotransform.");
    }

    force_wgs84(maskDs);

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "ALL_TOUCHED", "TRUE");
    const int bandList[1] = {1};
    OGRLayerH layers[1] = {layerH};
    const double burnValues[1] = {1.0};

    const CPLErr err = GDALRasterizeLayers(
        maskDs, 1, const_cast<int*>(bandList),
        1, layers, nullptr, nullptr,
        const_cast<double*>(burnValues), opts, nullptr, nullptr);
    CSLDestroy(opts);
    GDALClose(shpDs);
    if (err != CE_None) { GDALClose(maskDs); throw std::runtime_error("Rasterization failed."); }

    mask.data.assign(mask.xSize * (size_t)mask.ySize, 0);
    if (maskDs->GetRasterBand(1)->RasterIO(
            GF_Read, 0, 0, mask.xSize, mask.ySize,
            mask.data.data(), mask.xSize, mask.ySize,
            GDT_Byte, 0, 0, nullptr) != CE_None) {
        GDALClose(maskDs);
        throw std::runtime_error("Failed to read mask data.");
    }
    GDALClose(maskDs);
    return mask;
}

size_t cc_seas_rtrees(const std::vector<Record>& in,
                      std::vector<Record>& out,
                      const RasterMask& mask)
{
    const auto t0 = std::chrono::steady_clock::now();
    out.clear();
    out.reserve(in.size());

#ifdef _OPENMP
    #pragma omp parallel
    {
        std::vector<Record> local;
        local.reserve(in.size() / (omp_get_num_threads() ? omp_get_num_threads() : 1) + 1);

        #pragma omp for nowait
        for (int i = 0; i < (int)in.size(); ++i) {
            const auto& r = in[i];
            if (mask.is_land(r.lon, r.lat)) local.push_back(r);
        }
        #pragma omp critical
        out.insert(out.end(), local.begin(), local.end());
    }
#else
    for (const auto& r : in) if (mask.is_land(r.lon, r.lat)) out.push_back(r);
#endif

    const auto t1 = std::chrono::steady_clock::now();
    std::cout << "[cc_seas_rtrees] Kept " << out.size() << " of " << in.size()
              << " (land points) in "
              << std::chrono::duration<double>(t1 - t0).count() << " s.\n";
    return out.size();
}
