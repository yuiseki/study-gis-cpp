#include <iostream>
#include <string>

#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "cpl_conv.h" // CPLFree

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: 004 <path_to_geojson>\n";
        return 2;
    }
    const std::string path = argv[1];

    GDALAllRegister();

    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)
    );
    if (!ds) {
        std::cerr << "Failed to open dataset: " << path << "\n";
        return 1;
    }

    const int layer_count = ds->GetLayerCount();
    std::cout << "dataset: " << path << "\n";
    std::cout << "layers : " << layer_count << "\n\n";

    for (int i = 0; i < layer_count; i++) {
        OGRLayer* layer = ds->GetLayer(i);
        if (!layer) continue;

        std::cout << "== Layer[" << i << "] ==\n";
        std::cout << "name        : " << layer->GetName() << "\n";
        std::cout << "featureCount : " << layer->GetFeatureCount(TRUE) << "\n";

        const OGRwkbGeometryType gt = layer->GetGeomType();
        std::cout << "geomType    : " << OGRGeometryTypeToName(gt) << "\n";

        // SRS (WKT)
        const OGRSpatialReference* srs = layer->GetSpatialRef(); // owned by layer :contentReference[oaicite:5]{index=5}
        if (srs) {
            char* wkt = nullptr;
            srs->exportToPrettyWkt(&wkt);
            std::cout << "SRS(WKT)    :\n" << (wkt ? wkt : "(null)") << "\n";
            CPLFree(wkt);
        } else {
            std::cout << "SRS(WKT)    : (none)\n";
        }

        // Extent (Envelope)
        OGREnvelope env;
        if (layer->GetExtent(&env, TRUE) == OGRERR_NONE) {
            std::cout << "extent      : "
                      << "minX=" << env.MinX << ", minY=" << env.MinY
                      << ", maxX=" << env.MaxX << ", maxY=" << env.MaxY
                      << "\n";
        } else {
            std::cout << "extent      : (unknown)\n";
        }

        std::cout << "\n";
    }

    GDALClose(ds);
    return 0;
}
