#include <cmath>
#include <iostream>
#include <string>

#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "cpl_vsi.h"
#include "cpl_conv.h"   // CPLFree

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/PolygonArea.hpp>

#include <nlohmann/json.hpp>
using nlohmann::json;

struct Metrics {
    bool ok = false;
    double area_m2 = 0.0;               // outer - holes (abs-based)
    double perimeter_outer_m = 0.0;     // outer ring only
    double perimeter_all_rings_m = 0.0; // outer + holes
};

static void die(const std::string& msg, int code = 2) {
    std::cerr << msg << "\n";
    std::exit(code);
}

static double wrap_lon(double lon_deg) {
    lon_deg = std::fmod(lon_deg + 180.0, 360.0);
    if (lon_deg < 0) lon_deg += 360.0;
    return lon_deg - 180.0;
}

static bool nearly_equal(double a, double b, double eps = 1e-12) {
    return std::abs(a - b) <= eps;
}

static Metrics ring_metrics(
    const GeographicLib::Geodesic& geod,
    const OGRLinearRing* ring
) {
    Metrics m;
    if (!ring) return m;

    const int n = ring->getNumPoints();
    if (n < 4) return m; // LinearRingは通常「閉じた」座標列

    // 末尾が先頭と同一点なら重複を避ける
    int useN = n;
    if (n >= 2 &&
        nearly_equal(ring->getX(0), ring->getX(n - 1)) &&
        nearly_equal(ring->getY(0), ring->getY(n - 1))) {
        useN = n - 1;
    }
    if (useN < 3) return m;

    GeographicLib::PolygonArea poly(geod);
    for (int i = 0; i < useN; i++) {
        const double lon = wrap_lon(ring->getX(i)); // GeoJSON: x=lon :contentReference[oaicite:6]{index=6}
        const double lat = ring->getY(i);           // GeoJSON: y=lat :contentReference[oaicite:7]{index=7}
        poly.AddPoint(lat, lon);                    // AddPoint(lat, lon) :contentReference[oaicite:8]{index=8}
    }

    double per = 0.0, area = 0.0;
    poly.Compute(false, true, per, area);           // 周長は閉曲線を含む :contentReference[oaicite:9]{index=9}

    m.ok = true;
    m.perimeter_outer_m = per;      // 呼び出し側で使い分け
    m.area_m2 = std::abs(area);
    return m;
}

static Metrics polygon_metrics(
    const GeographicLib::Geodesic& geod,
    const OGRPolygon* poly
) {
    Metrics out;
    if (!poly) return out;

    const OGRLinearRing* ext = poly->getExteriorRing();
    Metrics extm = ring_metrics(geod, ext);
    if (!extm.ok) return out;

    double area = extm.area_m2;
    double per_outer = extm.perimeter_outer_m;
    double per_all = extm.perimeter_outer_m;

    const int holes = poly->getNumInteriorRings();
    for (int i = 0; i < holes; i++) {
        const OGRLinearRing* in = poly->getInteriorRing(i);
        Metrics hm = ring_metrics(geod, in);
        if (!hm.ok) continue;
        area -= hm.area_m2;
        per_all += hm.perimeter_outer_m;
    }

    if (area < 0) area = 0.0;

    out.ok = true;
    out.area_m2 = area;
    out.perimeter_outer_m = per_outer;
    out.perimeter_all_rings_m = per_all;
    return out;
}

static Metrics geometry_metrics(
    const GeographicLib::Geodesic& geod,
    const OGRGeometry* g
) {
    Metrics out;
    if (!g) return out;

    const auto t = wkbFlatten(g->getGeometryType());
    if (t == wkbPolygon) {
        return polygon_metrics(geod, g->toPolygon());
    }
    if (t == wkbMultiPolygon) {
        const auto* mp = g->toMultiPolygon();
        if (!mp) return out;

        double area = 0.0, per_outer = 0.0, per_all = 0.0;
        for (int i = 0; i < mp->getNumGeometries(); i++) {
            const OGRGeometry* gi = mp->getGeometryRef(i);
            Metrics pm = polygon_metrics(geod, gi ? gi->toPolygon() : nullptr);
            if (!pm.ok) continue;
            area += pm.area_m2;
            per_outer += pm.perimeter_outer_m;
            per_all += pm.perimeter_all_rings_m;
        }
        out.ok = true;
        out.area_m2 = area;
        out.perimeter_outer_m = per_outer;
        out.perimeter_all_rings_m = per_all;
        return out;
    }

    // Polygon/MultiPolygon以外は対象外
    return out;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: 005 <input.geojson> [output.geojson|-]\n";
        return 2;
    }
    const std::string inPath = argv[1];
    const std::string outPath = (argc >= 3) ? argv[2] : "-";

    GDALAllRegister();

    GDALDataset* inDs = static_cast<GDALDataset*>(
        GDALOpenEx(inPath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)
    );
    if (!inDs) die("Failed to open input: " + inPath);

    // WGS84楕円体（測地計算）
    const GeographicLib::Geodesic geod(
        GeographicLib::Constants::WGS84_a(),
        GeographicLib::Constants::WGS84_f()
    );

    // 出力先：ファイル or /vsimem（stdout用）
    const bool toStdout = (outPath == "-" || outPath.empty());
    const std::string realOut = toStdout ? "/vsimem/out_005.geojson" : outPath;

    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GeoJSON");
    if (!drv) die("GeoJSON driver not found (should be built-in)."); // :contentReference[oaicite:10]{index=10}

    // 既存を消して作り直したい場合（任意）
    // drv->Delete(realOut.c_str()); // 失敗しても無視されることがあります

    GDALDataset* outDs = drv->Create(realOut.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!outDs) die("Failed to create output: " + realOut);

    json summary = {
        {"input", inPath},
        {"output", toStdout ? "(stdout via /vsimem/)" : realOut},
        {"features_total", 0},
        {"features_written", 0},
        {"features_with_metrics", 0},
        {"total_area_m2", 0.0},
        {"total_perimeter_outer_m", 0.0},
        {"total_perimeter_all_rings_m", 0.0}
    };

    const int layerCount = inDs->GetLayerCount();
    for (int li = 0; li < layerCount; li++) {
        OGRLayer* inLy = inDs->GetLayer(li);
        if (!inLy) continue;

        const OGRSpatialReference* srs = inLy->GetSpatialRef();
        if (srs && !srs->IsGeographic()) {
            std::cerr << "[warn] Layer SRS is not geographic. This tool assumes lon/lat degrees.\n";
        }

        OGRFeatureDefn* inDefn = inLy->GetLayerDefn();
        const auto geomType = inDefn->GetGeomType();

        OGRLayer* outLy = outDs->CreateLayer(inLy->GetName(), srs, geomType, nullptr);
        if (!outLy) die("Failed to CreateLayer()");

        // 既存フィールドをコピー
        for (int fi = 0; fi < inDefn->GetFieldCount(); fi++) {
            OGRFieldDefn* fdef = inDefn->GetFieldDefn(fi);
            if (outLy->CreateField(fdef) != OGRERR_NONE) { // CreateField :contentReference[oaicite:11]{index=11}
                die("Failed to CreateField() while copying schema");
            }
        }

        // 追加フィールド
        {
            OGRFieldDefn fArea("geod_area_m2", OFTReal);
            OGRFieldDefn fPerO("geod_perim_outer_m", OFTReal);
            OGRFieldDefn fPerA("geod_perim_all_m", OFTReal);

            if (outLy->CreateField(&fArea) != OGRERR_NONE) die("Failed to add geod_area_m2");
            if (outLy->CreateField(&fPerO) != OGRERR_NONE) die("Failed to add geod_perim_outer_m");
            if (outLy->CreateField(&fPerA) != OGRERR_NONE) die("Failed to add geod_perim_all_m");
        }

        inLy->ResetReading();
        while (OGRFeature* f = inLy->GetNextFeature()) {
            summary["features_total"] = summary["features_total"].get<int>() + 1;

            OGRFeature* outF = OGRFeature::CreateFeature(outLy->GetLayerDefn());
            outF->SetFrom(f, TRUE); // 既存属性＋ジオメトリをコピー

            const OGRGeometry* g = f->GetGeometryRef();
            Metrics m = geometry_metrics(geod, g);
            if (m.ok) {
                outF->SetField("geod_area_m2", m.area_m2);
                outF->SetField("geod_perim_outer_m", m.perimeter_outer_m);
                outF->SetField("geod_perim_all_m", m.perimeter_all_rings_m);

                summary["features_with_metrics"] = summary["features_with_metrics"].get<int>() + 1;
                summary["total_area_m2"] = summary["total_area_m2"].get<double>() + m.area_m2;
                summary["total_perimeter_outer_m"] =
                    summary["total_perimeter_outer_m"].get<double>() + m.perimeter_outer_m;
                summary["total_perimeter_all_rings_m"] =
                    summary["total_perimeter_all_rings_m"].get<double>() + m.perimeter_all_rings_m;
            }

            if (outLy->CreateFeature(outF) != OGRERR_NONE) { // CreateFeature :contentReference[oaicite:12]{index=12}
                die("Failed to CreateFeature()");
            }
            summary["features_written"] = summary["features_written"].get<int>() + 1;

            OGRFeature::DestroyFeature(outF);
            OGRFeature::DestroyFeature(f);
        }
    }

    GDALClose(inDs);
    GDALClose(outDs);

    // stdoutモード：/vsimem から内容を取り出して表示（C++ならでは）:contentReference[oaicite:13]{index=13}
    if (toStdout) {
        vsi_l_offset nSize = 0;
        GByte* p = VSIGetMemFileBuffer(realOut.c_str(), &nSize, TRUE);
        if (!p) die("VSIGetMemFileBuffer() failed.");
        std::cout.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(nSize));
        CPLFree(p);
    }

    // サマリはstderrへ（stdoutがGeoJSONでも混ざらないように）
    std::cerr << summary.dump(2) << "\n";

    return 0;
}
