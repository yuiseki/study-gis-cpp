#include <cstdlib>
#include <iostream>
#include <string>

#include <proj.h>
#include <geos_c.h>
#include <nlohmann/json.hpp>

using nlohmann::json;

// GEOSMessageHandler_r: void (*)(const char* message, void* userdata)
static void geos_notice(const char* /*msg*/, void* /*userdata*/) {}
static void geos_error(const char* /*msg*/, void* /*userdata*/) {}

struct Args {
    double lon = 139.767125;   // 東京駅あたり
    double lat = 35.681236;
    double radius_m = 1000.0;
    int quad_segs = 8;         // バッファ円弧の細かさ
};

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string k = argv[i];
        auto need = [&](const char* name) {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                std::exit(2);
            }
            return std::string(argv[++i]);
        };

        if (k == "--lon") a.lon = std::stod(need("--lon"));
        else if (k == "--lat") a.lat = std::stod(need("--lat"));
        else if (k == "--radius") a.radius_m = std::stod(need("--radius"));
        else if (k == "--quad") a.quad_segs = std::stoi(need("--quad"));
        else if (k == "--help") {
            std::cout << "Usage: 001 [--lon <deg>] [--lat <deg>] [--radius <m>] [--quad <n>]\n";
            std::exit(0);
        }
    }
    return a;
}

// PROJ: CRS->CRS 変換を作り、可視化向け（lon,lat / easting,northing）に正規化する
static PJ* make_norm_transform(PJ_CONTEXT* ctx, const char* src, const char* dst) {
    PJ* P = proj_create_crs_to_crs(ctx, src, dst, nullptr);
    if (!P) return nullptr;
    PJ* N = proj_normalize_for_visualization(ctx, P);
    proj_destroy(P);
    return N;
}

// GEOSリング（LinearRing）の座標列を取り出して、EPSG:4326 (lon,lat) に戻してGeoJSON座標配列にする
static json ring_to_geojson_coords(GEOSContextHandle_t gctx, const GEOSGeometry* ring, PJ* proj_4326_to_3857_norm) {
    const GEOSCoordSequence* cs = GEOSGeom_getCoordSeq_r(gctx, ring);
    if (!cs) return json::array();

    unsigned int n = 0;
    if (!GEOSCoordSeq_getSize_r(gctx, cs, &n)) return json::array();

    json coords = json::array();
    // json自体に reserve() はないので、内部の array_t 参照を取得して reserve する
    coords.get_ref<json::array_t&>().reserve(n);

    for (unsigned int i = 0; i < n; i++) {
        double x = 0.0, y = 0.0;
        GEOSCoordSeq_getX_r(gctx, cs, i, &x);
        GEOSCoordSeq_getY_r(gctx, cs, i, &y);

        // 逆変換: EPSG:3857 -> EPSG:4326
        PJ_COORD c = proj_coord(x, y, 0, 0);
        PJ_COORD out = proj_trans(proj_4326_to_3857_norm, PJ_INV, c);

        // 正規化済みなので out.xy.x = lon(deg), out.xy.y = lat(deg) として扱う
        coords.push_back(json::array({out.xy.x, out.xy.y}));
    }

    return coords;
}

static json polygon_to_geojson(GEOSContextHandle_t gctx, const GEOSGeometry* poly, PJ* proj_norm) {
    const GEOSGeometry* ext = GEOSGetExteriorRing_r(gctx, poly); // 内部参照（freeしない）
    json rings = json::array();
    rings.push_back(ring_to_geojson_coords(gctx, ext, proj_norm));

    int holes = GEOSGetNumInteriorRings_r(gctx, poly);
    for (int i = 0; i < holes; i++) {
        const GEOSGeometry* in = GEOSGetInteriorRingN_r(gctx, poly, i); // 内部参照
        rings.push_back(ring_to_geojson_coords(gctx, in, proj_norm));
    }

    return json{
        {"type", "Polygon"},
        {"coordinates", rings}
    };
}

int main(int argc, char** argv) {
    const Args a = parse_args(argc, argv);

    // ---- PROJ ----
    PJ_CONTEXT* pj_ctx = proj_context_create();
    if (!pj_ctx) {
        std::cerr << "Failed to create PROJ context\n";
        return 1;
    }

    // EPSG:4326 -> EPSG:3857（正規化済み）
    PJ* T = make_norm_transform(pj_ctx, "EPSG:4326", "EPSG:3857");
    if (!T) {
        std::cerr << "proj_create_crs_to_crs / normalize failed\n";
        proj_context_destroy(pj_ctx);
        return 1;
    }

    // forward: lon,lat(deg) -> x,y(m)
    PJ_COORD c = proj_coord(a.lon, a.lat, 0, 0);
    PJ_COORD m = proj_trans(T, PJ_FWD, c);
    const double x = m.xy.x;
    const double y = m.xy.y;

    // ---- GEOS ----
    GEOSContextHandle_t gctx = GEOS_init_r();

    // userData引数が必要（nullptrでOK）
    GEOSContext_setNoticeMessageHandler_r(gctx, geos_notice, nullptr);
    GEOSContext_setErrorMessageHandler_r(gctx, geos_error, nullptr);

    GEOSCoordSequence* cs = GEOSCoordSeq_create_r(gctx, 1, 2);
    GEOSCoordSeq_setX_r(gctx, cs, 0, x);
    GEOSCoordSeq_setY_r(gctx, cs, 0, y);
    GEOSGeometry* pt = GEOSGeom_createPoint_r(gctx, cs);

    GEOSGeometry* buf = GEOSBuffer_r(gctx, pt, a.radius_m, a.quad_segs);
    GEOSGeom_destroy_r(gctx, pt);

    if (!buf) {
        std::cerr << "GEOSBuffer_r failed\n";
        GEOS_finish_r(gctx);
        proj_destroy(T);
        proj_context_destroy(pj_ctx);
        return 1;
    }

    // ---- GeoJSON化（EPSG:4326へ戻す）----
    json geom;
    const int type = GEOSGeomTypeId_r(gctx, buf);

    if (type == GEOS_POLYGON) {
        geom = polygon_to_geojson(gctx, buf, T);
    } else if (type == GEOS_MULTIPOLYGON) {
        json polys = json::array();
        const int n = GEOSGetNumGeometries_r(gctx, buf);
        polys.get_ref<json::array_t&>().reserve(static_cast<size_t>(n));

        for (int i = 0; i < n; i++) {
            const GEOSGeometry* p = GEOSGetGeometryN_r(gctx, buf, i); // 内部参照
            polys.push_back(polygon_to_geojson(gctx, p, T)["coordinates"]);
        }

        geom = json{{"type", "MultiPolygon"}, {"coordinates", polys}};
    } else {
        std::cerr << "Unexpected GEOS result type id: " << type << "\n";
        GEOSGeom_destroy_r(gctx, buf);
        GEOS_finish_r(gctx);
        proj_destroy(T);
        proj_context_destroy(pj_ctx);
        return 1;
    }

    json fc = {
        {"type", "FeatureCollection"},
        {"features", json::array({
            {
                {"type", "Feature"},
                {"properties", {
                    {"src", "EPSG:4326"},
                    {"radius_m", a.radius_m},
                    {"quad_segs", a.quad_segs}
                }},
                {"geometry", geom}
            }
        })}
    };

    std::cout << fc.dump(2) << "\n";

    GEOSGeom_destroy_r(gctx, buf);
    GEOS_finish_r(gctx);

    proj_destroy(T);
    proj_context_destroy(pj_ctx);

    return 0;
}
