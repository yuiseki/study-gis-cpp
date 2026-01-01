#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <geos_c.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>

#include <nlohmann/json.hpp>

using nlohmann::json;
namespace bg = boost::geometry;

// ------------------------------------------------------------
// 008: GEOS vs Boost.Geometry polygon boolean (PLANAR)
//
// 目的（比較軸を明確に）
//   同じ入力（平面ポリゴン）に対して、GEOS と Boost.Geometry で
//   union / intersection / difference / symdiff を実行し、結果を比較する。
//
// 前提（重要）
//   - ここでの比較は「平面（2Dユークリッド）」前提です。
//     EPSG:4326（lon/lat degrees）のまま “面積・距離” を評価する用途ではありません。
//     006でEPSG:3857等の平面CRSへ投影した上で使うのが王道です。
//   - GEOSもBoost.Geometryも、原則として「座標はただの(x,y)」として扱います。
//     測地（楕円体）上の面積・距離とは別物です（007/005の世界）。
//
// ライブラリ実装の雰囲気（コメント内に残す）
//
// [GEOS]
//   - JTS系（OverlayNG）の流れ：線分同士の交差（noding）を作り、planar graph から結果面を組み立てる。
//   - 精度モデルに応じて noder を切り替える設計があり、固定精度＋snap-rounding で robust を狙える。
//     浮動小数（フル精度）では例外（TopologyException）になり得る、という思想。
//   - 008では C API を使用し、WKT入出力と boolean を素直に叩く。
//
// [Boost.Geometry]
//   - overlay は「adapted Weiler–Atherton」系の traversal（境界を辿って結果リングを構成する）と説明される。
//   - テンプレート・ヘッダオンリーで高速に回せる一方、退化ケース（ほぼ接する/極細スリバーなど）で差が出やすい。
//   - 008では strategy の議論をせず、まず “素の結果” を観察する。
//     （将来、丸め・スナップ等の前処理を入れて差を見るのも面白い）
//
// 入力
//   --a <WKT or file>   : 1つ目のポリゴン（POLYGON / MULTIPOLYGON）
//   --b <WKT or file>   : 2つ目
//   --op <intersection|union|difference|symdiff>
//   [--emit-wkt]        : 結果WKTをJSONへ含める（長いので任意）
//
// 出力
//   JSON to stdout
// ------------------------------------------------------------

struct Args {
    std::string a;
    std::string b;
    std::string op = "intersection";
    bool emit_wkt = false;
};

static void die(const std::string& msg, int code = 2) {
    std::cerr << msg << "\n";
    std::exit(code);
}

static Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; i++) {
        std::string k = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) die(std::string("Missing value for ") + name);
            return std::string(argv[++i]);
        };

        if (k == "--a") args.a = need("--a");
        else if (k == "--b") args.b = need("--b");
        else if (k == "--op") args.op = need("--op");
        else if (k == "--emit-wkt") args.emit_wkt = true;
        else if (k == "--help") {
            std::cout
                << "Usage: 008 --a <WKT|file> --b <WKT|file> --op <intersection|union|difference|symdiff> [--emit-wkt]\n";
            std::exit(0);
        } else {
            die("Unknown option: " + k);
        }
    }
    if (args.a.empty() || args.b.empty()) {
        die("Both --a and --b are required. Use --help for usage.");
    }
    return args;
}

static std::optional<std::string> read_file_if_exists(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) return std::nullopt;
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

static std::string load_wkt_or_file(const std::string& s) {
    if (auto txt = read_file_if_exists(s)) return *txt;
    return s; // treat as inline WKT
}

// ---------- GEOS helpers ----------
static void geos_notice(const char* msg, void*) { (void)msg; }
static void geos_error(const char* msg, void*) { (void)msg; }

static double geos_area(GEOSContextHandle_t gctx, const GEOSGeometry* g) {
    double a = 0.0;
    if (!g) return std::nan("");
    if (!GEOSArea_r(gctx, g, &a)) return std::nan("");
    return a;
}

static double geos_boundary_length(GEOSContextHandle_t gctx, const GEOSGeometry* g) {
    if (!g) return std::nan("");
    GEOSGeometry* b = GEOSBoundary_r(gctx, g);
    if (!b) return std::nan("");
    double len = 0.0;
    const int ok = GEOSLength_r(gctx, b, &len);
    GEOSGeom_destroy_r(gctx, b);
    return ok ? len : std::nan("");
}

static uint64_t geos_count_coords_ring(GEOSContextHandle_t gctx, const GEOSGeometry* ring) {
    const GEOSCoordSequence* cs = GEOSGeom_getCoordSeq_r(gctx, ring);
    if (!cs) return 0;
    unsigned int n = 0;
    if (!GEOSCoordSeq_getSize_r(gctx, cs, &n)) return 0;
    return static_cast<uint64_t>(n);
}

static uint64_t geos_count_coords_polygon(GEOSContextHandle_t gctx, const GEOSGeometry* poly) {
    if (!poly) return 0;
    uint64_t total = 0;

    const GEOSGeometry* ext = GEOSGetExteriorRing_r(gctx, poly);
    total += geos_count_coords_ring(gctx, ext);

    const int holes = GEOSGetNumInteriorRings_r(gctx, poly);
    for (int i = 0; i < holes; i++) {
        const GEOSGeometry* in = GEOSGetInteriorRingN_r(gctx, poly, i);
        total += geos_count_coords_ring(gctx, in);
    }
    return total;
}

static uint64_t geos_vertex_count(GEOSContextHandle_t gctx, const GEOSGeometry* g) {
    if (!g) return 0;
    const int t = GEOSGeomTypeId_r(gctx, g);

    if (t == GEOS_POLYGON) {
        return geos_count_coords_polygon(gctx, g);
    } else if (t == GEOS_MULTIPOLYGON || t == GEOS_GEOMETRYCOLLECTION) {
        const int n = GEOSGetNumGeometries_r(gctx, g);
        uint64_t total = 0;
        for (int i = 0; i < n; i++) {
            const GEOSGeometry* gi = GEOSGetGeometryN_r(gctx, g, i);
            total += geos_vertex_count(gctx, gi);
        }
        return total;
    }
    // 今回は Polygon/MultiPolygon 前提
    return 0;
}

static std::string geos_to_wkt(GEOSContextHandle_t gctx, const GEOSGeometry* g) {
    if (!g) return "";
    GEOSWKTWriter* w = GEOSWKTWriter_create_r(gctx);
    if (!w) return "";
    GEOSWKTWriter_setRoundingPrecision_r(gctx, w, 10);
    char* s = GEOSWKTWriter_write_r(gctx, w, g);
    GEOSWKTWriter_destroy_r(gctx, w);
    std::string out = s ? std::string(s) : "";
    GEOSFree_r(gctx, s);
    return out;
}

// ---------- Boost.Geometry helpers ----------
using Pt = bg::model::d2::point_xy<double>;
using Poly = bg::model::polygon<Pt>;
using MP = bg::model::multi_polygon<Poly>;

static MP boost_read_multipolygon(const std::string& wkt) {
    // POLYGON / MULTIPOLYGON の両方に対応したいので、まず POLYGON を試してから MULTI を試す
    MP mp;

    try {
        Poly p;
        bg::read_wkt(wkt, p);
        bg::correct(p);
        mp.push_back(p);
        return mp;
    } catch (...) {
        // ignore
    }

    // MULTIPOLYGON を読む
    bg::read_wkt(wkt, mp);
    for (auto& p : mp) bg::correct(p);
    return mp;
}

static uint64_t boost_vertex_count(const MP& mp) {
    uint64_t total = 0;
    for (const auto& p : mp) {
        total += static_cast<uint64_t>(p.outer().size());
        for (const auto& in : p.inners()) total += static_cast<uint64_t>(in.size());
    }
    return total;
}

static std::string boost_to_wkt(const MP& mp) {
    std::ostringstream oss;
    oss << bg::wkt(mp);
    return oss.str();
}

// ---------- Result structs ----------
struct LibResult {
    bool ok = true;
    std::string error;

    bool valid = false;
    std::string valid_reason;

    double area = std::nan("");
    double perimeter_all_rings = std::nan("");

    uint64_t vertices = 0;
    double time_ms = std::nan("");

    std::string wkt; // optional
};

static json to_json(const LibResult& r) {
    json j;
    j["ok"] = r.ok;
    if (!r.ok) {
        j["error"] = r.error;
        return j;
    }
    j["valid"] = r.valid;
    if (!r.valid_reason.empty()) j["valid_reason"] = r.valid_reason;
    j["area"] = r.area;
    j["perimeter_all_rings"] = r.perimeter_all_rings;
    j["vertices"] = r.vertices;
    j["time_ms"] = r.time_ms;
    if (!r.wkt.empty()) j["wkt"] = r.wkt;
    return j;
}

static bool is_supported_op(const std::string& op) {
    return op == "intersection" || op == "union" || op == "difference" || op == "symdiff";
}

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    if (!is_supported_op(args.op)) {
        die("Unsupported --op. Use intersection|union|difference|symdiff");
    }

    const std::string wktA = load_wkt_or_file(args.a);
    const std::string wktB = load_wkt_or_file(args.b);

    json out;
    out["assumption"] = "All coordinates are planar (cartesian). Reproject lon/lat to a planar CRS before using.";
    out["op"] = args.op;

    // ---------------- GEOS path ----------------
    LibResult geosRes;
    {
        GEOSContextHandle_t gctx = GEOS_init_r();
        GEOSContext_setNoticeMessageHandler_r(gctx, geos_notice, nullptr);
        GEOSContext_setErrorMessageHandler_r(gctx, geos_error, nullptr);

        GEOSWKTReader* rdr = GEOSWKTReader_create_r(gctx);
        if (!rdr) {
            geosRes.ok = false;
            geosRes.error = "GEOSWKTReader_create_r failed";
        } else {
            GEOSGeometry* ga = GEOSWKTReader_read_r(gctx, rdr, wktA.c_str());
            GEOSGeometry* gb = GEOSWKTReader_read_r(gctx, rdr, wktB.c_str());

            if (!ga || !gb) {
                geosRes.ok = false;
                geosRes.error = "Failed to read WKT (GEOS). Ensure POLYGON/MULTIPOLYGON.";
            } else {
                // validity (input) - 参考として見る
                // （結果の validity は別途）
                // GEOSisValid_r は OGC validity 判定
                // ※ 退化ケースでは GEOS 側がより堅牢に扱えることが多いが、浮動小数精度では例外もあり得る（TopologyException 等）
                //    それが GEOS の「精度モデルとrobustness」を学べるポイント。
                const auto t0 = std::chrono::steady_clock::now();

                GEOSGeometry* gr = nullptr;
                if (args.op == "intersection") gr = GEOSIntersection_r(gctx, ga, gb);
                else if (args.op == "union") gr = GEOSUnion_r(gctx, ga, gb);
                else if (args.op == "difference") gr = GEOSDifference_r(gctx, ga, gb);
                else if (args.op == "symdiff") gr = GEOSSymDifference_r(gctx, ga, gb);

                const auto t1 = std::chrono::steady_clock::now();
                geosRes.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

                if (!gr) {
                    geosRes.ok = false;
                    geosRes.error = "GEOS boolean operation returned NULL (possible robustness/TopologyException case).";
                } else {
                    geosRes.valid = (GEOSisValid_r(gctx, gr) == 1);
                    geosRes.area = geos_area(gctx, gr);
                    geosRes.perimeter_all_rings = geos_boundary_length(gctx, gr);
                    geosRes.vertices = geos_vertex_count(gctx, gr);
                    if (args.emit_wkt) geosRes.wkt = geos_to_wkt(gctx, gr);

                    GEOSGeom_destroy_r(gctx, gr);
                }

                GEOSGeom_destroy_r(gctx, ga);
                GEOSGeom_destroy_r(gctx, gb);
            }

            GEOSWKTReader_destroy_r(gctx, rdr);
        }

        GEOS_finish_r(gctx);
    }

    // ---------------- Boost.Geometry path ----------------
    LibResult boostRes;
    {
        try {
            MP A = boost_read_multipolygon(wktA);
            MP B = boost_read_multipolygon(wktB);

            std::string reason;
            boostRes.valid = bg::is_valid(A, reason) && bg::is_valid(B, reason);
            if (!boostRes.valid) boostRes.valid_reason = reason;

            MP R;

            const auto t0 = std::chrono::steady_clock::now();
            if (args.op == "intersection") {
                bg::intersection(A, B, R); // set-theoretic intersection
            } else if (args.op == "union") {
                bg::union_(A, B, R);       // set-theoretic union
            } else if (args.op == "difference") {
                bg::difference(A, B, R);   // A \ B
            } else if (args.op == "symdiff") {
                bg::sym_difference(A, B, R);
            }
            const auto t1 = std::chrono::steady_clock::now();
            boostRes.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            // result validity
            std::string rreason;
            boostRes.valid = bg::is_valid(R, rreason);
            if (!boostRes.valid) boostRes.valid_reason = rreason;

            boostRes.area = bg::area(R);
            // perimeter: Boostの perimeter は外周＋内周の扱いがバージョンやジオメトリ型で差が出る可能性があるため、
            // ここでは「比較の一貫性」を優先して GEOS 側と同様に “境界長（全リング）相当” の値として perimeter() を採用する。
            // 必要になったら、outer/inner別に明示的に合算する実装へ発展させる。
            boostRes.perimeter_all_rings = bg::perimeter(R);

            boostRes.vertices = boost_vertex_count(R);
            if (args.emit_wkt) boostRes.wkt = boost_to_wkt(R);

        } catch (const std::exception& e) {
            boostRes.ok = false;
            boostRes.error = std::string("Boost.Geometry error: ") + e.what();
        } catch (...) {
            boostRes.ok = false;
            boostRes.error = "Boost.Geometry error: unknown exception";
        }
    }

    // ---------------- Compare ----------------
    json diff;
    if (geosRes.ok && boostRes.ok) {
        diff["abs_area_diff"] = std::abs(geosRes.area - boostRes.area);
        diff["abs_perimeter_diff"] = std::abs(geosRes.perimeter_all_rings - boostRes.perimeter_all_rings);
        diff["vertices_diff"] = (geosRes.vertices > boostRes.vertices)
            ? (geosRes.vertices - boostRes.vertices)
            : (boostRes.vertices - geosRes.vertices);
        diff["both_valid"] = (geosRes.valid && boostRes.valid);
    } else {
        diff["note"] = "diff is partial because one side failed.";
    }

    out["geos"] = to_json(geosRes);
    out["boost_geometry"] = to_json(boostRes);
    out["diff"] = diff;

    std::cout << out.dump(2) << "\n";
    return (geosRes.ok && boostRes.ok) ? 0 : 1;
}
