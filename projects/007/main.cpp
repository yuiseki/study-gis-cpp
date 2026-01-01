#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <proj.h>

// GeographicLib: 基準となる測地距離（楕円体上の最短経路）を計算する
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geodesic.hpp>

// Boost.Geometry: 距離計算を strategy を差し替えて比較する
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/srs/spheroid.hpp>
#include <boost/geometry/strategies/geographic/distance_andoyer.hpp>
#include <boost/geometry/strategies/geographic/distance_thomas.hpp>
#include <boost/geometry/strategies/geographic/distance_vincenty.hpp>
#include <boost/geometry/strategies/geographic/distance_karney.hpp>
#include <boost/geometry/strategies/spherical/distance_haversine.hpp>

#include <nlohmann/json.hpp>

using nlohmann::json;

namespace bg = boost::geometry;

// -------------------------------
// 007: Distance strategies compare
//
// 目的:
//   同じ2点 (lon/lat, deg) に対し、複数の距離モデル/strategy を一括で計算し JSON で出力する。
//   - GeographicLib (WGS84, 楕円体上の測地線距離) を「基準（真値寄り）」として扱う
//   - Boost.Geometry: geographic strategies (andoyer / thomas / vincenty / karney) で比較
//   - Boost.Geometry: spherical haversine (球) も比較
//   - PROJ: EPSG:4326 -> EPSG:3857 (WebMercator) に投影し、ユークリッド距離（平面）を出す
//
// キーワード整理（コメントに解説を残します）:
//
// [GeographicLib の s12]
//   Geodesic::Inverse() の出力 s12 は「点1→点2 を結ぶ測地線（最短経路）に沿った距離（m）」。
//   ここでは基準距離 geodesic_s12_m として使います。
//   ※ 大文字 S12 は “面積” 系の量で別物。小文字 s12 が距離です。
//
// [Boost.Geometry の strategy]
//   distance(p1, p2, strategy) で計算法を差し替え可能。
//   同じ geographic (lon/lat) でも strategy によって誤差特性や安定性が変わる。
//   - Andoyer / Thomas / Vincenty / Karney: いずれも回転楕円体（spheroid）上の距離
//   - Haversine: 完全な球を仮定して距離（半径を指定）
//
// [proj_webmerc_euclid_m]
//   EPSG:4326 を EPSG:3857（WebMercator / Pseudo-Mercator）に投影し、平面上の距離を hypot で計算。
//   WebMercator は可視化向けで、距離・面積の厳密さは保証しません（特に緯度が高いほど歪む）。
//
// [JSON の数値精度は足りる？]
//   nlohmann_json の浮動小数はデフォルトで double。
//   今回の用途（距離m・誤差比較）なら十分です。極端に厳密な桁が必要なら将来文字列化等を検討。
// -------------------------------

struct Args {
    double lon1 = 139.767125;  // 東京駅あたり
    double lat1 = 35.681236;
    double lon2 = 135.502253;  // 大阪駅あたり
    double lat2 = 34.702485;
};

static void die(const std::string& msg, int code = 2) {
    std::cerr << msg << "\n";
    std::exit(code);
}

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string k = argv[i];

        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) die(std::string("Missing value for ") + name);
            return std::string(argv[++i]);
        };

        if (k == "--lon1") a.lon1 = std::stod(need("--lon1"));
        else if (k == "--lat1") a.lat1 = std::stod(need("--lat1"));
        else if (k == "--lon2") a.lon2 = std::stod(need("--lon2"));
        else if (k == "--lat2") a.lat2 = std::stod(need("--lat2"));
        else if (k == "--help") {
            std::cout
                << "Usage:\n"
                << "  007 --lon1 <deg> --lat1 <deg> --lon2 <deg> --lat2 <deg>\n"
                << "\n"
                << "Output:\n"
                << "  JSON to stdout.\n";
            std::exit(0);
        } else {
            die("Unknown option: " + k);
        }
    }
    return a;
}

struct MethodResult {
    std::string id;      // machine-friendly key
    std::string family;  // e.g. "GeographicLib", "Boost.Geometry", "PROJ"
    std::string note;    // short explanation
    bool ok = true;
    double distance_m = std::numeric_limits<double>::quiet_NaN();
    std::string error;

    // reference comparison (vs GeographicLib s12)
    double abs_error_m = std::numeric_limits<double>::quiet_NaN();
    double rel_error = std::numeric_limits<double>::quiet_NaN(); // abs_error / ref
};

static json to_json(const MethodResult& r) {
    json j;
    j["id"] = r.id;
    j["family"] = r.family;
    j["note"] = r.note;
    j["ok"] = r.ok;
    if (r.ok) {
        j["distance_m"] = r.distance_m;
        j["abs_error_m"] = r.abs_error_m;
        j["rel_error"] = r.rel_error;
    } else {
        j["error"] = r.error;
    }
    return j;
}

static void attach_errors(MethodResult& r, double ref_m) {
    if (!r.ok || !std::isfinite(r.distance_m) || !std::isfinite(ref_m) || ref_m == 0.0) return;
    r.abs_error_m = std::abs(r.distance_m - ref_m);
    r.rel_error = r.abs_error_m / ref_m;
}

static PJ* make_norm_transform(PJ_CONTEXT* ctx, const char* src, const char* dst) {
    PJ* P = proj_create_crs_to_crs(ctx, src, dst, nullptr);
    if (!P) return nullptr;

    // 軸順を「伝統的なGIS順（lon,lat / easting,northing）」に正規化して扱いやすくする
    PJ* N = proj_normalize_for_visualization(ctx, P);
    proj_destroy(P);
    return N;
}

int main(int argc, char** argv) {
    const Args a = parse_args(argc, argv);

    // ---- GeographicLib（基準）: WGS84楕円体上の測地線距離 s12 ----
    MethodResult ref;
    ref.id = "geographiclib_s12";
    ref.family = "GeographicLib";
    ref.note =
        "WGS84 ellipsoid geodesic distance (Inverse). s12 = distance along the shortest geodesic [m].";

    double s12_m = std::numeric_limits<double>::quiet_NaN();
    double azi1_deg = std::numeric_limits<double>::quiet_NaN();
    double azi2_deg = std::numeric_limits<double>::quiet_NaN();

    try {
        const auto& geod = GeographicLib::Geodesic::WGS84();
        // Inverse の引数は (lat, lon) の順、結果 s12 は meters、azi は degrees
        geod.Inverse(a.lat1, a.lon1, a.lat2, a.lon2, s12_m, azi1_deg, azi2_deg);
        ref.distance_m = s12_m;
    } catch (const std::exception& e) {
        ref.ok = false;
        ref.error = e.what();
    }

    // ---- Boost.Geometry の準備 ----
    // Boost.Geometry の点は (x, y) を (lon, lat) として使うのが通例（geographic<degree>）
    using GeoPoint = bg::model::point<double, 2, bg::cs::geographic<bg::degree>>;
    GeoPoint p1(a.lon1, a.lat1);
    GeoPoint p2(a.lon2, a.lat2);

    // strategy に与える spheroid (a, f) は WGS84 を使う
    const double WGS84_a = GeographicLib::Constants::WGS84_a();
    const double WGS84_f = GeographicLib::Constants::WGS84_f();
    const double WGS84_b = WGS84_a * (1.0 - WGS84_f); // 短半径 b
    bg::srs::spheroid<double> spheroid(WGS84_a, WGS84_b);

    std::vector<MethodResult> methods;

    // 1) Boost.Geometry geographic strategies (ellipsoid)
    {
        MethodResult r;
        r.id = "boost_geographic_andoyer";
        r.family = "Boost.Geometry";
        r.note = "Ellipsoidal distance on WGS84 spheroid using Andoyer approximation.";
        try {
            bg::strategy::distance::andoyer<bg::srs::spheroid<double>> strat(spheroid);
            r.distance_m = bg::distance(p1, p2, strat);
        } catch (const std::exception& e) {
            r.ok = false; r.error = e.what();
        }
        methods.push_back(r);
    }
    {
        MethodResult r;
        r.id = "boost_geographic_thomas";
        r.family = "Boost.Geometry";
        r.note = "Ellipsoidal distance on WGS84 spheroid using Thomas formula (series-based).";
        try {
            bg::strategy::distance::thomas<bg::srs::spheroid<double>> strat(spheroid);
            r.distance_m = bg::distance(p1, p2, strat);
        } catch (const std::exception& e) {
            r.ok = false; r.error = e.what();
        }
        methods.push_back(r);
    }
    {
        MethodResult r;
        r.id = "boost_geographic_vincenty";
        r.family = "Boost.Geometry";
        r.note = "Ellipsoidal distance on WGS84 spheroid using Vincenty inverse method (iterative).";
        try {
            bg::strategy::distance::vincenty<bg::srs::spheroid<double>> strat(spheroid);
            r.distance_m = bg::distance(p1, p2, strat);
        } catch (const std::exception& e) {
            r.ok = false; r.error = e.what();
        }
        methods.push_back(r);
    }
    {
        MethodResult r;
        r.id = "boost_geographic_karney";
        r.family = "Boost.Geometry";
        r.note = "Ellipsoidal distance on WGS84 spheroid using Karney (2011) inverse geodesic algorithm.";
        try {
            bg::strategy::distance::karney<bg::srs::spheroid<double>> strat(spheroid);
            r.distance_m = bg::distance(p1, p2, strat);
        } catch (const std::exception& e) {
            r.ok = false; r.error = e.what();
        }
        methods.push_back(r);
    }

    // 2) Boost.Geometry spherical haversine (perfect sphere)
    //    球半径は WGS84_a を入れています（簡単な比較用。球モデルは本質的に近似です）
    {
        MethodResult r;
        r.id = "boost_spherical_haversine";
        r.family = "Boost.Geometry";
        r.note = "Spherical distance using haversine on a perfect sphere (radius = WGS84_a).";
        try {
            using SphPoint = bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>>;
            SphPoint s1(a.lon1, a.lat1);
            SphPoint s2(a.lon2, a.lat2);

            bg::strategy::distance::haversine<double> strat(WGS84_a);
            r.distance_m = bg::distance(s1, s2, strat);
        } catch (const std::exception& e) {
            r.ok = false; r.error = e.what();
        }
        methods.push_back(r);
    }

    // 3) PROJ: EPSG:4326 -> EPSG:3857 に投影して平面ユークリッド距離
    {
        MethodResult r;
        r.id = "proj_webmerc_euclid";
        r.family = "PROJ";
        r.note =
            "Project to EPSG:3857 (WebMercator) then compute Euclidean distance in meters (NOT geodesic).";
        PJ_CONTEXT* ctx = proj_context_create();
        PJ* T = nullptr;
        try {
            if (!ctx) throw std::runtime_error("proj_context_create failed");
            T = make_norm_transform(ctx, "EPSG:4326", "EPSG:3857");
            if (!T) throw std::runtime_error("proj_create_crs_to_crs / proj_normalize_for_visualization failed");

            PJ_COORD c1 = proj_coord(a.lon1, a.lat1, 0, 0);
            PJ_COORD c2 = proj_coord(a.lon2, a.lat2, 0, 0);
            PJ_COORD m1 = proj_trans(T, PJ_FWD, c1);
            PJ_COORD m2 = proj_trans(T, PJ_FWD, c2);

            const double dx = m2.xy.x - m1.xy.x;
            const double dy = m2.xy.y - m1.xy.y;
            r.distance_m = std::hypot(dx, dy);
        } catch (const std::exception& e) {
            r.ok = false; r.error = e.what();
        }

        if (T) proj_destroy(T);
        if (ctx) proj_context_destroy(ctx);
        methods.push_back(r);
    }

    // ---- 誤差計算（基準 = GeographicLib s12）----
    for (auto& m : methods) attach_errors(m, ref.distance_m);

    // ---- 出力 JSON ----
    json out;
    out["input"] = {
        {"lon1", a.lon1}, {"lat1", a.lat1},
        {"lon2", a.lon2}, {"lat2", a.lat2},
        {"units", {{"angles", "degrees"}, {"distance", "meters"}}},
        {"assumption", "Input lon/lat are in EPSG:4326 (WGS84)."}
    };

    out["reference"] = {
        {"id", ref.id},
        {"family", ref.family},
        {"note", ref.note},
        {"ok", ref.ok}
    };

    if (ref.ok) {
        out["reference"]["distance_m"] = ref.distance_m;
        out["reference"]["azi1_deg"] = azi1_deg;
        out["reference"]["azi2_deg"] = azi2_deg;
    } else {
        out["reference"]["error"] = ref.error;
    }

    json arr = json::array();
    for (const auto& m : methods) arr.push_back(to_json(m));
    out["methods"] = arr;

    // 便利な集計
    out["summary"] = {
        {"ref_distance_m", ref.distance_m},
        {"methods_count", methods.size()},
        {"note",
         "Compare multiple distance models/strategies. Use GeographicLib s12 as baseline; PROJ WebMercator+Euclid is included as a 'distortion grows' example."}
    };

    std::cout << out.dump(2) << "\n";
    return ref.ok ? 0 : 1;
}
