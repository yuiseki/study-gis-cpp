#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/PolygonArea.hpp>

using nlohmann::json;

struct Args {
    std::string path;
    bool pretty = true;
    int feature_index = -1; // -1 = all
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

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string k = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) die(std::string("Missing value for ") + name);
            return std::string(argv[++i]);
        };

        if (k == "--help") {
            std::cout
                << "Usage: 003 <geojson_path> [--feature <i>] [--compact]\n"
                << "  --feature <i>   : compute only feature i (0-based) in a FeatureCollection\n"
                << "  --compact       : compact JSON output\n";
            std::exit(0);
        } else if (k == "--feature") {
            a.feature_index = std::stoi(need("--feature"));
        } else if (k == "--compact") {
            a.pretty = false;
        } else if (!k.empty() && k[0] == '-') {
            die("Unknown option: " + k);
        } else {
            a.path = k;
        }
    }
    if (a.path.empty()) die("GeoJSON file path is required. Try --help");
    return a;
}

static json load_json_file(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) die("Failed to open file: " + path);
    json j;
    try {
        ifs >> j;
    } catch (const std::exception& e) {
        die(std::string("Failed to parse JSON: ") + e.what());
    }
    return j;
}

struct RingMetrics {
    double perimeter_m = 0.0;
    double area_m2_abs = 0.0;
    int vertices_used = 0;
};

static bool nearly_equal(double a, double b, double eps = 1e-12) {
    return std::abs(a - b) <= eps;
}

// ring: [[lon,lat,...], [lon,lat,...], ...]
// returns abs(area) and perimeter (both on ellipsoid)
static RingMetrics compute_ring_geodesic(
    const GeographicLib::Geodesic& geod,
    const json& ring_coords
) {
    if (!ring_coords.is_array() || ring_coords.size() < 4) {
        die("Invalid LinearRing: must be an array with >= 4 positions");
    }

    // Extract positions and drop closing point if it's a duplicate of the first.
    std::vector<std::pair<double,double>> pts;
    pts.reserve(ring_coords.size());

    for (const auto& pos : ring_coords) {
        if (!pos.is_array() || pos.size() < 2) die("Invalid position: expected [lon, lat, ...]");
        double lon = pos.at(0).get<double>();
        double lat = pos.at(1).get<double>();
        lon = wrap_lon(lon);
        pts.push_back({lon, lat});
    }

    if (pts.size() >= 2) {
        const auto& first = pts.front();
        const auto& last  = pts.back();
        if (nearly_equal(first.first, last.first) && nearly_equal(first.second, last.second)) {
            pts.pop_back(); // avoid double edge
        }
    }

    if (pts.size() < 3) die("Invalid LinearRing after de-dup: need >= 3 distinct vertices");

    GeographicLib::PolygonArea poly(geod /*polyline=false by default*/);

    for (const auto& [lon, lat] : pts) {
        // PolygonArea expects AddPoint(lat, lon)
        poly.AddPoint(lat, lon);
    }

    double perimeter = 0.0, area = 0.0;
    poly.Compute(false, true, perimeter, area);

    RingMetrics rm;
    rm.perimeter_m = perimeter;
    rm.area_m2_abs = std::abs(area);
    rm.vertices_used = static_cast<int>(pts.size());
    return rm;
}

// polygon_coords: [ ring0, ring1, ... ]
static json compute_polygon_metrics(
    const GeographicLib::Geodesic& geod,
    const json& polygon_coords
) {
    if (!polygon_coords.is_array() || polygon_coords.empty()) {
        die("Invalid Polygon coordinates: expected non-empty array of rings");
    }

    // Outer ring
    RingMetrics outer = compute_ring_geodesic(geod, polygon_coords.at(0));

    double area = outer.area_m2_abs;
    double per_outer = outer.perimeter_m;
    double per_all = outer.perimeter_m;

    int holes = 0;
    int vertices = outer.vertices_used;

    // Holes
    for (size_t i = 1; i < polygon_coords.size(); i++) {
        RingMetrics hole = compute_ring_geodesic(geod, polygon_coords.at(i));
        area -= hole.area_m2_abs;
        per_all += hole.perimeter_m;
        holes++;
        vertices += hole.vertices_used;
    }

    if (area < 0) area = 0; // guard tiny numeric issues

    return json{
        {"area_m2", area},
        {"perimeter_m_outer", per_outer},
        {"perimeter_m_all_rings", per_all},
        {"holes", holes},
        {"vertices_used_total", vertices}
    };
}

static void accumulate_totals(json& total, const json& poly) {
    total["polygons"] = total["polygons"].get<int>() + 1;
    total["area_m2"] = total["area_m2"].get<double>() + poly["area_m2"].get<double>();
    total["perimeter_m_outer"] =
        total["perimeter_m_outer"].get<double>() + poly["perimeter_m_outer"].get<double>();
    total["perimeter_m_all_rings"] =
        total["perimeter_m_all_rings"].get<double>() + poly["perimeter_m_all_rings"].get<double>();
}

static json compute_geometry_metrics(const GeographicLib::Geodesic& geod, const json& geom) {
    if (!geom.is_object()) die("Geometry must be an object");
    const std::string type = geom.value("type", "");
    if (type.empty()) die("Geometry missing 'type'");

    json out;
    out["type"] = type;
    out["polygons"] = json::array();

    json total = {
        {"polygons", 0},
        {"area_m2", 0.0},
        {"perimeter_m_outer", 0.0},
        {"perimeter_m_all_rings", 0.0}
    };

    if (type == "Polygon") {
        const json& coords = geom.at("coordinates");
        json poly = compute_polygon_metrics(geod, coords);
        out["polygons"].push_back(poly);
        accumulate_totals(total, poly);

    } else if (type == "MultiPolygon") {
        const json& coords = geom.at("coordinates");
        if (!coords.is_array()) die("MultiPolygon coordinates must be an array");
        for (const auto& poly_coords : coords) {
            json poly = compute_polygon_metrics(geod, poly_coords);
            out["polygons"].push_back(poly);
            accumulate_totals(total, poly);
        }
    } else {
        out["warning"] = "Unsupported geometry type for area/perimeter (only Polygon/MultiPolygon).";
    }

    out["total"] = total;
    return out;
}

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    const json root = load_json_file(args.path);

    const GeographicLib::Geodesic geod(
        GeographicLib::Constants::WGS84_a(),
        GeographicLib::Constants::WGS84_f()
    );

    json result;
    result["input"] = args.path;
    result["assumption"] = "GeoJSON coordinates are [lon, lat] in EPSG:4326 (RFC 7946).";
    result["features"] = json::array();

    auto push_feature_result = [&](int idx, const json& geom, const json& props) {
        json item;
        item["feature_index"] = idx;
        item["properties"] = props;
        item["geometry_metrics"] = compute_geometry_metrics(geod, geom);
        result["features"].push_back(item);
    };

    const std::string top_type = root.value("type", "");

    if (top_type == "FeatureCollection") {
        const json& feats = root.at("features");
        if (!feats.is_array()) die("FeatureCollection.features must be an array");

        for (size_t i = 0; i < feats.size(); i++) {
            if (args.feature_index >= 0 && static_cast<int>(i) != args.feature_index) continue;

            const json& f = feats.at(i);
            if (!f.is_object() || f.value("type", "") != "Feature") continue;
            if (!f.contains("geometry") || f.at("geometry").is_null()) continue;

            const json& geom = f.at("geometry");
            const json props = f.value("properties", json::object());
            push_feature_result(static_cast<int>(i), geom, props);
        }

    } else if (top_type == "Feature") {
        if (!root.contains("geometry") || root.at("geometry").is_null()) die("Feature has null geometry");
        const json& geom = root.at("geometry");
        const json props = root.value("properties", json::object());
        push_feature_result(0, geom, props);

    } else if (!top_type.empty()) {
        // geometry object directly
        push_feature_result(0, root, json::object());

    } else {
        die("Unknown top-level JSON (missing 'type')");
    }

    // Aggregate totals across features (Polygon/MultiPolygon only)
    json grand = {
        {"features_count", static_cast<int>(result["features"].size())},
        {"polygons", 0},
        {"area_m2", 0.0},
        {"perimeter_m_outer", 0.0},
        {"perimeter_m_all_rings", 0.0}
    };

    for (const auto& f : result["features"]) {
        const auto& gm = f["geometry_metrics"];
        if (!gm.contains("total")) continue;
        const auto& t = gm["total"];
        grand["polygons"] = grand["polygons"].get<int>() + t["polygons"].get<int>();
        grand["area_m2"] = grand["area_m2"].get<double>() + t["area_m2"].get<double>();
        grand["perimeter_m_outer"] =
            grand["perimeter_m_outer"].get<double>() + t["perimeter_m_outer"].get<double>();
        grand["perimeter_m_all_rings"] =
            grand["perimeter_m_all_rings"].get<double>() + t["perimeter_m_all_rings"].get<double>();
    }

    result["grand_total"] = grand;

    if (args.pretty) std::cout << result.dump(2) << "\n";
    else std::cout << result.dump() << "\n";

    return 0;
}
