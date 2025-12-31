#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/PolygonArea.hpp>

#include <nlohmann/json.hpp>

using nlohmann::json;

struct Args {
    double lon = 139.767125;   // 東京駅あたり
    double lat = 35.681236;
    double radius_m = 1000.0;
    int steps = 360;           // 方位角分割（多いほど滑らか）
};

// [-180, 180) に正規化（GeoJSONビューアの癖対策として最低限）
static double wrap_lon(double lon_deg) {
    lon_deg = std::fmod(lon_deg + 180.0, 360.0);
    if (lon_deg < 0) lon_deg += 360.0;
    return lon_deg - 180.0;
}

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
        else if (k == "--steps") a.steps = std::stoi(need("--steps"));
        else if (k == "--help") {
            std::cout << "Usage: 002 [--lon <deg>] [--lat <deg>] [--radius <m>] [--steps <n>]\n";
            std::exit(0);
        }
    }

    if (a.steps < 8) a.steps = 8;
    return a;
}

int main(int argc, char** argv) {
    const Args a = parse_args(argc, argv);

    // WGS84 楕円体で測地計算
    const GeographicLib::Geodesic geod(
        GeographicLib::Constants::WGS84_a(),
        GeographicLib::Constants::WGS84_f()
    );

    // 測地多角形の周長・面積を計算（polyline=false）
    GeographicLib::PolygonArea poly(geod);

    json ring = json::array();
    ring.get_ref<json::array_t&>().reserve(static_cast<size_t>(a.steps) + 1);

    double first_lon = 0.0, first_lat = 0.0;

    for (int i = 0; i < a.steps; i++) {
        const double azi = (360.0 * i) / a.steps; // deg
        double lat2 = 0.0, lon2 = 0.0;

        // direct: (lat1, lon1, azi1, s12) -> (lat2, lon2)
        geod.Direct(a.lat, a.lon, azi, a.radius_m, lat2, lon2);
        lon2 = wrap_lon(lon2);

        if (i == 0) {
            first_lon = lon2;
            first_lat = lat2;
        }

        // GeoJSONは [lon, lat]
        ring.push_back(json::array({lon2, lat2}));

        // PolygonAreaは (lat, lon)
        poly.AddPoint(lat2, lon2);
    }

    // GeoJSONのPolygon ringは閉じる（最初の点を末尾にもう一度）
    ring.push_back(json::array({first_lon, first_lat}));

    // reverse=false, sign=true（向きが逆でも「残りの地球」扱いにせず符号付きで返す）
    double perimeter_m = 0.0, area_m2_signed = 0.0;
    (void)poly.Compute(false, true, perimeter_m, area_m2_signed);
    const double area_m2_abs = std::abs(area_m2_signed);

    json geom = {
        {"type", "Polygon"},
        {"coordinates", json::array({ring})}
    };

    json fc = {
        {"type", "FeatureCollection"},
        {"features", json::array({
            {
                {"type", "Feature"},
                {"properties", {
                    {"method", "GeographicLib geodesic circle (WGS84)"},
                    {"center_lon", a.lon},
                    {"center_lat", a.lat},
                    {"radius_m", a.radius_m},
                    {"steps", a.steps},
                    {"perimeter_m", perimeter_m},
                    {"area_m2_abs", area_m2_abs},
                    {"area_m2_signed", area_m2_signed}
                }},
                {"geometry", geom}
            }
        })}
    };

    std::cout << fc.dump(2) << "\n";
    return 0;
}
