#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/number_utils.h>

#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2polygon.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using json = nlohmann::json;

// -------------------------
// CGAL Nef_S2 typedefs
// -------------------------
using RT = CGAL::Exact_integer;
using Kernel = CGAL::Homogeneous<RT>;
using NefS2 = CGAL::Nef_polyhedron_S2<Kernel>;
using Sphere_point = NefS2::Sphere_point;
using Sphere_circle = NefS2::Sphere_circle;
using Object_handle = NefS2::Object_handle;
using SFace_const_handle = NefS2::SFace_const_handle;

// -------------------------
// 小さなユーティリティ
// -------------------------
static double deg2rad(double d) { return d * M_PI / 180.0; }
static double rad2deg(double r) { return r * 180.0 / M_PI; }

struct Args {
  std::string op = "union"; // union|intersection|difference|xor
  std::string a_wkt;
  std::string b_wkt;
  int samples = 20000;
  uint64_t seed = 1;
  double expand_deg = 1.0;  // サンプリング矩形を少し広げる
  bool emit_loops = true;   // S2結果のループ座標をJSONに含める
};

// -------------------------
// ざっくり引数パース
// -------------------------
static Args parse_args(int argc, char** argv) {
  Args a;
  auto need = [&](int& i, const char* name) -> std::string {
    if (i + 1 >= argc) {
      std::cerr << "Missing value for " << name << "\n";
      std::exit(2);
    }
    return std::string(argv[++i]);
  };

  for (int i = 1; i < argc; i++) {
    std::string k = argv[i];
    if (k == "--op") a.op = need(i, "--op");
    else if (k == "--a") a.a_wkt = need(i, "--a");
    else if (k == "--b") a.b_wkt = need(i, "--b");
    else if (k == "--samples") a.samples = std::stoi(need(i, "--samples"));
    else if (k == "--seed") a.seed = static_cast<uint64_t>(std::stoull(need(i, "--seed")));
    else if (k == "--expand-deg") a.expand_deg = std::stod(need(i, "--expand-deg"));
    else if (k == "--no-emit-loops") a.emit_loops = false;
    else if (k == "--help") {
      std::cout
        << "Usage:\n"
        << "  009 --op union|intersection|difference|xor \\\n"
        << "      --a 'POLYGON((lon lat, ...))' \\\n"
        << "      --b 'POLYGON((lon lat, ...))' \\\n"
        << "      [--samples N] [--seed S] [--expand-deg D] [--no-emit-loops]\n\n"
        << "Assumption:\n"
        << "  - Coordinates are lon/lat degrees on the sphere (not planar).\n"
        << "  - For CGAL(Nef_S2) path in this minimal demo, polygons should be spherical-convex\n"
        << "    and contained within a hemisphere (so we can represent them as intersection of hemispheres).\n";
      std::exit(0);
    } else {
      std::cerr << "Unknown option: " << k << "\n";
      std::exit(2);
    }
  }

  if (a.a_wkt.empty() || a.b_wkt.empty()) {
    std::cerr << "Both --a and --b are required.\n";
    std::exit(2);
  }
  if (a.samples <= 0) a.samples = 1;
  return a;
}

// -------------------------
// WKT POLYGON((x y, ...)) の最小パーサ
// - 穴・MULTIなどは未対応（009は最小比較用）
// -------------------------
static std::vector<std::pair<double,double>> parse_wkt_polygon_lonlat(const std::string& wkt) {
  // かなり雑ですが学習プロジェクト用途なら十分
  auto up = wkt;
  // trim
  auto ltrim = [](std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char c){ return !std::isspace(c); }));
  };
  auto rtrim = [](std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c){ return !std::isspace(c); }).base(), s.end());
  };
  ltrim(up); rtrim(up);

  // 先頭がPOLYGONか確認（大文字小文字はざっくり）
  std::string head = up.substr(0, std::min<size_t>(7, up.size()));
  for (auto& c : head) c = std::toupper(static_cast<unsigned char>(c));
  if (head != "POLYGON") {
    std::cerr << "Only POLYGON is supported in this minimal parser.\n";
    std::exit(2);
  }

  auto p0 = up.find("((");
  auto p1 = up.rfind("))");
  if (p0 == std::string::npos || p1 == std::string::npos || p1 <= p0 + 1) {
    std::cerr << "Invalid POLYGON WKT.\n";
    std::exit(2);
  }
  std::string body = up.substr(p0 + 2, p1 - (p0 + 2));

  std::vector<std::pair<double,double>> pts;
  std::stringstream ss(body);
  std::string token;
  while (std::getline(ss, token, ',')) {
    ltrim(token); rtrim(token);
    if (token.empty()) continue;

    std::stringstream ps(token);
    double x=0.0, y=0.0;
    ps >> x >> y;
    if (!ps) {
      std::cerr << "Failed to parse coordinate: " << token << "\n";
      std::exit(2);
    }
    pts.emplace_back(x, y); // (lon, lat)
  }

  if (pts.size() < 4) {
    std::cerr << "Polygon ring must have at least 4 points (including closure).\n";
    std::exit(2);
  }

  // 閉じていなければ閉じる
  if (pts.front().first != pts.back().first || pts.front().second != pts.back().second) {
    pts.push_back(pts.front());
  }

  // S2Loopは最後=最初の重複頂点を好まないので、内部では落とす
  // ここではリングとしては閉じたまま返し、S2側で落とします。
  return pts;
}

// -------------------------
// S2 polygon construction
// - S2は「球面」前提。S2Loopの内側は左手側（CCW）と定義。 :contentReference[oaicite:3]{index=3}
// - コンストラクタは vector<unique_ptr<S2Loop>> を受け取れる形が一般的（Issue例）。 :contentReference[oaicite:4]{index=4}
// -------------------------
static void make_s2_polygon_from_lonlat_ring(const std::vector<std::pair<double,double>>& ring_lonlat_closed, S2Polygon* poly_out) {
  std::vector<S2Point> verts;
  verts.reserve(ring_lonlat_closed.size());

  // 最後の閉じ点は落とす（S2Loopは暗黙に最後->最初を結ぶ）
  for (size_t i = 0; i + 1 < ring_lonlat_closed.size(); i++) {
    const auto [lon, lat] = ring_lonlat_closed[i];
    S2Point p(S2LatLng::FromDegrees(lat, lon));
    verts.push_back(p);
  }

  auto loop = std::make_unique<S2Loop>(verts);
  loop->Normalize(); // 面積 <= 2π 側を内側に（「小さい方の多角形」へ寄せる）

  poly_out->Init(std::move(loop));
}

static json s2_polygon_loops_json(const S2Polygon& poly) {
  json out = json::array();
  for (int i = 0; i < poly.num_loops(); i++) {
    const S2Loop* loop = poly.loop(i);
    json ring = json::array();
    for (int j = 0; j < loop->num_vertices(); j++) {
      S2LatLng ll(loop->vertex(j));
      ring.push_back(json::array({ ll.lng().degrees(), ll.lat().degrees() })); // [lon,lat]
    }
    // GeoJSON風に閉じる
    if (loop->num_vertices() > 0) {
      S2LatLng ll0(loop->vertex(0));
      ring.push_back(json::array({ ll0.lng().degrees(), ll0.lat().degrees() }));
    }
    out.push_back(json{
      {"is_hole", loop->is_hole()},
      {"vertices", loop->num_vertices()},
      {"coordinates", ring}
    });
  }
  return out;
}

// 球面上の境界長（単位：ラジアン；単位球の弧長）
static double s2_polygon_boundary_length_rad(const S2Polygon& poly) {
  double sum = 0.0;
  for (int i = 0; i < poly.num_loops(); i++) {
    const S2Loop* loop = poly.loop(i);
    const int n = loop->num_vertices();
    for (int j = 0; j < n; j++) {
      const S2Point& a = loop->vertex(j);
      const S2Point& b = loop->vertex((j + 1) % n);
      sum += a.Angle(b); // 球面上の中心角
    }
  }
  return sum;
}

// -------------------------
// CGAL Nef_S2 construction (最小版)
// 「球面凸ポリゴン（半球内）」を、各辺が作る“半球”の交差として表現する。
// - Nef_S2は「球面上の大円弧（shalfedge）や面（sface）を持つ複体」でBoolean演算を行う。 :contentReference[oaicite:5]{index=5}
//
// ※ 一般の非凸球面ポリゴンはこの方法では表現できません（そこを攻めるのは009の次段階で）。
// -------------------------
static Sphere_point lonlat_to_sphere_point_approx(double lon_deg, double lat_deg) {
  const double lon = deg2rad(lon_deg);
  const double lat = deg2rad(lat_deg);
  const double x = std::cos(lat) * std::cos(lon);
  const double y = std::cos(lat) * std::sin(lon);
  const double z = std::sin(lat);

  // Homogeneous<Exact_integer> に入れるため整数へ丸め（方向ベクトルとして使うのでスケールは無関係）
  constexpr long long SCALE = 1000000000LL;
  long long ix = llround(x * SCALE);
  long long iy = llround(y * SCALE);
  long long iz = llround(z * SCALE);

  // 万一ゼロベクトル回避
  if (ix == 0 && iy == 0 && iz == 0) ix = 1;
  return Sphere_point(ix, iy, iz);
}

static std::array<double,3> to_dvec(const Sphere_point& p) {
  // Homogeneous kernel: p.x(), p.y(), p.z() are RT
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  double z = CGAL::to_double(p.z());
  // 同次座標のwはここでは無視（Sphere_pointは方向とみなす）
  double norm = std::sqrt(x*x + y*y + z*z);
  if (norm == 0) return {0,0,0};
  return {x/norm, y/norm, z/norm};
}

static double dot3(const std::array<double,3>& a, const std::array<double,3>& b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static NefS2 make_cgal_nef_convex_from_lonlat_ring(const std::vector<std::pair<double,double>>& ring_lonlat_closed) {
  // 内点（近似）：頂点方向ベクトルの平均
  double sx=0, sy=0, sz=0;
  std::vector<Sphere_point> v;
  v.reserve(ring_lonlat_closed.size());

  for (size_t i = 0; i + 1 < ring_lonlat_closed.size(); i++) {
    const auto [lon, lat] = ring_lonlat_closed[i];
    Sphere_point sp = lonlat_to_sphere_point_approx(lon, lat);
    auto dv = to_dvec(sp);
    sx += dv[0]; sy += dv[1]; sz += dv[2];
    v.push_back(sp);
  }
  double n = std::sqrt(sx*sx + sy*sy + sz*sz);
  if (n == 0) { sx = 1; sy = 0; sz = 0; n = 1; }
  sx /= n; sy /= n; sz /= n;

  // centroid方向をSphere_point（近似）に
  constexpr long long SCALE = 1000000000LL;
  Sphere_point centroid(llround(sx*SCALE), llround(sy*SCALE), llround(sz*SCALE));
  auto centroid_d = to_dvec(centroid);

  // 全空間から開始して半球で切っていく
  NefS2 nef(NefS2::COMPLETE);

  const size_t m = v.size();
  for (size_t i = 0; i < m; i++) {
    const Sphere_point& p = v[i];
    const Sphere_point& q = v[(i + 1) % m];

    Sphere_circle c(p, q);
    // cの“左半球”の極（orthogonal_pole）が、centroidと同じ側になるように向きを揃える
    Sphere_point pole = c.orthogonal_pole();
    auto pole_d = to_dvec(pole);
    if (dot3(pole_d, centroid_d) < 0) {
      c = c.opposite();
    }

    NefS2 half(c, NefS2::INCLUDED);
    nef = nef * half; // intersection
  }

  return nef;
}

// CGAL Nef_S2の点包含（境界ヒットは確率0なので最小実装では雑に扱う）
static bool cgal_nef_contains(const NefS2& nef, const Sphere_point& p) {
  Object_handle h = nef.locate(p);

  SFace_const_handle sf;
  if (CGAL::assign(sf, h)) {
    // sfaceには選択マーク（集合に含まれるか）がある :contentReference[oaicite:6]{index=6}
    return sf->mark();
  }
  // エッジ/頂点上など：ここでは「含まれる」と扱う（サンプル点ではほぼ起きない）
  return true;
}

// -------------------------
// サンプリング矩形（lat/lon）を作る（今回は日付変更線跨ぎは扱わない最小版）
// -------------------------
static void compute_lonlat_bbox(
  const std::vector<std::pair<double,double>>& a,
  const std::vector<std::pair<double,double>>& b,
  double expand_deg,
  double& lon_min, double& lon_max,
  double& lat_min, double& lat_max
) {
  lon_min = +std::numeric_limits<double>::infinity();
  lon_max = -std::numeric_limits<double>::infinity();
  lat_min = +std::numeric_limits<double>::infinity();
  lat_max = -std::numeric_limits<double>::infinity();

  auto upd = [&](const std::vector<std::pair<double,double>>& ring) {
    for (size_t i = 0; i + 1 < ring.size(); i++) {
      lon_min = std::min(lon_min, ring[i].first);
      lon_max = std::max(lon_max, ring[i].first);
      lat_min = std::min(lat_min, ring[i].second);
      lat_max = std::max(lat_max, ring[i].second);
    }
  };
  upd(a); upd(b);

  lon_min -= expand_deg;
  lon_max += expand_deg;
  lat_min -= expand_deg;
  lat_max += expand_deg;

  lon_min = std::max(lon_min, -180.0);
  lon_max = std::min(lon_max, +180.0);
  lat_min = std::max(lat_min, -90.0);
  lat_max = std::min(lat_max, +90.0);
}

// 矩形の球面面積（単位球の立体角 = steradian）
static double rect_area_steradians(double lon_min_deg, double lon_max_deg, double lat_min_deg, double lat_max_deg) {
  double lon1 = deg2rad(lon_min_deg);
  double lon2 = deg2rad(lon_max_deg);
  double lat1 = deg2rad(lat_min_deg);
  double lat2 = deg2rad(lat_max_deg);
  // 単位球上の緯度経度矩形の面積: Δλ * (sin φ2 - sin φ1)
  return std::max(0.0, (lon2 - lon1) * (std::sin(lat2) - std::sin(lat1)));
}

int main(int argc, char** argv) {
  const Args args = parse_args(argc, argv);

  // 入力WKT -> (lon,lat)リング
  auto a_ring = parse_wkt_polygon_lonlat(args.a_wkt);
  auto b_ring = parse_wkt_polygon_lonlat(args.b_wkt);

  // S2 polygons
  S2Polygon s2a, s2b;
  make_s2_polygon_from_lonlat_ring(a_ring, &s2a);
  make_s2_polygon_from_lonlat_ring(b_ring, &s2b);

  // CGAL Nef polygons (convex+hemisphere assumption)
  NefS2 ca = make_cgal_nef_convex_from_lonlat_ring(a_ring);
  NefS2 cb = make_cgal_nef_convex_from_lonlat_ring(b_ring);

  // Boolean op
  auto t0 = std::chrono::steady_clock::now();

  S2Polygon s2out;
  if (args.op == "union") s2out.InitToUnion(&s2a, &s2b);
  else if (args.op == "intersection") s2out.InitToIntersection(&s2a, &s2b);
  else if (args.op == "difference") s2out.InitToDifference(&s2a, &s2b); // A - B
  else if (args.op == "xor") {
    S2Polygon u, i;
    u.InitToUnion(&s2a, &s2b);
    i.InitToIntersection(&s2a, &s2b);
    s2out.InitToDifference(&u, &i); // (A ∪ B) \ (A ∩ B)
  } else {
    std::cerr << "Unknown --op: " << args.op << "\n";
    return 2;
  }

  NefS2 cgalout;
  if (args.op == "union") cgalout = ca + cb;
  else if (args.op == "intersection") cgalout = ca * cb;
  else if (args.op == "difference") cgalout = ca - cb; // A - B
  else if (args.op == "xor") {
    NefS2 u = ca + cb;
    NefS2 i = ca * cb;
    cgalout = u - i;
  }

  auto t1 = std::chrono::steady_clock::now();
  double op_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  // サンプリング矩形
  double lon_min, lon_max, lat_min, lat_max;
  compute_lonlat_bbox(a_ring, b_ring, args.expand_deg, lon_min, lon_max, lat_min, lat_max);
  const double rect_sr = rect_area_steradians(lon_min, lon_max, lat_min, lat_max);

  // サンプリング（緯度のsinを一様にすると矩形上で面積一様）
  std::mt19937_64 rng(args.seed);
  std::uniform_real_distribution<double> u_lon(lon_min, lon_max);
  std::uniform_real_distribution<double> u_sin(std::sin(deg2rad(lat_min)), std::sin(deg2rad(lat_max)));

  auto t2 = std::chrono::steady_clock::now();

  int in_s2 = 0;
  int in_cgal = 0;
  int mismatch = 0;

  for (int i = 0; i < args.samples; i++) {
    const double lon = u_lon(rng);
    const double s = u_sin(rng);
    const double lat = rad2deg(std::asin(std::max(-1.0, std::min(1.0, s))));

    S2Point p(S2LatLng::FromDegrees(lat, lon));
    bool a = s2out.Contains(p);

    Sphere_point q = lonlat_to_sphere_point_approx(lon, lat);
    bool b = cgal_nef_contains(cgalout, q);

    if (a) in_s2++;
    if (b) in_cgal++;
    if (a != b) mismatch++;
  }

  auto t3 = std::chrono::steady_clock::now();
  double sample_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

  // 面積推定（矩形内の比率 * 矩形面積）
  const double s2_area_est_sr = rect_sr * (static_cast<double>(in_s2) / args.samples);
  const double cgal_area_est_sr = rect_sr * (static_cast<double>(in_cgal) / args.samples);

  // S2は真の（単位球）面積を返せる
  const double s2_area_sr = s2out.GetArea();
  const double s2_boundary_rad = s2_polygon_boundary_length_rad(s2out);

  // 便宜的に“地球半径”を決めてm換算（球モデル）。値は比較のためのスケールとして扱ってください。
  constexpr double EARTH_RADIUS_M = 6371008.8;
  const double s2_area_m2 = s2_area_sr * EARTH_RADIUS_M * EARTH_RADIUS_M;
  const double s2_boundary_m = s2_boundary_rad * EARTH_RADIUS_M;
  const double s2_area_est_m2 = s2_area_est_sr * EARTH_RADIUS_M * EARTH_RADIUS_M;
  const double cgal_area_est_m2 = cgal_area_est_sr * EARTH_RADIUS_M * EARTH_RADIUS_M;

  json out;
  out["assumption"] = "Input coordinates are lon/lat degrees on a sphere. (This is NOT planar WKT semantics.)";
  out["note"] = "CGAL path in this minimal demo builds polygons as intersection of hemispheres: requires spherical-convex polygons within a hemisphere.";
  out["op"] = args.op;
  out["timing_ms"] = {
    {"boolean_op_total", op_ms},
    {"sampling_total", sample_ms}
  };
  out["input"] = {
    {"a_wkt", args.a_wkt},
    {"b_wkt", args.b_wkt},
    {"samples", args.samples},
    {"seed", args.seed},
    {"sampling_bbox_deg", {
      {"lon_min", lon_min}, {"lon_max", lon_max},
      {"lat_min", lat_min}, {"lat_max", lat_max}
    }},
    {"sampling_bbox_area_sr", rect_sr}
  };

  out["s2"] = {
    {"ok", true},
    {"valid", s2out.IsValid()},
    {"num_loops", s2out.num_loops()},
    {"num_vertices", s2out.num_vertices()},
    {"area_sr", s2_area_sr},
    {"area_m2_sphere", s2_area_m2},
    {"boundary_length_rad", s2_boundary_rad},
    {"boundary_length_m_sphere", s2_boundary_m}
  };

  if (args.emit_loops) {
    out["s2"]["loops"] = s2_polygon_loops_json(s2out);
  }

  out["cgal_nef_s2"] = {
    {"ok", true},
    {"contains_via_locate_mark", true}
  };

  out["comparison"] = {
    {"containment_mismatch", {
      {"count", mismatch},
      {"rate", static_cast<double>(mismatch) / args.samples}
    }},
    {"area_estimate_from_sampling", {
      {"s2_sr", s2_area_est_sr},
      {"cgal_sr", cgal_area_est_sr},
      {"abs_diff_sr", std::abs(s2_area_est_sr - cgal_area_est_sr)},
      {"s2_m2_sphere", s2_area_est_m2},
      {"cgal_m2_sphere", cgal_area_est_m2},
      {"abs_diff_m2_sphere", std::abs(s2_area_est_m2 - cgal_area_est_m2)}
    }}
  };

  std::cout << out.dump(2) << "\n";
  return 0;
}
