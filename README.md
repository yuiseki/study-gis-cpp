# study-gis-cpp

C++ で GIS 周辺の基礎的な処理を小さなサンプルに分けて検証するためのリポジトリです。
各 `projects/00x` は独立した最小構成で、CMake でビルドできます。

## Projects

- `projects/001` — **PROJ + GEOS + nlohmann_json**
  - EPSG:4326 の点を EPSG:3857 に投影し、メートル単位のバッファを作成。GeoJSON を標準出力。
- `projects/002` — **GeographicLib + nlohmann_json**
  - WGS84 楕円体上の測地円を多角形化し、周長・面積付き GeoJSON を出力。
- `projects/003` — **GeographicLib + nlohmann_json**
  - GeoJSON（Polygon/MultiPolygon）から測地周長・面積を計算し、JSON で結果を出力。
- `projects/004` — **GDAL**
  - GeoJSON のレイヤ情報・SRS・Extent などのメタ情報を表示。
- `projects/005` — **GDAL + GeographicLib + nlohmann_json**
  - GeoJSON を読み込み、Polygon/MultiPolygon の測地面積・周長を属性として付与して出力。
- `projects/006` — **GDAL + PROJ**
  - ベクタ GeoJSON を再投影して書き出し（必要なら --src で元SRSを指定）。
  - 軸順の扱いを切り替えたり、PROJ/GDAL の結果を1点で簡易比較できる。
- `projects/007` — **GeographicLib + Boost.Geometry + PROJ + nlohmann_json**
  - 同じ2点の距離を複数の計算戦略で比較し、基準（GeographicLib s12）との差分を JSON で出力。
- `projects/008` — **GEOS + Boost.Geometry + nlohmann_json**
  - 平面ポリゴンの boolean（union / intersection / difference / symdiff）結果を両ライブラリで比較。
  - 結果の面積・境界長・頂点数・実行時間などを JSON で出力。
- `projects/009` — **S2 + CGAL (Nef_S2) + nlohmann_json**
  - 球面ポリゴン boolean を S2 と CGAL で比較し、Monte Carlo で包含/面積の差を推定。
  - CGAL 側は「球面凸ポリゴン + 半球内」前提の最小実装。

## Build
各プロジェクトは個別にビルドします（C++17）。

```bash
cmake -S projects/001 -B build/001
cmake --build build/001
```

実行例（001）:

```bash
./build/001/001 --lon 139.767125 --lat 35.681236 --radius 1000 --quad 8
```

他のプロジェクトも `001` を `002`〜`009` に置き換えて同様にビルドします。

## Tests
自動テストは `projects/008` と `projects/009` にあります（CTest）。

```bash
cmake -S projects/008 -B build/008 -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/008 -j
ctest --test-dir build/008 --output-on-failure
```

`projects/009` も同様に `008` を `009` に置き換えて実行してください。

## Notes
- 依存ライブラリは各プロジェクトの `CMakeLists.txt` に記載しています。
- 入力/出力は EPSG:4326（lon/lat 度）前提のものが多いです。各 `projects/00x/README.md` を参照してください。
- `datasets/` は入力データ置き場、`forks/` は参照資料置き場です。
