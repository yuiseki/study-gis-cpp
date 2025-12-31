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

## Build
各プロジェクトは個別にビルドします。

```bash
cmake -S projects/001 -B build/001
cmake --build build/001
```

## Notes
- 依存ライブラリは各プロジェクトの `CMakeLists.txt` に記載しています。
- 入力/出力は EPSG:4326（lon/lat 度）前提のものが多いです。各 `projects/00x/README.md` を参照してください。
