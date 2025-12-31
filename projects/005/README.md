# projects/005

## About
GDAL で GeoJSON を読み込み、Polygon/MultiPolygon の測地面積・周長を
GeographicLib で計算して新しい GeoJSON に書き出します。
計算結果は属性（`geod_area_m2` など）として追加されます。

## Dependencies
- GDAL
- GeographicLib
- nlohmann_json

## Build & Run
```bash
cmake -S projects/005 -B build/005
cmake --build build/005
./build/005/005 input.geojson output.geojson
```

stdout に出す場合:
```bash
./build/005/005 input.geojson - > out.geojson
```

## Notes
- 入力は EPSG:4326（lon/lat 度）前提です。
- Polygon/MultiPolygon 以外は計算対象外です。
- 追加フィールド: `geod_area_m2`, `geod_perim_outer_m`, `geod_perim_all_m`
