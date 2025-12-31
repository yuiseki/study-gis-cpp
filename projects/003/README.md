# projects/003

## About
GeoJSON（Polygon / MultiPolygon）を読み込み、
GeographicLib で測地周長・面積を計算して JSON で出力します。
FeatureCollection / Feature / Geometry の各入力形式に対応します。

## Dependencies
- GeographicLib
- nlohmann_json

## Build & Run
```bash
cmake -S projects/003 -B build/003
cmake --build build/003
./build/003/003 path/to/input.geojson > metrics.json
```

## Options
- `--feature <i>`: FeatureCollection 内の特定 index のみ計算（0-based）
- `--compact`: JSON を1行で出力

## Notes
- 入力は EPSG:4326（[lon, lat]）前提です。
- Polygon/MultiPolygon 以外は警告のみ出力します。
