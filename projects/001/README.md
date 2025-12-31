# projects/001

## About

PROJ + GEOS + nlohmann_json で、

- EPSG:4326（lon/lat）入力
- PROJでEPSG:3857（m）へ投影
- GEOSで 半径Rメートル のバッファ（ポリゴン）生成
- PROJでEPSG:4326へ戻す
- nlohmann_jsonで GeoJSON FeatureCollection を標準出力

## Build & Run
```bash
cmake -S projects/001 -B build/001
cmake --build build/001
./build/001/001 --lon 139.767125 --lat 35.681236 --radius 1000 --quad 8 > webmerc.geojson
```

## Options
- `--lon <deg>`: 経度（度）
- `--lat <deg>`: 緯度（度）
- `--radius <m>`: バッファ半径（メートル）
- `--quad <n>`: 円弧分割数（大きいほど滑らか）
