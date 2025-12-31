# projects/002

## About
GeographicLib（WGS84）で測地円（lon/lat中心・半径m）を多角形化し、
GeoJSON FeatureCollection を標準出力します。周長・面積も計算します。

## Dependencies
- GeographicLib
- nlohmann_json

## Build & Run
```bash
cmake -S projects/002 -B build/002
cmake --build build/002
./build/002/002 --lon 139.767125 --lat 35.681236 --radius 1000 --steps 360 > geodesic.geojson
```

## Notes
- `--steps` は方位角の分割数（多いほど滑らか）。最低 8 にクランプされます。
- 出力は EPSG:4326 の GeoJSON（[lon, lat]）。
