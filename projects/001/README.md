# projects/001

## About

PROJ + GEOS + nlohmann_json で、

- EPSG:4326（lon/lat）入力
- PROJでEPSG:3857（m）へ投影
- GEOSで 半径Rメートル のバッファ（ポリゴン）生成
- PROJでEPSG:4326へ戻す
- nlohmann_jsonで GeoJSON FeatureCollection を標準出力

