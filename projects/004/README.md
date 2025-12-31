# projects/004

## About
GDAL/OGR で GeoJSON を開き、レイヤ情報・SRS・Extent などの概要を表示します。

## Dependencies
- GDAL

## Build & Run
```bash
cmake -S projects/004 -B build/004
cmake --build build/004
./build/004/004 path/to/input.geojson
```

## Notes
- 出力は標準出力のみ。GeoJSONの内容は変更しません。
