# projects/008 tests

These tests run the `008` CLI via CTest and validate JSON outputs.

## Run
From projects/008:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Notes

time_ms is intentionally ignored (non-deterministic).

All tests assume planar/cartesian coordinates (not lon/lat degrees).

## 7) 実行方法

```bash
cd projects/008
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```
