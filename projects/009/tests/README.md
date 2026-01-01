# projects/009 tests

These tests run `009` via CTest and validate JSON output.

## Build & Run
From projects/009:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

### Labels

fast: convex, hemisphere-friendly cases (expect mismatch=0)

experimental: intentionally breaks CGAL demo assumptions (expect mismatch > 0)

Run only fast tests:

```bash
ctest --test-dir build -L fast --output-on-failure
```

Run only experimental tests:

```bash
ctest --test-dir build -L experimental --output-on-failure
```

### Note

`009` uses Monte Carlo containment checks. Keep sample counts small for CI/ctest speed.
