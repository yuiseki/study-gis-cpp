# Repository Guidelines

## Project Structure & Module Organization
- `projects/00x/` contains small, self-contained C++ examples. Each subproject has its own `CMakeLists.txt` and `main.cpp`.
- `projects/00x/*.geojson` are sample outputs or fixtures for that example.
- `datasets/` is reserved for input data used by examples.
- `forks/` holds vendored or reference material.

## Build, Test, and Development Commands
Each project builds independently with CMake (C++17). Example for project 001:

```bash
cmake -S projects/001 -B build/001
cmake --build build/001
./build/001/001 --lon 139.767125 --lat 35.681236 --radius 1000 --quad 8
```

Other projects follow the same pattern (replace `001` with `002`–`005`).

## Coding Style & Naming Conventions
- C++17 is required (`CMAKE_CXX_STANDARD 17`).
- Indentation: 4 spaces; no tabs.
- Prefer clear, small functions and explicit error checks (see `projects/001/main.cpp`).
- Executable names match the project folder number (e.g., `add_executable(003 ...)`).

## Testing Guidelines
- No automated test framework is configured.
- Validate by running the executable and checking GeoJSON output or comparing with sample files under `projects/00x/`.
- If you add tests, document how to run them in this file and keep them per-project.

## Commit & Pull Request Guidelines
- Commit messages are short, imperative, and scoped (e.g., “Add projects/005”).
- PRs should include: a brief summary, the project number affected, build/run commands used, and sample output (e.g., GeoJSON snippet or file path). Screenshots are optional unless visual output is involved.

## Dependencies & Configuration Notes
- Project-specific dependencies are declared in each `CMakeLists.txt`:
  - `001`: PROJ, GEOS, nlohmann_json
  - `002`/`003`: GeographicLib, nlohmann_json
  - `004`: GDAL
  - `005`: GDAL, GeographicLib, nlohmann_json
- Many dependencies are discovered via CMake config packages or `pkg-config`; ensure development headers are installed.
