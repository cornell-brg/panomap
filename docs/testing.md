# Testing Infrastructure

This project uses CTest + doctest for unit tests, with optional FetchContent to vendor doctest when not installed.

## Layout
- Tests live under `tests/` with a single CMake target `piru_tests`.
- Minimal samples: `tests/test_main.cpp` (doctest entry) and `tests/test_metrics.cpp`.

## Configuration
- CMake options:
  - `BUILD_TESTING` (default ON): enable tests.
  - `PIRU_FETCH_DOCTEST` (default ON): FetchContent doctest if not found on the system.
- `tests/CMakeLists.txt` will:
  - Locate or fetch doctest (`doctest::doctest` target).
  - Build `piru_tests` and register tests via `doctest_discover_tests` when available (falls back to a single ctest invocation otherwise).
  - Create a `check` custom target to run ctest.

## Running tests
```bash
cmake -S piru -B piru/build -DBUILD_TESTING=ON
cmake --build piru/build
ctest --test-dir piru/build --output-on-failure   # or: cmake --build piru/build --target check
# if your generator supports it: make -C piru/build test
```

## Adding tests
- Add new `*_test.cpp` files to `tests/` and list them in `tests/CMakeLists.txt` under the `piru_tests` sources.
- Use doctest TEST_CASE/SCENARIO macros; include headers from `include/` and generated headers via `${CMAKE_BINARY_DIR}/generated`.
- Keep tests self-contained and fast; prefer unit-level coverage for library code.
- Example:
  ```cpp
  #include <doctest/doctest.h>
  #include "util/metrics.hpp"

  TEST_CASE("realtime monotonic") {
      auto t1 = piru::realtime();
      auto t2 = piru::realtime();
      CHECK(t2 >= t1);
  }
  ```
