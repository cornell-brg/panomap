# Concurrency Abstraction

This project hides the threading backend behind a small interface so implementations can be swapped without touching call sites.

## Interface
- Location: `include/concurrency/executor.hpp`
- `Executor` exposes `parallel_for(begin, end, grain, fn)` and `max_concurrency()`.
- `TaskGroup` supports lightweight task spawning + `wait()`.
- Factory helpers: `make_executor(int threads = -1)` and `make_task_group()`. A negative/zero thread hint lets the backend pick automatically.

## How it works
- `parallel_for` partitions `[begin, end)` into blocks of size `grain` and schedules each block on the backend (TBB `parallel_for` or a serial loop). Tune `grain` to balance overhead vs. work per task.
- `TaskGroup` wraps `tbb::task_group` when TBB is enabled; otherwise tasks run inline.
- `max_concurrency()` reports the backend’s concurrency (TBB arena or 1 for serial).

## Backends
- Default: oneTBB (`PIRU_USE_TBB=ON`), fetched via FetchContent by default. Build fails if requested but not found.
- Fallback: serial executor/task group when TBB is disabled or unavailable.
- Backend wiring lives in `src/concurrency/executor.cpp`; a TBB task_arena/task_group drives parallel work when enabled.

## Usage
- Include `concurrency/executor.hpp`, create an executor via `make_executor`, and call `parallel_for` with a sensible grain size.
- Example:
  ```cpp
  #include "concurrency/executor.hpp"
  #include <vector>

  void add_vectors() {
      auto exec = piru::concurrency::make_executor(4);  // hint 4 threads
      const std::size_t n = 1'000'000;
      std::vector<double> a(n, 1.0), b(n, 2.0), c(n, 0.0);
      const std::size_t grain = 4096;
      exec->parallel_for(0, n, grain, [&](std::size_t i) {
          const std::size_t end = std::min(i + grain, n);
          for (std::size_t j = i; j < end; ++j) {
              c[j] = a[j] + b[j];
          }
      });
  }
  ```
- The `mt-test` subcommand demonstrates parallel vector addition; adjust `-t/--threads` to see scaling effects.

## CMake Options
- `PIRU_USE_TBB` (default ON): enable the oneTBB backend.
- `PIRU_FETCH_TBB` (default ON): FetchContent oneTBB if not found; configuration fails if TBB is requested and unavailable.
