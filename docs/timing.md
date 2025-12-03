# Timing Instrumentation

This repository includes a lightweight timing collector to profile code blocks without a full profiler. Timing is off by default and opt-in per command.

## How it works
- Start/stop pairs accumulate time per label: `timing::start("label"); ... timing::stop("label");`.
- Labels are thread-safe and accumulate wall time across multiple calls; use outside tight inner loops.
- Reporting prints each label with total wall time and auto-scaled units (µs/ms/s).

## Enabling profiling
- Commands with profiling flags: `index-vg`, `index-dbg`, `map`, `mt-test` (use `--profile`/`-p`).
- Example:
  ```bash
  ./piru/build/piru index-vg --reads R --graph G --profile
  ```
- When enabled, the command wraps its work with the profiling macros and prints the timing summary to stderr before the footer.

## Adding instrumentation
- Include `util/timing.hpp` and wrap work with `timing::start("label"); ... timing::stop("label");` or the macros `PIRU_PROFILE_START/STOP(enabled, label)`.
- Keep labels stable and descriptive (e.g., `load-reads`, `build-index`, `map-batch`); repeated calls accumulate totals.
- Only emit `timing::report(std::cerr)` when profiling is enabled to avoid noise in normal runs.
- You can start/stop the same label multiple times; totals accumulate.
- Example:
  ```cpp
  #include "util/timing.hpp"
  #include <thread>

  void do_work(bool profile) {
      PIRU_PROFILE_START(profile, "work");
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
      PIRU_PROFILE_STOP(profile, "work");

      // Accumulate more time under the same label.
      PIRU_PROFILE_START(profile, "work");
      std::this_thread::sleep_for(std::chrono::milliseconds(25));
      PIRU_PROFILE_STOP(profile, "work");
      if (profile) piru::timing::report(std::cerr);
  }
  ```

## Notes
- Wall-clock time is recorded per label; use the footer for overall real/CPU/peak RSS.
- Minimal overhead when disabled; avoid wrapping ultra-tight inner loops unless measuring them explicitly.
