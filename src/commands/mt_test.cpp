#include "commands/mt_test.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <thread>
#include <vector>

#include "cli/parse.hpp"
#include "concurrency/executor.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"

namespace {

// mt-test: simple parallel vector addition to exercise the executor and timing.
// Options:
//   -t/--threads: number of parallel tasks/threads (default 8)
//   -n/--size:    vector length (default 1,000,000)
//   -p/--profile: emit timing summary
// The command fills A=1.0, B=2.0 and computes C=A+B, then logs a checksum.

std::size_t parse_size(const std::map<std::string, std::string>& values, const std::string& key,
                       std::size_t default_value) {
  auto it = values.find(key);
  if (it == values.end()) return default_value;
  try {
    return static_cast<std::size_t>(std::stoul(it->second));
  } catch (...) {
    return default_value;
  }
}

}  // namespace

int handle_mt_test(const std::vector<std::string>& args) {
  std::size_t threads = 8;
  std::size_t size = 1'000'000;
  bool profile = false;

  piru::cli::Parsed parsed;
  piru::cli::ParseConfig config;
  config.usage = "Usage: piru mt-test [options]";
  config.options = {
      {'h', "help", false, "Show this help message"},
      {'t', "threads", true, "Number of threads/tasks (default 8)"},
      {'n', "size", true, "Vector length (default 1,000,000)"},
      {'p', "profile", false, "Emit timing profile (tree)"},
  };
  config.on_error = [](const std::string&) { std::cerr << "mt-test: invalid option\n"; };

  if (!piru::cli::parse_args(args, config, parsed)) {
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.values.count("help")) {
    piru::cli::print_help(config, std::cout);
    return 0;
  }

  if (auto it = parsed.values.find("threads"); it != parsed.values.end()) {
    threads = parse_size({{"threads", it->second}}, "threads", threads);
  }
  if (auto it = parsed.values.find("size"); it != parsed.values.end()) {
    size = parse_size({{"size", it->second}}, "size", size);
  }
  profile = parsed.values.count("profile") > 0;

  auto exec = piru::concurrency::make_executor(static_cast<int>(threads));

  PIRU_PROFILE_START(profile, "mt-test");

  std::vector<double> a(size, 1.0);
  std::vector<double> b(size, 2.0);
  std::vector<double> c(size, 0.0);

  const std::size_t grain = 4096;
  exec->parallel_for(0, size, grain, [&](std::size_t i) {
    const std::size_t end = std::min(i + grain, size);
    for (std::size_t j = i; j < end; ++j) {
      c[j] = a[j] + b[j];
    }
  });

  // Simple checksum to ensure work is kept.
  double checksum = 0.0;
  for (double v : c) checksum += v;

  LOG_INFO("mt-test completed (threads=" + std::to_string(threads) +
           ", size=" + std::to_string(size) + ", backend=" + exec->backend_name() +
           ", max_concurrency=" + std::to_string(exec->max_concurrency()) +
           ", checksum=" + std::to_string(checksum) + ").");

  PIRU_PROFILE_STOP(profile, "mt-test");
  if (profile) piru::timing::report(std::cerr);
  return 0;
}
