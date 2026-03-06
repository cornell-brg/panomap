#include "util/timing.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <unordered_map>

#include "util/metrics.hpp"

namespace piru::timing {
namespace {

struct Entry {
  double wall_ms{0.0};
  bool running{false};
  std::chrono::steady_clock::time_point start_wall;
};

class Collector {
public:
  static Collector& instance() {
    static Collector c;
    return c;
  }

  void start(const std::string& label) {
    std::lock_guard<std::mutex> lock(mu_);
    auto& e = entries_[label];
    if (e.running) return;
    e.running = true;
    e.start_wall = std::chrono::steady_clock::now();
  }

  void stop(const std::string& label) {
    std::lock_guard<std::mutex> lock(mu_);
    auto it = entries_.find(label);
    if (it == entries_.end() || !it->second.running) return;
    auto end_wall = std::chrono::steady_clock::now();
    it->second.wall_ms += std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
                              end_wall - it->second.start_wall)
                              .count();
    it->second.running = false;
    has_data_ = true;
  }

  void print(std::ostream& os) {
    std::lock_guard<std::mutex> lock(mu_);
    if (!has_data_) return;
    os << "\n[timing]\n";
    for (const auto& [label, e] : entries_) {
      auto [val, unit] = pretty_time(e.wall_ms);
      os << "  " << label << ": " << std::fixed << std::setprecision(3) << val << " " << unit
         << "\n";
    }
  }

private:
  static std::pair<double, const char*> pretty_time(double ms) {
    if (ms >= 1000.0) {
      return {ms / 1000.0, "s"};
    }
    if (ms >= 1.0) {
      return {ms, "ms"};
    }
    return {ms * 1000.0, "µs"};
  }

  std::unordered_map<std::string, Entry> entries_;
  std::mutex mu_;
  bool has_data_{false};
};

}  // namespace

void start(const std::string& label) { Collector::instance().start(label); }

void stop(const std::string& label) { Collector::instance().stop(label); }

void report(std::ostream& os) { Collector::instance().print(os); }

}  // namespace piru::timing
