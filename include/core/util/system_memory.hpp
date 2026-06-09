#pragma once

#include <cstddef>
#include <string>

namespace panomap::util {

// Query system memory from /proc/meminfo (Linux).
// Returns 0 on failure or unsupported platform.
std::size_t getAvailableMemoryBytes();

// Compute a memory budget: min(max_memory, fraction * available).
// If max_memory == 0, uses fraction * available.
// Default fraction = 0.8 (leave 20% headroom for OS).
std::size_t computeMemoryBudget(std::size_t max_memory = 0, double fraction = 0.8);

// Format bytes as human-readable string (e.g., "12.3 GB").
std::string formatBytes(std::size_t bytes);

// Format large number with commas (e.g., 1234567 -> "1,234,567").
std::string formatCount(std::size_t count);

// Format seconds as human-readable duration (e.g., "2m 14s", "0.3s").
std::string formatDuration(double seconds);

// Print a progress bar to stderr: [████████░░░░░░░░] 42/100
// Updates in-place (carriage return). Call with done=total for final newline.
void printProgress(std::size_t done, std::size_t total, const std::string& label = "");

}  // namespace panomap::util
