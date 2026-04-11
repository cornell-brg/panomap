// SPDX-License-Identifier: MIT
/*
 * slow5_provider.hpp - ReadProvider implementation backed by SLOW5 files
 */

#pragma once

#include <string>
#include <vector>

extern "C" {
#include "slow5/slow5.h"
}

#include "io/reads/read_provider.hpp"

namespace piru::io {

class Slow5Provider : public ReadProvider {
public:
  explicit Slow5Provider(std::string filename);
  ~Slow5Provider() override;

  bool get_next(RawRead& read) override;
  void reset() override;
  std::string get_format_name() const override { return "slow5"; }

  // Expose input files (used for per-file output routing).
  const std::vector<std::string>& input_files() const { return filenames_; }

private:
  void open();
  void close();
  bool advance_to_next_file();

  static std::vector<std::string> collect_input_files(const std::string& path);

  std::string input_source_;
  std::vector<std::string> filenames_;
  std::size_t current_file_index_{0};
  slow5_file_t* fp_{nullptr};
};

}  // namespace piru::io
