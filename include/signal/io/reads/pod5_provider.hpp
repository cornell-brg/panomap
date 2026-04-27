// SPDX-License-Identifier: MIT
// POD5 read provider stub.

#pragma once

#include <memory>
#include <string>

#include "signal/io/reads/read_provider.hpp"

namespace piru::io {

class Pod5Provider : public ReadProvider {
public:
  explicit Pod5Provider(const std::string& path);
  ~Pod5Provider() override = default;

  bool get_next(RawRead& read) override;
  void reset() override;
  std::string get_format_name() const override;

private:
  std::string path_;
  bool warned_{false};
};

using Pod5ProviderPtr = std::unique_ptr<Pod5Provider>;

}  // namespace piru::io
