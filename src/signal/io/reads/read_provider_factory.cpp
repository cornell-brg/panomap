#include "signal/io/reads/read_provider_factory.hpp"

#include <algorithm>
#include <filesystem>

#include "signal/io/reads/pod5_provider.hpp"
#include "signal/io/reads/slow5_provider.hpp"
#include "core/util/logging.hpp"

namespace piru::io {

namespace {

std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}

std::string extension_of(const std::string& path) {
  const auto pos = path.find_last_of('.');
  if (pos == std::string::npos) return "";
  return to_lower(path.substr(pos + 1));
}

}  // namespace

ReadProviderPtr make_read_provider(const std::string& path) {
  // Directories: assume slow5/blow5 (Slow5Provider walks the tree).
  std::error_code ec;
  if (std::filesystem::is_directory(path, ec)) {
    return std::make_unique<Slow5Provider>(path);
  }

  const std::string ext = extension_of(path);
  if (ext == "slow5" || ext == "blow5") {
    return std::make_unique<Slow5Provider>(path);
  }
  if (ext == "pod5") {
    LOG_ERROR("POD5 reading not yet supported; please convert to slow5/blow5 (" + path + ")");
    return nullptr;
  }
  return nullptr;
}

}  // namespace piru::io
