#include "signal/io/reads/pod5_provider.hpp"

#include "core/util/logging.hpp"

namespace piru::io {

Pod5Provider::Pod5Provider(const std::string& path) : path_(path) {}

bool Pod5Provider::get_next(RawRead& read) {
  (void)read;
  if (!warned_) {
    LOG_ERROR("POD5 reading not supported yet; convert to slow5/blow5 ('" + path_ + "')");
    warned_ = true;
  }
  return false;
}

void Pod5Provider::reset() {
  // No-op stub.
}

std::string Pod5Provider::get_format_name() const { return "pod5"; }

}  // namespace piru::io
