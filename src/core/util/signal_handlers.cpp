#include "core/util/signal_handlers.hpp"

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <unistd.h>

#include "version.hpp"

namespace panomap {
namespace {

void segv_handler(int signum) {
  // Minimal async-signal-safe logging.
  const char header[] = "\n\n[panic] Caught SIGSEGV; dumping backtrace.\n";
  [[maybe_unused]] ssize_t _ = write(STDERR_FILENO, header, sizeof(header) - 1);

  const char version_prefix[] = "[panic] Version: ";
  [[maybe_unused]] ssize_t _2 = write(STDERR_FILENO, version_prefix, sizeof(version_prefix) - 1);
  [[maybe_unused]] ssize_t _3 =
      write(STDERR_FILENO, PIRU_VERSION_STRING, sizeof(PIRU_VERSION_STRING) - 1);
  const char nl[] = "\n";
  [[maybe_unused]] ssize_t _4 = write(STDERR_FILENO, nl, sizeof(nl) - 1);

  void* trace[64];
  int size = backtrace(trace, static_cast<int>(sizeof(trace) / sizeof(void*)));
  backtrace_symbols_fd(trace, size, STDERR_FILENO);

  const char tail[] = "\n[panic] Aborting due to segmentation fault.\n";
  [[maybe_unused]] ssize_t _5 = write(STDERR_FILENO, tail, sizeof(tail) - 1);
  std::_Exit(signum);
}

}  // namespace

void install_signal_handlers() { std::signal(SIGSEGV, segv_handler); }

}  // namespace panomap
