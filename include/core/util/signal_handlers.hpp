#pragma once

#include <csignal>

namespace panomap {

// Install handlers for fatal signals (currently SIGSEGV) to aid debugging.
void install_signal_handlers();

}  // namespace panomap
