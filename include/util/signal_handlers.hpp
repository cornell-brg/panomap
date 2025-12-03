#pragma once

#include <csignal>

namespace piru {

// Install handlers for fatal signals (currently SIGSEGV) to aid debugging.
void install_signal_handlers();

}  // namespace piru
