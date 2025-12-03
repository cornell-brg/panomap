#pragma once

#include <iosfwd>
#include <string>

namespace piru::timing {

// Start/stop timers by label. Multiple start/stop pairs accumulate.
void start(const std::string& label);
void stop(const std::string& label);

// Print collected timing summary (wall + CPU) to the given stream.
void report(std::ostream& os);

}  // namespace piru::timing

#define PIRU_PROFILE_START(enabled, label) \
    do {                                   \
        if (enabled) piru::timing::start(label); \
    } while (0)
#define PIRU_PROFILE_STOP(enabled, label) \
    do {                                  \
        if (enabled) piru::timing::stop(label); \
    } while (0)
