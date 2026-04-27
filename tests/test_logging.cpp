// SPDX-License-Identifier: MIT
#include <doctest/doctest.h>
#include <sstream>
#include <string>

#include "core/util/logging.hpp"

namespace {

bool contains_ansi_escape(const std::string& s) { return s.find("\033[") != std::string::npos; }

struct EnvVarGuard {
  explicit EnvVarGuard(const char* name) : name(name) {
    const char* current = std::getenv(name);
    if (current) {
      had_value = true;
      value = current;
    }
  }

  ~EnvVarGuard() {
    if (had_value) {
      ::setenv(name, value.c_str(), 1);
    } else {
      ::unsetenv(name);
    }
  }

  const char* name;
  bool had_value{false};
  std::string value;
};

}  // namespace

TEST_CASE("logger auto mode disables colors for non-terminal streams") {
  std::ostringstream capture;

  piru::logger.set_output_stream(capture);
  piru::logger.set_auto_colors();
  piru::logger.set_show_file_line(false);
  piru::logger.info("plain log");

  CHECK_FALSE(contains_ansi_escape(capture.str()));

  piru::logger.set_output_stream(std::cerr);
  piru::logger.set_auto_colors();
  piru::logger.set_show_file_line(true);
}

TEST_CASE("logger explicit color toggle overrides auto mode") {
  std::ostringstream capture;
  EnvVarGuard no_color("NO_COLOR");
  EnvVarGuard force_color("FORCE_COLOR");
  EnvVarGuard clicolor("CLICOLOR");
  EnvVarGuard clicolor_force("CLICOLOR_FORCE");

  ::unsetenv("NO_COLOR");
  ::unsetenv("FORCE_COLOR");
  ::unsetenv("CLICOLOR_FORCE");
  ::unsetenv("CLICOLOR");

  piru::logger.set_output_stream(capture);
  piru::logger.set_show_file_line(false);

  piru::logger.set_enable_colors(true);
  piru::logger.warn("forced color");
  CHECK(contains_ansi_escape(capture.str()));

  capture.str("");
  capture.clear();

  piru::logger.set_enable_colors(false);
  piru::logger.warn("forced plain");
  CHECK_FALSE(contains_ansi_escape(capture.str()));

  piru::logger.set_output_stream(std::cerr);
  piru::logger.set_auto_colors();
  piru::logger.set_show_file_line(true);
}
