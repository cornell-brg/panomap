// SPDX-License-Identifier: MIT
/*
 * logging.hpp
 * Author: Elton PJ Shih <ps2229@cornell.edu>
 */

#pragma once

#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>

namespace panomap {

enum class LogLevel { TRACE = 0, DEBUG = 1, INFO = 2, WARN = 3, ERROR = 4, FATAL = 5 };

class Logger {
private:
  enum class ColorMode { AUTO, ALWAYS, NEVER };

  static Logger* instance_;
  static std::mutex mutex_;

  LogLevel current_level_;
  std::ostream* output_stream_;
  bool show_timestamp_;
  bool show_thread_id_;
  bool show_file_line_;
  ColorMode color_mode_;

  bool env_var_is_set(const char* name) const {
    const char* value = std::getenv(name);
    return value && *value != '\0';
  }

  bool env_var_is_zero(const char* name) const {
    const char* value = std::getenv(name);
    return value && std::string(value) == "0";
  }

  bool stream_is_terminal() const {
    if (output_stream_ == &std::cout) return ::isatty(STDOUT_FILENO);
    if (output_stream_ == &std::cerr) return ::isatty(STDERR_FILENO);
    return false;
  }

  bool should_enable_colors() const {
    if (env_var_is_set("NO_COLOR")) return false;
    if (env_var_is_set("FORCE_COLOR") || env_var_is_set("CLICOLOR_FORCE")) return true;
    if (env_var_is_zero("CLICOLOR")) return false;

    if (color_mode_ == ColorMode::ALWAYS) return true;
    if (color_mode_ == ColorMode::NEVER) return false;

    const char* term = std::getenv("TERM");
    if (!stream_is_terminal()) return false;
    if (!term || std::string(term) == "dumb") return false;
    return true;
  }

  Logger()
      : current_level_(LogLevel::INFO),
        output_stream_(&std::cerr),
        show_timestamp_(false),
        show_thread_id_(false),
        show_file_line_(true),
        color_mode_(ColorMode::AUTO) {}

  std::string get_timestamp() const {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    ss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return ss.str();
  }

  std::string get_thread_id() const {
    std::stringstream ss;
    ss << std::this_thread::get_id();
    return ss.str();
  }

  // ANSI color codes
  static constexpr const char* RESET = "\033[0m";
  static constexpr const char* BOLD = "\033[1m";
  static constexpr const char* RED = "\033[31m";
  static constexpr const char* GREEN = "\033[32m";
  static constexpr const char* YELLOW = "\033[33m";
  static constexpr const char* BLUE = "\033[34m";
  static constexpr const char* MAGENTA = "\033[35m";
  static constexpr const char* CYAN = "\033[36m";
  static constexpr const char* WHITE = "\033[37m";
  static constexpr const char* GRAY = "\033[90m";

  std::string get_color_code(LogLevel level) const {
    if (!should_enable_colors()) return "";

    switch (level) {
      case LogLevel::TRACE:
        return GRAY;
      case LogLevel::DEBUG:
        return CYAN;
      case LogLevel::INFO:
        return GREEN;
      case LogLevel::WARN:
        return YELLOW;
      case LogLevel::ERROR:
        return RED;
      case LogLevel::FATAL:
        return BOLD + std::string(RED);
      default:
        return "";
    }
  }

  std::string get_bold_color_code(LogLevel level) const {
    if (!should_enable_colors()) return "";

    switch (level) {
      case LogLevel::TRACE:
        return BOLD + std::string(GRAY);
      case LogLevel::DEBUG:
        return BOLD + std::string(CYAN);
      case LogLevel::INFO:
        return BOLD + std::string(GREEN);
      case LogLevel::WARN:
        return BOLD + std::string(YELLOW);
      case LogLevel::ERROR:
        return BOLD + std::string(RED);
      case LogLevel::FATAL:
        return BOLD + std::string(RED);
      default:
        return BOLD;
    }
  }

  std::string level_to_string(LogLevel level) const {
    std::string bold_color = get_bold_color_code(level);
    std::string reset = should_enable_colors() ? RESET : "";

    switch (level) {
      case LogLevel::TRACE:
        return bold_color + "TRAC" + reset;
      case LogLevel::DEBUG:
        return bold_color + "DBUG" + reset;
      case LogLevel::INFO:
        return bold_color + "INFO" + reset;
      case LogLevel::WARN:
        return bold_color + "WARN" + reset;
      case LogLevel::ERROR:
        return bold_color + "ERRO" + reset;
      case LogLevel::FATAL:
        return bold_color + "FATL" + reset;
      default:
        return bold_color + "UNKN" + reset;
    }
  }

  void write_log(LogLevel level, const std::string& message, const char* file = nullptr,
                 int line = 0) {
    if (level < current_level_) return;

    std::lock_guard<std::mutex> lock(mutex_);
    std::stringstream ss;
    const bool colors_enabled = should_enable_colors();

    if (show_timestamp_) ss << "[" << get_timestamp() << "] ";

    if (show_file_line_ && file && line > 0) {
      const char* filename = std::strrchr(file, '/');
      filename = filename ? filename + 1 : file;
      ss << "[" << level_to_string(level) << "::" << filename << ":" << line << "] ";
    } else {
      ss << "[" << level_to_string(level) << "] ";
    }

    if (show_thread_id_) ss << "[T:" << get_thread_id() << "] ";

    // Color the message text
    std::string message_color = get_color_code(level);
    std::string reset = colors_enabled ? RESET : "";
    ss << message_color << message << reset << std::endl;

    *output_stream_ << ss.str();
    output_stream_->flush();
  }

public:
  static Logger& get_instance() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!instance_) {
      instance_ = new Logger();
    }
    return *instance_;
  }

  void set_level(LogLevel level) {
    std::lock_guard<std::mutex> lock(mutex_);
    current_level_ = level;
  }

  LogLevel get_level() const { return current_level_; }

  void set_output_stream(std::ostream& stream) {
    std::lock_guard<std::mutex> lock(mutex_);
    output_stream_ = &stream;
  }

  void set_show_timestamp(bool show) {
    std::lock_guard<std::mutex> lock(mutex_);
    show_timestamp_ = show;
  }

  void set_show_thread_id(bool show) {
    std::lock_guard<std::mutex> lock(mutex_);
    show_thread_id_ = show;
  }

  void set_show_file_line(bool show) {
    std::lock_guard<std::mutex> lock(mutex_);
    show_file_line_ = show;
  }

  void set_enable_colors(bool enable) {
    std::lock_guard<std::mutex> lock(mutex_);
    color_mode_ = enable ? ColorMode::ALWAYS : ColorMode::NEVER;
  }

  void set_auto_colors() {
    std::lock_guard<std::mutex> lock(mutex_);
    color_mode_ = ColorMode::AUTO;
  }

  void trace(const std::string& message, const char* file = nullptr, int line = 0) {
    write_log(LogLevel::TRACE, message, file, line);
  }

  void debug(const std::string& message, const char* file = nullptr, int line = 0) {
    write_log(LogLevel::DEBUG, message, file, line);
  }

  void info(const std::string& message, const char* file = nullptr, int line = 0) {
    write_log(LogLevel::INFO, message, file, line);
  }

  void warn(const std::string& message, const char* file = nullptr, int line = 0) {
    write_log(LogLevel::WARN, message, file, line);
  }

  void error(const std::string& message, const char* file = nullptr, int line = 0) {
    write_log(LogLevel::ERROR, message, file, line);
  }

  void fatal(const std::string& message, const char* file = nullptr, int line = 0) {
    write_log(LogLevel::FATAL, message, file, line);
    std::exit(1);
  }
};

// Static member initialization (inline to make it header-only)
inline Logger* Logger::instance_ = nullptr;
inline std::mutex Logger::mutex_;

// Global logger instance
inline Logger& logger = Logger::get_instance();

// Convenience macros for logging
#define LOG_TRACE(msg) panomap::logger.trace(msg, __FILE__, __LINE__)
#define LOG_DEBUG(msg) panomap::logger.debug(msg, __FILE__, __LINE__)
#define LOG_INFO(msg) panomap::logger.info(msg, __FILE__, __LINE__)
#define LOG_WARN(msg) panomap::logger.warn(msg, __FILE__, __LINE__)
#define LOG_ERROR(msg) panomap::logger.error(msg, __FILE__, __LINE__)
#define LOG_FATAL(msg) panomap::logger.fatal(msg, __FILE__, __LINE__)

// Shorter macros for common use
#define PANOMAP_TRACE(msg) LOG_TRACE(msg)
#define PANOMAP_DEBUG(msg) LOG_DEBUG(msg)
#define PANOMAP_INFO(msg) LOG_INFO(msg)
#define PANOMAP_WARN(msg) LOG_WARN(msg)
#define PANOMAP_ERROR(msg) LOG_ERROR(msg)
#define PANOMAP_FATAL(msg) LOG_FATAL(msg)

}  // namespace panomap
