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

namespace piru {

enum class LogLevel {
  TRACE = 0,
  DEBUG = 1,
  INFO  = 2,
  WARN  = 3,
  ERROR = 4,
  FATAL = 5
};

class Logger {
private:
  static Logger *instance_;
  static std::mutex mutex_;

  LogLevel current_level_;
  std::ostream *output_stream_;
  bool show_timestamp_;
  bool show_thread_id_;
  bool show_file_line_;
  bool enable_colors_;

  bool is_terminal_color_supported() const {
    // Check if output is to a terminal and supports colors
    if (output_stream_ == &std::cout || output_stream_ == &std::cerr) {
      const char *term = std::getenv("TERM");
      if (term && std::string(term) != "dumb") {
        return true;
      }
    }
    return false;
  }

  Logger()
      : current_level_(LogLevel::INFO), output_stream_(&std::cerr),
        show_timestamp_(false), show_thread_id_(false), show_file_line_(true),
        enable_colors_(true) {
    // Auto-detect color support
    enable_colors_ = is_terminal_color_supported();
  }

  std::string get_timestamp() const {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
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
  static constexpr const char *RESET   = "\033[0m";
  static constexpr const char *BOLD    = "\033[1m";
  static constexpr const char *RED     = "\033[31m";
  static constexpr const char *GREEN   = "\033[32m";
  static constexpr const char *YELLOW  = "\033[33m";
  static constexpr const char *BLUE    = "\033[34m";
  static constexpr const char *MAGENTA = "\033[35m";
  static constexpr const char *CYAN    = "\033[36m";
  static constexpr const char *WHITE   = "\033[37m";
  static constexpr const char *GRAY    = "\033[90m";

  std::string get_color_code(LogLevel level) const {
    if (!enable_colors_) return "";

    switch (level) {
    case LogLevel::TRACE: return GRAY;
    case LogLevel::DEBUG: return CYAN;
    case LogLevel::INFO:  return GREEN;
    case LogLevel::WARN:  return YELLOW;
    case LogLevel::ERROR: return RED;
    case LogLevel::FATAL: return BOLD + std::string(RED);
    default: return "";
    }
  }

  std::string get_bold_color_code(LogLevel level) const {
    if (!enable_colors_) return "";

    switch (level) {
    case LogLevel::TRACE: return BOLD + std::string(GRAY);
    case LogLevel::DEBUG: return BOLD + std::string(CYAN);
    case LogLevel::INFO:  return BOLD + std::string(GREEN);
    case LogLevel::WARN:  return BOLD + std::string(YELLOW);
    case LogLevel::ERROR: return BOLD + std::string(RED);
    case LogLevel::FATAL: return BOLD + std::string(RED);
    default: return BOLD;
    }
  }

  std::string level_to_string(LogLevel level) const {
    std::string bold_color = get_bold_color_code(level);
    std::string reset = enable_colors_ ? RESET : "";

    switch (level) {
    case LogLevel::TRACE: return bold_color + "TRAC" + reset;
    case LogLevel::DEBUG: return bold_color + "DBUG" + reset;
    case LogLevel::INFO:  return bold_color + "INFO" + reset;
    case LogLevel::WARN:  return bold_color + "WARN" + reset;
    case LogLevel::ERROR: return bold_color + "ERRO" + reset;
    case LogLevel::FATAL: return bold_color + "FATL" + reset;
    default: return bold_color + "UNKN" + reset;
    }
  }

  void write_log(LogLevel level, const std::string &message,
                 const char *file = nullptr, int line = 0) {
    if (level < current_level_) return;

    std::lock_guard<std::mutex> lock(mutex_);
    std::stringstream ss;

    if (show_timestamp_) ss << "[" << get_timestamp() << "] ";

    if (show_file_line_ && file && line > 0) {
      const char *filename = std::strrchr(file, '/');
      filename = filename ? filename + 1 : file;
      ss << "[" << level_to_string(level) << "::" << filename << ":" << line << "] ";
    } else {
      ss << "[" << level_to_string(level) << "] ";
    }

    if (show_thread_id_) ss << "[T:" << get_thread_id() << "] ";

    // Color the message text
    std::string message_color = get_color_code(level);
    std::string reset = enable_colors_ ? RESET : "";
    ss << message_color << message << reset << std::endl;

    *output_stream_ << ss.str();
    output_stream_->flush();
  }

public:
  static Logger &get_instance() {
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

  void set_output_stream(std::ostream &stream) {
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
    enable_colors_ = enable;
  }

  void trace(const std::string &message, const char *file = nullptr, int line = 0) {
    write_log(LogLevel::TRACE, message, file, line);
  }

  void debug(const std::string &message, const char *file = nullptr, int line = 0) {
    write_log(LogLevel::DEBUG, message, file, line);
  }

  void info(const std::string &message, const char *file = nullptr, int line = 0) {
    write_log(LogLevel::INFO, message, file, line);
  }

  void warn(const std::string &message, const char *file = nullptr, int line = 0) {
    write_log(LogLevel::WARN, message, file, line);
  }

  void error(const std::string &message, const char *file = nullptr, int line = 0) {
    write_log(LogLevel::ERROR, message, file, line);
  }

  void fatal(const std::string &message, const char *file = nullptr, int line = 0) {
    write_log(LogLevel::FATAL, message, file, line);
    std::exit(1);
  }
};

// Static member initialization (inline to make it header-only)
inline Logger *Logger::instance_ = nullptr;
inline std::mutex Logger::mutex_;

// Global logger instance
inline Logger &logger = Logger::get_instance();

// Convenience macros for logging
#define LOG_TRACE(msg) piru::logger.trace(msg, __FILE__, __LINE__)
#define LOG_DEBUG(msg) piru::logger.debug(msg, __FILE__, __LINE__)
#define LOG_INFO(msg)  piru::logger.info(msg, __FILE__, __LINE__)
#define LOG_WARN(msg)  piru::logger.warn(msg, __FILE__, __LINE__)
#define LOG_ERROR(msg) piru::logger.error(msg, __FILE__, __LINE__)
#define LOG_FATAL(msg) piru::logger.fatal(msg, __FILE__, __LINE__)

// Shorter macros for common use
#define PIRU_TRACE(msg) LOG_TRACE(msg)
#define PIRU_DEBUG(msg) LOG_DEBUG(msg)
#define PIRU_INFO(msg)  LOG_INFO(msg)
#define PIRU_WARN(msg)  LOG_WARN(msg)
#define PIRU_ERROR(msg) LOG_ERROR(msg)
#define PIRU_FATAL(msg) LOG_FATAL(msg)

}  // namespace piru
