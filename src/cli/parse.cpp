#include "cli/parse.hpp"

#include <sstream>

namespace piru::cli {

namespace {

std::string build_optstring(const std::vector<Option>& opts) {
  std::string optstring;
  for (const auto& o : opts) {
    if (o.short_opt == 0) continue;
    optstring.push_back(o.short_opt);
    if (o.requires_arg) optstring.push_back(':');
  }
  return optstring;
}

std::vector<::option> build_longopts(const std::vector<Option>& opts) {
  std::vector<::option> longs;
  longs.reserve(opts.size() + 1);
  for (const auto& o : opts) {
    if (o.long_opt.empty()) continue;
    // For long-only options (short_opt == 0), getopt_long returns 0,
    // so we use a unique non-zero value to distinguish them.
    // We'll use negative indices starting from -1 for long-only options.
    int val = o.short_opt;
    if (val == 0) {
      // Assign unique negative value based on position in vector
      val = -1 - static_cast<int>(&o - opts.data());
    }
    longs.push_back(
        {o.long_opt.c_str(), o.requires_arg ? required_argument : no_argument, nullptr, val});
  }
  longs.push_back({nullptr, 0, nullptr, 0});
  return longs;
}

}  // namespace

bool parse_args(const std::vector<std::string>& args, const ParseConfig& config, Parsed& out) {
  const std::string optstring = build_optstring(config.options);
  const auto longopts = build_longopts(config.options);

  // Build argv
  std::vector<std::string> argv_storage;
  argv_storage.reserve(args.size() + 1);
  argv_storage.push_back("cmd");
  argv_storage.insert(argv_storage.end(), args.begin(), args.end());

  std::vector<char*> argv_ptrs;
  argv_ptrs.reserve(argv_storage.size());
  for (auto& s : argv_storage) argv_ptrs.push_back(s.data());
  int argc = static_cast<int>(argv_ptrs.size());

  opterr = 0;
  optind = 0;
  int long_index = 0;
  int c = -1;
  while ((c = ::getopt_long(argc, argv_ptrs.data(), optstring.c_str(), longopts.data(),
                            &long_index)) != -1) {
    if (c == '?') {
      if (config.on_error) config.on_error("Unknown or invalid option");
      return false;
    }

    // Find the matching option
    std::vector<Option>::const_iterator it;
    if (c < 0) {
      // Negative value means this is a long-only option
      // Decode the index: c = -1 - index
      int index = -1 - c;
      if (index >= 0 && index < static_cast<int>(config.options.size())) {
        it = config.options.begin() + index;
      } else {
        if (config.on_error) config.on_error("Internal error: invalid option index");
        return false;
      }
    } else {
      // Regular option with short name
      it = std::find_if(config.options.begin(), config.options.end(),
                        [&](const Option& o) { return o.short_opt == c; });
    }

    if (it == config.options.end()) {
      if (config.on_error) config.on_error("Unknown option");
      return false;
    }
    const std::string key = it->long_opt.empty() ? std::string(1, it->short_opt) : it->long_opt;
    if (it->requires_arg) {
      out.values[key] = optarg ? std::string(optarg) : "";
    } else {
      out.values[key] = "true";
    }
  }

  for (int i = optind; i < argc; ++i) {
    out.positionals.emplace_back(argv_ptrs[i]);
  }
  return true;
}

void print_help(const ParseConfig& config, std::ostream& os) {
  os << config.usage << "\n";
  if (!config.positional_help.empty()) {
    os << "\nRequired arguments:\n";
    for (const auto& line : config.positional_help) {
      os << "  " << line << "\n";
    }
  }
  os << "\nOptions:\n";
  for (const auto& o : config.options) {
    // Skip hidden options (no help text, no short opt)
    if (o.help.empty() && !o.short_opt && !o.long_opt.empty()) continue;
    std::ostringstream line;
    line << "  ";
    if (o.short_opt) line << "-" << o.short_opt;
    if (o.short_opt && !o.long_opt.empty()) line << ", ";
    if (!o.long_opt.empty()) line << "--" << o.long_opt;
    os << line.str();
    if (!o.help.empty()) {
      if (line.str().size() < 24) {
        os << std::string(24 - line.str().size(), ' ');
      } else {
        os << " ";
      }
      os << o.help;
    }
    os << "\n";
  }
}

}  // namespace piru::cli
