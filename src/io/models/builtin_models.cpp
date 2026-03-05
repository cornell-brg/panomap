#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "io/models/builtin_model_common.hpp"
#include "io/models/builtin_model_r10.hpp"
#include "io/models/builtin_model_r9.hpp"
#include "io/models/model_factory.hpp"
#include "util/logging.hpp"

namespace piru::io {

namespace {

using Table = std::unordered_map<std::string, double>;
using TablePtr = std::shared_ptr<const Table>;

class StaticKmerModel : public KmerModel {
public:
    StaticKmerModel(std::string model_name, int kmer_size, TablePtr table)
        : name_(std::move(model_name)), k_(kmer_size), table_(std::move(table)) {}

    std::string name() const override { return name_; }
    int k() const override { return k_; }

    bool lookup(const std::string& kmer, double& mean) const override {
        const auto& table = *table_;
        auto it = table.find(kmer);
        if (it == table.end()) return false;
        mean = it->second;
        return true;
    }

private:
    std::string name_;
    int k_;
    TablePtr table_;
};

Table make_table(const generated::ModelEntry* entries, std::size_t count) {
    Table table;
    table.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const auto& entry = entries[i];
        table.emplace(entry.kmer, entry.mean);
    }
    return table;
}

TablePtr r9_table() {
    static const Table table = make_table(generated::kR9Entries, generated::kR9Count);
    static const TablePtr ptr(&table, [](const Table*) {});
    return ptr;
}

TablePtr r10_table() {
    static const Table table = make_table(generated::kR10Entries, generated::kR10Count);
    static const TablePtr ptr(&table, [](const Table*) {});
    return ptr;
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

bool is_r9_name(const std::string& name) {
    const std::string lower = to_lower(name);
    return lower == "r9.4" || lower == "r9_4" || lower == "r9.4_450bps" || lower == "r9_4_450bps" ||
           lower == "r9";
}

bool is_r10_name(const std::string& name) {
    const std::string lower = to_lower(name);
    return lower == "r10.4" || lower == "r10_4" || lower == "r10.4_400bps" ||
           lower == "r10_4_400bps" || lower == "r10";
}

ModelPtr make_r9() {
    return std::make_unique<StaticKmerModel>("r9.4_450bps", generated::kR9K, r9_table());
}

ModelPtr make_r10() {
    return std::make_unique<StaticKmerModel>("r10.4_400bps", generated::kR10K, r10_table());
}

}  // namespace

ModelPtr load_builtin_model(const std::string& name) {
    if (is_r9_name(name)) {
        return make_r9();
    }
    if (is_r10_name(name)) {
        return make_r10();
    }
    LOG_ERROR("Unknown built-in model: " + name);
    return nullptr;
}

ModelPtr load_model_from_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        LOG_ERROR("Failed to open model file: " + path);
        return nullptr;
    }

    std::string line;
    // Consume header if present.
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string first;
        iss >> first;
        if (first == "kmer" || first == "#kmer") {
            // header; skip and read next line
            continue;
        }
        // We already read the first token of the first data line.
        std::vector<std::string> tokens;
        tokens.push_back(first);
        std::string tok;
        while (iss >> tok) tokens.push_back(tok);

        const std::size_t num_cols = tokens.size();
        if (num_cols < 2) {
            LOG_ERROR("Could not parse model line in " + path + ": " + line);
            return nullptr;
        }

        const int kmer_size = static_cast<int>(tokens[0].size());
        Table table;
        table.reserve(4096);  // heuristic

        auto add_entry = [&](const std::vector<std::string>& cols) {
            if (cols.size() < 2) return;
            const std::string& kmer = cols[0];
            double mean = 0.0;
            try {
                mean = std::stod(cols[1]);
            } catch (const std::exception&) {
                return;
            }
            table.emplace(kmer, mean);
        };

        add_entry(tokens);
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream data(line);
            std::vector<std::string> cols;
            std::string c;
            while (data >> c) cols.push_back(c);
            add_entry(cols);
        }

        if (table.empty()) {
            LOG_ERROR("No entries parsed from model file: " + path);
            return nullptr;
        }

        const std::string model_name = std::filesystem::path(path).stem().string();
        return std::make_unique<StaticKmerModel>(model_name, kmer_size,
                                                 std::make_shared<Table>(std::move(table)));
    }

    LOG_ERROR("Model file empty: " + path);
    return nullptr;
}

}  // namespace piru::io
