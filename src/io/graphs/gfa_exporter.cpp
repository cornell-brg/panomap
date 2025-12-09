#include "io/graphs/gfa_exporter.hpp"

#include <fstream>
#include <iostream>

#include "io/graphs/graph.hpp"

namespace piru {

void GfaExporter::dumpImportedGraph(const io::ImportedGraph& graph, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Error: could not open file " << path << std::endl;
        return;
    }

    out << "H\tVN:Z:1.0" << std::endl;

    for (const auto& node : graph.nodes) {
        out << "S\t" << node.id << "\t" << node.sequence << "\tLN:i:" << node.sequence.length() << std::endl;
    }

    for (const auto& edge : graph.edges) {
        out << "L\t" << edge.from << "\t" << (edge.from_reverse ? "-" : "+") << "\t"
            << edge.to << "\t" << (edge.to_reverse ? "-" : "+") << "\t"
            << edge.overlap << "M" << std::endl;
    }

    for (const auto& path : graph.paths) {
        out << "P\t" << path.name << "\t";
        for (size_t i = 0; i < path.steps.size(); ++i) {
            out << path.steps[i].segment_id << (path.steps[i].is_reverse ? "-" : "+");
            if (i < path.steps.size() - 1) {
                out << ",";
            }
        }
        out << "\t";
        if (!path.overlaps.empty()) {
            for (size_t i = 0; i < path.overlaps.size(); ++i) {
                out << path.overlaps[i];
                if (i < path.overlaps.size() - 1) {
                    out << ",";
                }
            }
        } else {
            out << "*";
        }
        out << std::endl;
    }
}

void GfaExporter::dumpAlnGraph(const index::AlnGraph& graph, const std::string& path,
                               AlnGraphDumpMode mode, const void* data) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Error: could not open file " << path << std::endl;
        return;
    }

    out << "H\tVN:Z:1.0" << std::endl;

    if (mode == AlnGraphDumpMode::Bases) {
        for (size_t i = 0; i < graph.nodeCount(); ++i) {
            const auto& node = graph.node(i);
            out << "S\t" << node.id << "\t" << node.sequence << "\tLN:i:" << node.sequence.length();
            if (!node.original_id.empty()) {
                out << "\toi:Z:" << node.original_id << (node.is_reverse ? "-" : "+");
            }
            if (node.chain_id) {
                out << "\tci:i:" << *node.chain_id;
            }
            if (node.linear_position) {
                out << "\tlc:i:" << *node.linear_position;
            }
            out << std::endl;
        }
    } else if (mode == AlnGraphDumpMode::RawSignal) {
        const auto* signal_data = static_cast<const std::vector<std::vector<float>>*>(data);
        if (!signal_data) {
            return;
        }
        for (size_t i = 0; i < graph.nodeCount(); ++i) {
            const auto& node = graph.node(i);
            out << "S\t" << node.id << "\t";
            const auto& signals = (*signal_data)[i];
            for (size_t j = 0; j < signals.size(); ++j) {
                out << signals[j];
                if (j < signals.size() - 1) {
                    out << ",";
                }
            }
            out << "\tLN:i:" << signals.size() << "\tst:Z:raw_signal";
            if (!node.original_id.empty()) {
                out << "\toi:Z:" << node.original_id << (node.is_reverse ? "-" : "+");
            }
            if (node.chain_id) {
                out << "\tci:i:" << *node.chain_id;
            }
            if (node.linear_position) {
                out << "\tlc:i:" << *node.linear_position;
            }
            out << std::endl;
        }
    } else if (mode == AlnGraphDumpMode::FuzzyQuantized) {
        const auto* signal_data = static_cast<const std::vector<std::vector<int16_t>>*>(data);
        if (!signal_data) {
            return;
        }
        for (size_t i = 0; i < graph.nodeCount(); ++i) {
            const auto& node = graph.node(i);
            out << "S\t" << node.id << "\t";
            const auto& signals = (*signal_data)[i];
            for (size_t j = 0; j < signals.size(); ++j) {
                out << signals[j];
                if (j < signals.size() - 1) {
                    out << ",";
                }
            }
            out << "\tLN:i:" << signals.size() << "\tst:Z:fuzzy_quant";
            if (!node.original_id.empty()) {
                out << "\toi:Z:" << node.original_id << (node.is_reverse ? "-" : "+");
            }
            if (node.chain_id) {
                out << "\tci:i:" << *node.chain_id;
            }
            if (node.linear_position) {
                out << "\tlc:i:" << *node.linear_position;
            }
            out << std::endl;
        }
    } else if (mode == AlnGraphDumpMode::AlnQuantized) {
        const auto* signal_data = static_cast<const std::vector<std::vector<int16_t>>*>(data);
        if (!signal_data) {
            return;
        }
        for (size_t i = 0; i < graph.nodeCount(); ++i) {
            const auto& node = graph.node(i);
            out << "S\t" << node.id << "\t";
            const auto& signals = (*signal_data)[i];
            for (size_t j = 0; j < signals.size(); ++j) {
                out << signals[j];
                if (j < signals.size() - 1) {
                    out << ",";
                }
            }
            out << "\tLN:i:" << signals.size() << "\tst:Z:aln_quant";
            if (!node.original_id.empty()) {
                out << "\toi:Z:" << node.original_id << (node.is_reverse ? "-" : "+");
            }
            if (node.chain_id) {
                out << "\tci:i:" << *node.chain_id;
            }
            if (node.linear_position) {
                out << "\tlc:i:" << *node.linear_position;
            }
            out << std::endl;
        }
    }

    for (const auto& edge : graph.edges()) {
        out << "L\t" << edge.from << "\t+\t" << edge.to << "\t+\t"
            << edge.overlap_bases << "M" << std::endl;
    }
    
    for (const auto& path : graph.paths()) {
        out << "P\t" << path.name << "\t";
        for (size_t i = 0; i < path.steps.size(); ++i) {
            const auto& step = path.steps[i];
            out << step.node_id << (step.is_reverse ? "-" : "+");
            if (i < path.steps.size() - 1) {
                out << ",";
            }
        }
        out << "\t";
        if (!path.overlaps.empty()) {
            for (size_t i = 0; i < path.overlaps.size(); ++i) {
                out << path.overlaps[i] << "M";
                if (i < path.overlaps.size() - 1) {
                    out << ",";
                }
            }
        } else {
            out << "*";
        }
        out << std::endl;
    }
}

}  // namespace piru
