#include "io/graphs/gfa_exporter.hpp"

#include <fstream>
#include <iostream>
#include <variant>

#include "io/graphs/graph.hpp"
#include "signal/signal_types.hpp"

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

void GfaExporter::dumpAlnGraph(
    const index::AlnGraph& graph, const std::string& path,
    AlnGraphDumpMode mode, const void* data,
    const std::map<std::string, std::string>& header_tags) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Error: could not open file " << path << std::endl;
        return;
    }

    out << "H\tVN:Z:1.0";
    for (const auto& tag : header_tags) {
        out << "\t" << tag.first << ":" << tag.second;
    }
    out << std::endl;

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
        const auto* signal_data = static_cast<const std::vector<signal::AlignmentQuantizedSignal>*>(data);
        if (!signal_data) {
            return;
        }
        for (size_t i = 0; i < graph.nodeCount(); ++i) {
            const auto& node = graph.node(i);
            out << "S\t" << node.id << "\t";
            
            const auto& aln_sig = (*signal_data)[i];
            if (std::holds_alternative<std::vector<int16_t>>(aln_sig.data)) {
                const auto& signals = std::get<std::vector<int16_t>>(aln_sig.data);
                for (size_t j = 0; j < signals.size(); ++j) {
                    out << signals[j];
                    if (j < signals.size() - 1) {
                        out << ",";
                    }
                }
                out << "\tLN:i:" << signals.size();
            } else if (std::holds_alternative<std::vector<int8_t>>(aln_sig.data)) {
                const auto& signals = std::get<std::vector<int8_t>>(aln_sig.data);
                for (size_t j = 0; j < signals.size(); ++j) {
                    out << static_cast<int>(signals[j]);
                    if (j < signals.size() - 1) {
                        out << ",";
                    }
                }
                out << "\tLN:i:" << signals.size();
            } else if (std::holds_alternative<std::vector<float>>(aln_sig.data)) {
                const auto& signals = std::get<std::vector<float>>(aln_sig.data);
                for (size_t j = 0; j < signals.size(); ++j) {
                    out << signals[j];
                    if (j < signals.size() - 1) {
                        out << ",";
                    }
                }
                out << "\tLN:i:" << signals.size();
            } else {
                out << "\tLN:i:0";
            }
            
            out << "\tst:Z:aln_quant";
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
            << edge.overlap_bases << "M";
        const auto& from_node = graph.node(edge.from);
        const auto& to_node = graph.node(edge.to);
        if (!from_node.original_id.empty()) {
            out << "\tof:Z:" << from_node.original_id << (from_node.is_reverse ? "-" : "+");
        }
        if (!to_node.original_id.empty()) {
            out << "\tot:Z:" << to_node.original_id << (to_node.is_reverse ? "-" : "+");
        }
        out << std::endl;
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
