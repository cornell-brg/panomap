#include "io/graphs/vg_loader.hpp"

#include <exception>
#include <unordered_set>
#include <string>

#ifdef PIRU_HAS_LIBVGIO
#include <vg/io/basic_stream.hpp>
#include <vg/vg.pb.h>
#endif

#include "util/logging.hpp"

namespace piru::io {

VgLoader::VgLoader(std::string path) : path_(std::move(path)) {}

bool VgLoader::load(ImportedGraph& graph) {
    graph.clear();

#ifdef PIRU_HAS_LIBVGIO
    try {
        vg::Graph vg_graph = vg::io::inputStream(path_);
        for (const auto& node : vg_graph.node()) {
            ImportedGraphNode gnode;
            gnode.id = std::to_string(node.id());
            gnode.sequence = node.sequence();
            graph.add_node(std::move(gnode));
        }
        for (const auto& edge : vg_graph.edge()) {
            ImportedGraphEdge gedge;
            gedge.from = std::to_string(edge.from());
            gedge.to = std::to_string(edge.to());
            gedge.from_reverse = edge.from_start();
            // to_end=true means arriving at end of to node (reverse orientation in GFA).
            gedge.to_reverse = edge.to_end();
            gedge.overlap = std::to_string(edge.overlap());  // 0 when unspecified.
            if (edge.overlap() != 0) {
                gedge.overlap_bases = static_cast<std::size_t>(edge.overlap());
            }
            graph.add_edge(std::move(gedge));
        }
        std::unordered_map<std::string, std::size_t> path_indices;
        for (const auto& path : vg_graph.path()) {
            ImportedPath imported_path;
            imported_path.name = path.name();
            for (const auto& mapping : path.mapping()) {
                if (!mapping.has_position()) continue;
                const auto& pos = mapping.position();
                if (pos.node_id() == 0) continue;
                ImportedPathStep step;
                step.segment_id = std::to_string(pos.node_id());
                step.is_reverse = pos.is_reverse();
                imported_path.steps.push_back(std::move(step));
            }
            const auto it = path_indices.find(imported_path.name);
            if (it == path_indices.end()) {
                path_indices[imported_path.name] = graph.paths.size();
                graph.add_path(std::move(imported_path));
            } else {
                auto& existing = graph.paths[it->second];
                existing.steps.insert(existing.steps.end(), imported_path.steps.begin(),
                                      imported_path.steps.end());
            }
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to parse vg file '" + path_ + "': " + e.what());
        return false;
    } catch (...) {
        LOG_ERROR("Unknown error parsing vg file '" + path_ + "'");
        return false;
    }

    if (graph.nodes.empty() && graph.edges.empty()) {
        LOG_WARN("Parsed vg file '" + path_ + "' but found no nodes/edges");
    }

    return true;
#else
    if (!warned_) {
        LOG_ERROR("libvgio not linked; vg loading unavailable for '" + path_ + "'");
        warned_ = true;
    }
    return false;
#endif
}

std::string VgLoader::get_format_name() const { return "vg"; }

}  // namespace piru::io
