// SPDX-License-Identifier: MIT
#include "index/pseudo_linearize.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <queue>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace piru::index {
namespace {

struct Frame {
    std::size_t node;
    std::size_t next_child{0};
};

bool hasSelfLoop(const AlnGraph& graph, std::size_t node) {
    const auto& outs = graph.outgoing(node);
    return std::find(outs.begin(), outs.end(), node) != outs.end();
}

std::optional<std::size_t> findSuperbubble(const AlnGraph& graph, std::size_t start,
                                           const std::vector<bool>& ignorable_tip);
std::vector<std::size_t> traverseBubbleInterior(const AlnGraph& graph, std::size_t start,
                                                std::size_t end,
                                                const std::vector<bool>& ignorable_tip);

}  // namespace

SccResult computeScc(const AlnGraph& graph) {
    const std::size_t n = graph.nodeCount();
    SccResult res;
    res.component.assign(n, 0);
    res.in_scc.assign(n, false);

    std::vector<int> index(n, -1);
    std::vector<int> lowlink(n, 0);
    std::vector<bool> on_stack(n, false);
    std::vector<std::size_t> stack;
    stack.reserve(n);

    int index_counter = 0;
    std::size_t comp_counter = 0;

    for (std::size_t start = 0; start < n; ++start) {
        if (index[start] != -1) continue;

        std::vector<Frame> call_stack;
        call_stack.push_back({start, 0});

        while (!call_stack.empty()) {
            Frame& frame = call_stack.back();
            const std::size_t v = frame.node;

            if (index[v] == -1) {
                index[v] = lowlink[v] = index_counter++;
                stack.push_back(v);
                on_stack[v] = true;
            }

            const auto& outs = graph.outgoing(v);
            bool descended = false;
            while (frame.next_child < outs.size()) {
                const std::size_t w = outs[frame.next_child++];
                if (index[w] == -1) {
                    call_stack.push_back({w, 0});
                    descended = true;
                    break;
                }
                if (on_stack[w]) {
                    lowlink[v] = std::min(lowlink[v], index[w]);
                }
            }
            if (descended) continue;

            // Finished exploring v.
            call_stack.pop_back();

            // If v is a root node, pop the stack and generate an SCC.
            if (lowlink[v] == index[v]) {
                std::vector<std::size_t> component_nodes;
                while (!stack.empty()) {
                    const std::size_t w = stack.back();
                    stack.pop_back();
                    on_stack[w] = false;
                    res.component[w] = comp_counter;
                    component_nodes.push_back(w);
                    if (w == v) break;
                }

                const bool nontrivial = component_nodes.size() > 1 ||
                                        (component_nodes.size() == 1 && hasSelfLoop(graph, v));
                if (nontrivial) {
                    for (const auto w : component_nodes) {
                        res.in_scc[w] = true;
                    }
                }
                ++comp_counter;
            }

            // Propagate lowlink to parent.
            if (!call_stack.empty()) {
                const std::size_t parent = call_stack.back().node;
                lowlink[parent] = std::min(lowlink[parent], lowlink[v]);
            }
        }
    }

    res.components = comp_counter;

    // Reverse component numbering so edges go forward (higher component for successors).
    // Tarjan assigns components in reverse topological order; we need forward order
    // for tip folding and superbubble chaining.
    for (std::size_t v = 0; v < n; ++v) {
        res.component[v] = (comp_counter - 1) - res.component[v];
    }

    // Validation: edges should go from <= to >= component numbers (topological order)
    // This is an invariant check, not user-facing validation
    for (std::size_t v = 0; v < n; ++v) {
        for (const auto u : graph.outgoing(v)) {
            // Edge v → u should have component[v] <= component[u]
            if (res.component[u] < res.component[v]) {
                // This should never happen if Tarjan's algorithm is correct
                throw std::logic_error(
                    "SCC component ordering violation: edge " + std::to_string(v) +
                    " → " + std::to_string(u) + " goes from component " +
                    std::to_string(res.component[v]) + " to " +
                    std::to_string(res.component[u]));
            }
        }
    }

    return res;
}

std::size_t UnionFind::find(std::size_t x) {
    // Path compression: make every node point directly to root
    if (parent[x] != x) {
        parent[x] = find(parent[x]);
    }
    return parent[x];
}

void UnionFind::unite(std::size_t x, std::size_t y) {
    // Union by rank: attach smaller tree under root of larger tree
    x = find(x);
    y = find(y);
    if (x == y) return;

    if (rank[x] < rank[y]) {
        parent[x] = y;
    } else if (rank[x] > rank[y]) {
        parent[y] = x;
    } else {
        parent[y] = x;
        rank[x]++;
    }
}

TipFoldingResult chainTips(const AlnGraph& graph, const SccResult& scc) {
    const std::size_t n = graph.nodeCount();
    TipFoldingResult result(n);

    if (n == 0) return result;

    // Track which components are tips.
    std::vector<bool> is_forward_tip(scc.components, true);
    std::vector<bool> is_backward_tip(scc.components, true);

    // Order nodes by component number for efficient component-wise iteration.
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return scc.component[a] < scc.component[b];
    });

    // Forward tip detection (traverse components high → low, leaf-first).
    for (std::size_t ind = n - 1; ind < n; --ind) {  // Unsigned wrap-around idiom
        const std::size_t i = order[ind];
        const std::size_t ci = scc.component[i];

        if (!is_forward_tip[ci]) continue;

        for (const auto succ : graph.outgoing(i)) {
            const std::size_t cn = scc.component[succ];

            // If has edge to same component, not a tip
            if (cn == ci) {
                is_forward_tip[ci] = false;
                break;
            }
            // If has edge to non-tip component, not a tip
            if (!is_forward_tip[cn]) {
                is_forward_tip[ci] = false;
                break;
            }
        }
    }

    // Chain forward tips
    for (std::size_t ind = n - 1; ind < n; --ind) {
        const std::size_t i = order[ind];
        const std::size_t ci = scc.component[i];

        if (!is_forward_tip[ci]) continue;

        for (const auto succ : graph.outgoing(i)) {
            const std::size_t cn = scc.component[succ];
            if (!is_forward_tip[cn]) continue;
            if (i == succ) continue;  // Skip self-loops

            result.uf.unite(i, succ);
        }
    }

    // Backward tip detection (traverse components low → high).
    for (std::size_t ind = 0; ind < n; ++ind) {
        const std::size_t i = order[ind];
        const std::size_t ci = scc.component[i];

        if (!is_backward_tip[ci]) continue;

        for (const auto pred : graph.incoming(i)) {
            const std::size_t cn = scc.component[pred];

            // If has edge from same component, not a tip
            if (cn == ci) {
                is_backward_tip[ci] = false;
                break;
            }
            // If has edge from non-tip component, not a tip
            if (!is_backward_tip[cn]) {
                is_backward_tip[ci] = false;
                break;
            }
        }
    }

    // Chain backward tips
    for (std::size_t ind = 0; ind < n; ++ind) {
        const std::size_t i = order[ind];
        const std::size_t ci = scc.component[i];

        if (!is_backward_tip[ci]) continue;

        for (const auto pred : graph.incoming(i)) {
            const std::size_t cn = scc.component[pred];
            if (!is_backward_tip[cn]) continue;
            if (i == pred) continue;  // Skip self-loops

            result.uf.unite(i, pred);
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ci = scc.component[i];
        if (is_forward_tip[ci] || is_backward_tip[ci]) {
            result.ignorable_tip[i] = true;
        }
    }

    // Path compression for stable representatives
    for (std::size_t i = 0; i < n; ++i) {
        result.uf.find(i);
    }

    return result;
}

void chainCycles(const AlnGraph& graph, TipFoldingResult& tip_result) {
    const std::size_t n = graph.nodeCount();
    if (n == 0) return;

    constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();
    constexpr std::size_t MULTIPLE = std::numeric_limits<std::size_t>::max() - 1;

    for (std::size_t i = 0; i < n; ++i) {
        // Find unique forward neighbor (skip tips and self-loops)
        std::size_t unique_fw = NONE;
        for (const auto u : graph.outgoing(i)) {
            if (tip_result.ignorable_tip[u]) continue;  // Skip tips
            if (u == i) continue;  // Skip self-loops

            if (unique_fw == NONE) {
                unique_fw = u;
            } else if (u != unique_fw) {
                unique_fw = MULTIPLE;  // Multiple neighbors
                break;
            }
        }

        // Find unique backward neighbor (skip tips and self-loops)
        std::size_t unique_bw = NONE;
        for (const auto u : graph.incoming(i)) {
            if (tip_result.ignorable_tip[u]) continue;  // Skip tips
            if (u == i) continue;  // Skip self-loops

            if (unique_bw == NONE) {
                unique_bw = u;
            } else if (u != unique_bw) {
                unique_bw = MULTIPLE;  // Multiple neighbors
                break;
            }
        }

        // Must be exactly one neighbor in each direction and they must be the same
        if (unique_fw != unique_bw) continue;
        if (unique_fw == NONE) continue;      // No neighbors
        if (unique_fw == MULTIPLE) continue;  // Multiple neighbors
        if (unique_bw == NONE) continue;
        if (unique_bw == MULTIPLE) continue;

        // Found a 2-cycle: i ⇄ unique_fw
        tip_result.ignorable_tip[i] = true;
        tip_result.uf.unite(i, unique_fw);
    }
}

namespace {

// Detect a superbubble starting at `start`; return exit node if found.
std::optional<std::size_t> findSuperbubble(const AlnGraph& graph,
                                            std::size_t start,
                                            const std::vector<bool>& ignorable_tip) {
    // Onodera et al. 2013 superbubble detection algorithm
    std::vector<std::size_t> S;
    S.push_back(start);

    std::unordered_set<std::size_t> visited;
    std::unordered_set<std::size_t> seen;
    seen.insert(start);

    while (!S.empty()) {
        const std::size_t v = S.back();
        S.pop_back();

        if (seen.count(v) == 1) seen.erase(v);
        if (visited.count(v) == 1) continue;
        visited.insert(v);

        // If v has no outgoing edges, not a valid superbubble
        if (graph.outgoing(v).empty()) {
            return std::nullopt;
        }

        for (const auto u : graph.outgoing(v)) {
            if (ignorable_tip[u]) continue;  // Skip tips
            if (u == v) continue;            // Skip self-loops
            if (u == start) return std::nullopt;  // Cycle back to start

            if (visited.count(u) == 0) {
                seen.insert(u);

                // Check if u has any unvisited predecessors (excluding self-loops and tips)
                bool has_nonvisited_parent = false;
                for (const auto w : graph.incoming(u)) {
                    if (w == u) continue;            // Skip self-loops
                    if (ignorable_tip[w]) continue;  // Skip tips
                    if (visited.count(w) == 0) {
                        has_nonvisited_parent = true;
                        break;
                    }
                }

                // Only add to stack if all predecessors have been visited
                if (!has_nonvisited_parent) {
                    S.push_back(u);
                }
            }
        }

        // Check if we found the exit node (all paths reconverge)
        if (S.size() == 1 && seen.size() == 1 && seen.count(S[0]) == 1) {
            const std::size_t t = S.back();

            // Verify no cycle back to start from exit
            for (const auto u : graph.outgoing(t)) {
                if (u == start) {
                    return std::nullopt;
                }
            }

            return t;  // Found valid superbubble exit
        }
    }

    return std::nullopt;  // No superbubble found
}

// Helper: Traverse all nodes reachable from start without passing through end
// Returns set of node IDs
//
// Legacy reference: alignment_graph.cpp lines 1237-1246 (BFS traversal in _chain_bubble_ga)
std::vector<std::size_t> traverseBubbleInterior(const AlnGraph& graph,
                                                 std::size_t start,
                                                 std::size_t end,
                                                 const std::vector<bool>& ignorable_tip) {
    // BFS/DFS traversal from start, stopping at end (don't traverse past exit)
    std::unordered_set<std::size_t> visited;
    std::vector<std::size_t> stack;
    stack.push_back(start);

    while (!stack.empty()) {
        const std::size_t node = stack.back();
        stack.pop_back();

        if (visited.count(node) == 1) continue;
        if (ignorable_tip[node]) continue;

        visited.insert(node);

        for (const auto succ : graph.outgoing(node)) {
            if (visited.count(succ) == 1) continue;
            if (succ == end) continue;  // Don't traverse past exit

            stack.push_back(succ);
        }
    }

    // Convert set to vector
    return std::vector<std::size_t>(visited.begin(), visited.end());
}

}  // namespace

SuperbubbleResult chainSuperbubbles(const AlnGraph& graph,
                                     const SccResult& scc,
                                     const TipFoldingResult& tip_result) {
    SuperbubbleResult result(tip_result.uf, tip_result.ignorable_tip);

    const std::size_t n = graph.nodeCount();
    if (n == 0) return result;

    for (std::size_t start = 0; start < n; ++start) {
        // Skip nodes in SCCs (cycles)
        if (scc.in_scc[start]) continue;

        // Skip tips
        if (result.ignorable_tip[start]) continue;

        // Detect superbubble from start
        auto end_opt = findSuperbubble(graph, start, result.ignorable_tip);
        if (!end_opt.has_value()) continue;

        const std::size_t end = end_opt.value();

        // Verify topological constraint: end must be strictly downstream
        if (scc.component[end] <= scc.component[start]) {
            continue;  // Not a valid forward superbubble
        }

        // Union start and end
        result.uf.unite(start, end);

        // Traverse bubble interior and union all nodes with start
        const auto interior = traverseBubbleInterior(graph, start, end, result.ignorable_tip);
        for (const auto node : interior) {
            result.uf.unite(start, node);
        }
    }

    return result;
}

}  // namespace piru::index
