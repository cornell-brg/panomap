/**
 * panolayout.cpp
 *
 * Standalone PG-SGD 1D pangenome layout tool. Reads a GFA (P or W lines),
 * computes a 1D coordinate per node via path-guided SGD, and writes:
 *   - <prefix>.coords.tsv      one row per node: id, canon_start, canon_end, cc, length
 *   - <prefix>.path_steps.tsv  one row per step: path, idx, node, orient, cum_bp, c_start, c_end
 *
 * Algorithm ports piru's src/core/index/sort_1d.cpp without the piru
 * graph dependency, so it can be hacked on in isolation: cost function,
 * sampling distribution, initialization, post-processing, etc.
 *
 * Build (via parent CMake):  cmake --build build --target panolayout
 *
 * SPDX-License-Identifier: MIT
 */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <future>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "dirty_zipfian_int_distribution.h"

/* ---------------------------- Data ---------------------------- */

struct Node {
  std::uint32_t length;  // bp
};

struct PathStep {
  std::uint32_t node_id;
  std::uint8_t orient;     // 0 = '+', 1 = '-'
  std::uint64_t cum_bp;    // cumulative bp at START of this node on the path
};

struct Path {
  std::string name;
  std::vector<PathStep> steps;
  std::uint64_t total_bp{0};
  std::uint64_t strain_bp_offset{0};  // W-line start coord (0 for P-lines).
                                      // Lets multiple W-line fragments for one
                                      // strain share a single bp axis when plotted.
};

struct Graph {
  std::vector<Node> nodes;
  std::unordered_map<std::string, std::uint32_t> name_to_id;
  std::vector<Path> paths;
};

/* ----------------------- Parsing helpers ----------------------- */

static inline void split_tabs(const std::string& s, std::vector<std::string>& out) {
  out.clear();
  std::size_t a = 0;
  for (std::size_t i = 0; i <= s.size(); ++i) {
    if (i == s.size() || s[i] == '\t') {
      out.emplace_back(s, a, i - a);
      a = i + 1;
    }
  }
}

static std::uint32_t intern_node(Graph& g, const std::string& name, std::uint32_t length) {
  auto it = g.name_to_id.find(name);
  if (it != g.name_to_id.end()) return it->second;
  std::uint32_t id = static_cast<std::uint32_t>(g.nodes.size());
  g.nodes.push_back({length});
  g.name_to_id.emplace(name, id);
  return id;
}

static void parse_p_steps(const std::string& steps_str, std::vector<std::pair<std::string, char>>& out) {
  out.clear();
  std::size_t a = 0;
  for (std::size_t i = 0; i <= steps_str.size(); ++i) {
    if (i == steps_str.size() || steps_str[i] == ',') {
      if (i > a + 1) {
        char orient = steps_str[i - 1];
        out.emplace_back(steps_str.substr(a, i - 1 - a), orient);
      }
      a = i + 1;
    }
  }
}

/* W-line walk format: >id<id>id... (> = +, < = -). */
static void parse_w_walk(const std::string& walk, std::vector<std::pair<std::string, char>>& out) {
  out.clear();
  std::size_t i = 0;
  while (i < walk.size()) {
    char c = walk[i];
    if (c != '>' && c != '<') { ++i; continue; }
    char orient = (c == '>') ? '+' : '-';
    std::size_t j = i + 1;
    while (j < walk.size() && walk[j] != '>' && walk[j] != '<') ++j;
    out.emplace_back(walk.substr(i + 1, j - i - 1), orient);
    i = j;
  }
}

static bool parse_gfa(const std::string& path, Graph& g) {
  std::ifstream f(path);
  if (!f) { std::fprintf(stderr, "ERR: cannot open %s\n", path.c_str()); return false; }

  std::string line;
  std::vector<std::string> fields;
  std::vector<std::pair<std::string, char>> raw_steps;
  std::uint64_t n_s = 0, n_p = 0, n_w = 0;

  // Two passes: first collect all S-lines, then P/W (some GFAs are unsorted).
  std::vector<std::pair<std::string, std::uint32_t>> seg_pending;  // (name, length)
  struct PendingPath { std::string name; std::string body; char mode; std::uint64_t start; };
  std::vector<PendingPath> path_pending;

  while (std::getline(f, line)) {
    if (line.empty()) continue;
    char rec = line[0];
    if (rec == 'S') {
      split_tabs(line, fields);
      if (fields.size() < 3) continue;
      seg_pending.emplace_back(fields[1], static_cast<std::uint32_t>(fields[2].size()));
      ++n_s;
    } else if (rec == 'P') {
      split_tabs(line, fields);
      if (fields.size() < 3) continue;
      path_pending.push_back({fields[1], fields[2], 'P', 0});
      ++n_p;
    } else if (rec == 'W') {
      split_tabs(line, fields);
      if (fields.size() < 7) continue;
      std::string pname = fields[1] + "#" + fields[2] + "#" + fields[3];
      std::uint64_t start = std::strtoull(fields[4].c_str(), nullptr, 10);
      path_pending.push_back({pname, fields[6], 'W', start});
      ++n_w;
    }
  }

  // Materialize nodes in S-line order so IDs are stable.
  g.nodes.reserve(seg_pending.size());
  for (auto& [name, len] : seg_pending) intern_node(g, name, len);

  // Materialize paths.
  g.paths.reserve(path_pending.size());
  for (auto& pp : path_pending) {
    Path p;
    p.name = pp.name;
    p.strain_bp_offset = pp.start;
    if (pp.mode == 'P') parse_p_steps(pp.body, raw_steps);
    else parse_w_walk(pp.body, raw_steps);

    std::uint64_t cum = 0;
    p.steps.reserve(raw_steps.size());
    for (auto& [nname, orient] : raw_steps) {
      auto it = g.name_to_id.find(nname);
      if (it == g.name_to_id.end()) {
        std::fprintf(stderr, "WARN: path %s references unknown segment %s; skipping step\n",
                     pp.name.c_str(), nname.c_str());
        continue;
      }
      std::uint32_t nid = it->second;
      p.steps.push_back({nid, static_cast<std::uint8_t>(orient == '+' ? 0 : 1), cum});
      cum += g.nodes[nid].length;
    }
    p.total_bp = cum;
    g.paths.push_back(std::move(p));
  }

  std::fprintf(stderr, "[gfa] %s: %llu segments, %llu P-lines, %llu W-lines, %zu paths\n",
               path.c_str(),
               static_cast<unsigned long long>(n_s),
               static_cast<unsigned long long>(n_p),
               static_cast<unsigned long long>(n_w),
               g.paths.size());
  return true;
}

/* ----------------------- Connected components ----------------------- */
/* Path-based: two nodes share a cc iff there is a path that visits both. */

struct UnionFind {
  std::vector<std::uint32_t> parent;
  std::vector<std::uint32_t> rank_;
  void init(std::size_t n) { parent.resize(n); rank_.assign(n, 0); for (std::size_t i=0;i<n;++i) parent[i] = static_cast<std::uint32_t>(i); }
  std::uint32_t find(std::uint32_t x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
  }
  void unite(std::uint32_t a, std::uint32_t b) {
    a = find(a); b = find(b);
    if (a == b) return;
    if (rank_[a] < rank_[b]) std::swap(a, b);
    parent[b] = a;
    if (rank_[a] == rank_[b]) ++rank_[a];
  }
};

static std::vector<std::uint32_t> compute_cc(const Graph& g) {
  UnionFind uf;
  uf.init(g.nodes.size());
  for (const auto& p : g.paths) {
    for (std::size_t i = 1; i < p.steps.size(); ++i) {
      uf.unite(p.steps[i-1].node_id, p.steps[i].node_id);
    }
  }
  // Re-number cc ids to be dense 0..K-1
  std::unordered_map<std::uint32_t, std::uint32_t> remap;
  std::vector<std::uint32_t> cc(g.nodes.size());
  std::uint32_t next = 0;
  for (std::uint32_t i = 0; i < g.nodes.size(); ++i) {
    std::uint32_t r = uf.find(i);
    auto it = remap.find(r);
    if (it == remap.end()) { remap[r] = next; cc[i] = next; ++next; }
    else cc[i] = it->second;
  }
  std::fprintf(stderr, "[cc]  %u components\n", next);
  return cc;
}

/* ----------------------- PG-SGD ----------------------- */

struct Config {
  std::uint64_t iter_max{100};
  double theta{0.99};
  double first_cooling_frac{0.5};
  std::uint64_t space{0};        // 0 = use longest path
  std::uint64_t space_max{1000};
  std::uint64_t seed{42};
  std::size_t num_threads{1};
  double eps{0.01};              // min learning rate
  bool signed_loss{false};       // experimental: use signed bp distance
  std::string init{"node-order"}; // "node-order" (piru-style flat concat) or "path-progress"
};

// Per-iteration learning rate schedule (matches odgi's path_sgd_schedule).
static std::vector<double> compute_lr_schedule(std::uint64_t iter_max, double eta_max, double eps) {
  std::vector<double> etas(iter_max);
  if (iter_max == 0) return etas;
  double lambda = std::log(eta_max / eps) / (iter_max - 1);
  for (std::uint64_t i = 0; i < iter_max; ++i) {
    etas[i] = eta_max * std::exp(-lambda * i);
  }
  return etas;
}

/* The X[] array holds 2 coords per node: X[2*id] = node start, X[2*id+1] = node end. */
static void init_coords_node_order(const Graph& g, std::vector<double>& X) {
  // Piru-style: flat concat by node-ID order. Every node starts where the
  // previous node ended; every node lives on a single pre-flattened line at
  // init. SGD then warps the line toward path-distance consistency.
  X.assign(2 * g.nodes.size(), 0.0);
  double cum = 0.0;
  for (std::uint32_t i = 0; i < g.nodes.size(); ++i) {
    X[2*i]     = cum;
    X[2*i + 1] = cum + static_cast<double>(g.nodes[i].length);
    cum += static_cast<double>(g.nodes[i].length);
  }
}

static void init_coords_path_progress(const Graph& g, std::vector<double>& X) {
  // First-path-to-visit-wins; node coords carry their path-local cum_bp.
  // Unique-to-late-path segments end up offset because their cum_bp restarts
  // from 0 while shared anchors are pinned at the first-path's cum_bp.
  X.assign(2 * g.nodes.size(), 0.0);
  std::vector<bool> placed(g.nodes.size(), false);
  for (const auto& p : g.paths) {
    for (const auto& s : p.steps) {
      if (placed[s.node_id]) continue;
      double a = static_cast<double>(s.cum_bp);
      double b = a + static_cast<double>(g.nodes[s.node_id].length);
      if (s.orient == 1) std::swap(a, b);
      X[2*s.node_id]     = a;
      X[2*s.node_id + 1] = b;
      placed[s.node_id] = true;
    }
  }
}

static void run_sgd(const Graph& g, const Config& cfg, std::vector<double>& X) {
  // Step index for sampling
  struct StepRec { std::uint32_t node_id; std::int64_t ref_pos; std::uint32_t length; };
  std::vector<std::vector<StepRec>> path_index(g.paths.size());
  std::vector<std::size_t> valid_paths;
  std::uint64_t total_steps = 0;
  std::size_t max_path_steps = 0;
  for (std::size_t i = 0; i < g.paths.size(); ++i) {
    const auto& p = g.paths[i];
    if (p.steps.size() < 2) continue;
    auto& v = path_index[i];
    v.reserve(p.steps.size());
    for (const auto& s : p.steps)
      v.push_back({s.node_id, static_cast<std::int64_t>(s.cum_bp), g.nodes[s.node_id].length});
    valid_paths.push_back(i);
    total_steps += v.size();
    max_path_steps = std::max(max_path_steps, v.size());
  }
  if (valid_paths.empty()) { std::fprintf(stderr, "[sgd] no valid paths\n"); return; }

  std::uint64_t space = cfg.space ? cfg.space : max_path_steps;
  std::uint64_t space_max = std::min(cfg.space_max, space);

  // Precompute zeta values for Zipfian distributions
  std::vector<double> zetas(space_max + 2 + (space - space_max) / space_max + 1, 0.0);
  for (std::uint64_t i = 1; i <= space_max + 1; ++i) {
    double z = 0.0;
    for (std::uint64_t k = 1; k <= i; ++k) z += 1.0 / std::pow(static_cast<double>(k), cfg.theta);
    zetas[i] = z;
  }

  double eta_max = static_cast<double>(max_path_steps) * static_cast<double>(max_path_steps);
  auto etas = compute_lr_schedule(cfg.iter_max, eta_max, cfg.eps);
  std::uint64_t first_cooling = static_cast<std::uint64_t>(cfg.first_cooling_frac * cfg.iter_max);

  // Cumulative step counts for weighted path sampling
  std::vector<std::uint64_t> cum_steps;
  cum_steps.reserve(valid_paths.size());
  std::uint64_t cum = 0;
  for (auto p : valid_paths) { cum += path_index[p].size(); cum_steps.push_back(cum); }

  std::fprintf(stderr, "[sgd] %llu valid paths, %llu total steps, max=%zu, space=%llu, theta=%.3f\n",
               static_cast<unsigned long long>(valid_paths.size()),
               static_cast<unsigned long long>(total_steps), max_path_steps,
               static_cast<unsigned long long>(space), cfg.theta);
  std::fprintf(stderr, "[sgd] iter_max=%llu, threads=%zu, first_cooling=%llu, signed=%s\n",
               static_cast<unsigned long long>(cfg.iter_max), cfg.num_threads,
               static_cast<unsigned long long>(first_cooling),
               cfg.signed_loss ? "yes" : "no");

  auto one_update = [&](std::mt19937_64& rng, double eta, double adj_theta, bool cooling) {
    std::uniform_int_distribution<int> flip(0, 1);
    std::uniform_int_distribution<std::uint64_t> global_step(0, total_steps - 1);

    std::uint64_t gi = global_step(rng);
    auto it = std::lower_bound(cum_steps.begin(), cum_steps.end(), gi + 1);
    std::size_t pr = static_cast<std::size_t>(it - cum_steps.begin());
    std::size_t pid = valid_paths[pr];
    const auto& steps = path_index[pid];
    std::size_t off = (pr > 0) ? cum_steps[pr - 1] : 0;
    std::size_t a = static_cast<std::size_t>(gi - off);
    std::size_t n = steps.size();

    std::size_t b;
    if (cooling || flip(rng)) {  // zipfian
      bool back = (a > 0 && flip(rng)) || a == n - 1;
      std::uint64_t js = back ? std::min(space, static_cast<std::uint64_t>(a))
                              : std::min(space, static_cast<std::uint64_t>(n - a - 1));
      if (js == 0) return;
      std::uint64_t zi = js;
      if (js > space_max) zi = space_max + 1 + (js - space_max) / space_max;
      if (zi >= zetas.size()) zi = zetas.size() - 1;
      dirtyzipf::dirty_zipfian_int_distribution<std::uint64_t>::param_type zp(1, js, adj_theta, zetas[zi]);
      dirtyzipf::dirty_zipfian_int_distribution<std::uint64_t> z(zp);
      std::uint64_t step = z(rng);
      b = back ? a - step : a + step;
    } else {
      std::uniform_int_distribution<std::size_t> uni(0, n - 1);
      b = uni(rng);
    }
    if (b >= n || b == a) return;

    std::uint32_t ni = steps[a].node_id, nj = steps[b].node_id;
    if (ni == nj) return;

    std::uint32_t oi = flip(rng), oj = flip(rng);
    double pa = static_cast<double>(steps[a].ref_pos) + (oi ? steps[a].length : 0);
    double pb = static_cast<double>(steps[b].ref_pos) + (oj ? steps[b].length : 0);

    double xa = X[2*ni + oi], xb = X[2*nj + oj];
    double dx = xa - xb;
    if (dx == 0) dx = 1e-9;
    double mag = std::abs(dx);

    double target;
    if (cfg.signed_loss) {
      // Signed: target = pa - pb (carries direction). Want dx == target.
      double signed_target = pa - pb;
      double err = dx - signed_target;
      double w_ij = 1.0 / std::max(1.0, std::abs(signed_target));
      double mu = std::min(1.0, eta * w_ij);
      double delta_full = mu * err;
      X[2*ni + oi] -= 0.5 * delta_full;
      X[2*nj + oj] += 0.5 * delta_full;
      return;
    } else {
      target = std::abs(pa - pb);
      if (target == 0) return;
      double w_ij = 1.0 / target;
      double mu = std::min(1.0, eta * w_ij);
      double delta = mu * (mag - target) / 2.0;
      double r = delta / mag;
      X[2*ni + oi] -= r * dx;
      X[2*nj + oj] += r * dx;
    }
  };

  std::size_t nthr = std::max<std::size_t>(1, cfg.num_threads);

  for (std::uint64_t it = 0; it < cfg.iter_max; ++it) {
    double eta = etas[it];
    bool cooling = it > first_cooling;
    // Piru/odgi pattern: during cooling, drop theta toward 0 (near-uniform sampling)
    // to allow broader exploration and escape from local minima found in the early
    // (Zipf-concentrated) phase.
    double adj_theta = cooling ? 0.001 : cfg.theta;

    if (nthr <= 1) {
      std::mt19937_64 rng(cfg.seed + it * 0x100000001ULL);
      for (std::uint64_t s = 0; s < total_steps; ++s) one_update(rng, eta, adj_theta, cooling);
    } else {
      std::vector<std::future<void>> futs;
      futs.reserve(nthr);
      std::uint64_t per = total_steps / nthr;
      for (std::size_t t = 0; t < nthr; ++t) {
        std::uint64_t lo = t * per;
        std::uint64_t hi = (t + 1 == nthr) ? total_steps : (t + 1) * per;
        futs.emplace_back(std::async(std::launch::async, [&, lo, hi, t]() {
          std::mt19937_64 rng(cfg.seed + it * 0x100000001ULL + t * 0x9E3779B97F4A7C15ULL);
          for (std::uint64_t s = lo; s < hi; ++s) one_update(rng, eta, adj_theta, cooling);
        }));
      }
      for (auto& f : futs) f.get();
    }
    if ((it + 1) % std::max<std::uint64_t>(1, cfg.iter_max / 10) == 0) {
      std::fprintf(stderr, "[sgd] iter=%llu/%llu eta=%.3g cool=%d\n",
                   static_cast<unsigned long long>(it + 1),
                   static_cast<unsigned long long>(cfg.iter_max), eta, cooling ? 1 : 0);
    }
  }
}

/* ----------------------- Output ----------------------- */

static void write_coords_tsv(const Graph& g, const std::vector<double>& X,
                             const std::vector<std::uint32_t>& cc, const std::string& path) {
  std::ofstream out(path);
  out << "node_id\tcanon_start\tcanon_end\tcomponent_id\tlength\n";
  for (std::uint32_t i = 0; i < g.nodes.size(); ++i) {
    out << i << '\t' << X[2*i] << '\t' << X[2*i + 1] << '\t' << cc[i] << '\t' << g.nodes[i].length << '\n';
  }
  std::fprintf(stderr, "[out] wrote %s (%zu nodes)\n", path.c_str(), g.nodes.size());
}

static void write_path_steps_tsv(const Graph& g, const std::vector<double>& X, const std::string& path) {
  std::ofstream out(path);
  // cum_bp column = strain-wide bp position = strain_bp_offset + within-fragment cum_bp.
  // For P-lines (strain_bp_offset=0) this is just cum_bp. For W-lines, fragments that
  // share a path_name (same strain) now lie on a single bp axis with gaps where MC
  // clipped, so the plotter sees one diagonal per strain.
  out << "path_name\tstep_idx\tnode_id\torient\tcum_bp\tcoord_start\tcoord_end\n";
  for (const auto& p : g.paths) {
    for (std::size_t i = 0; i < p.steps.size(); ++i) {
      const auto& s = p.steps[i];
      out << p.name << '\t' << i << '\t' << s.node_id << '\t'
          << (s.orient == 0 ? '+' : '-') << '\t' << (p.strain_bp_offset + s.cum_bp) << '\t'
          << X[2*s.node_id] << '\t' << X[2*s.node_id + 1] << '\n';
    }
  }
  std::fprintf(stderr, "[out] wrote %s\n", path.c_str());
}

/* ----------------------- CLI ----------------------- */

static void usage() {
  std::fprintf(stderr,
    "Usage: panolayout <graph.gfa> -o <prefix> [options]\n"
    "Options:\n"
    "  -o <prefix>           Output prefix; writes <prefix>.coords.tsv and <prefix>.path_steps.tsv\n"
    "  --iter N              SGD iterations (default 100)\n"
    "  --threads T           Threads (default 1)\n"
    "  --theta F             Zipfian theta (default 0.99)\n"
    "  --space N             Max jump distance in path steps (0 = longest path; default 0)\n"
    "  --first-cooling F     Fraction of iters before cooling kicks in (default 0.5)\n"
    "  --seed S              RNG seed (default 42)\n"
    "  --signed-loss         EXPERIMENTAL: use signed bp distance in loss (orientation-aware)\n"
    "  --init MODE           Init: 'node-order' (piru-style flat concat; default) or 'path-progress'\n"
  );
}

int main(int argc, char** argv) {
  if (argc < 2) { usage(); return 1; }
  std::string gfa, prefix;
  Config cfg;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](int k){ if (i + k >= argc) { usage(); std::exit(1); } };
    if (a == "-h" || a == "--help") { usage(); return 0; }
    else if (a == "-o") { need(1); prefix = argv[++i]; }
    else if (a == "--iter") { need(1); cfg.iter_max = std::stoull(argv[++i]); }
    else if (a == "--threads") { need(1); cfg.num_threads = std::stoull(argv[++i]); }
    else if (a == "--theta") { need(1); cfg.theta = std::stod(argv[++i]); }
    else if (a == "--space") { need(1); cfg.space = std::stoull(argv[++i]); }
    else if (a == "--first-cooling") { need(1); cfg.first_cooling_frac = std::stod(argv[++i]); }
    else if (a == "--seed") { need(1); cfg.seed = std::stoull(argv[++i]); }
    else if (a == "--signed-loss") { cfg.signed_loss = true; }
    else if (a == "--init") { need(1); cfg.init = argv[++i]; }
    else if (a[0] != '-' && gfa.empty()) { gfa = a; }
    else { std::fprintf(stderr, "ERR: unknown arg '%s'\n", a.c_str()); usage(); return 1; }
  }
  if (gfa.empty() || prefix.empty()) { usage(); return 1; }

  Graph g;
  if (!parse_gfa(gfa, g)) return 1;
  if (g.nodes.empty()) { std::fprintf(stderr, "ERR: no nodes parsed\n"); return 1; }

  auto cc = compute_cc(g);
  std::vector<double> X;
  if (cfg.init == "node-order") init_coords_node_order(g, X);
  else if (cfg.init == "path-progress") init_coords_path_progress(g, X);
  else { std::fprintf(stderr, "ERR: unknown --init '%s'\n", cfg.init.c_str()); return 1; }
  std::fprintf(stderr, "[init] %s\n", cfg.init.c_str());
  run_sgd(g, cfg, X);

  // Piru-style post-SGD: zero-base each component (per-cc min subtracted from
  // all nodes in that cc). Layout structure unchanged, global offset removed.
  std::uint32_t ncc = 0;
  for (auto v : cc) ncc = std::max(ncc, v + 1);
  std::vector<double> comp_min(ncc, std::numeric_limits<double>::max());
  for (std::uint32_t i = 0; i < g.nodes.size(); ++i) {
    comp_min[cc[i]] = std::min(comp_min[cc[i]], std::min(X[2*i], X[2*i+1]));
  }
  for (std::uint32_t i = 0; i < g.nodes.size(); ++i) {
    X[2*i]     -= comp_min[cc[i]];
    X[2*i + 1] -= comp_min[cc[i]];
  }
  std::fprintf(stderr, "[post] zero-based %u component(s)\n", ncc);

  write_coords_tsv(g, X, cc, prefix + ".coords.tsv");
  write_path_steps_tsv(g, X, prefix + ".path_steps.tsv");
  return 0;
}
