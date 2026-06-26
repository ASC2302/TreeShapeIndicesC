// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <string>
#include <functional>
#include <cctype>
#include <cstdlib>

using namespace Rcpp;

using Node = int;
using BranchLength = double;
using Abundance = double;

struct BranchInfo {
  Node startNode;
  Node endNode;
  BranchLength branchLength; // (possibly shortened) Branch_length in R
  BranchLength x;            // distance from reference node to end of this branch (possibly shifted)
};

struct TISResult {
  double T = 0.0;
  std::vector<double> S;
  std::vector<double> x;
};

struct SIAResult {
  std::vector<double> S;
  std::vector<double> x;
  std::vector<std::vector<double>> abundances;
  std::vector<std::vector<int>> nodes;
};

struct DJSeriesResult {
  bool individual = false;
  std::vector<double> x;
  std::vector<double> In;
  std::vector<double> D1;
  std::vector<double> J1;
  std::vector<double> D0;
};

struct IntegralResult {
  bool individual = false;
  double value = 0.0;
  std::array<double, 3> values{{0.0, 0.0, 0.0}};
};

static inline bool approx_equal(double a, double b, double tol = 1e-12) {
  double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
  return std::fabs(a - b) <= tol * scale;
}

static inline std::string to_upper_copy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return static_cast<char>(std::toupper(c)); });
  return s;
}

static inline bool is_null_or_false(const Rcpp::RObject& x) {
  if (x.isNULL()) return true;
  if (TYPEOF(x) == LGLSXP && Rf_length(x) > 0) {
    return (Rcpp::as<Rcpp::LogicalVector>(x)[0] == FALSE);
  }
  return false;
}

static inline bool parse_int_string(const std::string& s, int& out) {
  if (s.empty()) return false;
  char* end = nullptr;
  long v = std::strtol(s.c_str(), &end, 10);
  if (end == s.c_str() || *end != '\0') return false;
  if (v < 0 || v > std::numeric_limits<int>::max()) return false;
  out = static_cast<int>(v);
  return true;
}

// Numeric node abundances indexed by node number. This replaces the repeated
// int -> string conversion and unordered_map lookup used in the original code.
struct AbundanceLookup {
  std::vector<double> value; // 1-indexed when node labels are standard phylo ids

  inline double get(int node) const {
    if (node < 0 || node >= static_cast<int>(value.size())) return 0.0;
    return value[node];
  }
};

static inline AbundanceLookup preprocess_abundances_vec(const NumericVector& abundances) {
  AbundanceLookup out;

  Rcpp::RObject nm_obj = abundances.attr("names");
  bool has_names = !nm_obj.isNULL();

  if (!has_names) {
    out.value.assign(abundances.size() + 1, 0.0);
    for (int i = 0; i < abundances.size(); ++i) out.value[i + 1] = abundances[i];
    return out;
  }

  CharacterVector nms = Rcpp::as<CharacterVector>(nm_obj);
  int max_node = 0;
  std::vector<int> ids(abundances.size(), -1);
  for (int i = 0; i < abundances.size(); ++i) {
    Rcpp::String si = nms[i];
    if (si == NA_STRING) continue;
    int id = -1;
    if (parse_int_string(std::string(si), id)) {
      ids[i] = id;
      max_node = std::max(max_node, id);
    }
  }

  out.value.assign(max_node + 1, 0.0);
  for (int i = 0; i < abundances.size(); ++i) {
    int id = ids[i];
    if (id >= 0 && id < static_cast<int>(out.value.size())) out.value[id] = abundances[i];
  }
  return out;
}

struct TreeCache {
  std::vector<std::vector<std::pair<int,double>>> children;
  std::vector<int> parent;
  std::vector<double> len_to_parent;
  int root = NA_INTEGER;
  int max_node = 0;
};

static inline TreeCache build_cache(const List& tree) {
  IntegerMatrix edge = tree["edge"];
  NumericVector el = tree["edge.length"];

  TreeCache tc;

  int max_node = 0;
  for (int i = 0; i < edge.nrow(); ++i) {
    max_node = std::max(max_node, edge(i,0));
    max_node = std::max(max_node, edge(i,1));
  }
  tc.max_node = max_node;
  tc.children.assign(max_node + 1, std::vector<std::pair<int,double>>());
  tc.parent.assign(max_node + 1, NA_INTEGER);
  tc.len_to_parent.assign(max_node + 1, 0.0);
  std::vector<char> appears_as_child(max_node + 1, 0);
  std::vector<char> appears_as_parent(max_node + 1, 0);

  for (int i = 0; i < edge.nrow(); ++i) {
    int p = edge(i,0);
    int c = edge(i,1);
    double l = el[i];

    if (p >= 0 && p <= max_node && c >= 0 && c <= max_node) {
      tc.children[p].push_back({c,l});
      tc.parent[c] = p;
      tc.len_to_parent[c] = l;
      appears_as_child[c] = 1;
      appears_as_parent[p] = 1;
    }
  }

  int root = NA_INTEGER;
  for (int i = 0; i <= max_node; ++i) {
    if (appears_as_parent[i] && !appears_as_child[i]) { root = i; break; }
  }
  if (root == NA_INTEGER && edge.nrow() > 0) root = edge(0,0);
  tc.root = root;

  return tc;
}

static inline double distance_cpp(const TreeCache& tc, int top_node, int bottom_node) {
  double dist = 0.0;
  int cur = bottom_node;
  while (cur != top_node) {
    if (cur < 0 || cur >= static_cast<int>(tc.parent.size())) break;
    int p = tc.parent[cur];
    if (p == NA_INTEGER) break; // not reachable
    dist += tc.len_to_parent[cur];
    cur = p;
  }
  return dist;
}

// R-equivalent: all descendant edges in subtree(root=node), with x = distance(node, endNode).
static inline std::vector<BranchInfo> descendant_branches_with_x(const TreeCache& tc, int node) {
  std::vector<BranchInfo> out;
  if (node < 0 || node >= static_cast<int>(tc.children.size())) return out;

  out.reserve(std::max(0, tc.max_node));
  std::vector<std::pair<int,double>> stack;
  stack.reserve(std::max(1, tc.max_node));
  stack.push_back({node, 0.0}); // (current node, distance from root-of-subtree to this node)

  while (!stack.empty()) {
    int u = stack.back().first;
    double du = stack.back().second;
    stack.pop_back();

    if (u < 0 || u >= static_cast<int>(tc.children.size())) continue;
    const auto& ch = tc.children[u];
    for (const auto& pr : ch) {
      int v = pr.first;
      double len = pr.second;
      double dv = du + len;
      out.push_back({u, v, len, dv});
      stack.push_back({v, dv});
    }
  }
  return out;
}

static inline DataFrame tis_to_dataframe(const TISResult& r) {
  return DataFrame::create(_["S_i"] = r.S, _["x"] = r.x);
}

static inline TISResult dataframe_to_tis(const DataFrame& df) {
  NumericVector S = df["S_i"];
  NumericVector X = df["x"];
  TISResult out;
  out.S.assign(S.begin(), S.end());
  out.x.assign(X.begin(), X.end());
  out.T = 0.0;
  return out;
}

static inline NumericVector make_named_numeric(const std::vector<double>& vals,
                                               const std::vector<int>& nodes) {
  NumericVector out(vals.size());
  CharacterVector nms(vals.size());
  for (int i = 0; i < static_cast<int>(vals.size()); ++i) {
    out[i] = vals[i];
    nms[i] = std::to_string(nodes[i]);
  }
  out.names() = nms;
  return out;
}

static inline List sia_to_list(const SIAResult& r) {
  List abund_list;
  for (int i = 0; i < static_cast<int>(r.abundances.size()); ++i) {
    abund_list.push_back(make_named_numeric(r.abundances[i], r.nodes[i]));
  }
  DataFrame DF_S_i = DataFrame::create(_["S_i"] = r.S, _["x"] = r.x);
  return List::create(DF_S_i, abund_list);
}

static inline DataFrame series_to_individual_df(const std::vector<double>& values,
                                                const std::vector<double>& x,
                                                const std::string& value_name) {
  if (value_name == "In_i_a") return DataFrame::create(_["In_i_a"] = values, _["x"] = x);
  return DataFrame::create(_["In"] = values, _["x"] = x);
}

static inline TISResult compute_T_i_S_i_cached(const TreeCache& tc,
                                               int node,
                                               const AbundanceLookup& abund) {
  TISResult out;
  if (node < 0 || node >= static_cast<int>(tc.children.size())) return out;

  std::vector<std::pair<int,double>> desc; // (endNode, branchLength)
  const auto& ch = tc.children[node];
  desc.reserve(ch.size());
  for (const auto& pr : ch) desc.push_back(pr);

  if (desc.empty()) return out;

  std::sort(desc.begin(), desc.end(), [](const auto& a, const auto& b){
    return a.second < b.second;
  });

  out.S.reserve(desc.size());
  out.x.reserve(desc.size());

  double remaining_sum = 0.0;
  for (const auto& pr : desc) remaining_sum += abund.get(pr.first);

  double prev_x = 0.0;
  int i = 0;
  while (i < static_cast<int>(desc.size())) {
    double len = desc[i].second;
    out.S.push_back(remaining_sum);
    out.x.push_back(len);
    out.T += remaining_sum * (len - prev_x);
    prev_x = len;

    int j = i;
    while (j < static_cast<int>(desc.size()) && approx_equal(desc[j].second, len)) {
      remaining_sum -= abund.get(desc[j].first);
      ++j;
    }
    i = j;
  }

  return out;
}

// [[Rcpp::export]]
List compute_T_i_S_i(List tree, int node, NumericVector abundances) {
  TreeCache tc = build_cache(tree);
  AbundanceLookup abund = preprocess_abundances_vec(abundances);
  TISResult r = compute_T_i_S_i_cached(tc, node, abund);
  // Match R's list(T_i, DF_S_i)
  return List::create(r.T, tis_to_dataframe(r));
}

static inline SIAResult calculate_S_i_a_star_cached(const TreeCache& tc,
                                                    int node,
                                                    const AbundanceLookup& abund) {
  SIAResult out;
  std::vector<BranchInfo> br = descendant_branches_with_x(tc, node);
  for (auto &b : br) b.x = b.branchLength; // star: x = Branch_length

  const int n = static_cast<int>(br.size());
  out.S.reserve(n);
  out.x.reserve(n);
  out.abundances.reserve(n);
  out.nodes.reserve(n);

  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int lhs, int rhs){
    return br[lhs].x < br[rhs].x;
  });

  std::vector<char> alive(n, 1);
  double remaining_sum = 0.0;
  for (const auto& b : br) remaining_sum += abund.get(b.endNode);

  int i = 0;
  int remaining_count = n;
  while (i < n) {
    double min_x = br[order[i]].x;

    // Preserve the original remaining-branch order in the returned abundance vectors.
    std::vector<double> av;
    std::vector<int> nd;
    av.reserve(remaining_count);
    nd.reserve(remaining_count);
    for (int k = 0; k < n; ++k) {
      if (!alive[k]) continue;
      av.push_back(abund.get(br[k].endNode));
      nd.push_back(br[k].endNode);
    }

    out.S.push_back(remaining_sum);
    out.x.push_back(min_x);
    out.abundances.push_back(std::move(av));
    out.nodes.push_back(std::move(nd));

    int j = i;
    while (j < n && approx_equal(br[order[j]].x, min_x)) {
      int idx = order[j];
      if (alive[idx]) {
        remaining_sum -= abund.get(br[idx].endNode);
        alive[idx] = 0;
        --remaining_count;
      }
      ++j;
    }
    i = j;
  }

  return out;
}

// [[Rcpp::export]]
List calculate_S_i_a_star(List tree, int node, NumericVector abundances) {
  TreeCache tc = build_cache(tree);
  AbundanceLookup abund = preprocess_abundances_vec(abundances);
  return sia_to_list(calculate_S_i_a_star_cached(tc, node, abund));
}

static inline SIAResult calculate_S_i_a_cached(const TreeCache& tc,
                                               int node,
                                               const AbundanceLookup& abund,
                                               int curr_ancestor,
                                               double h,
                                               double l_i) {
  SIAResult out;

  std::vector<BranchInfo> edges = descendant_branches_with_x(tc, (node == curr_ancestor ? node : curr_ancestor));

  if (node != curr_ancestor) {
    // Filter x > h
    edges.erase(std::remove_if(edges.begin(), edges.end(),
      [&](const BranchInfo& b){ return !(b.x > h + 1e-12); }), edges.end());

    // Correct branch lengths for edges starting before node: (x - Branch_length) < h
    for (auto &b : edges) {
      double startDist = b.x - b.branchLength;
      if (startDist < h - 1e-12) {
        b.branchLength = b.x - h;
      }
      b.x = b.x - h; // shift x to be measured from node
    }
  }

  out.S.reserve(edges.size());
  out.x.reserve(edges.size());
  out.abundances.reserve(edges.size());
  out.nodes.reserve(edges.size());

  double x_global = 0.0;
  while (!edges.empty() && x_global < l_i - 1e-12) {
    std::vector<int> active_idx;
    active_idx.reserve(edges.size());
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
      if (approx_equal(edges[i].x, edges[i].branchLength)) active_idx.push_back(i);
    }
    if (active_idx.empty()) break; // should not happen in well-formed rooted phylo

    double min_x_active = std::numeric_limits<double>::infinity();
    for (int idx : active_idx) min_x_active = std::min(min_x_active, edges[idx].x);

    std::vector<double> av;
    std::vector<int> nd;
    av.reserve(active_idx.size());
    nd.reserve(active_idx.size());
    double Ssum = 0.0;
    for (int idx : active_idx) {
      double a = abund.get(edges[idx].endNode);
      av.push_back(a);
      nd.push_back(edges[idx].endNode);
      Ssum += a;
    }

    out.S.push_back(Ssum);
    out.x.push_back(min_x_active + x_global);
    out.abundances.push_back(std::move(av));
    out.nodes.push_back(std::move(nd));

    double prev_x = min_x_active;
    x_global += prev_x;

    // indices to delete: ALL edges with x == min_x_active (R does not restrict to active)
    std::vector<char> del(edges.size(), 0);
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
      if (approx_equal(edges[i].x, min_x_active)) del[i] = 1;
    }

    // shift x for all edges
    for (auto &b : edges) b.x -= prev_x;

    // shorten Branch_length only for active edges
    for (int idx : active_idx) edges[idx].branchLength -= prev_x;

    std::vector<BranchInfo> kept;
    kept.reserve(edges.size());
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
      if (!del[i]) kept.push_back(edges[i]);
    }
    edges.swap(kept);
  }

  return out;
}

// [[Rcpp::export]]
List calculate_S_i_a(List tree, int node, NumericVector abundances, int curr_ancestor, double h, double l_i) {
  TreeCache tc = build_cache(tree);
  AbundanceLookup abund = preprocess_abundances_vec(abundances);
  return sia_to_list(calculate_S_i_a_cached(tc, node, abund, curr_ancestor, h, l_i));
}

static inline DJSeriesResult calculate_DJ_i_a_cached(const TreeCache& tc,
                                                     int node,
                                                     const AbundanceLookup& abund,
                                                     int curr_ancestor,
                                                     double h,
                                                     double l_i,
                                                     std::string index_letter,
                                                     int q,
                                                     bool individual) {
  index_letter = to_upper_copy(index_letter);

  SIAResult sres = calculate_S_i_a_cached(tc, node, abund, curr_ancestor, h, l_i);

  DJSeriesResult out;
  out.individual = individual;
  out.x = sres.x;
  const int n = static_cast<int>(sres.abundances.size());

  if (individual) {
    out.In.assign(n, 0.0);
    for (int k = 0; k < n; ++k) {
      const std::vector<double>& av = sres.abundances[k];
      int m = static_cast<int>(av.size());
      double S = sres.S[k];

      if (index_letter == "J") {
        if (m == 1) out.In[k] = 1.0;
        else {
          double log_m = std::log(static_cast<double>(m));
          double sum = 0.0;
          for (double a : av) {
            double p = a / S;
            sum += -p * (std::log(p) / log_m);
          }
          out.In[k] = sum;
        }
      } else { // D
        if (q == 1) {
          double sum = 0.0;
          for (double a : av) {
            double p = a / S;
            sum += -p * std::log(p);
          }
          out.In[k] = sum;
        } else { // q == 0
          out.In[k] = std::log(static_cast<double>(m));
        }
      }
    }
    return out;
  }

  out.D1.assign(n, 0.0);
  out.J1.assign(n, 0.0);
  out.D0.assign(n, 0.0);
  for (int k = 0; k < n; ++k) {
    const std::vector<double>& av = sres.abundances[k];
    int m = static_cast<int>(av.size());
    double S = sres.S[k];

    double d1 = 0.0;
    for (double a : av) {
      double p = a / S;
      d1 += -p * std::log(p);
    }
    out.D1[k] = d1;
    out.D0[k] = std::log(static_cast<double>(m));

    if (m == 1) out.J1[k] = 1.0;
    else {
      double log_m = std::log(static_cast<double>(m));
      double sum = 0.0;
      for (double a : av) {
        double p = a / S;
        sum += -p * (std::log(p) / log_m);
      }
      out.J1[k] = sum;
    }
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::RObject calculate_DJ_i_a(List tree, int node, NumericVector abundances,
                               int curr_ancestor, double h, double l_i,
                               std::string index_letter, int q, bool individual) {

  TreeCache tc = build_cache(tree);
  AbundanceLookup abund = preprocess_abundances_vec(abundances);
  DJSeriesResult r = calculate_DJ_i_a_cached(tc, node, abund, curr_ancestor, h, l_i,
                                             index_letter, q, individual);

  if (individual) {
    return DataFrame::create(_["In"] = r.In, _["x"] = r.x);
  }

  return List::create(
    _["1DN"] = DataFrame::create(_["D1"] = r.D1, _["x"] = r.x),
    _["1JN"] = DataFrame::create(_["J1"] = r.J1, _["x"] = r.x),
    _["0DN"] = DataFrame::create(_["D0"] = r.D0, _["x"] = r.x)
  );
}

static inline List calculate_EJM_i_a_from_sia(const SIAResult& sres,
                                              std::string index_letter,
                                              int q,
                                              bool individual) {
  index_letter = to_upper_copy(index_letter);
  const int n = static_cast<int>(sres.abundances.size());

  auto term = [](double a, double b, double log_base){
    if (a == 0.0 || b == 0.0) return 0.0;
    double p = a / b;
    return -p * (std::log(p) / log_base);
  };

  if (!individual) {
    std::vector<double> E(n), J(n), M(n);
    for (int k = 0; k < n; ++k) {
      const std::vector<double>& av = sres.abundances[k];
      int m = static_cast<int>(av.size());
      double S = sres.S[k];

      M[k] = std::log(static_cast<double>(m));

      double e = 0.0;
      for (double a : av) e += term(a, S, 1.0);
      E[k] = e;

      if (m == 1) J[k] = 1.0;
      else {
        double log_m = std::log(static_cast<double>(m));
        double j = 0.0;
        for (double a : av) j += term(a, S, log_m);
        J[k] = j;
      }
    }
    return List::create(
      _["1DN"] = DataFrame::create(_["E_i_a"] = E, _["x"] = sres.x),
      _["1JN"] = DataFrame::create(_["J_i_a"] = J, _["x"] = sres.x),
      _["0DN"] = DataFrame::create(_["M_i_a"] = M, _["x"] = sres.x)
    );
  }

  std::vector<double> In(n);
  for (int k = 0; k < n; ++k) {
    const std::vector<double>& av = sres.abundances[k];
    int m = static_cast<int>(av.size());
    double S = sres.S[k];

    if (index_letter == "D" && q == 0) In[k] = std::log(static_cast<double>(m));
    else if (index_letter == "D" && q == 1) {
      double e = 0.0;
      for (double a : av) e += term(a, S, 1.0);
      In[k] = e;
    } else { // J
      if (m == 1) In[k] = 1.0;
      else {
        double log_m = std::log(static_cast<double>(m));
        double j = 0.0;
        for (double a : av) j += term(a, S, log_m);
        In[k] = j;
      }
    }
  }
  return List::create(DataFrame::create(_["In_i_a"] = In, _["x"] = sres.x));
}

// [[Rcpp::export]]
List calculate_EJM_i_a(List tree, int node, NumericVector abundances,
                       int curr_ancestor, double h, double l_i,
                       std::string index_letter, int q, bool individual) {

  TreeCache tc = build_cache(tree);
  AbundanceLookup abund = preprocess_abundances_vec(abundances);
  SIAResult sres = calculate_S_i_a_cached(tc, node, abund, curr_ancestor, h, l_i);
  List res = calculate_EJM_i_a_from_sia(sres, index_letter, q, individual);
  if (individual) return Rcpp::as<DataFrame>(res[0]);
  return res;
}

static inline double max_x(const std::vector<double>& x) {
  double out = 0.0;
  for (double v : x) out = std::max(out, v);
  return out;
}

static inline double integral_step_vectors(const TISResult& S_i,
                                           const std::vector<double>& vx,
                                           const std::vector<double>& vv,
                                           double l_i) {
  double prev = 0.0;
  int rs = 0, ri = 0;
  double acc = 0.0;

  while (rs < static_cast<int>(S_i.x.size()) && ri < static_cast<int>(vx.size())) {
    double x_s = S_i.x[rs];
    double x_i = vx[ri];
    double v_s = S_i.S[rs];
    double v_i = vv[ri];

    if (x_s < x_i - 1e-12) {
      acc += v_s * v_i * (x_s - prev);
      prev = x_s;
      rs++;
    } else if (approx_equal(x_s, x_i)) {
      acc += v_s * v_i * (x_s - prev);
      prev = x_s;
      rs++; ri++;
    } else { // x_i < x_s
      double right = (x_i < l_i ? x_i : l_i);
      acc += v_s * v_i * (right - prev);
      prev = x_i;
      ri++;
    }
  }
  return acc;
}

static inline IntegralResult calculate_integral_cached(const TreeCache& tc,
                                                       int node,
                                                       int curr_ancestor,
                                                       const TISResult& S_i,
                                                       const AbundanceLookup& abund,
                                                       std::string index_letter,
                                                       int q,
                                                       bool individual) {
  IntegralResult out;
  out.individual = individual;

  double l_i = max_x(S_i.x);

  if (approx_equal(l_i, 0.0)) return out;

  double h = distance_cpp(tc, curr_ancestor, node);

  double d_parent;
  if (curr_ancestor == tc.root) {
    d_parent = l_i;
  } else if (curr_ancestor >= 0 && curr_ancestor < static_cast<int>(tc.len_to_parent.size())) {
    d_parent = tc.len_to_parent[curr_ancestor] + h;
  } else {
    d_parent = l_i;
  }

  // int_2
  if (l_i <= h + 1e-12) return out;

  double int_2 = 0.0;
  double prev_x = h;

  int start_row = -1;
  for (int i = 0; i < static_cast<int>(S_i.x.size()); ++i) {
    if (S_i.x[i] > h + 1e-12) { start_row = i; break; }
  }
  if (start_row == -1) return out;

  for (int r = start_row; r < static_cast<int>(S_i.x.size()); ++r) {
    double x_s = S_i.x[r];
    double value_s = S_i.S[r];
    if (x_s > d_parent + 1e-12) {
      int_2 += value_s * (d_parent - prev_x);
      break;
    } else {
      int_2 += value_s * (x_s - prev_x);
      prev_x = x_s;
    }
  }

  if (approx_equal(int_2, 0.0)) return out;

  DJSeriesResult sum_i_a = calculate_DJ_i_a_cached(tc, node, abund, curr_ancestor, h, l_i,
                                                   index_letter, q, individual);

  if (individual) {
    out.value = integral_step_vectors(S_i, sum_i_a.x, sum_i_a.In, l_i) * int_2;
  } else {
    out.values[0] = integral_step_vectors(S_i, sum_i_a.x, sum_i_a.D1, l_i) * int_2;
    out.values[1] = integral_step_vectors(S_i, sum_i_a.x, sum_i_a.J1, l_i) * int_2;
    out.values[2] = integral_step_vectors(S_i, sum_i_a.x, sum_i_a.D0, l_i) * int_2;
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::RObject calculate_integral(List tree, int node, int curr_ancestor, DataFrame S_i,
                                 NumericVector abundances, std::string index_letter,
                                 int q, bool individual) {

  TreeCache tc = build_cache(tree);
  AbundanceLookup abund = preprocess_abundances_vec(abundances);
  TISResult s = dataframe_to_tis(S_i);

  IntegralResult r = calculate_integral_cached(tc, node, curr_ancestor, s, abund,
                                               index_letter, q, individual);
  if (individual) return wrap(r.value);
  return NumericVector::create(r.values[0], r.values[1], r.values[2]);
}

// ---------- Helpers for high-level indices (long_star, node, all_indices) ----------

static inline double sum_edge_lengths(const List& tree) {
  NumericVector el = tree["edge.length"];
  double s = 0.0;
  for (double v : el) s += v;
  return s;
}

// Compute subtree-total abundances for every node, following abundance_phylo() in the R file.
// - If node_abundances is absent/false: tips get 1/Ntip, internal nodes get 0.
// - If node_abundances is a data.frame with columns (names, values): those are the *direct* node abundances;
//   we then convert to subtree totals.
static inline NumericVector abundance_phylo_cpp_cached(const List& tree,
                                                       const Rcpp::RObject& node_abundances,
                                                       const TreeCache& tc) {
  int max_node = tc.max_node;
  if (max_node <= 0) return NumericVector();

  CharacterVector tip_labels;
  CharacterVector node_labels;
  int Ntip = 0, Nnode = 0;

  if (tree.containsElementNamed("tip.label")) {
    tip_labels = tree["tip.label"];
    Ntip = tip_labels.size();
  }
  if (tree.containsElementNamed("node.label")) {
    node_labels = tree["node.label"];
    Nnode = node_labels.size();
  } else if (tree.containsElementNamed("Nnode")) {
    Nnode = Rcpp::as<int>(tree["Nnode"]);
    node_labels = CharacterVector(Nnode);
    for (int i = 0; i < Nnode; ++i) node_labels[i] = std::to_string(Ntip + 1 + i);
  }

  std::vector<double> direct(max_node + 1, 0.0); // 1-indexed
  if (is_null_or_false(node_abundances)) {
    // Default: leaves equally abundant, internal nodes zero
    if (Ntip <= 0) {
      // Fallback: infer tips as nodes with no children.
      int inferred = 0;
      for (int i = 1; i <= max_node; ++i) {
        if (i < static_cast<int>(tc.children.size()) && tc.children[i].empty()) inferred++;
      }
      Ntip = inferred;
    }
    double leaf = (Ntip > 0) ? (1.0 / static_cast<double>(Ntip)) : 0.0;
    for (int i = 1; i <= max_node; ++i) {
      direct[i] = (i <= Ntip) ? leaf : 0.0;
    }
  } else {
    // Map label -> value
    DataFrame df = Rcpp::as<DataFrame>(node_abundances);
    if (df.size() < 2) stop("node_abundances must have at least 2 columns: names and values.");
    CharacterVector nm = df[0];
    NumericVector val = df[1];
    std::unordered_map<std::string, double> label_value;
    label_value.reserve(nm.size() * 2 + 1);
    for (int i = 0; i < nm.size(); ++i) {
      Rcpp::String si = nm[i];
      if (si == NA_STRING) continue;
      label_value[std::string(si)] = val[i];
    }

    // For each node number, find its label (tip.label / node.label) and pull direct abundance.
    for (int i = 1; i <= max_node; ++i) {
      std::string label;
      if (i <= Ntip) {
        label = (Ntip > 0) ? std::string(tip_labels[i-1]) : std::to_string(i);
      } else {
        int idx = i - Ntip - 1;
        if (idx >= 0 && idx < node_labels.size()) label = std::string(node_labels[idx]);
        else label = std::to_string(i);
      }
      auto it = label_value.find(label);
      direct[i] = (it == label_value.end()) ? 0.0 : it->second;
    }
  }

  // Postorder recursion to get subtree totals
  std::vector<double> total(max_node + 1, std::numeric_limits<double>::quiet_NaN());
  std::function<double(int)> dfs = [&](int u) -> double {
    if (u <= 0 || u > max_node) return 0.0;
    double &memo = total[u];
    if (!std::isnan(memo)) return memo;

    double s = direct[u];
    if (u < static_cast<int>(tc.children.size())) {
      for (const auto &ch : tc.children[u]) s += dfs(ch.first);
    }
    memo = s;
    return memo;
  };

  if (tc.root != NA_INTEGER) dfs(tc.root);

  NumericVector out(max_node);
  CharacterVector nms(max_node);
  for (int i = 1; i <= max_node; ++i) {
    nms[i-1] = std::to_string(i);
    double v = total[i];
    out[i-1] = std::isnan(v) ? 0.0 : v;
  }
  out.names() = nms;
  return out;
}

static inline NumericVector abundance_phylo_cpp(const List& tree, const Rcpp::RObject& node_abundances) {
  TreeCache tc = build_cache(tree);
  return abundance_phylo_cpp_cached(tree, node_abundances, tc);
}

// R-equivalent: calculate_DJ_i()
static inline Rcpp::RObject calculate_DJ_i_cpp_impl_cached(const TreeCache& tc,
                                                           int node,
                                                           const AbundanceLookup& abund,
                                                           std::string index_letter,
                                                           int q,
                                                           bool individual) {
  index_letter = to_upper_copy(index_letter);

  // Ancestors of node (excluding node), then node itself (order doesn't matter for summation)
  std::vector<int> ancestors;
  int cur = node;
  while (cur >= 0 && cur < static_cast<int>(tc.parent.size()) && tc.parent[cur] != NA_INTEGER) {
    cur = tc.parent[cur];
    ancestors.push_back(cur);
  }
  ancestors.push_back(node);

  TISResult S_i = compute_T_i_S_i_cached(tc, node, abund);
  double T_i = S_i.T;

  if (approx_equal(T_i, 0.0)) {
    if (individual) return Rcpp::wrap(0.0);
    return Rcpp::NumericVector::create(0.0, 0.0, 0.0);
  }

  if (individual) {
    double acc = 0.0;
    for (int a : ancestors) {
      IntegralResult tmp = calculate_integral_cached(tc, node, a, S_i, abund, index_letter, q, true);
      acc += tmp.value;
    }
    return Rcpp::wrap(acc / T_i);
  }

  std::array<double, 3> acc{{0.0, 0.0, 0.0}};
  for (int a : ancestors) {
    IntegralResult tmp = calculate_integral_cached(tc, node, a, S_i, abund, index_letter, q, false);
    acc[0] += tmp.values[0];
    acc[1] += tmp.values[1];
    acc[2] += tmp.values[2];
  }

  return Rcpp::NumericVector::create(acc[0] / T_i, acc[1] / T_i, acc[2] / T_i);
}

static inline Rcpp::RObject calculate_DJ_i_cpp_impl(const List& tree,
                                                    int node,
                                                    const NumericVector& abundances,
                                                    std::string index_letter,
                                                    int q,
                                                    bool individual) {
  TreeCache tc = build_cache(tree);
  AbundanceLookup abund = preprocess_abundances_vec(abundances);
  return calculate_DJ_i_cpp_impl_cached(tc, node, abund, index_letter, q, individual);
}

static inline Rcpp::RObject node_cpp_impl(const List& tree, const Rcpp::RObject& node_abundances,
                                          std::string index_letter, int q, bool individual) {

  index_letter = to_upper_copy(index_letter);

  CharacterVector tip_labels = tree.containsElementNamed("tip.label") ? CharacterVector(tree["tip.label"]) : CharacterVector();
  int Ntip = tip_labels.size();
  int Nnode = tree.containsElementNamed("Nnode") ? Rcpp::as<int>(tree["Nnode"]) : 0;

  // Linear tree (only one tip)
  if (Ntip == 1) {
    if (individual) return Rcpp::wrap(1.0);
    return List::create(_["D1N"] = 1.0, _["J1N"] = 1.0, _["D0N"] = 1.0);
  }

  TreeCache tc = build_cache(tree);
  NumericVector abund_vec = abundance_phylo_cpp_cached(tree, node_abundances, tc);
  AbundanceLookup abund = preprocess_abundances_vec(abund_vec);

  IntegerMatrix edge = tree["edge"];
  NumericVector el = tree["edge.length"];

  // T = sum_{edges} a_child * length(edge)
  double T = 0.0;
  for (int i = 0; i < edge.nrow(); ++i) {
    int child = edge(i,1);
    T += abund.get(child) * el[i];
  }
  if (approx_equal(T, 0.0)) {
    if (individual) return Rcpp::wrap(0.0);
    return List::create(_["D1N"] = 1.0, _["J1N"] = 0.0, _["D0N"] = 1.0);
  }

  // Internal node numbers: (Ntip+1) .. (Ntip+Nnode)
  std::vector<int> nodes;
  nodes.reserve(std::max(0, Nnode));
  for (int v = Ntip + 1; v <= Ntip + Nnode; ++v) nodes.push_back(v);

  if (!individual) {
    double sumD1 = 0.0, sumJ1 = 0.0, sumD0 = 0.0;
    for (int v : nodes) {
      Rcpp::RObject tmp = calculate_DJ_i_cpp_impl_cached(tc, v, abund, "D", 1, false);
      NumericVector dj = Rcpp::as<NumericVector>(tmp);
      if (dj.size() == 3) {
        sumD1 += dj[0];
        sumJ1 += dj[1];
        sumD0 += dj[2];
      }
    }
    double D1N_log = sumD1 / T;
    double J1N     = sumJ1 / T;
    double D0N_log = sumD0 / T;

    return List::create(_["D1N"] = std::exp(D1N_log),
                        _["J1N"] = J1N,
                        _["D0N"] = std::exp(D0N_log));
  }

  double acc = 0.0;
  for (int v : nodes) {
    Rcpp::RObject tmp = calculate_DJ_i_cpp_impl_cached(tc, v, abund, index_letter, q, true);
    acc += Rcpp::as<double>(tmp);
  }
  double mean = acc / T;
  if (index_letter == "J") return Rcpp::wrap(mean);
  return Rcpp::wrap(std::exp(mean));
}

static inline Rcpp::RObject long_star_cpp_impl(const List& tree, const Rcpp::RObject& node_abundances,
                                               std::string mean_type, std::string index_letter,
                                               int q, bool individual) {

  mean_type = to_upper_copy(mean_type);
  index_letter = to_upper_copy(index_letter);

  CharacterVector tip_labels = tree.containsElementNamed("tip.label") ? CharacterVector(tree["tip.label"]) : CharacterVector();
  int Ntip = tip_labels.size();
  int Nnode = tree.containsElementNamed("Nnode") ? Rcpp::as<int>(tree["Nnode"]) : 0;

  // Linear tree (only one tip)
  if (Ntip == 1 && mean_type == "LONGITUDINAL") {
    if (individual) return Rcpp::wrap(1.0);
    return List::create(_["D0L"] = 1.0, _["D1L"] = 1.0, _["J1L"] = 1.0);
  }
  if (Ntip == 1 && mean_type == "STAR" && is_null_or_false(node_abundances)) {
    // Match R's special-case shortcut when no abundances are provided
    if (individual) return Rcpp::wrap(1.0);
    return List::create(_["D0S"] = static_cast<double>(Nnode),
                        _["D1S"] = static_cast<double>(Nnode),
                        _["J1S"] = 1.0);
  }

  TreeCache tc = build_cache(tree);
  int root = tc.root;
  if (root == NA_INTEGER) stop("Could not determine root.");

  NumericVector abund_vec = abundance_phylo_cpp_cached(tree, node_abundances, tc);
  AbundanceLookup abund = preprocess_abundances_vec(abund_vec);

  SIAResult sres;
  if (mean_type == "LONGITUDINAL") {
    double L = sum_edge_lengths(tree);
    sres = calculate_S_i_a_cached(tc, root, abund, root, 0.0, L);
  } else if (mean_type == "STAR") {
    sres = calculate_S_i_a_star_cached(tc, root, abund);
  } else {
    stop("mean_type must be 'STAR' or 'LONGITUDINAL'.");
  }

  int n = static_cast<int>(sres.x.size());
  if (n == 0) {
    if (individual) return Rcpp::wrap(0.0);
    if (mean_type == "STAR") return List::create(_["D0S"]=1.0, _["D1S"]=1.0, _["J1S"]=0.0);
    return List::create(_["D0L"]=1.0, _["D1L"]=1.0, _["J1L"]=0.0);
  }

  std::vector<double> region(n);
  region[0] = sres.x[0];
  for (int i = 1; i < n; ++i) region[i] = sres.x[i] - sres.x[i-1];

  double T_S_sum = 0.0;
  for (int i = 0; i < n; ++i) T_S_sum += sres.S[i] * region[i];
  if (approx_equal(T_S_sum, 0.0)) {
    if (individual) return Rcpp::wrap(0.0);
    if (mean_type == "STAR") return List::create(_["D0S"]=1.0, _["D1S"]=1.0, _["J1S"]=0.0);
    return List::create(_["D0L"]=1.0, _["D1L"]=1.0, _["J1L"]=0.0);
  }

  auto contrib_entropy = [](const std::vector<double>& av, double Ssum) {
    // sum(-a * log(a/S))
    double out = 0.0;
    for (double a : av) {
      if (a <= 0.0 || Ssum <= 0.0) continue;
      double p = a / Ssum;
      out += -a * std::log(p);
    }
    return out;
  };

  auto contrib_J = [](const std::vector<double>& av, double Ssum) {
    // sum(-a * log(a/S, base=m)) = sum(-a * log(a/S)/log(m))
    int m = static_cast<int>(av.size());
    if (m <= 1) return Ssum; // matches R: 1 * S * (region length later)
    double denom = std::log(static_cast<double>(m));
    double out = 0.0;
    for (double a : av) {
      if (a <= 0.0 || Ssum <= 0.0) continue;
      double p = a / Ssum;
      out += -a * (std::log(p) / denom);
    }
    return out;
  };

  if (individual) {
    double acc = 0.0;
    for (int i = 0; i < n; ++i) {
      const std::vector<double>& av = sres.abundances[i];
      int m = static_cast<int>(av.size());
      double Ssum = sres.S[i];
      double r = region[i];

      if (index_letter == "J") {
        double c = (m <= 1) ? Ssum : contrib_J(av, Ssum);
        acc += c * r;
      } else {
        if (q == 1) {
          acc += contrib_entropy(av, Ssum) * r;
        } else { // q==0
          acc += Ssum * r * std::log(static_cast<double>(m));
        }
      }
    }
    double mean = acc / T_S_sum;
    if (index_letter == "J") return Rcpp::wrap(mean);
    return Rcpp::wrap(std::exp(mean));
  }

  double acc_h0 = 0.0, acc_h1 = 0.0, acc_j1 = 0.0;
  for (int i = 0; i < n; ++i) {
    const std::vector<double>& av = sres.abundances[i];
    int m = static_cast<int>(av.size());
    double Ssum = sres.S[i];
    double r = region[i];

    acc_h1 += contrib_entropy(av, Ssum) * r;
    acc_h0 += Ssum * r * std::log(static_cast<double>(m));
    acc_j1 += ((m <= 1) ? Ssum : contrib_J(av, Ssum)) * r;
  }

  double D0 = std::exp(acc_h0 / T_S_sum);
  double D1 = std::exp(acc_h1 / T_S_sum);
  double J1 = acc_j1 / T_S_sum;

  if (mean_type == "STAR") {
    return List::create(_["D0S"] = D0, _["D1S"] = D1, _["J1S"] = J1);
  }
  return List::create(_["D0L"] = D0, _["D1L"] = D1, _["J1L"] = J1);
}

// Convert file/path/phylo -> phylo list using the user's R helper read_convert() when needed.
static inline List ensure_phylo(SEXP file) {
  Rcpp::RObject obj(file);
  if (obj.inherits("phylo")) return Rcpp::as<List>(obj);

  stop("Input must be a phylo object at the C++ level. Use the exported R wrapper, which calls read_convert() first.");
}

// [[Rcpp::export]]
Rcpp::RObject calculate_DJ_i(List tree, int node, NumericVector abundances,
                             std::string index_letter = "D", int q = 1, bool individual = false) {
  return calculate_DJ_i_cpp_impl(tree, node, abundances, index_letter, q, individual);
}

// [[Rcpp::export]]
Rcpp::RObject node_cpp(SEXP file, Rcpp::RObject node_abundances = R_NilValue,
                   std::string index_letter = "D", int q = 1, bool individual = false) {
  List tree = ensure_phylo(file);
  return node_cpp_impl(tree, node_abundances, index_letter, q, individual);
}

// [[Rcpp::export]]
Rcpp::RObject long_star_cpp(SEXP file, Rcpp::RObject node_abundances = R_NilValue,
                        std::string mean_type = "Star", std::string index_letter = "D",
                        int q = 1, bool individual = false) {
  List tree = ensure_phylo(file);
  return long_star_cpp_impl(tree, node_abundances, mean_type, index_letter, q, individual);
}

// [[Rcpp::export]]
List all_indices_cpp(SEXP file, Rcpp::RObject node_abundances = R_NilValue) {
  List tree = ensure_phylo(file);

  List node_res = Rcpp::as<List>(node_cpp_impl(tree, node_abundances, "D", 1, false));
  List star_res = Rcpp::as<List>(long_star_cpp_impl(tree, node_abundances, "STAR", "D", 0, false));
  List long_res = Rcpp::as<List>(long_star_cpp_impl(tree, node_abundances, "LONGITUDINAL", "D", 0, false));

  return List::create(
    _["D0N"] = node_res["D0N"],
    _["D1N"] = node_res["D1N"],
    _["J1N"] = node_res["J1N"],
    _["D0S"] = star_res["D0S"],
    _["D1S"] = star_res["D1S"],
    _["J1S"] = star_res["J1S"],
    _["D0L"] = long_res["D0L"],
    _["D1L"] = long_res["D1L"],
    _["J1L"] = long_res["J1L"]
  );
}
