#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;

struct Edge {
    int to;
    double res_cap;
    int rev_idx;
};

class FlowNetwork {
  public:
    int N;
    vector<vector<Edge>> graph;

    FlowNetwork(int n) : N{n}, graph{(size_t) n} {}

    void addEdge(int u, int v, double cap) {
        graph[u].emplace_back(Edge{v, cap, (int) graph[v].size()});
        graph[v].emplace_back(Edge{u, 0.0, (int) graph[u].size() - 1});
    }
};

double dinic(int s, int t, FlowNetwork& fn) {
    int N = fn.N;
    vector<int> level(N), ptr(N);

    const double INF = 1e18;

    auto bfs = [&]() -> bool {
        fill(level.begin(), level.end(), -1);
        queue<int> q;
        q.push(s);
        level[s] = 0;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (auto& e : fn.graph[u]) {
                if (level[e.to] == -1 && e.res_cap > 1e-9) {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                }
            }
        }
        return level[t] != -1;
    };

    function<double(int, double)> dfs = [&](int u, double pushed) -> double {
        if (pushed < 1e-9)
            return 0.0;
        if (u == t)
            return pushed;

        for (int& cid = ptr[u]; cid < (int) fn.graph[u].size(); cid++) {
            Edge& e = fn.graph[u][cid];

            if (level[e.to] != level[u] + 1 || e.res_cap < 1e-9)
                continue;

            double tr = dfs(e.to, min(pushed, e.res_cap));
            if (tr < 1e-9)
                continue;

            e.res_cap -= tr;
            fn.graph[e.to][e.rev_idx].res_cap += tr;
            return tr;
        }
        return 0.0;
    };

    double flow = 0.0;

    while (bfs()) {
        fill(ptr.begin(), ptr.end(), 0);
        while (double pushed = dfs(s, INF)) {
            flow += pushed;
        }
    }

    return flow;
}

double edmondsKarp(int s, int t, FlowNetwork& fn) {
    double maxflow = 0.0;
    while (true) {
        vector<int> parent(fn.N, -1);
        vector<int> parentEdge(fn.N, -1);
        queue<int> q;
        q.push(s);
        parent[s] = s;

        while (!q.empty() && parent[t] == -1) {
            int u = q.front();
            q.pop();
            for (int i = 0; i < (int) fn.graph[u].size(); i++) {
                Edge& e = fn.graph[u][i];
                if (parent[e.to] == -1 && e.res_cap > 1e-9) {
                    parent[e.to] = u;
                    parentEdge[e.to] = i;
                    q.push(e.to);
                }
            }
        }
        if (parent[t] == -1)
            break;

        double flow = 1e18;
        for (int v = t; v != s; v = parent[v]) {
            Edge& e = fn.graph[parent[v]][parentEdge[v]];
            flow = min(flow, e.res_cap);
        }
        for (int v = t; v != s; v = parent[v]) {
            Edge& e = fn.graph[parent[v]][parentEdge[v]];
            e.res_cap -= flow;
            fn.graph[v][e.rev_idx].res_cap += flow;
        }
        maxflow += flow;
    }
    return maxflow;
}

vector<bool> minCutReachable(int s, FlowNetwork& fn) {
    vector<bool> visited(fn.N, false);
    queue<int> q;
    q.push(s);
    visited[s] = true;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (auto& e : fn.graph[u]) {
            if (!visited[e.to] && e.res_cap > 1e-9) {
                visited[e.to] = true;
                q.push(e.to);
            }
        }
    }
    return visited;
}

// triangle enumeration
vector<array<int, 3>> getTriangles(int n, const vector<vector<int>>& adj) {
    vector<array<int, 3>> triangles;
    for (int u = 0; u < n; u++) {
        for (int v : adj[u])
            if (v > u) {
                int i = 0, j = 0;
                int su = adj[u].size(), sv = adj[v].size();
                while (i < su && j < sv) {
                    if (adj[u][i] == adj[v][j]) {
                        int w = adj[u][i];
                        if (w > v)
                            triangles.push_back({u, v, w});
                        i++;
                        j++;
                    } else if (adj[u][i] < adj[v][j])
                        i++;
                    else
                        j++;
                }
            }
    }
    return triangles;
}

//(k, Psi) core decomposition
vector<int> cliqueCoreDecomposition(int n,
                                    const vector<vector<int>>& /*adj*/,
                                    const vector<array<int, 3>>& triangles,
                                    double& best_residual_density) {
    int T = triangles.size();

    vector<int> deg(n, 0);
    vector<vector<int>> tri_of(n);
    for (int t = 0; t < T; t++) {
        for (int x : {triangles[t][0], triangles[t][1], triangles[t][2]}) {
            deg[x]++;
            tri_of[x].push_back(t);
        }
    }

    int max_deg = 0;
    for (int v = 0; v < n; v++)
        max_deg = max(max_deg, deg[v]);

    vector<vector<int>> bin(max_deg + 1);
    vector<int> pos(n, -1);
    for (int v = 0; v < n; v++) {
        pos[v] = bin[deg[v]].size();
        bin[deg[v]].push_back(v);
    }

    vector<bool> removed_v(n, false);
    vector<bool> dead_t(T, false);
    vector<int> core(n, 0);

    auto moveBin = [&](int v, int old_d, int new_d) {
        int p = pos[v];
        int last = bin[old_d].back();
        bin[old_d][p] = last;
        pos[last] = p;
        bin[old_d].pop_back();
        if ((int) bin.size() <= new_d)
            bin.resize(new_d + 1);
        pos[v] = bin[new_d].size();
        bin[new_d].push_back(v);
    };

    int remaining_v = n;
    int remaining_t = T;
    best_residual_density = (remaining_v > 0)
                                ? (double) remaining_t / remaining_v
                                : 0.0;

    int processed = 0;
    int scan = 0;
    while (processed < n) {
        while (scan < (int) bin.size() && bin[scan].empty())
            scan++;
        if (scan >= (int) bin.size())
            break;

        int v = bin[scan].back();
        bin[scan].pop_back();
        pos[v] = -1;

        int deg_v = deg[v];
        core[v] = deg_v;
        removed_v[v] = true;
        processed++;
        remaining_v--;

        for (int tid : tri_of[v]) {
            if (dead_t[tid])
                continue;
            dead_t[tid] = true;
            remaining_t--;
            for (int x : {triangles[tid][0], triangles[tid][1], triangles[tid][2]}) {
                if (x == v || removed_v[x])
                    continue;
                // Only decrement strictly-larger neighbours (Alg.3 line 8).
                if (deg[x] > deg_v) {
                    int old_d = deg[x];
                    int new_d = old_d - 1;
                    deg[x] = new_d;
                    moveBin(x, old_d, new_d);
                    if (new_d < scan)
                        scan = new_d;
                }
            }
        }

        if (remaining_v > 0) {
            double d = (double) remaining_t / remaining_v;
            if (d > best_residual_density)
                best_residual_density = d;
        }
    }

    return core;
}

// Connected components of the subgraph induced by a vertex mask.
vector<vector<int>> connectedComponents(const vector<vector<int>>& adj,
                                        const vector<bool>& in_sub) {
    int n = adj.size();
    vector<bool> visited(n, false);
    vector<vector<int>> comps;
    for (int v = 0; v < n; v++) {
        if (!in_sub[v] || visited[v])
            continue;
        vector<int> comp;
        queue<int> q;
        q.push(v);
        visited[v] = true;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            comp.push_back(u);
            for (int w : adj[u]) {
                if (in_sub[w] && !visited[w]) {
                    visited[w] = true;
                    q.push(w);
                }
            }
        }
        comps.push_back(std::move(comp));
    }
    return comps;
}

// Number of triangles fully contained in the vertex mask.
int countTrianglesInMask(const vector<bool>& mask,
                         const vector<array<int, 3>>& triangles) {
    int cnt = 0;
    for (auto& t : triangles)
        if (mask[t[0]] && mask[t[1]] && mask[t[2]])
            cnt++;
    return cnt;
}

double triangleDensityOfMask(const vector<bool>& mask,
                             const vector<array<int, 3>>& triangles,
                             int sz) {
    if (sz <= 0)
        return 0.0;
    return (double) countTrianglesInMask(mask, triangles) / sz;
}

// Build the flow network restricted to a subgraph C.
FlowNetwork buildFlowNetworkOnSubgraph(
    const vector<int>& comp,
    const vector<int>& g2l,
    const vector<vector<int>>& adj,
    const vector<array<int, 3>>& triangles,
    double alpha,
    int& s_out, int& t_out) {
    static const double INF = 1e18;
    int nc = comp.size();

    // collect edges within C and assign local edge IDs
    map<pair<int, int>, int> edge_id;
    int eid = 0;
    for (int gu : comp) {
        for (int gv : adj[gu]) {
            if (gu < gv && g2l[gv] != -1) {
                edge_id[{gu, gv}] = eid++;
            }
        }
    }

    // collect triangles within C and per-vertex (in-component) triangle degree
    vector<array<int, 3>> local_tris;
    local_tris.reserve(triangles.size());
    vector<int> tri_deg(nc, 0);
    for (auto& tri : triangles) {
        int a = tri[0], b = tri[1], c = tri[2];
        if (g2l[a] != -1 && g2l[b] != -1 && g2l[c] != -1) {
            local_tris.push_back(tri);
            tri_deg[g2l[a]]++;
            tri_deg[g2l[b]]++;
            tri_deg[g2l[c]]++;
        }
    }

    int N = nc + eid;
    int s = N, t = N + 1;
    s_out = s;
    t_out = t;
    FlowNetwork fn(N + 2);

    // s -> v with cap = tri_deg(v) within C
    for (int i = 0; i < nc; i++)
        fn.addEdge(s, i, tri_deg[i]);
    // v -> t with cap = |V_Psi| * alpha = 3 * alpha
    for (int i = 0; i < nc; i++)
        fn.addEdge(i, t, 3.0 * alpha);
    // edge-node -> its two vertices with cap +inf
    for (auto& kv : edge_id) {
        int gu = kv.first.first, gv = kv.first.second;
        int e_node = nc + kv.second;
        fn.addEdge(e_node, g2l[gu], INF);
        fn.addEdge(e_node, g2l[gv], INF);
    }
    // for each triangle (a,b,c): a -> e_{bc}, b -> e_{ac}, c -> e_{ab}, cap 1
    auto mmp = [](int x, int y) {
        return pair<int, int>{min(x, y), max(x, y)};
    };
    for (auto& tri : local_tris) {
        int a = tri[0], b = tri[1], c = tri[2];
        int e_ab = nc + edge_id[mmp(a, b)];
        int e_bc = nc + edge_id[mmp(b, c)];
        int e_ac = nc + edge_id[mmp(a, c)];
        fn.addEdge(g2l[a], e_bc, 1);
        fn.addEdge(g2l[b], e_ac, 1);
        fn.addEdge(g2l[c], e_ab, 1);
    }
    return fn;
}

// Convenience: rebuild g2l for a (possibly shrunken) component.
void refreshG2L(const vector<int>& comp, vector<int>& g2l) {
    fill(g2l.begin(), g2l.end(), -1);
    for (int i = 0; i < (int) comp.size(); i++)
        g2l[comp[i]] = i;
}

// Algorithm 4: CoreExact
vector<int> findDensestSubgraphCoreExact(const vector<vector<int>>& adj) {
    int n = adj.size();
    vector<array<int, 3>> triangles = getTriangles(n, adj);

    // core decomposition (Algorithm 3) ----
    // Side effect: rho_prime = best residual density seen during peeling.
    double rho_prime = 0.0;
    vector<int> core = cliqueCoreDecomposition(n, adj, triangles, rho_prime);

    int kmax = 0;
    for (int v = 0; v < n; v++)
        kmax = max(kmax, core[v]);

    if (kmax == 0)
        return {}; // no triangles at all

    // locate (k'', Psi)-core via Pruning1 + Pruning2 ----
    // Pruning1: lower bound on rho_opt is rho'  =>  k' = ceil(rho')
    auto ceilEps = [](double x) {
        return (int) ceil(x - 1e-9);
    };
    int k_prime = max(0, ceilEps(rho_prime));

    // (k', Psi)-core
    vector<bool> in_kprime(n, false);
    for (int v = 0; v < n; v++)
        if (core[v] >= k_prime)
            in_kprime[v] = true;

    // Pruning2: tighten lower bound to max density across components
    auto comps_kprime = connectedComponents(adj, in_kprime);
    double rho_dprime = rho_prime;
    vector<int> seed_D;
    double seed_density = -1.0;
    for (auto& comp : comps_kprime) {
        vector<bool> mask(n, false);
        for (int v : comp)
            mask[v] = true;
        double d = triangleDensityOfMask(mask, triangles, comp.size());
        if (d > rho_dprime)
            rho_dprime = d;
        if (d > seed_density) {
            seed_density = d;
            seed_D = comp;
        }
    }
    int k_dprime = max(k_prime, ceilEps(rho_dprime));

    // (k'', Psi)-core, then its connected components -> the set C
    vector<bool> in_kdprime(n, false);
    for (int v = 0; v < n; v++)
        if (core[v] >= k_dprime)
            in_kdprime[v] = true;
    auto C = connectedComponents(adj, in_kdprime);

    // initialise
    double l = rho_dprime;
    double u = (double) kmax;
    vector<int> D = seed_D;
    double best_density = seed_density;

    vector<int> g2l(n, -1);

    // loop over connected components of (k'', Psi)-core
    for (auto comp : C) { // copy: we will shrink it
        // if our current lower bound l exceeds k'', we can
        // safely restrict to the (ceil(l), Psi)-core inside this comp.
        if (l > (double) k_dprime + 1e-12) {
            int kc = ceilEps(l);
            vector<int> filt;
            for (int v : comp)
                if (core[v] >= kc)
                    filt.push_back(v);
            comp = std::move(filt);
        }
        if (comp.empty())
            continue;

        refreshG2L(comp, g2l);

        // feasibility probe at alpha = l
        int sx, tx;
        FlowNetwork fn0 = buildFlowNetworkOnSubgraph(
            comp, g2l, adj, triangles, l, sx, tx);
        dinic(sx, tx, fn0);
        vector<bool> reach0 = minCutReachable(sx, fn0);
        bool nonempty = false;
        for (int i = 0; i < (int) comp.size(); i++)
            if (reach0[i]) {
                nonempty = true;
                break;
            }
        if (!nonempty)
            continue; // skip this component

        // Initial best subgraph from this component = the feasible side at l
        vector<int> U;
        for (int i = 0; i < (int) comp.size(); i++)
            if (reach0[i])
                U.push_back(comp[i]);

        // Lines 10-19: binary search inside this component
        int Vc = comp.size();
        while (u - l >= 1.0 / ((double) Vc * (Vc - 1) + 1e-12)) {
            double alpha = (l + u) / 2.0;

            refreshG2L(comp, g2l);
            int sa, ta;
            FlowNetwork fnA = buildFlowNetworkOnSubgraph(
                comp, g2l, adj, triangles, alpha, sa, ta);
            dinic(sa, ta, fnA);
            vector<bool> reachA = minCutReachable(sa, fnA);

            bool any = false;
            vector<int> U_new;
            for (int i = 0; i < (int) comp.size(); i++) {
                if (reachA[i]) {
                    any = true;
                    U_new.push_back(comp[i]);
                }
            }

            if (!any) {
                // No subgraph of density >= alpha in C; tighten upper bound.
                u = alpha;
            } else {
                // Found one; tighten lower bound and (Pruning1 mid-search)
                // shrink C to the (ceil(alpha), Psi)-core if alpha pushed
                // past the previous integer ceiling.
                int old_ceil = ceilEps(l);
                if (alpha > (double) old_ceil + 1e-12) {
                    int kc = ceilEps(alpha);
                    vector<int> filt;
                    for (int v : comp)
                        if (core[v] >= kc)
                            filt.push_back(v);
                    comp = std::move(filt);
                    Vc = comp.size();
                    if (Vc < 2) {
                        l = alpha;
                        U = U_new;
                        break;
                    }
                }
                l = alpha;
                U = std::move(U_new);
            }
        }

        // keep the best subgraph across all components
        if (!U.empty()) {
            vector<bool> umask(n, false);
            for (int v : U)
                umask[v] = true;
            double dU = triangleDensityOfMask(umask, triangles, U.size());
            if (dU > best_density) {
                best_density = dU;
                D = U;
            }
        }
    }

    return D;
}
vector<vector<int>> buildGraphFromInput(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file " << filename << endl;
        exit(1);
    }

    set<pair<int, int>> unique_edges;
    int u, v;
    int maxV = 0;
    string line;

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        stringstream ss(line);
        if (ss >> u >> v) {
            if (u == v)
                continue;
            if (u > v)
                swap(u, v); // Normalize to ensure undirected uniqueness
            unique_edges.insert({u, v});
            maxV = max({u, v, maxV});
        }
    }
    int numVertices = maxV + 1;
    cout << "Graph info: " << numVertices << " vertices, " << unique_edges.size() << " unique edges\n";

    vector<vector<int>> graph(numVertices);

    for (const auto& edge : unique_edges) {
        int u = edge.first, v = edge.second;
        graph[u].push_back(v);
        graph[v].push_back(u);
    }

    // sort adj lists once for triangle intersection finding
    for (int i = 0; i < numVertices; i++)
        sort(graph[i].begin(), graph[i].end());

    return graph;
}

void display_memory_usage(ofstream& out) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    // cout << "Peak memory usage: " << usage.ru_maxrss / (1024 * 1024) << "MB" << endl; // for macos, unit is in bytes
    out << "Peak memory usage: " << usage.ru_maxrss / 1024 << "MB" << endl; // for linux, unit is in kb, uncomment this line
}

int main(int argc, char* argv[]) {

    string filename = "email-Enron.txt";
    if (argc > 1)
        filename = argv[1];

    vector<vector<int>> adj = buildGraphFromInput(filename);

    auto start = chrono::high_resolution_clock::now();
    vector<int> DS = findDensestSubgraphCoreExact(adj);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> diff = end - start;
    cout << "Densest Subgraph found in " << diff.count() << " seconds." << endl;

    // get only neighbours in DS
    unordered_set<int> inDS(DS.begin(), DS.end());

    // output
    string outfile;

    if (argc > 2) {
        outfile = argv[2];
    }

    ofstream out(outfile);

    int triangle_count = 0;

    display_memory_usage(out);

    for (const int& u : DS) {
        for (int v : adj[u])
            if (v > u && inDS.count(v)) {
                int i = 0, j = 0;
                const auto& au = adj[u];
                const auto& av = adj[v];

                while (i < au.size() && j < av.size()) {
                    if (au[i] == av[j]) {
                        int w = au[i];
                        if (w > v && inDS.count(w)) {
                            triangle_count++;
                        }
                        i++;
                        j++;
                    } else if (au[i] < av[j])
                        i++;
                    else
                        j++;
                }
            }
    }

    double density = DS.empty() ? 0.0 : (double) triangle_count / DS.size();

    out << "Clique density of densest subgraph: " << density << endl;
    out << "Densest Subgraph found in " << diff.count() << " seconds.\n";
    out << "Nodes which are part of densest subgraph: \n";
    for (const int& node : DS) {
        out << node << " ";
    }
    out << "\n";
    out.close();
    cout << "Results written to " << outfile << endl;

    return 0;
}