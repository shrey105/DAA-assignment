#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <sys/resource.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

struct Edge {
    int to;
    double res_cap;
    int rev_idx;

    Edge(int to, double res_cap, int rev_idx) : to(to), res_cap(res_cap), rev_idx(rev_idx) {}
};

class FlowNetwork {
  public:
    int N;
    vector<vector<Edge>> graph;

    FlowNetwork(int n) : N{n}, graph{(size_t) n} {}

    void addEdge(int u, int v, double cap) {
        graph[u].emplace_back(v, cap, (int) graph[v].size());
        graph[v].emplace_back(u, 0.0, (int) graph[u].size() - 1);
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
        vector<int> parent(fn.N, -1);     // node
        vector<int> parentEdge(fn.N, -1); // edge index

        queue<int> q{};
        q.push(s);
        parent[s] = s;

        // BFS for augmenting path
        while (!q.empty() && parent[t] == -1) {
            int u = q.front();
            q.pop();

            for (int i = 0; i < fn.graph[u].size(); i++) {
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

        // bottleneck edge
        double flow = 1e18;
        int v = t;
        while (v != s) {
            int u = parent[v];
            Edge& e = fn.graph[u][parentEdge[v]];
            flow = min(flow, e.res_cap);
            v = u;
        }

        v = t;
        while (v != s) {
            int u = parent[v];
            Edge& e = fn.graph[u][parentEdge[v]];
            e.res_cap -= flow;
            fn.graph[v][e.rev_idx].res_cap += flow;
            v = u;
        }

        maxflow += flow;
    }
    return maxflow;
}

vector<array<int, 3>> getTriangles(int n, const vector<vector<int>>& adj) {
    vector<array<int, 3>> triangles{};

    for (int u = 0; u < n; u++) {
        for (const int& v : adj[u])
            if (v > u) {
                int i = 0, j = 0;
                int su = adj[u].size(), sv = adj[v].size();
                while (i < su && j < sv) {
                    if (adj[u][i] == adj[v][j]) {
                        int w = adj[u][i];
                        if (w > v) {
                            triangles.push_back({u, v, w});
                        }
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

inline long long edgeKey(int u, int v) {
    if (u > v)
        swap(u, v);
    return ((long long) u << 32) | v;
}

inline pair<int, int> decodeEdge(long long key) {
    int u = key >> 32;
    int v = key & 0xffffffff;
    return {u, v};
}

FlowNetwork buildFlowNetwork(int n, const vector<vector<int>>& adj, const unordered_map<long long, int>& edge_ids, const vector<array<int, 3>>& triangles, const vector<int>& tri_deg, double alpha) {
    static const double INF = 1e18;

    const int edges = edge_ids.size();
    int N = n + edges;
    int s = N, t = N + 1;
    FlowNetwork fn(N + 2);

    // s to v
    for (int v = 0; v < n; v++) {
        fn.addEdge(s, v, tri_deg[v]);
    }

    // v to t
    for (int v = 0; v < n; v++) {
        double cap = 3 * alpha;
        fn.addEdge(v, t, cap);
    }

    // edges-nodes -> vertices
    for (const auto& p : edge_ids) {
        auto [u, v] = decodeEdge(p.first);
        int e_node = n + p.second;
        fn.addEdge(e_node, u, INF);
        fn.addEdge(e_node, v, INF);
    }

    // triangle connections
    const auto min_max_pair = [&](int u, int v) -> long long {
        return edgeKey(min(u, v), max(u, v));
    };

    for (const auto& tri : triangles) {
        int u = tri[0], v = tri[1], w = tri[2];
        int e_uv = n + edge_ids.at(min_max_pair(u, v));
        int e_vw = n + edge_ids.at(min_max_pair(v, w));
        int e_uw = n + edge_ids.at(min_max_pair(u, w));

        fn.addEdge(u, e_vw, 1);
        fn.addEdge(v, e_uw, 1);
        fn.addEdge(w, e_uv, 1);
    }

    return fn;
}

vector<bool> minCutReachable(int s, FlowNetwork& fn) {
    vector<bool> visited(fn.N, false);

    queue<int> q{};
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

vector<int> findDensestSubgraph(const vector<vector<int>>& adj) {
    int n = adj.size();

    // compute triangle degrees only once
    vector<array<int, 3>> triangles = getTriangles(n, adj);
    vector<int> tri_deg(n, 0);
    for (const auto& tri : triangles) {
        tri_deg[tri[0]]++;
        tri_deg[tri[1]]++;
        tri_deg[tri[2]]++;
    }

    // map edge to ids
    unordered_map<long long, int> edge_ids{};
    int eid = 0;
    for (int u = 0; u < n; u++)
        for (const int& v : adj[u])
            if (u < v)
                edge_ids[edgeKey(u, v)] = eid++;

    double u = 0, l = 0;
    for (int i = 0; i < n; i++) {
        u = max(u, (double) tri_deg[i]);
    }

    vector<int> densestSubgraph{};

    double thresh = max(1e-6, 1.0 / (1.0 * n * (n - 1)));
    while (u - l > thresh) {
        double alpha = l + (u - l) / 2.0;
        // cout << "Building flow network with alpha = " << alpha << endl;

        FlowNetwork fn = buildFlowNetwork(n, adj, edge_ids, triangles, tri_deg, alpha);
        int s = fn.N - 2, t = fn.N - 1;

        // double maxflow = edmondsKarp(s, t, fn);
        double maxflow = dinic(s, t, fn);
        // cout << "Max flow = " << maxflow << endl;

        vector<bool> reachable = minCutReachable(s, fn);

        vector<int> S{};
        // exclude s and t
        for (int i = 0; i < n; i++)
            if (reachable[i])
                S.push_back(i);

        if (S.empty())
            u = alpha;
        else {
            l = alpha;
            densestSubgraph = S;
        }
    }
    return densestSubgraph;
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

void display_memory_usage(ostream& out) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    // out << "Peak memory usage: " << usage.ru_maxrss / (1024 * 1024) << "MB" << endl; // for macos, unit is in bytes
    out << "Peak memory usage: " << usage.ru_maxrss / 1024 << "MB" << endl; // for linux, unit is in kb, uncomment this line
}

int main(int argc, char* argv[]) {

    string filename = "../datasets/email-Enron.txt";
    if (argc > 1)
        filename = argv[1];

    vector<vector<int>> adj = buildGraphFromInput(filename);

    ostream* out = &cout;
    ofstream file;

    if (argc > 2) {
        file.open(argv[2]);
        out = &file;
    }

    auto start = chrono::high_resolution_clock::now();
    vector<int> DS = findDensestSubgraph(adj);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> diff = end - start;
    *out << "Densest Subgraph found in " << diff.count() << " seconds." << endl;

    display_memory_usage(*out);

    // get only neighbours in DS
    unordered_set<int> inDS(DS.begin(), DS.end());

    // calc density
    int triangle_count = 0;

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

    // output
    *out << "Clique density of densest subgraph: " << density << endl;
    *out << "Nodes which are part of densest subgraph: \n";
    for (const int& node : DS) {
        *out << node << " ";
    }
    *out << endl;

    return 0;
}