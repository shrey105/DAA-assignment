#include <chrono>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct Graph {
    int n;
    int m;
    vector<vector<int>> adj;

    Graph(int n) : n(n), m(0), adj(n) {}

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
        m++;
    }
};

struct DensestSubgraph {
    double density;
    vector<int> nodes;
};

DensestSubgraph greedy_plus_plus(const Graph& G, int T) {
    const int n = G.n;

    vector<long long> load(n, 0LL);

    DensestSubgraph best;
    best.density = (n > 0) ? (double) G.m / n : 0.0;
    best.nodes.resize(n);
    for (int i = 0; i < n; i++)
        best.nodes[i] = i;

    vector<int> cur_deg(n);
    vector<bool> alive(n);
    vector<int> peeling_order;
    peeling_order.reserve(n);

    for (int iter = 1; iter <= T; ++iter) {
        // initialize
        for (int v = 0; v < n; ++v) {
            cur_deg[v] = (int) G.adj[v].size();
            alive[v] = true;
        }

        peeling_order.clear();

        int H_size = n;
        long long H_edges = G.m;

        using PQEntry = pair<long long, int>;
        priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>> pq;

        for (int v = 0; v < n; ++v)
            pq.push({load[v] + cur_deg[v], v});

        vector<long long> new_load(n, 0LL);

        double best_iter_density = best.density;
        int best_cut = -1; // number of removed nodes at best point

        int removed = 0;

        while (!pq.empty()) {
            auto [key, u] = pq.top();
            pq.pop();

            if (!alive[u])
                continue;
            if (key != load[u] + cur_deg[u]) {
                pq.push({load[u] + cur_deg[u], u});
                continue;
            }

            new_load[u] = load[u] + cur_deg[u];

            alive[u] = false;
            peeling_order.push_back(u);
            removed++;

            H_edges -= cur_deg[u];
            H_size--;

            for (int nb : G.adj[u]) {
                if (!alive[nb])
                    continue;
                cur_deg[nb]--;
                pq.push({load[nb] + cur_deg[nb], nb});
            }

            if (H_size > 0) {
                double rho = (double) H_edges / H_size;
                if (rho > best_iter_density) {
                    best_iter_density = rho;
                    best_cut = removed;
                }
            }
        }

        // reconstruct candidate ONCE
        vector<bool> removed_flag(n, false);
        for (int i = 0; i < best_cut; ++i)
            removed_flag[peeling_order[i]] = true;

        vector<int> candidate;
        candidate.reserve(n);
        for (int v = 0; v < n; ++v)
            if (!removed_flag[v])
                candidate.push_back(v);

        if (best_iter_density > best.density) {
            best.density = best_iter_density;
            best.nodes.swap(candidate);
        }

        load.swap(new_load);
    }

    return best;
}

Graph read_graph(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file " << filename << "\n";
        exit(1);
    }

    vector<pair<int, int>> edges;
    int u, v, maxV = 0;
    string line;

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        stringstream ss(line);
        if (ss >> u >> v) {
            if (u == v)
                continue;
            if (u > v)
                swap(u, v); // normalize
            edges.emplace_back(u, v);
            maxV = max(maxV, max(u, v));
        }
    }

    // sort + deduplicate (O(m log m), but cache-friendly and much faster than set)
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());

    Graph G(maxV + 1);

    for (auto& [a, b] : edges)
        G.add_edge(a, b);

    return G;
}

int main(int argc, char* argv[]) {
    cin.tie(nullptr);

    Graph G(0);
    int T = 10;

    if (argc >= 3) {
        G = read_graph(argv[1]);
        T = atoi(argv[2]);
    } else if (argc == 2) {
        G = read_graph(argv[1]);
    }

    cout << "Graph: n=" << G.n << "  m=" << G.m
         << "  iterations=" << T << "\n\n";

    auto start = chrono::high_resolution_clock::now();
    DensestSubgraph result = greedy_plus_plus(G, T);
    auto end = chrono::high_resolution_clock::now();

    double elapsed_ms =
        chrono::duration<double, milli>(end - start).count();

    cout << "Densest subgraph found:\n";
    cout << "  Density : " << fixed << setprecision(6) << result.density << "\n";
    cout << "  Size    : " << result.nodes.size() << " nodes\n";

    if (G.n <= 50 && !result.nodes.empty()) {
        cout << "  Nodes   : ";
        for (int v : result.nodes)
            cout << v << " ";
        cout << "\n";
    }

    cout << "Runtime : " << fixed << setprecision(3) << elapsed_ms << " ms\n";
    return 0;
}