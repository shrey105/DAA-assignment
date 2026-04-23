#include <chrono>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <queue>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct Graph {
    int n;                   // number of vertices
    int m;                   // number of edges
    vector<vector<int>> adj; // adjacency list

    Graph(int n, int m) : n(n), m(m), adj(n) {}

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

struct DensestSubgraph {
    double density;
    vector<int> nodes; // nodes in the densest subgraph found
};

DensestSubgraph greedy_plus_plus(const Graph& G, int T) {
    const int n = G.n;

    vector<long long> load(n, 0LL);

    DensestSubgraph best;
    best.density = (G.n > 0) ? (double) G.m / G.n : 0.0;
    best.nodes.resize(G.n);
    best.nodes.clear();
    for (int i = 0; i < n; i++) {
        best.nodes.push_back(i);
    }

    for (int iter = 1; iter <= T; ++iter) {
        vector<int> cur_deg(n);
        for (int v = 0; v < n; ++v)
            cur_deg[v] = static_cast<int>(G.adj[v].size());

        // alive[v] = true if v has not been peeled yet
        vector<bool> alive(n, true);

        // peeling_order[i] = vertex removed at step i
        vector<int> peeling_order;
        peeling_order.reserve(n);

        // Track which vertices are still in H for density computation
        int H_size = n;
        long long H_edges = G.m;

        // pq ordered by (load + cur_deg)
        using PQEntry = pair<long long, int>;
        priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>> pq;

        for (int v = 0; v < n; ++v)
            pq.push({load[v] + cur_deg[v], v});

        // new_load[v] = load accumulated in this iteration
        vector<long long> new_load(n, 0LL);

        while (!pq.empty()) {
            // Find vertex u with minimum (load + cur_deg)
            auto [key, u] = pq.top();
            pq.pop();

            // Skip stale entries (lazy deletion)
            if (!alive[u])
                continue;
            if (key != load[u] + cur_deg[u]) {
                pq.push({load[u] + cur_deg[u], u});
                continue;
            }

            // Record the load assigned in this iteration
            new_load[u] = load[u] + cur_deg[u];

            // Remove u and its adjacent edges from H
            alive[u] = false;
            peeling_order.push_back(u);

            H_edges -= cur_deg[u]; // edges leaving with u
            H_size -= 1;

            for (int nb : G.adj[u]) {
                if (!alive[nb])
                    continue;
                cur_deg[nb]--;
                // Push updated entry when key for nb changes
                pq.push({load[nb] + cur_deg[nb], nb});
            }

            // Check density of remaining H
            if (H_size > 0) {
                double rho = static_cast<double>(H_edges) / H_size;
                if (rho > best.density) {
                    best.density = rho;
                    best.nodes.clear();
                    for (int v = 0; v < n; ++v)
                        if (alive[v])
                            best.nodes.push_back(v);
                }
            }
        }

        // Update load vector for the next iteration
        load = new_load;
    }

    return best;
}

Graph read_graph(const string& filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Error: cannot open file " << filename << "\n";
        exit(1);
    }

    int n, m;
    fin >> n >> m;

    Graph G(n, m);
    for (int i = 0; i < m; ++i) {
        int u, v;
        fin >> u >> v;
        if (u != v) // skip self-loops
            G.add_edge(u, v);
    }
    return G;
}

Graph make_example_graph() {
    // Vertices 0..1 on left of bipartite, 2..11 on right, 12..15 form K_4
    const int n = 16;
    const int m_dummy = 0; // will be overwritten
    Graph G(n, m_dummy);
    G.m = 0;

    // K_{2,10}
    for (int l = 0; l < 2; ++l)
        for (int r = 2; r < 12; ++r) {
            G.add_edge(l, r);
            G.m++;
        }

    // K_4
    for (int u = 12; u < 16; ++u)
        for (int v = u + 1; v < 16; ++v) {
            G.add_edge(u, v);
            G.m++;
        }

    return G;
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Graph G(0, 0);
    int T = 10; // default iteration count

    if (argc >= 3) {
        G = read_graph(argv[1]);
        T = atoi(argv[2]);
    } else if (argc == 2) {
        G = read_graph(argv[1]);
    } else {
        cout << "No input file given — running on built-in example graph.\n\n";
        G = make_example_graph();
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

    // Print nodes only for small graphs
    if (G.n <= 50 && !result.nodes.empty()) {
        cout << "  Nodes   : ";
        for (int v : result.nodes)
            cout << v << " ";
        cout << "\n";
    }

    cout << "\nRuntime : " << fixed << setprecision(3) << elapsed_ms << " ms\n";
    return 0;
}