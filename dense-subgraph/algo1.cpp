#include <bits/stdc++.h>
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

    FlowNetwork(int n) : N{n}, graph{n} {}

    void addEdge(int u, int v, double cap) {
        graph[u].emplace_back(v, cap, graph[v].size());
        graph[v].emplace_back(u, 0.0, graph[u].size() - 1);
    }
};

double edmondsKarp(int s, int t, FlowNetwork& fn) {
    double maxflow = 0.0;

    while (true) {
        vector<int> parent(fn.N, -1); // node
        vector<int> parentEdge(fn.N, -1); // edge index

        queue<int> q {};
        q.push(s);
        parent[s] = s;

        // BFS for augmenting path
        while (!q.empty()) {
            int u = q.front(); q.pop();

            for (int i = 0; i < fn.graph[u].size(); i++) {
                Edge &e = fn.graph[u][i];

                if (parent[e.to] == -1 && e.res_cap > 1e-9) {
                    parent[e.to] = u;
                    parentEdge[e.to] = i;
                    q.push(e.to);
                }
            }
        }
        if (parent[t] == -1) break;

        // bottleneck edge
        double flow = 1e9;
        int v = t;
        while (v != s) {
            int u = parent[v];
            Edge &e = fn.graph[u][parentEdge[v]];
            flow = min(flow, e.res_cap);
            v = u;
        }

        v = t;
        while (v != s) {
            int u = parent[v];
            Edge &e = fn.graph[u][parentEdge[v]];
            e.res_cap -= flow;
            fn.graph[v][e.rev_idx].res_cap += flow;
            v = u;
        }

        maxflow += flow;
    }
    return maxflow;
} 

FlowNetwork buildFlowNetwork(int n, const vector<vector<int>>& adj, double alpha) {
    int s = n, t = n + 1;
    FlowNetwork fn(n + 2);

    int m = 0;
    vector<int> degree(n, 0);

    for (int u = 0; u < n; u++) {
        degree[u] = adj[u].size();
        m += degree[u];
    }
    m /= 2;

    // s to v
    for (int v = 0; v < n; v++) {
        fn.addEdge(s, v, m);
    }

    // v to t
    for (int v = 0; v < n; v++) {
        double cap = m + 2 * alpha - degree[v];
        fn.addEdge(v, t, cap);
    }

    // edges
    for (int u = 0; u < n; u++) {
        for (int v: adj[u]) {
            if (u < v) {
                fn.addEdge(u, v, 1);
                fn.addEdge(v, u, 1);
            }
        }
    }
    return fn;
}

vector<bool> minCutReachable(int s, FlowNetwork& fn) {
    vector<bool> visited(fn.N, false);

    queue<int> q {};
    q.push(s);
    visited[s] = true;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto &e: fn.graph[u]) {
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
    double u = 0, l = 0;
    for (int i = 0; i < n; i++) {
        u = max(u, (double) adj[i].size());
    }

    vector<int> densestSubgraph {};

    while (u - l > 1.0 / (n * (n - 1))) {
        double alpha = (l + u) / 2.0;
        
        FlowNetwork fn = buildFlowNetwork(n, adj, alpha);
        int s = n, t = n + 1;

        double maxflow = edmondsKarp(s, t, fn);

        vector<bool> reachable = minCutReachable(s, fn);

        vector<int> S {};
        // exclude s and t
        for (int i = 0; i < n; i++)
            if (reachable[i]) S.push_back(i);
        
        if (S.empty())
            u = alpha;
        else {
            l = alpha;
            densestSubgraph = S;
        }
    }
    return densestSubgraph;
}

vector<vector<int>> buildGraphFromInput(FILE *file);

int main() {

    return 0;
}