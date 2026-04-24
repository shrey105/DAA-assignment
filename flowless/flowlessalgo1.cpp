#include <bits/stdc++.h>
using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
//  SNAP Dataset Loader
// ─────────────────────────────────────────────────────────────────────────────
vector<pair<int,int>> load_snap(
    const string&  filepath,
    map<int,int>&  id_to_new,
    vector<int>&   new_to_id)
{
    ifstream fin(filepath);
    if (!fin.is_open()) {
        cerr << "ERROR: cannot open file: " << filepath << "\n";
        exit(1);
    }

    vector<pair<int,int>> edges;
    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream ss(line);
        int u, v;
        if (!(ss >> u >> v)) continue;
        if (!id_to_new.count(u)) {
            int nid = (int)id_to_new.size();
            id_to_new[u] = nid;
            new_to_id.push_back(u);
        }
        if (!id_to_new.count(v)) {
            int nid = (int)id_to_new.size();
            id_to_new[v] = nid;
            new_to_id.push_back(v);
        }
        int nu = id_to_new[u], nv = id_to_new[v];
        if (nu != nv) edges.push_back({nu, nv});
    }
    return edges;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Charikar's Greedy Peeling — Densest Subgraph (Fixed)
// ─────────────────────────────────────────────────────────────────────────────
vector<int> greedyDensestSubgraph(int n, vector<pair<int,int>>& edges) {
    vector<set<int>> adj(n);
    for (auto& [u, v] : edges) {
        if (u != v) { // Prevent self-loops from skewing density
            adj[u].insert(v);
            adj[v].insert(u);
        }
    }

    // Calculate the exact number of unique, undirected edges
    int curEdges = 0;
    for (int i = 0; i < n; i++) {
        curEdges += adj[i].size();
    }
    curEdges /= 2;

    vector<bool> active(n, true);
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<>> pq;
    
    for (int v = 0; v < n; v++) {
        pq.push({(int)adj[v].size(), v});
    }

    int curNodes = n;
    double bestDensity = (curNodes > 0) ? (double)curEdges / curNodes : 0.0;
    int bestCutoff = 0; 
    vector<int> removalOrder;
    removalOrder.reserve(n);

    while (curNodes > 0) {
        int u = -1;
        
        // Standard Dijkstra-style lazy extraction
        while (!pq.empty()) {
            auto [deg, v] = pq.top();
            pq.pop();
            if (!active[v]) continue;
            
            // If the degree matches the current actual degree, we found our minimum
            if (deg == (int)adj[v].size()) {
                u = v;
                break;
            }
            // If deg != adj[v].size(), it's an outdated entry. Just ignore it.
        }

        if (u == -1) break; // Priority queue exhausted

        active[u] = false;
        removalOrder.push_back(u);

        for (int nb : adj[u]) {
            if (active[nb]) {
                adj[nb].erase(u);
                curEdges--;
                
                // Push the newly decreased degree into the min-heap immediately
                pq.push({(int)adj[nb].size(), nb}); 
            }
        }
        
        adj[u].clear();
        curNodes--;

        // Track best density AND the cutoff index
        if (curNodes > 0) {
            double density = (double)curEdges / curNodes;
            if (density > bestDensity) {
                bestDensity = density;
                bestCutoff = (int)removalOrder.size();
            }
        }
    }

    // Build result: vertices NOT in the first bestCutoff removals
    set<int> removed(removalOrder.begin(), removalOrder.begin() + bestCutoff);
    vector<int> result;
    for (int v = 0; v < n; v++) {
        if (!removed.count(v)) {
            result.push_back(v);
        }
    }

    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    string filepath = "Email-Enron.txt";   // default
    if (argc >= 2) filepath = argv[1];

    cout << "Loading dataset : " << filepath << endl;

    map<int,int> id_to_new;
    vector<int>  new_to_id;

    auto t0 = chrono::high_resolution_clock::now();
    vector<pair<int,int>> edges = load_snap(filepath, id_to_new, new_to_id);
    auto t1 = chrono::high_resolution_clock::now();

    int N = (int)id_to_new.size();
    cout << "Vertices loaded : " << N << endl;
    cout << "Edges loaded    : " << edges.size() << endl;
    cout << "Load time       : "
         << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count()
         << " ms\n";

    cout << "\nRunning greedy densest subgraph..." << endl;
    auto t2 = chrono::high_resolution_clock::now();
    vector<int> densest = greedyDensestSubgraph(N, edges);
    auto t3 = chrono::high_resolution_clock::now();

    cout << "Algorithm time  : "
         << chrono::duration_cast<chrono::milliseconds>(t3 - t2).count()
         << " ms\n";

    // Compute density of result
    set<int> inResult(densest.begin(), densest.end());
    int subEdges = 0;
    for (auto& [u, v] : edges)
        if (inResult.count(u) && inResult.count(v))
            subEdges++;

    // Divide by 2 because each undirected edge is counted twice if (u,v) and (v,u) exist, 
    // but the input might just have (u,v). Let's do a precise check:
    // To be safe and consistent with the algorithm's density logic, we'll build an exact sub-graph edge count.
    subEdges = 0;
    set<pair<int,int>> unique_edges;
    for (auto& [u, v] : edges) {
        if (inResult.count(u) && inResult.count(v)) {
            int mn = min(u, v);
            int mx = max(u, v);
            unique_edges.insert({mn, mx});
        }
    }
    subEdges = unique_edges.size();

    double density = (densest.size() > 0) ? (double)subEdges / densest.size() : 0.0;

    cout << "\n=== Densest Subgraph Results ===" << endl;
    cout << "Vertices in densest subgraph : " << densest.size() << endl;
    cout << "Edges in densest subgraph    : " << subEdges << endl;
    cout << "Density (edges/vertices)     : " << fixed << setprecision(4) << density << endl;

    // Save results
    string out_file = filepath + "_densest_results.txt";
    ofstream fout(out_file);
    fout << "# Densest Subgraph Results\n";
    fout << "# Dataset: " << filepath << "\n";
    fout << "# Density: " << density << "\n";
    fout << "# Format: original_node_id\n";
    for (int v : densest)
        fout << new_to_id[v] << "\n";
    fout.close();
    cout << "Results saved to: " << out_file << endl;

    return 0;
}