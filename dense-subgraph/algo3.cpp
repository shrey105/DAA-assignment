#include <bits/stdc++.h>
using namespace std;


using Graph  = vector<set<int>>;
using Clique = vector<int>;

void enumerate_cliques(
    const Graph&    G,
    const set<int>& active,
    int             sz,
    vector<Clique>& out)
{
    vector<int> verts(active.begin(), active.end());
    int n = (int)verts.size();

    function<void(int, Clique&)> bt = [&](int start, Clique& cur) {
        if ((int)cur.size() == sz) { out.push_back(cur); return; }
        int need = sz - (int)cur.size();
        for (int i = start; i <= n - need; ++i) {
            int v = verts[i];
            bool ok = true;
            for (int u : cur)
                if (!G[u].count(v)) { ok = false; break; }
            if (!ok) continue;
            cur.push_back(v);
            bt(i + 1, cur);
            cur.pop_back();
        }
    };

    Clique cur;
    bt(0, cur);
}

void compute_clique_degrees(
    const Graph&         G,
    const set<int>&      active,
    int                  clique_size,
    vector<int>&         cdeg,
    vector<Clique>&      cliques,
    vector<vector<int>>& v2cliques)
{
    enumerate_cliques(G, active, clique_size, cliques);
    for (int ci = 0; ci < (int)cliques.size(); ++ci)
        for (int v : cliques[ci]) {
            cdeg[v]++;
            v2cliques[v].push_back(ci);
        }
}

vector<int> kpsi_core_decomposition(const Graph& adj, int clique_size = 3)
{
    int N = (int)adj.size();
    vector<int> core(N, 0);
    Graph G = adj;

    set<int> active;
    for (int v = 0; v < N; ++v) active.insert(v);

    vector<int>          cdeg(N, 0);
    vector<Clique>       cliques;
    vector<vector<int>>  v2cliques(N);

    compute_clique_degrees(G, active, clique_size, cdeg, cliques, v2cliques);

    set<pair<int,int>> sorted_verts;
    for (int v = 0; v < N; ++v)
        sorted_verts.insert({cdeg[v], v});

    vector<bool> clique_alive(cliques.size(), true);

    while (!sorted_verts.empty()) {
        auto it = sorted_verts.begin();
        auto [deg_v, v] = *it;
        sorted_verts.erase(it);
        active.erase(v);
        core[v] = cdeg[v];

        for (int ci : v2cliques[v]) {
            if (!clique_alive[ci]) continue;
            clique_alive[ci] = false;
            for (int u : cliques[ci]) {
                if (u == v || !active.count(u)) continue;
                if (cdeg[u] > core[v]) {
                    sorted_verts.erase({cdeg[u], u});
                    cdeg[u]--;
                    sorted_verts.insert({cdeg[u], u});
                }
            }
        }

        for (int nbr : G[v]) G[nbr].erase(v);
        G[v].clear();
    }

    return core;
}

Graph load_snap(
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
        edges.push_back({u, v});
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
    }

    int N = (int)id_to_new.size();
    Graph G(N);
    for (auto [u, v] : edges) {
        int nu = id_to_new[u], nv = id_to_new[v];
        if (nu == nv) continue;
        G[nu].insert(nv);
        G[nv].insert(nu);
    }
    return G;
}


void print_results(
    const vector<int>& core,
    const vector<int>& new_to_id,
    int                clique_size,
    bool               verbose = false)
{
    int N = (int)core.size();
    int max_core = *max_element(core.begin(), core.end());

    cout << "\n=== (k, Psi)-Core Decomposition Results ===" << endl;
    cout << "Clique size (Psi) : " << clique_size << endl;
    cout << "Total vertices    : " << N << endl;
    cout << "Max clique-core # : " << max_core << endl;

    map<int,int> dist;
    for (int v = 0; v < N; ++v) dist[core[v]]++;

    cout << "\nCore-number distribution:" << endl;
    cout << "  Core#  |  # Vertices" << endl;
    cout << "---------|-------------" << endl;
    for (auto [k, cnt] : dist)
        cout << "    " << setw(4) << k << "   |   " << cnt << endl;

    if (verbose) {
        cout << "\nPer-vertex (internal_id -> original_id : core#):" << endl;
        for (int v = 0; v < N; ++v)
            cout << "  " << v << " -> " << new_to_id[v]
                 << " : " << core[v] << endl;
    }
}

int main(int argc, char* argv[])
{
    string filepath    = "email-Enron.txt";
    int    clique_size = 3;

    if (argc >= 2) filepath    = argv[1];
    if (argc >= 3) clique_size = stoi(argv[2]);

    if (clique_size < 2) {
        cerr << "ERROR: clique_size must be >= 2\n";
        return 1;
    }

    cout << "Loading dataset : " << filepath << endl;
    cout << "Clique size     : " << clique_size << endl;

    map<int,int> id_to_new;
    vector<int>  new_to_id;

    auto t0 = chrono::high_resolution_clock::now();
    Graph G = load_snap(filepath, id_to_new, new_to_id);
    auto t1 = chrono::high_resolution_clock::now();

    int N = (int)G.size();
    long long E = 0;
    for (int v = 0; v < N; ++v) E += G[v].size();
    E /= 2;

    cout << "Vertices loaded : " << N << endl;
    cout << "Edges loaded    : " << E << endl;
    cout << "Load time       : "
         << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count()
         << " ms\n";

    cout << "\nRunning (k, Psi)-core decomposition..." << endl;
    auto t2 = chrono::high_resolution_clock::now();
    vector<int> core = kpsi_core_decomposition(G, clique_size);
    auto t3 = chrono::high_resolution_clock::now();

    cout << "Algorithm time  : "
         << chrono::duration_cast<chrono::milliseconds>(t3 - t2).count()
         << " ms\n";

    // verbose=false: just show distribution; set true for per-vertex output
    print_results(core, new_to_id, clique_size, /*verbose=*/false);

    // Save results to file
    string out_file = filepath + "_core_results.txt";
    ofstream fout(out_file);
    fout << "# (k,Psi)-Core Decomposition Results\n";
    fout << "# Dataset: "     << filepath    << "\n";
    fout << "# Clique size: " << clique_size << "\n";
    fout << "# Format: original_node_id  clique_core_number\n";
    for (int v = 0; v < N; ++v)
        fout << new_to_id[v] << "\t" << core[v] << "\n";
    fout.close();
    cout << "\nResults saved to: " << out_file << endl;

    return 0;
}