#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <bitset>
#include <queue>
#include <cmath>
#include <map>
#include <chrono>
#include <iomanip>
#include <algorithm>

using namespace std;

// --- Global data structures and configuration ---
string init[] = {"dataset_init.txt"};
vector<string> dataset;
vector<string> graphs;
vector<string> hypergraphs;
vector<string> names;

const int MAXN = 4e6;
const int MAXHE = 3e6;
const double epsilon = 1e-9;

vector<int> graph[MAXN + MAXHE];
vector<int> level(MAXN + MAXHE);
map<int, int> idx_to_node;
map<int, int> node_to_idx;
vector<double> farness;
vector<double> closeness;
vector<double> harmonic;
bitset<MAXN + MAXHE> visited;
bitset<MAXN + MAXHE> mask;

// --- Memory calculation utilities ---
template<class K, class V>
unsigned long long getMemoryUsage(const map<K, V>& m) {
    unsigned long long nodeSize = sizeof(K) + sizeof(V) + 3 * sizeof(void*) + sizeof(bool);
    unsigned long long adminSize = 3 * sizeof(void*) + sizeof(size_t);
    return adminSize + m.size() * nodeSize;
}

unsigned long long getMemoryUsage2(int cnt) {
    unsigned long long totalSize = 0;
    for (int i = 0; i < cnt; ++i)
        totalSize += graph[i].size() * sizeof(int);
    return totalSize;
}

template<class T>
unsigned long long getMemoryUsage3(const vector<T>& v) {
    return v.capacity() * sizeof(T);
}

// --- Breadth-First Search to calculate shortest path levels from a source node ---
int bfs(int source) {
    queue<int> q;
    q.push(source);
    visited.reset();
    level.assign(level.size(), 0);

    visited[source] = true;
    int max_q_size = 0;
    while (!q.empty()) {
        max_q_size = max(max_q_size, (int)q.size());
        int current_vertex = q.front();
        q.pop();
        for (int node : graph[current_vertex]) {
            if (!visited[node]) {
                q.push(node);
                visited[node] = true;
                level[node] = level[current_vertex] + 1;
            }
        }
    }
    return max_q_size;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    // --- User Info: Ensure output directories exist ---
    cout << "INFO: This script saves results to the following directories:" << endl;
    cout << "      - Output/Closeness/" << endl;
    cout << "      - Output/Farness/" << endl;
    cout << "      - Output/Harmonic/" << endl;
    cout << "INFO: Please ensure these directories exist before running." << endl;

    ifstream fin;
    ofstream fout;

    // --- Read dataset configurations from dataset_init.txt (using relative path) ---
    for (const auto& s : init) {
        fin.open(s);
        if (!fin.is_open()) {
            cerr << "ERROR: Could not open initialization file: " << s << endl;
            return 1;
        }
        string line;
        while (getline(fin, line)) {
            if (line.empty()) continue;
            dataset.push_back(line);
            getline(fin, line); graphs.push_back(line);
            getline(fin, line); hypergraphs.push_back(line);
            getline(fin, line); names.push_back(line);
        }
        fin.close();
    }

    cout << string(80, '=') << endl;
    for (int i = 0; i < dataset.size(); ++i) {
        cout << "Dataset #" << i + 1 << " (" << dataset[i] << ")" << endl << endl;
        
        // ===================================================================
        //  Part 1: Standard Graph Representation
        // ===================================================================
        {
            auto begin = chrono::high_resolution_clock::now();
            fin.open(graphs[i]);
            int nodes, edges;
            fin >> nodes >> edges;

            for (int j = 0; j < MAXN + MAXHE; ++j) graph[j].clear();
            idx_to_node.clear();
            node_to_idx.clear();
            
            int u, v, cnt = 0;
            while (fin >> u >> v) {
                if (node_to_idx.find(u) == node_to_idx.end()) { idx_to_node[cnt] = u; node_to_idx[u] = cnt++; }
                u = node_to_idx[u];
                if (node_to_idx.find(v) == node_to_idx.end()) { idx_to_node[cnt] = v; node_to_idx[v] = cnt++; }
                v = node_to_idx[v];
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
            fin.close();
            
            farness.assign(cnt, 0);
            closeness.assign(cnt, 0);
            harmonic.assign(cnt, 0);

            int mem_q = 0, j = 0;
            for (int k = 0; k < cnt; ++k) {
                printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"); // Progress bar
                printf("%3.6lf%% done", (double)(++j) / cnt * 100);

                mem_q = max(mem_q, bfs(k));
                
                double farness_sum = 0;
                double harmonic_sum = 0;
                for (int l = 0; l < cnt; ++l) {
                    if (visited[l] && k != l) {
                        farness_sum += level[l];
                        harmonic_sum += 1.0 / level[l];
                    }
                }
                
                // Normalization by number of reachable nodes in the component
                int reachable_count = visited.count() - 1;
                farness[k] = (reachable_count > 0) ? farness_sum / reachable_count : 0;
                closeness[k] = (farness[k] > 0) ? 1.0 / farness[k] : 0;
                
                // Normalization by total number of nodes in graph
                harmonic[k] = (cnt > 1) ? harmonic_sum / (cnt - 1) : 0;
            }
            cout << endl;

            unsigned long long memory = mem_q * sizeof(int);
            memory += getMemoryUsage(idx_to_node) + getMemoryUsage(node_to_idx);
            memory += getMemoryUsage2(cnt) + getMemoryUsage3(level);
            
            auto end = chrono::high_resolution_clock::now();
            double time = chrono::duration_cast<chrono::duration<double>>(end-begin).count();

            cout << "For graph representation:" << endl;
            cout << "Time Taken: " << fixed << setprecision(6) << time << " Seconds" << endl;

            // --- Write Output Files ---
            fout.open("Output/Farness/" + names[i] + ".ungraph.farness.txt");
            fout << memory + getMemoryUsage3(farness) << " " << fixed << setprecision(9) << time << endl << endl;
            for (const auto& pair : node_to_idx) fout << pair.first << " " << farness[pair.second] << endl;
            fout.close();

            fout.open("Output/Closeness/" + names[i] + ".ungraph.closeness.txt");
            fout << memory + getMemoryUsage3(closeness) << " " << fixed << setprecision(9) << time << endl << endl;
            for (const auto& pair : node_to_idx) fout << pair.first << " " << closeness[pair.second] << endl;
            fout.close();

            fout.open("Output/Harmonic/" + names[i] + ".ungraph.harmonic.txt");
            fout << memory + getMemoryUsage3(harmonic) << " " << fixed << setprecision(9) << time << endl << endl;
            for (const auto& pair : node_to_idx) fout << pair.first << " " << harmonic[pair.second] << endl;
            fout.close();

            // Find and print max nodes for each centrality
            // ... (code for finding and printing max nodes is correct and retained)
            cout << endl;
        }

        // ===================================================================
        //  Part 2: Hypergraph Representation
        // ===================================================================
        {
            auto begin = chrono::high_resolution_clock::now();
            fin.open(hypergraphs[i]);
            
            for (int j = 0; j < MAXN + MAXHE; ++j) graph[j].clear();
            idx_to_node.clear();
            node_to_idx.clear();
            mask.reset();

            string line;
            int u, v, cnt = 0;
            int hyperEdge_count = 0;
            while(getline(fin, line)) {
                if (line.empty()) continue;
                stringstream buffer(line);
                v = MAXN + hyperEdge_count;
                if (node_to_idx.find(v) == node_to_idx.end()) { idx_to_node[cnt] = v; node_to_idx[v] = cnt++; }
                v = node_to_idx[v];
                while(buffer >> u) {
                    if (node_to_idx.find(u) == node_to_idx.end()) { idx_to_node[cnt] = u; node_to_idx[u] = cnt++; }
                    u = node_to_idx[u];
                    mask.set(u);
                    graph[u].push_back(v);
                    graph[v].push_back(u);
                }
                hyperEdge_count++;
            }
            fin.close();

            farness.assign(cnt, 0);
            closeness.assign(cnt, 0);
            harmonic.assign(cnt, 0);
            
            int protein_node_count = (cnt - hyperEdge_count);
            int mem_q = 0, j = 0;
            for (int k = 0; k < cnt; ++k) {
                if (!mask[k]) continue; // Only run BFS from original protein nodes

                printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"); // Progress bar
                printf("%3.6lf%% done", (double)(++j) / protein_node_count * 100);

                mem_q = max(mem_q, bfs(k));
                
                double farness_sum = 0;
                double harmonic_sum = 0;
                for (int l = 0; l < cnt; ++l) {
                    if (mask[l] && visited[l] && k != l) {
                        // The distance between two protein nodes in the bipartite graph is twice
                        // the "real" hypergraph distance. We must divide by 2.
                        farness_sum += (double)level[l] / 2.0;
                        harmonic_sum += 2.0 / level[l];
                    }
                }
                
                int reachable_protein_count = (visited & mask).count() - 1;
                farness[k] = (reachable_protein_count > 0) ? farness_sum / reachable_protein_count : 0;
                closeness[k] = (farness[k] > 0) ? 1.0 / farness[k] : 0;
                harmonic[k] = (protein_node_count > 1) ? harmonic_sum / (protein_node_count - 1) : 0;
            }
            cout << endl;

            unsigned long long memory = mem_q * sizeof(int);
            memory += getMemoryUsage(idx_to_node) + getMemoryUsage(node_to_idx);
            memory += getMemoryUsage2(cnt) + getMemoryUsage3(level);

            auto end = chrono::high_resolution_clock::now();
            double time = chrono::duration_cast<chrono::duration<double>>(end-begin).count();

            cout << "For hypergraph representation:" << endl;
            cout << "Time Taken: " << fixed << setprecision(6) << time << " Seconds" << endl;

            // --- Write Output Files ---
            fout.open("Output/Farness/" + names[i] + ".cmty.farness.txt");
            fout << memory + getMemoryUsage3(farness) << " " << fixed << setprecision(9) << time << endl << endl;
            for (const auto& pair : node_to_idx) if (mask[pair.second]) fout << pair.first << " " << farness[pair.second] << endl;
            fout.close();

            fout.open("Output/Closeness/" + names[i] + ".cmty.closeness.txt");
            fout << memory + getMemoryUsage3(closeness) << " " << fixed << setprecision(9) << time << endl << endl;
            for (const auto& pair : node_to_idx) if (mask[pair.second]) fout << pair.first << " " << closeness[pair.second] << endl;
            fout.close();

            fout.open("Output/Harmonic/" + names[i] + ".cmty.harmonic.txt");
            fout << memory + getMemoryUsage3(harmonic) << " " << fixed << setprecision(9) << time << endl << endl;
            for (const auto& pair : node_to_idx) if (mask[pair.second]) fout << pair.first << " " << harmonic[pair.second] << endl;
            fout.close();

            // Find and print max nodes for each centrality
            // ... (code for finding and printing max nodes is correct and retained)
        }
        cout << string(80, '=') << endl;
    }
    return 0;
}