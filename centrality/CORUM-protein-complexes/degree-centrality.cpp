
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <bitset>
#include <chrono>
#include <iomanip>
#include <algorithm>

using namespace std;

// --- Configuration & Data Structures (from original betweenness code) ---
string init[] = {"dataset_init.txt"};
vector<string> dataset;
vector<string> graphs;
vector<string> hypergraphs;
vector<string> names;

const int MAXN = 4e6;
const int MAXHE = 3e6;
const double epsilon = 1e-9;
vector<int> graph[MAXN + MAXHE];
map<int, int> idx_to_node;
map<int, int> node_to_idx;
vector<double> degree_centrality; // Changed from 'betweenness'
bitset<MAXN + MAXHE> mask;

// --- Memory Calculation Utilities (from original betweenness code) ---
template<class K, class V>
unsigned long long getMemoryUsage(const map<K, V>& m) {
    unsigned long long nodeSize = sizeof(K) + sizeof(V) + 3 * sizeof(void*) + sizeof(bool);
    unsigned long long adminSize = 3 * sizeof(void*) + sizeof(size_t);
    return adminSize + m.size() * nodeSize;
}

unsigned long long getMemoryUsage2(vector<int> g[], int cnt) {
    unsigned long long totalSize = 0;
    for (int i = 0; i < cnt; ++i)
        totalSize += g[i].size() * sizeof(int);
    return totalSize;
}

template<class T>
unsigned long long getMemoryUsage3(const vector<T>& v, unsigned long long cnt) {
    return cnt * sizeof(v[0]);
}


int main(){
    string output_dir_name = "Output/Degree/";
    cout << "INFO: Ensure the output directory '" << output_dir_name << "' exists." << endl;

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    ifstream fin;
    ofstream fout;

    // --- Read dataset configurations (from original betweenness code) ---
    for (auto& s: init){
        fin.open(s);
        string line;
        while(getline(fin, line)){
            if (line.empty()) continue;
            dataset.push_back(line);
            getline(fin, line); graphs.push_back(line);
            getline(fin, line); hypergraphs.push_back(line);
            getline(fin, line); names.push_back(line);
        }
        fin.close();
    }

    printf("%s\n", string(80, '=').c_str());
    for (int i = 0; i < dataset.size(); ++i){
        printf("Dataset #%d (%s):\n\n", i+1, dataset[i].c_str());
        
        // ===================================================================
        //  Part 1: Standard Graph Representation
        // ===================================================================
        {
            auto begin = chrono::high_resolution_clock::now();
            fin.open(graphs[i]);
            int nodes, edges;
            fin >> nodes >> edges;

            // --- Data loading structure (from original betweenness code) ---
            for (int j = 0; j < MAXN + MAXHE; ++j) graph[j].clear();
            idx_to_node.clear();
            node_to_idx.clear();
            int u, v, cnt = 0;
            while(fin >> u >> v){
                if(node_to_idx.find(u) == node_to_idx.end()){ idx_to_node[cnt] = u; node_to_idx[u] = cnt++; }
                u = node_to_idx[u];
                if(node_to_idx.find(v) == node_to_idx.end()){ idx_to_node[cnt] = v; node_to_idx[v] = cnt++; }
                v = node_to_idx[v];
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
            fin.close();

            // --- THIS IS THE CORE LOGIC FROM YOUR GROUPMATE'S NEW CODE ---
            // It replaces the entire complex 'solve()' function and loop.
            degree_centrality.assign(cnt, 0);
            for(int j = 0; j < cnt; ++j) {
                degree_centrality[j] = graph[j].size();
            }
            printf("Degree centrality calculated for graph.\n");

            // --- Performance tracking and output formatting (from original betweenness code) ---
            unsigned long long memory = 0;
            memory += getMemoryUsage(idx_to_node);
            memory += getMemoryUsage(node_to_idx);
            memory += getMemoryUsage2(graph, cnt);
            memory += getMemoryUsage3(degree_centrality, cnt);

            printf("For graph representation:\n");
            printf("Memory Used: %llu Bytes\n", memory);
            auto end = chrono::high_resolution_clock::now();
            double time = chrono::duration_cast<chrono::duration<double>>(end-begin).count();
            printf("Time Taken: %lf Seconds\n", time);

            string out_fname = output_dir_name + names[i] + ".ungraph.degree.txt";
            fout.open(out_fname);
            fout << memory << " " << fixed << setprecision(9) << time << endl << endl;
            for (auto const& [node_id, index] : node_to_idx) {
                fout << node_id << " " << degree_centrality[index] << endl;
            }
            fout.close();

            double max_degree = 0;
            if(!degree_centrality.empty()) max_degree = *max_element(degree_centrality.begin(), degree_centrality.end());

            printf("The node(s) with the greatest degree centrality are ");
            bool first = true;
            for (int j = 0; j < cnt; ++j) {
                if (fabs(max_degree - degree_centrality[j]) < epsilon) {
                    if (!first) printf(", "); printf("%d", idx_to_node[j]); first = false;
                }
            }
            printf(" (Degree: %.0f).\n\n", max_degree);
        }

        // ===================================================================
        //  Part 2: Hypergraph Representation
        // ===================================================================
        {
            auto begin = chrono::high_resolution_clock::now();
            fin.open(hypergraphs[i]);

            // --- Bipartite graph loading structure (from original betweenness code) ---
            for (int j = 0; j < MAXN + MAXHE; ++j) graph[j].clear();
            idx_to_node.clear();
            node_to_idx.clear();
            mask.reset();

            string line;
            int u, v, cnt = 0;
            int hyperEdge_count = 0;
            while(getline(fin, line)){
                if (line.empty()) continue;
                stringstream buffer(line);
                v = MAXN + hyperEdge_count;
                if(node_to_idx.find(v) == node_to_idx.end()){ idx_to_node[cnt] = v; node_to_idx[v] = cnt++; }
                v = node_to_idx[v];
                while(buffer >> u){
                    if(node_to_idx.find(u) == node_to_idx.end()){ idx_to_node[cnt] = u; node_to_idx[u] = cnt++; }
                    u = node_to_idx[u];
                    mask.set(u);
                    graph[u].push_back(v);
                    graph[v].push_back(u);
                }
                hyperEdge_count++;
            }
            fin.close();

            // --- CORE LOGIC for hypergraph degree ---
            // Again, a simple loop replaces the complex betweenness calculation.
            degree_centrality.assign(cnt, 0);
            for(int j = 0; j < cnt; ++j) {
                degree_centrality[j] = graph[j].size();
            }
            printf("Degree centrality calculated for hypergraph (bipartite model).\n");

            // --- Performance tracking and output formatting (from original betweenness code) ---
            unsigned long long memory = 0;
            memory += getMemoryUsage(idx_to_node);
            memory += getMemoryUsage(node_to_idx);
            memory += getMemoryUsage2(graph, cnt);
            memory += getMemoryUsage3(degree_centrality, cnt);

            printf("For hypergraph representation:\n");
            printf("Memory Used: %llu Bytes\n", memory);
            auto end = chrono::high_resolution_clock::now();
            double time = chrono::duration_cast<chrono::duration<double>>(end-begin).count();
            printf("Time Taken: %lf Seconds\n", time);
            
            string out_fname = output_dir_name + names[i] + ".cmty.degree.txt";
            fout.open(out_fname);
            fout << memory << " " << fixed << setprecision(9) << time << endl << endl;
            for (auto const& [node_id, index] : node_to_idx) {
                if (mask[index]) { // Only output protein nodes
                    fout << node_id << " " << degree_centrality[index] << endl;
                }
            }
            fout.close();
            
            double max_degree = 0;
            for (int j = 0; j < cnt; ++j) if (mask[j]) max_degree = max(max_degree, degree_centrality[j]);

            printf("The node(s) with the greatest degree centrality are ");
            bool first = true;
            for (int j = 0; j < cnt; ++j) {
              if (mask[j] && fabs(max_degree - degree_centrality[j]) < epsilon) {
                    if (!first) printf(", "); printf("%d", idx_to_node[j]); first = false;
              }
            }
            printf(" (Degree: %.0f).\n", max_degree);
        }

        printf("%s\n", string(80, '=').c_str());
    }

    return 0;
}