#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <queue>
#include <map>
#include <cmath>
#include <bitset>
#include <chrono>
#include <iomanip>
using namespace std;
// #define endl '\n'

string init[] = {"dataset_init.txt"};
vector<string> dataset;
vector<string> graphs;
vector<string> hypergraphs;
vector<string> names;

const int MAXN = 4e6;
const int MAXHE = 3e6;
const double epsilon = 1e-9;
vector<int> graph[MAXN + MAXHE];
vector<int> predecessor[MAXN + MAXHE];
vector<int> level;
vector<double> sigma;
map<int, int> idx_to_node;
map<int, int> node_to_idx;
vector<double> betweenness;
bitset<MAXN + MAXHE> visited;
bitset<MAXN + MAXHE> mask;

template<class K, class V>
unsigned long long getMemoryUsage(const map<K, V>& m) {
  // Node structure estimation for Red-Black trees
  // With the key and the value, 2 child pointers, 1 parent pointer, and 1 color bit
  unsigned long long nodeSize = sizeof(K) + sizeof(V) + 3 * sizeof(void*) + sizeof(bool);
  // Administrative overhead assuming 3 pointers + 1 size_t
  unsigned long long adminSize = 3 * sizeof(void*) + sizeof(size_t);
  unsigned long long totalSize = adminSize + m.size() * nodeSize;

  return totalSize;
}

unsigned long long getMemoryUsage2(vector<int> g[], int cnt) {
  unsigned long long totalSize = 0;
  for (int i = 0; i < cnt; ++i)
    totalSize += g[i].size();
  totalSize *= sizeof(int);
  return totalSize;
}

template<class T>
unsigned long long getMemoryUsage3(const vector<T>& v, unsigned long long cnt) {
  return cnt * sizeof(v[0]);
}

// Brandes Algorithm
unsigned long long solve(int source, int cnt){
    stack<int> S;
    level.resize(cnt);
    queue<int> Q;
    sigma.assign(cnt, 0);
    for (int i = 0; i < cnt; ++i)
        predecessor[i].clear();
    visited.reset();
    unsigned long long mem = 0;

    sigma[source] = 1;
    level[source] = 0;
    Q.push(source);
    visited[source] = true;
    while(!Q.empty()){
        mem = max(mem, (unsigned long long)Q.size());
        int current_vertex = Q.front();
        Q.pop();
        S.push(current_vertex);
        for(auto node: graph[current_vertex]){
            if(!visited[node]) {
                Q.push(node);
                level[node] = level[current_vertex] +1;
            }
            visited[node] = true;
            // Update sigma value of node if path through current_vertex was shortest
            if(level[node] == level[current_vertex]+1){
                sigma[node] = sigma[node] + sigma[current_vertex];
                predecessor[node].push_back(current_vertex);
            }
        }
    }
    mem *= sizeof(int);
    mem += S.size() * sizeof(int);

    vector<double> delta(cnt, 0);
    // Pop elements out of the stack, starting from terminal node
    // work backward frontier b frontier, computing delta values
    while(!S.empty()) {
        int v = S.top();
        S.pop();
        for(auto& u: predecessor[v])
            delta[u] += (sigma[u]/sigma[v]) * (mask[v]+delta[v]);
        if(v != source)
            betweenness[v] += delta[v]/2;
    }

    mem += getMemoryUsage3(level, cnt);
    mem += getMemoryUsage3(sigma, cnt);
    mem += getMemoryUsage2(predecessor, cnt);
    mem += getMemoryUsage3(delta, cnt);
    return mem;
}

int main(){
    ios_base::sync_with_stdio(false);
    // cin.tie(NULL);
    // cin.exceptions(cin.failbit);

    ifstream fin;
    ofstream fout;
    for (auto& s: init){
        fin.open(s);
        string line;
        while(getline(fin, line)){
            if (line.empty())
                continue;
            dataset.push_back(line);
            getline(fin, line);
            graphs.push_back(line);
            getline(fin, line);
            hypergraphs.push_back(line);
            getline(fin, line);
            names.push_back(line);
        }
        fin.close();
    }

    printf("%s\n", string(80, '=').c_str());
    for (int i = 0; i < dataset.size(); ++i){
        printf("Dataset #%d (%s):\n\n", i+1, dataset[i].c_str());
        int nodes, edges;

        unsigned long long memory = 0;
        auto begin = chrono::high_resolution_clock::now();

        // Graph
        fin.open(graphs[i]);

        fin>>nodes>>edges;

        for (int i = 0; i < MAXN + MAXHE; ++i)
        graph[i].clear();
        idx_to_node.clear();
        node_to_idx.clear();
        betweenness.assign(MAXN + MAXHE, 0);

        int u, v, cnt = 0;
        mask.reset();
        while(fin>>u>>v){
            if(!node_to_idx.count(u)){
                idx_to_node[cnt] = u;
                node_to_idx[u] = cnt++;
            }
            u = node_to_idx[u];
            if(!node_to_idx.count(v)){
                idx_to_node[cnt] = v;
                node_to_idx[v] = cnt++;
            }
            v = node_to_idx[v];
            mask[u] = mask[v] = true;
            graph[u].push_back(v);
            graph[v].push_back(u);
        }

        unsigned long long mem = 0;
        int j = 0;
        for (int i = 0; i < cnt; ++i){
            double temp = double(++j)/cnt;

            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("%3.6lf%% done", temp*100);

            mem = max(mem, solve(i, cnt));
        }
        printf("\n");

        memory += mem;
        memory += getMemoryUsage(idx_to_node);
        memory += getMemoryUsage(node_to_idx);
        memory += getMemoryUsage2(graph, cnt);
        memory += getMemoryUsage3(betweenness, cnt);

        printf("For graph representation:\n");
        printf("Memory Used: %lld Bytes\n", memory);
        auto end = chrono::high_resolution_clock::now();
        double time = chrono::duration_cast<chrono::duration<double>>(end-begin).count();
        printf("Time Taken: %lf Seconds\n", time);

        double mx;
        fout.open("Output/Betweenness/" + names[i] + ".ungraph.betweenness.txt");
        fout<<memory<<" "<<fixed<<setprecision(9)<<time<<endl<<endl;
        for (auto& [i, j]: node_to_idx)
            fout<<i<<" "<<betweenness[j]<<endl;
        fout.close();

        mx = 0;
        for (int i = 0; i < cnt; ++i)
          mx = max(mx, betweenness[i]);

        printf("The node(s) with the greatest betweenness are ");
        for (int i = 0; i <= cnt; ++i)
          if (fabs(mx - betweenness[i]) < epsilon)
            printf("%d, ", idx_to_node[i]);
        printf("\b\b.\n");

        fin.close();






        // Hypergraph
        int hyperEdge = 0;
        fin.open(hypergraphs[i]);

        memory = 0;
        begin = chrono::high_resolution_clock::now();

        for (int i = 0; i < MAXN + MAXHE; ++i)
            graph[i].clear();
        idx_to_node.clear();
        node_to_idx.clear();
        betweenness.assign(MAXN + MAXHE, 0);

        string line;
        cnt = 0;
        mask.reset();
        while(getline(fin, line)){
            stringstream buffer(line);
            v = MAXN + hyperEdge;
            idx_to_node[cnt] = v;
            node_to_idx[v] = cnt++;
            v = node_to_idx[v];
            while(buffer>>u){
                if(!node_to_idx.count(u)){
                    idx_to_node[cnt] = u;
                    node_to_idx[u] = cnt++;
                }
                u = node_to_idx[u];
                mask.set(u);
                graph[u].push_back(v);
                graph[v].push_back(u);
          }
          ++hyperEdge;
        }

        mem = 0;
        j = 0;
        for (int i = 0; i < cnt; ++i){
            if (!mask[i])
                continue;

            double temp = double(++j)/(cnt-hyperEdge);

            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("%3.6lf%% done", temp*100);

            mem = max(mem, solve(i, cnt));
        }
        printf("\n");

        memory += mem;
        memory += getMemoryUsage(idx_to_node);
        memory += getMemoryUsage(node_to_idx);
        memory += getMemoryUsage2(graph, cnt);
        memory += getMemoryUsage3(betweenness, cnt);

        printf("For hypergraph representation:\n");
        printf("Memory Used: %lld Bytes\n", memory);
        end = chrono::high_resolution_clock::now();
        time = chrono::duration_cast<chrono::duration<double>>(end-begin).count();
        printf("Time Taken: %lf Seconds\n", time);

        fout.open("Output/Betweenness/" + names[i] + ".cmty.betweenness.txt");
        fout<<memory<<" "<<fixed<<setprecision(9)<<time<<endl<<endl;
        for (auto& [i, j]: node_to_idx)
            if (i < MAXN)
                fout<<i<<" "<<betweenness[j]<<endl;
        fout.close();

        mx = 0;
        for (int i = 0; i < cnt; ++i)
          mx = max(mx, mask[i]*betweenness[i]);

        printf("The node(s) with the greatest betweenness are ");
        for (int i = 0; i < cnt; ++i)
          if (mask[i] && fabs(mx - betweenness[i]) < epsilon)
            printf("%d, ", idx_to_node[i]);
        printf("\b\b.\n");

        fin.close();


        printf("%s\n", string(80, '=').c_str());
    }

    return 0;

}
