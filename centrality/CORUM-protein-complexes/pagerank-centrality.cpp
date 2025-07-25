#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <ctype.h>
#include <limits>

using namespace std;

// --- Global constants and file management vectors (from our established structure) ---
const int MAXN = 4e6; // Used to distinguish protein nodes from hyperedge nodes
const double epsilon = 1e-9;
string init[] = {"dataset_init.txt"};
vector<string> dataset;
vector<string> graphs;
vector<string> hypergraphs;
vector<string> names;

// --- Provided PageRank Implementation (wrapped in a namespace) ---
namespace pageRank {
    // --- Default parameters for the PageRank algorithm ---
    const double DEFAULT_ALPHA = 0.85;
    const double DEFAULT_CONVERGENCE = 0.00001;
    const unsigned long DEFAULT_MAX_ITERATIONS = 10000;
    const bool DEFAULT_NUMERIC = true; // Our input files use numeric IDs
    const bool DEFAULT_UNDIRECTED = true;

    /*
     * A PageRank calculator class.
     * Original source: https://github.com/louridas/pagerank/tree/master
     * Modified to fit the project's needs (memory calculation, etc.).
     */
    class Table {
    private:
        double alpha;
        double convergence;
        unsigned long max_iterations;
        bool undirected;
        bool hypergraph;

        vector<vector<size_t>> rows;
        vector<size_t> num_outgoing;
        map<string, size_t> nodes_to_idx;
        map<size_t, string> idx_to_nodes;
        vector<double> pr;

        void trim(string &str) {
            size_t startpos = str.find_first_not_of(" \t");
            if (string::npos == startpos) {
                str = "";
            } else {
                str = str.substr(startpos, str.find_last_not_of(" \t") - startpos + 1);
            }
        }

        template <class Vector, class T>
        bool insert_into_vector(Vector& v, const T& t) {
            typename Vector::iterator i = lower_bound(v.begin(), v.end(), t);
            if (i == v.end() || t < *i) {
                v.insert(i, t);
                return true;
            }
            return false;
        }

        void reset() {
            num_outgoing.clear();
            rows.clear();
            nodes_to_idx.clear();
            idx_to_nodes.clear();
            pr.clear();
        }

        size_t insert_mapping(const string &key) {
            auto it = nodes_to_idx.find(key);
            if (it != nodes_to_idx.end()) {
                return it->second;
            }
            size_t index = nodes_to_idx.size();
            nodes_to_idx[key] = index;
            idx_to_nodes[index] = key;
            return index;
        }

        void add_arc(size_t from, size_t to) {
            size_t max_dim = max(from, to);
            if (rows.size() <= max_dim) {
                rows.resize(max_dim + 1);
                num_outgoing.resize(max_dim + 1);
            }
            if (insert_into_vector(rows[to], from)) {
                num_outgoing[from]++;
            }
        }

    public:
        Table(double a = DEFAULT_ALPHA, double c = DEFAULT_CONVERGENCE, size_t i = DEFAULT_MAX_ITERATIONS)
            : alpha(a), convergence(c), max_iterations(i), undirected(DEFAULT_UNDIRECTED), hypergraph(false) {}

        void set_undirected(bool u) { undirected = u; }
        void set_hypergraph(bool h) { hypergraph = h; }

        int read_file(const string &filename) {
            reset();
            ifstream infile(filename.c_str());
            if (!infile) {
                cerr << "Error: Cannot open file " << filename << endl;
                return 1;
            }

            string line;
            if (!hypergraph) {
                getline(infile, line); // Skip header line in .graph files
            }

            size_t linenum = 0;
            while (getline(infile, line)) {
                if (line.empty()) continue;
                stringstream buffer(line);

                if (hypergraph) {
                    string hyperedge_node_str = to_string(MAXN + linenum);
                    size_t hyperedge_idx = insert_mapping(hyperedge_node_str);
                    string protein_node_str;
                    while (buffer >> protein_node_str) {
                        size_t protein_idx = insert_mapping(protein_node_str);
                        add_arc(protein_idx, hyperedge_idx);
                        add_arc(hyperedge_idx, protein_idx);
                    }
                } else {
                    string from_str, to_str;
                    buffer >> from_str >> to_str;
                    size_t from_idx = insert_mapping(from_str);
                    size_t to_idx = insert_mapping(to_str);
                    add_arc(from_idx, to_idx);
                    if (undirected) {
                        add_arc(to_idx, from_idx);
                    }
                }
                linenum++;
            }
            return 0;
        }

        void calculate_pagerank() {
            size_t num_rows = rows.size();
            if (num_rows == 0) return;

            pr.assign(num_rows, 1.0 / num_rows);
            vector<double> old_pr;
            double diff = 1;
            unsigned long num_iterations = 0;

            while (diff > convergence && num_iterations < max_iterations) {
                old_pr = pr;
                double dangling_pr_sum = 0;
                for (size_t k = 0; k < num_rows; ++k) {
                    if (num_outgoing[k] == 0) {
                        dangling_pr_sum += old_pr[k];
                    }
                }

                for (size_t i = 0; i < num_rows; ++i) {
                    double h = 0.0;
                    for (size_t incoming_node : rows[i]) {
                        h += old_pr[incoming_node] / num_outgoing[incoming_node];
                    }
                    pr[i] = alpha * (h + dangling_pr_sum / num_rows) + (1.0 - alpha) / num_rows;
                }

                diff = 0;
                for (size_t i = 0; i < num_rows; ++i) {
                    diff += fabs(pr[i] - old_pr[i]);
                }
                num_iterations++;
            }
        }
        
        // --- Accessor methods to get results ---
        const vector<double>& get_pagerank_vector() const { return pr; }
        const map<size_t, string>& get_idx_to_node_map() const { return idx_to_nodes; }

        // --- New method to calculate memory usage of this class instance ---
        unsigned long long getMemoryUsage() const {
            unsigned long long totalSize = 0;
            // Memory for vector<vector<size_t>> rows
            for(const auto& row : rows) {
                totalSize += row.capacity() * sizeof(size_t);
            }
            totalSize += rows.capacity() * sizeof(vector<size_t>);
            // Other vectors
            totalSize += num_outgoing.capacity() * sizeof(size_t);
            totalSize += pr.capacity() * sizeof(double);
            // Memory for maps (estimation)
            totalSize += nodes_to_idx.size() * (sizeof(string) + sizeof(size_t) + 3 * sizeof(void*));
            totalSize += idx_to_nodes.size() * (sizeof(size_t) + sizeof(string) + 3 * sizeof(void*));
            return totalSize;
        }
    };
} // namespace pageRank

int main() {
    string output_dir_name = "Output/PageRank/";
    cout << "INFO: This script will save results to '" << output_dir_name << "'." << endl;
    cout << "INFO: Please ensure this directory exists." << endl;

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    ifstream fin;
    ofstream fout;

    for (const auto& s : init) {
        fin.open(s); // Uses relative path
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
            pageRank::Table calculator;
            calculator.set_undirected(true);
            calculator.set_hypergraph(false);

            auto begin = chrono::high_resolution_clock::now();
            calculator.read_file(graphs[i]); // Using relative path from init file
            calculator.calculate_pagerank();
            auto end = chrono::high_resolution_clock::now();
            double time = chrono::duration_cast<chrono::duration<double>>(end - begin).count();
            
            unsigned long long memory = calculator.getMemoryUsage();

            cout << "For graph representation:" << endl;
            cout << "Memory Used: " << memory << " Bytes" << endl;
            cout << "Time Taken: " << fixed << setprecision(9) << time << " Seconds" << endl;

            string out_fname = output_dir_name + names[i] + ".ungraph.pagerank.txt";
            fout.open(out_fname);
            fout << memory << " " << fixed << setprecision(9) << time << endl << endl;

            const auto& ranks = calculator.get_pagerank_vector();
            const auto& mapping = calculator.get_idx_to_node_map();
            double max_rank = 0.0;
            
            for (size_t j = 0; j < ranks.size(); ++j) {
                const string& node_name = mapping.at(j);
                fout << node_name << " " << ranks[j] << endl;
                if (ranks[j] > max_rank) {
                    max_rank = ranks[j];
                }
            }
            fout.close();

            cout << "The node(s) with the greatest PageRank are ";
            bool first = true;
            for (size_t j = 0; j < ranks.size(); ++j) {
                if (fabs(ranks[j] - max_rank) < epsilon) {
                    if (!first) cout << ", ";
                    cout << mapping.at(j);
                    first = false;
                }
            }
            cout << ".\n" << endl;
        }
        
        // ===================================================================
        //  Part 2: Hypergraph Representation
        // ===================================================================
        {
            pageRank::Table calculator;
            calculator.set_undirected(true); // Our bipartite model is undirected
            calculator.set_hypergraph(true);

            auto begin = chrono::high_resolution_clock::now();
            calculator.read_file(hypergraphs[i]);
            calculator.calculate_pagerank();
            auto end = chrono::high_resolution_clock::now();
            double time = chrono::duration_cast<chrono::duration<double>>(end - begin).count();

            unsigned long long memory = calculator.getMemoryUsage();
            
            cout << "For hypergraph representation:" << endl;
            cout << "Memory Used: " << memory << " Bytes" << endl;
            cout << "Time Taken: " << fixed << setprecision(9) << time << " Seconds" << endl;

            string out_fname = output_dir_name + names[i] + ".cmty.pagerank.txt";
            fout.open(out_fname);
            fout << memory << " " << fixed << setprecision(9) << time << endl << endl;

            const auto& ranks = calculator.get_pagerank_vector();
            const auto& mapping = calculator.get_idx_to_node_map();
            double max_rank = 0.0;

            for (size_t j = 0; j < ranks.size(); ++j) {
                const string& node_name = mapping.at(j);
                // Only consider original protein nodes for output and max rank
                if (stoll(node_name) < MAXN) {
                    fout << node_name << " " << ranks[j] << endl;
                    if (ranks[j] > max_rank) {
                        max_rank = ranks[j];
                    }
                }
            }
            fout.close();

            cout << "The node(s) with the greatest PageRank are ";
            bool first = true;
            for (size_t j = 0; j < ranks.size(); ++j) {
                 const string& node_name = mapping.at(j);
                if (stoll(node_name) < MAXN && fabs(ranks[j] - max_rank) < epsilon) {
                    if (!first) cout << ", ";
                    cout << node_name;
                    first = false;
                }
            }
            cout << ".\n" << endl;
        }

        cout << string(80, '=') << endl;
    }

    return 0;
}