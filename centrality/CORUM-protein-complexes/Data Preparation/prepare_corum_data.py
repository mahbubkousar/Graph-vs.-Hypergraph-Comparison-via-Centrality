import os
import itertools

# --- Configuration ---
# Input file from the CORUM database
CORUM_INPUT_FILE = 'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Data Preparation\gene_set_library_crisp.gmt'

# Base name for the output files
DATASET_BASE_NAME = 'corum_dataset'

# --- Output filenames with desired .graph and .hypergraph extensions ---
GRAPH_OUTPUT_FILE = f'{DATASET_BASE_NAME}.graph'
HYPERGRAPH_OUTPUT_FILE = f'{DATASET_BASE_NAME}.hypergraph'
DATASET_INIT_FILE = 'dataset_init.txt' # Configuration file for your C++ program

# The C++ code saves results here, so the script ensures it exists
OUTPUT_DIR_FOR_C_CODE = 'Output/Betweenness'

# --- Data Structures to hold processed data ---
gene_to_id = {}
next_gene_id = 0
all_complexes_gene_ids = []

# --- 1. Parse the Gene Set Library and Create Gene-ID Mapping ---
print(f"Reading input file: {CORUM_INPUT_FILE}...")
try:
    with open(CORUM_INPUT_FILE, 'r') as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                print(f"Warning: Skipping line {line_num + 1} with fewer than 3 columns.")
                continue

            current_complex_gene_symbols = parts[2:]
            current_complex_gene_ids = []

            for gene_symbol in current_complex_gene_symbols:
                if gene_symbol not in gene_to_id:
                    gene_to_id[gene_symbol] = next_gene_id
                    next_gene_id += 1
                current_complex_gene_ids.append(gene_to_id[gene_symbol])

            if current_complex_gene_ids:
                all_complexes_gene_ids.append(current_complex_gene_ids)
except FileNotFoundError:
    print(f"ERROR: The input file '{CORUM_INPUT_FILE}' was not found.")
    print("Please make sure the file is in the same directory as the script.")
    exit()


print(f"Processed {len(all_complexes_gene_ids)} complexes and found {len(gene_to_id)} unique genes.")

# --- 2. Generate .graph File (Pairwise Edges) ---
print(f"Generating graph file: {GRAPH_OUTPUT_FILE}...")
pairwise_edges = set()
for complex_gene_ids in all_complexes_gene_ids:
    if len(complex_gene_ids) < 2:
        continue
    for u, v in itertools.combinations(complex_gene_ids, 2):
        pairwise_edges.add(tuple(sorted((u, v))))

num_nodes_graph = len(gene_to_id)
num_edges_graph = len(pairwise_edges)

with open(GRAPH_OUTPUT_FILE, 'w') as f:
    f.write(f"{num_nodes_graph} {num_edges_graph}\n")
    for u, v in sorted(list(pairwise_edges)):
        f.write(f"{u} {v}\n")

print(f"Graph file generated: {num_nodes_graph} nodes, {num_edges_graph} edges.")

# --- 3. Generate .hypergraph File (Complex Membership) ---
print(f"Generating hypergraph file: {HYPERGRAPH_OUTPUT_FILE}...")
with open(HYPERGRAPH_OUTPUT_FILE, 'w') as f:
    for complex_gene_ids in all_complexes_gene_ids:
        f.write(" ".join(map(str, complex_gene_ids)) + "\n")

print(f"Hypergraph file generated: {len(all_complexes_gene_ids)} hyperedges.")


# --- 4. Create dataset_init.txt for C++ program ---
print(f"Creating C++ initialization file: {DATASET_INIT_FILE}...")
os.makedirs(OUTPUT_DIR_FOR_C_CODE, exist_ok=True)

with open(DATASET_INIT_FILE, 'w') as f:
    f.write("CORUM Protein Complexes Dataset\n")
    # Point the C++ program to the new filenames
    f.write(f"{GRAPH_OUTPUT_FILE}\n")
    f.write(f"{HYPERGRAPH_OUTPUT_FILE}\n")
    # This is the base name your C++ program will use for its output files
    f.write(f"{DATASET_BASE_NAME}\n")

print(f"All required files for the C++ program have been generated.")
print("\n--- Next Steps ---")
print(f"1. Place '{GRAPH_OUTPUT_FILE}', '{HYPERGRAPH_OUTPUT_FILE}', and '{DATASET_INIT_FILE}' in the same directory as your C++ executable.")
print("2. Compile and run your C++ code. It will now read the correctly named files.")