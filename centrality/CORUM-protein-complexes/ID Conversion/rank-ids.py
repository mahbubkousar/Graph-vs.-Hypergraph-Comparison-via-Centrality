import os

# --- 1. CONFIGURATION ---

# List of the GSEA-ready files you want to sort.
# This list is generated by the 'convert_to_gsea.py' script.
# Add or remove files from this list as needed.
FILES_TO_CONVERT = [
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Betweenness\corum_dataset.cmty.betweenness.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Betweenness\corum_dataset.ungraph.betweenness.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Degree\corum_dataset.cmty.degree.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Degree\corum_dataset.ungraph.degree.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\PageRank\corum_dataset.ungraph.pagerank.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\PageRank\corum_dataset.cmty.pagerank.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Closeness\corum_dataset.ungraph.closeness.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Closeness\corum_dataset.cmty.closeness.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Harmonic\corum_dataset.ungraph.harmonic.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Harmonic\corum_dataset.cmty.harmonic.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Farness\corum_dataset.cmty.farness.GSEA.rnk',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Farness\corum_dataset.ungraph.farness.GSEA.rnk'
]

# --- 2. PROCESS AND SORT EACH FILE ---

print("Starting to sort GSEA files...")
print("-" * 40)

for input_filepath in FILES_TO_CONVERT:
    
    # Define the new output filename. We'll add '.ranked' to it.
    # e.g., '...betweenness.GSEA.rnk' -> '...betweenness.GSEA.ranked.rnk'
    path_parts = os.path.splitext(input_filepath)
    output_filepath = f"{path_parts[0]}.ranked.rnk"

    print(f"  -> Processing '{input_filepath}'")

    try:
        data_to_sort = []
        with open(input_filepath, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) != 2:
                    print(f"      Warning: Skipping malformed line: '{line}'")
                    continue
                
                gene_symbol = parts[0]
                try:
                    # Convert score to a float for correct numeric sorting
                    score = float(parts[1])
                    data_to_sort.append([gene_symbol, score])
                except ValueError:
                    print(f"      Warning: Could not convert score to number on line: '{line}'")
        
        # --- The core sorting logic ---
        # We sort the list of lists based on the second element (the score),
        # in descending (reverse) order.
        data_to_sort.sort(key=lambda x: x[1], reverse=True)
        
        # --- Write the sorted data to the new file ---
        with open(output_filepath, 'w') as outfile:
            for item in data_to_sort:
                gene_symbol = item[0]
                score = item[1]
                outfile.write(f"{gene_symbol}\t{score}\n")

        print(f"      Success! Created sorted file: '{output_filepath}'")

    except FileNotFoundError:
        print(f"      ERROR: Input file not found: '{input_filepath}'. Skipping.")
        continue

print("-" * 40)
print("Ranking complete.")