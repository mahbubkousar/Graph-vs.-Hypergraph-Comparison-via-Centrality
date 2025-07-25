import os

# --- 1. CONFIGURATION ---

ORIGINAL_CORUM_FILE = 'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Data Preparation\gene_set_library_crisp.gmt'

FILES_TO_CONVERT = [
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Betweenness\corum_dataset.cmty.betweenness.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Betweenness\corum_dataset.ungraph.betweenness.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Degree\corum_dataset.cmty.degree.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Degree\corum_dataset.ungraph.degree.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\PageRank\corum_dataset.ungraph.pagerank.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\PageRank\corum_dataset.cmty.pagerank.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Closeness\corum_dataset.ungraph.closeness.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Closeness\corum_dataset.cmty.closeness.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Harmonic\corum_dataset.ungraph.harmonic.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Harmonic\corum_dataset.cmty.harmonic.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Farness\corum_dataset.cmty.farness.txt',
    'F:\Projects\graphs-and-hypergraphs\centrality\CORUM-protein-complexes\Output\Farness\corum_dataset.ungraph.farness.txt'
]

# --- 2. RE-GENERATE THE GENE MAPPING ---

print("Step 1: Re-generating the gene-to-ID mapping...")

gene_to_id = {}
next_gene_id = 0

try:
    with open(ORIGINAL_CORUM_FILE, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            # The same logic as the initial preparation script
            gene_symbols = parts[2:]
            for symbol in gene_symbols:
                if symbol not in gene_to_id:
                    gene_to_id[symbol] = next_gene_id
                    next_gene_id += 1
except FileNotFoundError:
    print(f"FATAL ERROR: The original CORUM file '{ORIGINAL_CORUM_FILE}' was not found.")
    print("This file is required to map integer IDs back to gene names.")
    exit()

# --- Create the reverse mapping (ID -> Gene) which we need for conversion ---
id_to_gene = {v: k for k, v in gene_to_id.items()}

print(f"Mapping successfully generated. Found {len(id_to_gene)} unique genes.")
print("-" * 40)


# --- 3. PROCESS EACH CENTRALITY FILE ---

print("Step 2: Converting centrality files to GSEA format...")

for input_filepath in FILES_TO_CONVERT:
    
    # Create the new output filename with a .GSEA.rnk extension
    # e.g., 'corum_dataset.ungraph.betweenness.GSEA.rnk'
    path_parts = os.path.splitext(input_filepath)
    output_filepath = f"{path_parts[0]}.GSEA.rnk"

    print(f"  -> Converting '{input_filepath}' to '{output_filepath}'")

    try:
        with open(input_filepath, 'r') as infile, open(output_filepath, 'w') as outfile:
            # The first two lines of the centrality files are header info (memory, time)
            # which we need to skip.
            infile.readline() # Skip memory/time line
            infile.readline() # Skip blank line

            for line in infile:
                line = line.strip()
                if not line:
                    continue
                
                # The line is formatted as: 'integer_id  centrality_score'
                parts = line.split()
                if len(parts) != 2:
                    continue

                try:
                    integer_id = int(parts[0])
                    score = parts[1]
                    
                    # Look up the integer ID to get the gene symbol
                    gene_symbol = id_to_gene[integer_id]
                    
                    # Write the new GSEA-compatible line to the output file
                    outfile.write(f"{gene_symbol}\t{score}\n")
                
                except (ValueError, KeyError) as e:
                    print(f"      Warning: Skipping malformed or unmappable line: '{line}' ({e})")

    except FileNotFoundError:
        print(f"      ERROR: Input file not found: '{input_filepath}'. Skipping.")
        continue

print("-" * 40)
print("Conversion complete. The new '.GSEA.rnk' files are ready for GSEA analysis.")