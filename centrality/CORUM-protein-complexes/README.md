# Centrality Measures of Protein Complex using Graph and Hypergraph Representations

This project analyzes the structure of protein complexes from the CORUM dataset by calculating betweenness centrality on both graph and hypergraph models.

## Dataset: CORUM Protein Complexes

- **Link:** [The CORUM Database](http://mips.helmholtz-muenchen.de/corum/)
- **Citation:** Ruepp, A et al. (2010) CORUM: the comprehensive resource of mammalian protein complexes--2009. Nucleic Acids Res. 38:D497-501.

### Description

The **CORUM (the comprehensive resource of mammalian protein complexes)** database is a collection of experimentally verified protein complexes from mammalian organisms. It is manually curated from scientific literature by expert biologists. This high-quality dataset provides information on which proteins associate to form functional units (complexes).

For this project, we use the **Gene Set Library** file (`gene_set_library_crisp.gmt`), which lists each protein complex and its constituent member proteins (represented by their gene symbols).

## Data Conversion for Analysis

To analyze the network structure, the raw `gene_set_library_crisp.gmt` file is processed into two distinct representations using the `prepare_corum_data.py` script. This script maps each unique gene symbol to an integer ID to prepare it for the C++ analysis program.

### 1. Graph Representation (`.graph` file)

The standard graph model represents **pairwise interactions**. We assume that if two proteins are part of the same complex, they interact.

- **Nodes:** Each unique protein (gene) is a node in the graph.
- **Edges:** An edge is created between every pair of proteins that appear together in the same complex. For a complex with N proteins, this process generates N\*(N-1)/2 edges, forming a "clique".
- **Format (`corum_dataset.graph`):**
  1.  The first line contains the total number of unique nodes and the total number of unique edges.
  2.  Each subsequent line represents an edge, containing the integer IDs of the two connected nodes, separated by a space (`node1 node2`).

This representation is useful for standard network algorithms but loses the higher-order information of the multi-protein complex.

### 2. Hypergraph Representation (`.hypergraph` file)

The hypergraph model is a more natural fit for protein complex data, as it can represent interactions involving more than two members simultaneously.

- **Vertices (Nodes):** Each unique protein (gene) is a vertex.
- **Hyperedges:** Each protein complex itself is a single **hyperedge**. A hyperedge is simply the set of all vertices (proteins) that belong to that complex.
- **Format (`corum_dataset.hypergraph`):**
  1.  The file is a direct translation of the original dataset.
  2.  Each line represents a single hyperedge (a complex).
  3.  The line contains the integer IDs of all member proteins, separated by spaces (`node1 node2 node3 ...`).

This representation is used by our C++ program, which internally converts it to a bipartite graph (protein nodes vs. complex nodes) to calculate centrality, preserving the group-level information.
