# Centrality Calculation

This folder contains centrality calculation algorithms for graphs and hypergraphs, specifically designed to analyze protein complex networks from the CORUM database.

## Acknowledgments

The centrality calculation codes in this directory were written by **Azim** and **Sahil**.

## Project Structure

### CORUM-protein-complexes/
Main implementation for protein complex analysis using CORUM database:

#### Data Preparation/
- `prepare_corum_data.py` - Converts CORUM GMT files to graph/hypergraph formats
- `gene_set_library_crisp.gmt` - Input gene set library from CORUM database

#### ID Conversion/
- `convert-id.py` - Converts numerical IDs back to gene symbols for interpretability
- `rank-ids.py` - Ranks genes by centrality scores for GSEA analysis

#### Core Algorithms (C++)
- `betweenness-centrality.cpp` - Betweenness centrality using Brandes algorithm
- `degree-centrality.cpp` - Degree centrality calculation
- `pagerank-centrality.cpp` - PageRank centrality implementation
- `closeness-farness-harmonic-centrality.cpp` - Multiple distance-based centrality measures

#### Output/
Organized results by centrality measure (Betweenness, Closeness, Degree, Farness, Harmonic, PageRank):
- `.txt` files - Raw centrality scores
- `.rnk` files - Gene symbols with scores for GSEA
- `.ranked.rnk` files - Ranked gene lists for enrichment analysis
- Separate results for community (`cmty`) and ungraphed versions

### testing/
Contains test implementations with smaller social network datasets for algorithm validation.

## Centrality Measures

- **Betweenness** - Measures node importance based on shortest path betweenness using Brandes algorithm
- **Degree** - Simple node degree centrality for both graphs and hypergraphs
- **PageRank** - Google's PageRank algorithm adapted for protein networks
- **Closeness** - Average distance to all other nodes in the network
- **Farness** - Sum of distances to all other nodes (inverse of closeness)
- **Harmonic** - Harmonic mean of distances to other nodes

## Workflow

1. **Data Preparation**: Convert CORUM GMT files to graph/hypergraph format
2. **Centrality Calculation**: Run C++ algorithms on processed data
3. **ID Conversion**: Convert results back to gene symbols
4. **GSEA Preparation**: Generate ranked files for Gene Set Enrichment Analysis

## Graph Types

- **Standard Graph** (`ungraph`): Traditional protein-protein interaction networks
- **Community Graph** (`cmty`): Networks with community structure preserved
