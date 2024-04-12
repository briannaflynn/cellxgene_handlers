# Aggregate data processing and normalization

The following tables are samples, full size tables available on the Marcotte pods from the CELLxGENE directory in SCRATCH here: ```CELLxGENE/collections/aggregated_summary_tables/```. See cellxgene manuscript and explorer documentation for additional details. 

**Removal of Duplicate Cells**

Some data on CellxGene is duplicated due to independent submissions, for example metaanalysis vs original data. All data submitted on Discover is curated to indicate whether any cell
is the primary data. Only cells demarcated as primary data are included in the processing steps
below.

**Removal of Cells Based on Sequencing Assay**

Only cells from sequencing assays that measure gene expression and don't require gene-length
normalization are included ( see table 1 of cellxgene manuscript )

**Data Normalization**

Read counts are normalized using the ln(CPTT+1) transformation of raw counts.
Normalized matrices from multiple datasets of the same tissue are concatenated along the gene
axis.

```
import numpy as np

# Example raw read counts for a single gene across multiple cells
raw_read_counts = np.array([100, 500, 0, 200])

# Total counts per cell for each of the cells (example values)
total_counts_per_cell = np.array([10000, 20000, 5000, 15000])

# Calculate CPTT (counts per ten thousand)
cptt = (raw_read_counts / total_counts_per_cell) * 10000

# Apply the log transformation
transformed_counts = np.log(cptt + 1)

print("Transformed counts:", transformed_counts)
```

**Removal of Noisy Ultra-low Expression Values**

After applying normalization, any gene/cell combination counts less or equal than 2 are set to
missing data. This allows for removal of noise due to ultra-lowly expressed genes and provides
a cleaner visualization.

_____________________________________________________________________________________________

# Note about cell ontology

The cell ontology is a hierarchical tree structure that represents the relationship between cell
types. For example, the cell type "B cell" is a descendant of the cell type "lymphocyte". For a
particular cell type, the cell ontology is used to sum up the expression values and cell counts of
cells labeled as that cell type as well as those labeled as its descendants. In the aforementioned
example, the average expression of "lymphocyte" would include "B cells" and all its other
descendants.
This rollup operation accounts for the fact that different datasets may have labeled their cells
with different levels of granularity. It provides a more robust measure of the average expression
of low-granularity cell type terms, such as "secretory cell" or "lymphocyte".
