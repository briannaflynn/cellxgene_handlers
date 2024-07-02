#!/usr/bin/env python3
#Date: June 4, 2024
#Author: Muyoung Lee
#Description: Calculate percentile ranks for each annotation, like tissue, cell types. Then, make a table to be used as an metadata in Seurat package.
#Usage: [THIS SCRIPT] [metadata] [what's the name of the column for transcript counts?] [annotation]

import sys
import pandas as pd

col_nCounts = sys.argv[2]
col_annot = sys.argv[3]

df = pd.read_csv(sys.argv[1], sep="\t", header=0, index_col=0)
df["pct_rank"] = df.groupby(col_annot)[col_nCounts].rank(pct=True, ascending=False) * 100

output = sys.argv[1].replace("meta", "perc_rank")
df["pct_rank"].to_csv(output, sep="\t")
