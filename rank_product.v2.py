#!/usr/bin/env python3
#Date: May 21, 2024
#Author: Muyoung Lee
#Description: Calculate the rank product of genes based on Spearman correlations from each case.
#Usage: [THIS SCRIPT]

import glob
import pandas as pd

target_dict = {"ENSG00000138073": "PREB",
				"ENSG00000152700": "SAR1B",
				"ENSG00000064835": "PIT1",
				"ENSG00000150527": "MIA2",
				"ENSG00000164081": "TEX264",
				"ENSG00000008018": "PSMB1",
				"ENSG00000187535": "IFT140",
				"ENSG00000118965": "IFT121",
				"ENSG00000123607": "IFT139",
				"ENSG00000125686": "MED1",
				"ENSG00000127463": "EMC1"}

with open("random_gene_sets.txt") as RANDOM:
	count = 0
	for line in RANDOM:
		count += 1
		gene1 = line.strip().split(",")[0]
		target_dict[gene1] = f"Random_set{count}_gene1"
	
for target in target_dict:
	file_list = glob.glob(f"human_cell_landscape/*.{target}_spearman.txt")
	file_list += glob.glob(f"Tabula_Sapiens/*.{target}_spearman.txt")
	file_list += glob.glob(f"human_cell_atlas_of_fetal_gene_exp/*.{target}_spearman.txt")
	if len(file_list) < 6:
		continue
	flag_gene_rank = {}
	gene_flag = {}
	for file_ in file_list:
		collection= file_.split("/")[0]
		sampling = file_.split("/")[1].split(".")[0]
		flag = f"{sampling}@{collection}"
		if flag not in flag_gene_rank:
			flag_gene_rank[flag] = {}
		df = pd.read_csv(file_, sep="\t", names=["Gene", "Corr"], skiprows=1, index_col=0)
		df["pct_rank"] = df["Corr"].rank(pct=True, ascending=False)
		for gene in df.index:
			pct_rank = df.loc[gene, "pct_rank"]
			if not pd.isna(pct_rank):
				flag_gene_rank[flag][gene] = pct_rank * 100
				if gene not in gene_flag:
					gene_flag[gene] = []
				gene_flag[gene].append(flag)
	
	gene_RP = {}
	for gene in gene_flag:
	#	if len(gene_flag[gene]) > (len(file_list)/2):
		if len(gene_flag[gene]) > 5:
			rank_product = 1
			for flag in gene_flag[gene]:
				rank_product *= flag_gene_rank[flag][gene]
			rank_product = rank_product ** (1/len(gene_flag[gene]))
			gene_RP[gene] = rank_product
	
	gene_RP = dict(sorted(gene_RP.items(), key=lambda item: item[1]))

	OUT = open(f"{target_dict[target]}.co-expressed_genes.rank_product.tsv", "w")
	print("gene", "RP_of_perctile_ranks", "n(Subsamples)", sep="\t", file=OUT)
	for gene in gene_RP:
		print(gene, gene_RP[gene], len(gene_flag[gene]), sep="\t", file=OUT)
	OUT.close()
