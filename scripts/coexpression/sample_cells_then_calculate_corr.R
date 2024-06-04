#!/usr/bin/env Rscript
#Date: May 20, 2024; May 23, 2024; June 4, 2024
#Author: Muyoung Lee
#Description: Collect various types of cells from the RDS file. Calculate correlations between the gene of interest and all the other genes.
#Usage: [THIS SCIRPT]

library(Seurat)
'%!in%' <- function(x,y)!('%in%'(x,y))

total <- readRDS("b84def55-a776-4aa4-a9a6-7aab8b973086.rds")
nCountRNA = "nCounts_RNA_UMIs"
pct_rank_cut = 80
zero_percent_cut = 30

gene_of_interest <- c("ENSG00000138073", "ENSG00000152700", "ENSG00000064835",
	"ENSG00000150527", "ENSG00000164081", "ENSG00000008018", "ENSG00000187535",
	"ENSG00000118965", "ENSG00000123607", "ENSG00000125686", "ENSG00000127463")
# PREB: ENSG00000138073
# SAR1B: ENSG00000152700
# PIT1 (POU1F1): ENSG00000064835
# MIA2: ENSG00000150527
# TEX264: ENSG00000164081
# PSMA1: ENSG00000129084 # Deleted
# PSMB1: ENSG00000008018
# IFT140: ENSG00000187535
# IFT121: ENSG00000118965
# IFT139: ENSG00000123607
# MED1: ENSG00000125686
# EMC1: ENSG00000127463

random <- read.csv("../random_gene_sets.txt", header=F, sep=",")[, 1]
target_genes <- c(gene_of_interest, random)

main_func <- function(annot, sample_size, my_obj){
	print(paste(annot, sample_size, sep=" "))

	Idents(my_obj) <- annot
	sub_obj <- subset(my_obj, downsample=sample_size)

	meta_file = paste(sample_size, annot, "meta", "txt", sep=".")
	write.table(sub_obj@meta.data, file=meta_file, quote=F, sep="\t")
	system(paste("../calculate_perc_per_annot.py", meta_file, nCountRNA, annot, sep=" "))
	perc_rank_file = paste(sample_size, annot, "perc_rank", "txt", sep=".")
	perc_ranks <- read.table(file=perc_rank_file, sep="\t", header=T)
	sub_obj <- AddMetaData(sub_obj, perc_ranks, col.name="pct_rank")
	sub_obj <- subset(sub_obj, subset= pct_rank <= pct_rank_cut)

	exp_matrix <- as.matrix(sub_obj[["RNA"]]$data)
	rm(sub_obj)

	total_genes <- rownames(exp_matrix)
	num_cells <- ncol(exp_matrix)
	num_zero_cells <- apply(exp_matrix, 1, function(x){sum(x==0)})
	output_genes = c()
	output_spearmans = c()

	for (target in target_genes){
		if (target %in% total_genes){
			target_exp <- as.numeric(exp_matrix[target, ])
			if (num_zero_cells[target] > (num_cells * zero_percent_cut / 100)){
				print(paste(target, " expression is zero more than ", zero_percent_cut, "% of cells in this subset", sep=""))
			}
			else{
				for (gene in total_genes){
					if (num_zero_cells[gene] <= (num_cells * zero_percent_cut / 100)){			
						spearman <- cor(target_exp, exp_matrix[gene, ], method="spearman")
						append(output_genes, gene)
						append(output_spearmans, spearman)
					}
				}
			}
		}
		else{
			print(paste(target, "is not in the data", sep=" "))
		}
	}

	output_df = data.frame(output_genes, output_spearmans)
	write.table(output_df, file="temp.txt", quote=F, sep="\t")
	output_name = paste(sample_size, annot, "spearman", "txt", sep=".")
	cmd <- paste("tail -n +2 temp.txt | sort -k2 -gr >", output_name, sep=" ") 
	print(cmd)
	system(cmd)
	#system("rm temp.txt")
}

### tissue_in_publication ###
main_func("tissue_in_publication", 2000, total)

###### tissue ######
subset2 <- subset(total, subset = tissue %!in% c("uterus", "lacrimal gland", "tongue"))
main_func("tissue", 500, subset2)

###### cell_type ######
blacklist <- c("type B pancreatic cell", "gut endothelial cell", "plasmacytoid dendritic cell", "bronchial smooth muscle cell", "tracheal goblet cell", "pancreatic A cell", "intestinal crypt stem cell of small intestine", "pancreatic PP cell", "mucus secreting cell", "tongue muscle cell", "surface ectodermal cell", "epithelial cell of lacrimal sac", "retinal pigment epithelial cell", "blood vessel endothelial cell", "salivary gland cell", "immature natural killer cell", "adipocyte", "duodenum glandular cell", "intrahepatic cholangiocyte", "serous cell of epithelium of trachea", "mesothelial cell", "Langerhans cell", "retinal bipolar neuron", "liver dendritic cell", "DN4 thymocyte", "ciliated epithelial cell", "erythroid lineage cell", "pulmonary ionocyte", "mature conventional dendritic cell", "pigmented ciliary epithelial cell", "serous cell of epithelium of bronchus", "plasmablast", "sperm", "cell of skeletal muscle", "Schwann cell", "pancreatic D cell", "double-positive, alpha-beta thymocyte", "retina horizontal cell", "retinal ganglion cell", "myeloid dendritic cell", "unknown")
subset3 <- subset(total, subset = cell_type %!in% blacklist)
main_func("cell_type", 100, subset3)

### cell_ontology_class ###
blacklist2 <- c("gut endothelial cell", "pancreatic beta cell", "plasmacytoid dendritic cell", "bronchial smooth muscle cell", "tracheal goblet cell", "limbal stromal cell", "intestinal crypt stem cell of small intestine", "pancreatic pp cell", "pancreatic alpha cell", "mucus secreting cell", "tongue muscle cell", "ocular surface cell", "epithelial cell of lacrimal sac", "retinal pigment epithelial cell", "bronchial vessel endothelial cell", "salivary gland cell", "lymphatic endothelial cell", "immature natural killer cell", "duodenum glandular cell", "adipocyte", "intrahepatic cholangiocyte", "serous cell of epithelium of trachea", "mesothelial cell", "langerhans cell", "liver dendritic cell", "retinal bipolar neuron", "dn4 thymocyte", "ciliated epithelial cell", "erythroid lineage cell", "pulmonary ionocyte", "ciliary body", "mature conventional dendritic cell", "serous cell of epithelium of bronchus", "plasmablast", "sperm", "cell of skeletal muscle", "schwann cell", "pancreatic delta cell", "respiratory mucous cell", "double-positive, alpha-beta thymocyte", "retina horizontal cell", "retinal ganglion cell", "myeloid dendritic cell")
subset4 <- subset(total, subset = cell_ontology_class %!in% blacklist2)
main_func("cell_ontology_class", 100, subset4)
