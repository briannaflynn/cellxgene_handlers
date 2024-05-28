#!/usr/bin/env Rscript
#Date: May 20, 2024; May 23, 2024
#Author: Muyoung Lee
#Description: Collect various types of cells from the RDS file. Calculate correlations between the gene of interest and all the other genes.
#Usage: [THIS SCIRPT]

library(Seurat)
'%!in%' <- function(x,y)!('%in%'(x,y))

total <- readRDS("b84def55-a776-4aa4-a9a6-7aab8b973086.rds")

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

### tissue_in_publication ###
Idents(total) <- "tissue_in_publication"
sub_obj <- subset(total, downsample=2000)
write.table(sub_obj@meta.data, file="2000_cells_per_tissue_in_pub.meta.txt", quote=F, sep="\t")

exp_matrix <- as.matrix(sub_obj[["RNA"]]$data)
rm(sub_obj)
total_genes <- rownames(exp_matrix)

for (gene in target_genes){
	if (gene %in% total_genes){
		print(gene)
		target <- as.numeric(exp_matrix[gene, ])
		spearman <- apply(exp_matrix, 1, function(x){cor(target, x, method="spearman")})
		write.table(spearman, file="temp.txt", quote=F, sep="\t")
		cmd <- paste("tail -n +2 temp.txt | sort -k2 -gr > human_2000_cells_per_tissue_in_pub.", gene, "_spearman.txt", sep="") 
		system(cmd)
		system("rm temp.txt")
	}
}

### tissue ###
Idents(total) <- "tissue"
sub_obj <- subset(total, subset = tissue %!in% c("uterus", "lacrimal gland", "tongue"))
sub_obj <- subset(sub_obj, downsample=500)
write.table(sub_obj@meta.data, file="500_cells_per_tissue.meta.txt", quote=F, sep="\t")

exp_matrix <- as.matrix(sub_obj[["RNA"]]$data)
rm(sub_obj)
total_genes <- rownames(exp_matrix)

for (gene in target_genes){
	if (gene %in% total_genes){
		print(gene)
		target <- as.numeric(exp_matrix[gene, ])
		spearman <- apply(exp_matrix, 1, function(x){cor(target, x, method="spearman")})
		write.table(spearman, file="temp.txt", quote=F, sep="\t")
		cmd <- paste("tail -n +2 temp.txt | sort -k2 -gr > human_500_cells_per_tissue.", gene, "_spearman.txt", sep="") 
		system(cmd)
		system("rm temp.txt")
	}
}

### cell_type ###
blacklist <- c("type B pancreatic cell", "gut endothelial cell", "plasmacytoid dendritic cell", "bronchial smooth muscle cell", "tracheal goblet cell", "pancreatic A cell", "intestinal crypt stem cell of small intestine", "pancreatic PP cell", "mucus secreting cell", "tongue muscle cell", "surface ectodermal cell", "epithelial cell of lacrimal sac", "retinal pigment epithelial cell", "blood vessel endothelial cell", "salivary gland cell", "immature natural killer cell", "adipocyte", "duodenum glandular cell", "intrahepatic cholangiocyte", "serous cell of epithelium of trachea", "mesothelial cell", "Langerhans cell", "retinal bipolar neuron", "liver dendritic cell", "DN4 thymocyte", "ciliated epithelial cell", "erythroid lineage cell", "pulmonary ionocyte", "mature conventional dendritic cell", "pigmented ciliary epithelial cell", "serous cell of epithelium of bronchus", "plasmablast", "sperm", "cell of skeletal muscle", "Schwann cell", "pancreatic D cell", "double-positive, alpha-beta thymocyte", "retina horizontal cell", "retinal ganglion cell", "myeloid dendritic cell", "unknown")

Idents(total) <- "cell_type"
sub_obj <- subset(total, subset = cell_type %!in% blacklist)
sub_obj <- subset(sub_obj, downsample=100)
write.table(sub_obj@meta.data, file="100_cells_per_cell_type.meta.txt", quote=F, sep="\t")

exp_matrix <- as.matrix(sub_obj[["RNA"]]$data)
rm(sub_obj)
total_genes <- rownames(exp_matrix)

for (gene in target_genes){
	if (gene %in% total_genes){
		print(gene)
		target <- as.numeric(exp_matrix[gene, ])
		spearman <- apply(exp_matrix, 1, function(x){cor(target, x, method="spearman")})
		write.table(spearman, file="temp.txt", quote=F, sep="\t")
		cmd <- paste("tail -n +2 temp.txt | sort -k2 -gr > human_100_cells_per_cell_type.", gene, "_spearman.txt", sep="") 
		system(cmd)
		system("rm temp.txt")
	}
}

### cell_ontology_class ###
blacklist2 <- c("gut endothelial cell", "pancreatic beta cell", "plasmacytoid dendritic cell", "bronchial smooth muscle cell", "tracheal goblet cell", "limbal stromal cell", "intestinal crypt stem cell of small intestine", "pancreatic pp cell", "pancreatic alpha cell", "mucus secreting cell", "tongue muscle cell", "ocular surface cell", "epithelial cell of lacrimal sac", "retinal pigment epithelial cell", "bronchial vessel endothelial cell", "salivary gland cell", "lymphatic endothelial cell", "immature natural killer cell", "duodenum glandular cell", "adipocyte", "intrahepatic cholangiocyte", "serous cell of epithelium of trachea", "mesothelial cell", "langerhans cell", "liver dendritic cell", "retinal bipolar neuron", "dn4 thymocyte", "ciliated epithelial cell", "erythroid lineage cell", "pulmonary ionocyte", "ciliary body", "mature conventional dendritic cell", "serous cell of epithelium of bronchus", "plasmablast", "sperm", "cell of skeletal muscle", "schwann cell", "pancreatic delta cell", "respiratory mucous cell", "double-positive, alpha-beta thymocyte", "retina horizontal cell", "retinal ganglion cell", "myeloid dendritic cell")

Idents(total) <- "cell_ontology_class"
sub_obj <- subset(total, subset = cell_ontology_class %!in% blacklist2)
sub_obj <- subset(sub_obj, downsample=100)
write.table(sub_obj@meta.data, file="100_cells_per_cell_ontology.meta.txt",
			quote=F, sep="\t")

exp_matrix <- as.matrix(sub_obj[["RNA"]]$data)
rm(sub_obj)
total_genes <- rownames(exp_matrix)

for (gene in target_genes){
	if (gene %in% total_genes){
		print(gene)
		target <- as.numeric(exp_matrix[gene, ])
		spearman <- apply(exp_matrix, 1, function(x){cor(target, x, method="spearman")})
		write.table(spearman, file="temp.txt", quote=F, sep="\t")
		cmd <- paste("tail -n +2 temp.txt | sort -k2 -gr > human_100_cells_per_cell_ontology.", gene, "_spearman.txt", sep="") 
		system(cmd)
		system("rm temp.txt")
	}
}
