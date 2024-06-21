#!/usr/bin/env python3
#Date: June 18, 2024
#Author: Muyoung Lee
#Description: Draw ROC and PR curves of sorted gene list with well-known biological complexes
#Usage: [THIS SCRIPT] [SORTED GENE LIST]

import sys
from sklearn import metrics
import matplotlib.pyplot as plt

# Functions
def cal_TPR_FPR_AUROC(sorted_gene_list, positive_ctrl_set):
	total = len(sorted_gene_list)
	actual_positive = positive_ctrl_set & set(sorted_gene_list)
	dy_TPR = 1 / len(actual_positive)
	dx_FPR = 1 / (total - len(actual_positive))
	calculated_values = []
	x, y, AUC = 0, 0, 0
	for gene in sorted_gene_list:
		dAUC = 0
		if gene in actual_positive:
			y += 1
		else:
			x += 1
			dAUC = dx_FPR * y * dy_TPR
		AUC += dAUC
		calculated_values.append((x * dx_FPR, y * dy_TPR, AUC))
	remainder_AUC = (1 - x * dx_FPR) * (y * dy_TPR + 1) * 0.5
	total_AUC = AUC + remainder_AUC
	calculated_values.append((1.0, 1.0, total_AUC))
	return calculated_values

def cal_pre_rec_AUPRC(sorted_gene_list, positive_ctrl_set):
	total = len(sorted_gene_list)
	actual_positive = positive_ctrl_set & set(sorted_gene_list)
	# x: recall = TP / Actual Positive
	# y: precision = TP / Predicted Positive
	AP = len(actual_positive)
	PP = 0
	TP = 0
	calculated_values = []
	for gene in sorted_gene_list:
		PP += 1
		if gene in actual_positive:
			TP += 1
		calculated_values.append((TP/AP, TP/PP))
	recall = [x[0] for x in calculated_values]
	precision = [x[1] for x in calculated_values]
	AUC = metrics.auc(recall, precision)
	return (recall, precision, AUC)

def draw_curve_ROC(coordinates_for_graph, complex_name):
	x = []
	y = []
	for values in coordinates_for_graph:
		x.append(values[0])
		y.append(values[1])
	total_AUC = coordinates_for_graph[-1][-1]
	plt.plot(x, y, linewidth=1, alpha=0.6, label=f"{complex_name}_{total_AUC:0.3f}")
	plt.xlim([0, 1])
	plt.ylim([0, 1])

def draw_curve_PRC(coordinates_for_graph, complex_name):
	x = coordinates_for_graph[0]
	y = coordinates_for_graph[1]
	total_AUC = coordinates_for_graph[2]
	plt.plot(x, y, linewidth=1, alpha=0.6, label=f"{complex_name}_{total_AUC:0.3f}")
	plt.xlim([0, 1])
	plt.ylim([0, 1])

# Main
def main():
	geneset_dict = {"IFT": set(), "20S": set(), "MED": set(), "EMC": set()}
	with open("known_complexes.txt") as REF:
		REF.readline()
		for line in REF:
			nickname, name, EnsID, _, gene_type, chromosome = line.strip().split("\t")
			if gene_type == "protein_coding":
				if chromosome == "X" or chromosome.isdigit():
					if nickname.startswith("IFT"):
						geneset_dict["IFT"].add(EnsID)
					elif nickname.startswith("PSM"):
						geneset_dict["20S"].add(EnsID)

	with open("EMC_complex_genes.txt") as REF2:
		REF2.readline()
		for line in REF2:
			EnsID, _, gene_type, chromosome = line.strip().split("\t")
			if gene_type == "protein_coding":
				if chromosome == "X" or chromosome.isdigit():
					geneset_dict["EMC"].add(EnsID)

	with open("mediator_complex_genes.txt") as REF3:
		REF3.readline()
		for line in REF3:
			EnsID, _, gene_type, chromosome = line.strip().split("\t")
			if gene_type == "protein_coding":
				if chromosome == "X" or chromosome.isdigit():
					geneset_dict["MED"].add(EnsID)
	
	with open("random_gene_sets.txt") as RANDOM:
		count = 0
		for line in RANDOM:
			count += 1
			EnsID_list = line.strip().split(",")
			geneset_dict[f"RANDOM{count}"] = set(EnsID_list)
			
	# input gene list
	target_gene = sys.argv[1].split("/")[-1].split(".")[0]
	location = "/".join(sys.argv[1].split("/")[:-1])
	if target_gene.startswith("RANDOM"):
		target_gene += "_gene1"
	sorted_gene_list = []
	with open(sys.argv[1]) as INPUT:
		INPUT.readline()
		for line in INPUT:
			gene = line.strip().split("\t")[0]
			sorted_gene_list.append(gene)

	print("Query_gene\tGeneset\tAUROC")	
	for known_complex in geneset_dict:
		values = cal_TPR_FPR_AUROC(sorted_gene_list, geneset_dict[known_complex])
		draw_curve_ROC(values, known_complex)
		print(target_gene, known_complex, f"{values[-1][-1]:0.3f}", sep="\t")
	
	plt.plot([0,1], [0,1], color="gray", linewidth=1, linestyle="dashed")
	plt.legend(loc="lower right")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.title(f"When gene are sorted\n by corr with {target_gene}")
	plt.savefig(f"{location}/{target_gene}.ROC.png", dpi=300)
	plt.clf()	
	
	print("Query_gene\tGeneset\tAUPRC")	
	for known_complex in geneset_dict:
		values = cal_pre_rec_AUPRC(sorted_gene_list, geneset_dict[known_complex])
		draw_curve_PRC(values, known_complex)
		print(target_gene, known_complex, f"{values[-1]:0.3f}", sep="\t")
	
	plt.plot([0,1], [0,1], color="gray", linewidth=1, linestyle="dashed")
	plt.legend(loc="lower right")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.title(f"When gene are sorted\n by corr with {target_gene}")
	plt.savefig(f"{location}/{target_gene}.PR.png", dpi=300)
	plt.clf()
	
if __name__ == "__main__":
	main()