import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
from collections import Counter
from typing import Optional, Union, List, Dict
from datetime import datetime
from synthetic_data_generator import generate_synthetic_anndata
from cxg_logger import *
from adata_preprocessor import Preprocessor
from anndata import AnnData

class WGCNA:
    def __init__(
        self,
        data: AnnData,
        power: int = 6,
        num_modules: int = 5,
    ):
        self.data: AnnData = data
        self.power: int = power
        self.num_modules: int = num_modules
        self.correlation_matrix: Optional[pd.DataFrame] = None
        self.adjacency_matrix: Optional[np.ndarray] = None
        self.TOM: Optional[np.ndarray] = None
        self.dissimilarity_TOM: Optional[np.ndarray] = None
        self.linkage_matrix: Optional[np.ndarray] = None
        self.modules: Optional[np.ndarray] = None

        self.dataframe: pd.DataFrame = None

    def export_to_dataframe(self) -> pd.DataFrame:
        # Convert sparse matrix to dense
        data_dense = self.data.X.toarray()
        
        df = pd.DataFrame(data_dense, index=self.data.obs['cell_id'], columns=self.data.var['gene_id'])
        
        # Add observation metadata 
        for col in self.data.obs.columns:
            df[col] = self.data.obs[col].values
        
        return df

    def preprocess_data(self) -> None:
        
        filter_preprocessor = Preprocessor(
            use_key=None,
            filter_gene_by_counts=10,
            gene_filter=False,
            filter_cell_by_counts=200,
            normalize_total=True,
            result_normed_key="X_normed",
            log1p=True,
            result_log1p_key="X_log1p",
            subset_hvg=True,
            hvg_use_key=None,
            hvg_flavor="seurat_v3",
            execute_filter_genes='on',
            execute_filter_cells='on',
            execute_normalize_total='on',
            execute_log1p='on',
            execute_subset_hvg='on',
            execute_binning='off'
        )

        # writes the changes in place for memory efficiency
        filter_preprocessor(self.data)

        self.dataframe = self.export_to_dataframe()

        # Outlier detection
        # commenting out for now ...

        # distance_matrix = pairwise_distances(self.dataframe, metric='euclidean')
        # linkage_matrix = linkage(squareform(distance_matrix), method='average')
        # max_dist = np.percentile(linkage_matrix[:, 2], 90)
        # clusters = fcluster(linkage_matrix, t=max_dist, criterion='distance')
        # cluster_counts = Counter(clusters)
        # outlier_clusters = [cluster for cluster, count in cluster_counts.items() if count < (len(self.data) / 100)]
        # outliers = [i for i, cluster in enumerate(clusters) if cluster in outlier_clusters]
        # self.dataframe = self.dataframe.drop(self.dataframe.index[outliers])

    def compute_correlation_matrix(self):
        self.correlation_matrix = self.dataframe.corr(method='spearmann')

    def compute_adjacency_matrix(self):
        def adjacency_function(correlation_matrix, power):
            return np.power(correlation_matrix, power)
        
        self.adjacency_matrix = adjacency_function(self.correlation_matrix, self.power)

    def compute_TOM(self):
        def TOM_similarity(adjacency_matrix):
            N = adjacency_matrix.shape[0]
            TOM = np.zeros((N, N))
            for i in range(N):
                for j in range(N):
                    TOM[i, j] = (np.sum(np.minimum(adjacency_matrix[i], adjacency_matrix[j])) + adjacency_matrix[i, j]) / (
                        np.sum(adjacency_matrix[i]) + np.sum(adjacency_matrix[j]) - adjacency_matrix[i, j]
                    )
            return TOM
        
        self.TOM = TOM_similarity(self.adjacency_matrix)
        self.dissimilarity_TOM = 1 - self.TOM

    def perform_clustering(self):
        self.linkage_matrix = linkage(squareform(self.dissimilarity_TOM), method='average')

    def identify_modules(self):
        def cut_tree_to_modules(linkage_matrix, num_modules):
            clusters = cut_tree(linkage_matrix, n_clusters=num_modules)
            return clusters.reshape(-1)

        self.modules = cut_tree_to_modules(self.linkage_matrix, self.num_modules)
        self.dataframe['Module'] = self.modules

    def save_results(self, filename='gene_modules.csv'):
        self.dataframe.to_csv(filename)
        print(f"Module assignments saved to '{filename}'")

    def run(self):
        current_time = datetime.now().strftime('%Y%m%d_%H%M%S')
        print(f'Starting analysis at {current_time}')
        self.preprocess_data()
        print(f"Data preprocessed\n\n{self.dataframe}")
        self.compute_correlation_matrix()
        print(f"Correlation matrix computed")
        self.compute_adjacency_matrix()
        print("Adjacency matrix computed")
        self.compute_TOM()
        print("Finished computing topological overlap matrix, TOM similarity for how similar two genes are in terms of network connectivity, TOM dissimilarity 1 - TOM for clustering")
        self.perform_clustering()
        print("Clustering and identifying modules")
        self.identify_modules()
        self.save_results(filename=f'{current_time}_gene_modules.csv')

    def plot_dendrogram(self):
        plt.figure(figsize=(10, 7))
        dendrogram(self.linkage_matrix)
        plt.title('Gene Clustering Dendrogram')
        plt.xlabel('Gene index')
        plt.ylabel('Distance')
        plt.show()


# example gene expression data for testing
synthetic_adata = generate_synthetic_anndata(num_cells=50, num_genes=1000, sparsity = 0.8, observations={'cell_type':['neuron'] * 50, 'disease':['normal']* 50, 'organism': ['human'] * 50})

# initialize and run WGCNA
# should update with kwarg dict to pass in when initializing Preprocessor
wgcna = WGCNA(synthetic_adata, power=6, num_modules=5)
wgcna.run()

# separate from run, plot dendrogram
# wgcna.plot_dendrogram()