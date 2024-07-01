import anndata
import numpy as np
import scipy.sparse as sp 

def generate_synthetic_anndata(num_cells:int, num_genes:int, sparsity:float=0.1, observations:dict=None) -> anndata.AnnData:
    
    # create random data
    data = np.random.poisson(1, (num_cells, num_genes))
    # make sparse
    mask = np.random.binomial(1, 1 - sparsity, data.shape).astype(bool)
    data[~mask] = 0
    # make sparse matrix data type
    data_sparse = sp.csr_matrix(data)

    # create anndata object
    adata = anndata.AnnData(X=data_sparse)

    adata.obs['cell_id'] = [f'cell_{i}' for i in range(num_cells)]
    adata.var['gene_id'] = [f'gene_{i}' for i in range(num_genes)]

    # add metadata? optional values are lists that should have correct matching dimensions with num_cells
    if observations:
        for ok, ov in observations.items():
            assert num_cells == len(ov)
            adata.obs[ok] = ov

    return adata

