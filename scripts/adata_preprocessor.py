from typing import Dict, Optional, Union, Tuple, List
import numpy as np
from scipy.sparse import issparse
import scanpy as sc
from scanpy.get import _get_obs_rep, _set_obs_rep
from anndata import AnnData
from cxg_logger import *
from datetime import datetime

# date and time as logfile suffix 
current_time = datetime.now().strftime('%Y%m%d_%H%M%S')
log_filename = f'../logs/preprocess/preprocessor_{current_time}.log'

# Set up logger to log to both console and file
logger = setup_logger('cxg_logger', log_filename)

with StreamToLogger(logger):
    print("This will be logged to both the console and the output log file with the current date and time suffix.")

class Preprocessor:
    """
    This class contains methods for preprocessing the CellXGene data prior to downstream analysis. 
    Mainly, includes methods for outlier removal (filter_gene_by_counts and filter_cell_by_counts), normalization of raw expression values (normalize_total), log1p transformation (log1p), 
    highly variable gene subsetting (subset_hvg) and correction for batch and scale differences (binning).
    """

    def __init__(
        self,
        use_key: Optional[str] = None,
        filter_gene_by_counts: Union[int, bool] = False,
        gene_filter: Union[Dict[str, Union[str, List[str]]], bool] = False,
        genes_to_keep: Optional[List[str]] = None,
        filter_cell_by_counts: Union[int, bool] = False,
        normalize_total: Union[float, bool] = 1e4,
        result_normed_key: Optional[str] = "X_normed",
        log1p: bool = False,
        result_log1p_key: str = "X_log1p",
        subset_hvg: Union[int, bool] = False,
        hvg_use_key: Optional[str] = None,
        hvg_flavor: str = "seurat_v3",
        binning: Optional[int] = None,
        result_binned_key: str = "X_binned",
        even_binning: bool = False,
        filter_observations: dict = None,
        execute_filter_genes: str = 'on',
        execute_filter_cells: str = 'on',
        execute_normalize_total: str = 'on',
        execute_log1p: str = 'on',
        execute_subset_hvg: str = 'on',
        execute_binning: str = 'on'
    ):
        """
        Initializes the Preprocessor with the specified settings.
        Args:
            (Refer to previous docstring for argument descriptions)
        """
        self.use_key = use_key
        self.filter_gene_by_counts = filter_gene_by_counts
        self.gene_filter = gene_filter
        self.genes_to_keep = genes_to_keep
        self.filter_cell_by_counts = filter_cell_by_counts
        self.normalize_total = normalize_total
        self.result_normed_key = result_normed_key
        self.log1p = log1p
        self.result_log1p_key = result_log1p_key
        self.subset_hvg = subset_hvg
        self.hvg_use_key = hvg_use_key
        self.hvg_flavor = hvg_flavor
        self.binning = binning
        self.result_binned_key = result_binned_key
        self.even_binning = even_binning
        self.filter_observations=filter_observations

        self.execute_filter_genes = execute_filter_genes
        self.execute_filter_cells = execute_filter_cells
        self.execute_normalize_total = execute_normalize_total
        self.execute_log1p = execute_log1p
        self.execute_subset_hvg = execute_subset_hvg
        self.execute_binning = execute_binning

    def __call__(self, adata: AnnData, batch_key: Optional[str] = None) -> Dict:
        """
        Processes the AnnData object based on initialized configurations.
        Args:
            adata: AnnData object to preprocess.
            batch_key: Optional batch information key for HVG selection.
        """
        key_to_process = self._resolve_key_to_process(adata)
        is_logged = self.check_logged(adata, obs_key=key_to_process)
        
        # if dictionary provided, filter out rows by specified observation column row pairs
        if self.filter_observations:
            self.filter_obs(adata, self.filter_observations)

        # new: determine what method to run via on or off parameter
        if self.execute_filter_genes == 'on':
            self._filter_genes(adata)
        if self.execute_filter_cells == 'on':
            self._filter_cells(adata)
        if self.execute_normalize_total == 'on':
            key_to_process = self._normalize_total(adata, key_to_process)
        if self.execute_log1p == 'on':
            key_to_process = self._apply_log1p(adata, key_to_process, is_logged)
        if self.execute_subset_hvg == 'on':
            self._subset_hvg(adata, batch_key)
        if self.execute_binning == 'on':
            self._apply_binning(adata, key_to_process)

    def _resolve_key_to_process(self, adata: AnnData) -> Optional[str]:
        return None if self.use_key == "X" else self.use_key

    def _filter_genes(self, adata: AnnData) -> None:
        
        if self.gene_filter:
            logger.info("Filtering genes by provided dictionary, where key is a column (check adata.var.columns for details), and value is a string or list to build query...")

            query = np.ones(adata.shape[1], dtype=bool)
            for key, value in self.gene_filter.items():
                if key not in adata.var:
                    raise ValueError(f"Column '{key}' not found in adata.var.")
                if isinstance(value, list):
                    condition = adata.var[key].isin(value)
                else:
                    condition = adata.var[key] == value
                query = query & condition
            adata._inplace_subset_var(query)
        elif self.genes_to_keep:
            logger.info("Filtering genes to keep only the specified genes...")
            genes_to_keep_set = set(self.genes_to_keep)
            genes_mask = [gene in genes_to_keep_set for gene in adata.var_names]
            adata._inplace_subset_var(genes_mask)

        if self.filter_gene_by_counts:
            logger.info("Filtering genes by counts...")
            sc.pp.filter_genes(adata, min_counts=self.filter_gene_by_counts if isinstance(self.filter_gene_by_counts, int) else None)


    def _filter_cells(self, adata: AnnData) -> None:
        if isinstance(self.filter_cell_by_counts, int):
            logger.info("Filtering cells by counts...")
            sc.pp.filter_cells(adata, min_counts=self.filter_cell_by_counts)

    def _normalize_total(self, adata: AnnData, key_to_process: Optional[str]) -> Optional[str]:
        if self.normalize_total:
            logger.info("Normalizing total counts...")
            normed_data = sc.pp.normalize_total(adata, target_sum=self.normalize_total if isinstance(self.normalize_total, float) else None, layer=key_to_process, inplace=False)["X"]
            key_to_process = self.result_normed_key or key_to_process
            _set_obs_rep(adata, normed_data, layer=key_to_process)
        return key_to_process

    def _apply_log1p(self, adata: AnnData, key_to_process: Optional[str], is_logged: bool) -> Optional[str]:
        if self.log1p:
            logger.info("Log1p transforming...")
            if is_logged:
                logger.warning(
                    "The input data appears to be already log1p transformed. "
                    "Set `log1p=False` to avoid double log1p transformation."
                )
            if self.result_log1p_key:
                _set_obs_rep(
                    adata,
                    _get_obs_rep(adata, layer=key_to_process),
                    layer=self.result_log1p_key,
                )
                key_to_process = self.result_log1p_key
            sc.pp.log1p(adata, layer=key_to_process)
        return key_to_process

    def _subset_hvg(self, adata: AnnData, batch_key: Optional[str] = None) -> None:
        if self.subset_hvg:
            logger.info("Subsetting highly variable genes...")
            if batch_key is None:
                logger.warning(
                    "No batch_key is provided. Using all cells for HVG selection."
                )
            sc.pp.highly_variable_genes(
                adata,
                layer=self.hvg_use_key,
                n_top_genes=self.subset_hvg if isinstance(self.subset_hvg, int) else None,
                batch_key=batch_key,
                flavor=self.hvg_flavor,
                subset=True,
            )

    def _apply_binning(self, adata: AnnData, key_to_process: Optional[str]) -> None:
        if self.binning:
            logger.info("Binning data...")
            if not isinstance(self.binning, int):
                raise ValueError(f"Binning argument must be an integer, but got {self.binning}.")
            binned_rows, bin_edges = self._digitize_all(adata, key_to_process)
            adata.layers[self.result_binned_key] = np.stack(binned_rows)
            adata.obsm["bin_edges"] = np.stack(bin_edges)

    def _digitize_all(self, adata: AnnData, layer: Optional[str]) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """
        Digitize all rows in the provided AnnData layer.

        Args:
            adata: AnnData object containing the data to bin.
            layer: Layer name to access the data from AnnData.

        Returns:
            A tuple of lists containing binned rows and bin edges.
        """
        n_bins = self.binning
        layer_data = _get_obs_rep(adata, layer=layer)
        layer_data = layer_data.A if issparse(layer_data) else layer_data

        binned_rows = []
        bin_edges = []
        for row in layer_data:
            non_zero_ids = row.nonzero()
            non_zero_row = row[non_zero_ids]
            bins = np.quantile(non_zero_row, np.linspace(0, 1, n_bins - 1))
            non_zero_digits = self._digitize(non_zero_row, bins)
            binned_row = np.zeros_like(row, dtype=np.int64)
            binned_row[non_zero_ids] = non_zero_digits
            binned_rows.append(binned_row)
            bin_edges.append(np.concatenate([[0], bins]))

        return binned_rows, bin_edges

    def _digitize(self, x: np.ndarray, bins: np.ndarray) -> np.ndarray:
        """
        Digitize the data into bins, distributing uniformly when bin edges have the same values.

        Args:
            x: Data to digitize.
            bins: Bin edges for digitization.

        Returns:
            Digitized data.
        """
        assert x.ndim == 1 and bins.ndim == 1
        left_digits = np.digitize(x, bins)
        right_digits = np.digitize(x, bins, right=True)
        rands = np.random.rand(len(x))
        digits = rands * (right_digits - left_digits) + left_digits
        return np.ceil(digits).astype(np.int64)
    
    def filter_obs(self, adata: AnnData, filters: Dict[str, Union[str, List[str]]]) -> None:
        """
        Filters the rows of an AnnData object based on key-value pairs in adata.obs and updates the adata object in place.

        Args:
            adata: AnnData object to filter.
            filters: Dictionary where keys are column names in adata.obs and values are the desired values to filter by. 
                    Values can be a single string or a list of strings.

        Returns:
            None
        """
        query = np.ones(adata.shape[0], dtype=bool)

        for key, value in filters.items():
            if key not in adata.obs:
                raise ValueError(f"Column '{key}' not found in adata.obs.")

            if isinstance(value, list):
                condition = adata.obs[key].isin(value)
            else:
                condition = adata.obs[key] == value

            query = query & condition

        adata._inplace_subset_obs(query)

    def check_logged(self, adata: AnnData, obs_key: Optional[str] = None) -> bool:
        """
        Checks if the data is already log1p transformed.

        Args:
            adata: AnnData object to check.
            obs_key: Key to access the data in the AnnData object.

        Returns:
            True if data is already log1p transformed, False otherwise.
        """
        data = _get_obs_rep(adata, layer=obs_key)
        max_, min_ = data.max(), data.min()
        if max_ > 30 or min_ < 0:
            return False

        non_zero_min = data[data > 0].min()
        return non_zero_min < 1
