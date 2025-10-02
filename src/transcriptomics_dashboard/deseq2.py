"""DESeq2 wrapper using rpy2 for differential expression analysis."""

import logging
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    RPY2_AVAILABLE = True
except ImportError:
    RPY2_AVAILABLE = False
    logging.warning("rpy2 not available. DESeq2 analysis will not work.")


logger = logging.getLogger(__name__)


class DESeq2Error(Exception):
    """Exception for DESeq2-related errors."""
    pass


class DESeq2Wrapper:
    """Wrapper for DESeq2 differential expression analysis."""
    
    def __init__(self):
        """Initialize DESeq2 wrapper and check R environment."""
        if not RPY2_AVAILABLE:
            raise DESeq2Error("rpy2 is not installed. Please install it with: pip install rpy2")
        
        self._check_r_packages()
        self._load_r_packages()
    
    def _check_r_packages(self):
        """Check if required R packages are installed."""
        required_packages = ['DESeq2', 'BiocParallel']
        
        utils = importr('utils')
        base = importr('base')
        
        installed = base.rownames(utils.installed_packages())
        
        missing = []
        for pkg in required_packages:
            if pkg not in installed:
                missing.append(pkg)
        
        if missing:
            error_msg = (
                f"Required R packages not found: {', '.join(missing)}\n"
                "Please install them in R using:\n"
                "  if (!require('BiocManager', quietly = TRUE))\n"
                "      install.packages('BiocManager')\n"
                f"  BiocManager::install(c({', '.join(f\"'{p}'\" for p in missing)}))"
            )
            raise DESeq2Error(error_msg)
    
    def _load_r_packages(self):
        """Load required R packages."""
        try:
            self.deseq2 = importr('DESeq2')
            self.biocparallel = importr('BiocParallel')
            self.base = importr('base')
            logger.info("Successfully loaded DESeq2 and dependencies")
        except Exception as e:
            raise DESeq2Error(f"Failed to load R packages: {str(e)}")
    
    def _convert_to_r_matrix(self, df: pd.DataFrame) -> ro.Matrix:
        """Convert pandas DataFrame to R matrix."""
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_df = ro.conversion.py2rpy(df)
        
        # Convert to matrix
        r_matrix = self.base.as_matrix(r_df)
        
        # Set row and column names
        r_matrix.rownames = ro.StrVector(df.index)
        r_matrix.colnames = ro.StrVector(df.columns)
        
        return r_matrix
    
    def _convert_to_r_dataframe(self, df: pd.DataFrame) -> ro.DataFrame:
        """Convert pandas DataFrame to R DataFrame."""
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_df = ro.conversion.py2rpy(df)
        return r_df
    
    def _convert_from_r_dataframe(self, r_df) -> pd.DataFrame:
        """Convert R DataFrame to pandas DataFrame."""
        with localconverter(ro.default_converter + pandas2ri.converter):
            pd_df = ro.conversion.rpy2py(r_df)
        return pd_df
    
    def create_deseq_dataset(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        design: str = "~condition"
    ):
        """
        Create DESeqDataSet object.
        
        Args:
            counts: Count matrix (genes x samples)
            metadata: Sample metadata
            design: Design formula (e.g., "~condition")
            
        Returns:
            DESeqDataSet R object
        """
        logger.info(f"Creating DESeqDataSet with design: {design}")
        
        # Ensure sample order matches
        metadata = metadata.loc[counts.columns]
        
        # Convert to R objects
        count_matrix = self._convert_to_r_matrix(counts)
        col_data = self._convert_to_r_dataframe(metadata)
        
        # Create DESeqDataSet
        try:
            dds = self.deseq2.DESeqDataSetFromMatrix(
                countData=count_matrix,
                colData=col_data,
                design=ro.Formula(design)
            )
            logger.info(f"Created DESeqDataSet with {counts.shape[0]} genes and {counts.shape[1]} samples")
            return dds
        except Exception as e:
            raise DESeq2Error(f"Failed to create DESeqDataSet: {str(e)}")
    
    def run_deseq(self, dds, parallel: bool = True, n_threads: int = 1):
        """
        Run DESeq2 analysis.
        
        Args:
            dds: DESeqDataSet object
            parallel: Use parallel processing
            n_threads: Number of threads for parallel processing
            
        Returns:
            DESeqDataSet with results
        """
        logger.info("Running DESeq2 analysis...")
        
        try:
            if parallel and n_threads > 1:
                # Set up parallel processing
                self.biocparallel.register(
                    self.biocparallel.MulticoreParam(workers=n_threads)
                )
            
            # Run DESeq
            dds = self.deseq2.DESeq(dds)
            logger.info("DESeq2 analysis completed successfully")
            return dds
            
        except Exception as e:
            error_msg = str(e)
            if "rank deficient" in error_msg.lower():
                raise DESeq2Error(
                    "Design matrix is rank deficient. This usually means:\n"
                    "  - Too few replicates per condition\n"
                    "  - Perfect correlation between variables\n"
                    "  - All samples have identical values for a variable"
                )
            else:
                raise DESeq2Error(f"DESeq2 analysis failed: {error_msg}")
    
    def get_results(
        self,
        dds,
        contrast: Optional[List[str]] = None,
        alpha: float = 0.05,
        lfc_threshold: float = 0
    ) -> pd.DataFrame:
        """
        Extract results from DESeq analysis.
        
        Args:
            dds: DESeqDataSet with results
            contrast: Contrast specification [factor, numerator, denominator]
            alpha: FDR threshold for independent filtering
            lfc_threshold: Log2 fold change threshold
            
        Returns:
            DataFrame with DE results
        """
        logger.info(f"Extracting results (alpha={alpha}, lfcThreshold={lfc_threshold})")
        
        try:
            # Get results
            if contrast:
                res = self.deseq2.results(
                    dds,
                    contrast=ro.StrVector(contrast),
                    alpha=alpha,
                    lfcThreshold=lfc_threshold
                )
            else:
                res = self.deseq2.results(
                    dds,
                    alpha=alpha,
                    lfcThreshold=lfc_threshold
                )
            
            # Convert to DataFrame
            res_df = self._convert_from_r_dataframe(self.base.as_data_frame(res))
            
            # Add gene names from row names
            with localconverter(ro.default_converter + pandas2ri.converter):
                gene_names = list(self.base.rownames(res))
            res_df.insert(0, 'gene', gene_names)
            
            # Standardize column names
            res_df.columns = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
            
            # Sort by adjusted p-value
            res_df = res_df.sort_values('padj', na_position='last')
            
            # Add significance flags
            res_df['significant'] = (
                (res_df['padj'] < alpha) & 
                (res_df['log2FoldChange'].abs() > lfc_threshold)
            )
            
            res_df['direction'] = 'not_sig'
            res_df.loc[
                (res_df['padj'] < alpha) & (res_df['log2FoldChange'] > lfc_threshold),
                'direction'
            ] = 'up'
            res_df.loc[
                (res_df['padj'] < alpha) & (res_df['log2FoldChange'] < -lfc_threshold),
                'direction'
            ] = 'down'
            
            # Count significant genes
            n_up = (res_df['direction'] == 'up').sum()
            n_down = (res_df['direction'] == 'down').sum()
            logger.info(f"Found {n_up} up-regulated and {n_down} down-regulated genes")
            
            return res_df
            
        except Exception as e:
            raise DESeq2Error(f"Failed to extract results: {str(e)}")
    
    def get_normalized_counts(self, dds) -> pd.DataFrame:
        """
        Get normalized counts from DESeqDataSet.
        
        Args:
            dds: DESeqDataSet object
            
        Returns:
            DataFrame with normalized counts
        """
        try:
            norm_counts = self.deseq2.counts(dds, normalized=True)
            norm_df = self._convert_from_r_dataframe(self.base.as_data_frame(norm_counts))
            
            # Add gene names
            with localconverter(ro.default_converter + pandas2ri.converter):
                gene_names = list(self.base.rownames(norm_counts))
            norm_df.index = gene_names
            
            return norm_df
        except Exception as e:
            raise DESeq2Error(f"Failed to get normalized counts: {str(e)}")
    
    def get_vst_counts(self, dds) -> pd.DataFrame:
        """
        Get variance-stabilized transformation of counts.
        
        Args:
            dds: DESeqDataSet object
            
        Returns:
            DataFrame with VST-transformed counts
        """
        try:
            vst = self.deseq2.vst(dds, blind=False)
            vst_matrix = self.deseq2.assay(vst)
            vst_df = self._convert_from_r_dataframe(self.base.as_data_frame(vst_matrix))
            
            # Add gene names
            with localconverter(ro.default_converter + pandas2ri.converter):
                gene_names = list(self.base.rownames(vst_matrix))
            vst_df.index = gene_names
            
            return vst_df
        except Exception as e:
            raise DESeq2Error(f"Failed to get VST counts: {str(e)}")


def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str,
    contrast: Optional[Tuple[str, str]] = None,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0,
    n_threads: int = 1
) -> Dict:
    """
    Run complete DESeq2 analysis pipeline.
    
    Args:
        counts: Count matrix (genes x samples)
        metadata: Sample metadata with condition column
        condition_col: Name of condition column in metadata
        contrast: Tuple of (numerator_condition, denominator_condition)
        fdr_threshold: FDR threshold
        lfc_threshold: Log2 fold change threshold
        n_threads: Number of threads for parallel processing
        
    Returns:
        Dictionary containing:
            - results: DE results DataFrame
            - normalized_counts: Normalized count matrix
            - vst_counts: VST-transformed counts
            - dds: DESeqDataSet object (for further analysis)
    """
    wrapper = DESeq2Wrapper()
    
    # Create design formula
    design = f"~{condition_col}"
    
    # Create DESeqDataSet
    dds = wrapper.create_deseq_dataset(counts, metadata, design)
    
    # Run analysis
    dds = wrapper.run_deseq(dds, parallel=(n_threads > 1), n_threads=n_threads)
    
    # Build contrast specification
    contrast_spec = None
    if contrast:
        contrast_spec = [condition_col, contrast[0], contrast[1]]
    
    # Get results
    results = wrapper.get_results(
        dds,
        contrast=contrast_spec,
        alpha=fdr_threshold,
        lfc_threshold=lfc_threshold
    )
    
    # Get normalized counts
    normalized_counts = wrapper.get_normalized_counts(dds)
    
    # Get VST counts
    vst_counts = wrapper.get_vst_counts(dds)
    
    return {
        'results': results,
        'normalized_counts': normalized_counts,
        'vst_counts': vst_counts,
        'dds': dds
    }


if __name__ == "__main__":
    # Example usage
    import numpy as np
    
    # Create sample data
    np.random.seed(42)
    n_genes = 1000
    n_samples = 6
    
    # Simulate counts with some differential expression
    counts = pd.DataFrame(
        np.random.negative_binomial(10, 0.1, (n_genes, n_samples)),
        index=[f"Gene_{i}" for i in range(n_genes)],
        columns=[f"Sample_{i}" for i in range(n_samples)]
    )
    
    # Add differential expression to first 50 genes
    counts.iloc[:50, 3:] = counts.iloc[:50, 3:] * 3
    
    metadata = pd.DataFrame({
        'condition': ['control'] * 3 + ['treatment'] * 3
    }, index=counts.columns)
    
    # Run analysis
    try:
        result = run_deseq2(
            counts,
            metadata,
            condition_col='condition',
            contrast=('treatment', 'control'),
            fdr_threshold=0.05,
            lfc_threshold=1.0
        )
        
        print("\nTop 10 differentially expressed genes:")
        print(result['results'].head(10)[['gene', 'log2FoldChange', 'padj', 'direction']])
        
    except DESeq2Error as e:
        print(f"Error: {e}")