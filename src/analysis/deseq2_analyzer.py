"""
DESeq2 Differential Expression Analyzer
Provides interface to R/DESeq2 and Python-based alternatives for DE analysis.
"""

import os
import subprocess
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
from loguru import logger
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import nbinom
from sklearn.preprocessing import StandardScaler
import warnings

warnings.filterwarnings('ignore')


class DESeq2Analyzer:
    """
    Differential expression analyzer using DESeq2 and Python alternatives.
    """
    
    def __init__(self, output_dir: str = "results/differential_expression"):
        """
        Initialize DESeq2 analyzer.
        
        Args:
            output_dir: Directory to store DE analysis results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.r_available = self._check_r_available()
        self.deseq2_available = self._check_deseq2_available() if self.r_available else False
        
    def _check_r_available(self) -> bool:
        """Check if R is available on the system."""
        try:
            result = subprocess.run(['R', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def _check_deseq2_available(self) -> bool:
        """Check if DESeq2 is installed in R."""
        try:
            r_script = '''
            if (!require("DESeq2", quietly = TRUE)) {
                quit(status = 1)
            } else {
                quit(status = 0)
            }
            '''
            result = subprocess.run(['R', '--slave', '--no-restore'], 
                                  input=r_script, text=True, 
                                  capture_output=True, timeout=30)
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def analyze_differential_expression(self,
                                      count_matrix: pd.DataFrame,
                                      sample_metadata: pd.DataFrame,
                                      control_condition: str = "control",
                                      treatment_condition: str = "treatment",
                                      padj_threshold: float = 0.05,
                                      lfc_threshold: float = 1.0,
                                      use_deseq2: bool = True) -> Dict:
        """
        Perform differential expression analysis.
        
        Args:
            count_matrix: Gene count matrix (genes as rows, samples as columns)
            sample_metadata: Sample metadata with condition information
            control_condition: Name of control condition
            treatment_condition: Name of treatment condition
            padj_threshold: Adjusted p-value threshold
            lfc_threshold: Log2 fold change threshold
            use_deseq2: Whether to use DESeq2 (if available) or Python alternative
            
        Returns:
            Dictionary containing DE results and statistics
        """
        logger.info("Starting differential expression analysis")
        
        # Validate inputs
        self._validate_inputs(count_matrix, sample_metadata)
        
        # Try DESeq2 first if requested and available
        if use_deseq2 and self.deseq2_available:
            logger.info("Using DESeq2 for differential expression analysis")
            try:
                results = self._run_deseq2(count_matrix, sample_metadata, 
                                         control_condition, treatment_condition)
            except Exception as e:
                logger.warning(f"DESeq2 failed: {e}. Falling back to Python implementation.")
                results = self._run_python_de(count_matrix, sample_metadata,
                                            control_condition, treatment_condition)
        else:
            logger.info("Using Python implementation for differential expression analysis")
            results = self._run_python_de(count_matrix, sample_metadata,
                                        control_condition, treatment_condition)
        
        # Filter significant genes
        significant_genes = self._filter_significant_genes(
            results, padj_threshold, lfc_threshold
        )
        
        # Create summary statistics
        summary_stats = self._create_summary_stats(results, significant_genes)
        
        # Save results
        self._save_results(results, significant_genes, summary_stats)
        
        return {
            'all_results': results,
            'significant_genes': significant_genes,
            'summary_stats': summary_stats,
            'method_used': 'DESeq2' if (use_deseq2 and self.deseq2_available) else 'Python'
        }
    
    def _validate_inputs(self, count_matrix: pd.DataFrame, sample_metadata: pd.DataFrame):
        """Validate input data."""
        if count_matrix.empty:
            raise ValueError("Count matrix is empty")
        
        if sample_metadata.empty:
            raise ValueError("Sample metadata is empty")
        
        if 'condition' not in sample_metadata.columns:
            raise ValueError("Sample metadata must contain 'condition' column")
        
        # Check if sample names match
        common_samples = set(count_matrix.columns) & set(sample_metadata.index)
        if len(common_samples) != len(count_matrix.columns):
            logger.warning("Not all samples in count matrix have metadata")
    
    def _run_deseq2(self, count_matrix: pd.DataFrame, sample_metadata: pd.DataFrame,
                   control_condition: str, treatment_condition: str) -> pd.DataFrame:
        """Run DESeq2 analysis using R."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            
            # Save count matrix and metadata
            count_file = temp_dir / "counts.csv"
            metadata_file = temp_dir / "metadata.csv"
            results_file = temp_dir / "deseq2_results.csv"
            
            count_matrix.to_csv(count_file)
            sample_metadata.to_csv(metadata_file)
            
            # Create R script
            r_script = f'''
            library(DESeq2)
            library(readr)
            
            # Load data
            counts <- read.csv("{count_file}", row.names=1)
            metadata <- read.csv("{metadata_file}", row.names=1)
            
            # Ensure sample order matches
            metadata <- metadata[colnames(counts), , drop=FALSE]
            
            # Filter low count genes
            keep <- rowSums(counts >= 10) >= 3
            counts <- counts[keep, ]
            
            # Create DESeq2 object
            dds <- DESeqDataSetFromMatrix(countData = counts,
                                        colData = metadata,
                                        design = ~ condition)
            
            # Set reference level
            dds$condition <- relevel(dds$condition, ref = "{control_condition}")
            
            # Run DESeq2
            dds <- DESeq(dds)
            
            # Get results
            res <- results(dds, contrast = c("condition", "{treatment_condition}", "{control_condition}"))
            
            # Convert to data frame and save
            res_df <- as.data.frame(res)
            res_df$gene <- rownames(res_df)
            write.csv(res_df, "{results_file}", row.names=FALSE)
            '''
            
            # Run R script
            result = subprocess.run(['R', '--slave', '--no-restore'], 
                                  input=r_script, text=True, 
                                  capture_output=True, timeout=300)
            
            if result.returncode != 0:
                raise RuntimeError(f"DESeq2 failed: {result.stderr}")
            
            # Load results
            results_df = pd.read_csv(results_file)
            results_df.set_index('gene', inplace=True)
            
            return results_df
    
    def _run_python_de(self, count_matrix: pd.DataFrame, sample_metadata: pd.DataFrame,
                      control_condition: str, treatment_condition: str) -> pd.DataFrame:
        """Run differential expression analysis using Python (edgeR-like approach)."""
        logger.info("Running Python-based differential expression analysis")
        
        # Filter samples
        control_samples = sample_metadata[sample_metadata['condition'] == control_condition].index
        treatment_samples = sample_metadata[sample_metadata['condition'] == treatment_condition].index
        
        all_samples = list(control_samples) + list(treatment_samples)
        count_matrix_filtered = count_matrix[all_samples]
        
        # Filter low count genes
        min_count = 10
        min_samples = 3
        keep_genes = (count_matrix_filtered >= min_count).sum(axis=1) >= min_samples
        count_matrix_filtered = count_matrix_filtered[keep_genes]
        
        logger.info(f"Analyzing {len(count_matrix_filtered)} genes across {len(all_samples)} samples")
        
        # Calculate library sizes and normalization factors
        lib_sizes = count_matrix_filtered.sum(axis=0)
        norm_factors = lib_sizes / np.median(lib_sizes)
        
        results_list = []
        
        for gene in count_matrix_filtered.index:
            try:
                control_counts = count_matrix_filtered.loc[gene, control_samples].values
                treatment_counts = count_matrix_filtered.loc[gene, treatment_samples].values
                
                # Normalize counts
                control_norm = control_counts / norm_factors[control_samples]
                treatment_norm = treatment_counts / norm_factors[treatment_samples]
                
                # Calculate basic statistics
                control_mean = np.mean(control_norm)
                treatment_mean = np.mean(treatment_norm)
                
                # Avoid division by zero
                if control_mean == 0:
                    control_mean = 0.1
                if treatment_mean == 0:
                    treatment_mean = 0.1
                
                log2_fold_change = np.log2(treatment_mean / control_mean)
                base_mean = np.mean([control_mean, treatment_mean])
                
                # Statistical test (Welch's t-test on log-transformed data)
                control_log = np.log2(control_norm + 1)
                treatment_log = np.log2(treatment_norm + 1)
                
                t_stat, p_value = stats.ttest_ind(treatment_log, control_log, equal_var=False)
                
                results_list.append({
                    'gene': gene,
                    'baseMean': base_mean,
                    'log2FoldChange': log2_fold_change,
                    'stat': t_stat,
                    'pvalue': p_value
                })
                
            except Exception as e:
                logger.warning(f"Error analyzing gene {gene}: {e}")
                continue
        
        # Create results dataframe
        results_df = pd.DataFrame(results_list)
        results_df.set_index('gene', inplace=True)
        
        # Adjust p-values using Benjamini-Hochberg
        from statsmodels.stats.multitest import multipletests
        
        valid_pvals = ~np.isnan(results_df['pvalue'])
        adjusted_pvals = np.full(len(results_df), np.nan)
        
        if valid_pvals.sum() > 0:
            _, adj_pvals, _, _ = multipletests(
                results_df.loc[valid_pvals, 'pvalue'], 
                method='fdr_bh'
            )
            adjusted_pvals[valid_pvals] = adj_pvals
        
        results_df['padj'] = adjusted_pvals
        
        return results_df
    
    def _filter_significant_genes(self, results: pd.DataFrame, 
                                padj_threshold: float, 
                                lfc_threshold: float) -> pd.DataFrame:
        """Filter significantly differentially expressed genes."""
        significant = (
            (results['padj'] < padj_threshold) & 
            (results['log2FoldChange'].abs() > lfc_threshold) &
            (~results['padj'].isna())
        )
        
        return results[significant].copy()
    
    def _create_summary_stats(self, all_results: pd.DataFrame, 
                            significant_genes: pd.DataFrame) -> Dict:
        """Create summary statistics."""
        upregulated = significant_genes[significant_genes['log2FoldChange'] > 0]
        downregulated = significant_genes[significant_genes['log2FoldChange'] < 0]
        
        return {
            'total_genes_tested': len(all_results),
            'total_significant': len(significant_genes),
            'upregulated': len(upregulated),
            'downregulated': len(downregulated),
            'percent_significant': (len(significant_genes) / len(all_results)) * 100
        }
    
    def _save_results(self, all_results: pd.DataFrame, 
                     significant_genes: pd.DataFrame, 
                     summary_stats: Dict):
        """Save analysis results to files."""
        # Save all results
        all_results.to_csv(self.output_dir / "all_de_results.csv")
        
        # Save significant genes
        significant_genes.to_csv(self.output_dir / "significant_genes.csv")
        
        # Save summary statistics
        summary_df = pd.DataFrame([summary_stats])
        summary_df.to_csv(self.output_dir / "summary_statistics.csv", index=False)
        
        logger.info(f"Results saved to {self.output_dir}")
    
    def create_de_plots(self, results: pd.DataFrame, 
                       padj_threshold: float = 0.05,
                       lfc_threshold: float = 1.0) -> Dict[str, str]:
        """
        Create differential expression plots.
        
        Args:
            results: DE results dataframe
            padj_threshold: Adjusted p-value threshold
            lfc_threshold: Log fold change threshold
            
        Returns:
            Dictionary of plot file paths
        """
        plot_paths = {}
        
        # Filter out genes with missing values
        clean_results = results.dropna(subset=['log2FoldChange', 'padj'])
        
        if len(clean_results) == 0:
            logger.warning("No valid results for plotting")
            return plot_paths
        
        # Volcano plot
        plt.figure(figsize=(10, 8))
        
        # Define significance
        significant = (
            (clean_results['padj'] < padj_threshold) & 
            (clean_results['log2FoldChange'].abs() > lfc_threshold)
        )
        
        # Plot non-significant genes
        plt.scatter(clean_results.loc[~significant, 'log2FoldChange'],
                   -np.log10(clean_results.loc[~significant, 'padj']),
                   alpha=0.6, color='lightgray', s=20, label='Non-significant')
        
        # Plot significant genes
        if significant.sum() > 0:
            sig_up = significant & (clean_results['log2FoldChange'] > 0)
            sig_down = significant & (clean_results['log2FoldChange'] < 0)
            
            if sig_up.sum() > 0:
                plt.scatter(clean_results.loc[sig_up, 'log2FoldChange'],
                           -np.log10(clean_results.loc[sig_up, 'padj']),
                           alpha=0.8, color='red', s=30, label=f'Upregulated ({sig_up.sum()})')
            
            if sig_down.sum() > 0:
                plt.scatter(clean_results.loc[sig_down, 'log2FoldChange'],
                           -np.log10(clean_results.loc[sig_down, 'padj']),
                           alpha=0.8, color='blue', s=30, label=f'Downregulated ({sig_down.sum()})')
        
        # Add threshold lines
        plt.axhline(y=-np.log10(padj_threshold), color='black', linestyle='--', alpha=0.7)
        plt.axvline(x=lfc_threshold, color='black', linestyle='--', alpha=0.7)
        plt.axvline(x=-lfc_threshold, color='black', linestyle='--', alpha=0.7)
        
        plt.xlabel('log₂(Fold Change)', fontsize=12)
        plt.ylabel('-log₁₀(Adjusted p-value)', fontsize=12)
        plt.title('Volcano Plot - Differential Expression', fontsize=14, fontweight='bold')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        volcano_path = self.output_dir / "volcano_plot.png"
        plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['volcano_plot'] = str(volcano_path)
        
        # MA plot
        plt.figure(figsize=(10, 8))
        
        plt.scatter(clean_results.loc[~significant, 'baseMean'],
                   clean_results.loc[~significant, 'log2FoldChange'],
                   alpha=0.6, color='lightgray', s=20, label='Non-significant')
        
        if significant.sum() > 0:
            plt.scatter(clean_results.loc[significant, 'baseMean'],
                       clean_results.loc[significant, 'log2FoldChange'],
                       alpha=0.8, color='red', s=30, label=f'Significant ({significant.sum()})')
        
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        plt.axhline(y=lfc_threshold, color='black', linestyle='--', alpha=0.7)
        plt.axhline(y=-lfc_threshold, color='black', linestyle='--', alpha=0.7)
        
        plt.xlabel('Mean Expression', fontsize=12)
        plt.ylabel('log₂(Fold Change)', fontsize=12)
        plt.title('MA Plot - Differential Expression', fontsize=14, fontweight='bold')
        plt.xscale('log')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        ma_path = self.output_dir / "ma_plot.png"
        plt.savefig(ma_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['ma_plot'] = str(ma_path)
        
        logger.info(f"Created DE plots: {list(plot_paths.keys())}")
        return plot_paths
    
    def create_mock_count_data(self, n_genes: int = 1000, n_samples: int = 6) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Create mock count data for demonstration.
        
        Args:
            n_genes: Number of genes
            n_samples: Number of samples
            
        Returns:
            Tuple of (count_matrix, sample_metadata)
        """
        np.random.seed(42)
        
        # Create sample metadata
        conditions = ['control'] * (n_samples // 2) + ['treatment'] * (n_samples // 2)
        sample_names = [f"Sample_{i+1}" for i in range(n_samples)]
        
        sample_metadata = pd.DataFrame({
            'condition': conditions
        }, index=sample_names)
        
        # Create count matrix
        gene_names = [f"Gene_{i+1:04d}" for i in range(n_genes)]
        
        # Generate realistic count data using negative binomial distribution
        count_data = np.random.negative_binomial(100, 0.1, size=(n_genes, n_samples)).astype(int)
        
        # Add some differential expression
        de_genes = n_genes // 10  # 10% DE genes
        for i in range(de_genes):
            # Increase expression in treatment samples
            treatment_cols = n_samples // 2
            multiplier = np.random.uniform(2, 5)
            count_data[i, treatment_cols:] = (count_data[i, treatment_cols:] * multiplier).astype(int)
        
        count_matrix = pd.DataFrame(count_data, 
                                   index=gene_names, 
                                   columns=sample_names)
        
        return count_matrix, sample_metadata


def install_deseq2():
    """Install DESeq2 in R (helper function)."""
    r_script = '''
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("DESeq2")
    '''
    
    try:
        result = subprocess.run(['R', '--slave', '--no-restore'], 
                              input=r_script, text=True, 
                              capture_output=True, timeout=600)
        
        if result.returncode == 0:
            print("DESeq2 installed successfully!")
            return True
        else:
            print(f"Failed to install DESeq2: {result.stderr}")
            return False
    except Exception as e:
        print(f"Error installing DESeq2: {e}")
        return False 