"""
Python-based Differential Expression Analyzer
Pure Python implementation for RNA-seq differential expression analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from loguru import logger
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import quantile_transform
from statsmodels.stats.multitest import multipletests
import warnings

warnings.filterwarnings('ignore')


class PythonDEAnalyzer:
    """
    Pure Python differential expression analyzer.
    Implements methods similar to edgeR/limma for RNA-seq analysis.
    """
    
    def __init__(self, output_dir: str = "results/differential_expression"):
        """
        Initialize Python DE analyzer.
        
        Args:
            output_dir: Directory to store DE analysis results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def analyze_differential_expression(self,
                                      count_matrix: pd.DataFrame,
                                      sample_metadata: pd.DataFrame,
                                      control_condition: str = "control",
                                      treatment_condition: str = "treatment",
                                      method: str = "edger_like",
                                      padj_threshold: float = 0.05,
                                      lfc_threshold: float = 1.0) -> Dict:
        """
        Perform differential expression analysis using Python methods.
        
        Args:
            count_matrix: Gene count matrix (genes as rows, samples as columns)
            sample_metadata: Sample metadata with condition information
            control_condition: Name of control condition
            treatment_condition: Name of treatment condition
            method: Analysis method ('edger_like', 'limma_like', 'simple_ttest')
            padj_threshold: Adjusted p-value threshold
            lfc_threshold: Log2 fold change threshold
            
        Returns:
            Dictionary containing DE results and statistics
        """
        logger.info(f"Starting Python DE analysis using {method} method")
        
        # Validate inputs
        self._validate_inputs(count_matrix, sample_metadata)
        
        # Choose analysis method
        if method == "edger_like":
            results = self._run_edger_like_analysis(
                count_matrix, sample_metadata, control_condition, treatment_condition
            )
        elif method == "limma_like":
            results = self._run_limma_like_analysis(
                count_matrix, sample_metadata, control_condition, treatment_condition
            )
        elif method == "simple_ttest":
            results = self._run_simple_ttest(
                count_matrix, sample_metadata, control_condition, treatment_condition
            )
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Filter significant genes
        significant_genes = self._filter_significant_genes(
            results, padj_threshold, lfc_threshold
        )
        
        # Create summary statistics
        summary_stats = self._create_summary_stats(results, significant_genes)
        
        # Save results
        self._save_results(results, significant_genes, summary_stats, method)
        
        return {
            'all_results': results,
            'significant_genes': significant_genes,
            'summary_stats': summary_stats,
            'method_used': f'Python_{method}'
        }
    
    def _validate_inputs(self, count_matrix: pd.DataFrame, sample_metadata: pd.DataFrame):
        """Validate input data."""
        if count_matrix.empty:
            raise ValueError("Count matrix is empty")
        
        if sample_metadata.empty:
            raise ValueError("Sample metadata is empty")
        
        if 'condition' not in sample_metadata.columns:
            raise ValueError("Sample metadata must contain 'condition' column")
    
    def _run_edger_like_analysis(self, count_matrix: pd.DataFrame, 
                                sample_metadata: pd.DataFrame,
                                control_condition: str, 
                                treatment_condition: str) -> pd.DataFrame:
        """
        Run edgeR-like analysis using negative binomial distribution assumptions.
        """
        logger.info("Running edgeR-like analysis")
        
        # Filter samples
        control_samples = sample_metadata[sample_metadata['condition'] == control_condition].index
        treatment_samples = sample_metadata[sample_metadata['condition'] == treatment_condition].index
        
        all_samples = list(control_samples) + list(treatment_samples)
        count_matrix_filtered = count_matrix[all_samples]
        
        # Filter low count genes
        min_count = 10
        min_samples = max(2, len(all_samples) // 4)
        keep_genes = (count_matrix_filtered >= min_count).sum(axis=1) >= min_samples
        count_matrix_filtered = count_matrix_filtered[keep_genes]
        
        logger.info(f"Analyzing {len(count_matrix_filtered)} genes across {len(all_samples)} samples")
        
        # Calculate normalization factors (TMM-like)
        norm_factors = self._calculate_tmm_factors(count_matrix_filtered)
        
        results_list = []
        
        for gene in count_matrix_filtered.index:
            try:
                control_counts = count_matrix_filtered.loc[gene, control_samples].values
                treatment_counts = count_matrix_filtered.loc[gene, treatment_samples].values
                
                # Apply normalization
                control_norm = control_counts / norm_factors[control_samples]
                treatment_norm = treatment_counts / norm_factors[treatment_samples]
                
                # Calculate means
                control_mean = np.mean(control_norm)
                treatment_mean = np.mean(treatment_norm)
                base_mean = np.mean([control_mean, treatment_mean])
                
                # Calculate log2 fold change
                if control_mean == 0:
                    control_mean = 0.1
                if treatment_mean == 0:
                    treatment_mean = 0.1
                
                log2_fold_change = np.log2(treatment_mean / control_mean)
                
                # Estimate dispersion and perform test
                dispersion = self._estimate_dispersion(control_norm, treatment_norm)
                
                # Negative binomial test (approximated using normal distribution)
                control_var = control_mean + dispersion * (control_mean ** 2)
                treatment_var = treatment_mean + dispersion * (treatment_mean ** 2)
                
                pooled_se = np.sqrt(control_var / len(control_norm) + treatment_var / len(treatment_norm))
                
                if pooled_se > 0:
                    z_stat = log2_fold_change / (pooled_se / np.log(2))
                    p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))
                else:
                    z_stat = 0
                    p_value = 1.0
                
                results_list.append({
                    'gene': gene,
                    'baseMean': base_mean,
                    'log2FoldChange': log2_fold_change,
                    'stat': z_stat,
                    'pvalue': p_value
                })
                
            except Exception as e:
                logger.warning(f"Error analyzing gene {gene}: {e}")
                continue
        
        # Create results dataframe
        results_df = pd.DataFrame(results_list)
        results_df.set_index('gene', inplace=True)
        
        # Adjust p-values
        results_df = self._adjust_pvalues(results_df)
        
        return results_df
    
    def _run_limma_like_analysis(self, count_matrix: pd.DataFrame,
                                sample_metadata: pd.DataFrame,
                                control_condition: str,
                                treatment_condition: str) -> pd.DataFrame:
        """
        Run limma-like analysis on log-transformed data.
        """
        logger.info("Running limma-like analysis")
        
        # Filter samples
        control_samples = sample_metadata[sample_metadata['condition'] == control_condition].index
        treatment_samples = sample_metadata[sample_metadata['condition'] == treatment_condition].index
        
        all_samples = list(control_samples) + list(treatment_samples)
        count_matrix_filtered = count_matrix[all_samples]
        
        # Filter low count genes
        min_count = 10
        min_samples = max(2, len(all_samples) // 4)
        keep_genes = (count_matrix_filtered >= min_count).sum(axis=1) >= min_samples
        count_matrix_filtered = count_matrix_filtered[keep_genes]
        
        # Normalize and log-transform
        norm_factors = self._calculate_tmm_factors(count_matrix_filtered)
        normalized_counts = count_matrix_filtered.astype(float).div(norm_factors, axis=1)
        
        # Log2 transform with pseudocount
        log_counts = np.log2(normalized_counts + 1)
        
        results_list = []
        
        for gene in log_counts.index:
            try:
                control_values = log_counts.loc[gene, control_samples].values
                treatment_values = log_counts.loc[gene, treatment_samples].values
                
                # Calculate statistics
                control_mean = np.mean(control_values)
                treatment_mean = np.mean(treatment_values)
                log2_fold_change = treatment_mean - control_mean
                base_mean = np.mean([control_mean, treatment_mean])
                
                # Moderated t-test (simplified version)
                pooled_std = np.sqrt(
                    ((len(control_values) - 1) * np.var(control_values, ddof=1) +
                     (len(treatment_values) - 1) * np.var(treatment_values, ddof=1)) /
                    (len(control_values) + len(treatment_values) - 2)
                )
                
                if pooled_std > 0:
                    t_stat = log2_fold_change / (pooled_std * np.sqrt(1/len(control_values) + 1/len(treatment_values)))
                    df = len(control_values) + len(treatment_values) - 2
                    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df))
                else:
                    t_stat = 0
                    p_value = 1.0
                
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
        
        # Adjust p-values
        results_df = self._adjust_pvalues(results_df)
        
        return results_df
    
    def _run_simple_ttest(self, count_matrix: pd.DataFrame,
                         sample_metadata: pd.DataFrame,
                         control_condition: str,
                         treatment_condition: str) -> pd.DataFrame:
        """
        Run simple t-test on log-transformed normalized counts.
        """
        logger.info("Running simple t-test analysis")
        
        # Filter samples
        control_samples = sample_metadata[sample_metadata['condition'] == control_condition].index
        treatment_samples = sample_metadata[sample_metadata['condition'] == treatment_condition].index
        
        all_samples = list(control_samples) + list(treatment_samples)
        count_matrix_filtered = count_matrix[all_samples]
        
        # Filter low count genes
        min_count = 5
        min_samples = max(2, len(all_samples) // 4)
        keep_genes = (count_matrix_filtered >= min_count).sum(axis=1) >= min_samples
        count_matrix_filtered = count_matrix_filtered[keep_genes]
        
        # Simple normalization (library size)
        lib_sizes = count_matrix_filtered.sum(axis=0).astype(float)
        normalized_counts = count_matrix_filtered.astype(float).div(lib_sizes, axis=1) * 1e6  # CPM
        
        # Log2 transform
        log_counts = np.log2(normalized_counts + 1)
        
        results_list = []
        
        for gene in log_counts.index:
            try:
                control_values = log_counts.loc[gene, control_samples].values
                treatment_values = log_counts.loc[gene, treatment_samples].values
                
                # Calculate basic statistics
                control_mean = np.mean(control_values)
                treatment_mean = np.mean(treatment_values)
                log2_fold_change = treatment_mean - control_mean
                base_mean = np.mean([control_mean, treatment_mean])
                
                # Welch's t-test
                t_stat, p_value = stats.ttest_ind(treatment_values, control_values, equal_var=False)
                
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
        
        # Adjust p-values
        results_df = self._adjust_pvalues(results_df)
        
        return results_df
    
    def _calculate_tmm_factors(self, count_matrix: pd.DataFrame) -> pd.Series:
        """
        Calculate TMM (Trimmed Mean of M-values) normalization factors.
        Simplified implementation.
        """
        # Calculate library sizes
        lib_sizes = count_matrix.sum(axis=0).astype(float)
        
        # Use median library size as reference
        ref_lib_size = float(np.median(lib_sizes))
        
        # Calculate scaling factors
        scaling_factors = lib_sizes / ref_lib_size
        
        return scaling_factors.astype(float)
    
    def _estimate_dispersion(self, control_counts: np.ndarray, 
                           treatment_counts: np.ndarray) -> float:
        """
        Estimate dispersion parameter for negative binomial distribution.
        """
        all_counts = np.concatenate([control_counts, treatment_counts])
        
        if len(all_counts) < 3:
            return 0.1  # Default dispersion
        
        mean_count = np.mean(all_counts)
        var_count = np.var(all_counts)
        
        if mean_count > 0 and var_count > mean_count:
            # Estimate dispersion from mean-variance relationship
            dispersion = (var_count - mean_count) / (mean_count ** 2)
            return max(0.01, min(dispersion, 10))  # Bound dispersion
        else:
            return 0.1
    
    def _adjust_pvalues(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """Adjust p-values for multiple testing."""
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
            'percent_significant': (len(significant_genes) / len(all_results)) * 100 if len(all_results) > 0 else 0
        }
    
    def _save_results(self, all_results: pd.DataFrame, 
                     significant_genes: pd.DataFrame, 
                     summary_stats: Dict,
                     method: str):
        """Save analysis results to files."""
        # Save all results
        all_results.to_csv(self.output_dir / f"python_{method}_all_results.csv")
        
        # Save significant genes
        significant_genes.to_csv(self.output_dir / f"python_{method}_significant_genes.csv")
        
        # Save summary statistics
        summary_df = pd.DataFrame([summary_stats])
        summary_df.to_csv(self.output_dir / f"python_{method}_summary.csv", index=False)
        
        logger.info(f"Python {method} results saved to {self.output_dir}")
    
    def compare_methods(self, count_matrix: pd.DataFrame,
                       sample_metadata: pd.DataFrame,
                       control_condition: str = "control",
                       treatment_condition: str = "treatment") -> Dict:
        """
        Compare different DE analysis methods.
        
        Returns:
            Dictionary containing results from all methods
        """
        logger.info("Comparing DE analysis methods")
        
        methods = ['edger_like', 'limma_like', 'simple_ttest']
        all_method_results = {}
        
        for method in methods:
            try:
                logger.info(f"Running {method} analysis")
                result = self.analyze_differential_expression(
                    count_matrix, sample_metadata, 
                    control_condition, treatment_condition,
                    method=method
                )
                all_method_results[method] = result
            except Exception as e:
                logger.error(f"Error running {method}: {e}")
                continue
        
        # Create comparison summary
        comparison_summary = self._create_method_comparison(all_method_results)
        
        return {
            'method_results': all_method_results,
            'comparison_summary': comparison_summary
        }
    
    def _create_method_comparison(self, all_method_results: Dict) -> pd.DataFrame:
        """Create a comparison table of different methods."""
        comparison_data = []
        
        for method, results in all_method_results.items():
            stats = results['summary_stats']
            comparison_data.append({
                'Method': method,
                'Total_Genes_Tested': stats['total_genes_tested'],
                'Significant_Genes': stats['total_significant'],
                'Upregulated': stats['upregulated'],
                'Downregulated': stats['downregulated'],
                'Percent_Significant': stats['percent_significant']
            })
        
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df.to_csv(self.output_dir / "method_comparison.csv", index=False)
        
        return comparison_df 