"""
Plot generator for RNA-seq analysis.
Creates various plots for RNA-seq data visualization.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from loguru import logger
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots


class PlotGenerator:
    """Generate various plots for RNA-seq analysis."""
    
    def __init__(self, output_dir: str = "results/plots"):
        """
        Initialize plot generator.
        
        Args:
            output_dir: Directory to store generated plots
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set plotting style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
        
    def create_quality_control_plots(self, qc_data: pd.DataFrame) -> Dict[str, str]:
        """
        Create quality control plots.
        
        Args:
            qc_data: DataFrame with QC metrics
            
        Returns:
            Dictionary mapping plot names to file paths
        """
        plot_paths = {}
        
        # 1. Per base quality plot
        plt.figure(figsize=(12, 6))
        plt.title("Per Base Sequence Quality")
        plt.xlabel("Position in read")
        plt.ylabel("Quality score")
        plt.axhline(y=20, color='red', linestyle='--', alpha=0.7, label='Quality threshold')
        plt.legend()
        plot_path = self.output_dir / "per_base_quality.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['per_base_quality'] = str(plot_path)
        
        # 2. Per sequence quality distribution
        plt.figure(figsize=(10, 6))
        plt.title("Per Sequence Quality Distribution")
        plt.xlabel("Quality Score")
        plt.ylabel("Count")
        plot_path = self.output_dir / "per_sequence_quality.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['per_sequence_quality'] = str(plot_path)
        
        # 3. GC content distribution
        plt.figure(figsize=(10, 6))
        plt.title("GC Content Distribution")
        plt.xlabel("GC Content (%)")
        plt.ylabel("Count")
        plot_path = self.output_dir / "gc_content.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['gc_content'] = str(plot_path)
        
        return plot_paths
    
    def create_alignment_plots(self, alignment_stats: Dict) -> Dict[str, str]:
        """
        Create alignment statistics plots.
        
        Args:
            alignment_stats: Dictionary with alignment statistics
            
        Returns:
            Dictionary mapping plot names to file paths
        """
        plot_paths = {}
        
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(alignment_stats, orient='index')
        
        # 1. Mapping rate bar plot
        plt.figure(figsize=(12, 6))
        plt.bar(df.index, df['mapped_percent'])
        plt.title("Mapping Rate by Sample")
        plt.xlabel("Sample")
        plt.ylabel("Mapping Rate (%)")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plot_path = self.output_dir / "mapping_rate.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['mapping_rate'] = str(plot_path)
        
        # 2. Total reads vs mapped reads
        plt.figure(figsize=(10, 6))
        plt.scatter(df['total_reads'], df['mapped_reads'], alpha=0.7)
        plt.plot([df['total_reads'].min(), df['total_reads'].max()], 
                [df['total_reads'].min(), df['total_reads'].max()], 
                'r--', alpha=0.5, label='Perfect mapping')
        plt.title("Total Reads vs Mapped Reads")
        plt.xlabel("Total Reads")
        plt.ylabel("Mapped Reads")
        plt.legend()
        plot_path = self.output_dir / "reads_comparison.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['reads_comparison'] = str(plot_path)
        
        return plot_paths
    
    def create_expression_plots(self, expression_data: pd.DataFrame) -> Dict[str, str]:
        """
        Create gene expression plots.
        
        Args:
            expression_data: DataFrame with gene expression data
            
        Returns:
            Dictionary mapping plot names to file paths
        """
        plot_paths = {}
        
        # 1. Expression distribution
        plt.figure(figsize=(12, 6))
        plt.hist(np.log2(expression_data.values.flatten() + 1), bins=50, alpha=0.7)
        plt.title("Gene Expression Distribution (log2)")
        plt.xlabel("log2(Expression + 1)")
        plt.ylabel("Frequency")
        plot_path = self.output_dir / "expression_distribution.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['expression_distribution'] = str(plot_path)
        
        # 2. Sample correlation heatmap
        plt.figure(figsize=(10, 8))
        correlation_matrix = expression_data.corr()
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0)
        plt.title("Sample Correlation Matrix")
        plot_path = self.output_dir / "sample_correlation.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['sample_correlation'] = str(plot_path)
        
        # 3. PCA plot
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        
        # Prepare data for PCA
        scaled_data = StandardScaler().fit_transform(expression_data.T)
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(scaled_data)
        
        plt.figure(figsize=(10, 8))
        plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7)
        for i, sample in enumerate(expression_data.columns):
            plt.annotate(sample, (pca_result[i, 0], pca_result[i, 1]))
        plt.title("PCA of Gene Expression Data")
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
        plot_path = self.output_dir / "pca_plot.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['pca_plot'] = str(plot_path)
        
        return plot_paths
    
    def create_differential_expression_plots(self, 
                                           de_results: pd.DataFrame,
                                           logfc_threshold: float = 1.0,
                                           pvalue_threshold: float = 0.05) -> Dict[str, str]:
        """
        Create differential expression analysis plots.
        
        Args:
            de_results: DataFrame with differential expression results
            logfc_threshold: Log fold change threshold
            pvalue_threshold: P-value threshold
            
        Returns:
            Dictionary mapping plot names to file paths
        """
        plot_paths = {}
        
        # 1. Volcano plot
        plt.figure(figsize=(10, 8))
        
        # Color points based on significance
        significant = (de_results['log2FoldChange'].abs() > logfc_threshold) & \
                     (de_results['padj'] < pvalue_threshold)
        
        plt.scatter(de_results['log2FoldChange'], 
                   -np.log10(de_results['padj']), 
                   alpha=0.6, c=significant, cmap='viridis')
        
        plt.axhline(-np.log10(pvalue_threshold), color='red', linestyle='--', alpha=0.7)
        plt.axvline(logfc_threshold, color='red', linestyle='--', alpha=0.7)
        plt.axvline(-logfc_threshold, color='red', linestyle='--', alpha=0.7)
        
        plt.title("Volcano Plot")
        plt.xlabel("log2 Fold Change")
        plt.ylabel("-log10(adjusted p-value)")
        plt.colorbar(label='Significant')
        
        plot_path = self.output_dir / "volcano_plot.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['volcano_plot'] = str(plot_path)
        
        # 2. MA plot
        plt.figure(figsize=(10, 8))
        plt.scatter(de_results['baseMean'], 
                   de_results['log2FoldChange'], 
                   alpha=0.6, c=significant, cmap='viridis')
        
        plt.axhline(0, color='black', linestyle='-', alpha=0.5)
        plt.axhline(logfc_threshold, color='red', linestyle='--', alpha=0.7)
        plt.axhline(-logfc_threshold, color='red', linestyle='--', alpha=0.7)
        
        plt.title("MA Plot")
        plt.xlabel("Mean Expression")
        plt.ylabel("log2 Fold Change")
        plt.xscale('log')
        plt.colorbar(label='Significant')
        
        plot_path = self.output_dir / "ma_plot.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['ma_plot'] = str(plot_path)
        
        # 3. Heatmap of top differentially expressed genes
        if significant.sum() > 0:
            top_genes = de_results[significant].nlargest(50, 'baseMean')
            
            plt.figure(figsize=(12, 10))
            # This would require the original expression data
            # For now, create a placeholder
            plt.title("Heatmap of Top Differentially Expressed Genes")
            plot_path = self.output_dir / "de_heatmap.png"
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            plot_paths['de_heatmap'] = str(plot_path)
        
        return plot_paths
    
    def create_interactive_plots(self, data: pd.DataFrame, plot_type: str = "scatter") -> str:
        """
        Create interactive plots using Plotly.
        
        Args:
            data: DataFrame with data to plot
            plot_type: Type of plot to create
            
        Returns:
            Path to HTML file with interactive plot
        """
        if plot_type == "scatter":
            fig = px.scatter(data, x=data.columns[0], y=data.columns[1], 
                           title=f"{plot_type.title()} Plot")
        elif plot_type == "heatmap":
            fig = px.imshow(data, title=f"{plot_type.title()} Plot")
        elif plot_type == "box":
            fig = px.box(data, title=f"{plot_type.title()} Plot")
        else:
            fig = px.scatter(data, title="Interactive Plot")
        
        plot_path = self.output_dir / f"interactive_{plot_type}.html"
        fig.write_html(str(plot_path))
        
        return str(plot_path)
    
    def create_summary_report(self, plot_paths: Dict[str, str]) -> str:
        """
        Create a summary report with all plots.
        
        Args:
            plot_paths: Dictionary of plot file paths
            
        Returns:
            Path to generated report
        """
        report_path = self.output_dir / "visualization_report.html"
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>RNA-seq Visualization Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .plot {{ text-align: center; margin: 20px 0; }}
                img {{ max-width: 100%; height: auto; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>RNA-seq Visualization Report</h1>
                <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
        """
        
        # Group plots by category
        categories = {
            'Quality Control': [k for k in plot_paths.keys() if 'quality' in k.lower() or 'gc' in k.lower()],
            'Alignment': [k for k in plot_paths.keys() if 'mapping' in k.lower() or 'reads' in k.lower()],
            'Expression': [k for k in plot_paths.keys() if 'expression' in k.lower() or 'correlation' in k.lower() or 'pca' in k.lower()],
            'Differential Expression': [k for k in plot_paths.keys() if 'volcano' in k.lower() or 'ma' in k.lower() or 'de' in k.lower()]
        }
        
        for category, plots in categories.items():
            if plots:
                html_content += f"""
                <div class="section">
                    <h2>{category}</h2>
                """
                
                for plot_name in plots:
                    plot_path = plot_paths[plot_name]
                    if plot_path.endswith('.html'):
                        # Interactive plot
                        html_content += f"""
                        <div class="plot">
                            <h3>{plot_name.replace('_', ' ').title()}</h3>
                            <iframe src="{plot_path}" width="100%" height="600px" frameborder="0"></iframe>
                        </div>
                        """
                    else:
                        # Static plot
                        html_content += f"""
                        <div class="plot">
                            <h3>{plot_name.replace('_', ' ').title()}</h3>
                            <img src="{plot_path}" alt="{plot_name}">
                        </div>
                        """
                
                html_content += "</div>"
        
        html_content += """
        </body>
        </html>
        """
        
        with open(report_path, 'w') as f:
            f.write(html_content)
            
        logger.info(f"Visualization report generated: {report_path}")
        return str(report_path) 