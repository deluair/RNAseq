"""
Demonstration script for running a full RNA-seq analysis pipeline,
including differential expression and RNA splicing analysis.

This script showcases how to use the various analysis modules of the platform
in a cohesive workflow.

Author: RNA-seq Analysis Platform
"""

import os
import pandas as pd
from loguru import logger
import json
import numpy as np
import base64
from pathlib import Path

# Import analysis modules
from src.analysis.deseq2_analyzer import DESeq2Analyzer
from src.analysis.splicing_analyzer import SplicingAnalyzer, create_mock_bam_files

def image_to_base64(img_path):
    """Converts an image file to a base64 string for HTML embedding."""
    if not img_path or not Path(img_path).exists():
        return None
    with open(img_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def generate_combined_html_report(de_results, splicing_results, de_dir, splicing_dir, output_dir):
    """Generates a single HTML report summarizing all analysis results."""
    
    report_path = Path(output_dir) / "full_analysis_report.html"
    de_dir = Path(de_dir)
    splicing_dir = Path(splicing_dir)

    # Get base64 encoded images
    volcano_b64 = image_to_base64(de_dir / "volcano_plot.png")
    ma_b64 = image_to_base64(de_dir / "ma_plot.png")
    splicing_overview_b64 = image_to_base64(splicing_results.get("overview_plot"))
    diff_splicing_b64 = image_to_base64(splicing_results.get("differential_plot"))

    # Start HTML content
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Full RNA-Seq Analysis Report</title>
        <style>
            body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; margin: 0 auto; max-width: 1200px; padding: 20px; color: #333; }
            h1, h2, h3 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
            h1 { text-align: center; }
            .container { display: flex; flex-wrap: wrap; gap: 20px; }
            .card { background: #fff; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); padding: 20px; flex: 1; min-width: 400px; }
            .metric-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px; text-align: center; }
            .metric { background-color: #f8f9fa; padding: 15px; border-radius: 5px; }
            .metric .value { font-size: 2em; font-weight: bold; color: #3498db; }
            .metric .label { font-size: 0.9em; color: #555; }
            table { border-collapse: collapse; width: 100%; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
            th { background-color: #3498db; color: white; }
            .plot { text-align: center; margin: 20px 0; }
            img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; }
            .footer { text-align: center; margin-top: 40px; color: #777; font-size: 0.9em; }
        </style>
    </head>
    <body>
        <h1>Full RNA-Seq Analysis Report</h1>
    """

    # Section 1: Differential Expression
    html += """
        <h2>Differential Gene Expression Analysis</h2>
        <div class="card">
            <h3>Summary Statistics</h3>
    """
    if de_results and de_results.get('summary_stats'):
        stats = de_results['summary_stats']
        html += f"""
            <div class="metric-grid">
                <div class="metric"><div class="value">{stats['total_genes_tested']:,}</div><div class="label">Genes Tested</div></div>
                <div class="metric"><div class="value">{stats['total_significant']:,}</div><div class="label">Significant</div></div>
                <div class="metric"><div class="value" style="color: #e74c3c;">{stats['upregulated']:,}</div><div class="label">Upregulated</div></div>
                <div class="metric"><div class="value" style="color: #2980b9;">{stats['downregulated']:,}</div><div class="label">Downregulated</div></div>
            </div>
        """
    else:
        html += "<p>No differential expression results available.</p>"
    
    html += "<h3>Plots</h3><div class='container'>"
    if volcano_b64:
        html += f"<div class='plot card'><h4>Volcano Plot</h4><img src='data:image/png;base64,{volcano_b64}'></div>"
    if ma_b64:
        html += f"<div class='plot card'><h4>MA Plot</h4><img src='data:image/png;base64,{ma_b64}'></div>"
    html += "</div>"
    
    html += """
        <h3>Downloads</h3>
        <ul>
            <li><a href='./differential_expression/all_de_results.csv' download>Full Results (CSV)</a></li>
            <li><a href='./differential_expression/significant_genes.csv' download>Significant Genes (CSV)</a></li>
        </ul>
    </div>
    """

    # Section 2: RNA Splicing
    html += """
        <h2>RNA Splicing Analysis</h2>
        <div class="card">
            <h3>Summary Statistics</h3>
    """
    if splicing_results and splicing_results.get('summary'):
        s_stats = splicing_results['summary']
        html += f"""
            <div class="metric-grid">
                <div class="metric"><div class="value">{s_stats['total_junctions']:,}</div><div class="label">Total Junctions</div></div>
                <div class="metric"><div class="value">{s_stats['as_events_count']['SE']:,}</div><div class="label">Skipped Exons</div></div>
                <div class="metric"><div class="value">{s_stats['as_events_count']['RI']:,}</div><div class="label">Retained Introns</div></div>
                <div class="metric"><div class="value">{sum(s_stats['as_events_count'].values()):,}</div><div class="label">Total AS Events</div></div>
            </div>
        """
    else:
        html += "<p>No splicing analysis results available.</p>"
        
    html += "<h3>Plots</h3><div class='container'>"
    if splicing_overview_b64:
        html += f"<div class='plot card'><h4>Splicing Overview</h4><img src='data:image/png;base64,{splicing_overview_b64}'></div>"
    if diff_splicing_b64:
        html += f"<div class='plot card'><h4>Differential Splicing</h4><img src='data:image/png;base64,{diff_splicing_b64}'></div>"
    html += "</div>"

    html += """
        <h3>Downloads & Links</h3>
        <ul>
            <li><a href='./splicing_analysis/splicing_report.html'>Detailed Splicing Report (HTML)</a></li>
            <li><a href='./splicing_analysis/splice_junctions.csv' download>Splice Junctions (CSV)</a></li>
            <li><a href='./splicing_analysis/splicing_efficiency.csv' download>Splicing Efficiency (CSV)</a></li>
        </ul>
    </div>
    """
    
    html += """
        <div class="footer">
            <p>Report generated by the RNA-Seq Analysis Platform.</p>
        </div>
    </body>
    </html>
    """
    
    with open(report_path, "w") as f:
        f.write(html)
        
    return report_path

def run_full_analysis_pipeline():
    """
    Executes a full, mock analysis pipeline from start to finish.
    """
    logger.info("Starting full RNA-seq analysis pipeline demonstration")
    
    # ==========================================================================
    # 1. SETUP & MOCK DATA PREPARATION
    # ==========================================================================
    logger.info("Step 1: Setting up analysis and preparing mock data")
    
    # Define output directory for this full analysis run
    output_dir = "results/full_analysis_example"
    os.makedirs(output_dir, exist_ok=True)
    
    # Define sample information
    n_samples_per_condition = 3
    total_samples = n_samples_per_condition * 2
    
    condition1_samples = [f'sample_{i+1}' for i in range(n_samples_per_condition)]
    condition2_samples = [f'sample_{i+n_samples_per_condition+1}' for i in range(n_samples_per_condition)]
    all_samples = condition1_samples + condition2_samples
    
    # Create mock alignment (BAM) files for splicing analysis
    mock_bam_dir = f"{output_dir}/mock_alignment"
    mock_bam_files = create_mock_bam_files(output_dir=mock_bam_dir, n_samples=total_samples)
    
    # Create realistic mock count data for differential expression analysis
    n_genes = 2000
    np.random.seed(42) # for reproducibility

    # Generate base counts from a negative binomial distribution
    base_counts = np.random.negative_binomial(p=0.1, n=100, size=(n_genes, total_samples))

    # Introduce differential expression for a subset of genes
    n_de_genes = int(n_genes * 0.15) # 15% DE genes
    de_indices = np.random.choice(n_genes, n_de_genes, replace=False)

    n_upregulated = int(n_de_genes * 0.7) # 70% of DE genes are upregulated
    up_indices = de_indices[:n_upregulated]
    down_indices = de_indices[n_upregulated:]
    
    # Define fold changes
    up_fold_changes = np.random.uniform(2.5, 5, size=n_upregulated)
    down_fold_changes = np.random.uniform(0.1, 0.4, size=len(down_indices))

    # Apply fold changes to treatment samples
    treatment_slice = slice(n_samples_per_condition, total_samples)
    
    # Upregulate genes
    base_counts[up_indices, treatment_slice] = (base_counts[up_indices, treatment_slice] * up_fold_changes[:, np.newaxis]).astype(int)
    
    # Downregulate genes
    base_counts[down_indices, treatment_slice] = (base_counts[down_indices, treatment_slice] * down_fold_changes[:, np.newaxis]).astype(int)

    mock_counts = pd.DataFrame(
        data=base_counts,
        index=[f'gene_{i+1}' for i in range(n_genes)],
        columns=all_samples
    )
    mock_counts_file = f"{output_dir}/mock_counts.csv"
    mock_counts.to_csv(mock_counts_file)
    
    logger.info(f"Mock data created in: {output_dir}")
    logger.info(f"Mock count matrix file: {mock_counts_file}")
    logger.info(f"Mock BAM files created in: {mock_bam_dir}")

    # ==========================================================================
    # 2. DIFFERENTIAL EXPRESSION ANALYSIS
    # ==========================================================================
    logger.info("Step 2: Running Differential Expression Analysis")
    
    de_output_dir = f"{output_dir}/differential_expression"
    de_analyzer = DESeq2Analyzer(output_dir=de_output_dir)
    
    # Check R/DESeq2 status
    deseq2_status = {
        "r_available": de_analyzer.r_available,
        "deseq2_available": de_analyzer.deseq2_available
    }
    logger.info(f"DESeq2 status: {deseq2_status}")
    
    # Define analysis parameters
    de_params = {
        'method': 'python_limma',  # Use a Python-based method
        'condition1': 'condition1',
        'condition2': 'condition2',
        'condition1_samples': condition1_samples,
        'condition2_samples': condition2_samples,
        'alpha': 0.05
    }
    
    # Run DE analysis
    de_results = de_analyzer.analyze_differential_expression(
        count_matrix=mock_counts, 
        sample_metadata=pd.DataFrame(index=all_samples, data={'condition': [de_params['condition1']]*n_samples_per_condition + [de_params['condition2']]*n_samples_per_condition}),
        control_condition=de_params['condition1'],
        treatment_condition=de_params['condition2'],
        padj_threshold=de_params['alpha']
    )
    
    if de_results and not de_results['all_results'].empty:
        logger.info("Differential expression analysis completed successfully")
        
        # Generate plots
        de_analyzer.create_de_plots(de_results['all_results'])
        
        logger.info(f"DE results and plots saved to: {de_output_dir}")
        
        # Display top DE genes
        top_de_genes = de_results['significant_genes'].head()
        print("\n--- Top 5 Differentially Expressed Genes ---")
        print(top_de_genes)
        print("------------------------------------------\n")
    else:
        logger.error("Differential expression analysis failed or produced no results.")

    # ==========================================================================
    # 3. RNA SPLICING ANALYSIS
    # ==========================================================================
    logger.info("Step 3: Running RNA Splicing Analysis")
    
    splicing_output_dir = f"{output_dir}/splicing_analysis"
    splicing_analyzer = SplicingAnalyzer(output_dir=splicing_output_dir)
    
    # Run the complete splicing analysis pipeline
    splicing_results = splicing_analyzer.run_complete_analysis(
        bam_files=mock_bam_files,
        condition1_samples=[os.path.basename(f).replace('_Aligned.sortedByCoord.out.bam', '') for f in mock_bam_files[:n_samples_per_condition]],
        condition2_samples=[os.path.basename(f).replace('_Aligned.sortedByCoord.out.bam', '') for f in mock_bam_files[n_samples_per_condition:]]
    )
    
    logger.info("RNA splicing analysis completed successfully")
    logger.info(f"Splicing analysis results saved to: {splicing_output_dir}")
    
    # Display summary of splicing analysis
    print("\n--- Splicing Analysis Summary ---")
    print(json.dumps(splicing_results.get('summary', {}), indent=2))
    print("---------------------------------\n")

    # ==========================================================================
    # 4. GENERATE COMBINED HTML REPORT
    # ==========================================================================
    logger.info("Step 4: Generating combined HTML report")
    combined_report_path = generate_combined_html_report(
        de_results, splicing_results, de_output_dir, splicing_output_dir, output_dir
    )
    logger.info(f"Combined HTML report generated at: {combined_report_path}")


    # ==========================================================================
    # 5. FINAL SUMMARY
    # ==========================================================================
    logger.info("Step 5: Generating final summary")
    
    print("\n========================================================")
    print("      FULL RNA-SEQ ANALYSIS PIPELINE COMPLETE")
    print("========================================================")
    print(f"\nAll results for this demonstration run are located in:\n  {os.path.abspath(output_dir)}\n")
    print("Key outputs include:")
    print(f"  - Combined HTML Report: {os.path.abspath(combined_report_path)}")
    print(f"  - Differential Expression Results: {os.path.abspath(de_output_dir)}")
    print(f"    - DE Results Table: {os.path.join(os.path.abspath(de_output_dir), 'all_de_results.csv')}")
    print(f"    - Volcano Plot: {os.path.join(os.path.abspath(de_output_dir), 'volcano_plot.png')}")
    print(f"    - MA Plot: {os.path.join(os.path.abspath(de_output_dir), 'ma_plot.png')}")
    print("\n")
    print(f"  - Splicing Analysis Results: {os.path.abspath(splicing_output_dir)}")
    print(f"    - Splicing HTML Report: {os.path.abspath(splicing_results.get('report', ''))}")
    print(f"    - Splicing Overview Plot: {os.path.abspath(splicing_results.get('overview_plot', ''))}")
    if splicing_results.get('differential_plot'):
        print(f"    - Differential Splicing Plot: {os.path.abspath(splicing_results.get('differential_plot'))}")
    print("\nTo view the results, please check the files listed above.")
    print("========================================================\n")

if __name__ == "__main__":
    # Configure logger
    log_file = "logs/full_analysis_demo.log"
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger.add(log_file, rotation="10 MB", level="INFO")
    
    run_full_analysis_pipeline() 