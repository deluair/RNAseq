#!/usr/bin/env python3
"""
RNA-seq Analysis Platform Demo
Demonstrates the key features of the platform.
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from loguru import logger

# Add src to path
sys.path.append(str(Path(__file__).parent / 'src'))

from data_retrieval import NCBIDownloader
from preprocessing import QualityControl
from alignment import STARAligner
from visualization import PlotGenerator


def create_demo_data():
    """Create demo data for visualization."""
    logger.info("Creating demo data")
    
    # Create demo expression data
    np.random.seed(42)
    genes = [f"Gene_{i:04d}" for i in range(1000)]
    samples = ["Control_1", "Control_2", "Control_3", "Treatment_1", "Treatment_2", "Treatment_3"]
    
    # Generate expression data
    expression_data = pd.DataFrame(
        np.random.poisson(100, size=(len(genes), len(samples))),
        index=genes,
        columns=samples
    )
    
    # Add some differential expression
    treatment_genes = genes[:100]
    for gene in treatment_genes:
        expression_data.loc[gene, ["Treatment_1", "Treatment_2", "Treatment_3"]] *= 2
    
    # Save demo data
    demo_dir = Path("demo_data")
    demo_dir.mkdir(exist_ok=True)
    
    expression_data.to_csv(demo_dir / "expression_data.csv")
    
    # Create demo DE results
    de_results = pd.DataFrame({
        'gene': genes,
        'baseMean': expression_data.mean(axis=1),
        'log2FoldChange': np.random.normal(0, 1, len(genes)),
        'padj': np.random.uniform(0, 1, len(genes))
    })
    
    # Make some genes significant
    de_results.iloc[:100, de_results.columns.get_loc('log2FoldChange')] = np.random.normal(2, 0.5, 100)
    de_results.iloc[:100, de_results.columns.get_loc('padj')] = np.random.uniform(0, 0.05, 100)
    
    de_results.to_csv(demo_dir / "de_results.csv", index=False)
    
    return expression_data, de_results


def demo_visualization():
    """Demonstrate visualization capabilities."""
    logger.info("Running visualization demo")
    
    # Create demo data
    expression_data, de_results = create_demo_data()
    
    # Initialize plot generator
    plot_generator = PlotGenerator(output_dir="demo_results/plots")
    
    # Create various plots
    plot_paths = {}
    
    # Expression plots
    expression_plots = plot_generator.create_expression_plots(expression_data)
    plot_paths.update(expression_plots)
    
    # DE plots
    de_plots = plot_generator.create_differential_expression_plots(de_results)
    plot_paths.update(de_plots)
    
    # Create summary report
    summary_report = plot_generator.create_summary_report(plot_paths)
    
    logger.info(f"Demo plots created: {len(plot_paths)}")
    logger.info(f"Summary report: {summary_report}")
    
    return plot_paths, summary_report


def demo_ncbi_search():
    """Demonstrate NCBI data search."""
    logger.info("Running NCBI search demo")
    
    downloader = NCBIDownloader()
    
    # Search for soybean datasets
    datasets = downloader.get_soybean_datasets(limit=5)
    
    logger.info(f"Found {len(datasets)} soybean datasets")
    
    for dataset in datasets:
        logger.info(f"  - {dataset['run_id']}: {dataset['title']}")
    
    return datasets


def demo_quality_control():
    """Demonstrate quality control functionality."""
    logger.info("Running quality control demo")
    
    qc = QualityControl(output_dir="demo_results/qc")
    
    # Create demo QC data
    qc_data = pd.DataFrame({
        'sample': ['Sample_1', 'Sample_2', 'Sample_3'],
        'Total Sequences': [1000000, 950000, 1050000],
        'Sequences flagged as poor quality': [50000, 45000, 55000],
        'GC content': [45.2, 44.8, 45.5]
    })
    
    # Create QC plots
    plot_paths = qc.create_qc_summary_plots(qc_data)
    
    # Generate QC report
    qc_report = qc.generate_qc_report(qc_data, plot_paths)
    
    logger.info(f"QC demo completed. Report: {qc_report}")
    
    return plot_paths, qc_report


def demo_alignment():
    """Demonstrate alignment functionality."""
    logger.info("Running alignment demo")
    
    # Create demo alignment stats
    alignment_stats = {
        'Sample_1': {
            'total_reads': 1000000,
            'mapped_reads': 950000,
            'mapped_percent': 95.0
        },
        'Sample_2': {
            'total_reads': 950000,
            'mapped_reads': 900000,
            'mapped_percent': 94.7
        },
        'Sample_3': {
            'total_reads': 1050000,
            'mapped_reads': 1000000,
            'mapped_percent': 95.2
        }
    }
    
    # Initialize aligner for demo
    aligner = STARAligner(genome_dir="genome_index", output_dir="demo_results/alignment")
    
    # Generate alignment report
    alignment_report = aligner.generate_alignment_report(alignment_stats)
    
    logger.info(f"Alignment demo completed. Report: {alignment_report}")
    
    return alignment_stats, alignment_report


def main():
    """Run the complete demo."""
    logger.info("Starting RNA-seq Analysis Platform Demo")
    
    # Create demo directories
    demo_dirs = ["demo_data", "demo_results", "demo_results/plots", 
                "demo_results/qc", "demo_results/alignment"]
    for dir_path in demo_dirs:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    try:
        # Demo 1: NCBI Search
        logger.info("\n" + "="*50)
        logger.info("DEMO 1: NCBI Data Search")
        logger.info("="*50)
        datasets = demo_ncbi_search()
        
        # Demo 2: Quality Control
        logger.info("\n" + "="*50)
        logger.info("DEMO 2: Quality Control")
        logger.info("="*50)
        qc_plots, qc_report = demo_quality_control()
        
        # Demo 3: Alignment
        logger.info("\n" + "="*50)
        logger.info("DEMO 3: Alignment")
        logger.info("="*50)
        alignment_stats, alignment_report = demo_alignment()
        
        # Demo 4: Visualization
        logger.info("\n" + "="*50)
        logger.info("DEMO 4: Visualization")
        logger.info("="*50)
        viz_plots, viz_report = demo_visualization()
        
        # Summary
        logger.info("\n" + "="*50)
        logger.info("DEMO SUMMARY")
        logger.info("="*50)
        logger.info(f"✓ NCBI datasets found: {len(datasets)}")
        logger.info(f"✓ QC plots created: {len(qc_plots)}")
        logger.info(f"✓ Alignment samples processed: {len(alignment_stats)}")
        logger.info(f"✓ Visualization plots created: {len(viz_plots)}")
        logger.info(f"✓ Reports generated: 3")
        
        logger.info("\nDemo completed successfully!")
        logger.info("Check the demo_results/ directory for generated files.")
        
    except Exception as e:
        logger.error(f"Demo failed: {e}")
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 