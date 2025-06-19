#!/usr/bin/env python3
"""
RNA-seq Analysis Pipeline
Command-line interface for RNA-seq data analysis.
"""

import os
import sys
import yaml
import argparse
import subprocess
from pathlib import Path
from loguru import logger
import pandas as pd

# Add src to path
sys.path.append(str(Path(__file__).parent / 'src'))

from data_retrieval import NCBIDownloader
from preprocessing import QualityControl
from alignment import STARAligner
from visualization import PlotGenerator


class RNAseqPipeline:
    """RNA-seq analysis pipeline."""
    
    def __init__(self, config_file: str):
        """
        Initialize the pipeline.
        
        Args:
            config_file: Path to configuration file
        """
        self.config = self.load_config(config_file)
        self.setup_logging()
        
    def load_config(self, config_file: str) -> dict:
        """Load configuration from YAML file."""
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def setup_logging(self):
        """Setup logging configuration."""
        log_dir = Path("logs")
        log_dir.mkdir(exist_ok=True)
        
        logger.add(
            log_dir / "pipeline.log",
            rotation="1 day",
            retention="7 days",
            level="INFO"
        )
        
    def run_pipeline(self):
        """Run the complete RNA-seq analysis pipeline."""
        logger.info("Starting RNA-seq analysis pipeline")
        
        pipeline_steps = self.config['pipeline']['steps']
        
        for step in pipeline_steps:
            logger.info(f"Running step: {step}")
            
            try:
                if step == "data_retrieval":
                    self.run_data_retrieval()
                elif step == "quality_control":
                    self.run_quality_control()
                elif step == "alignment":
                    self.run_alignment()
                elif step == "quantification":
                    self.run_quantification()
                elif step == "differential_expression":
                    self.run_differential_expression()
                elif step == "visualization":
                    self.run_visualization()
                else:
                    logger.warning(f"Unknown pipeline step: {step}")
                    
            except Exception as e:
                logger.error(f"Error in step {step}: {e}")
                if not self.config['pipeline'].get('continue_on_error', False):
                    raise
                    
        logger.info("Pipeline completed successfully")
        
    def run_data_retrieval(self):
        """Run data retrieval step."""
        logger.info("Starting data retrieval")
        
        config = self.config['data_retrieval']
        downloader = NCBIDownloader(output_dir=config['output_directory'])
        
        all_datasets = []
        for query in config['search_queries']:
            logger.info(f"Searching for: {query}")
            datasets = downloader.search_sra(
                query, 
                max_results=config['max_results_per_query']
            )
            all_datasets.extend(datasets)
            
        # Remove duplicates
        unique_datasets = {}
        for dataset in all_datasets:
            unique_datasets[dataset['run_id']] = dataset
            
        datasets = list(unique_datasets.values())
        
        if datasets:
            # Download datasets
            run_ids = [d['run_id'] for d in datasets]
            downloaded_files = downloader.download_sra_data(run_ids)
            
            # Create sample sheet
            sample_sheet = downloader.create_sample_sheet(datasets)
            
            logger.info(f"Downloaded {len(downloaded_files)} files")
            logger.info(f"Sample sheet created: {sample_sheet}")
        else:
            logger.warning("No datasets found")
            
    def run_quality_control(self):
        """Run quality control step."""
        logger.info("Starting quality control")
        
        config = self.config['quality_control']
        qc = QualityControl(output_dir=config['output_directory'])
        
        # Find FASTQ files
        data_dir = Path(self.config['data_retrieval']['output_directory'])
        fastq_files = list(data_dir.rglob("*.fastq*"))
        
        if not fastq_files:
            logger.warning("No FASTQ files found for quality control")
            return
            
        logger.info(f"Found {len(fastq_files)} FASTQ files")
        
        # Run FastQC
        qc_outputs = qc.run_fastqc(fastq_files, threads=config['threads'])
        
        if qc_outputs:
            # Run MultiQC
            multiqc_report = qc.run_multiqc(qc_outputs)
            
            # Parse results and create plots
            qc_df = qc.parse_fastqc_results(qc_outputs)
            plot_paths = qc.create_qc_summary_plots(qc_df)
            qc_report = qc.generate_qc_report(qc_df, plot_paths)
            
            logger.info(f"Quality control completed. Report: {qc_report}")
        else:
            logger.warning("No QC outputs generated")
            
    def run_alignment(self):
        """Run alignment step."""
        logger.info("Starting alignment")
        
        config = self.config['alignment']
        genome_dir = config['genome_directory']
        
        # Check if genome index exists
        if not Path(genome_dir).exists():
            logger.info("Building genome index")
            self.build_genome_index()
            
        aligner = STARAligner(
            genome_dir=genome_dir,
            output_dir=config['output_directory'],
            threads=config['threads'],
            memory_limit=config['memory_limit']
        )
        
        # Find FASTQ files
        data_dir = Path(self.config['data_retrieval']['output_directory'])
        fastq_files = list(data_dir.rglob("*.fastq*"))
        
        if not fastq_files:
            logger.warning("No FASTQ files found for alignment")
            return
            
        # Get sample names from configuration
        samples = self.config['samples']
        sample_names = [s['name'] for s in samples]
        
        # Determine if paired-end or single-end
        if len(fastq_files) == len(sample_names):
            # Single-end
            bam_files = aligner.align_single_reads(fastq_files, sample_names)
        elif len(fastq_files) == 2 * len(sample_names):
            # Paired-end
            read1_files = fastq_files[:len(sample_names)]
            read2_files = fastq_files[len(sample_names):]
            bam_files = aligner.align_paired_reads(read1_files, read2_files, sample_names)
        else:
            logger.error("Mismatch between number of files and samples")
            return
            
        if bam_files:
            # Create alignment summary
            alignment_summary = aligner.create_alignment_summary(bam_files)
            alignment_report = aligner.generate_alignment_report(alignment_summary)
            
            logger.info(f"Alignment completed. Report: {alignment_report}")
        else:
            logger.warning("No BAM files generated")
            
    def build_genome_index(self):
        """Build STAR genome index."""
        logger.info("Building STAR genome index")
        
        ref_config = self.config['reference_genome']
        genome_fasta = ref_config['genome_fasta']
        gtf_file = ref_config['gtf_file']
        index_dir = ref_config['index_directory']
        
        # Check if reference files exist
        if not Path(genome_fasta).exists():
            logger.error(f"Genome FASTA file not found: {genome_fasta}")
            return False
            
        if not Path(gtf_file).exists():
            logger.error(f"GTF file not found: {gtf_file}")
            return False
            
        aligner = STARAligner(genome_dir=index_dir)
        success = aligner.build_genome_index(genome_fasta, gtf_file, index_dir)
        
        if success:
            logger.info("Genome index built successfully")
        else:
            logger.error("Failed to build genome index")
            
        return success
        
    def run_quantification(self):
        """Run quantification step."""
        logger.info("Starting quantification")
        
        config = self.config['quantification']
        
        # Find BAM files
        alignment_dir = Path(self.config['alignment']['output_directory'])
        bam_files = list(alignment_dir.rglob("*.bam"))
        
        if not bam_files:
            logger.warning("No BAM files found for quantification")
            return
            
        # Run featureCounts
        gtf_file = self.config['reference_genome']['gtf_file']
        output_file = Path(config['output_directory']) / "counts.txt"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            'featureCounts',
            '-T', str(config['threads']),
            '-t', config['feature_type'],
            '-g', config['attribute'],
            '-o', str(output_file),
            '-a', gtf_file
        ] + bam_files
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info(f"Quantification completed: {output_file}")
            else:
                logger.error(f"Quantification failed: {result.stderr}")
                
        except Exception as e:
            logger.error(f"Error running quantification: {e}")
            
    def run_differential_expression(self):
        """Run differential expression analysis."""
        logger.info("Starting differential expression analysis")
        
        config = self.config['differential_expression']
        
        # This would typically involve R and DESeq2
        # For now, create a placeholder
        logger.info("Differential expression analysis would be implemented here")
        
    def run_visualization(self):
        """Run visualization step."""
        logger.info("Starting visualization")
        
        config = self.config['visualization']
        plot_generator = PlotGenerator(output_dir=config['output_directory'])
        
        # Create various plots based on available data
        plot_paths = {}
        
        # Quality control plots
        qc_dir = Path(self.config['quality_control']['output_directory'])
        if qc_dir.exists():
            # This would require parsing QC data
            logger.info("Creating QC plots")
            
        # Alignment plots
        alignment_dir = Path(self.config['alignment']['output_directory'])
        if alignment_dir.exists():
            # This would require parsing alignment data
            logger.info("Creating alignment plots")
            
        # Expression plots
        quantification_dir = Path(self.config['quantification']['output_directory'])
        counts_file = quantification_dir / "counts.txt"
        if counts_file.exists():
            logger.info("Creating expression plots")
            
        # Generate summary report
        if plot_paths:
            summary_report = plot_generator.create_summary_report(plot_paths)
            logger.info(f"Visualization completed. Report: {summary_report}")
        else:
            logger.info("No plots generated")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="RNA-seq Analysis Pipeline")
    parser.add_argument(
        '--config', 
        type=str, 
        default='config/soybean_analysis.yaml',
        help='Configuration file path'
    )
    parser.add_argument(
        '--step',
        type=str,
        choices=['data_retrieval', 'quality_control', 'alignment', 
                'quantification', 'differential_expression', 'visualization'],
        help='Run specific pipeline step'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without executing'
    )
    
    args = parser.parse_args()
    
    # Check if config file exists
    if not Path(args.config).exists():
        logger.error(f"Configuration file not found: {args.config}")
        sys.exit(1)
        
    # Initialize pipeline
    pipeline = RNAseqPipeline(args.config)
    
    if args.dry_run:
        logger.info("Dry run mode - no actual execution")
        logger.info(f"Configuration: {args.config}")
        if args.step:
            logger.info(f"Would run step: {args.step}")
        else:
            logger.info("Would run all pipeline steps")
        return
        
    # Run pipeline
    try:
        if args.step:
            # Run specific step
            if args.step == "data_retrieval":
                pipeline.run_data_retrieval()
            elif args.step == "quality_control":
                pipeline.run_quality_control()
            elif args.step == "alignment":
                pipeline.run_alignment()
            elif args.step == "quantification":
                pipeline.run_quantification()
            elif args.step == "differential_expression":
                pipeline.run_differential_expression()
            elif args.step == "visualization":
                pipeline.run_visualization()
        else:
            # Run complete pipeline
            pipeline.run_pipeline()
            
    except KeyboardInterrupt:
        logger.info("Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 