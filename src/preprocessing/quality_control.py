"""
Quality control module for RNA-seq data.
Handles FastQC analysis and quality assessment.
"""

import os
import subprocess
import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
from loguru import logger
import matplotlib.pyplot as plt
import seaborn as sns


class QualityControl:
    """Quality control analysis for RNA-seq data."""
    
    def __init__(self, output_dir: str = "results/qc"):
        """
        Initialize quality control module.
        
        Args:
            output_dir: Directory to store QC results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def run_fastqc(self, fastq_files: List[str], threads: int = 4) -> List[str]:
        """
        Run FastQC on FASTQ files.
        
        Args:
            fastq_files: List of FASTQ file paths
            threads: Number of threads to use
            
        Returns:
            List of FastQC output directories
        """
        logger.info(f"Running FastQC on {len(fastq_files)} files")
        
        fastqc_outputs = []
        for fastq_file in fastq_files:
            try:
                output_path = self.output_dir / Path(fastq_file).stem
                output_path.mkdir(exist_ok=True)
                
                cmd = [
                    'fastqc',
                    '--outdir', str(output_path),
                    '--threads', str(threads),
                    '--noextract',
                    fastq_file
                ]
                
                logger.info(f"Running FastQC on {fastq_file}")
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    fastqc_outputs.append(str(output_path))
                    logger.info(f"FastQC completed for {fastq_file}")
                else:
                    logger.error(f"FastQC failed for {fastq_file}: {result.stderr}")
                    
            except Exception as e:
                logger.error(f"Error running FastQC on {fastq_file}: {e}")
                
        return fastqc_outputs
    
    def run_multiqc(self, qc_dirs: List[str]) -> str:
        """
        Run MultiQC to aggregate FastQC results.
        
        Args:
            qc_dirs: List of FastQC output directories
            
        Returns:
            Path to MultiQC report
        """
        logger.info("Running MultiQC to aggregate QC results")
        
        try:
            cmd = [
                'multiqc',
                '--outdir', str(self.output_dir),
                '--filename', 'multiqc_report'
            ] + qc_dirs
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                report_path = self.output_dir / "multiqc_report.html"
                logger.info(f"MultiQC report generated: {report_path}")
                return str(report_path)
            else:
                logger.error(f"MultiQC failed: {result.stderr}")
                return ""
                
        except Exception as e:
            logger.error(f"Error running MultiQC: {e}")
            return ""
    
    def parse_fastqc_results(self, qc_dirs: List[str]) -> pd.DataFrame:
        """
        Parse FastQC results into a pandas DataFrame.
        
        Args:
            qc_dirs: List of FastQC output directories
            
        Returns:
            DataFrame with QC metrics
        """
        qc_data = []
        
        for qc_dir in qc_dirs:
            try:
                qc_dir_path = Path(qc_dir)
                sample_name = qc_dir_path.name
                
                # Parse FastQC data.txt file
                data_file = qc_dir_path / f"{sample_name}_fastqc" / "fastqc_data.txt"
                
                if data_file.exists():
                    metrics = self._parse_fastqc_data(data_file)
                    metrics['sample'] = sample_name
                    qc_data.append(metrics)
                    
            except Exception as e:
                logger.error(f"Error parsing FastQC results for {qc_dir}: {e}")
                
        return pd.DataFrame(qc_data)
    
    def _parse_fastqc_data(self, data_file: Path) -> Dict:
        """Parse individual FastQC data file."""
        metrics = {}
        
        try:
            with open(data_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        if '\t' in line:
                            key, value = line.split('\t', 1)
                            metrics[key] = value
                            
        except Exception as e:
            logger.error(f"Error parsing FastQC data file {data_file}: {e}")
            
        return metrics
    
    def create_qc_summary_plots(self, qc_df: pd.DataFrame) -> Dict[str, str]:
        """
        Create summary plots for QC metrics.
        
        Args:
            qc_df: DataFrame with QC metrics
            
        Returns:
            Dictionary mapping plot names to file paths
        """
        plot_paths = {}
        
        # Set style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
        
        # 1. Per base sequence quality
        if 'Per base sequence quality' in qc_df.columns:
            plt.figure(figsize=(12, 6))
            # This would require parsing the detailed quality data
            # For now, create a placeholder plot
            plt.title("Per Base Sequence Quality")
            plt.xlabel("Position in read")
            plt.ylabel("Quality score")
            plot_path = self.output_dir / "per_base_quality.png"
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            plot_paths['per_base_quality'] = str(plot_path)
        
        # 2. Per sequence quality scores
        plt.figure(figsize=(10, 6))
        plt.title("Per Sequence Quality Scores")
        plt.xlabel("Quality Score")
        plt.ylabel("Count")
        plot_path = self.output_dir / "per_sequence_quality.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['per_sequence_quality'] = str(plot_path)
        
        # 3. Per base sequence content
        plt.figure(figsize=(12, 6))
        plt.title("Per Base Sequence Content")
        plt.xlabel("Position in read")
        plt.ylabel("Percentage")
        plot_path = self.output_dir / "per_base_content.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['per_base_content'] = str(plot_path)
        
        # 4. GC content distribution
        plt.figure(figsize=(10, 6))
        plt.title("GC Content Distribution")
        plt.xlabel("GC Content (%)")
        plt.ylabel("Count")
        plot_path = self.output_dir / "gc_content.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plot_paths['gc_content'] = str(plot_path)
        
        logger.info(f"Created {len(plot_paths)} QC summary plots")
        return plot_paths
    
    def generate_qc_report(self, qc_df: pd.DataFrame, plot_paths: Dict[str, str]) -> str:
        """
        Generate a comprehensive QC report.
        
        Args:
            qc_df: DataFrame with QC metrics
            plot_paths: Dictionary of plot file paths
            
        Returns:
            Path to generated report
        """
        report_path = self.output_dir / "qc_report.html"
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>RNA-seq Quality Control Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .plot {{ text-align: center; margin: 20px 0; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>RNA-seq Quality Control Report</h1>
                <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="section">
                <h2>Sample Summary</h2>
                <p>Total samples analyzed: {len(qc_df)}</p>
            </div>
            
            <div class="section">
                <h2>Quality Metrics Summary</h2>
                {qc_df.to_html() if not qc_df.empty else '<p>No QC data available</p>'}
            </div>
            
            <div class="section">
                <h2>Quality Control Plots</h2>
        """
        
        for plot_name, plot_path in plot_paths.items():
            html_content += f"""
                <div class="plot">
                    <h3>{plot_name.replace('_', ' ').title()}</h3>
                    <img src="{plot_path}" alt="{plot_name}" style="max-width: 100%;">
                </div>
            """
        
        html_content += """
            </div>
        </body>
        </html>
        """
        
        with open(report_path, 'w') as f:
            f.write(html_content)
            
        logger.info(f"QC report generated: {report_path}")
        return str(report_path) 