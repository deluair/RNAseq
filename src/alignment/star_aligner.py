"""
STAR aligner module for RNA-seq analysis.
Handles RNA-seq read alignment using STAR aligner.
"""

import os
import subprocess
import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from loguru import logger
import pandas as pd


class STARAligner:
    """STAR aligner for RNA-seq data."""
    
    def __init__(self, 
                 genome_dir: str,
                 output_dir: str = "results/alignment",
                 threads: int = 4,
                 memory_limit: str = "32G"):
        """
        Initialize STAR aligner.
        
        Args:
            genome_dir: Directory containing STAR genome index
            output_dir: Directory to store alignment results
            threads: Number of threads to use
            memory_limit: Memory limit for STAR
        """
        self.genome_dir = Path(genome_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        self.memory_limit = memory_limit
        
        # Check if STAR is available
        self._check_star_installation()
        
    def _check_star_installation(self):
        """Check if STAR is installed and accessible."""
        try:
            result = subprocess.run(['STAR', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"STAR version: {result.stdout.strip()}")
            else:
                logger.warning("STAR not found or not working properly")
        except FileNotFoundError:
            logger.error("STAR not found. Please install STAR aligner.")
            
    def build_genome_index(self, 
                          genome_fasta: str,
                          gtf_file: str,
                          index_dir: str = "genome_index") -> bool:
        """
        Build STAR genome index.
        
        Args:
            genome_fasta: Path to genome FASTA file
            gtf_file: Path to GTF annotation file
            index_dir: Directory to store genome index
            
        Returns:
            True if successful, False otherwise
        """
        logger.info("Building STAR genome index")
        
        index_path = Path(index_dir)
        index_path.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            'STAR',
            '--runMode', 'genomeGenerate',
            '--genomeDir', str(index_path),
            '--genomeFastaFiles', genome_fasta,
            '--sjdbGTFfile', gtf_file,
            '--runThreadN', str(self.threads),
            '--genomeSAindexNbases', '14'  # For large genomes
        ]
        
        try:
            logger.info("Running STAR genome index generation...")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info(f"Genome index built successfully in {index_path}")
                self.genome_dir = index_path
                return True
            else:
                logger.error(f"STAR genome index generation failed: {result.stderr}")
                return False
                
        except Exception as e:
            logger.error(f"Error building genome index: {e}")
            return False
    
    def align_paired_reads(self, 
                          read1_files: List[str],
                          read2_files: List[str],
                          sample_names: List[str]) -> List[str]:
        """
        Align paired-end RNA-seq reads.
        
        Args:
            read1_files: List of R1 FASTQ files
            read2_files: List of R2 FASTQ files
            sample_names: List of sample names
            
        Returns:
            List of BAM file paths
        """
        logger.info(f"Aligning {len(sample_names)} samples with STAR")
        
        bam_files = []
        
        for i, (r1, r2, sample) in enumerate(zip(read1_files, read2_files, sample_names)):
            try:
                sample_output_dir = self.output_dir / sample
                sample_output_dir.mkdir(exist_ok=True)
                
                cmd = [
                    'STAR',
                    '--genomeDir', str(self.genome_dir),
                    '--readFilesIn', r1, r2,
                    '--outFileNamePrefix', str(sample_output_dir / f"{sample}_"),
                    '--runThreadN', str(self.threads),
                    '--outSAMtype', 'BAM', 'SortedByCoordinate',
                    '--outBAMsortingThreadN', str(self.threads),
                    '--limitBAMsortRAM', self.memory_limit,
                    '--outSAMunmapped', 'Within',
                    '--outSAMattributes', 'Standard',
                    '--quantMode', 'GeneCounts',
                    '--twopassMode', 'Basic'
                ]
                
                # Handle compressed files
                if r1.endswith('.gz'):
                    cmd.extend(['--readFilesCommand', 'zcat'])
                
                logger.info(f"Aligning sample {sample} ({i+1}/{len(sample_names)})")
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    bam_file = sample_output_dir / f"{sample}_Aligned.sortedByCoord.out.bam"
                    if bam_file.exists():
                        bam_files.append(str(bam_file))
                        logger.info(f"Successfully aligned {sample}")
                    else:
                        logger.error(f"BAM file not found for {sample}")
                else:
                    logger.error(f"STAR alignment failed for {sample}: {result.stderr}")
                    
            except Exception as e:
                logger.error(f"Error aligning {sample}: {e}")
                
        return bam_files
    
    def align_single_reads(self, 
                          fastq_files: List[str],
                          sample_names: List[str]) -> List[str]:
        """
        Align single-end RNA-seq reads.
        
        Args:
            fastq_files: List of FASTQ files
            sample_names: List of sample names
            
        Returns:
            List of BAM file paths
        """
        logger.info(f"Aligning {len(sample_names)} single-end samples with STAR")
        
        bam_files = []
        
        for i, (fastq_file, sample) in enumerate(zip(fastq_files, sample_names)):
            try:
                sample_output_dir = self.output_dir / sample
                sample_output_dir.mkdir(exist_ok=True)
                
                cmd = [
                    'STAR',
                    '--genomeDir', str(self.genome_dir),
                    '--readFilesIn', fastq_file,
                    '--outFileNamePrefix', str(sample_output_dir / f"{sample}_"),
                    '--runThreadN', str(self.threads),
                    '--outSAMtype', 'BAM', 'SortedByCoordinate',
                    '--outBAMsortingThreadN', str(self.threads),
                    '--limitBAMsortRAM', self.memory_limit,
                    '--outSAMunmapped', 'Within',
                    '--outSAMattributes', 'Standard',
                    '--quantMode', 'GeneCounts'
                ]
                
                # Handle compressed files
                if fastq_file.endswith('.gz'):
                    cmd.extend(['--readFilesCommand', 'zcat'])
                
                logger.info(f"Aligning sample {sample} ({i+1}/{len(sample_names)})")
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    bam_file = sample_output_dir / f"{sample}_Aligned.sortedByCoord.out.bam"
                    if bam_file.exists():
                        bam_files.append(str(bam_file))
                        logger.info(f"Successfully aligned {sample}")
                    else:
                        logger.error(f"BAM file not found for {sample}")
                else:
                    logger.error(f"STAR alignment failed for {sample}: {result.stderr}")
                    
            except Exception as e:
                logger.error(f"Error aligning {sample}: {e}")
                
        return bam_files
    
    def create_alignment_summary(self, bam_files: List[str]) -> Dict:
        """
        Create alignment summary statistics.
        
        Args:
            bam_files: List of BAM file paths
            
        Returns:
            Dictionary with alignment statistics
        """
        logger.info("Creating alignment summary")
        
        summary = {}
        
        for bam_file in bam_files:
            try:
                sample_name = Path(bam_file).parent.name
                
                # Use samtools to get alignment statistics
                cmd = ['samtools', 'flagstat', bam_file]
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    stats = self._parse_samtools_flagstat(result.stdout)
                    summary[sample_name] = stats
                    logger.info(f"Alignment stats for {sample_name}: {stats['mapped_percent']:.1f}% mapped")
                else:
                    logger.error(f"Failed to get stats for {bam_file}")
                    
            except Exception as e:
                logger.error(f"Error getting alignment stats for {bam_file}: {e}")
                
        return summary
    
    def _parse_samtools_flagstat(self, flagstat_output: str) -> Dict:
        """Parse samtools flagstat output."""
        stats = {}
        
        for line in flagstat_output.strip().split('\n'):
            if 'total' in line:
                stats['total_reads'] = int(line.split()[0])
            elif 'mapped' in line and 'mapped (' in line:
                stats['mapped_reads'] = int(line.split()[0])
                stats['mapped_percent'] = float(line.split('(')[1].split('%')[0])
            elif 'paired in sequencing' in line:
                stats['paired_reads'] = int(line.split()[0])
                
        return stats
    
    def generate_alignment_report(self, summary: Dict) -> str:
        """
        Generate alignment report.
        
        Args:
            summary: Alignment summary statistics
            
        Returns:
            Path to generated report
        """
        report_path = self.output_dir / "alignment_report.html"
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>RNA-seq Alignment Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>RNA-seq Alignment Report</h1>
                <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="section">
                <h2>Alignment Summary</h2>
                <table>
                    <tr>
                        <th>Sample</th>
                        <th>Total Reads</th>
                        <th>Mapped Reads</th>
                        <th>Mapping Rate (%)</th>
                    </tr>
        """
        
        for sample, stats in summary.items():
            html_content += f"""
                    <tr>
                        <td>{sample}</td>
                        <td>{stats.get('total_reads', 'N/A')}</td>
                        <td>{stats.get('mapped_reads', 'N/A')}</td>
                        <td>{stats.get('mapped_percent', 'N/A'):.1f}</td>
                    </tr>
            """
        
        html_content += """
                </table>
            </div>
        </body>
        </html>
        """
        
        with open(report_path, 'w') as f:
            f.write(html_content)
            
        logger.info(f"Alignment report generated: {report_path}")
        return str(report_path) 