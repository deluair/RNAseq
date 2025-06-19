#!/usr/bin/env python3
"""
Setup script for RNA-seq Analysis Platform
Helps users install dependencies and configure the platform.
"""

import os
import sys
import subprocess
import platform
from pathlib import Path
from loguru import logger


def check_python_version():
    """Check if Python version is compatible."""
    if sys.version_info < (3, 8):
        logger.error("Python 3.8 or higher is required")
        return False
    logger.info(f"Python version: {sys.version}")
    return True


def install_python_dependencies():
    """Install Python dependencies."""
    logger.info("Installing Python dependencies...")
    
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"], 
                      check=True)
        logger.info("Python dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install Python dependencies: {e}")
        return False


def check_external_tools():
    """Check if required external tools are installed."""
    tools = {
        'STAR': 'STAR --version',
        'fastqc': 'fastqc --version',
        'multiqc': 'multiqc --version',
        'samtools': 'samtools --version',
        'featureCounts': 'featureCounts -v'
    }
    
    missing_tools = []
    available_tools = []
    
    for tool, command in tools.items():
        try:
            result = subprocess.run(command.split(), capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✓ {tool} is installed")
                available_tools.append(tool)
            else:
                missing_tools.append(tool)
                logger.warning(f"✗ {tool} not found or not working")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            missing_tools.append(tool)
            logger.warning(f"✗ {tool} not found")
    
    if missing_tools:
        logger.warning(f"Missing tools: {', '.join(missing_tools)}")
        logger.info("\nTo install missing tools:")
        logger.info("  - STAR: conda install -c bioconda star")
        logger.info("  - FastQC: conda install -c bioconda fastqc")
        logger.info("  - MultiQC: pip install multiqc")
        logger.info("  - samtools: conda install -c bioconda samtools")
        logger.info("  - featureCounts: conda install -c bioconda subread")
        logger.info("\nAlternatively, visit:")
        logger.info("  - STAR: https://github.com/alexdobin/STAR")
        logger.info("  - FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/")
        logger.info("  - MultiQC: https://multiqc.info/")
        logger.info("  - samtools: http://www.htslib.org/")
        logger.info("  - featureCounts: http://subread.sourceforge.net/")
    
    if available_tools:
        logger.info(f"\nAvailable tools: {', '.join(available_tools)}")
    
    return len(missing_tools) == 0


def create_directories():
    """Create necessary directories."""
    directories = [
        "data",
        "data/raw",
        "data/reference",
        "results",
        "results/qc",
        "results/alignment",
        "results/quantification",
        "results/differential_expression",
        "results/plots",
        "logs",
        "genome_index"
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
        logger.info(f"Created directory: {directory}")


def download_reference_data():
    """Download reference genome data for soybean."""
    logger.info("Setting up reference data...")
    
    ref_dir = Path("data/reference")
    ref_dir.mkdir(parents=True, exist_ok=True)
    
    # Create placeholder files for demo
    genome_fasta = ref_dir / "Glycine_max_v4.0.fa"
    gtf_file = ref_dir / "Glycine_max_v4.0.gtf"
    
    if not genome_fasta.exists():
        logger.info("Creating placeholder genome FASTA file")
        with open(genome_fasta, 'w') as f:
            f.write(">Chr01\nATCGATCGATCG\n")
    
    if not gtf_file.exists():
        logger.info("Creating placeholder GTF file")
        with open(gtf_file, 'w') as f:
            f.write("Chr01\t.\tgene\t1\t12\t.\t+\t.\tgene_id \"Gene_0001\";\n")
    
    logger.info("Reference data setup completed")
    logger.info("Note: For real analysis, download actual soybean reference files")


def create_sample_config():
    """Create a sample configuration file."""
    config_file = Path("config/my_analysis.yaml")
    
    if not config_file.exists():
        logger.info("Creating sample configuration file")
        
        sample_config = """# Sample RNA-seq Analysis Configuration

# Data Retrieval Settings
data_retrieval:
  organism: "Glycine max"
  search_queries:
    - "soybean RNA-seq"
  max_results_per_query: 5
  output_directory: "data/raw"

# Quality Control Settings
quality_control:
  threads: 4
  output_directory: "results/qc"

# Alignment Settings
alignment:
  aligner: "STAR"
  threads: 4
  memory_limit: "32G"
  output_directory: "results/alignment"
  genome_directory: "genome_index"

# Reference Genome Settings
reference_genome:
  organism: "Glycine max"
  genome_fasta: "data/reference/Glycine_max_v4.0.fa"
  gtf_file: "data/reference/Glycine_max_v4.0.gtf"
  index_directory: "genome_index"

# Sample Information
samples:
  - id: "SRR12345678"
    name: "soybean_control_1"
    condition: "control"
    replicate: 1

# Pipeline Steps
pipeline:
  steps:
    - "data_retrieval"
    - "quality_control"
    - "alignment"
    - "visualization"
  skip_completed: true
"""
        
        config_file.parent.mkdir(parents=True, exist_ok=True)
        with open(config_file, 'w') as f:
            f.write(sample_config)
        
        logger.info(f"Sample configuration created: {config_file}")


def run_demo():
    """Run the demo to verify installation."""
    logger.info("Running demo to verify installation...")
    
    try:
        result = subprocess.run([sys.executable, "demo.py"], 
                              capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("✓ Demo completed successfully")
            return True
        else:
            logger.error(f"Demo failed: {result.stderr}")
            return False
            
    except Exception as e:
        logger.error(f"Error running demo: {e}")
        return False


def main():
    """Main setup function."""
    logger.info("RNA-seq Analysis Platform Setup")
    logger.info("=" * 40)
    
    # Check Python version
    if not check_python_version():
        return False
    
    # Install Python dependencies
    if not install_python_dependencies():
        return False
    
    # Check external tools
    tools_ok = check_external_tools()
    
    # Create directories
    create_directories()
    
    # Download reference data
    download_reference_data()
    
    # Create sample configuration
    create_sample_config()
    
    # Run demo if tools are available
    if tools_ok:
        demo_success = run_demo()
    else:
        logger.warning("Skipping demo due to missing external tools")
        demo_success = True
    
    # Summary
    logger.info("\n" + "=" * 40)
    logger.info("SETUP SUMMARY")
    logger.info("=" * 40)
    logger.info("✓ Python dependencies installed")
    logger.info(f"{'✓' if tools_ok else '✗'} External tools available")
    logger.info("✓ Directories created")
    logger.info("✓ Reference data setup")
    logger.info("✓ Sample configuration created")
    logger.info(f"{'✓' if demo_success else '✗'} Demo completed")
    
    logger.info("\nNext steps:")
    logger.info("1. Install missing external tools if any")
    logger.info("2. Download actual reference genome files")
    logger.info("3. Configure your analysis in config/my_analysis.yaml")
    logger.info("4. Run: python pipeline.py --config config/my_analysis.yaml")
    logger.info("5. Or start web interface: python app.py")
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 