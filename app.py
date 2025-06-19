"""
RNA-seq Analysis Platform Web Application.
Provides a web interface for RNA-seq data analysis.
"""

import os
import json
import argparse
import pandas as pd
from pathlib import Path
from flask import Flask, render_template, request, jsonify, send_file, redirect, url_for
from flask_cors import CORS
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.utils
from loguru import logger

# Import our modules
from src.data_retrieval import NCBIDownloader
from src.preprocessing import QualityControl
from src.alignment import STARAligner
from src.visualization import PlotGenerator

app = Flask(__name__)
CORS(app)

# Configure logging
logger.add("logs/app.log", rotation="1 day", retention="7 days")

# Global variables to store analysis state
current_analysis = {
    'datasets': [],
    'qc_results': None,
    'alignment_results': None,
    'expression_data': None,
    'de_results': None
}

@app.route('/')
def index():
    """Main page of the RNA-seq analysis platform."""
    return render_template('index.html')

@app.route('/api/search_datasets', methods=['POST'])
def search_datasets():
    """Search for RNA-seq datasets on NCBI SRA."""
    try:
        data = request.get_json()
        query = data.get('query', 'soybean RNA-seq')
        max_results = data.get('max_results', 10)
        
        downloader = NCBIDownloader()
        datasets = downloader.search_sra(query, max_results)
        
        current_analysis['datasets'] = datasets
        
        return jsonify({
            'success': True,
            'datasets': datasets,
            'count': len(datasets)
        })
        
    except Exception as e:
        logger.error(f"Error searching datasets: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/download_datasets', methods=['POST'])
def download_datasets():
    """Download selected datasets from NCBI SRA with improved error handling."""
    try:
        data = request.get_json()
        run_ids = data.get('run_ids', [])
        quick_download = data.get('quick_download', True)  # Default to quick download
        
        if not run_ids:
            return jsonify({
                'success': False,
                'error': 'No run IDs provided'
            }), 400
        
        downloader = NCBIDownloader()
        
        # Use quick download for demonstration purposes
        if quick_download:
            downloaded_files = downloader.download_datasets_quick(run_ids, max_reads=5000)
            download_type = "Quick download (5,000 reads per sample)"
        else:
            downloaded_files = downloader.download_sra_data(run_ids)
            download_type = "Full download"
        
        # Store downloaded files in current analysis
        current_analysis['downloaded_files'] = downloaded_files
        current_analysis['datasets'] = run_ids
        
        return jsonify({
            'success': True,
            'downloaded_files': downloaded_files,
            'count': len(downloaded_files),
            'download_type': download_type,
            'run_ids': run_ids,
            'message': f'Successfully downloaded {len(downloaded_files)} files for {len(run_ids)} datasets'
        })
        
    except Exception as e:
        logger.error(f"Error downloading datasets: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/download_datasets_full', methods=['POST'])
def download_datasets_full():
    """Download full datasets (may take longer)."""
    try:
        data = request.get_json()
        run_ids = data.get('run_ids', [])
        
        if not run_ids:
            return jsonify({
                'success': False,
                'error': 'No run IDs provided'
            }), 400
        
        downloader = NCBIDownloader()
        downloaded_files = downloader.download_sra_data(run_ids)
        
        # Store downloaded files in current analysis
        current_analysis['downloaded_files'] = downloaded_files
        current_analysis['datasets'] = run_ids
        
        return jsonify({
            'success': True,
            'downloaded_files': downloaded_files,
            'count': len(downloaded_files),
            'download_type': "Full download",
            'run_ids': run_ids,
            'message': f'Successfully downloaded {len(downloaded_files)} files for {len(run_ids)} datasets'
        })
        
    except Exception as e:
        logger.error(f"Error downloading full datasets: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/run_quality_control', methods=['POST'])
def run_quality_control():
    """Run quality control analysis on downloaded data."""
    try:
        data = request.get_json()
        fastq_files = data.get('fastq_files', [])
        
        if not fastq_files:
            return jsonify({
                'success': False,
                'error': 'No FASTQ files provided'
            }), 400
        
        qc = QualityControl()
        qc_outputs = qc.run_fastqc(fastq_files)
        
        if qc_outputs:
            multiqc_report = qc.run_multiqc(qc_outputs)
            qc_df = qc.parse_fastqc_results(qc_outputs)
            plot_paths = qc.create_qc_summary_plots(qc_df)
            qc_report = qc.generate_qc_report(qc_df, plot_paths)
            
            current_analysis['qc_results'] = {
                'qc_outputs': qc_outputs,
                'multiqc_report': multiqc_report,
                'qc_data': qc_df.to_dict() if not qc_df.empty else {},
                'plot_paths': plot_paths,
                'qc_report': qc_report
            }
        
        return jsonify({
            'success': True,
            'qc_outputs': qc_outputs,
            'multiqc_report': multiqc_report if qc_outputs else None
        })
        
    except Exception as e:
        logger.error(f"Error running quality control: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/run_alignment', methods=['POST'])
def run_alignment():
    """Run RNA-seq alignment using STAR."""
    try:
        data = request.get_json()
        fastq_files = data.get('fastq_files', [])
        sample_names = data.get('sample_names', [])
        genome_dir = data.get('genome_dir', 'genome_index')
        
        if not fastq_files or not sample_names:
            return jsonify({
                'success': False,
                'error': 'FASTQ files and sample names required'
            }), 400
        
        aligner = STARAligner(genome_dir=genome_dir)
        
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
            return jsonify({
                'success': False,
                'error': 'Mismatch between number of files and samples'
            }), 400
        
        if bam_files:
            alignment_summary = aligner.create_alignment_summary(bam_files)
            alignment_report = aligner.generate_alignment_report(alignment_summary)
            
            current_analysis['alignment_results'] = {
                'bam_files': bam_files,
                'alignment_summary': alignment_summary,
                'alignment_report': alignment_report
            }
        
        return jsonify({
            'success': True,
            'bam_files': bam_files,
            'alignment_summary': alignment_summary if bam_files else {}
        })
        
    except Exception as e:
        logger.error(f"Error running alignment: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/create_visualizations', methods=['POST'])
def create_visualizations():
    """Create visualizations for the analysis results."""
    try:
        plot_generator = PlotGenerator()
        plot_paths = {}
        
        # Create QC plots if available
        if current_analysis.get('qc_results'):
            qc_data = pd.DataFrame(current_analysis['qc_results']['qc_data'])
            if not qc_data.empty:
                qc_plots = plot_generator.create_quality_control_plots(qc_data)
                plot_paths.update(qc_plots)
        
        # Create alignment plots if available
        if current_analysis.get('alignment_results'):
            alignment_stats = current_analysis['alignment_results']['alignment_summary']
            if alignment_stats:
                alignment_plots = plot_generator.create_alignment_plots(alignment_stats)
                plot_paths.update(alignment_plots)
        
        # Create expression plots if available
        if current_analysis.get('expression_data'):
            expression_data = pd.DataFrame(current_analysis['expression_data'])
            if not expression_data.empty:
                expression_plots = plot_generator.create_expression_plots(expression_data)
                plot_paths.update(expression_plots)
        
        # Create differential expression plots if available
        if current_analysis.get('de_results'):
            de_results = pd.DataFrame(current_analysis['de_results'])
            if not de_results.empty:
                de_plots = plot_generator.create_differential_expression_plots(de_results)
                plot_paths.update(de_plots)
        
        # Generate summary report
        if plot_paths:
            summary_report = plot_generator.create_summary_report(plot_paths)
        else:
            summary_report = None
        
        return jsonify({
            'success': True,
            'plot_paths': plot_paths,
            'summary_report': summary_report
        })
        
    except Exception as e:
        logger.error(f"Error creating visualizations: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/get_analysis_status')
def get_analysis_status():
    """Get current analysis status."""
    return jsonify({
        'datasets_count': len(current_analysis.get('datasets', [])),
        'qc_completed': current_analysis.get('qc_results') is not None,
        'alignment_completed': current_analysis.get('alignment_results') is not None,
        'expression_data_available': current_analysis.get('expression_data') is not None,
        'de_results_available': current_analysis.get('de_results') is not None
    })

@app.route('/api/get_soybean_datasets')
def get_soybean_datasets():
    """Get soybean-specific datasets."""
    try:
        downloader = NCBIDownloader()
        datasets = downloader.get_soybean_datasets(limit=20)
        
        return jsonify({
            'success': True,
            'datasets': datasets
        })
        
    except Exception as e:
        logger.error(f"Error getting soybean datasets: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/reports/<path:filename>')
def serve_report(filename):
    """Serve generated reports."""
    report_path = Path('results') / filename
    if report_path.exists():
        return send_file(str(report_path))
    else:
        return "Report not found", 404

@app.route('/plots/<path:filename>')
def serve_plot(filename):
    """Serve generated plots."""
    plot_path = Path('results/plots') / filename
    if plot_path.exists():
        return send_file(str(plot_path))
    else:
        return "Plot not found", 404

@app.route('/api/upload_expression_data', methods=['POST'])
def upload_expression_data():
    """Upload expression data file."""
    try:
        if 'file' not in request.files:
            return jsonify({
                'success': False,
                'error': 'No file uploaded'
            }), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({
                'success': False,
                'error': 'No file selected'
            }), 400
        
        # Read the file
        if file.filename.endswith('.csv'):
            df = pd.read_csv(file)
        elif file.filename.endswith('.tsv'):
            df = pd.read_csv(file, sep='\t')
        elif file.filename.endswith('.xlsx'):
            df = pd.read_excel(file)
        else:
            return jsonify({
                'success': False,
                'error': 'Unsupported file format'
            }), 400
        
        # Store the data
        current_analysis['expression_data'] = df.to_dict()
        
        return jsonify({
            'success': True,
            'rows': len(df),
            'columns': len(df.columns),
            'columns': list(df.columns)
        })
        
    except Exception as e:
        logger.error(f"Error uploading expression data: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/upload_de_results', methods=['POST'])
def upload_de_results():
    """Upload differential expression results."""
    try:
        if 'file' not in request.files:
            return jsonify({
                'success': False,
                'error': 'No file uploaded'
            }), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({
                'success': False,
                'error': 'No file selected'
            }), 400
        
        # Read the file
        if file.filename.endswith('.csv'):
            df = pd.read_csv(file)
        elif file.filename.endswith('.tsv'):
            df = pd.read_csv(file, sep='\t')
        elif file.filename.endswith('.xlsx'):
            df = pd.read_excel(file)
        else:
            return jsonify({
                'success': False,
                'error': 'Unsupported file format'
            }), 400
        
        # Store the data
        current_analysis['de_results'] = df.to_dict()
        
        return jsonify({
            'success': True,
            'rows': len(df),
            'columns': list(df.columns)
        })
        
    except Exception as e:
        logger.error(f"Error uploading DE results: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/run_differential_expression', methods=['POST'])
def run_differential_expression():
    """Run differential expression analysis using DESeq2 or Python alternatives."""
    try:
        data = request.get_json()
        use_mock_data = data.get('use_mock_data', True)
        method = data.get('method', 'auto')  # 'deseq2', 'python', or 'auto'
        control_condition = data.get('control_condition', 'control')
        treatment_condition = data.get('treatment_condition', 'treatment')
        
        # Import the analysis modules
        from src.analysis import DESeq2Analyzer, PythonDEAnalyzer
        
        if use_mock_data:
            # Create mock data for demonstration
            logger.info("Creating mock count data for DE analysis")
            deseq2_analyzer = DESeq2Analyzer()
            count_matrix, sample_metadata = deseq2_analyzer.create_mock_count_data(
                n_genes=2000, n_samples=6
            )
        else:
            # Try to use uploaded data
            count_file = data.get('count_file')
            metadata_file = data.get('metadata_file')
            
            if not count_file or not metadata_file:
                return jsonify({
                    'success': False,
                    'error': 'Count matrix and metadata files required when not using mock data'
                }), 400
            
            # Load user data (this would need file upload handling)
            # For now, return error if no mock data
            return jsonify({
                'success': False,
                'error': 'File upload not implemented yet. Please use mock data.'
            }), 400
        
        # Run differential expression analysis
        if method == 'deseq2' or method == 'auto':
            analyzer = DESeq2Analyzer()
            results = analyzer.analyze_differential_expression(
                count_matrix, sample_metadata,
                control_condition, treatment_condition,
                use_deseq2=True
            )
        else:
            analyzer = PythonDEAnalyzer()
            results = analyzer.analyze_differential_expression(
                count_matrix, sample_metadata,
                control_condition, treatment_condition,
                method='edger_like'
            )
        
        # Create plots
        if method == 'deseq2' or method == 'auto':
            plot_paths = analyzer.create_de_plots(results['all_results'])
        else:
            # For Python analyzer, we need to add plotting method
            plot_paths = {}
        
        # Store results in current analysis
        current_analysis['de_results'] = results['all_results'].to_dict()
        current_analysis['de_significant'] = results['significant_genes'].to_dict()
        current_analysis['de_summary'] = results['summary_stats']
        current_analysis['de_plots'] = plot_paths
        
        return jsonify({
            'success': True,
            'method_used': results['method_used'],
            'summary_stats': results['summary_stats'],
            'plot_paths': plot_paths,
            'significant_genes_count': len(results['significant_genes']),
            'total_genes_tested': len(results['all_results'])
        })
        
    except Exception as e:
        logger.error(f"Error running differential expression analysis: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/get_de_results')
def get_de_results():
    """Get differential expression analysis results."""
    try:
        if 'de_results' not in current_analysis:
            return jsonify({
                'success': False,
                'error': 'No DE analysis results available'
            }), 404
        
        # Convert to DataFrame for easier processing
        import pandas as pd
        de_results = pd.DataFrame.from_dict(current_analysis['de_results'], orient='index')
        
        # Get top significant genes
        if 'padj' in de_results.columns and 'log2FoldChange' in de_results.columns:
            significant = de_results[
                (de_results['padj'] < 0.05) & 
                (de_results['log2FoldChange'].abs() > 1) &
                (~de_results['padj'].isna())
            ]
            top_up = significant[significant['log2FoldChange'] > 0].nlargest(10, 'log2FoldChange')
            top_down = significant[significant['log2FoldChange'] < 0].nsmallest(10, 'log2FoldChange')
        else:
            top_up = pd.DataFrame()
            top_down = pd.DataFrame()
        
        return jsonify({
            'success': True,
            'summary': current_analysis.get('de_summary', {}),
            'top_upregulated': top_up.to_dict('records') if not top_up.empty else [],
            'top_downregulated': top_down.to_dict('records') if not top_down.empty else [],
            'plot_paths': current_analysis.get('de_plots', {})
        })
        
    except Exception as e:
        logger.error(f"Error getting DE results: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/compare_de_methods', methods=['POST'])
def compare_de_methods():
    """Compare different DE analysis methods."""
    try:
        from src.analysis import DESeq2Analyzer, PythonDEAnalyzer
        
        # Create mock data
        deseq2_analyzer = DESeq2Analyzer()
        count_matrix, sample_metadata = deseq2_analyzer.create_mock_count_data(
            n_genes=1000, n_samples=6
        )
        
        # Run DESeq2 analysis
        try:
            deseq2_results = deseq2_analyzer.analyze_differential_expression(
                count_matrix, sample_metadata, use_deseq2=True
            )
            deseq2_available = True
        except Exception as e:
            logger.warning(f"DESeq2 analysis failed: {e}")
            deseq2_results = deseq2_analyzer.analyze_differential_expression(
                count_matrix, sample_metadata, use_deseq2=False
            )
            deseq2_available = False
        
        # Run Python analysis comparison
        python_analyzer = PythonDEAnalyzer()
        python_comparison = python_analyzer.compare_methods(
            count_matrix, sample_metadata
        )
        
        # Combine results
        all_results = {
            'deseq2_results': deseq2_results,
            'python_comparison': python_comparison,
            'deseq2_available': deseq2_available
        }
        
        return jsonify({
            'success': True,
            'comparison_results': {
                'deseq2_method': deseq2_results['method_used'],
                'deseq2_summary': deseq2_results['summary_stats'],
                'python_methods': python_comparison['comparison_summary'].to_dict('records')
            },
            'deseq2_available': deseq2_available
        })
        
    except Exception as e:
        logger.error(f"Error comparing DE methods: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/check_deseq2_status')
def check_deseq2_status():
    """Check if R and DESeq2 are available."""
    try:
        from src.analysis import DESeq2Analyzer
        
        analyzer = DESeq2Analyzer()
        
        return jsonify({
            'r_available': analyzer.r_available,
            'deseq2_available': analyzer.deseq2_available,
            'status': 'ready' if analyzer.deseq2_available else 'python_only'
        })
        
    except Exception as e:
        logger.error(f"Error checking DESeq2 status: {e}")
        return jsonify({
            'r_available': False,
            'deseq2_available': False,
            'status': 'error',
            'error': str(e)
        }), 500

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RNA-seq Analysis Platform')
    parser.add_argument('--port', type=int, default=5000, help='Port to run the web server on')
    args = parser.parse_args()

    # Create necessary directories
    os.makedirs('logs', exist_ok=True)
    os.makedirs('templates', exist_ok=True)
    os.makedirs('static', exist_ok=True)
    
    # Run the application
    app.run(debug=True, host='0.0.0.0', port=args.port) 