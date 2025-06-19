"""
RNA Splicing Analysis Module

This module provides comprehensive RNA splicing analysis capabilities including:
- Splice junction detection and quantification
- Alternative splicing event identification (SE, MXE, A3SS, A5SS, RI)
- Differential splicing analysis
- Splicing efficiency calculations
- Visualization of splicing patterns

Author: RNA-seq Analysis Platform
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from collections import defaultdict, Counter
import warnings
from typing import Dict, List, Tuple, Optional
import json
from loguru import logger

class SplicingAnalyzer:
    """
    Comprehensive RNA splicing analysis class
    """
    
    def __init__(self, output_dir: str = "results/splicing"):
        """
        Initialize splicing analyzer
        
        Args:
            output_dir: Directory to save results
        """
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(f"{self.output_dir}/plots", exist_ok=True)
        
        # Alternative splicing event types
        self.as_event_types = {
            'SE': 'Skipped Exon',
            'MXE': 'Mutually Exclusive Exons', 
            'A3SS': 'Alternative 3\' Splice Site',
            'A5SS': 'Alternative 5\' Splice Site',
            'RI': 'Retained Intron'
        }
        
        logger.info(f"SplicingAnalyzer initialized with output directory: {output_dir}")
    
    def detect_splice_junctions(self, bam_files: List[str], gtf_file: str = None) -> pd.DataFrame:
        """
        Detect splice junctions from BAM files
        
        Args:
            bam_files: List of BAM file paths
            gtf_file: GTF annotation file (optional)
            
        Returns:
            DataFrame with splice junction information
        """
        logger.info("Detecting splice junctions from alignment files")
        
        # Mock splice junction data for demonstration
        # In real implementation, would parse BAM files using pysam
        junctions = []
        
        for i, bam_file in enumerate(bam_files):
            sample_name = os.path.basename(bam_file).replace('.bam', '')
            
            # Generate mock junction data
            n_junctions = np.random.randint(500, 1500)
            
            for j in range(n_junctions):
                # Mock chromosomal coordinates
                chr_name = f"chr{np.random.randint(1, 23)}"
                start = np.random.randint(1000000, 50000000)
                end = start + np.random.randint(100, 10000)
                
                # Mock read counts
                unique_reads = np.random.poisson(10)
                multi_reads = np.random.poisson(2)
                
                # Junction type
                junction_type = np.random.choice(['known', 'novel'], p=[0.8, 0.2])
                
                junctions.append({
                    'sample': sample_name,
                    'chromosome': chr_name,
                    'start': start,
                    'end': end,
                    'strand': np.random.choice(['+', '-']),
                    'unique_reads': unique_reads,
                    'multi_reads': multi_reads,
                    'total_reads': unique_reads + multi_reads,
                    'junction_type': junction_type,
                    'junction_id': f"{chr_name}:{start}-{end}"
                })
        
        df_junctions = pd.DataFrame(junctions)
        
        # Save results
        output_file = f"{self.output_dir}/splice_junctions.csv"
        df_junctions.to_csv(output_file, index=False)
        logger.info(f"Splice junctions saved to {output_file}")
        
        return df_junctions
    
    def identify_alternative_splicing_events(self, junction_data: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Identify alternative splicing events from junction data
        
        Args:
            junction_data: DataFrame with splice junction information
            
        Returns:
            Dictionary of DataFrames for each AS event type
        """
        logger.info("Identifying alternative splicing events")
        
        as_events = {}
        
        for event_type in self.as_event_types.keys():
            events = self._generate_mock_as_events(junction_data, event_type)
            as_events[event_type] = events
            
            # Save individual event files
            output_file = f"{self.output_dir}/as_events_{event_type}.csv"
            events.to_csv(output_file, index=False)
            logger.info(f"{event_type} events saved to {output_file}")
        
        return as_events
    
    def _generate_mock_as_events(self, junction_data: pd.DataFrame, event_type: str) -> pd.DataFrame:
        """Generate mock alternative splicing events for demonstration"""
        
        samples = junction_data['sample'].unique()
        n_events = np.random.randint(50, 200)
        
        events = []
        for i in range(n_events):
            # Mock event coordinates
            chr_name = f"chr{np.random.randint(1, 23)}"
            start = np.random.randint(1000000, 50000000)
            
            # Event-specific coordinates
            if event_type == 'SE':  # Skipped Exon
                exon_start = start + 1000
                exon_end = exon_start + np.random.randint(50, 500)
                upstream_end = exon_start - 50
                downstream_start = exon_end + 50
                
                event_coords = {
                    'upstream_exon_end': upstream_end,
                    'exon_start': exon_start,
                    'exon_end': exon_end,
                    'downstream_exon_start': downstream_start
                }
            
            elif event_type == 'A3SS':  # Alternative 3' SS
                event_coords = {
                    'long_exon_start': start,
                    'short_exon_start': start + np.random.randint(50, 200),
                    'common_end': start + 1000
                }
            
            elif event_type == 'A5SS':  # Alternative 5' SS  
                event_coords = {
                    'common_start': start,
                    'short_exon_end': start + 500,
                    'long_exon_end': start + 500 + np.random.randint(50, 200)
                }
            
            elif event_type == 'MXE':  # Mutually Exclusive Exons
                event_coords = {
                    'upstream_end': start,
                    'exon1_start': start + 100,
                    'exon1_end': start + 300,
                    'exon2_start': start + 400,
                    'exon2_end': start + 600,
                    'downstream_start': start + 700
                }
            
            else:  # RI - Retained Intron
                event_coords = {
                    'exon1_end': start,
                    'intron_start': start + 1,
                    'intron_end': start + np.random.randint(100, 2000),
                    'exon2_start': start + np.random.randint(100, 2000) + 1
                }
            
            # Generate PSI values for each sample
            base_event = {
                'event_id': f"{event_type}_{chr_name}_{start}_{i}",
                'event_type': event_type,
                'chromosome': chr_name,
                'strand': np.random.choice(['+', '-']),
                'gene_name': f"GENE_{i:04d}",
                **event_coords
            }
            
            for sample in samples:
                # PSI (Percent Spliced In) values
                psi_value = np.random.beta(2, 2)  # Values between 0 and 1
                inclusion_reads = np.random.poisson(20)
                exclusion_reads = np.random.poisson(15)
                
                event = base_event.copy()
                event.update({
                    'sample': sample,
                    'psi': psi_value,
                    'inclusion_reads': inclusion_reads,
                    'exclusion_reads': exclusion_reads,
                    'total_reads': inclusion_reads + exclusion_reads,
                    'confidence': 'high' if inclusion_reads + exclusion_reads > 10 else 'low'
                })
                
                events.append(event)
        
        return pd.DataFrame(events)
    
    def calculate_differential_splicing(self, as_events: Dict[str, pd.DataFrame], 
                                      condition1_samples: List[str],
                                      condition2_samples: List[str]) -> Dict[str, pd.DataFrame]:
        """
        Calculate differential splicing between conditions
        
        Args:
            as_events: Dictionary of alternative splicing events
            condition1_samples: Sample names for condition 1
            condition2_samples: Sample names for condition 2
            
        Returns:
            Dictionary of differential splicing results
        """
        logger.info("Calculating differential splicing between conditions")
        
        diff_splicing = {}
        
        for event_type, events_df in as_events.items():
            logger.info(f"Analyzing differential splicing for {event_type} events")
            
            # Group events by event_id
            diff_events = []
            
            for event_id in events_df['event_id'].unique():
                event_data = events_df[events_df['event_id'] == event_id]
                
                # Get PSI values for each condition
                cond1_psi = event_data[event_data['sample'].isin(condition1_samples)]['psi'].values
                cond2_psi = event_data[event_data['sample'].isin(condition2_samples)]['psi'].values
                
                if len(cond1_psi) > 0 and len(cond2_psi) > 0:
                    # Calculate statistics
                    mean_psi_cond1 = np.mean(cond1_psi)
                    mean_psi_cond2 = np.mean(cond2_psi)
                    delta_psi = mean_psi_cond2 - mean_psi_cond1
                    
                    # Statistical test
                    if len(cond1_psi) > 1 and len(cond2_psi) > 1:
                        statistic, p_value = stats.ttest_ind(cond1_psi, cond2_psi)
                    else:
                        statistic, p_value = 0, 1.0
                    
                    # Get event metadata
                    event_info = event_data.iloc[0]
                    
                    diff_events.append({
                        'event_id': event_id,
                        'event_type': event_type,
                        'chromosome': event_info['chromosome'],
                        'strand': event_info['strand'],
                        'gene_name': event_info['gene_name'],
                        'mean_psi_cond1': mean_psi_cond1,
                        'mean_psi_cond2': mean_psi_cond2,
                        'delta_psi': delta_psi,
                        'p_value': p_value,
                        'statistic': statistic,
                        'n_samples_cond1': len(cond1_psi),
                        'n_samples_cond2': len(cond2_psi)
                    })
            
            if diff_events:
                diff_df = pd.DataFrame(diff_events)
                
                # Multiple testing correction
                from scipy.stats import false_discovery_control
                diff_df['p_adj'] = false_discovery_control(diff_df['p_value'])
                
                # Add significance labels
                diff_df['significant'] = (diff_df['p_adj'] < 0.05) & (abs(diff_df['delta_psi']) > 0.1)
                
                # Sort by significance
                diff_df = diff_df.sort_values('p_adj')
                
                diff_splicing[event_type] = diff_df
                
                # Save results
                output_file = f"{self.output_dir}/differential_splicing_{event_type}.csv"
                diff_df.to_csv(output_file, index=False)
                logger.info(f"Differential splicing results for {event_type} saved to {output_file}")
        
        return diff_splicing
    
    def calculate_splicing_efficiency(self, junction_data: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate splicing efficiency metrics
        
        Args:
            junction_data: DataFrame with splice junction information
            
        Returns:
            DataFrame with splicing efficiency metrics
        """
        logger.info("Calculating splicing efficiency metrics")
        
        efficiency_data = []
        
        # Group by sample
        for sample in junction_data['sample'].unique():
            sample_data = junction_data[junction_data['sample'] == sample]
            
            # Calculate various efficiency metrics
            total_junctions = len(sample_data)
            total_reads = sample_data['total_reads'].sum()
            novel_junctions = len(sample_data[sample_data['junction_type'] == 'novel'])
            known_junctions = len(sample_data[sample_data['junction_type'] == 'known'])
            
            # Splicing efficiency (reads per junction)
            splicing_efficiency = total_reads / total_junctions if total_junctions > 0 else 0
            
            # Novel junction rate
            novel_rate = novel_junctions / total_junctions if total_junctions > 0 else 0
            
            efficiency_data.append({
                'sample': sample,
                'total_junctions': total_junctions,
                'total_reads': total_reads,
                'known_junctions': known_junctions,
                'novel_junctions': novel_junctions,
                'splicing_efficiency': splicing_efficiency,
                'novel_junction_rate': novel_rate,
                'mean_reads_per_junction': sample_data['total_reads'].mean(),
                'median_reads_per_junction': sample_data['total_reads'].median()
            })
        
        efficiency_df = pd.DataFrame(efficiency_data)
        
        # Save results
        output_file = f"{self.output_dir}/splicing_efficiency.csv"
        efficiency_df.to_csv(output_file, index=False)
        logger.info(f"Splicing efficiency results saved to {output_file}")
        
        return efficiency_df
    
    def plot_splicing_overview(self, junction_data: pd.DataFrame, as_events: Dict[str, pd.DataFrame]) -> str:
        """
        Create comprehensive splicing overview plots
        
        Args:
            junction_data: Splice junction data
            as_events: Alternative splicing events data
            
        Returns:
            Path to saved plot file
        """
        logger.info("Creating splicing overview plots")
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('RNA Splicing Analysis Overview', fontsize=16, fontweight='bold')
        
        # 1. Junction types distribution
        ax1 = axes[0, 0]
        junction_counts = junction_data.groupby(['sample', 'junction_type']).size().unstack(fill_value=0)
        junction_counts.plot(kind='bar', ax=ax1, color=['skyblue', 'lightcoral'])
        ax1.set_title('Splice Junction Types by Sample')
        ax1.set_xlabel('Sample')
        ax1.set_ylabel('Number of Junctions')
        ax1.legend(title='Junction Type')
        plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45)
        
        # 2. Read distribution per junction
        ax2 = axes[0, 1]
        ax2.hist(junction_data['total_reads'], bins=50, alpha=0.7, color='lightgreen')
        ax2.set_title('Distribution of Reads per Junction')
        ax2.set_xlabel('Number of Reads')
        ax2.set_ylabel('Frequency')
        ax2.axvline(junction_data['total_reads'].median(), color='red', linestyle='--', 
                   label=f'Median: {junction_data["total_reads"].median():.1f}')
        ax2.legend()
        
        # 3. Alternative splicing events count
        ax3 = axes[0, 2]
        as_counts = {event_type: len(df['event_id'].unique()) for event_type, df in as_events.items()}
        bars = ax3.bar(as_counts.keys(), as_counts.values(), 
                      color=['#FF9999', '#66B2FF', '#99FF99', '#FFCC99', '#FF99CC'])
        ax3.set_title('Alternative Splicing Events Count')
        ax3.set_xlabel('Event Type')
        ax3.set_ylabel('Number of Events')
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom')
        
        # 4. PSI distribution by event type
        ax4 = axes[1, 0]
        psi_data = []
        event_types = []
        for event_type, df in as_events.items():
            psi_data.extend(df['psi'].values)
            event_types.extend([event_type] * len(df))
        
        psi_df = pd.DataFrame({'PSI': psi_data, 'Event_Type': event_types})
        sns.boxplot(data=psi_df, x='Event_Type', y='PSI', ax=ax4)
        ax4.set_title('PSI Distribution by Event Type')
        ax4.set_xlabel('Event Type')
        ax4.set_ylabel('PSI Value')
        plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45)
        
        # 5. Junction reads by chromosome
        ax5 = axes[1, 1]
        chr_reads = junction_data.groupby('chromosome')['total_reads'].sum().sort_values(ascending=False)
        top_chrs = chr_reads.head(10)
        ax5.bar(range(len(top_chrs)), top_chrs.values, color='lightblue')
        ax5.set_title('Total Junction Reads by Chromosome (Top 10)')
        ax5.set_xlabel('Chromosome')
        ax5.set_ylabel('Total Reads')
        ax5.set_xticks(range(len(top_chrs)))
        ax5.set_xticklabels(top_chrs.index, rotation=45)
        
        # 6. Sample correlation based on PSI values
        ax6 = axes[1, 2]
        # Create sample correlation matrix based on SE events
        if 'SE' in as_events:
            se_events = as_events['SE']
            psi_matrix = se_events.pivot(index='event_id', columns='sample', values='psi')
            psi_matrix = psi_matrix.fillna(0)
            
            if psi_matrix.shape[1] > 1:
                correlation_matrix = psi_matrix.corr()
                sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, ax=ax6)
                ax6.set_title('Sample Correlation (PSI values)')
            else:
                ax6.text(0.5, 0.5, 'Insufficient samples\nfor correlation', 
                        ha='center', va='center', transform=ax6.transAxes)
                ax6.set_title('Sample Correlation')
        
        plt.tight_layout()
        
        # Save plot
        plot_file = f"{self.output_dir}/plots/splicing_overview.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Splicing overview plot saved to {plot_file}")
        return plot_file
    
    def plot_differential_splicing(self, diff_splicing: Dict[str, pd.DataFrame]) -> str:
        """
        Create differential splicing volcano plots
        
        Args:
            diff_splicing: Differential splicing results
            
        Returns:
            Path to saved plot file
        """
        logger.info("Creating differential splicing plots")
        
        n_event_types = len(diff_splicing)
        if n_event_types == 0:
            logger.warning("No differential splicing data to plot")
            return None
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        fig.suptitle('Differential Splicing Analysis', fontsize=16, fontweight='bold')
        
        for i, (event_type, diff_df) in enumerate(diff_splicing.items()):
            if i >= 6:  # Maximum 6 subplots
                break
                
            ax = axes[i]
            
            # Volcano plot
            if not diff_df.empty:
                # Prepare data
                x = diff_df['delta_psi']
                y = -np.log10(diff_df['p_adj'].replace(0, 1e-300))  # Avoid log(0)
                
                # Color by significance
                colors = ['red' if sig else 'gray' for sig in diff_df['significant']]
                
                scatter = ax.scatter(x, y, c=colors, alpha=0.6, s=30)
                
                # Add significance thresholds
                ax.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p_adj = 0.05')
                ax.axvline(0.1, color='black', linestyle=':', alpha=0.5, label='ΔΨ = 0.1')
                ax.axvline(-0.1, color='black', linestyle=':', alpha=0.5)
                
                ax.set_xlabel('ΔΨ (Condition2 - Condition1)')
                ax.set_ylabel('-log10(p_adj)')
                ax.set_title(f'{self.as_event_types[event_type]} Events\n'
                           f'({len(diff_df[diff_df["significant"]])} significant)')
                
                # Add text for significant events count
                n_sig = len(diff_df[diff_df['significant']])
                ax.text(0.05, 0.95, f'Significant: {n_sig}', transform=ax.transAxes, 
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
            else:
                ax.text(0.5, 0.5, f'No {event_type} events', ha='center', va='center', 
                       transform=ax.transAxes)
                ax.set_title(f'{self.as_event_types[event_type]} Events')
        
        # Hide unused subplots
        for i in range(len(diff_splicing), 6):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        # Save plot
        plot_file = f"{self.output_dir}/plots/differential_splicing.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Differential splicing plot saved to {plot_file}")
        return plot_file
    
    def generate_splicing_report(self, junction_data: pd.DataFrame, 
                               as_events: Dict[str, pd.DataFrame],
                               diff_splicing: Dict[str, pd.DataFrame] = None,
                               efficiency_data: pd.DataFrame = None) -> str:
        """
        Generate comprehensive splicing analysis report
        
        Args:
            junction_data: Splice junction data
            as_events: Alternative splicing events
            diff_splicing: Differential splicing results (optional)
            efficiency_data: Splicing efficiency data (optional)
            
        Returns:
            Path to HTML report
        """
        logger.info("Generating splicing analysis report")
        
        # Calculate summary statistics
        n_samples = junction_data['sample'].nunique()
        total_junctions = len(junction_data)
        novel_junctions = len(junction_data[junction_data['junction_type'] == 'novel'])
        
        as_summary = {}
        for event_type, df in as_events.items():
            as_summary[event_type] = {
                'total_events': len(df['event_id'].unique()),
                'total_observations': len(df),
                'mean_psi': df['psi'].mean(),
                'high_confidence': len(df[df['confidence'] == 'high'])
            }
        
        # HTML report template
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>RNA Splicing Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                h1, h2, h3 {{ color: #2c3e50; }}
                .summary-box {{ background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin: 20px 0; }}
                table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
                th {{ background-color: #3498db; color: white; }}
                .metric {{ display: inline-block; margin: 10px; padding: 15px; background-color: #ecf0f1; border-radius: 5px; }}
                .significant {{ color: #e74c3c; font-weight: bold; }}
            </style>
        </head>
        <body>
            <h1>RNA Splicing Analysis Report</h1>
            
            <div class="summary-box">
                <h2>Analysis Summary</h2>
                <div class="metric">
                    <strong>Samples Analyzed:</strong> {n_samples}
                </div>
                <div class="metric">
                    <strong>Total Splice Junctions:</strong> {total_junctions:,}
                </div>
                <div class="metric">
                    <strong>Novel Junctions:</strong> {novel_junctions:,} ({novel_junctions/total_junctions*100:.1f}%)
                </div>
            </div>
            
            <h2>Splice Junction Analysis</h2>
            <p>Detected {total_junctions:,} splice junctions across {n_samples} samples.</p>
            
            <h3>Junction Statistics by Sample</h3>
            <table>
                <tr>
                    <th>Sample</th>
                    <th>Total Junctions</th>
                    <th>Known Junctions</th>
                    <th>Novel Junctions</th>
                    <th>Total Reads</th>
                </tr>
        """
        
        # Add sample statistics
        for sample in junction_data['sample'].unique():
            sample_data = junction_data[junction_data['sample'] == sample]
            total_junc = len(sample_data)
            known_junc = len(sample_data[sample_data['junction_type'] == 'known'])
            novel_junc = len(sample_data[sample_data['junction_type'] == 'novel'])
            total_reads = sample_data['total_reads'].sum()
            
            html_content += f"""
                <tr>
                    <td>{sample}</td>
                    <td>{total_junc:,}</td>
                    <td>{known_junc:,}</td>
                    <td>{novel_junc:,}</td>
                    <td>{total_reads:,}</td>
                </tr>
            """
        
        html_content += """
            </table>
            
            <h2>Alternative Splicing Events</h2>
        """
        
        # Add AS events summary
        for event_type, summary in as_summary.items():
            event_name = self.as_event_types[event_type]
            html_content += f"""
            <h3>{event_name} ({event_type})</h3>
            <div class="summary-box">
                <div class="metric">
                    <strong>Unique Events:</strong> {summary['total_events']:,}
                </div>
                <div class="metric">
                    <strong>Total Observations:</strong> {summary['total_observations']:,}
                </div>
                <div class="metric">
                    <strong>Mean PSI:</strong> {summary['mean_psi']:.3f}
                </div>
                <div class="metric">
                    <strong>High Confidence:</strong> {summary['high_confidence']:,}
                </div>
            </div>
            """
        
        # Add differential splicing if available
        if diff_splicing:
            html_content += "<h2>Differential Splicing Analysis</h2>"
            
            for event_type, diff_df in diff_splicing.items():
                if not diff_df.empty:
                    n_significant = len(diff_df[diff_df['significant']])
                    event_name = self.as_event_types[event_type]
                    
                    html_content += f"""
                    <h3>{event_name} - Differential Events</h3>
                    <div class="summary-box">
                        <div class="metric">
                            <strong>Total Events Tested:</strong> {len(diff_df):,}
                        </div>
                        <div class="metric significant">
                            <strong>Significant Events:</strong> {n_significant:,}
                        </div>
                        <div class="metric">
                            <strong>Median ΔΨ:</strong> {diff_df['delta_psi'].median():.3f}
                        </div>
                    </div>
                    """
        
        # Add efficiency data if available
        if efficiency_data is not None:
            html_content += f"""
            <h2>Splicing Efficiency</h2>
            <table>
                <tr>
                    <th>Sample</th>
                    <th>Splicing Efficiency</th>
                    <th>Novel Junction Rate</th>
                    <th>Mean Reads/Junction</th>
                </tr>
            """
            
            for _, row in efficiency_data.iterrows():
                html_content += f"""
                <tr>
                    <td>{row['sample']}</td>
                    <td>{row['splicing_efficiency']:.2f}</td>
                    <td>{row['novel_junction_rate']:.1%}</td>
                    <td>{row['mean_reads_per_junction']:.1f}</td>
                </tr>
                """
            
            html_content += "</table>"
        
        html_content += """
            <h2>Files Generated</h2>
            <ul>
                <li>splice_junctions.csv - Splice junction data</li>
                <li>as_events_*.csv - Alternative splicing events by type</li>
                <li>splicing_efficiency.csv - Splicing efficiency metrics</li>
                <li>plots/splicing_overview.png - Overview visualization</li>
        """
        
        if diff_splicing:
            html_content += "<li>differential_splicing_*.csv - Differential splicing results</li>"
            html_content += "<li>plots/differential_splicing.png - Differential splicing plots</li>"
        
        html_content += """
            </ul>
            
            <p><em>Report generated by RNA-seq Analysis Platform</em></p>
        </body>
        </html>
        """
        
        # Save report
        report_file = f"{self.output_dir}/splicing_report.html"
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Splicing analysis report saved to {report_file}")
        return report_file
    
    def run_complete_analysis(self, bam_files: List[str], 
                            condition1_samples: List[str] = None,
                            condition2_samples: List[str] = None,
                            gtf_file: str = None) -> Dict:
        """
        Run complete splicing analysis pipeline
        
        Args:
            bam_files: List of BAM files
            condition1_samples: Samples for condition 1 (for differential analysis)
            condition2_samples: Samples for condition 2 (for differential analysis)
            gtf_file: GTF annotation file
            
        Returns:
            Dictionary with all analysis results
        """
        logger.info("Running complete splicing analysis pipeline")
        
        results = {}
        
        # 1. Detect splice junctions
        junction_data = self.detect_splice_junctions(bam_files, gtf_file)
        results['junctions'] = junction_data
        
        # 2. Identify alternative splicing events
        as_events = self.identify_alternative_splicing_events(junction_data)
        results['as_events'] = as_events
        
        # 3. Calculate splicing efficiency
        efficiency_data = self.calculate_splicing_efficiency(junction_data)
        results['efficiency'] = efficiency_data
        
        # 4. Differential splicing analysis (if conditions provided)
        diff_splicing = None
        if condition1_samples and condition2_samples:
            diff_splicing = self.calculate_differential_splicing(
                as_events, condition1_samples, condition2_samples
            )
            results['differential'] = diff_splicing
        
        # 5. Generate visualizations
        overview_plot = self.plot_splicing_overview(junction_data, as_events)
        results['overview_plot'] = overview_plot
        
        if diff_splicing:
            diff_plot = self.plot_differential_splicing(diff_splicing)
            results['differential_plot'] = diff_plot
        
        # 6. Generate comprehensive report
        report_file = self.generate_splicing_report(
            junction_data, as_events, diff_splicing, efficiency_data
        )
        results['report'] = report_file
        
        logger.info("Complete splicing analysis finished successfully")
        
        # Save analysis summary
        summary_file = f"{self.output_dir}/analysis_summary.json"
        summary = {
            'n_samples': junction_data['sample'].nunique(),
            'total_junctions': len(junction_data),
            'as_events_count': {et: len(df['event_id'].unique()) for et, df in as_events.items()},
            'differential_analysis': diff_splicing is not None,
            'files_generated': [
                'splice_junctions.csv',
                'splicing_efficiency.csv',
                'splicing_report.html',
                'plots/splicing_overview.png'
            ]
        }
        
        if diff_splicing:
            summary['significant_events'] = {
                et: len(df[df['significant']]) for et, df in diff_splicing.items()
            }
            summary['files_generated'].extend([
                'plots/differential_splicing.png'
            ] + [f'differential_splicing_{et}.csv' for et in diff_splicing.keys()])
        
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        results['summary'] = summary
        
        return results


def create_mock_bam_files(output_dir: str = "data/alignment", n_samples: int = 6) -> List[str]:
    """
    Create mock BAM file paths for demonstration
    
    Args:
        output_dir: Directory where BAM files would be located
        n_samples: Number of sample files to create
        
    Returns:
        List of BAM file paths
    """
    os.makedirs(output_dir, exist_ok=True)
    
    bam_files = []
    for i in range(n_samples):
        bam_file = f"{output_dir}/sample_{i+1}_Aligned.sortedByCoord.out.bam"
        bam_files.append(bam_file)
        
        # Create empty file for demonstration
        with open(bam_file, 'w') as f:
            f.write("# Mock BAM file for splicing analysis\n")
    
    return bam_files


# Example usage
if __name__ == "__main__":
    # Initialize analyzer
    analyzer = SplicingAnalyzer()
    
    # Create mock BAM files
    bam_files = create_mock_bam_files()
    
    # Run complete analysis
    results = analyzer.run_complete_analysis(
        bam_files=bam_files,
        condition1_samples=['sample_1', 'sample_2', 'sample_3'],
        condition2_samples=['sample_4', 'sample_5', 'sample_6']
    )
    
    print("Splicing analysis completed!")
    print(f"Results saved to: {analyzer.output_dir}")
    print(f"Report available at: {results['report']}") 