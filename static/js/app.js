// RNA-seq Analysis Platform JavaScript

class RNAseqPlatform {
    constructor() {
        this.selectedDatasets = [];
        this.analysisStatus = {
            data: false,
            qc: false,
            alignment: false,
            visualization: false
        };
        this.init();
    }

    init() {
        this.bindEvents();
        this.updateStatus();
        this.loadSoybeanDatasets();
    }

    bindEvents() {
        // Search form
        document.getElementById('searchForm').addEventListener('submit', (e) => {
            e.preventDefault();
            this.searchDatasets();
        });

        // Download button
        document.getElementById('downloadBtn').addEventListener('click', () => {
            this.downloadDatasets();
        });

        // QC form
        document.getElementById('qcForm').addEventListener('submit', (e) => {
            e.preventDefault();
            this.runQualityControl();
        });

        // Alignment form
        document.getElementById('alignmentForm').addEventListener('submit', (e) => {
            e.preventDefault();
            this.runAlignment();
        });

        // Upload form
        document.getElementById('uploadForm').addEventListener('submit', (e) => {
            e.preventDefault();
            this.uploadFiles();
        });

        // Generate plots button
        document.getElementById('generatePlotsBtn').addEventListener('click', () => {
            this.generateVisualizations();
        });

        // Tab change events
        document.querySelectorAll('[data-bs-toggle="tab"]').forEach(tab => {
            tab.addEventListener('shown.bs.tab', (e) => {
                this.onTabChange(e.target.getAttribute('data-bs-target'));
            });
        });
    }

    async searchDatasets() {
        const query = document.getElementById('searchQuery').value;
        const maxResults = document.getElementById('maxResults').value;

        this.showLoading('Searching NCBI SRA...');

        try {
            const response = await fetch('/api/search_datasets', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    query: query,
                    max_results: parseInt(maxResults)
                })
            });

            const data = await response.json();

            if (data.success) {
                this.displayDatasets(data.datasets);
                this.showAlert(`Found ${data.count} datasets`, 'success');
            } else {
                this.showAlert(data.error, 'danger');
            }
        } catch (error) {
            this.showAlert('Error searching datasets: ' + error.message, 'danger');
        } finally {
            this.hideLoading();
        }
    }

    async loadSoybeanDatasets() {
        try {
            const response = await fetch('/api/get_soybean_datasets');
            const data = await response.json();

            if (data.success) {
                this.displayDatasets(data.datasets);
            }
        } catch (error) {
            console.error('Error loading soybean datasets:', error);
        }
    }

    displayDatasets(datasets) {
        const container = document.getElementById('datasetsList');
        
        if (datasets.length === 0) {
            container.innerHTML = '<div class="text-center text-muted">No datasets found</div>';
            return;
        }

        container.innerHTML = datasets.map(dataset => `
            <div class="dataset-card" data-run-id="${dataset.run_id}">
                <div class="form-check">
                    <input class="form-check-input dataset-checkbox" type="checkbox" 
                           value="${dataset.run_id}" id="dataset_${dataset.run_id}">
                    <label class="form-check-label" for="dataset_${dataset.run_id}">
                        <h6>${dataset.title}</h6>
                        <div class="mb-2">
                            <span class="badge bg-primary">${dataset.organism}</span>
                            <span class="badge bg-secondary">${dataset.layout}</span>
                        </div>
                        <small class="text-muted">Run ID: ${dataset.run_id}</small>
                    </label>
                </div>
            </div>
        `).join('');

        // Bind checkbox events
        document.querySelectorAll('.dataset-checkbox').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                this.toggleDatasetSelection(e.target.value, e.target.checked);
            });
        });
    }

    toggleDatasetSelection(runId, selected) {
        if (selected) {
            this.selectedDatasets.push(runId);
        } else {
            this.selectedDatasets = this.selectedDatasets.filter(id => id !== runId);
        }

        const downloadBtn = document.getElementById('downloadBtn');
        downloadBtn.disabled = this.selectedDatasets.length === 0;
    }

    async downloadDatasets() {
        if (this.selectedDatasets.length === 0) {
            this.showAlert('Please select datasets to download', 'warning');
            return;
        }

        this.showLoading('Downloading datasets...');
        this.showProgress();

        try {
            const response = await fetch('/api/download_datasets', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    run_ids: this.selectedDatasets
                })
            });

            const data = await response.json();

            if (data.success) {
                this.showAlert(`Successfully downloaded ${data.count} files`, 'success');
                this.analysisStatus.data = true;
                this.updateStatus();
                this.updatePipelineStep(1, true);
            } else {
                this.showAlert(data.error, 'danger');
            }
        } catch (error) {
            this.showAlert('Error downloading datasets: ' + error.message, 'danger');
        } finally {
            this.hideLoading();
            this.hideProgress();
        }
    }

    async runQualityControl() {
        const fastqFiles = Array.from(document.getElementById('fastqFiles').selectedOptions)
            .map(option => option.value);
        const threads = document.getElementById('qcThreads').value;

        if (fastqFiles.length === 0) {
            this.showAlert('Please select FASTQ files for quality control', 'warning');
            return;
        }

        this.showLoading('Running quality control...');

        try {
            const response = await fetch('/api/run_quality_control', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    fastq_files: fastqFiles,
                    threads: parseInt(threads)
                })
            });

            const data = await response.json();

            if (data.success) {
                this.showAlert('Quality control completed successfully', 'success');
                this.analysisStatus.qc = true;
                this.updateStatus();
                this.updatePipelineStep(2, true);
                this.displayQCResults(data);
            } else {
                this.showAlert(data.error, 'danger');
            }
        } catch (error) {
            this.showAlert('Error running quality control: ' + error.message, 'danger');
        } finally {
            this.hideLoading();
        }
    }

    async runAlignment() {
        const alignmentFiles = Array.from(document.getElementById('alignmentFiles').selectedOptions)
            .map(option => option.value);
        const sampleNames = document.getElementById('sampleNames').value.split('\n').filter(name => name.trim());
        const genomeDir = document.getElementById('genomeDir').value;

        if (alignmentFiles.length === 0 || sampleNames.length === 0) {
            this.showAlert('Please provide alignment files and sample names', 'warning');
            return;
        }

        this.showLoading('Running alignment...');

        try {
            const response = await fetch('/api/run_alignment', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    fastq_files: alignmentFiles,
                    sample_names: sampleNames,
                    genome_dir: genomeDir
                })
            });

            const data = await response.json();

            if (data.success) {
                this.showAlert('Alignment completed successfully', 'success');
                this.analysisStatus.alignment = true;
                this.updateStatus();
                this.updatePipelineStep(3, true);
                this.displayAlignmentResults(data);
            } else {
                this.showAlert(data.error, 'danger');
            }
        } catch (error) {
            this.showAlert('Error running alignment: ' + error.message, 'danger');
        } finally {
            this.hideLoading();
        }
    }

    async uploadFiles() {
        const expressionFile = document.getElementById('expressionFile').files[0];
        const deFile = document.getElementById('deFile').files[0];

        if (!expressionFile && !deFile) {
            this.showAlert('Please select at least one file to upload', 'warning');
            return;
        }

        this.showLoading('Uploading files...');

        try {
            const formData = new FormData();
            if (expressionFile) formData.append('file', expressionFile);
            if (deFile) formData.append('file', deFile);

            const response = await fetch('/api/upload_expression_data', {
                method: 'POST',
                body: formData
            });

            const data = await response.json();

            if (data.success) {
                this.showAlert('Files uploaded successfully', 'success');
                this.analysisStatus.visualization = true;
                this.updateStatus();
            } else {
                this.showAlert(data.error, 'danger');
            }
        } catch (error) {
            this.showAlert('Error uploading files: ' + error.message, 'danger');
        } finally {
            this.hideLoading();
        }
    }

    async generateVisualizations() {
        this.showLoading('Generating visualizations...');

        try {
            const response = await fetch('/api/create_visualizations', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                }
            });

            const data = await response.json();

            if (data.success) {
                this.showAlert('Visualizations generated successfully', 'success');
                this.updatePipelineStep(4, true);
                this.displayVisualizationResults(data);
            } else {
                this.showAlert(data.error, 'danger');
            }
        } catch (error) {
            this.showAlert('Error generating visualizations: ' + error.message, 'danger');
        } finally {
            this.hideLoading();
        }
    }

    displayQCResults(data) {
        const container = document.getElementById('qcResults');
        container.innerHTML = `
            <div class="alert alert-success">
                <i class="fas fa-check-circle me-2"></i>
                Quality control completed successfully
            </div>
            <div class="result-item">
                <h6>QC Outputs</h6>
                <p>Generated ${data.qc_outputs.length} QC reports</p>
                ${data.multiqc_report ? `<a href="${data.multiqc_report}" target="_blank" class="btn btn-sm btn-primary">View MultiQC Report</a>` : ''}
            </div>
        `;
    }

    displayAlignmentResults(data) {
        const container = document.getElementById('alignmentResults');
        const summary = data.alignment_summary;
        
        let summaryHtml = '';
        if (summary && Object.keys(summary).length > 0) {
            summaryHtml = Object.entries(summary).map(([sample, stats]) => `
                <div class="mb-2">
                    <strong>${sample}:</strong> ${stats.mapped_percent?.toFixed(1) || 'N/A'}% mapped
                </div>
            `).join('');
        }

        container.innerHTML = `
            <div class="alert alert-success">
                <i class="fas fa-check-circle me-2"></i>
                Alignment completed successfully
            </div>
            <div class="result-item">
                <h6>Alignment Summary</h6>
                ${summaryHtml}
            </div>
        `;
    }

    displayVisualizationResults(data) {
        const container = document.getElementById('plotResults');
        const plotPaths = data.plot_paths;
        
        if (plotPaths && Object.keys(plotPaths).length > 0) {
            const plotsHtml = Object.entries(plotPaths).map(([name, path]) => `
                <div class="result-item">
                    <h6>${name.replace(/_/g, ' ').toUpperCase()}</h6>
                    <a href="${path}" target="_blank" class="btn btn-sm btn-primary">View Plot</a>
                </div>
            `).join('');

            container.innerHTML = `
                <div class="alert alert-success">
                    <i class="fas fa-check-circle me-2"></i>
                    Visualizations generated successfully
                </div>
                ${plotsHtml}
                ${data.summary_report ? `
                    <div class="result-item">
                        <h6>Summary Report</h6>
                        <a href="${data.summary_report}" target="_blank" class="btn btn-sm btn-primary">View Report</a>
                    </div>
                ` : ''}
            `;
        } else {
            container.innerHTML = `
                <div class="alert alert-warning">
                    <i class="fas fa-exclamation-triangle me-2"></i>
                    No visualization data available
                </div>
            `;
        }
    }

    updateStatus() {
        document.getElementById('dataStatus').textContent = this.analysisStatus.data ? 'Completed' : 'Not Started';
        document.getElementById('qcStatus').textContent = this.analysisStatus.qc ? 'Completed' : 'Not Started';
        document.getElementById('alignmentStatus').textContent = this.analysisStatus.alignment ? 'Completed' : 'Not Started';
        document.getElementById('vizStatus').textContent = this.analysisStatus.visualization ? 'Completed' : 'Not Started';
    }

    updatePipelineStep(step, completed = false) {
        const stepElement = document.getElementById(`step${step}`);
        if (stepElement) {
            stepElement.classList.remove('active');
            if (completed) {
                stepElement.classList.add('completed');
            } else {
                stepElement.classList.add('active');
            }
        }
    }

    onTabChange(targetId) {
        // Update active pipeline step based on tab
        const stepMap = {
            '#data': 1,
            '#qc': 2,
            '#alignment': 3,
            '#visualization': 4,
            '#results': 4
        };

        const step = stepMap[targetId];
        if (step) {
            // Reset all steps
            document.querySelectorAll('.pipeline-step').forEach(el => {
                el.classList.remove('active');
            });
            
            // Activate current step
            this.updatePipelineStep(step, false);
        }
    }

    showLoading(message = 'Processing...') {
        document.getElementById('loadingMessage').textContent = message;
        const modal = new bootstrap.Modal(document.getElementById('loadingModal'));
        modal.show();
    }

    hideLoading() {
        const modal = bootstrap.Modal.getInstance(document.getElementById('loadingModal'));
        if (modal) {
            modal.hide();
        }
    }

    showProgress() {
        const progress = document.getElementById('downloadProgress');
        progress.style.display = 'block';
        
        // Simulate progress
        const progressBar = progress.querySelector('.progress-bar');
        let width = 0;
        const interval = setInterval(() => {
            if (width >= 90) {
                clearInterval(interval);
            } else {
                width += 10;
                progressBar.style.width = width + '%';
            }
        }, 500);
    }

    hideProgress() {
        const progress = document.getElementById('downloadProgress');
        progress.style.display = 'none';
        progress.querySelector('.progress-bar').style.width = '0%';
    }

    showAlert(message, type = 'info') {
        const alertDiv = document.createElement('div');
        alertDiv.className = `alert alert-${type} alert-dismissible fade show`;
        alertDiv.innerHTML = `
            ${message}
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        `;

        // Insert at the top of the container
        const container = document.querySelector('.container');
        container.insertBefore(alertDiv, container.firstChild);

        // Auto-dismiss after 5 seconds
        setTimeout(() => {
            if (alertDiv.parentNode) {
                alertDiv.remove();
            }
        }, 5000);
    }
}

// Initialize the platform when the page loads
document.addEventListener('DOMContentLoaded', () => {
    window.rnaseqPlatform = new RNAseqPlatform();
}); 