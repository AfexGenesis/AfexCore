// ===== GENE COMPARATOR - DIAMOND-INSPIRED TOOL =====
// Direct Python script execution from Electron renderer

class PythonGeneComparator {
    constructor() {
        this.queryFile = null;
        this.subjectFile = null;
        this.currentResults = null;
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUploads();
        console.log('🧬 Python Gene Comparator initialized');
    }

    bindEvents() {
        // Compare button
        const compareBtn = document.getElementById('gene-compare-btn');
        if (compareBtn) {
            compareBtn.addEventListener('click', () => this.performComparison());
        }

        // Clear button
        const clearBtn = document.getElementById('gene-clear-btn');
        if (clearBtn) {
            clearBtn.addEventListener('click', () => this.clearAll());
        }

        // Copy buttons
        const copyAlignmentBtn = document.getElementById('copy-alignment-btn');
        if (copyAlignmentBtn) {
            copyAlignmentBtn.addEventListener('click', () => this.copyAlignment());
        }

        const copyStatsBtn = document.getElementById('copy-stats-btn');
        if (copyStatsBtn) {
            copyStatsBtn.addEventListener('click', () => this.copyStatistics());
        }

        // Download report button
        const downloadBtn = document.getElementById('download-report-btn');
        if (downloadBtn) {
            downloadBtn.addEventListener('click', () => this.downloadReport());
        }

        // Input validation
        const queryInput = document.getElementById('query-sequence');
        const subjectInput = document.getElementById('subject-sequence');
        
        if (queryInput) {
            queryInput.addEventListener('input', () => this.updateSequenceInfo('query'));
        }
        if (subjectInput) {
            subjectInput.addEventListener('input', () => this.updateSequenceInfo('subject'));
        }
    }

    setupFileUploads() {
        // Query file upload
        this.setupFileUpload('query-upload-area', 'query-file-input', (file) => {
            this.queryFile = file;
            this.handleFileUpload(file, 'query');
        });

        // Subject file upload
        this.setupFileUpload('subject-upload-area', 'subject-file-input', (file) => {
            this.subjectFile = file;
            this.handleFileUpload(file, 'subject');
        });
    }

    setupFileUpload(areaId, inputId, callback) {
        const uploadArea = document.getElementById(areaId);
        const fileInput = document.getElementById(inputId);

        if (!uploadArea || !fileInput) return;

        // Click to upload
        uploadArea.addEventListener('click', () => fileInput.click());

        // File input change
        fileInput.addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (file) {
                callback(file);
            }
        });

        // Drag and drop
        uploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            uploadArea.classList.add('border-emerald-500/70');
        });

        uploadArea.addEventListener('dragleave', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('border-emerald-500/70');
        });

        uploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('border-emerald-500/70');
            
            const files = e.dataTransfer.files;
            if (files.length > 0) {
                callback(files[0]);
            }
        });
    }

    async handleFileUpload(file, type) {
        try {
            // Check file size (warn for files > 10MB)
            const maxSize = 10 * 1024 * 1024; // 10MB
            if (file.size > maxSize) {
                const proceed = confirm(`File "${file.name}" is ${this.formatFileSize(file.size)}. Large files may take longer to process. Continue?`);
                if (!proceed) return;
            }

            // Show loading
            this.showFileLoading(type, true);

            // Read file content
            const fileContent = await this.readFileContent(file);
            
            // Update textarea with file content
            const textarea = document.getElementById(`${type}-sequence`);
            if (textarea) {
                textarea.value = fileContent;
                this.updateSequenceInfo(type);
            }

            // Update upload status
            this.updateUploadStatus(type, file.name, true);
            this.showNotification(`File "${file.name}" loaded successfully`, 'success');

        } catch (error) {
            console.error('File upload error:', error);
            this.updateUploadStatus(type, 'Upload failed', false);
            this.showNotification('Failed to load file. Please try again.', 'error');
        } finally {
            this.showFileLoading(type, false);
        }
    }

    readFileContent(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = (e) => resolve(e.target.result);
            reader.onerror = (e) => reject(e);
            reader.readAsText(file);
        });
    }

    updateSequenceInfo(type) {
        const textarea = document.getElementById(`${type}-sequence`);
        const infoElement = document.getElementById(`${type}-info`);
        
        if (textarea && infoElement) {
            const sequence = textarea.value.replace(/\s/g, '').replace(/>/g, '');
            const length = sequence.length;
            const gcContent = this.calculateGCContent(sequence);
            
            infoElement.innerHTML = `
                <span class="text-emerald-400">Length:</span> ${length} bp | 
                <span class="text-emerald-400">GC:</span> ${gcContent}%
            `;
        }
    }

    calculateGCContent(sequence) {
        const cleanSeq = sequence.toUpperCase().replace(/[^ATGC]/g, '');
        if (cleanSeq.length === 0) return 0;
        
        const gcCount = (cleanSeq.match(/[GC]/g) || []).length;
        return ((gcCount / cleanSeq.length) * 100).toFixed(1);
    }

    updateUploadStatus(type, fileName, success) {
        const statusElement = document.getElementById(`${type}-upload-status`);
        if (statusElement) {
            statusElement.innerHTML = success 
                ? `<span class="text-emerald-400">✓ ${fileName}</span>`
                : `<span class="text-red-400">✗ ${fileName}</span>`;
        }
    }

    showFileLoading(type, show) {
        const uploadArea = document.getElementById(`${type}-upload-area`);
        if (uploadArea) {
            if (show) {
                uploadArea.classList.add('opacity-50');
                uploadArea.innerHTML = '<div class="text-center">📄 Loading file...</div>';
            } else {
                uploadArea.classList.remove('opacity-50');
                uploadArea.innerHTML = `
                    <div class="text-center">
                        <div class="text-2xl mb-2">📁</div>
                        <div class="text-sm text-gray-300">Click or drag ${type} sequence file</div>
                        <div class="text-xs text-gray-400 mt-1">FASTA, GenBank, or text</div>
                    </div>
                `;
            }
        }
    }

    async performComparison() {
        try {
            // Get sequences
            const querySequence = document.getElementById('query-sequence')?.value.trim();
            const subjectSequence = document.getElementById('subject-sequence')?.value.trim();

            if (!querySequence || !subjectSequence) {
                this.showNotification('Please provide both query and subject sequences', 'error');
                return;
            }

            // Show loading state
            this.setLoadingState(true);
            this.hideResults();

            // Call Python script
            const result = await this.callPythonScript(querySequence, subjectSequence);

            if (result.success) {
                this.currentResults = result;
                this.displayResults(result);
                this.showNotification('Gene comparison completed successfully!', 'success');
            } else {
                throw new Error(result.error || 'Comparison failed');
            }

        } catch (error) {
            console.error('Comparison error:', error);
            this.showNotification(`Comparison failed: ${error.message}`, 'error');
        } finally {
            this.setLoadingState(false);
        }
    }

    async callPythonScript(querySequence, subjectSequence) {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Python script (without sequences)
            const scriptPath = path.join(__dirname, 'assets', 'gene-comparator.py');
            const args = [
                scriptPath,
                '--query', 'STDIN_QUERY',
                '--subject', 'STDIN_SUBJECT',
                '--format', 'auto'
            ];

            console.log('🐍 Calling Python Gene Comparator:', scriptPath);

            // Spawn Python process
            const pythonProcess = spawn('python', args, {
                cwd: path.join(__dirname, '..'),
                stdio: ['pipe', 'pipe', 'pipe']
            });

            let stdout = '';
            let stderr = '';

            pythonProcess.stdout.on('data', (data) => {
                stdout += data.toString();
            });

            pythonProcess.stderr.on('data', (data) => {
                stderr += data.toString();
            });

            // Send sequences via stdin as JSON
            const inputData = JSON.stringify({
                query: querySequence,
                subject: subjectSequence
            });
            
            pythonProcess.stdin.write(inputData);
            pythonProcess.stdin.end();

            pythonProcess.on('close', (code) => {
                console.log(`🐍 Python process exited with code: ${code}`);
                
                if (code === 0) {
                    try {
                        const result = JSON.parse(stdout);
                        console.log('✅ Python result:', result);
                        resolve(result);
                    } catch (parseError) {
                        console.error('❌ Failed to parse Python output:', parseError);
                        console.error('Raw stdout:', stdout);
                        reject(new Error(`Failed to parse Python output: ${parseError.message}`));
                    }
                } else {
                    console.error('❌ Python script failed:', stderr);
                    reject(new Error(`Python script failed: ${stderr || 'Unknown error'}`));
                }
            });

            pythonProcess.on('error', (error) => {
                console.error('❌ Failed to start Python process:', error);
                reject(new Error(`Failed to start Python process: ${error.message}`));
            });
        });
    }

    displayResults(result) {
        // Show results container
        const resultsContainer = document.getElementById('comparison-results');
        if (resultsContainer) {
            resultsContainer.style.display = 'block';
        }

        // Display summary
        this.displaySummary(result);
        
        // Display statistics
        this.displayStatistics(result.statistics);
        
        // Display composition
        this.displayComposition(result.composition);
        
        // Display alignment
        this.displayAlignment(result.alignment_display);
    }

    displaySummary(result) {
        const summaryElement = document.getElementById('comparison-summary');
        if (!summaryElement) return;

        const stats = result.statistics;
        const summary = result.summary;
        
        // Determine color based on identity
        let identityColor = 'text-red-400';
        if (stats.identity_percent >= 80) identityColor = 'text-emerald-400';
        else if (stats.identity_percent >= 60) identityColor = 'text-yellow-400';

        summaryElement.innerHTML = `
            <div class="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                <div class="bg-slate-700/30 rounded-lg p-4 text-center">
                    <div class="text-2xl font-bold ${identityColor}">${stats.identity_percent}%</div>
                    <div class="text-sm text-gray-300">Identity</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-4 text-center">
                    <div class="text-2xl font-bold text-cyan-400">${stats.similarity_percent}%</div>
                    <div class="text-sm text-gray-300">Similarity</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-4 text-center">
                    <div class="text-2xl font-bold text-blue-400">${Math.round((stats.query_coverage + stats.subject_coverage) / 2)}%</div>
                    <div class="text-sm text-gray-300">Avg Coverage</div>
                </div>
            </div>
            <div class="bg-slate-800/50 rounded-lg p-4">
                <h4 class="text-lg font-semibold text-white mb-2">${summary.description}</h4>
                <p class="text-gray-300">${summary.functional_prediction}</p>
            </div>
        `;
    }

    displayStatistics(stats) {
        const statsElement = document.getElementById('detailed-statistics');
        if (!statsElement) return;

        statsElement.innerHTML = `
            <div class="grid grid-cols-2 md:grid-cols-4 gap-4">
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">Alignment Length</div>
                    <div class="text-lg font-semibold text-white">${stats.alignment_length}</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">Matches</div>
                    <div class="text-lg font-semibold text-emerald-400">${stats.matches}</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">Mismatches</div>
                    <div class="text-lg font-semibold text-red-400">${stats.mismatches}</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">Gaps</div>
                    <div class="text-lg font-semibold text-yellow-400">${stats.gaps}</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">Query Coverage</div>
                    <div class="text-lg font-semibold text-cyan-400">${stats.query_coverage}%</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">Subject Coverage</div>
                    <div class="text-lg font-semibold text-cyan-400">${stats.subject_coverage}%</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">Bit Score</div>
                    <div class="text-lg font-semibold text-blue-400">${stats.bit_score}</div>
                </div>
                <div class="bg-slate-700/30 rounded-lg p-3">
                    <div class="text-sm text-gray-400">E-value</div>
                    <div class="text-lg font-semibold text-purple-400">${stats.e_value}</div>
                </div>
            </div>
        `;
    }

    displayComposition(composition) {
        const compositionElement = document.getElementById('sequence-composition');
        if (!compositionElement) return;

        const createCompositionCard = (title, data, color) => `
            <div class="bg-slate-700/30 rounded-lg p-4">
                <h4 class="text-lg font-semibold ${color} mb-3">${title}</h4>
                <div class="grid grid-cols-2 gap-2 text-sm">
                    <div>Length: <span class="text-white font-mono">${data.length} bp</span></div>
                    <div>GC Content: <span class="text-white font-mono">${data.gc_content}%</span></div>
                    <div>A: <span class="text-white font-mono">${data.a_count}</span></div>
                    <div>T: <span class="text-white font-mono">${data.t_count}</span></div>
                    <div>G: <span class="text-white font-mono">${data.g_count}</span></div>
                    <div>C: <span class="text-white font-mono">${data.c_count}</span></div>
                </div>
            </div>
        `;

        compositionElement.innerHTML = `
            <div class="grid grid-cols-1 md:grid-cols-2 gap-4">
                ${createCompositionCard('Query Sequence', composition.query, 'text-emerald-400')}
                ${createCompositionCard('Subject Sequence', composition.subject, 'text-cyan-400')}
            </div>
        `;
    }

    displayAlignment(alignmentDisplay) {
        const alignmentElement = document.getElementById('alignment-display');
        if (!alignmentElement) return;

        alignmentElement.innerHTML = `
            <pre class="text-xs font-mono text-gray-300 whitespace-pre-wrap overflow-x-auto bg-slate-900/50 p-4 rounded-lg border border-slate-600/30">${alignmentDisplay}</pre>
        `;
    }

    hideResults() {
        const resultsContainer = document.getElementById('comparison-results');
        if (resultsContainer) {
            resultsContainer.style.display = 'none';
        }
    }

    clearAll() {
        // Clear textareas
        const queryInput = document.getElementById('query-sequence');
        const subjectInput = document.getElementById('subject-sequence');
        
        if (queryInput) queryInput.value = '';
        if (subjectInput) subjectInput.value = '';

        // Clear file references
        this.queryFile = null;
        this.subjectFile = null;
        this.currentResults = null;

        // Reset upload areas
        this.showFileLoading('query', false);
        this.showFileLoading('subject', false);

        // Clear upload status
        this.updateUploadStatus('query', '', false);
        this.updateUploadStatus('subject', '', false);

        // Clear sequence info
        this.updateSequenceInfo('query');
        this.updateSequenceInfo('subject');

        // Hide results
        this.hideResults();

        this.showNotification('All data cleared', 'success');
    }

    setLoadingState(isLoading) {
        const compareBtn = document.getElementById('gene-compare-btn');
        if (compareBtn) {
            if (isLoading) {
                compareBtn.disabled = true;
                compareBtn.innerHTML = '🔄 Comparing...';
                compareBtn.classList.add('opacity-50');
            } else {
                compareBtn.disabled = false;
                compareBtn.innerHTML = '🧬 Compare Genes';
                compareBtn.classList.remove('opacity-50');
            }
        }
    }

    async copyAlignment() {
        if (!this.currentResults?.alignment_display) {
            this.showNotification('No alignment to copy', 'error');
            return;
        }

        try {
            await navigator.clipboard.writeText(this.currentResults.alignment_display);
            this.showNotification('Alignment copied to clipboard!', 'success');
        } catch (error) {
            console.error('Copy failed:', error);
            this.showNotification('Failed to copy alignment', 'error');
        }
    }

    async copyStatistics() {
        if (!this.currentResults?.statistics) {
            this.showNotification('No statistics to copy', 'error');
            return;
        }

        const stats = this.currentResults.statistics;
        const statsText = `Gene Comparison Statistics:
Identity: ${stats.identity_percent}%
Similarity: ${stats.similarity_percent}%
Query Coverage: ${stats.query_coverage}%
Subject Coverage: ${stats.subject_coverage}%
Alignment Length: ${stats.alignment_length}
Matches: ${stats.matches}
Mismatches: ${stats.mismatches}
Gaps: ${stats.gaps}
Bit Score: ${stats.bit_score}
E-value: ${stats.e_value}`;

        try {
            await navigator.clipboard.writeText(statsText);
            this.showNotification('Statistics copied to clipboard!', 'success');
        } catch (error) {
            console.error('Copy failed:', error);
            this.showNotification('Failed to copy statistics', 'error');
        }
    }

    downloadReport() {
        if (!this.currentResults) {
            this.showNotification('No results to download', 'error');
            return;
        }

        const report = this.generateReport();
        const blob = new Blob([report], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        
        const a = document.createElement('a');
        a.href = url;
        a.download = `gene_comparison_report_${new Date().toISOString().slice(0, 10)}.txt`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);

        this.showNotification('Report downloaded successfully!', 'success');
    }

    generateReport() {
        const result = this.currentResults;
        const stats = result.statistics;
        const composition = result.composition;

        return `AFEX Gene Comparator Report
Generated: ${new Date().toLocaleString()}
========================================

SUMMARY
${result.summary.description}
${result.summary.functional_prediction}

ALIGNMENT STATISTICS
Identity: ${stats.identity_percent}%
Similarity: ${stats.similarity_percent}%
Query Coverage: ${stats.query_coverage}%
Subject Coverage: ${stats.subject_coverage}%
Alignment Length: ${stats.alignment_length}
Matches: ${stats.matches}
Mismatches: ${stats.mismatches}
Gaps: ${stats.gaps} (${stats.gap_opens} opens)
Bit Score: ${stats.bit_score}
E-value: ${stats.e_value}

SEQUENCE COMPOSITION
Query Sequence:
  Length: ${composition.query.length} bp
  GC Content: ${composition.query.gc_content}%
  A: ${composition.query.a_count}, T: ${composition.query.t_count}, G: ${composition.query.g_count}, C: ${composition.query.c_count}

Subject Sequence:
  Length: ${composition.subject.length} bp
  GC Content: ${composition.subject.gc_content}%
  A: ${composition.subject.a_count}, T: ${composition.subject.t_count}, G: ${composition.subject.g_count}, C: ${composition.subject.c_count}

ALIGNMENT
${result.alignment_display}

========================================
Generated by AFEX Gene Comparator
DIAMOND-inspired gene comparison tool`;
    }

    formatFileSize(bytes) {
        if (bytes === 0) return '0 Bytes';
        const k = 1024;
        const sizes = ['Bytes', 'KB', 'MB', 'GB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));
        return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
    }

    showNotification(message, type = 'info') {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `fixed top-4 right-4 z-50 px-6 py-3 rounded-lg shadow-lg transition-all duration-300 ${
            type === 'success' ? 'bg-emerald-600 text-white' :
            type === 'error' ? 'bg-red-600 text-white' :
            'bg-blue-600 text-white'
        }`;
        notification.textContent = message;

        document.body.appendChild(notification);

        // Remove after 3 seconds
        setTimeout(() => {
            notification.style.opacity = '0';
            setTimeout(() => {
                if (notification.parentNode) {
                    notification.parentNode.removeChild(notification);
                }
            }, 300);
        }, 3000);
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    // Only initialize if we're on the gene comparator page
    if (document.getElementById('genome-comparator-page')) {
        window.geneComparator = new PythonGeneComparator();
    }
});