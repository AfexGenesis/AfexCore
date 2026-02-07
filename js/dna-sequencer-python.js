// ===== DNA SEQUENCER - PYTHON INTEGRATION =====
// Advanced DNA sequence analysis and visualization

class DNASequencerPython {
    constructor() {
        this.currentSequence = null;
        this.currentResults = null;
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUpload();
        this.updateSequenceLength();
        this.updateAnalyzeButton();
        console.log('DNA Sequencer initialized');
    }

    bindEvents() {
        // Analyze button - we need to add this to the UI
        const analyzeBtn = document.getElementById('sequencer-analyze-btn');
        if (analyzeBtn) {
            analyzeBtn.addEventListener('click', () => this.performAnalysis());
        }

        // Clear button
        const clearBtn = document.getElementById('sequencer-clear-btn');
        if (clearBtn) {
            clearBtn.addEventListener('click', () => this.clearInput());
        }

        // Manual input textarea
        const manualInput = document.getElementById('sequencer-manual-input');
        if (manualInput) {
            manualInput.addEventListener('input', () => {
                this.updateSequenceLength();
                this.updateRealtimeGC();
                this.updateAnalyzeButton();
            });
            
            // Also add paste event for immediate display
            manualInput.addEventListener('paste', () => {
                setTimeout(() => {
                    this.updateSequenceLength();
                    this.updateRealtimeGC();
                    this.updateAnalyzeButton();
                }, 100);
            });
        }

        // Copy results button
        const copyBtn = document.getElementById('copy-sequencer-results');
        if (copyBtn) {
            copyBtn.addEventListener('click', () => this.copyResults());
        }

        // Edit Lock Mode button
        const editLockBtn = document.getElementById('edit-lock-mode');
        if (editLockBtn) {
            editLockBtn.addEventListener('click', () => this.toggleEditLockMode());
        }
    }

    setupFileUpload() {
        const fileInput = document.getElementById('sequencer-file-upload');
        const dropZone = fileInput?.parentElement;

        if (!fileInput || !dropZone) return;

        // File input change event
        fileInput.addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (file) {
                this.handleFileUpload(file);
            }
        });

        // Drag and drop events
        dropZone.addEventListener('dragover', (e) => {
            e.preventDefault();
            dropZone.classList.add('border-emerald-500/70');
        });

        dropZone.addEventListener('dragleave', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-emerald-500/70');
        });

        dropZone.addEventListener('drop', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-emerald-500/70');
            
            const files = e.dataTransfer.files;
            if (files.length > 0) {
                this.handleFileUpload(files[0]);
            }
        });
    }

    async handleFileUpload(file) {
        try {
            // Check for large files (>100MB)
            const maxSafeSize = 100 * 1024 * 1024; // 100MB
            const maxRecommendedSize = 5 * 1024 * 1024 * 1024; // 5GB recommended limit
            const maxAllowedSize = 10 * 1024 * 1024 * 1024; // 10GB absolute limit
            
            if (file.size > maxAllowedSize) {
                this.showNotification('File size exceeds maximum limit of 10GB', 'error');
                return;
            }

            if (file.size > maxSafeSize) {
                const proceed = await this.showLargeFileWarning(file, 'DNA sequence', maxRecommendedSize);
                if (!proceed) {
                    return;
                }
            }

            // Show loading for large files
            if (file.size > maxSafeSize) {
                this.showLoadingOverlay('Reading large sequence file...');
            }

            const fileContent = await this.readFileContent(file);
            const fileName = file.name.toLowerCase();
            let fileFormat = 'txt';

            // Determine file format
            if (fileName.endsWith('.fasta') || fileName.endsWith('.fa')) {
                fileFormat = 'fasta';
            } else if (fileName.endsWith('.gb') || fileName.endsWith('.genbank') || fileName.endsWith('.gbk')) {
                fileFormat = 'genbank';
            }

            // Parse file and load sequence
            const result = await this.callPythonScript('', fileContent, fileFormat, {});
            
            if (result.success) {
                // Load the sequence into the manual input area
                const manualInput = document.getElementById('sequencer-manual-input');
                if (manualInput && result.clean_sequence) {
                    manualInput.value = result.clean_sequence;
                    this.updateSequenceLength();
                    this.updateRealtimeGC();
                    this.updateAnalyzeButton();
                }

                // Automatically display the sequence
                this.displaySequence(result);
                
                this.hideLoadingOverlay();
                this.showNotification(`File "${file.name}" loaded successfully`, 'success');
            } else {
                throw new Error(result.error || 'Failed to parse file');
            }

        } catch (error) {
            console.error('File upload error:', error);
            this.hideLoadingOverlay();
            this.showNotification('Failed to load file. Please try again.', 'error');
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

    updateSequenceLength() {
        const manualInput = document.getElementById('sequencer-manual-input');
        const lengthDisplay = document.getElementById('sequencer-sequence-length');
        
        if (manualInput && lengthDisplay) {
            const sequence = manualInput.value.replace(/\s/g, '');
            lengthDisplay.textContent = `Length: ${sequence.length}`;
        }
    }

    updateRealtimeGC() {
        const manualInput = document.getElementById('sequencer-manual-input');
        const gcDisplay = document.getElementById('sequencer-input-gc');
        
        if (manualInput && gcDisplay) {
            const sequence = manualInput.value.replace(/\s/g, '').toUpperCase();
            if (sequence.length > 0) {
                const gcCount = (sequence.match(/[GC]/g) || []).length;
                const gcContent = ((gcCount / sequence.length) * 100).toFixed(1);
                gcDisplay.textContent = `${gcContent}%`;
            } else {
                gcDisplay.textContent = '0%';
            }
        }
    }

    updateAnalyzeButton() {
        const manualInput = document.getElementById('sequencer-manual-input');
        const analyzeBtn = document.getElementById('sequencer-analyze-btn');
        
        if (manualInput && analyzeBtn) {
            const sequence = manualInput.value.trim();
            analyzeBtn.disabled = !sequence;
        }
    }

    clearInput() {
        const manualInput = document.getElementById('sequencer-manual-input');
        if (manualInput) {
            manualInput.value = '';
            this.updateSequenceLength();
            this.updateRealtimeGC();
            this.updateAnalyzeButton();
            this.hideResults();
        }
    }

    async performAnalysis() {
        try {
            // Get input sequence
            const manualInput = document.getElementById('sequencer-manual-input');
            const sequence = manualInput?.value.trim();

            if (!sequence) {
                this.showNotification('Please enter a DNA sequence or upload a file', 'error');
                return;
            }

            // Get analysis options
            const options = this.getAnalysisOptions();

            // Show loading state
            this.setLoadingState(true);

            // Call Python script
            const result = await this.callPythonScript(sequence, null, null, options);

            if (result.success) {
                this.displaySequence(result);
                this.showNotification('Sequence analysis completed!', 'success');
            } else {
                throw new Error(result.error || 'Analysis failed');
            }

        } catch (error) {
            console.error('Analysis error:', error);
            this.showNotification(`Analysis failed: ${error.message}`, 'error');
        } finally {
            this.setLoadingState(false);
        }
    }

    getAnalysisOptions() {
        const options = {};

        // Reading frames
        const readingFrames = [];
        if (document.getElementById('reading-frame-1')?.checked) readingFrames.push('+1');
        if (document.getElementById('reading-frame-2')?.checked) readingFrames.push('+2');
        if (document.getElementById('reading-frame-3')?.checked) readingFrames.push('+3');
        if (document.getElementById('reading-frame-all')?.checked) {
            readingFrames.push('+1', '+2', '+3', '-1', '-2', '-3');
        }
        if (readingFrames.length > 0) {
            options.reading_frames = readingFrames;
        }

        // GC content analysis
        if (document.getElementById('analyze-gc-content')?.checked) {
            options.analyze_gc_content = true;
        }

        // Show reverse complement
        if (document.getElementById('show-reverse-complement')?.checked) {
            options.show_reverse_complement = true;
        }

        // Restriction enzyme finder
        const restrictionSelect = document.getElementById('restriction-enzymes-finder');
        if (restrictionSelect && restrictionSelect.value) {
            options.restriction_enzyme = restrictionSelect.value;
        }

        // Color scheme
        const colorScheme = document.getElementById('color-scheme');
        if (colorScheme && colorScheme.value) {
            options.color_scheme = colorScheme.value;
        }

        // Amino acids finder
        const aminoAcidsFinder = document.getElementById('amino-acids-finder');
        if (aminoAcidsFinder && aminoAcidsFinder.value) {
            options.amino_acids_finder = aminoAcidsFinder.value;
        }

        // GC Content Optimizer
        const gcOptimizer = document.getElementById('sequencer-gc-optimizer');
        if (gcOptimizer && gcOptimizer.value) {
            if (gcOptimizer.value === 'custom') {
                // Get custom slider value
                const gcSlider = document.getElementById('gc-target-slider');
                if (gcSlider) {
                    options.gc_optimizer = gcSlider.value; // Send the exact percentage
                }
            } else {
                options.gc_optimizer = gcOptimizer.value; // Send the range (e.g., "20-30")
            }
        }

        // Codon Usage Optimizer
        const codonOptimizer = document.getElementById('sequencer-codon-optimizer');
        if (codonOptimizer && codonOptimizer.value) {
            options.codon_optimizer = codonOptimizer.value;
        }

        return options;
    }

    async callPythonScript(sequence, fileContent = null, fileFormat = null, options = {}) {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Python script (without large content)
            const scriptPath = path.join(__dirname, 'assets', 'dna-sequencer.py');
            const args = [scriptPath];

            // Prepare input data as JSON to send via stdin
            const inputData = {
                sequence: sequence,
                fileContent: fileContent,
                fileFormat: fileFormat,
                options: options
            };

            console.log('Calling DNA Sequencer Python script:', scriptPath);

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

            // Send input data via stdin
            try {
                pythonProcess.stdin.write(JSON.stringify(inputData));
                pythonProcess.stdin.end();
            } catch (error) {
                console.error('Failed to write to Python stdin:', error);
                reject(new Error(`Failed to send data to Python: ${error.message}`));
                return;
            }

            pythonProcess.on('close', (code) => {
                console.log(`Python process exited with code: ${code}`);
                
                if (code === 0) {
                    try {
                        const result = JSON.parse(stdout);
                        console.log('Python result:', result);
                        resolve(result);
                    } catch (parseError) {
                        console.error('Failed to parse Python output:', parseError);
                        console.error('Raw stdout:', stdout);
                        reject(new Error(`Failed to parse Python output: ${parseError.message}`));
                    }
                } else {
                    console.error('Python script failed:', stderr);
                    reject(new Error(`Python script failed: ${stderr || 'Unknown error'}`));
                }
            });

            pythonProcess.on('error', (error) => {
                console.error('Failed to start Python process:', error);
                reject(new Error(`Failed to start Python process: ${error.message}`));
            });
        });
    }

    displaySequence(result) {
        // Store current results
        this.currentResults = result;
        this.currentSequence = result.clean_sequence;

        // Hide empty state
        const emptyState = document.getElementById('sequencer-empty-state');
        if (emptyState) {
            emptyState.style.display = 'none';
        }

        // Show results container
        const resultsContainer = document.getElementById('sequencer-results');
        if (resultsContainer) {
            resultsContainer.classList.remove('hidden');
            resultsContainer.style.display = 'block';
        }

        // Display formatted sequence
        const outputDiv = document.getElementById('sequencer-output');
        if (outputDiv && result.results.formatted_sequence) {
            let html = '<div class="space-y-1">';
            
            // Add sequence info header
            if (result.file_info) {
                html += `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-cyan-500/30">`;
                html += `<h4 class="text-cyan-400 font-medium mb-1">📄 ${result.file_info.id}</h4>`;
                html += `<p class="text-gray-400 text-sm">${result.file_info.description}</p>`;
                html += `<p class="text-gray-500 text-xs mt-1">Format: ${result.file_info.format.toUpperCase()}</p>`;
                html += `</div>`;
            }

            // Add sequence statistics
            const seqInfo = result.results.sequence_info;
            html += `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-emerald-500/30">`;
            html += `<h4 class="text-emerald-400 font-medium mb-2">📊 Sequence Statistics</h4>`;
            html += `<div class="grid grid-cols-2 md:grid-cols-4 gap-3 text-sm">`;
            html += `<div><span class="text-gray-400">Length:</span> <span class="text-white font-mono">${seqInfo.length.toLocaleString()} bp</span></div>`;
            html += `<div><span class="text-gray-400">GC Content:</span> <span class="text-white font-mono">${seqInfo.gc_content}%</span></div>`;
            html += `<div><span class="text-red-400">A:</span> <span class="text-white font-mono">${seqInfo.composition.A.toLocaleString()}</span></div>`;
            html += `<div><span class="text-blue-400">T:</span> <span class="text-white font-mono">${seqInfo.composition.T.toLocaleString()}</span></div>`;
            html += `<div><span class="text-green-400">G:</span> <span class="text-white font-mono">${seqInfo.composition.G.toLocaleString()}</span></div>`;
            html += `<div><span class="text-yellow-400">C:</span> <span class="text-white font-mono">${seqInfo.composition.C.toLocaleString()}</span></div>`;
            if (seqInfo.composition.N > 0) {
                html += `<div><span class="text-gray-400">N:</span> <span class="text-white font-mono">${seqInfo.composition.N.toLocaleString()}</span></div>`;
            }
            html += `</div></div>`;

            // Add formatted sequence lines
            html += `<div class="mb-4 w-full">`;
            
            // Check if reverse complement is being shown
            const showingReverseComplement = this.getAnalysisOptions().show_reverse_complement;
            const sequenceTitle = showingReverseComplement ? 
                `🧬 DNA Double Strand (5'→3' and 3'→5')` : 
                `🧬 DNA Sequence`;
            
            html += `<h4 class="text-cyan-400 font-medium mb-2">${sequenceTitle}</h4>`;
            html += `<div class="bg-black/30 rounded-lg p-4 font-mono text-sm leading-relaxed w-full max-w-full">`;
            
            // Use colored sequence if available, otherwise apply default coloring
            const sequenceLines = result.results.colored_sequence || result.results.formatted_sequence;
            const useCustomColors = result.results.colored_sequence;
            
            sequenceLines.forEach((line, index) => {
                if (useCustomColors) {
                    // Already colored by Python
                    if (showingReverseComplement) {
                        // Handle separated line numbers and sequences for perfect alignment
                        if (line.includes('|')) {
                            const [lineNum, sequence] = line.split('|');
                            if (lineNum.trim().match(/^\d+$/)) {
                                // Forward strand line
                                html += `<div class="mb-1 w-full break-all flex">
                                    <span class="text-gray-500 text-xs mr-2 flex-shrink-0">5'</span>
                                    <span class="text-gray-400 mr-2 flex-shrink-0 font-mono">${lineNum}:</span>
                                    <span class="font-mono">${sequence}</span>
                                    <span class="text-gray-500 text-xs ml-2 flex-shrink-0">3'</span>
                                </div>`;
                            } else {
                                // Reverse complement line
                                html += `<div class="mb-1 w-full break-all flex text-gray-300">
                                    <span class="text-gray-500 text-xs mr-2 flex-shrink-0">3'</span>
                                    <span class="text-transparent mr-2 flex-shrink-0 font-mono">${lineNum.replace(/\s/g, '0')}:</span>
                                    <span class="font-mono">${sequence}</span>
                                    <span class="text-gray-500 text-xs ml-2 flex-shrink-0">5'</span>
                                </div>`;
                            }
                        } else if (line.trim() === '') {
                            // Empty separator line
                            html += `<div class="mb-1 w-full">${line}</div>`;
                        } else {
                            // Fallback for other lines
                            html += `<div class="mb-1 w-full break-all">${line}</div>`;
                        }
                    } else {
                        // Regular single strand display
                        if (line.includes('|')) {
                            const [lineNum, sequence] = line.split('|');
                            html += `<div class="mb-1 w-full break-all">
                                <span class="text-gray-400 mr-2 font-mono">${lineNum}:</span>
                                <span class="font-mono">${sequence}</span>
                            </div>`;
                        } else {
                            html += `<div class="mb-1 w-full break-all">${line}</div>`;
                        }
                    }
                } else {
                    // Apply default coloring
                    let coloredLine = line.replace(/([ATGC])/g, (match) => {
                        switch(match) {
                            case 'A': return `<span class="text-red-400">${match}</span>`;
                            case 'T': return `<span class="text-blue-400">${match}</span>`;
                            case 'G': return `<span class="text-green-400">${match}</span>`;
                            case 'C': return `<span class="text-yellow-400">${match}</span>`;
                            default: return match;
                        }
                    });
                    
                    if (showingReverseComplement) {
                        // Handle separated line numbers and sequences for perfect alignment
                        if (line.includes('|')) {
                            const [lineNum, sequence] = line.split('|');
                            // Apply coloring to the sequence part
                            let coloredSequence = sequence.replace(/([ATGC])/g, (match) => {
                                switch(match) {
                                    case 'A': return `<span class="text-red-400">${match}</span>`;
                                    case 'T': return `<span class="text-blue-400">${match}</span>`;
                                    case 'G': return `<span class="text-green-400">${match}</span>`;
                                    case 'C': return `<span class="text-yellow-400">${match}</span>`;
                                    default: return match;
                                }
                            });
                            
                            if (lineNum.trim().match(/^\d+$/)) {
                                // Forward strand line
                                html += `<div class="mb-1 w-full break-all flex">
                                    <span class="text-gray-500 text-xs mr-2 flex-shrink-0">5'</span>
                                    <span class="text-gray-400 mr-2 flex-shrink-0 font-mono">${lineNum}:</span>
                                    <span class="font-mono">${coloredSequence}</span>
                                    <span class="text-gray-500 text-xs ml-2 flex-shrink-0">3'</span>
                                </div>`;
                            } else {
                                // Reverse complement line
                                html += `<div class="mb-1 w-full break-all flex text-gray-300">
                                    <span class="text-gray-500 text-xs mr-2 flex-shrink-0">3'</span>
                                    <span class="text-transparent mr-2 flex-shrink-0 font-mono">${lineNum.replace(/\s/g, '0')}:</span>
                                    <span class="font-mono">${coloredSequence}</span>
                                    <span class="text-gray-500 text-xs ml-2 flex-shrink-0">5'</span>
                                </div>`;
                            }
                        } else if (line.trim() === '') {
                            // Empty separator line
                            html += `<div class="mb-1 w-full">${line}</div>`;
                        } else {
                            // Fallback for other lines
                            html += `<div class="mb-1 w-full break-all">${coloredLine}</div>`;
                        }
                    } else {
                        // Regular single strand display
                        if (line.includes('|')) {
                            const [lineNum, sequence] = line.split('|');
                            // Apply coloring to the sequence part
                            let coloredSequence = sequence.replace(/([ATGC])/g, (match) => {
                                switch(match) {
                                    case 'A': return `<span class="text-red-400">${match}</span>`;
                                    case 'T': return `<span class="text-blue-400">${match}</span>`;
                                    case 'G': return `<span class="text-green-400">${match}</span>`;
                                    case 'C': return `<span class="text-yellow-400">${match}</span>`;
                                    default: return match;
                                }
                            });
                            html += `<div class="mb-1 w-full break-all">
                                <span class="text-gray-400 mr-2 font-mono">${lineNum}:</span>
                                <span class="font-mono">${coloredSequence}</span>
                            </div>`;
                        } else {
                            html += `<div class="mb-1 w-full break-all">${coloredLine}</div>`;
                        }
                    }
                }
            });
            
            html += `</div></div>`;

            html += '</div>';
            outputDiv.innerHTML = html;
        }

        // Display DETAILED ANALYSIS in bottom section
        this.displayDetailedAnalysis(result);

        // Update GC content display
        const gcDisplay = document.getElementById('sequencer-output-gc');
        if (gcDisplay && result.results.sequence_info) {
            gcDisplay.textContent = `${result.results.sequence_info.gc_content}%`;
        }

        // Update melting temperature display
        const tmDisplay = document.getElementById('sequencer-output-tm');
        if (tmDisplay && result.results.melting_temp) {
            // Use salt-adjusted Tm as it's most accurate for PCR
            const tm = result.results.melting_temp.salt_adjusted || result.results.melting_temp.basic || result.results.melting_temp.nearest_neighbor;
            tmDisplay.textContent = `${tm}°C`;
        }
    }

    displayDetailedAnalysis(result) {
        // Check if we have any detailed analysis results
        const hasDetailedResults = result.results.restriction_sites || 
                                 result.results.translations || 
                                 result.results.amino_acid_analysis || 
                                 result.results.gc_optimization || 
                                 result.results.codon_optimization || 
                                 result.results.color_scheme_used || 
                                 result.results.amino_acid_highlighted;

        const detailedContainer = document.getElementById('sequencer-detailed-results');
        const detailedOutput = document.getElementById('sequencer-detailed-output');

        if (!hasDetailedResults) {
            // Hide detailed section if no detailed results
            if (detailedContainer) {
                detailedContainer.classList.add('hidden');
            }
            return;
        }

        // Show detailed section
        if (detailedContainer) {
            detailedContainer.classList.remove('hidden');
        }

        if (detailedOutput) {
            let html = '';

            // Add all detailed analysis results
            if (result.results.gc_optimization) {
                html += this.formatGCOptimization(result.results.gc_optimization, result.results.sequence_optimized, result.results.original_sequence_info);
            }

            if (result.results.codon_optimization) {
                html += this.formatCodonOptimization(result.results.codon_optimization);
            }

            if (result.results.restriction_sites) {
                html += this.formatRestrictionSites(result.results.restriction_sites);
            }

            if (result.results.translations) {
                html += this.formatTranslations(result.results.translations);
            }

            if (result.results.amino_acid_analysis) {
                html += this.formatAminoAcidAnalysis(result.results.amino_acid_analysis);
            }

            if (result.results.color_scheme_used) {
                html += this.formatColorSchemeInfo(result.results.color_scheme_used);
            }

            if (result.results.amino_acid_highlighted) {
                html += this.formatAminoAcidHighlighting(result.results.amino_acid_highlighted);
            }

            detailedOutput.innerHTML = html;
        }
    }

    formatRestrictionSites(sites) {
        let html = `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-red-500/30">`;
        html += `<h4 class="text-red-400 font-medium mb-2">✂️ Restriction Enzyme Sites</h4>`;
        
        if (Object.keys(sites).length === 0) {
            html += `<div class="p-3 bg-red-900/20 border border-red-500/30 rounded-lg">`;
            html += `<p class="text-red-400 text-sm">❌ No restriction sites found for selected enzyme(s)</p>`;
            html += `</div>`;
        } else {
            Object.entries(sites).forEach(([enzyme, data]) => {
                if (data.count > 0) {
                    html += `<div class="mb-2 p-2 bg-green-900/20 border border-green-500/30 rounded">`;
                    html += `<span class="text-green-300 font-medium">✅ ${enzyme}</span> `;
                    html += `<span class="text-gray-400">(${data.recognition_site})</span><br>`;
                    html += `<span class="text-white">${data.count} sites at positions: </span>`;
                    html += `<span class="text-cyan-300 font-mono text-sm">${data.positions.join(', ')}</span>`;
                    html += `</div>`;
                } else {
                    html += `<div class="mb-2 p-2 bg-red-900/20 border border-red-500/30 rounded">`;
                    html += `<span class="text-red-300 font-medium">❌ ${enzyme}</span> `;
                    html += `<span class="text-gray-400">(${data.recognition_site})</span><br>`;
                    html += `<span class="text-red-400 text-sm">No ${enzyme} sites found in sequence</span>`;
                    html += `</div>`;
                }
            });
        }
        
        html += `</div>`;
        return html;
    }

    formatTranslations(translations) {
        let html = `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-blue-500/30">`;
        html += `<h4 class="text-blue-400 font-medium mb-2">🔄 Reading Frame Translations</h4>`;
        
        Object.entries(translations).forEach(([frame, data]) => {
            html += `<div class="mb-3 p-2 bg-black/30 rounded">`;
            html += `<h5 class="text-blue-300 font-medium mb-1">Frame ${frame}</h5>`;
            html += `<div class="text-xs font-mono text-gray-300 break-all mb-1">${data.protein}</div>`;
            html += `<div class="text-xs text-gray-500">Length: ${data.length} bp</div>`;
            html += `</div>`;
        });
        
        html += `</div>`;
        return html;
    }

    formatAminoAcidAnalysis(aminoAcidData) {
        let html = `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-emerald-500/30">`;
        html += `<h4 class="text-emerald-400 font-medium mb-2">🧪 Amino Acid Analysis</h4>`;
        
        if (Object.keys(aminoAcidData).length === 0) {
            html += `<p class="text-gray-400 text-sm">No amino acid matches found</p>`;
        } else {
            Object.entries(aminoAcidData).forEach(([aa, data]) => {
                html += `<div class="mb-3 p-2 bg-black/30 rounded">`;
                html += `<h5 class="text-emerald-300 font-medium mb-1">${data.amino_acid} - ${data.total_count} occurrences</h5>`;
                
                // Show codon usage
                if (data.codon_usage && Object.keys(data.codon_usage).length > 0) {
                    html += `<div class="text-xs text-gray-300 mb-2">Codon usage: `;
                    const codonUsage = Object.entries(data.codon_usage).map(([codon, count]) => `${codon}(${count})`);
                    html += codonUsage.join(', ');
                    html += `</div>`;
                }
                
                // Show first few positions
                if (data.positions && data.positions.length > 0) {
                    const firstPositions = data.positions.slice(0, 10).map(pos => Array.isArray(pos) ? pos[0] : pos);
                    html += `<div class="text-xs text-cyan-300">Positions: ${firstPositions.join(', ')}`;
                    if (data.positions.length > 10) {
                        html += ` (+${data.positions.length - 10} more)`;
                    }
                    html += `</div>`;
                }
                html += `</div>`;
            });
        }
        
        html += `</div>`;
        return html;
    }

    formatGCOptimization(gcOptData, sequenceOptimized, originalInfo) {
        let html = `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-indigo-500/30">`;
        html += `<h4 class="text-indigo-400 font-medium mb-2">⚡ GC Content Optimization</h4>`;
        
        if (sequenceOptimized && gcOptData.optimized_sequence && gcOptData.changes_made > 0) {
            // Show optimization results
            html += `<div class="mb-3 p-3 bg-green-900/20 border border-green-500/30 rounded-lg">`;
            html += `<h5 class="text-green-300 font-medium mb-2">🎯 Sequence Successfully Optimized!</h5>`;
            
            html += `<div class="grid grid-cols-2 md:grid-cols-4 gap-3 text-sm mb-2">`;
            html += `<div><span class="text-gray-400">Original GC:</span> <span class="text-red-300 font-mono">${originalInfo.gc_content.toFixed(1)}%</span></div>`;
            html += `<div><span class="text-gray-400">Optimized GC:</span> <span class="text-green-300 font-mono">${gcOptData.optimized_gc.toFixed(1)}%</span></div>`;
            html += `<div><span class="text-gray-400">Target:</span> <span class="text-white font-mono">${gcOptData.target_range}</span></div>`;
            html += `<div><span class="text-gray-400">Changes:</span> <span class="text-cyan-300 font-mono">${gcOptData.changes_made} codons</span></div>`;
            html += `</div>`;
            
            const improvement = gcOptData.improvement > 0 ? 'improved' : 'adjusted';
            const improvementColor = gcOptData.improvement > 0 ? 'text-green-400' : 'text-yellow-400';
            html += `<p class="text-sm ${improvementColor}">✅ GC content ${improvement} by ${Math.abs(gcOptData.improvement).toFixed(1)} percentage points</p>`;
            html += `<p class="text-xs text-gray-400 mt-1">The sequence above shows the optimized version with synonymous codons</p>`;
            html += `</div>`;
            
        } else if (gcOptData.message && gcOptData.changes_made === 0) {
            // Show "already optimal" message
            html += `<div class="mb-3 p-3 bg-blue-900/20 border border-blue-500/30 rounded-lg">`;
            html += `<h5 class="text-blue-300 font-medium mb-2">✅ GC Content Already Optimal!</h5>`;
            html += `<div class="grid grid-cols-2 md:grid-cols-3 gap-3 text-sm mb-2">`;
            html += `<div><span class="text-gray-400">Current GC:</span> <span class="text-green-300 font-mono">${gcOptData.current_gc.toFixed(1)}%</span></div>`;
            html += `<div><span class="text-gray-400">Target:</span> <span class="text-white font-mono">${gcOptData.target_range}</span></div>`;
            html += `<div><span class="text-gray-400">Status:</span> <span class="text-green-300">In range</span></div>`;
            html += `</div>`;
            html += `<p class="text-sm text-blue-300">✅ ${gcOptData.message}</p>`;
            html += `</div>`;
            
        } else if (gcOptData.error) {
            // Show error message
            html += `<div class="mb-3 p-3 bg-red-900/20 border border-red-500/30 rounded-lg">`;
            html += `<p class="text-red-400 text-sm">❌ ${gcOptData.error}</p>`;
            html += `</div>`;
            
        } else {
            // Show analysis only (no optimization performed)
            html += `<div class="grid grid-cols-2 gap-4 mb-3">`;
            html += `<div class="text-sm">`;
            html += `<span class="text-gray-400">Current GC:</span> `;
            html += `<span class="text-white font-mono">${gcOptData.current_gc}%</span>`;
            html += `</div>`;
            html += `<div class="text-sm">`;
            html += `<span class="text-gray-400">Target Range:</span> `;
            html += `<span class="text-white font-mono">${gcOptData.target_range}</span>`;
            html += `</div>`;
            html += `</div>`;
            
            const statusColor = gcOptData.optimization_needed ? 'text-yellow-400' : 'text-green-400';
            const statusIcon = gcOptData.optimization_needed ? '⚠️' : '✅';
            const statusText = gcOptData.optimization_needed ? 'Optimization recommended' : 'GC content is optimal';
            html += `<p class="text-sm ${statusColor} mb-2">${statusIcon} ${statusText}</p>`;
        }
        
        html += `</div>`;
        return html;
    }

    formatCodonOptimization(codonOptData) {
        let html = `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-purple-500/30">`;
        html += `<h4 class="text-purple-400 font-medium mb-2">🧬 Codon Usage Optimization</h4>`;
        
        if (codonOptData.error || !codonOptData.success) {
            html += `<p class="text-red-400 text-sm">${codonOptData.error || 'Codon optimization failed'}</p>`;
        } else {
            // Display optimization results using NEW format
            html += `<div class="grid grid-cols-2 md:grid-cols-4 gap-3 mb-3 text-sm">`;
            html += `<div><span class="text-gray-400">Organism:</span> <span class="text-white">${codonOptData.organism}</span></div>`;
            html += `<div><span class="text-gray-400">Optimized:</span> <span class="text-green-400 font-mono">${codonOptData.optimization_percentage || 0}%</span></div>`;
            html += `<div><span class="text-gray-400">Changes:</span> <span class="text-cyan-400 font-mono">${codonOptData.changes_made || 0}</span></div>`;
            html += `<div><span class="text-gray-400">Total Codons:</span> <span class="text-white font-mono">${codonOptData.total_codons || 0}</span></div>`;
            html += `</div>`;
            
            // Show optimization message
            if (codonOptData.message) {
                html += `<div class="mb-3 p-2 bg-green-900/30 rounded border border-green-500/30">`;
                html += `<div class="text-green-400 font-medium text-sm">✅ ${codonOptData.message}</div>`;
                html += `</div>`;
            }
            
            // Show protein verification
            if (codonOptData.protein_identical !== undefined) {
                const proteinStatus = codonOptData.protein_identical ? 
                    '<span class="text-green-400">✅ Preserved</span>' : 
                    '<span class="text-red-400">❌ Changed</span>';
                html += `<div class="mb-3 text-sm">`;
                html += `<span class="text-gray-400">Protein Sequence:</span> ${proteinStatus}`;
                html += `</div>`;
            }
            
            // Show GC content change if available
            if (codonOptData.gc_change !== undefined) {
                const gcChangeColor = codonOptData.gc_change > 0 ? 'text-green-400' : 
                                    codonOptData.gc_change < 0 ? 'text-red-400' : 'text-gray-400';
                html += `<div class="mb-3 text-sm">`;
                html += `<span class="text-gray-400">GC Change:</span> `;
                html += `<span class="${gcChangeColor} font-mono">${codonOptData.gc_change > 0 ? '+' : ''}${codonOptData.gc_change}%</span>`;
                html += `<span class="text-gray-400"> (${codonOptData.original_gc}% → ${codonOptData.optimized_gc}%)</span>`;
                html += `</div>`;
            }
            
            // Show sample codon changes
            if (codonOptData.codon_changes && codonOptData.codon_changes.length > 0) {
                html += `<div class="mb-3">`;
                html += `<div class="text-gray-400 text-sm mb-2">Sample Codon Changes:</div>`;
                html += `<div class="space-y-1">`;
                
                codonOptData.codon_changes.slice(0, 5).forEach(change => {
                    html += `<div class="text-xs bg-black/30 rounded p-2">`;
                    html += `<span class="text-cyan-400">${change.amino_acid}</span>: `;
                    html += `<span class="text-red-400 font-mono">${change.original_codon}</span> → `;
                    html += `<span class="text-green-400 font-mono">${change.optimized_codon}</span>`;
                    html += `<span class="text-gray-400 ml-2">(${change.improvement})</span>`;
                    html += `</div>`;
                });
                
                if (codonOptData.codon_changes.length > 5) {
                    html += `<div class="text-xs text-gray-400 mt-1">... and ${codonOptData.codon_changes.length - 5} more changes</div>`;
                }
                
                html += `</div></div>`;
            }
        }
        
        html += `</div>`;
        return html;
    }

    formatColorSchemeInfo(colorScheme) {
        let html = `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-yellow-500/30">`;
        html += `<h4 class="text-yellow-400 font-medium mb-2">🎨 Color Scheme Applied</h4>`;
        html += `<p class="text-gray-300 text-sm">Using <span class="text-white font-medium">${colorScheme}</span> color scheme for nucleotide visualization</p>`;
        html += `</div>`;
        return html;
    }

    formatAminoAcidHighlighting(targetAA) {
        let html = `<div class="mb-4 p-3 bg-slate-800/50 rounded-lg border border-yellow-300/30">`;
        html += `<h4 class="text-yellow-300 font-medium mb-2">🔍 Amino Acid Highlighting</h4>`;
        
        let aaName = '';
        switch(targetAA) {
            case 'M': aaName = 'Methionine (Start codon)'; break;
            case 'START': aaName = 'Start codons (ATG)'; break;
            case 'STOP': aaName = 'Stop codons (TAA, TAG, TGA)'; break;
            case 'ALL': aaName = 'All amino acids'; break;
            case 'A': aaName = 'Alanine'; break;
            case 'R': aaName = 'Arginine'; break;
            case 'N': aaName = 'Asparagine'; break;
            case 'D': aaName = 'Aspartic acid'; break;
            case 'C': aaName = 'Cysteine'; break;
            case 'E': aaName = 'Glutamic acid'; break;
            case 'Q': aaName = 'Glutamine'; break;
            case 'G': aaName = 'Glycine'; break;
            case 'H': aaName = 'Histidine'; break;
            case 'I': aaName = 'Isoleucine'; break;
            case 'L': aaName = 'Leucine'; break;
            case 'K': aaName = 'Lysine'; break;
            case 'F': aaName = 'Phenylalanine'; break;
            case 'P': aaName = 'Proline'; break;
            case 'S': aaName = 'Serine'; break;
            case 'T': aaName = 'Threonine'; break;
            case 'W': aaName = 'Tryptophan'; break;
            case 'Y': aaName = 'Tyrosine'; break;
            case 'V': aaName = 'Valine'; break;
            default: aaName = targetAA;
        }
        
        html += `<p class="text-gray-300 text-sm">Codons coding for <span class="text-white font-medium">${aaName}</span> are highlighted with `;
        html += `<span class="bg-yellow-300 text-black font-bold px-1 rounded">yellow background</span> in the sequence above</p>`;
        html += `</div>`;
        return html;
    }

    hideResults() {
        const emptyState = document.getElementById('sequencer-empty-state');
        const resultsContainer = document.getElementById('sequencer-results');
        const detailedContainer = document.getElementById('sequencer-detailed-results');
        
        if (emptyState) {
            emptyState.style.display = 'block';
        }
        
        if (resultsContainer) {
            resultsContainer.classList.add('hidden');
            resultsContainer.style.display = 'none';
        }
        
        // Also hide detailed results section
        if (detailedContainer) {
            detailedContainer.classList.add('hidden');
        }
    }

    setLoadingState(isLoading) {
        const analyzeBtn = document.getElementById('sequencer-analyze-btn');
        if (analyzeBtn) {
            if (isLoading) {
                analyzeBtn.disabled = true;
                analyzeBtn.innerHTML = '🔄 Analyzing...';
                analyzeBtn.classList.add('opacity-50');
            } else {
                analyzeBtn.disabled = false;
                analyzeBtn.innerHTML = '🧬 Analyze Sequence';
                analyzeBtn.classList.remove('opacity-50');
            }
        }
    }

    async copyResults() {
        if (!this.currentSequence) {
            this.showNotification('No sequence to copy', 'error');
            return;
        }

        try {
            await navigator.clipboard.writeText(this.currentSequence);
            this.showNotification('Sequence copied to clipboard!', 'success');
        } catch (error) {
            console.error('Copy failed:', error);
            this.showNotification('Failed to copy to clipboard', 'error');
        }
    }

    // Utility methods for notifications and overlays
    showNotification(message, type = 'info') {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `fixed top-4 right-4 z-50 p-4 rounded-lg shadow-lg transition-all duration-300 ${
            type === 'success' ? 'bg-green-600 text-white' :
            type === 'error' ? 'bg-red-600 text-white' :
            'bg-blue-600 text-white'
        }`;
        notification.textContent = message;
        
        document.body.appendChild(notification);
        
        // Auto remove after 3 seconds
        setTimeout(() => {
            notification.remove();
        }, 3000);
    }

    showLoadingOverlay(message) {
        const overlay = document.createElement('div');
        overlay.id = 'sequencer-loading-overlay';
        overlay.className = 'fixed inset-0 bg-black/80 backdrop-blur-sm z-50 flex items-center justify-center';
        overlay.innerHTML = `
            <div class="bg-slate-800 border border-cyan-500/30 rounded-xl p-6 text-center">
                <div class="animate-spin w-8 h-8 border-2 border-cyan-500 border-t-transparent rounded-full mx-auto mb-4"></div>
                <p class="text-white">${message}</p>
            </div>
        `;
        document.body.appendChild(overlay);
    }

    hideLoadingOverlay() {
        const overlay = document.getElementById('sequencer-loading-overlay');
        if (overlay) {
            overlay.remove();
        }
    }

    showLargeFileWarning(file, fileType, maxRecommendedSize) {
        return new Promise((resolve) => {
            const overlay = document.createElement('div');
            overlay.className = 'fixed inset-0 bg-black/80 backdrop-blur-sm z-50 flex items-center justify-center p-4';
            
            const modal = document.createElement('div');
            modal.className = 'bg-slate-800 border-2 border-red-500 rounded-xl p-6 max-w-2xl w-full mx-auto shadow-2xl';
            
            const fileSize = this.formatFileSize(file.size);
            const recommendedSize = this.formatFileSize(maxRecommendedSize);
            
            modal.innerHTML = `
                <div class="text-center">
                    <div class="text-6xl mb-4">⚠️</div>
                    <h3 class="text-2xl font-bold text-red-400 mb-4">Large File Warning</h3>
                    <div class="text-left space-y-3 mb-6">
                        <p class="text-gray-300">You're about to load a large ${fileType} file:</p>
                        <div class="bg-slate-900/50 rounded-lg p-3 font-mono text-sm">
                            <div class="text-cyan-400">File: ${file.name}</div>
                            <div class="text-yellow-400">Size: ${fileSize}</div>
                            <div class="text-gray-400">Recommended limit: ${recommendedSize}</div>
                        </div>
                        <div class="text-yellow-300 text-sm">
                            <p>⚡ Large files may take significant time to process and could impact performance.</p>
                            <p>💾 Consider using smaller files or more powerful hardware for optimal experience.</p>
                        </div>
                    </div>
                    <div class="flex gap-3 justify-center">
                        <button id="cancel-large-file" class="px-6 py-2 bg-gray-600 hover:bg-gray-700 text-white rounded-lg transition-colors">
                            Cancel
                        </button>
                        <button id="proceed-large-file" class="px-6 py-2 bg-red-600 hover:bg-red-700 text-white rounded-lg transition-colors">
                            Proceed Anyway
                        </button>
                    </div>
                </div>
            `;
            
            overlay.appendChild(modal);
            document.body.appendChild(overlay);
            
            // Event listeners
            modal.querySelector('#cancel-large-file').addEventListener('click', () => {
                overlay.remove();
                resolve(false);
            });
            
            modal.querySelector('#proceed-large-file').addEventListener('click', () => {
                overlay.remove();
                resolve(true);
            });
        });
    }

    formatFileSize(bytes) {
        if (bytes === 0) return '0 Bytes';
        const k = 1024;
        const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));
        return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
    }

    toggleEditLockMode() {
        const editLockBtn = document.getElementById('edit-lock-mode');
        const sequencerOutput = document.getElementById('sequencer-output');
        
        if (!editLockBtn || !sequencerOutput) {
            this.showNotification('Edit mode not available', 'error');
            return;
        }

        // Check if we're currently in edit mode
        const isEditMode = sequencerOutput.contentEditable === 'true';
        
        if (isEditMode) {
            // Exit edit mode
            sequencerOutput.contentEditable = 'false';
            sequencerOutput.style.border = '';
            sequencerOutput.style.backgroundColor = '';
            editLockBtn.innerHTML = '🔒 Edit Lock';
            editLockBtn.className = 'px-3 py-1 bg-orange-600/20 border border-orange-500/30 rounded text-orange-300 hover:bg-orange-600/30 transition-all text-sm';
            
            // Update the current sequence with edited content
            this.updateSequenceFromEdit();
            
            this.showNotification('Edit mode disabled - sequence locked', 'info');
        } else {
            // Enter edit mode
            if (!this.currentSequence) {
                this.showNotification('No sequence to edit - analyze a sequence first', 'error');
                return;
            }
            
            sequencerOutput.contentEditable = 'true';
            sequencerOutput.style.border = '2px solid #f59e0b';
            sequencerOutput.style.backgroundColor = 'rgba(245, 158, 11, 0.1)';
            editLockBtn.innerHTML = '🔓 Exit Edit';
            editLockBtn.className = 'px-3 py-1 bg-green-600/20 border border-green-500/30 rounded text-green-300 hover:bg-green-600/30 transition-all text-sm';
            
            this.showNotification('Edit mode enabled - you can now freely edit the DNA sequence!', 'success');
            
            // Focus on the editable area
            sequencerOutput.focus();
        }
    }

    updateSequenceFromEdit() {
        const sequencerOutput = document.getElementById('sequencer-output');
        if (!sequencerOutput) return;
        
        // Extract plain text from the edited content, removing HTML tags
        let editedText = sequencerOutput.innerText || sequencerOutput.textContent || '';
        
        // Clean up the text - remove non-DNA characters and normalize
        editedText = editedText.replace(/[^ATGCN\s]/gi, '').replace(/\s+/g, '').toUpperCase();
        
        if (editedText.length > 0) {
            this.currentSequence = editedText;
            
            // Update the GC content display
            this.updateGCContentFromSequence(editedText);
            
            this.showNotification(`Sequence updated! Length: ${editedText.length} bp`, 'success');
        }
    }

    updateGCContentFromSequence(sequence) {
        const gcOutput = document.getElementById('sequencer-output-gc');
        if (!gcOutput || !sequence) return;
        
        // Calculate GC content
        const gcCount = (sequence.match(/[GC]/g) || []).length;
        const gcContent = ((gcCount / sequence.length) * 100).toFixed(1);
        
        gcOutput.textContent = `${gcContent}%`;
    }
}

// Initialize DNA Sequencer when the page loads
document.addEventListener('DOMContentLoaded', () => {
    // Only initialize if we're on the DNA sequencer page
    if (document.getElementById('dna-sequencer-page')) {
        window.dnaSequencer = new DNASequencerPython();
    }
});

// Also initialize when navigating to the page
document.addEventListener('pageChanged', (e) => {
    if (e.detail.page === 'dna-sequencer') {
        if (!window.dnaSequencer) {
            window.dnaSequencer = new DNASequencerPython();
        }
    }
});