class CodonFinderPython {
    constructor() {
        this.selectedAminoAcid = null;
        this.currentSequence = '';
        this.codonData = null;
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUpload();
        this.setupProteinSearch();
    }

    bindEvents() {
        // Analyze button
        const analyzeBtn = document.getElementById('analyze-codons-btn');
        if (analyzeBtn) {
            analyzeBtn.addEventListener('click', () => this.performAnalysis());
        }

        // Clear button
        const clearBtn = document.getElementById('codon-clear-btn');
        if (clearBtn) {
            clearBtn.addEventListener('click', () => this.clearInput());
        }

        // DNA input textarea
        const dnaInput = document.getElementById('codon-dna-input');
        if (dnaInput) {
            dnaInput.addEventListener('input', () => this.updateSequenceLength());
        }

        // Protein selection
        const proteinItems = document.querySelectorAll('.protein-item');
        proteinItems.forEach(item => {
            item.addEventListener('click', () => this.selectAminoAcid(item));
        });
    }

    setupFileUpload() {
        const fileInput = document.getElementById('codon-file-upload');
        const dropZone = fileInput?.parentElement;

        if (fileInput && dropZone) {
            fileInput.addEventListener('change', (e) => {
                if (e.target.files.length > 0) {
                    this.handleFileUpload(e.target.files[0]);
                }
            });

            // Drag and drop
            dropZone.addEventListener('dragover', (e) => {
                e.preventDefault();
                dropZone.classList.add('border-emerald-500/50');
            });

            dropZone.addEventListener('dragleave', (e) => {
                e.preventDefault();
                dropZone.classList.remove('border-emerald-500/50');
            });

            dropZone.addEventListener('drop', (e) => {
                e.preventDefault();
                dropZone.classList.remove('border-emerald-500/50');
                
                if (e.dataTransfer.files.length > 0) {
                    this.handleFileUpload(e.dataTransfer.files[0]);
                }
            });
        }
    }

    setupProteinSearch() {
        const searchInput = document.getElementById('protein-search');
        if (searchInput) {
            searchInput.addEventListener('input', (e) => {
                this.filterProteins(e.target.value);
            });
        }
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

            const fileContent = await this.readFile(file);
            const fileName = file.name.toLowerCase();
            
            let fileFormat = 'txt';
            if (fileName.endsWith('.fasta') || fileName.endsWith('.fa')) {
                fileFormat = 'fasta';
            }

            // Parse file content and extract sequence
            let sequence = '';
            if (fileFormat === 'fasta') {
                const lines = fileContent.split('\n');
                const sequenceLines = lines.filter(line => !line.startsWith('>'));
                sequence = sequenceLines.join('').replace(/\s/g, '');
            } else {
                sequence = fileContent.replace(/\s/g, '');
            }

            // Put sequence in textarea
            const sequenceInput = document.getElementById('codon-dna-input');
            if (sequenceInput) {
                sequenceInput.value = sequence;
                this.updateSequenceLength();
            }

            this.hideLoadingOverlay();
            this.showNotification(`File "${file.name}" loaded successfully`, 'success');

        } catch (error) {
            console.error('File upload error:', error);
            this.hideLoadingOverlay();
            this.showNotification(`Failed to load file: ${error.message}`, 'error');
        }
    }

    readFile(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = (e) => resolve(e.target.result);
            reader.onerror = (e) => reject(new Error('Failed to read file'));
            reader.readAsText(file);
        });
    }

    updateSequenceLength() {
        const sequenceInput = document.getElementById('codon-dna-input');
        const lengthDisplay = document.getElementById('codon-sequence-length');
        
        if (sequenceInput && lengthDisplay) {
            const sequence = sequenceInput.value.replace(/\s/g, '');
            lengthDisplay.textContent = `Length: ${sequence.length} bp`;
            this.currentSequence = sequence;
        }
    }

    selectAminoAcid(item) {
        // Remove previous selection
        document.querySelectorAll('.protein-item').forEach(el => {
            el.classList.remove('ring-2', 'ring-purple-500', 'bg-purple-600/20');
        });

        // Add selection to clicked item
        item.classList.add('ring-2', 'ring-purple-500', 'bg-purple-600/20');

        // Get amino acid info
        const protein = item.dataset.protein;
        this.selectedAminoAcid = this.getAminoAcidCode(protein);

        // Update selected amino acid display
        this.updateSelectedAminoDisplay(protein);

        // Enable analyze button if we have sequence
        this.updateAnalyzeButton();

        this.showNotification(`Selected: ${this.getAminoAcidName(protein)}`, 'info');
    }

    getAminoAcidCode(protein) {
        const mapping = {
            'alanine': 'A',
            'arginine': 'R',
            'asparagine': 'N',
            'aspartic_acid': 'D',
            'cysteine': 'C',
            'glutamic_acid': 'E',
            'glutamine': 'Q',
            'glycine': 'G',
            'histidine': 'H',
            'isoleucine': 'I',
            'leucine': 'L',
            'lysine': 'K',
            'methionine': 'M',
            'phenylalanine': 'F',
            'proline': 'P',
            'serine': 'S',
            'threonine': 'T',
            'tryptophan': 'W',
            'tyrosine': 'Y',
            'valine': 'V',
            'stop': '*'
        };
        return mapping[protein] || protein.toUpperCase();
    }

    getAminoAcidName(protein) {
        const mapping = {
            'alanine': 'Alanine (Ala, A)',
            'arginine': 'Arginine (Arg, R)',
            'asparagine': 'Asparagine (Asn, N)',
            'aspartic_acid': 'Aspartic Acid (Asp, D)',
            'cysteine': 'Cysteine (Cys, C)',
            'glutamic_acid': 'Glutamic Acid (Glu, E)',
            'glutamine': 'Glutamine (Gln, Q)',
            'glycine': 'Glycine (Gly, G)',
            'histidine': 'Histidine (His, H)',
            'isoleucine': 'Isoleucine (Ile, I)',
            'leucine': 'Leucine (Leu, L)',
            'lysine': 'Lysine (Lys, K)',
            'methionine': 'Methionine (Met, M) - START',
            'phenylalanine': 'Phenylalanine (Phe, F)',
            'proline': 'Proline (Pro, P)',
            'serine': 'Serine (Ser, S)',
            'threonine': 'Threonine (Thr, T)',
            'tryptophan': 'Tryptophan (Trp, W)',
            'tyrosine': 'Tyrosine (Tyr, Y)',
            'valine': 'Valine (Val, V)',
            'stop': 'Stop Codons - STOP'
        };
        return mapping[protein] || protein;
    }

    updateSelectedAminoDisplay(protein) {
        const selectedInfo = document.getElementById('selected-amino-info');
        const selectedName = document.getElementById('selected-amino-name');
        
        if (selectedInfo && selectedName) {
            selectedInfo.classList.remove('hidden');
            selectedName.textContent = this.getAminoAcidName(protein);
        }
    }

    updateAnalyzeButton() {
        const analyzeBtn = document.getElementById('analyze-codons-btn');
        if (analyzeBtn) {
            const hasSequence = this.currentSequence.length > 0;
            const hasSelection = this.selectedAminoAcid !== null;
            
            if (hasSequence && hasSelection) {
                analyzeBtn.disabled = false;
                analyzeBtn.classList.remove('opacity-50');
                analyzeBtn.innerHTML = '🔍 Analyze Codons';
            } else {
                analyzeBtn.disabled = true;
                analyzeBtn.classList.add('opacity-50');
                if (!hasSequence && !hasSelection) {
                    analyzeBtn.innerHTML = '🔍 Enter sequence & select amino acid';
                } else if (!hasSequence) {
                    analyzeBtn.innerHTML = '🔍 Enter DNA sequence';
                } else {
                    analyzeBtn.innerHTML = '🔍 Select amino acid';
                }
            }
        }
    }

    filterProteins(searchTerm) {
        const proteinItems = document.querySelectorAll('.protein-item');
        const term = searchTerm.toLowerCase();

        proteinItems.forEach(item => {
            const text = item.textContent.toLowerCase();
            if (text.includes(term)) {
                item.style.display = 'block';
            } else {
                item.style.display = 'none';
            }
        });
    }

    clearInput() {
        const sequenceInput = document.getElementById('codon-dna-input');
        if (sequenceInput) {
            sequenceInput.value = '';
            this.currentSequence = '';
            this.updateSequenceLength();
            this.updateAnalyzeButton();
        }

        // Clear results
        this.hideResults();
        this.showNotification('Input cleared', 'info');
    }

    async performAnalysis() {
        try {
            if (!this.currentSequence || !this.selectedAminoAcid) {
                this.showNotification('Please enter a DNA sequence and select an amino acid', 'error');
                return;
            }

            // Show processing state
            this.setProcessingState(true);

            // Call Python script
            const result = await this.callPythonScript(this.currentSequence, this.selectedAminoAcid);

            if (result.success) {
                this.codonData = result;
                this.displayResults(result);
                this.showNotification(`Found ${result.statistics.target_count} ${result.amino_acid_info.symbol} codons!`, 'success');
            } else {
                throw new Error(result.error || 'Codon analysis failed');
            }

        } catch (error) {
            console.error('Analysis error:', error);
            this.showNotification(`Analysis failed: ${error.message}`, 'error');
        } finally {
            this.setProcessingState(false);
        }
    }

    async callPythonScript(sequence, aminoAcid) {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Python script
            const scriptPath = path.join(__dirname, 'assets', 'codon-finder.py');
            const args = [
                scriptPath,
                '--sequence', sequence,
                '--amino-acid', aminoAcid
            ];

            console.log('🔍 Calling Codon Finder:', scriptPath);

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

            pythonProcess.on('close', (code) => {
                console.log(`🔍 Python process exited with code: ${code}`);
                
                if (code === 0) {
                    try {
                        const result = JSON.parse(stdout);
                        console.log('✅ Codon Finder result:', result);
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
        // Hide empty state
        const emptyState = document.getElementById('codon-empty-state');
        if (emptyState) {
            emptyState.style.display = 'none';
        }

        // Show results
        const resultsDiv = document.getElementById('codon-results');
        if (resultsDiv) {
            resultsDiv.classList.remove('hidden');
            resultsDiv.style.display = 'block';
        }

        // Update highlighted sequence
        this.displayHighlightedSequence(result);

        // Update statistics
        this.updateStatistics(result.statistics);

        // Update codon details
        this.displayCodonDetails(result.codon_details);

        // Setup hover tooltips
        this.setupCodonTooltips();
    }

    displayHighlightedSequence(result) {
        const sequenceDiv = document.getElementById('highlighted-sequence');
        const codonCount = document.getElementById('codon-count');
        const sequenceLength = document.getElementById('sequence-length-display');

        if (sequenceDiv) {
            sequenceDiv.innerHTML = result.highlighted_sequence;
        }

        if (codonCount) {
            const count = result.statistics.target_count;
            codonCount.textContent = `${count} ${result.amino_acid_info.symbol} codon${count !== 1 ? 's' : ''} found`;
        }

        if (sequenceLength) {
            sequenceLength.textContent = `${result.statistics.sequence_length} bp`;
        }
    }

    updateStatistics(stats) {
        const totalCodons = document.getElementById('total-codons-stat');
        const targetCodons = document.getElementById('target-codons-stat');
        const frequency = document.getElementById('frequency-stat');

        if (totalCodons) totalCodons.textContent = stats.total_codons;
        if (targetCodons) targetCodons.textContent = stats.target_count;
        if (frequency) frequency.textContent = `${stats.frequency}%`;
    }

    displayCodonDetails(details) {
        const detailsDiv = document.getElementById('codon-details');
        if (!detailsDiv) return;

        if (details.length === 0) {
            detailsDiv.innerHTML = '<p class="text-gray-400 text-sm">No target codons found in sequence</p>';
            return;
        }

        let html = '<div class="space-y-2">';
        
        details.forEach((detail, index) => {
            html += `
                <div class="flex items-center justify-between bg-slate-900/30 rounded-lg p-3 border border-slate-600/20">
                    <div class="flex items-center gap-3">
                        <span class="w-6 h-6 bg-purple-600/20 border border-purple-500/30 rounded-full flex items-center justify-center text-purple-400 text-xs font-bold">${detail.index}</span>
                        <div>
                            <h6 class="text-white font-mono font-medium">${detail.codon}</h6>
                            <p class="text-gray-400 text-xs">${detail.amino_acid_name} (${detail.amino_acid_symbol})</p>
                        </div>
                    </div>
                    <div class="text-right">
                        <p class="text-gray-300 text-sm">Position ${detail.position}</p>
                        <p class="text-gray-500 text-xs capitalize">${detail.amino_acid_type}</p>
                    </div>
                </div>
            `;
        });
        
        html += '</div>';
        detailsDiv.innerHTML = html;
    }

    setupCodonTooltips() {
        const codonElements = document.querySelectorAll('[data-codon]');
        
        codonElements.forEach(element => {
            element.addEventListener('mouseenter', (e) => {
                this.showCodonTooltip(e);
            });
            
            element.addEventListener('mouseleave', () => {
                this.hideCodonTooltip();
            });
        });
    }

    showCodonTooltip(event) {
        const element = event.target;
        const codon = element.dataset.codon;
        const amino = element.dataset.amino;
        const position = element.dataset.position;

        // Remove existing tooltip
        this.hideCodonTooltip();

        // Create tooltip
        const tooltip = document.createElement('div');
        tooltip.id = 'codon-tooltip';
        tooltip.className = 'fixed bg-slate-800 border border-slate-600 rounded-lg p-3 text-sm z-50 shadow-xl';
        
        const aminoInfo = this.getAminoAcidInfo(amino);
        
        tooltip.innerHTML = `
            <div class="space-y-2">
                <div class="flex items-center gap-2">
                    <span class="font-mono font-bold text-white">${codon}</span>
                    <span class="text-gray-400">→</span>
                    <span class="text-cyan-400 font-medium">${aminoInfo.name}</span>
                </div>
                <div class="text-xs text-gray-400">
                    <p>Symbol: ${aminoInfo.symbol}</p>
                    <p>Position: ${parseInt(position) + 1}</p>
                    <p>Type: ${aminoInfo.type}</p>
                </div>
            </div>
        `;

        document.body.appendChild(tooltip);

        // Position tooltip
        const rect = element.getBoundingClientRect();
        tooltip.style.left = `${rect.left + rect.width / 2 - tooltip.offsetWidth / 2}px`;
        tooltip.style.top = `${rect.top - tooltip.offsetHeight - 8}px`;

        // Ensure tooltip stays on screen
        const tooltipRect = tooltip.getBoundingClientRect();
        if (tooltipRect.left < 0) {
            tooltip.style.left = '8px';
        }
        if (tooltipRect.right > window.innerWidth) {
            tooltip.style.left = `${window.innerWidth - tooltipRect.width - 8}px`;
        }
        if (tooltipRect.top < 0) {
            tooltip.style.top = `${rect.bottom + 8}px`;
        }
    }

    hideCodonTooltip() {
        const tooltip = document.getElementById('codon-tooltip');
        if (tooltip) {
            document.body.removeChild(tooltip);
        }
    }

    getAminoAcidInfo(code) {
        const info = {
            'A': { name: 'Alanine', symbol: 'Ala', type: 'nonpolar' },
            'R': { name: 'Arginine', symbol: 'Arg', type: 'basic' },
            'N': { name: 'Asparagine', symbol: 'Asn', type: 'polar' },
            'D': { name: 'Aspartic Acid', symbol: 'Asp', type: 'acidic' },
            'C': { name: 'Cysteine', symbol: 'Cys', type: 'polar' },
            'E': { name: 'Glutamic Acid', symbol: 'Glu', type: 'acidic' },
            'Q': { name: 'Glutamine', symbol: 'Gln', type: 'polar' },
            'G': { name: 'Glycine', symbol: 'Gly', type: 'nonpolar' },
            'H': { name: 'Histidine', symbol: 'His', type: 'basic' },
            'I': { name: 'Isoleucine', symbol: 'Ile', type: 'nonpolar' },
            'L': { name: 'Leucine', symbol: 'Leu', type: 'nonpolar' },
            'K': { name: 'Lysine', symbol: 'Lys', type: 'basic' },
            'M': { name: 'Methionine', symbol: 'Met', type: 'nonpolar' },
            'F': { name: 'Phenylalanine', symbol: 'Phe', type: 'nonpolar' },
            'P': { name: 'Proline', symbol: 'Pro', type: 'nonpolar' },
            'S': { name: 'Serine', symbol: 'Ser', type: 'polar' },
            'T': { name: 'Threonine', symbol: 'Thr', type: 'polar' },
            'W': { name: 'Tryptophan', symbol: 'Trp', type: 'nonpolar' },
            'Y': { name: 'Tyrosine', symbol: 'Tyr', type: 'polar' },
            'V': { name: 'Valine', symbol: 'Val', type: 'nonpolar' },
            '*': { name: 'Stop Codon', symbol: 'Stop', type: 'stop' }
        };
        
        return info[code] || { name: 'Unknown', symbol: 'X', type: 'unknown' };
    }

    hideResults() {
        const emptyState = document.getElementById('codon-empty-state');
        const resultsDiv = document.getElementById('codon-results');
        
        if (emptyState) {
            emptyState.style.display = 'block';
        }
        if (resultsDiv) {
            resultsDiv.classList.add('hidden');
            resultsDiv.style.display = 'none';
        }
    }

    setProcessingState(processing) {
        const analyzeBtn = document.getElementById('analyze-codons-btn');
        
        if (analyzeBtn) {
            if (processing) {
                analyzeBtn.disabled = true;
                analyzeBtn.innerHTML = '⏳ Analyzing...';
                analyzeBtn.classList.add('opacity-50');
            } else {
                this.updateAnalyzeButton();
            }
        }
    }

    // Large file warning system (copied from other tools)
    showLargeFileWarning(file, fileType, maxRecommendedSize = 5 * 1024 * 1024 * 1024) {
        return new Promise((resolve) => {
            // Create overlay
            const overlay = document.createElement('div');
            overlay.className = 'fixed inset-0 bg-black/80 backdrop-blur-sm z-50 flex items-center justify-center p-4';
            
            // Create warning modal
            const modal = document.createElement('div');
            modal.className = 'bg-slate-800 border-2 border-red-500 rounded-xl p-6 max-w-2xl w-full mx-auto shadow-2xl';
            
            const fileSize = this.formatFileSize(file.size);
            const isAboveRecommended = file.size > maxRecommendedSize; // Above 5GB
            const isVeryLarge = file.size > 7 * 1024 * 1024 * 1024; // Above 7GB
            
            modal.innerHTML = `
                <div class="text-center mb-6">
                    <div class="w-16 h-16 bg-red-600/20 rounded-full flex items-center justify-center mx-auto mb-4">
                        <span class="text-3xl">⚠️</span>
                    </div>
                    <h2 class="text-2xl font-bold text-red-400 mb-2">Large File Warning</h2>
                    <p class="text-gray-300">You're about to upload a large ${fileType} file</p>
                </div>
                
                <div class="bg-red-900/20 border border-red-500/30 rounded-lg p-4 mb-6">
                    <div class="flex items-center gap-3 mb-3">
                        <span class="text-red-400 text-xl">🚨</span>
                        <h3 class="text-red-400 font-semibold">File: ${file.name}</h3>
                    </div>
                    <p class="text-white font-mono text-lg mb-2">Size: ${fileSize}</p>
                    <div class="text-red-300 text-sm space-y-1">
                        <p>• This file is larger than the recommended 100MB limit</p>
                        ${isAboveRecommended ? '<p class="text-orange-400 font-semibold">• File exceeds 5GB recommended limit - requires powerful computer</p>' : ''}
                        <p>• Processing may take several minutes to complete</p>
                        <p>• Your computer may become slow or unresponsive</p>
                        ${isVeryLarge ? '<p class="text-red-400 font-semibold">• Risk of application crash or memory errors (7GB+)</p>' : ''}
                        <p>• RAM usage could exceed ${Math.ceil(file.size / (1024 * 1024 * 1024) * 2)}GB during processing</p>
                    </div>
                </div>
                
                ${isAboveRecommended ? `
                <div class="bg-orange-900/20 border border-orange-500/30 rounded-lg p-4 mb-6">
                    <h4 class="text-orange-400 font-semibold mb-2">🔥 High-Performance Computing Required:</h4>
                    <div class="text-orange-300 text-sm space-y-1">
                        <p>• Your file is above the 5GB recommended limit</p>
                        <p>• Maximum supported: 10GB (for powerful computers only)</p>
                        <p>• Requires high-end computer with excellent cooling</p>
                        <p>• Consider splitting large genomes into smaller chunks</p>
                    </div>
                </div>
                ` : ''}
                
                <div class="bg-yellow-900/20 border border-yellow-500/30 rounded-lg p-4 mb-6">
                    <h4 class="text-yellow-400 font-semibold mb-2">⚡ System Requirements:</h4>
                    <div class="text-yellow-300 text-sm space-y-1">
                        <p>• Minimum ${Math.ceil(file.size / (1024 * 1024 * 1024) * 3)}GB available RAM recommended</p>
                        <p>• ${isAboveRecommended ? 'High-end CPU (8+ cores recommended)' : 'Modern CPU with good performance'}</p>
                        <p>• Close other applications to free up memory</p>
                        <p>• Ensure stable power supply (processing may take ${isAboveRecommended ? '30-60' : '10-30'} minutes)</p>
                        <p>• Save your work in other applications before proceeding</p>
                        ${isAboveRecommended ? '<p class="text-orange-300 font-semibold">• Consider using a dedicated workstation for files this large</p>' : ''}
                    </div>
                </div>
                
                <div class="bg-slate-900/50 border border-slate-600/30 rounded-lg p-4 mb-6">
                    <h4 class="text-white font-semibold mb-3">Type "CONFIRM" to proceed with large file processing:</h4>
                    <input 
                        type="text" 
                        id="confirm-input" 
                        class="w-full bg-slate-800 border border-slate-600 rounded-lg px-4 py-3 text-white font-mono text-lg focus:border-red-500 focus:outline-none"
                        placeholder="Type CONFIRM here..."
                        autocomplete="off"
                        spellcheck="false"
                    >
                </div>
                
                <div class="flex gap-4">
                    <button 
                        id="cancel-upload" 
                        class="flex-1 px-6 py-3 bg-slate-600 hover:bg-slate-500 border border-slate-500 rounded-lg text-white font-medium transition-all"
                    >
                        ❌ Cancel Upload
                    </button>
                    <button 
                        id="proceed-upload" 
                        class="flex-1 px-6 py-3 bg-red-600/50 border border-red-500 rounded-lg text-red-300 font-medium transition-all disabled:opacity-30 disabled:cursor-not-allowed"
                        disabled
                    >
                        ⚠️ Proceed Anyway
                    </button>
                </div>
            `;
            
            overlay.appendChild(modal);
            document.body.appendChild(overlay);
            
            // Get elements
            const confirmInput = modal.querySelector('#confirm-input');
            const proceedBtn = modal.querySelector('#proceed-upload');
            const cancelBtn = modal.querySelector('#cancel-upload');
            
            // Handle confirmation input
            confirmInput.addEventListener('input', () => {
                const isValid = confirmInput.value.trim().toUpperCase() === 'CONFIRM';
                proceedBtn.disabled = !isValid;
                
                if (isValid) {
                    proceedBtn.classList.remove('bg-red-600/50', 'text-red-300');
                    proceedBtn.classList.add('bg-red-600', 'text-white', 'hover:bg-red-500');
                } else {
                    proceedBtn.classList.add('bg-red-600/50', 'text-red-300');
                    proceedBtn.classList.remove('bg-red-600', 'text-white', 'hover:bg-red-500');
                }
            });
            
            // Handle proceed
            proceedBtn.addEventListener('click', () => {
                if (confirmInput.value.trim().toUpperCase() === 'CONFIRM') {
                    document.body.removeChild(overlay);
                    resolve(true);
                }
            });
            
            // Handle cancel
            cancelBtn.addEventListener('click', () => {
                document.body.removeChild(overlay);
                resolve(false);
            });
            
            // Focus input
            setTimeout(() => confirmInput.focus(), 100);
        });
    }

    showLoadingOverlay(message) {
        // Remove existing overlay
        this.hideLoadingOverlay();
        
        const overlay = document.createElement('div');
        overlay.id = 'loading-overlay';
        overlay.className = 'fixed inset-0 bg-black/70 backdrop-blur-sm z-40 flex items-center justify-center';
        
        overlay.innerHTML = `
            <div class="bg-slate-800 border border-slate-600 rounded-xl p-8 text-center max-w-md mx-4">
                <div class="animate-spin w-12 h-12 border-4 border-blue-500/30 border-t-blue-500 rounded-full mx-auto mb-4"></div>
                <h3 class="text-white text-xl font-semibold mb-2">Processing Large File</h3>
                <p class="text-gray-300 mb-4">${message}</p>
                <div class="bg-slate-900/50 rounded-lg p-3">
                    <p class="text-yellow-400 text-sm">⚠️ Please wait, do not close the application</p>
                    <p class="text-gray-400 text-xs mt-1">This may take several minutes...</p>
                </div>
            </div>
        `;
        
        document.body.appendChild(overlay);
    }

    hideLoadingOverlay() {
        const overlay = document.getElementById('loading-overlay');
        if (overlay) {
            document.body.removeChild(overlay);
        }
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
        notification.className = `fixed top-4 right-4 px-6 py-3 rounded-lg text-white font-medium z-50 transition-all duration-300 transform translate-x-full`;
        
        // Set color based on type
        switch (type) {
            case 'success':
                notification.classList.add('bg-emerald-600');
                break;
            case 'error':
                notification.classList.add('bg-red-600');
                break;
            case 'warning':
                notification.classList.add('bg-yellow-600');
                break;
            default:
                notification.classList.add('bg-blue-600');
        }
        
        notification.textContent = message;
        document.body.appendChild(notification);
        
        // Animate in
        setTimeout(() => {
            notification.classList.remove('translate-x-full');
        }, 100);
        
        // Animate out and remove
        setTimeout(() => {
            notification.classList.add('translate-x-full');
            setTimeout(() => {
                document.body.removeChild(notification);
            }, 300);
        }, 3000);
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    if (document.getElementById('codon-finder-page')) {
        window.codonFinderPython = new CodonFinderPython();
    }
});