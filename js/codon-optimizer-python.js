class CodonOptimizerPython {
    constructor() {
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUpload();
        this.updateSequenceLength();
        this.updateGCContent();
    }

    bindEvents() {
        // Optimize button
        const optimizeBtn = document.getElementById('optimize-btn');
        if (optimizeBtn) {
            optimizeBtn.addEventListener('click', () => this.performOptimization());
        }

        // Clear button
        const clearBtn = document.getElementById('optimizer-clear-btn');
        if (clearBtn) {
            clearBtn.addEventListener('click', () => this.clearInput());
        }

        // Input textarea
        const sequenceInput = document.getElementById('optimizer-input');
        if (sequenceInput) {
            sequenceInput.addEventListener('input', () => {
                this.updateSequenceLength();
                this.updateGCContent();
            });
        }

        // GC optimizer dropdown
        const gcOptimizer = document.getElementById('gc-optimizer');
        if (gcOptimizer) {
            gcOptimizer.addEventListener('change', () => this.handleGCOptimizerChange());
        }

        // Custom GC target slider
        const gcSlider = document.getElementById('gc-target-slider');
        if (gcSlider) {
            gcSlider.addEventListener('input', () => this.updateGCSliderValue());
        }

        // Copy button
        const copyBtn = document.getElementById('copy-optimized');
        if (copyBtn) {
            copyBtn.addEventListener('click', () => this.copyToClipboard('optimized-result'));
        }

        // Initialize UI
        this.handleGCOptimizerChange();
        this.updateGCSliderValue();
    }

    setupFileUpload() {
        const fileInput = document.getElementById('optimizer-file-upload');
        const dropZone = fileInput?.parentElement;

        if (!fileInput || !dropZone) return;

        // File input change
        fileInput.addEventListener('change', (e) => {
            if (e.target.files.length > 0) {
                this.handleFileUpload(e.target.files[0]);
            }
        });

        // Drag and drop
        dropZone.addEventListener('dragover', (e) => {
            e.preventDefault();
            dropZone.classList.add('border-lime-500/50');
        });

        dropZone.addEventListener('dragleave', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-lime-500/50');
        });

        dropZone.addEventListener('drop', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-lime-500/50');
            
            if (e.dataTransfer.files.length > 0) {
                this.handleFileUpload(e.dataTransfer.files[0]);
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
            const sequenceInput = document.getElementById('optimizer-input');
            if (sequenceInput) {
                sequenceInput.value = sequence;
                this.updateSequenceLength();
                this.updateGCContent();
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

    handleGCOptimizerChange() {
        const gcOptimizer = document.getElementById('gc-optimizer');
        const customGCTarget = document.getElementById('custom-gc-target');
        
        if (gcOptimizer && customGCTarget) {
            if (gcOptimizer.value === 'custom') {
                customGCTarget.classList.remove('hidden');
            } else {
                customGCTarget.classList.add('hidden');
            }
        }
    }

    updateGCSliderValue() {
        const slider = document.getElementById('gc-target-slider');
        const valueDisplay = document.getElementById('gc-target-value');
        
        if (slider && valueDisplay) {
            valueDisplay.textContent = `${slider.value}%`;
        }
    }

    updateSequenceLength() {
        const sequenceInput = document.getElementById('optimizer-input');
        const lengthDisplay = document.getElementById('optimizer-sequence-length');
        
        if (sequenceInput && lengthDisplay) {
            const sequence = sequenceInput.value.replace(/\s/g, '');
            const cleanSequence = sequence.replace(/[^ATGC]/gi, '');
            lengthDisplay.textContent = `Length: ${cleanSequence.length} nucleotides`;
        }
    }

    updateGCContent() {
        const sequenceInput = document.getElementById('optimizer-input');
        const gcDisplay = document.getElementById('current-gc-content');
        
        if (sequenceInput && gcDisplay) {
            const sequence = sequenceInput.value.replace(/\s/g, '').toUpperCase();
            const cleanSequence = sequence.replace(/[^ATGC]/g, '');
            
            if (cleanSequence.length > 0) {
                const gcCount = (cleanSequence.match(/[GC]/g) || []).length;
                const gcPercent = ((gcCount / cleanSequence.length) * 100).toFixed(1);
                gcDisplay.textContent = `${gcPercent}%`;
            } else {
                gcDisplay.textContent = '0%';
            }
        }
    }

    clearInput() {
        const sequenceInput = document.getElementById('optimizer-input');
        if (sequenceInput) {
            sequenceInput.value = '';
            this.updateSequenceLength();
            this.updateGCContent();
        }
        this.hideAllOutputs();
    }

    async performOptimization() {
        try {
            // Get input values
            const sequenceInput = document.getElementById('optimizer-input');
            const sequence = sequenceInput?.value.trim();

            if (!sequence) {
                this.showNotification('Please enter a DNA sequence', 'error');
                return;
            }

            // Get optimization options
            const gcOptimizer = document.getElementById('gc-optimizer')?.value;
            const codonOptimizer = document.getElementById('codon-optimizer')?.value;
            const avoidRareCodons = document.getElementById('avoid-rare-codons')?.checked;
            const removeRestrictionSites = document.getElementById('remove-restriction-sites')?.checked;
            const optimizeFolding = document.getElementById('optimize-folding')?.checked;

            // Get custom GC target if selected
            let gcTarget = gcOptimizer;
            let customGC = null;
            if (gcOptimizer === 'custom') {
                const slider = document.getElementById('gc-target-slider');
                customGC = slider ? parseFloat(slider.value) : 50;
            }

            // Show loading state
            this.setLoadingState(true);

            // Call Python script
            const result = await this.callPythonScript(
                sequence, 
                gcTarget, 
                customGC,
                codonOptimizer, 
                avoidRareCodons, 
                removeRestrictionSites, 
                optimizeFolding
            );

            if (result.success) {
                this.displayResults(result);
                this.showNotification('Sequence optimization completed successfully!', 'success');
            } else {
                throw new Error(result.error || 'Optimization failed');
            }

        } catch (error) {
            console.error('Optimization error:', error);
            this.showNotification(`Optimization failed: ${error.message}`, 'error');
        } finally {
            this.setLoadingState(false);
        }
    }

    async callPythonScript(sequence, gcTarget, customGC, organism, avoidRareCodons, removeRestrictionSites, optimizeFolding) {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Python script
            const scriptPath = path.join(__dirname, 'assets', 'codon-optimizer.py');
            const args = [
                scriptPath,
                '--sequence', sequence
            ];

            // Add GC target
            if (gcTarget && gcTarget !== '') {
                if (gcTarget === 'custom' && customGC) {
                    args.push('--gc-target', 'custom');
                    args.push('--custom-gc', customGC.toString());
                } else {
                    args.push('--gc-target', gcTarget);
                }
            }

            // Add organism
            if (organism && organism !== '') {
                args.push('--organism', organism);
            }

            // Add advanced options
            if (avoidRareCodons) {
                args.push('--avoid-rare-codons');
            }
            if (removeRestrictionSites) {
                args.push('--remove-restriction-sites');
            }
            if (optimizeFolding) {
                args.push('--optimize-folding');
            }

            console.log('Calling Codon Optimizer:', scriptPath);
            console.log('Full command:', 'python', args.join(' '));

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
                console.log(`Python process exited with code: ${code}`);
                
                if (code === 0) {
                    try {
                        const result = JSON.parse(stdout);
                        console.log('Codon Optimizer result:', result);
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

    displayResults(result) {
        // Hide empty state
        const emptyState = document.getElementById('optimizer-empty-state');
        if (emptyState) {
            emptyState.style.display = 'none';
        }

        // Show optimized sequence
        const optimizedOutput = document.getElementById('optimized-output');
        const optimizedResult = document.getElementById('optimized-result');
        const outputGCContent = document.getElementById('output-gc-content');

        if (optimizedOutput && optimizedResult) {
            optimizedOutput.classList.remove('hidden');
            optimizedOutput.style.display = 'block';
            optimizedResult.textContent = result.optimized_sequence;
            
            if (outputGCContent) {
                outputGCContent.textContent = `${result.statistics.optimized_gc_content}%`;
            }
        }

        // Show statistics
        const statsOutput = document.getElementById('optimization-stats');
        if (statsOutput) {
            statsOutput.classList.remove('hidden');
            statsOutput.style.display = 'block';

            // Update individual statistics
            const gcContentEl = document.getElementById('gc-content');
            const caiScoreEl = document.getElementById('cai-score');
            const changesMadeEl = document.getElementById('changes-made');
            const efficiencyScoreEl = document.getElementById('efficiency-score');

            if (gcContentEl) gcContentEl.textContent = `${result.statistics.optimized_gc_content}%`;
            if (caiScoreEl) caiScoreEl.textContent = result.statistics.cai_score.toFixed(3);
            if (changesMadeEl) changesMadeEl.textContent = result.statistics.changes_made;
            if (efficiencyScoreEl) efficiencyScoreEl.textContent = `${result.statistics.efficiency_score}%`;
        }

        // Log result info
        console.log('Optimization Info:', result.optimization_info);
        console.log('Statistics:', result.statistics);
    }

    hideAllOutputs() {
        const outputElements = ['optimized-output', 'optimization-stats'];
        outputElements.forEach(id => {
            const element = document.getElementById(id);
            if (element) {
                element.classList.add('hidden');
                element.style.display = 'none';
            }
        });

        // Show empty state
        const emptyState = document.getElementById('optimizer-empty-state');
        if (emptyState) {
            emptyState.style.display = 'block';
        }
    }

    setLoadingState(loading) {
        const optimizeBtn = document.getElementById('optimize-btn');
        if (optimizeBtn) {
            if (loading) {
                optimizeBtn.disabled = true;
                optimizeBtn.innerHTML = '⚡ Optimizing...';
                optimizeBtn.classList.add('opacity-50');
            } else {
                optimizeBtn.disabled = false;
                optimizeBtn.innerHTML = '⚡ Optimize Sequence';
                optimizeBtn.classList.remove('opacity-50');
            }
        }
    }

    async copyToClipboard(elementId) {
        try {
            const element = document.getElementById(elementId);
            if (element) {
                await navigator.clipboard.writeText(element.textContent);
                this.showNotification('Copied to clipboard!', 'success');
            }
        } catch (error) {
            console.error('Copy failed:', error);
            this.showNotification('Failed to copy to clipboard', 'error');
        }
    }

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
                    <p class="text-gray-400 text-xs mt-2">This confirmation ensures you understand the risks</p>
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
            
            // Handle escape key
            const handleEscape = (e) => {
                if (e.key === 'Escape') {
                    document.body.removeChild(overlay);
                    document.removeEventListener('keydown', handleEscape);
                    resolve(false);
                }
            };
            document.addEventListener('keydown', handleEscape);
            
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
    if (document.getElementById('codon-optimizer-page')) {
        window.codonOptimizerPython = new CodonOptimizerPython();
    }
});