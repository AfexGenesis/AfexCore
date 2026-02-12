class BaseCounterPython {
    constructor() {
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUpload();
        this.updateSequenceLength();
        this.updateRealtimeGC();
    }

    bindEvents() {
        // Analyze button
        const analyzeBtn = document.getElementById('analyze-btn');
        if (analyzeBtn) {
            analyzeBtn.addEventListener('click', () => this.performAnalysis());
        }

        // Clear button
        const clearBtn = document.getElementById('counter-clear-btn');
        if (clearBtn) {
            clearBtn.addEventListener('click', () => this.clearInput());
        }

        // Input textarea - real-time updates
        const sequenceInput = document.getElementById('counter-input');
        if (sequenceInput) {
            sequenceInput.addEventListener('input', () => {
                this.updateSequenceLength();
                this.updateRealtimeGC();
            });
        }
    }

    setupFileUpload() {
        const fileInput = document.getElementById('counter-file-upload');
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
            dropZone.classList.add('border-indigo-500/50');
        });

        dropZone.addEventListener('dragleave', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-indigo-500/50');
        });

        dropZone.addEventListener('drop', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-indigo-500/50');
            
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
            const sequenceInput = document.getElementById('counter-input');
            if (sequenceInput) {
                sequenceInput.value = sequence;
                this.updateSequenceLength();
                this.updateRealtimeGC();
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
        const sequenceInput = document.getElementById('counter-input');
        const lengthDisplay = document.getElementById('counter-sequence-length');
        
        if (sequenceInput && lengthDisplay) {
            const sequence = sequenceInput.value.replace(/\s/g, '');
            const cleanSequence = sequence.replace(/[^ATGC]/gi, '');
            lengthDisplay.textContent = `Length: ${cleanSequence.length} bases`;
        }
    }

    updateRealtimeGC() {
        const sequenceInput = document.getElementById('counter-input');
        const gcDisplay = document.getElementById('realtime-gc-content');
        
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
        const sequenceInput = document.getElementById('counter-input');
        if (sequenceInput) {
            sequenceInput.value = '';
            this.updateSequenceLength();
            this.updateRealtimeGC();
        }
        this.hideAllResults();
    }

    async performAnalysis() {
        try {
            // Get input sequence
            const sequenceInput = document.getElementById('counter-input');
            const sequence = sequenceInput?.value.trim();

            if (!sequence) {
                this.showNotification('Please enter a DNA sequence', 'error');
                return;
            }

            // Show loading state
            this.setLoadingState(true);

            // Call Python script
            const result = await this.callPythonScript(sequence);

            if (result.success) {
                this.displayResults(result);
                this.showNotification('Base analysis completed successfully!', 'success');
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

    async callPythonScript(sequence) {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Python script
            const scriptPath = path.join(__dirname, '..', 'assets', 'base-counter.py');
            const args = [scriptPath];

            console.log('Calling Base Counter:', scriptPath);

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

            // Send input data via stdin as JSON
            const inputData = {
                sequence: sequence
            };

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
                        console.log('Base Counter result:', result);
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
        // Hide empty states
        const emptyState = document.getElementById('counter-empty-state');
        const emptySummary = document.getElementById('empty-summary');
        
        if (emptyState) emptyState.style.display = 'none';
        if (emptySummary) emptySummary.style.display = 'none';

        // Show base results
        const baseResults = document.getElementById('base-results');
        if (baseResults) {
            baseResults.classList.remove('hidden');
            baseResults.style.display = 'block';

            // Update individual base counts
            this.updateBaseCount('adenine', 'A', result.base_counts.A, result.base_percentages.A);
            this.updateBaseCount('thymine', 'T', result.base_counts.T, result.base_percentages.T);
            this.updateBaseCount('guanine', 'G', result.base_counts.G, result.base_percentages.G);
            this.updateBaseCount('cytosine', 'C', result.base_counts.C, result.base_percentages.C);

            // Update GC content highlight
            const finalGCContent = document.getElementById('final-gc-content');
            const gcClassification = document.getElementById('gc-classification');
            
            if (finalGCContent) finalGCContent.textContent = `${result.gc_content}%`;
            if (gcClassification) gcClassification.textContent = result.gc_classification;
        }

        // Show sequence summary
        const sequenceSummary = document.getElementById('sequence-summary');
        if (sequenceSummary) {
            sequenceSummary.classList.remove('hidden');
            sequenceSummary.style.display = 'block';

            // Update summary statistics
            const totalBasesCount = document.getElementById('total-bases-count');
            const summaryGCContent = document.getElementById('summary-gc-content');
            const atContent = document.getElementById('at-content');
            const meltingTemp = document.getElementById('melting-temp');

            if (totalBasesCount) totalBasesCount.textContent = result.total_bases;
            if (summaryGCContent) summaryGCContent.textContent = `${result.gc_content}%`;
            if (atContent) atContent.textContent = `${result.at_content}%`;
            if (meltingTemp) meltingTemp.textContent = `${result.melting_temperature}°C`;
        }

        // Log result info
        console.log('Base Analysis:', result);
    }

    updateBaseCount(baseName, baseLetter, count, percentage) {
        const countElement = document.getElementById(`${baseName}-count`);
        const percentElement = document.getElementById(`${baseName}-percent`);
        
        if (countElement) countElement.textContent = count;
        if (percentElement) percentElement.textContent = `${percentage}%`;
    }

    hideAllResults() {
        const resultElements = ['base-results', 'sequence-summary'];
        resultElements.forEach(id => {
            const element = document.getElementById(id);
            if (element) {
                element.classList.add('hidden');
                element.style.display = 'none';
            }
        });

        // Show empty states
        const emptyState = document.getElementById('counter-empty-state');
        const emptySummary = document.getElementById('empty-summary');
        
        if (emptyState) emptyState.style.display = 'block';
        if (emptySummary) emptySummary.style.display = 'block';
    }

    setLoadingState(loading) {
        const analyzeBtn = document.getElementById('analyze-btn');
        if (analyzeBtn) {
            if (loading) {
                analyzeBtn.disabled = true;
                analyzeBtn.innerHTML = '📊 Analyzing...';
                analyzeBtn.classList.add('opacity-50');
            } else {
                analyzeBtn.disabled = false;
                analyzeBtn.innerHTML = '📊 Analyze Sequence';
                analyzeBtn.classList.remove('opacity-50');
            }
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
    if (document.getElementById('base-counter-page')) {
        window.baseCounterPython = new BaseCounterPython();
    }
});