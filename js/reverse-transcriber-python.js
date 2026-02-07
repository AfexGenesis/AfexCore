// ===== REVERSE TRANSCRIBER - PURE PYTHON EXECUTION =====
// Direct Python script execution from Electron renderer

class PythonReverseTranscriber {
    constructor() {
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUpload();
        console.log('Python Reverse Transcriber initialized');
    }

    bindEvents() {
        // Reverse transcribe button
        const transcribeBtn = document.querySelector('#reverse-transcriber-page .bg-gradient-to-r');
        if (transcribeBtn) {
            transcribeBtn.addEventListener('click', () => this.performReverseTranscription());
        }

        // Clear button
        const clearBtn = document.getElementById('reverse-clear-btn');
        if (clearBtn) {
            clearBtn.addEventListener('click', () => this.clearInput());
        }

        // Input textarea
        const sequenceInput = document.getElementById('sequence-input');
        if (sequenceInput) {
            sequenceInput.addEventListener('input', () => this.updateSequenceLength());
        }

        // Input type radio buttons
        const inputTypeRadios = document.querySelectorAll('input[name="input-type"]');
        inputTypeRadios.forEach(radio => {
            radio.addEventListener('change', () => this.updateInputTypeUI());
        });

        // Copy buttons
        this.bindCopyButtons();
        
        // Initialize UI
        this.updateInputTypeUI();
        this.updateSequenceLength();
    }

    bindCopyButtons() {
        const copyButtons = [
            { id: 'copy-dna', resultId: 'dna-result' }
        ];

        copyButtons.forEach(({ id, resultId }) => {
            const button = document.getElementById(id);
            if (button) {
                button.addEventListener('click', () => this.copyToClipboard(resultId));
            }
        });
    }

    setupFileUpload() {
        const fileInput = document.getElementById('reverse-file-upload');
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
            dropZone.classList.add('border-purple-500/70');
        });

        dropZone.addEventListener('dragleave', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-purple-500/70');
        });

        dropZone.addEventListener('drop', (e) => {
            e.preventDefault();
            dropZone.classList.remove('border-purple-500/70');
            
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
                const proceed = await this.showLargeFileWarning(file, 'RNA sequence', maxRecommendedSize);
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
            }

            // Use Reverse Transcriber to parse the file and extract sequence
            const result = await this.callReverseTranscriberScript('dummy', 'mrna', ['dna'], fileContent, fileFormat);
            
            if (result.success) {
                // Extract the clean sequence and put it in the textarea
                const sequenceInput = document.getElementById('sequence-input');
                if (sequenceInput && result.original_sequence) {
                    sequenceInput.value = result.original_sequence;
                    this.updateSequenceLength();
                }
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

    updateInputTypeUI() {
        const selectedType = document.querySelector('input[name="input-type"]:checked')?.value;
        const sequenceInput = document.getElementById('sequence-input');
        
        if (sequenceInput && selectedType) {
            // Update placeholder based on input type
            const placeholders = {
                'mrna': 'Enter mRNA sequence (e.g., AUGCGAUCG...)',
                'trna': 'Enter tRNA sequence (e.g., UACGCUAGC...)',
                'protein': 'Enter protein sequence (e.g., MRS...)'
            };
            
            sequenceInput.placeholder = placeholders[selectedType] || 'Enter sequence...';
        }
    }

    updateSequenceLength() {
        const sequenceInput = document.getElementById('sequence-input');
        const lengthDisplay = document.getElementById('reverse-sequence-length');
        
        if (sequenceInput && lengthDisplay) {
            const sequence = sequenceInput.value.replace(/\s/g, '');
            const selectedType = document.querySelector('input[name="input-type"]:checked')?.value;
            
            let unit = 'characters';
            if (selectedType === 'mrna' || selectedType === 'trna') {
                unit = 'nucleotides';
            } else if (selectedType === 'protein') {
                unit = 'amino acids';
            }
            
            lengthDisplay.textContent = `Length: ${sequence.length} ${unit}`;
        }
    }

    clearInput() {
        const sequenceInput = document.getElementById('sequence-input');
        if (sequenceInput) {
            sequenceInput.value = '';
            this.updateSequenceLength();
            this.hideAllOutputs();
        }
    }

    async performReverseTranscription() {
        try {
            // Get input values
            const sequenceInput = document.getElementById('sequence-input');
            const sequence = sequenceInput?.value.trim();

            if (!sequence) {
                this.showNotification('Please enter a sequence', 'error');
                return;
            }

            // Get selected input type
            const inputType = document.querySelector('input[name="input-type"]:checked')?.value;
            if (!inputType) {
                this.showNotification('Please select an input type', 'error');
                return;
            }

            // Show loading state
            this.setLoadingState(true);

            // Call Reverse Transcriber Python script - always output DNA
            const result = await this.callReverseTranscriberScript(sequence, inputType, ['dna']);

            if (result.success) {
                // Display only DNA result
                this.displayResults({ dna: result.primary_dna }, result.sequence_info);
                this.showNotification('Reverse transcription completed successfully!', 'success');
            } else {
                throw new Error(result.error || 'Reverse transcription failed');
            }

        } catch (error) {
            console.error('Reverse transcription error:', error);
            this.showNotification(`Reverse transcription failed: ${error.message}`, 'error');
        } finally {
            this.setLoadingState(false);
        }
    }

    async callReverseTranscriberScript(sequence, inputType, outputTypes, fileContent = null, fileFormat = null) {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Reverse Transcriber Python script
            const scriptPath = path.join(__dirname, 'assets', 'reverse-transcriber.py');
            const args = [
                scriptPath,
                '--sequence', sequence,
                '--input-type', inputType,
                '--output-types', outputTypes.join(','),
                '--codon-usage', 'optimal'
            ];

            // Add file content if provided
            if (fileContent && fileFormat) {
                args.push('--file-content', fileContent);
                args.push('--file-format', fileFormat);
            }

            console.log('Calling Reverse Transcriber:', scriptPath);
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
                        console.log('Reverse Transcriber result:', result);
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

    displayResults(results, sequenceInfo) {
        // Hide empty state
        const emptyOutput = document.getElementById('reverse-empty-output');
        if (emptyOutput) {
            emptyOutput.style.display = 'none';
        }

        // Display each result type
        Object.keys(results).forEach(type => {
            const outputDiv = document.getElementById(`${type}-output`);
            const resultElement = document.getElementById(`${type}-result`);

            if (outputDiv && resultElement) {
                outputDiv.classList.remove('hidden');
                outputDiv.style.display = 'block';
                resultElement.textContent = results[type];
            }
        });

        // Log sequence info
        console.log('Sequence Info:', sequenceInfo);
    }

    hideAllOutputs() {
        const outputTypes = ['dna'];
        outputTypes.forEach(type => {
            const outputDiv = document.getElementById(`${type}-output`);
            if (outputDiv) {
                outputDiv.classList.add('hidden');
                outputDiv.style.display = 'none';
            }
        });

        // Show empty state
        const emptyOutput = document.getElementById('reverse-empty-output');
        if (emptyOutput) {
            emptyOutput.style.display = 'block';
        }
    }

    setLoadingState(isLoading) {
        const transcribeBtn = document.querySelector('#reverse-transcriber-page .bg-gradient-to-r');
        if (transcribeBtn) {
            if (isLoading) {
                transcribeBtn.disabled = true;
                transcribeBtn.innerHTML = '🔄 Processing...';
                transcribeBtn.classList.add('opacity-50');
            } else {
                transcribeBtn.disabled = false;
                transcribeBtn.innerHTML = '🔄 Reverse Transcribe';
                transcribeBtn.classList.remove('opacity-50');
            }
        }
    }

    async copyToClipboard(resultId) {
        const resultElement = document.getElementById(resultId);
        if (!resultElement) return;

        try {
            await navigator.clipboard.writeText(resultElement.textContent);
            this.showNotification('Copied to clipboard!', 'success');
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
                notification.classList.add('bg-purple-600');
                break;
            case 'error':
                notification.classList.add('bg-red-600');
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

        // Remove after 3 seconds
        setTimeout(() => {
            notification.classList.add('translate-x-full');
            setTimeout(() => {
                document.body.removeChild(notification);
            }, 300);
        }, 3000);
    }
}

// Initialize Python Reverse Transcriber when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    window.pythonReverseTranscriber = new PythonReverseTranscriber();
});