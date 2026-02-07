class ChromosomeSplitterPython {
    constructor() {
        this.fastaFile = null;
        this.chromosomeData = null;
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUploads();
    }

    bindEvents() {
        // Split button
        const splitBtn = document.getElementById('split-chromosomes-btn');
        if (splitBtn) {
            splitBtn.addEventListener('click', () => this.performSplit());
        }

        // Download all button
        const downloadAllBtn = document.getElementById('download-all-btn');
        if (downloadAllBtn) {
            downloadAllBtn.addEventListener('click', () => this.downloadAllChromosomes());
        }

        // Remove file button
        const removeFastaBtn = document.getElementById('remove-fasta');
        
        if (removeFastaBtn) {
            removeFastaBtn.addEventListener('click', () => this.removeFastaFile());
        }
    }

    setupFileUploads() {
        // FASTA file upload
        const fastaInput = document.getElementById('fasta-file-upload');
        const fastaDropZone = fastaInput?.parentElement;

        if (fastaInput && fastaDropZone) {
            fastaInput.addEventListener('change', (e) => {
                if (e.target.files.length > 0) {
                    this.handleFastaUpload(e.target.files[0]);
                }
            });

            // Drag and drop for FASTA
            fastaDropZone.addEventListener('dragover', (e) => {
                e.preventDefault();
                fastaDropZone.classList.add('border-rose-500/50');
            });

            fastaDropZone.addEventListener('dragleave', (e) => {
                e.preventDefault();
                fastaDropZone.classList.remove('border-rose-500/50');
            });

            fastaDropZone.addEventListener('drop', (e) => {
                e.preventDefault();
                fastaDropZone.classList.remove('border-rose-500/50');
                
                if (e.dataTransfer.files.length > 0) {
                    this.handleFastaUpload(e.dataTransfer.files[0]);
                }
            });
        }
    }

    async handleFastaUpload(file) {
        try {
            console.log('📁 Starting FASTA upload:', {
                name: file.name,
                size: file.size,
                type: file.type
            });

            // Validate file type
            const fileName = file.name.toLowerCase();
            if (!fileName.endsWith('.fasta') && !fileName.endsWith('.fa') && !fileName.endsWith('.fas')) {
                this.showNotification('Please upload a valid FASTA file (.fasta, .fa, .fas)', 'error');
                return;
            }

            // Show loading for large files
            if (file.size > 100 * 1024 * 1024) { // 100MB
                this.showLoadingOverlay('Reading FASTA file...');
            }

            console.log('Reading file...');
            const fileContent = await this.readFile(file);
            
            // Check if we got a temp file path instead of content (for very large files)
            if (typeof fileContent === 'string' && fileContent.startsWith('TEMP_FILE:')) {
                const tempFilePath = fileContent.substring('TEMP_FILE:'.length);
                console.log('Large file streamed to temp location:', tempFilePath);
                
                this.fastaFile = {
                    name: file.name,
                    size: file.size,
                    content: null, // No content in memory
                    tempFilePath: tempFilePath, // Path to temp file
                    isLargeFile: true
                };
            } else {
                console.log('File read successfully, content length:', fileContent.length);
                
                this.fastaFile = {
                    name: file.name,
                    size: file.size,
                    content: fileContent,
                    tempFilePath: null,
                    isLargeFile: false
                };
            }

            // Update UI
            this.updateFastaFileInfo();
            this.updateSplitButton();
            
            this.hideLoadingOverlay();
            this.showNotification(`FASTA file "${file.name}" loaded successfully`, 'success');

        } catch (error) {
            console.error('FASTA upload error:', error);
            this.hideLoadingOverlay();
            this.showNotification(`Failed to load FASTA file: ${error.message}`, 'error');
        }
    }

    readFile(file) {
        return new Promise((resolve, reject) => {
            if (!file) {
                reject(new Error('No file provided'));
                return;
            }
            
            if (file.size === 0) {
                reject(new Error('File is empty'));
                return;
            }
            
            // For very large files (>10GB), show error
            if (file.size > 10 * 1024 * 1024 * 1024) { // > 10GB
                reject(new Error(`File "${file.name}" is too large (${this.formatFileSize(file.size)}). Maximum supported size is 10GB.`));
                return;
            }

            // For large files (>1GB), stream to temp file instead of loading into memory
            if (file.size > 1 * 1024 * 1024 * 1024) { // > 1GB
                this.streamLargeFileToTemp(file, resolve, reject);
            } else {
                // For smaller files, read normally
                const reader = new FileReader();
                reader.onload = (e) => resolve(e.target.result);
                reader.onerror = (e) => {
                    console.error('FileReader error:', e);
                    reject(new Error(`Failed to read file: ${file.name}. File may be corrupted or inaccessible.`));
                };
                reader.readAsText(file);
            }
        });
    }

    async streamLargeFileToTemp(file, resolve, reject) {
        try {
            const fileSizeGB = (file.size / (1024 * 1024 * 1024)).toFixed(2);
            const fileSizeMB = (file.size / (1024 * 1024)).toFixed(1);
            
            console.log(`LARGE FASTA FILE DETECTED`);
            console.log(`File: ${file.name}`);
            console.log(`Size: ${fileSizeGB} GB (${fileSizeMB} MB)`);
            console.log(`Creating temporary copy for processing...`);
            
            // Show user notification about large file processing
            this.showNotification(
                `Processing large file (${fileSizeGB}GB). Creating temporary copy for reliable processing...`, 
                'info', 
                5000
            );
            
            // Create temporary file using Node.js
            const fs = require('fs');
            const path = require('path');
            const os = require('os');
            
            // Create temporary file
            const tempDir = os.tmpdir();
            // Truncate filename if too long to avoid ENAMETOOLONG
            const safeName = file.name.length > 50 ? file.name.substring(0, 50) + '...' + path.extname(file.name) : file.name;
            const tempFileName = `chromosome_splitter_${Date.now()}_${safeName}`;
            const tempFilePath = path.join(tempDir, tempFileName);
            
            console.log(`Creating temporary copy: ${tempFilePath}`);
            
            // Stream file to disk in chunks
            await this.streamFileToDisk(file, tempFilePath);
            
            console.log(`Reading temp file back in chunks to avoid 2GB limit...`);
            
            // For files >2GB, we can't load into JavaScript string - use streaming approach
            if (file.size > 2 * 1024 * 1024 * 1024) { // > 2GB
                console.log(`File is >2GB - using streaming approach without loading into memory`);
                
                // Return a special marker indicating this is a temp file
                // The processing functions will handle streaming directly from disk
                resolve(`TEMP_FILE:${tempFilePath}`);
                
                // Don't clean up temp file yet - it will be cleaned up after processing
                return;
            }
            
            // For files 1-2GB, try to read into memory
            this.readLargeFileInChunks(tempFilePath)
                .then(tempContent => {
                    // Clean up temp file
                    fs.unlink(tempFilePath, (unlinkErr) => {
                        if (unlinkErr) {
                            console.warn(`Could not delete temp file: ${tempFilePath}`);
                        } else {
                            console.log(`Cleaned up temp file: ${tempFilePath}`);
                        }
                    });
                    
                    console.log(`Large file read successfully in chunks, content length: ${tempContent.length}`);
                    resolve(tempContent);
                })
                .catch(err => {
                    // Clean up temp file on error
                    fs.unlink(tempFilePath, () => {});
                    reject(new Error(`Failed to read large temp file: ${err.message}`));
                });
            
        } catch (error) {
            console.error('Failed to stream large file:', error);
            reject(new Error(`Failed to process large file: ${error.message}`));
        }
    }

    streamFileToDisk(file, tempFilePath) {
        return new Promise((resolve, reject) => {
            const fs = require('fs');
            const chunkSize = 64 * 1024 * 1024; // 64MB chunks
            let offset = 0;
            
            // Create write stream
            const writeStream = fs.createWriteStream(tempFilePath);
            
            const writeNextChunk = () => {
                if (offset >= file.size) {
                    writeStream.end();
                    resolve();
                    return;
                }

                const chunk = file.slice(offset, offset + chunkSize);
                const reader = new FileReader();
                
                reader.onload = (e) => {
                    writeStream.write(e.target.result);
                    offset += chunkSize;
                    
                    // Show progress
                    const progress = Math.min(100, Math.round((offset / file.size) * 100));
                    console.log(`Creating temp copy: ${progress}%`);
                    
                    // Continue writing
                    setTimeout(writeNextChunk, 10);
                };
                
                reader.onerror = (e) => {
                    writeStream.destroy();
                    reject(new Error(`Failed to read chunk at position ${offset}`));
                };
                
                reader.readAsText(chunk);
            };

            writeStream.on('error', (error) => {
                reject(new Error(`Failed to write to temporary file: ${error.message}`));
            });

            writeNextChunk();
        });
    }

    readLargeFileInChunks(filePath) {
        return new Promise((resolve, reject) => {
            const fs = require('fs');
            const chunkSize = 64 * 1024 * 1024; // 64MB chunks
            let contentChunks = []; // Store chunks in array instead of concatenating
            let position = 0;
            
            // Get file size first
            fs.stat(filePath, (err, stats) => {
                if (err) {
                    reject(new Error(`Failed to get file stats: ${err.message}`));
                    return;
                }
                
                const fileSize = stats.size;
                console.log(`Reading ${(fileSize / (1024 * 1024 * 1024)).toFixed(2)}GB file in ${Math.ceil(fileSize / chunkSize)} chunks...`);
                
                const readNextChunk = () => {
                    if (position >= fileSize) {
                        // Join all chunks at the end - this is more memory efficient
                        console.log(`Joining ${contentChunks.length} chunks into final content...`);
                        try {
                            const finalContent = contentChunks.join('');
                            contentChunks = null; // Free memory
                            console.log(`Large file read successfully, final content length: ${finalContent.length}`);
                            resolve(finalContent);
                        } catch (joinError) {
                            console.error('Failed to join chunks - file too large for JavaScript string:', joinError);
                            reject(new Error(`File is too large to process in JavaScript. The file exceeds JavaScript's string length limit. Consider using a smaller file or processing it in smaller sections.`));
                        }
                        return;
                    }
                    
                    const currentChunkSize = Math.min(chunkSize, fileSize - position);
                    const buffer = Buffer.alloc(currentChunkSize);
                    
                    fs.open(filePath, 'r', (err, fd) => {
                        if (err) {
                            reject(new Error(`Failed to open file: ${err.message}`));
                            return;
                        }
                        
                        fs.read(fd, buffer, 0, currentChunkSize, position, (err, bytesRead) => {
                            fs.close(fd, () => {}); // Close file descriptor
                            
                            if (err) {
                                reject(new Error(`Failed to read chunk: ${err.message}`));
                                return;
                            }
                            
                            // Store chunk in array instead of concatenating immediately
                            contentChunks.push(buffer.toString('utf8', 0, bytesRead));
                            position += bytesRead;
                            
                            // Show progress
                            const progress = Math.min(100, Math.round((position / fileSize) * 100));
                            console.log(`Reading chunks: ${progress}%`);
                            
                            // Continue reading next chunk
                            setTimeout(readNextChunk, 10);
                        });
                    });
                };
                
                readNextChunk();
            });
        });
    }

    updateFastaFileInfo() {
        const fileInfo = document.getElementById('fasta-file-info');
        const fileName = document.getElementById('fasta-file-name');
        const fileSize = document.getElementById('fasta-file-size');

        if (fileInfo && fileName && fileSize && this.fastaFile) {
            fileInfo.classList.remove('hidden');
            fileName.textContent = `FASTA: ${this.fastaFile.name}`;
            fileSize.textContent = this.formatFileSize(this.fastaFile.size);
        }
    }

    updateSplitButton() {
        const splitBtn = document.getElementById('split-chromosomes-btn');
        if (splitBtn) {
            if (this.fastaFile) {
                splitBtn.disabled = false;
                splitBtn.classList.remove('opacity-50');
            } else {
                splitBtn.disabled = true;
                splitBtn.classList.add('opacity-50');
            }
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
        notification.className = `fixed top-4 right-4 p-4 rounded-lg shadow-lg z-50 max-w-md ${
            type === 'error' ? 'bg-red-600 text-white' :
            type === 'success' ? 'bg-green-600 text-white' :
            'bg-blue-600 text-white'
        }`;
        notification.textContent = message;

        document.body.appendChild(notification);

        // Auto remove after 5 seconds
        setTimeout(() => {
            if (notification.parentNode) {
                notification.parentNode.removeChild(notification);
            }
        }, 5000);
    }

    showLoadingOverlay(message) {
        let overlay = document.getElementById('loading-overlay');
        if (!overlay) {
            overlay = document.createElement('div');
            overlay.id = 'loading-overlay';
            overlay.className = 'fixed inset-0 bg-black/50 flex items-center justify-center z-50';
            overlay.innerHTML = `
                <div class="bg-slate-800 border border-blue-500/30 rounded-xl p-6 text-center">
                    <div class="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-500 mx-auto mb-4"></div>
                    <p class="text-white text-lg">${message}</p>
                </div>
            `;
            document.body.appendChild(overlay);
        }
    }

    hideLoadingOverlay() {
        const overlay = document.getElementById('loading-overlay');
        if (overlay) {
            overlay.remove();
        }
    }

    removeFastaFile() {
        // Clean up temp file if it exists
        if (this.fastaFile && this.fastaFile.isLargeFile && this.fastaFile.tempFilePath) {
            const fs = require('fs');
            fs.unlink(this.fastaFile.tempFilePath, (err) => {
                if (err) {
                    console.warn(`Could not delete temp file: ${this.fastaFile.tempFilePath}`);
                } else {
                    console.log(`Cleaned up temp file: ${this.fastaFile.tempFilePath}`);
                }
            });
        }
        
        // Clean up chromosome temp directory if it exists
        if (this.chromosomeData && this.chromosomeData.tempDir) {
            const fs = require('fs');
            fs.rmSync(this.chromosomeData.tempDir, { recursive: true, force: true });
            console.log(`Cleaned up chromosome temp directory: ${this.chromosomeData.tempDir}`);
        }
        
        this.fastaFile = null;
        this.chromosomeData = null;
        
        // Hide file info
        const fileInfo = document.getElementById('fasta-file-info');
        if (fileInfo) {
            fileInfo.classList.add('hidden');
        }
        
        // Clear file input
        const fastaInput = document.getElementById('fasta-file-upload');
        if (fastaInput) {
            fastaInput.value = '';
        }
        
        // Update buttons
        this.updateSplitButton();
        
        // Clear results
        const resultsContainer = document.getElementById('chromosome-results');
        if (resultsContainer) {
            resultsContainer.innerHTML = '';
            resultsContainer.classList.add('hidden');
        }
        
        this.showNotification('FASTA file removed', 'info');
    }

    cleanupTempFiles() {
        try {
            const fs = require('fs');
            
            // Clean up FASTA temp file
            if (this.fastaFile && this.fastaFile.isLargeFile && this.fastaFile.tempFilePath) {
                fs.unlink(this.fastaFile.tempFilePath, (err) => {
                    if (err) {
                        console.warn(`Could not delete temp file: ${this.fastaFile.tempFilePath}`);
                    } else {
                        console.log(`Cleaned up temp file: ${this.fastaFile.tempFilePath}`);
                    }
                });
            }
            
            // Clean up chromosome temp directory
            if (this.chromosomeData && this.chromosomeData.tempDir) {
                fs.rmSync(this.chromosomeData.tempDir, { recursive: true, force: true });
                console.log(`Cleaned up chromosome temp directory: ${this.chromosomeData.tempDir}`);
            }
        } catch (error) {
            console.warn('Error during temp file cleanup, if u see this pls fix this code ts gave me a headache its still not working:', error);
        }
    }

    async performSplit() {
        if (!this.fastaFile) {
            this.showNotification('Please upload a FASTA file first', 'error');
            return;
        }

        try {
            this.showLoadingOverlay('Splitting chromosomes...');
            
            let chromosomes;
            let tempDir = null;
            
            // Handle large files with streaming parser
            if (this.fastaFile.isLargeFile && this.fastaFile.tempFilePath) {
                console.log('Using streaming parser for large file...');
                const result = await this.parseFastaStreaming(this.fastaFile.tempFilePath);
                chromosomes = result.chromosomes;
                tempDir = result.tempDir;
            } else {
                // Use regular parser for smaller files
                chromosomes = this.parseFasta(this.fastaFile.content);
            }
            
            if (chromosomes.length === 0) {
                throw new Error('No valid chromosomes found in FASTA file');
            }

            this.chromosomeData = chromosomes;
            this.chromosomeData.tempDir = tempDir; // Store temp directory for cleanup
            
            console.log(`About to display ${chromosomes.length} chromosomes`);
            console.log('First few chromosomes:', chromosomes.slice(0, 3));
            
            this.displayResults(chromosomes);
            
            // Clean up temp files if they exist
            if (this.fastaFile.isLargeFile && this.fastaFile.tempFilePath) {
                const fs = require('fs');
                fs.unlink(this.fastaFile.tempFilePath, (err) => {
                    if (err) {
                        console.warn(`Could not delete temp file: ${this.fastaFile.tempFilePath}`);
                    } else {
                        console.log(`Cleaned up temp file: ${this.fastaFile.tempFilePath}`);
                    }
                });
                this.fastaFile.tempFilePath = null;
            }
            
            // Store temp directory for later cleanup (don't auto-delete)
            if (tempDir) {
                this.chromosomeData.tempDir = tempDir;
                console.log(`Temp directory stored for manual cleanup: ${tempDir}`);
                
                // Add cleanup on page unload
                window.addEventListener('beforeunload', () => {
                    this.cleanupTempFiles();
                });
            }
            
            this.hideLoadingOverlay();
            this.showNotification(`Successfully split ${chromosomes.length} chromosomes`, 'success');

        } catch (error) {
            console.error('Split error:', error);
            this.hideLoadingOverlay();
            this.showNotification(`Failed to split chromosomes: ${error.message}`, 'error');
        }
    }

    parseFasta(content) {
        const chromosomes = [];
        const lines = content.split('\n');
        let currentChromosome = null;
        let currentSequence = '';

        for (let line of lines) {
            line = line.trim();
            
            if (line.startsWith('>')) {
                // Save previous chromosome
                if (currentChromosome && currentSequence) {
                    chromosomes.push({
                        name: currentChromosome,
                        sequence: currentSequence
                    });
                }
                
                // Start new chromosome
                currentChromosome = line.substring(1).trim();
                currentSequence = '';
            } else if (line && currentChromosome) {
                currentSequence += line;
            }
        }

        // Save last chromosome
        if (currentChromosome && currentSequence) {
            chromosomes.push({
                name: currentChromosome,
                sequence: currentSequence
            });
        }

        return chromosomes;
    }

    parseFastaStreaming(filePath) {
        return new Promise((resolve, reject) => {
            const fs = require('fs');
            const readline = require('readline');
            const path = require('path');
            const os = require('os');
            
            console.log('Starting streaming FASTA parser with direct-to-disk writing...');
            
            const chromosomes = [];
            let currentChromosome = null;
            let currentWriteStream = null;
            let currentSequenceLength = 0;
            let lineCount = 0;
            let chromosomeCount = 0;
            
            // Create temp directory for chromosome files
            const tempDir = path.join(os.tmpdir(), `chromosome_split_${Date.now()}`);
            fs.mkdirSync(tempDir, { recursive: true });
            console.log(`Created temp directory: ${tempDir}`);
            
            // Create read stream
            const fileStream = fs.createReadStream(filePath, { encoding: 'utf8' });
            
            // Create readline interface
            const rl = readline.createInterface({
                input: fileStream,
                crlfDelay: Infinity // Handle Windows line endings
            });
            
            const finishCurrentChromosome = () => {
                if (currentChromosome && currentWriteStream) {
                    currentWriteStream.end();
                    chromosomes.push({
                        name: currentChromosome,
                        sequence: null, // No sequence in memory
                        sequenceLength: currentSequenceLength,
                        filePath: path.join(tempDir, `${currentChromosome.replace(/[<>:"/\\|?*]/g, '_')}.fasta`),
                        isLargeFile: true
                    });
                    console.log(`Wrote chromosome: ${currentChromosome} (${currentSequenceLength.toLocaleString()} bp) to disk`);
                    chromosomeCount++;
                }
            };
            
            rl.on('line', (line) => {
                lineCount++;
                
                // Show progress every 100k lines
                if (lineCount % 100000 === 0) {
                    console.log(`Processed ${lineCount.toLocaleString()} lines, ${chromosomeCount} chromosomes...`);
                }
                
                line = line.trim();
                
                if (line.startsWith('>')) {
                    // Finish previous chromosome
                    finishCurrentChromosome();
                    
                    // Start new chromosome
                    currentChromosome = line.substring(1).trim();
                    currentSequenceLength = 0;
                    
                    // Create new file for this chromosome
                    const safeFileName = currentChromosome.replace(/[<>:"/\\|?*]/g, '_');
                    const chromosomeFilePath = path.join(tempDir, `${safeFileName}.fasta`);
                    
                    currentWriteStream = fs.createWriteStream(chromosomeFilePath);
                    currentWriteStream.write(`>${currentChromosome}\n`);
                    
                    console.log(`🧬 Found chromosome: ${currentChromosome} -> ${safeFileName}.fasta`);
                } else if (line && currentChromosome && currentWriteStream) {
                    // Write sequence line directly to file
                    currentWriteStream.write(line + '\n');
                    currentSequenceLength += line.length;
                }
            });
            
            rl.on('close', () => {
                // Finish last chromosome
                finishCurrentChromosome();
                
                console.log(`Streaming parser completed! Found ${chromosomes.length} chromosomes from ${lineCount.toLocaleString()} lines`);
                console.log(`All chromosomes written to: ${tempDir}`);
                
                // Return object with chromosomes and temp directory info
                resolve({
                    chromosomes: chromosomes,
                    tempDir: tempDir
                });
            });
            
            rl.on('error', (error) => {
                console.error('Streaming parser error:', error);
                // Clean up temp directory on error
                fs.rmSync(tempDir, { recursive: true, force: true });
                reject(new Error(`Failed to parse FASTA file: ${error.message}`));
            });
            
            fileStream.on('error', (error) => {
                console.error('File stream error wtf?:', error);
                // Clean up temp directory on error
                fs.rmSync(tempDir, { recursive: true, force: true });
                reject(new Error(`Failed to read FASTA file: ${error.message}`));
            });
        });
    }

    displayResults(chromosomes) {
        console.log(`displayResults called with ${chromosomes.length} chromosomes`);
        
        const resultsContainer = document.getElementById('chromosome-results');
        if (!resultsContainer) {
            console.error('Results container not found!');
            return;
        }

        console.log('Building HTML for chromosomes...');
        resultsContainer.innerHTML = `
            <div class="bg-slate-800 border border-green-500/30 rounded-xl p-6">
                <h3 class="text-xl font-semibold text-white mb-4">Split Results</h3>
                <div class="space-y-3">
                    ${chromosomes.map((chr, index) => `
                        <div class="bg-slate-700 border border-slate-600 rounded-lg p-4">
                            <div class="flex justify-between items-center">
                                <div>
                                    <h4 class="text-white font-medium">${chr.name}</h4>
                                    <p class="text-gray-400 text-sm">Length: ${chr.isLargeFile ? chr.sequenceLength.toLocaleString() : chr.sequence.length.toLocaleString()} bp</p>
                                </div>
                                <button onclick="chromosomeSplitter.downloadChromosome(${index})" 
                                        class="px-4 py-2 bg-blue-600/20 border border-blue-500/30 rounded text-blue-300 hover:bg-blue-600/30 transition-all">
                                    Download
                                </button>
                            </div>
                        </div>
                    `).join('')}
                </div>
                <div class="mt-6 text-center">
                    <button id="download-all-btn" 
                            class="px-6 py-3 bg-green-600/20 border border-green-500/30 rounded-lg text-green-300 hover:bg-green-600/30 transition-all">
                        Download All Chromosomes
                    </button>
                </div>
            </div>
        `;

        console.log('HTML built and inserted into DOM');
        
        // Show the results container
        resultsContainer.classList.remove('hidden');
        console.log('Results container made visible');
        
        // Re-bind download all button
        const downloadAllBtn = document.getElementById('download-all-btn');
        if (downloadAllBtn) {
            downloadAllBtn.addEventListener('click', () => this.downloadAllChromosomes());
            console.log('Download all button re-bound');
        } else {
            console.error('Download all button not found after HTML insertion');
        }
    }

    async downloadChromosome(index) {
        if (!this.chromosomeData || !this.chromosomeData[index]) {
            this.showNotification('Chromosome data not found', 'error');
            return;
        }

        const chromosome = this.chromosomeData[index];
        
        try {
            let content;
            
            if (chromosome.isLargeFile && chromosome.filePath) {
                // Read from file for large chromosomes
                const fs = require('fs');
                content = fs.readFileSync(chromosome.filePath, 'utf8');
            } else {
                // Use in-memory sequence for smaller chromosomes
                content = `>${chromosome.name}\n${chromosome.sequence}`;
            }
            
            const filename = `${chromosome.name.replace(/[^a-zA-Z0-9]/g, '_')}.fasta`;
            this.downloadFile(content, filename);
            
        } catch (error) {
            console.error('Download error:', error);
            this.showNotification(`Failed to download chromosome: ${error.message}`, 'error');
        }
    }

    async downloadAllChromosomes() {
        if (!this.chromosomeData || this.chromosomeData.length === 0) {
            this.showNotification('No chromosome data available', 'error');
            return;
        }

        try {
            this.showLoadingOverlay('Preparing downloads...');

            // Create a zip file with all chromosomes
            for (let i = 0; i < this.chromosomeData.length; i++) {
                const chromosome = this.chromosomeData[i];
                const content = `>${chromosome.name}\n${chromosome.sequence}`;
                const filename = `${chromosome.name.replace(/[^a-zA-Z0-9]/g, '_')}.fasta`;
                
                // Small delay between downloads to prevent browser blocking
                setTimeout(() => {
                    this.downloadFile(content, filename);
                }, i * 100);
            }

            this.hideLoadingOverlay();
            this.showNotification(`Started download of ${this.chromosomeData.length} chromosome files`, 'success');

        } catch (error) {
            console.error('Download error:', error);
            this.hideLoadingOverlay();
            this.showNotification(`Failed to download chromosomes: ${error.message}`, 'error');
        }
    }

    downloadFile(content, filename) {
        const blob = new Blob([content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }
}

// Initialize the chromosome splitter
const chromosomeSplitter = new ChromosomeSplitterPython();