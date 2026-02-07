// ===== GENE EXTRACTOR - PURE PYTHON EXECUTION =====
// Direct Python script execution from Electron renderer

class PythonGeneExtractor {
    constructor() {
        this.fastaFile = null;
        this.gffFile = null;
        this.extractedGenes = [];
        this.currentResults = null;
        this.currentPage = 1;
        this.genesPerPage = 500; // Match server-side pagination
        this.init();
    }

    init() {
        this.bindEvents();
        this.setupFileUploads();
        console.log('Python Gene Extractor initialized');
    }

    bindEvents() {
        // Analyze button
        const analyzeBtn = document.getElementById('gene-analyze-btn');
        if (analyzeBtn) {
            analyzeBtn.addEventListener('click', () => this.performAnalysis());
        }

        // Clear button
        const clearBtn = document.getElementById('gene-clear-btn');
        if (clearBtn) {
            clearBtn.addEventListener('click', () => this.clearResults());
        }

        // Download button
        const downloadBtn = document.getElementById('gene-download-btn');
        if (downloadBtn) {
            downloadBtn.addEventListener('click', () => this.downloadAllGenes());
        }

        // Manual sequence input
        const manualInput = document.getElementById('gene-manual-sequence');
        if (manualInput) {
            manualInput.addEventListener('input', () => this.validateInput());
        }
    }

    setupFileUploads() {
        // FASTA/GenBank file upload
        this.setupFileUpload('gene-fasta-upload-area', 'gene-fasta-file', (file) => {
            this.fastaFile = file;
            this.updateUploadStatus('fasta', file.name);
            this.validateInput();
        });

        // GFF file upload
        this.setupFileUpload('gene-gff-upload-area', 'gene-gff-file', (file) => {
            this.gffFile = file;
            this.updateUploadStatus('gff', file.name);
            this.validateInput();
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
                this.handleFileUpload(file, callback);
            }
        });

        // Drag and drop
        uploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            uploadArea.classList.add('border-blue-500/70');
        });

        uploadArea.addEventListener('dragleave', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('border-blue-500/70');
        });

        uploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('border-blue-500/70');
            
            const files = e.dataTransfer.files;
            if (files.length > 0) {
                this.handleFileUpload(files[0], callback);
            }
        });
    }

    async handleFileUpload(file, callback) {
        try {
            // Check file size limits - warning at 5GB, max at 10GB
            const maxRecommendedSize = 5 * 1024 * 1024 * 1024; // 5GB
            const maxAllowedSize = 10 * 1024 * 1024 * 1024; // 10GB

            if (file.size > maxAllowedSize) {
                this.showNotification('File size exceeds maximum limit of 10GB', 'error');
                return;
            }

            if (file.size > maxRecommendedSize) {
                const proceed = await this.showLargeFileWarning(file, 'genome file', maxRecommendedSize);
                if (!proceed) {
                    return;
                }
            }

            callback(file);
            this.showNotification(`File "${file.name}" loaded successfully`, 'success');

        } catch (error) {
            console.error('File upload error:', error);
            this.showNotification('Failed to load file. Please try again.', 'error');
        }
    }

    updateUploadStatus(type, fileName) {
        const uploadArea = document.getElementById(`gene-${type}-upload-area`);
        if (!uploadArea) return;

        // Update visual feedback
        uploadArea.classList.remove('border-slate-600');
        uploadArea.classList.add('border-green-500/50', 'bg-green-500/10');
        
        const textElement = uploadArea.querySelector('p');
        if (textElement) {
            textElement.textContent = `✓ ${fileName}`;
            textElement.classList.add('text-green-400');
        }
    }

    validateInput() {
        const analyzeBtn = document.getElementById('gene-analyze-btn');
        const manualSequence = document.getElementById('gene-manual-sequence')?.value.trim();
        
        const hasInput = this.fastaFile || manualSequence;
        
        if (analyzeBtn) {
            if (hasInput) {
                analyzeBtn.disabled = false;
                analyzeBtn.classList.remove('opacity-50', 'cursor-not-allowed');
            } else {
                analyzeBtn.disabled = true;
                analyzeBtn.classList.add('opacity-50', 'cursor-not-allowed');
            }
        }
    }

    async performAnalysis() {
        try {
            const manualSequence = document.getElementById('gene-manual-sequence')?.value.trim();
            
            // Check if we have sequence input
            if (!this.fastaFile && !manualSequence) {
                this.showNotification('Please upload a FASTA/GenBank file or enter a sequence manually', 'error');
                return;
            }

            // Check for GFF file - show warning if missing
            if (!this.gffFile) {
                const proceed = await this.showGFFWarning();
                if (!proceed) {
                    return;
                }
            }

            // Show loading state
            this.setLoadingState(true);

            // Prepare file contents
            let fastaContent = null;
            let gffContent = null;

            if (this.fastaFile) {
                fastaContent = await this.readFileContent(this.fastaFile);
            }

            if (this.gffFile) {
                gffContent = await this.readFileContent(this.gffFile);
            }

            // Call Python script - TRY LOADING ALL GENES AT ONCE!
            const result = await this.callPythonScript(
                manualSequence || fastaContent,
                gffContent,
                this.fastaFile?.name || 'manual_sequence',
                1,    // page 1 (ignored when loadAll=true)
                500,  // genes per page (ignored when loadAll=true)
                true  // loadAll = true - EXPERIMENT TIME!
            );

            if (result.success) {
                
                this.extractedGenes = result.genes; // Store ALL genes in memory
                
                // But only DISPLAY the first page (500 genes)
                const totalGenes = result.genes.length;
                const totalPages = Math.ceil(totalGenes / this.genesPerPage);
                const firstPageGenes = result.genes.slice(0, this.genesPerPage);
                
                // Create proper pagination object for first page
                const firstPagePagination = {
                    total_genes: totalGenes,
                    genes_per_page: this.genesPerPage,
                    total_pages: totalPages,
                    current_page: 1,
                    showing_genes: firstPageGenes.length,
                    start_gene: 1,
                    end_gene: Math.min(this.genesPerPage, totalGenes),
                    has_next: totalPages > 1,
                    has_prev: false
                };
                
                // Store the current results for pagination
                this.currentResults = {
                    statistics: result.statistics,
                    pagination: firstPagePagination
                };
                
                // Display only first page
                this.displayResults(firstPageGenes, result.statistics, firstPagePagination);
                
                // Show appropriate success message based on extraction method
                if (result.method === 'FAST_GFF_scan') {
                    this.showNotification(
                        `🚀 ALL ${totalGenes.toLocaleString()} GENES LOADED! Scan took ${result.statistics.scan_time}s. Instant page switching enabled! 🚀`, 
                        'success', 
                        8000
                    );
                } else if (result.method === 'GFF_annotation') {
                    if (result.pagination && result.pagination.truncated) {
                        this.showNotification(
                            `🧬 REAL GENES EXTRACTED! Found ${result.total_genes_found.toLocaleString()} actual genes from GFF annotations, showing first ${result.genes.length} to prevent memory overflow.`, 
                            'success', 
                            8000
                        );
                    } else {
                        this.showNotification(
                            `🧬 SUCCESS! Extracted ${result.genes.length.toLocaleString()} REAL genes from GFF annotations!`, 
                            'success', 
                            6000
                        );
                    }
                } else if (result.method === 'ORF_detection') {
                    this.showNotification(
                        `⚠️ Found ${result.genes.length.toLocaleString()} potential genes using ORF detection. For REAL genes, please provide a GFF file!`, 
                        'warning', 
                        8000
                    );
                } else {
                    this.showNotification(`Successfully extracted ${result.genes.length} genes!`, 'success');
                }
            } else {
                throw new Error(result.error || 'Gene extraction failed');
            }

        } catch (error) {
            console.error('Analysis error:', error);
            const errorMessage = error.message || 'Unknown error occurred';
            this.showNotification(`Analysis failed: ${errorMessage}`, 'error');
        } finally {
            this.setLoadingState(false);
        }
    }

    async callPythonScript(sequenceData, gffData, fileName, page = 1, genesPerPage = 500, loadAll = false) {
        return new Promise((resolve, reject) => {
            try {
                // Check if we're in Electron environment
                if (typeof require === 'undefined') {
                    reject(new Error('Node.js integration not available. Please run in Electron environment.'));
                    return;
                }

                const { spawn } = require('child_process');
                const path = require('path');

                // Prepare arguments for Python script
                const scriptPath = path.join(__dirname, 'assets', 'gene-extractor.py');
                const args = [scriptPath];

                // Add sequence data (handle temp files)
                if (sequenceData) {
                    if (sequenceData.startsWith('TEMP_FILE:')) {
                        const tempFilePath = sequenceData.replace('TEMP_FILE:', '');
                        args.push('--sequence-file', tempFilePath);
                    } else {
                        // For non-temp files, check if data is too long for command line
                        if (sequenceData.length > 4000) {
                            // Create a temp file for very long sequence data
                            const fs = require('fs');
                            const os = require('os');
                            const tempDir = os.tmpdir();
                            const tempSeqFile = path.join(tempDir, `seq_${Date.now()}.tmp`);
                            fs.writeFileSync(tempSeqFile, sequenceData);
                            args.push('--sequence-file', tempSeqFile);
                            console.log('Created temp file for long sequence data:', tempSeqFile);
                        } else {
                            args.push('--sequence-data', sequenceData);
                        }
                    }
                }

                // Add GFF data (handle temp files)
                if (gffData) {
                    if (gffData.startsWith('TEMP_FILE:')) {
                        const tempFilePath = gffData.replace('TEMP_FILE:', '');
                        args.push('--gff-file', tempFilePath);
                    } else {
                        // For non-temp files, check if data is too long for command line
                        if (gffData.length > 4000) {
                            // Create a temp file for very long GFF data
                            const fs = require('fs');
                            const os = require('os');
                            const tempDir = os.tmpdir();
                            const tempGffFile = path.join(tempDir, `gff_${Date.now()}.tmp`);
                            fs.writeFileSync(tempGffFile, gffData);
                            args.push('--gff-file', tempGffFile);
                            console.log('Created temp file for long GFF data:', tempGffFile);
                        } else {
                            args.push('--gff-data', gffData);
                        }
                    }
                }

                // Add file name for context (truncate if too long)
                const safeFileName = fileName.length > 100 ? fileName.substring(0, 100) + '...' : fileName;
                args.push('--file-name', safeFileName);
                
                // Add pagination parameters
                args.push('--page', page.toString());
                args.push('--genes-per-page', genesPerPage.toString());
                
                // Add load-all flag if requested
                if (loadAll) {
                    args.push('--load-all');
                }

                console.log('Calling Python Gene Extractor:', scriptPath);
                console.log('Arguments:', args);
                console.log('Command line length:', JSON.stringify(args).length);
                
                // Check for potential ENAMETOOLONG issues
                const fullCommand = `python ${args.join(' ')}`;
                console.log('Full command length:', fullCommand.length);
                if (fullCommand.length > 8000) {
                    console.warn('Command line is very long, may cause ENAMETOOLONG error');
                }

                // Spawn Python process
                const pythonProcess = spawn('python', args, {
                    cwd: path.join(__dirname, '..'),
                    stdio: ['pipe', 'pipe', 'pipe']
                });

                // Handle spawn errors (like ENAMETOOLONG)
                pythonProcess.on('error', (error) => {
                    console.error('Python process spawn error:', error);
                    this.cleanupTempFiles(sequenceData, gffData);
                    reject(new Error(`Script execution error: ${error.code || error.message}`));
                });

                let stdoutChunks = [];
                let stdoutSize = 0;
                let stderr = '';
                const MAX_BUFFER_SIZE = 100 * 1024 * 1024; // 100MB limit

                pythonProcess.stdout.on('data', (data) => {
                    stdoutSize += data.length;
                    
                    // Check if we're exceeding safe buffer limits
                    if (stdoutSize > MAX_BUFFER_SIZE) {
                        console.error('Output too large, terminating process');
                        pythonProcess.kill('SIGTERM');
                        reject(new Error(`Output too large (>${MAX_BUFFER_SIZE / (1024*1024)}MB). The genome file produced too much data. Try using a smaller file or enable gene filtering.`));
                        return;
                    }
                    
                    stdoutChunks.push(data);
                });

                pythonProcess.stderr.on('data', (data) => {
                    const stderrChunk = data.toString();
                    stderr += stderrChunk;
                    
                    // Show progress from stderr in real-time
                    if (stderrChunk.includes('Processed') || stderrChunk.includes('Loaded sequence') || stderrChunk.includes('Found')) {
                        console.log('Progress:', stderrChunk.trim());
                    }
                });

                pythonProcess.on('close', (code) => {
                    console.log(`Python process exited with code: ${code}`);
                    console.log('Stdr output:', stderr);
                    
                    // Cleanup temporary files
                    this.cleanupTempFiles(sequenceData, gffData);
                    
                    if (code === 0) {
                        try {
                            // Combine chunks efficiently
                            const stdout = Buffer.concat(stdoutChunks).toString();
                            
                            if (!stdout.trim()) {
                                reject(new Error('Python script returned empty output'));
                                return;
                            }
                            
                            console.log(`Output size: ${(stdout.length / (1024*1024)).toFixed(2)}MB`);
                            
                            const result = JSON.parse(stdout);
                            console.log('Python result parsed successfully');
                            console.log(`Found ${result.genes?.length || 0} genes`);
                            resolve(result);
                        } catch (parseError) {
                            console.error(' Failed to fuck:', parseError);
                            reject(new Error(`Failed to parse Python output: ${parseError.message}`));
                        }
                    } else {
                        console.error('Python script failed:', stderr);
                        const stdout = Buffer.concat(stdoutChunks).toString();
                        reject(new Error(`Python script failed (exit code ${code}): ${stderr || stdout || 'Unknown error'}`));
                    }
                });

                pythonProcess.on('error', (error) => {
                    console.error('Failed to start Python process:', error);
                    reject(new Error(`Failed to start Python process: ${error.message}`));
                });

            } catch (error) {
                console.error('Error in callPythonScript:', error);
                reject(new Error(`Script execution error: ${error.message}`));
            }
        });
    }

    displayResults(genes, statistics, pagination = null) {
        // DON'T overwrite this.extractedGenes - it contains ALL genes!
        // Only store the current page results
        this.currentResults = { genes, statistics, pagination };
        
        // Hide placeholder
        const placeholder = document.getElementById('gene-results-placeholder');
        if (placeholder) {
            placeholder.style.display = 'none';
        }

        // Show results content
        const resultsContent = document.getElementById('gene-results-content');
        if (resultsContent) {
            resultsContent.classList.remove('hidden');
            resultsContent.innerHTML = '';

            // Add pagination info and Download All button
            this.createPaginationHeader(resultsContent, pagination);

            // Show all genes from current page (server-side pagination)
            // No need to slice - Python already sent us the right genes for this page
            genes.forEach((gene, index) => {
                const globalIndex = pagination ? (pagination.start_gene - 1 + index) : index;
                const geneCard = this.createGeneCard(gene, globalIndex);
                resultsContent.appendChild(geneCard);
            });

            // Add pagination controls
            this.createPaginationControls(resultsContent, pagination);
        }

        // Update statistics
        this.updateStatistics(statistics);

        // Enable download button
        const downloadBtn = document.getElementById('gene-download-btn');
        if (downloadBtn) {
            downloadBtn.disabled = false;
            downloadBtn.classList.remove('opacity-50', 'cursor-not-allowed');
        }
    }

    createPaginationHeader(container, pagination) {
        const header = document.createElement('div');
        header.className = 'bg-slate-800/50 border border-slate-600/50 rounded-lg p-4 mb-4';
        
        const totalGenes = pagination ? pagination.total_genes : this.extractedGenes.length;
        const totalPages = pagination ? pagination.total_pages : Math.ceil(totalGenes / this.genesPerPage);
        
        header.innerHTML = `
            <div class="flex items-center justify-between">
                <div class="flex items-center gap-4">
                    <div class="flex items-center gap-2">
                        <img src="emojis/dna.svg" class="w-5 h-5">
                        <h3 class="text-lg font-bold text-white">Extracted Genes</h3>
                    </div>
                    <div class="text-sm text-gray-300">
                        <span class="font-medium">${totalGenes.toLocaleString()}</span> genes found
                    </div>
                </div>
                <div class="flex items-center gap-3">
                    <div class="text-sm text-gray-400">
                        Page ${this.currentPage} of ${totalPages.toLocaleString()}
                    </div>
                    <button id="download-all-genes-btn" class="px-3 py-1.5 bg-green-600/20 border border-green-500/30 rounded-lg text-green-300 hover:bg-green-600/30 transition-all flex items-center gap-2 text-sm font-medium">
                        <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"></path>
                        </svg>
                        Download All ${totalGenes.toLocaleString()} Genes
                    </button>
                </div>
            </div>
        `;
        
        container.appendChild(header);
        
        // Bind Download All button
        const downloadAllBtn = header.querySelector('#download-all-genes-btn');
        if (downloadAllBtn) {
            downloadAllBtn.addEventListener('click', () => {
                // Add loading state
                downloadAllBtn.disabled = true;
                downloadAllBtn.innerHTML = `
                    <svg class="w-4 h-4 animate-spin" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15"></path>
                    </svg>
                    Preparing Download...
                `;
                
                this.downloadAllGenes().finally(() => {
                    // Reset button state
                    downloadAllBtn.disabled = false;
                    downloadAllBtn.innerHTML = `
                        <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"></path>
                        </svg>
                        Download All ${totalGenes.toLocaleString()} Genes
                    `;
                });
            });
        }
    }

    createPaginationControls(container, pagination) {
        const totalGenes = pagination ? pagination.total_genes : this.extractedGenes.length;
        const totalPages = pagination ? pagination.total_pages : Math.ceil(totalGenes / this.genesPerPage);
        
        if (totalPages <= 1) return; // No pagination needed
        
        const controls = document.createElement('div');
        controls.className = 'bg-slate-800/50 border border-slate-600/50 rounded-lg p-4 mt-4';
        
        const startGene = ((this.currentPage - 1) * this.genesPerPage) + 1;
        const endGene = Math.min(this.currentPage * this.genesPerPage, totalGenes);
        
        controls.innerHTML = `
            <div class="flex items-center justify-between">
                <div class="text-sm text-gray-400">
                    Showing genes ${startGene.toLocaleString()}-${endGene.toLocaleString()} of ${totalGenes.toLocaleString()}
                </div>
                <div class="flex items-center gap-2">
                    <button id="prev-page-btn" class="px-3 py-1 bg-slate-600/20 border border-slate-500/30 rounded text-slate-300 hover:bg-slate-600/30 transition-all disabled:opacity-50 disabled:cursor-not-allowed" ${this.currentPage === 1 ? 'disabled' : ''}>
                        Previous
                    </button>
                    <div class="flex items-center gap-1">
                        ${this.generatePageNumbers(totalPages)}
                    </div>
                    <button id="next-page-btn" class="px-3 py-1 bg-slate-600/20 border border-slate-500/30 rounded text-slate-300 hover:bg-slate-600/30 transition-all disabled:opacity-50 disabled:cursor-not-allowed" ${this.currentPage === totalPages ? 'disabled' : ''}>
                        Next
                    </button>
                </div>
            </div>
        `;
        
        container.appendChild(controls);
        
        // Bind pagination events
        const prevBtn = controls.querySelector('#prev-page-btn');
        const nextBtn = controls.querySelector('#next-page-btn');
        
        if (prevBtn) {
            prevBtn.addEventListener('click', () => {
                console.log('Previous button clicked');
                this.goToPage(this.currentPage - 1);
            });
        }
        
        if (nextBtn) {
            nextBtn.addEventListener('click', () => {
                console.log('Next button clicked');
                this.goToPage(this.currentPage + 1);
            });
        }
        
        // Bind page number buttons
        controls.querySelectorAll('.page-number-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const page = parseInt(e.target.dataset.page);
                console.log(`Page number ${page} button clicked`);
                this.goToPage(page);
            });
        });
    }

    generatePageNumbers(totalPages) {
        const current = this.currentPage;
        const delta = 2; // Show 2 pages before and after current
        const range = [];
        const rangeWithDots = [];
        
        for (let i = Math.max(2, current - delta); i <= Math.min(totalPages - 1, current + delta); i++) {
            range.push(i);
        }
        
        if (current - delta > 2) {
            rangeWithDots.push(1, '...');
        } else {
            rangeWithDots.push(1);
        }
        
        rangeWithDots.push(...range);
        
        if (current + delta < totalPages - 1) {
            rangeWithDots.push('...', totalPages);
        } else {
            rangeWithDots.push(totalPages);
        }
        
        return rangeWithDots.map(page => {
            if (page === '...') {
                return '<span class="px-2 py-1 text-gray-500">...</span>';
            } else {
                const isActive = page === current;
                return `<button class="page-number-btn px-2 py-1 rounded text-sm ${isActive ? 'bg-blue-600/30 border border-blue-500/50 text-blue-300' : 'text-gray-400 hover:text-white hover:bg-slate-600/30'} transition-all" data-page="${page}">${page}</button>`;
            }
        }).join('');
    }

    goToPage(page) {
        console.log(`goToPage called with page: ${page}`);
        
        // Since all genes are loaded in memory, just do client-side pagination
        const totalGenes = this.extractedGenes.length;
        const totalPages = Math.ceil(totalGenes / this.genesPerPage);
        
        console.log(`Total genes: ${totalGenes}, Total pages: ${totalPages}, Current page: ${this.currentPage}`);
        
        if (page < 1 || page > totalPages) {
            console.log(`Wongpage: ${page} (valid range: 1-${totalPages})`);
            return;
        }
        
        this.currentPage = page;
        
        // Calculate which genes to show for this page
        const startIdx = (page - 1) * this.genesPerPage;
        const endIdx = Math.min(startIdx + this.genesPerPage, totalGenes);
        const pageGenes = this.extractedGenes.slice(startIdx, endIdx);
        
        // Create fake pagination object for display
        const fakePagination = {
            total_genes: totalGenes,
            genes_per_page: this.genesPerPage,
            total_pages: totalPages,
            current_page: page,
            showing_genes: pageGenes.length,
            start_gene: startIdx + 1,
            end_gene: endIdx,
            has_next: page < totalPages,
            has_prev: page > 1
        };
        
        // Update display with current page genes
        this.displayResults(pageGenes, this.currentResults.statistics, fakePagination);
        
        // Show quick success message
        this.showNotification(
            `📄 Page ${page} - Showing genes ${startIdx + 1}-${endIdx} of ${totalGenes.toLocaleString()}`, 
            'success', 
            2000
        );
    }

    createGeneCard(gene, index) {
        const card = document.createElement('div');
        card.className = 'bg-slate-700/30 border border-slate-600/50 rounded-lg p-4 hover:border-blue-500/50 transition-all';
        
        card.innerHTML = `
            <div class="flex items-center justify-between mb-3">
                <div class="flex items-center gap-3">
                    <img src="emojis/dna.svg" class="w-5 h-5">
                    <h4 class="text-white font-medium">${gene.id || gene.name || `Gene_${index + 1}`}</h4>
                </div>
                <button class="extract-gene-btn px-3 py-1 bg-blue-600/20 border border-blue-500/30 rounded text-blue-300 hover:bg-blue-600/30 transition-all text-sm" data-gene-index="${index}">
                    Extract
                </button>
            </div>
            <div class="grid grid-cols-2 gap-4 text-sm">
                <div>
                    <span class="text-gray-400">Coordinates:</span>
                    <span class="text-white font-mono ml-2">${gene.start}-${gene.end}</span>
                </div>
                <div>
                    <span class="text-gray-400">Length:</span>
                    <span class="text-white font-mono ml-2">${gene.length} bp</span>
                </div>
                <div>
                    <span class="text-gray-400">Strand:</span>
                    <span class="text-white font-mono ml-2">${gene.strand}</span>
                </div>
                <div>
                    <span class="text-gray-400">Type:</span>
                    <span class="text-white font-mono ml-2">${gene.feature_type || 'gene'}</span>
                </div>
            </div>
        `;

        // Add extract button event
        const extractBtn = card.querySelector('.extract-gene-btn');
        extractBtn.addEventListener('click', () => this.extractSingleGene(index));

        return card;
    }

    updateStatistics(stats) {
        // Show stats panel
        const statsPanel = document.getElementById('gene-stats-panel');
        if (statsPanel) {
            statsPanel.classList.remove('hidden');
        }

        // Update individual stats
        const updates = [
            { id: 'gene-total-count', value: stats.total_genes },
            { id: 'gene-avg-length', value: `${stats.average_length} bp` },
            { id: 'gene-max-length', value: `${stats.longest_gene} bp` },
            { id: 'gene-min-length-stat', value: `${stats.shortest_gene} bp` }
        ];

        updates.forEach(({ id, value }) => {
            const element = document.getElementById(id);
            if (element) {
                element.textContent = value;
            }
        });
    }

    async extractSingleGene(geneIndex) {
        try {
            const gene = this.extractedGenes[geneIndex];
            if (!gene) return;

            // Create FASTA content for single gene
            const fastaContent = `>${gene.id || gene.name || `Gene_${geneIndex + 1}`} | ${gene.start}-${gene.end} | ${gene.strand} | ${gene.length}bp\n${gene.sequence}`;

            // Download the file
            this.downloadFile(fastaContent, `${gene.id || gene.name || `Gene_${geneIndex + 1}`}.fasta`);

            this.showNotification(`Gene "${gene.id || gene.name || `Gene_${geneIndex + 1}`}" extracted successfully!`, 'success');

        } catch (error) {
            console.error('Single gene extraction error:', error);
            this.showNotification('Failed to extract gene', 'error');
        }
    }

    async downloadAllGenes() {
        try {
            if (this.extractedGenes.length === 0) {
                this.showNotification('No genes to download', 'error');
                return;
            }

            const totalGenes = this.extractedGenes.length;
            
            // Show progress notification for large downloads
            if (totalGenes > 1000) {
                this.showNotification(`Preparing download of ${totalGenes.toLocaleString()} genes... This may take a moment.`, 'info', 5000);
            }

            // Create multi-FASTA content efficiently for large datasets
            const fastaChunks = [];
            const chunkSize = 1000; // Process 1000 genes at a time
            
            for (let i = 0; i < totalGenes; i += chunkSize) {
                const chunk = this.extractedGenes.slice(i, i + chunkSize);
                const chunkContent = chunk.map((gene, index) => {
                    const globalIndex = i + index;
                    return `>${gene.id || gene.name || `Gene_${globalIndex + 1}`} | ${gene.start}-${gene.end} | ${gene.strand} | ${gene.length}bp\n${gene.sequence}`;
                }).join('\n\n');
                
                fastaChunks.push(chunkContent);
                
                // Show progress for very large datasets
                if (totalGenes > 10000 && (i + chunkSize) % 10000 === 0) {
                    console.log(`Processed ${i + chunkSize} / ${totalGenes} genes...`);
                }
            }

            const fastaContent = fastaChunks.join('\n\n');
            
            // Generate filename with gene count
            const filename = `AnthuriumAmnicola_${totalGenes.toLocaleString()}_genes.fasta`;
            
            // Download the file
            this.downloadFile(fastaContent, filename);
            
            this.showNotification(`Successfully downloaded ${totalGenes.toLocaleString()} genes as ${filename}!`, 'success', 6000);

        } catch (error) {
            console.error('Download all genes error:', error);
            this.showNotification('Failed to download genes. File may be too large for browser memory.', 'error');
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

    clearResults() {
        // Reset file uploads
        this.fastaFile = null;
        this.gffFile = null;

        // Clear manual input
        const manualInput = document.getElementById('gene-manual-sequence');
        if (manualInput) {
            manualInput.value = '';
        }

        // Reset upload areas
        this.resetUploadArea('gene-fasta-upload-area', 'Upload FASTA/GenBank File');
        this.resetUploadArea('gene-gff-upload-area', 'Upload GFF/GFF3 Annotation File');

        // Hide results and show placeholder
        const resultsContent = document.getElementById('gene-results-content');
        if (resultsContent) {
            resultsContent.classList.add('hidden');
            resultsContent.innerHTML = '';
        }

        const placeholder = document.getElementById('gene-results-placeholder');
        if (placeholder) {
            placeholder.style.display = 'flex';
        }

        // Hide stats panel
        const statsPanel = document.getElementById('gene-stats-panel');
        if (statsPanel) {
            statsPanel.classList.add('hidden');
        }

        // Disable download button
        const downloadBtn = document.getElementById('gene-download-btn');
        if (downloadBtn) {
            downloadBtn.disabled = true;
            downloadBtn.classList.add('opacity-50', 'cursor-not-allowed');
        }

        // Reset extracted genes
        this.extractedGenes = [];

        // Validate input state
        this.validateInput();

        this.showNotification('Results cleared', 'success');
    }

    resetUploadArea(areaId, originalText) {
        const uploadArea = document.getElementById(areaId);
        if (!uploadArea) return;

        uploadArea.classList.remove('border-green-500/50', 'bg-green-500/10');
        uploadArea.classList.add('border-slate-600');

        const textElement = uploadArea.querySelector('p');
        if (textElement) {
            textElement.textContent = originalText;
            textElement.classList.remove('text-green-400');
        }
    }

    setLoadingState(isLoading) {
        const analyzeBtn = document.getElementById('gene-analyze-btn');
        if (analyzeBtn) {
            if (isLoading) {
                analyzeBtn.disabled = true;
                analyzeBtn.innerHTML = `
                    <img src="emojis/microscope.svg" class="w-5 h-5">
                    🔄 Analyzing Genome...
                `;
                analyzeBtn.classList.add('opacity-50');
            } else {
                analyzeBtn.disabled = false;
                analyzeBtn.innerHTML = `
                    <img src="emojis/microscope.svg" class="w-5 h-5">
                    Analyze & Extract Genes
                `;
                analyzeBtn.classList.remove('opacity-50');
                this.validateInput(); // Re-validate in case there's no input
            }
        }
    }

    showGFFWarning() {
        return new Promise((resolve) => {
            // Create overlay
            const overlay = document.createElement('div');
            overlay.className = 'fixed inset-0 bg-black/80 backdrop-blur-sm z-50 flex items-center justify-center p-4';
            
            // Create warning modal
            const modal = document.createElement('div');
            modal.className = 'bg-slate-800 border-2 border-yellow-500 rounded-xl p-6 max-w-2xl w-full mx-auto shadow-2xl';
            
            modal.innerHTML = `
                <div class="flex items-center gap-3 mb-4">
                    <img src="emojis/warning.svg" class="w-8 h-8">
                    <h3 class="text-xl font-bold text-yellow-400">Missing GFF/GFF3 Annotation File</h3>
                </div>
                <div class="space-y-4 mb-6">
                    <p class="text-gray-300">
                        🧬 <strong>For REAL gene extraction</strong>, you need a GFF/GFF3 annotation file! Without annotations, the system will fall back to ORF detection which finds <em>potential</em> genes, not actual annotated genes.
                    </p>
                    <div class="bg-red-500/10 border border-red-500/30 rounded-lg p-4">
                        <h4 class="text-red-400 font-medium mb-2">🚨 REAL vs PREDICTED Genes:</h4>
                        <div class="grid grid-cols-2 gap-4 text-sm">
                            <div>
                                <h5 class="text-green-400 font-medium">✅ WITH GFF (Real Genes):</h5>
                                <ul class="text-gray-300 space-y-1">
                                    <li>• Actual annotated genes</li>
                                    <li>• Precise gene boundaries</li>
                                    <li>• Gene names & functions</li>
                                    <li>• Strand information</li>
                                </ul>
                            </div>
                            <div>
                                <h5 class="text-yellow-400 font-medium">⚠️ WITHOUT GFF (ORF Detection):</h5>
                                <ul class="text-gray-300 space-y-1">
                                    <li>• Potential genes only</li>
                                    <li>• May include false positives</li>
                                    <li>• No functional annotation</li>
                                    <li>• Less accurate boundaries</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                    <p class="text-gray-400 text-sm">
                        💡 <strong>Tip:</strong> For your AnthuriumAmnicola genome, use the corresponding .gff file for real gene extraction!
                    </p>
                </div>
                <div class="flex gap-3 justify-end">
                    <button id="gff-warning-cancel" class="px-4 py-2 bg-slate-600/20 border border-slate-500/30 rounded-lg text-slate-300 hover:bg-slate-600/30 transition-all">
                        Cancel & Upload GFF
                    </button>
                    <button id="gff-warning-continue" class="px-4 py-2 bg-yellow-600/20 border border-yellow-500/30 rounded-lg text-yellow-300 hover:bg-yellow-600/30 transition-all">
                        Continue Anyway
                    </button>
                </div>
            `;

            overlay.appendChild(modal);
            document.body.appendChild(overlay);

            // Handle button clicks
            const cancelBtn = modal.querySelector('#gff-warning-cancel');
            const continueBtn = modal.querySelector('#gff-warning-continue');

            cancelBtn.addEventListener('click', () => {
                document.body.removeChild(overlay);
                resolve(false);
            });

            continueBtn.addEventListener('click', () => {
                document.body.removeChild(overlay);
                resolve(true);
            });

            // Close on overlay click
            overlay.addEventListener('click', (e) => {
                if (e.target === overlay) {
                    document.body.removeChild(overlay);
                    resolve(false);
                }
            });
        });
    }

    showLargeFileWarning(file, fileType, maxRecommendedSize) {
        return new Promise((resolve) => {
            // Create overlay
            const overlay = document.createElement('div');
            overlay.className = 'fixed inset-0 bg-black/80 backdrop-blur-sm z-50 flex items-center justify-center p-4';
            
            // Create warning modal
            const modal = document.createElement('div');
            modal.className = 'bg-slate-800 border-2 border-orange-500 rounded-xl p-6 max-w-2xl w-full mx-auto shadow-2xl';
            
            const fileSize = this.formatFileSize(file.size);
            const recommendedSize = this.formatFileSize(maxRecommendedSize);
            
            modal.innerHTML = `
                <div class="flex items-center gap-3 mb-4">
                    <img src="emojis/warning.svg" class="w-8 h-8">
                    <h3 class="text-xl font-bold text-orange-400">Large File Warning</h3>
                </div>
                <div class="space-y-4 mb-6">
                    <p class="text-gray-300">
                        You're about to process a large ${fileType}:
                    </p>
                    <div class="bg-orange-500/10 border border-orange-500/30 rounded-lg p-4">
                        <div class="grid grid-cols-2 gap-4 text-sm">
                            <div>
                                <span class="text-orange-400 font-medium">File Size:</span>
                                <span class="text-white ml-2">${fileSize}</span>
                            </div>
                            <div>
                                <span class="text-orange-400 font-medium">Recommended:</span>
                                <span class="text-white ml-2">< ${recommendedSize}</span>
                            </div>
                        </div>
                    </div>
                    <div class="bg-red-900/20 border border-red-500/30 rounded-lg p-4">
                        <h4 class="text-red-300 font-medium mb-2">⚠️ FULL GENOME ANALYSIS - RAM WARNING:</h4>
                        <ul class="text-red-200 text-sm space-y-2">
                            <li><strong>• RAM Usage: ~${fileSize} of RAM will be used</strong></li>
                            <li>• No chunking - entire genome loaded for complete gene extraction</li>
                            <li>• Processing time: Several minutes for multi-GB files</li>
                            <li>• May cause system slowdown if insufficient RAM</li>
                            <li>• Close other applications to free up memory</li>
                        </ul>
                    </div>
                    <div class="bg-blue-900/20 border border-blue-500/30 rounded-lg p-4">
                        <h4 class="text-blue-300 font-medium mb-2">✅ Benefits of Full Analysis:</h4>
                        <ul class="text-blue-200 text-sm space-y-1">
                            <li>• Complete genes (no artificial cuts)</li>
                            <li>• Accurate gene boundaries</li>
                            <li>• Full genome coverage</li>
                            <li>• Real gene sequences (not fragments)</li>
                        </ul>
                    </div>
                </div>
                <div class="flex gap-3 justify-end">
                    <button id="large-file-cancel" class="px-4 py-2 bg-slate-600/20 border border-slate-500/30 rounded-lg text-slate-300 hover:bg-slate-600/30 transition-all">
                        Cancel
                    </button>
                    <button id="large-file-continue" class="px-4 py-2 bg-green-600/20 border border-green-500/30 rounded-lg text-green-300 hover:bg-green-600/30 transition-all">
                        🧬 Analyze Full Genome
                    </button>
                </div>
            `;

            overlay.appendChild(modal);
            document.body.appendChild(overlay);

            // Handle button clicks
            const cancelBtn = modal.querySelector('#large-file-cancel');
            const continueBtn = modal.querySelector('#large-file-continue');

            cancelBtn.addEventListener('click', () => {
                document.body.removeChild(overlay);
                resolve(false);
            });

            continueBtn.addEventListener('click', () => {
                document.body.removeChild(overlay);
                resolve(true);
            });
        });
    }

    readFileContent(file) {
        return new Promise((resolve, reject) => {
            // For very large files (>10GB), show error
            if (file.size > 10 * 1024 * 1024 * 1024) { // > 10GB
                reject(new Error(`File "${file.name}" is too large (${this.formatFileSize(file.size)}). Maximum supported size is 10GB.`));
                return;
            }

            // For large files (>1GB), we need to stream to Python instead of loading into memory
            if (file.size > 1 * 1024 * 1024 * 1024) { // > 1GB
                this.streamLargeFileToPython(file, resolve, reject);
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

    async streamLargeFileToPython(file, resolve, reject) {
        try {
            const fileSizeGB = (file.size / (1024 * 1024 * 1024)).toFixed(2);
            const fileSizeMB = (file.size / (1024 * 1024)).toFixed(1);
            
            // Show RAM usage warning
            console.log(`FULL GENOME ANALYSIS MODE`);
            console.log(`File: ${file.name}`);
            console.log(`Size: ${fileSizeGB} GB (${fileSizeMB} MB)`);
            console.log(`RAM USAGE WARNING: This will use ~${fileSizeGB} GB of RAM`);
            console.log(`No chunking - analyzing entire genome for complete gene extraction`);
            
            // Show user notification about RAM usage
            this.showNotification(
                `RAM Warning: This ${fileSizeGB}GB file will use ~${fileSizeGB}GB of RAM for complete genome analysis. Make sure you have enough free memory!`, 
                'warning', 
                8000
            );
            
            // For large files, we'll save to a temporary file and pass the path to Python
            const fs = require('fs');
            const path = require('path');
            const os = require('os');
            
            // Create temporary file
            const tempDir = os.tmpdir();
            // Truncate filename if too long to avoid ENAMETOOLONG
            const safeName = file.name.length > 50 ? file.name.substring(0, 50) + '...' + path.extname(file.name) : file.name;
            const tempFileName = `gene_extractor_${Date.now()}_${safeName}`;
            const tempFilePath = path.join(tempDir, tempFileName);
            
            console.log(`Creating temporary copy: ${tempFilePath}`);
            console.log(`This creates an EXACT COPY of your file (no cuts/chunks)`);
            console.log(`Python will then read the ENTIRE temp file into RAM for complete analysis`);
            
            // Stream file to disk
            await this.streamFileToDisk(file, tempFilePath);
            
            // Return the temp file path instead of content
            resolve(`TEMP_FILE:${tempFilePath}`);
            
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

    cleanupTempFiles(sequenceData, gffData) {
        try {
            const fs = require('fs');
            
            // Clean up sequence temp file
            if (sequenceData && sequenceData.startsWith('TEMP_FILE:')) {
                const tempFilePath = sequenceData.replace('TEMP_FILE:', '');
                if (fs.existsSync(tempFilePath)) {
                    fs.unlinkSync(tempFilePath);
                    console.log(`Cleaned up temp file: ${tempFilePath}`);
                }
            }
            
            // Clean up GFF temp file
            if (gffData && gffData.startsWith('TEMP_FILE:')) {
                const tempFilePath = gffData.replace('TEMP_FILE:', '');
                if (fs.existsSync(tempFilePath)) {
                    fs.unlinkSync(tempFilePath);
                    console.log(`Cleaned up temp file: ${tempFilePath}`);
                }
            }
        } catch (error) {
            console.warn('Failed to cleanup temp file ok:', error.message);
        }
    }

    formatFileSize(bytes) {
        if (bytes === 0) return '0 Bytes';
        const k = 1024;
        const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));
        return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
    }

    showNotification(message, type = 'info') {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `fixed top-4 right-4 z-50 px-6 py-3 rounded-lg shadow-lg transition-all duration-300 transform translate-x-full`;
        
        // Set colors based on type
        switch (type) {
            case 'success':
                notification.className += ' bg-green-600/90 border border-green-500/50 text-green-100';
                break;
            case 'error':
                notification.className += ' bg-red-600/90 border border-red-500/50 text-red-100';
                break;
            case 'warning':
                notification.className += ' bg-yellow-600/90 border border-yellow-500/50 text-yellow-100';
                break;
            default:
                notification.className += ' bg-blue-600/90 border border-blue-500/50 text-blue-100';
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
                if (document.body.contains(notification)) {
                    document.body.removeChild(notification);
                }
            }, 300);
        }, 3000);
    }
}

// Initialize when page loads
document.addEventListener('DOMContentLoaded', () => {
    if (document.getElementById('gene-extractor-page')) {
        window.geneExtractor = new PythonGeneExtractor();
    }
});