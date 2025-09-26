// SNP Highlighter JavaScript
document.addEventListener('DOMContentLoaded', function() {
    // Elements
    const referenceFileUpload = document.getElementById('snp-reference-file-upload');
    const sampleFileUpload = document.getElementById('snp-sample-file-upload');
    const referenceInput = document.getElementById('snp-reference-input');
    const sampleInput = document.getElementById('snp-sample-input');
    const analyzeBtn = document.getElementById('analyze-snp-btn');
    
    // Stats elements
    const referenceLength = document.getElementById('snp-reference-length');
    const sampleLength = document.getElementById('snp-sample-length');
    const referenceGC = document.getElementById('snp-reference-gc');
    const sampleGC = document.getElementById('snp-sample-gc');
    
    // Clear buttons
    const referenceClearBtn = document.getElementById('snp-reference-clear-btn');
    const sampleClearBtn = document.getElementById('snp-sample-clear-btn');
    
    // Output elements
    const referenceResultsContent = document.getElementById('reference-results-content');
    const sampleResultsContent = document.getElementById('sample-results-content');
    
    // Load buttons
    const referenceLoadBtn = document.getElementById('reference-load-btn');
    const sampleLoadBtn = document.getElementById('sample-load-btn');
    
    // Edit Lock buttons
    const referenceEditLockBtn = document.getElementById('reference-edit-lock-mode');
    const sampleEditLockBtn = document.getElementById('sample-edit-lock-mode');
    
    // Results indicators
    const snpOutputCount = document.getElementById('snp-output-count');
    const snpOutputSimilarity = document.getElementById('snp-output-similarity');
    
    // Copy and download buttons
    const copySnpResults = document.getElementById('copy-snp-results');
    const downloadSnpResults = document.getElementById('download-snp-results');
    
    console.log('SNP Copy button:', copySnpResults);
    console.log('SNP Download button:', downloadSnpResults);

    // File upload handlers
    if (referenceFileUpload) {
        referenceFileUpload.addEventListener('change', function(e) {
            const file = e.target.files[0];
            if (file) {
                const reader = new FileReader();
                reader.onload = function(e) {
                    let content = e.target.result;
                    
                    // Parse different file formats
                    if (file.name.toLowerCase().includes('.fasta') || file.name.toLowerCase().includes('.fa')) {
                        content = parseFasta(content);
                    } else if (file.name.toLowerCase().includes('.gb') || file.name.toLowerCase().includes('.genbank')) {
                        content = parseGenBank(content);
                    }
                    
                    // Populate the manual input textarea
                    if (referenceInput) {
                        referenceInput.value = content.toUpperCase();
                        updateSequenceStats('reference');
                    }
                    
                    showNotification(`Reference file "${file.name}" loaded successfully!`, 'success');
                };
                reader.readAsText(file);
            }
        });
    }

    if (sampleFileUpload) {
        sampleFileUpload.addEventListener('change', function(e) {
            const file = e.target.files[0];
            if (file) {
                const reader = new FileReader();
                reader.onload = function(e) {
                    let content = e.target.result;
                    
                    // Parse different file formats
                    if (file.name.toLowerCase().includes('.fasta') || file.name.toLowerCase().includes('.fa')) {
                        content = parseFasta(content);
                    } else if (file.name.toLowerCase().includes('.gb') || file.name.toLowerCase().includes('.genbank')) {
                        content = parseGenBank(content);
                    }
                    
                    // Populate the manual input textarea
                    if (sampleInput) {
                        sampleInput.value = content.toUpperCase();
                        updateSequenceStats('sample');
                    }
                    
                    showNotification(`Sample file "${file.name}" loaded successfully!`, 'success');
                };
                reader.readAsText(file);
            }
        });
    }

    // DNA input handlers
    if (referenceInput) {
        referenceInput.addEventListener('input', function() {
            // Only allow valid DNA bases
            this.value = this.value.toUpperCase().replace(/[^ATGCNRYSWKMBDHV]/g, '');
            updateSequenceStats('reference');
        });
    }

    if (sampleInput) {
        sampleInput.addEventListener('input', function() {
            // Only allow valid DNA bases
            this.value = this.value.toUpperCase().replace(/[^ATGCNRYSWKMBDHV]/g, '');
            updateSequenceStats('sample');
        });
    }

    // Clear button handlers
    if (referenceClearBtn) {
        referenceClearBtn.addEventListener('click', function() {
            if (referenceInput) {
                referenceInput.value = '';
                updateSequenceStats('reference');
            }
            // Clear the results content too
            if (referenceResultsContent) {
                referenceResultsContent.innerHTML = `
                    <div class="text-center text-gray-500 mt-20">
                        <img src="emojis/bookmark.svg" class="w-12 h-12 mx-auto mb-3 opacity-50">
                        <p class="text-base">Reference Analysis</p>
                        <p class="text-xs mt-2">Reference sequence results will appear here</p>
                    </div>
                `;
            }
            showNotification('Reference sequence cleared', 'info');
        });
    }

    if (sampleClearBtn) {
        sampleClearBtn.addEventListener('click', function() {
            if (sampleInput) {
                sampleInput.value = '';
                updateSequenceStats('sample');
            }
            // Clear the results content too
            if (sampleResultsContent) {
                sampleResultsContent.innerHTML = `
                    <div class="text-center text-gray-500 mt-20">
                        <img src="emojis/test_tube.svg" class="w-12 h-12 mx-auto mb-3 opacity-50">
                        <p class="text-base">Sample Analysis</p>
                        <p class="text-xs mt-2">Sample sequence results will appear here</p>
                    </div>
                `;
            }
            showNotification('Sample sequence cleared', 'info');
        });
    }

    // Update sequence statistics
    function updateSequenceStats(type) {
        const input = type === 'reference' ? referenceInput : sampleInput;
        const lengthElement = type === 'reference' ? referenceLength : sampleLength;
        const gcElement = type === 'reference' ? referenceGC : sampleGC;
        
        if (!input || !lengthElement || !gcElement) return;
        
        const sequence = input.value.replace(/\s/g, '');
        const length = sequence.length;
        
        lengthElement.textContent = `Length: ${length} bp`;
        
        if (length > 0) {
            const gcCount = (sequence.match(/[GC]/g) || []).length;
            const gcPercent = ((gcCount / length) * 100).toFixed(1);
            gcElement.textContent = `${gcPercent}%`;
        } else {
            gcElement.textContent = '0%';
        }
    }

    // Reference Load button
    if (referenceLoadBtn) {
        referenceLoadBtn.addEventListener('click', function() {
            const sequence = referenceInput ? referenceInput.value.trim() : '';
            
            if (!sequence) {
                showNotification('Please enter a reference DNA sequence to load', 'error');
                return;
            }
            
            loadSequence(sequence, 'reference');
        });
    }

    // Sample Load button
    if (sampleLoadBtn) {
        sampleLoadBtn.addEventListener('click', function() {
            const sequence = sampleInput ? sampleInput.value.trim() : '';
            
            if (!sequence) {
                showNotification('Please enter a sample DNA sequence to load', 'error');
                return;
            }
            
            loadSequence(sequence, 'sample');
        });
    }

    // Edit Lock button functionality
    if (referenceEditLockBtn) {
        referenceEditLockBtn.addEventListener('click', function() {
            toggleEditLockMode('reference');
        });
    }

    if (sampleEditLockBtn) {
        sampleEditLockBtn.addEventListener('click', function() {
            toggleEditLockMode('sample');
        });
    }

    // Load sequence function
    async function loadSequence(sequence, type) {
        const loadBtn = type === 'reference' ? referenceLoadBtn : sampleLoadBtn;
        const resultsContent = type === 'reference' ? referenceResultsContent : sampleResultsContent;
        
        try {
            if (loadBtn) {
                loadBtn.disabled = true;
                loadBtn.textContent = '🔄 Loading...';
            }
            
            const result = await callPythonScript('load', sequence, type);
            
            if (result.success) {
                displayLoadedSequence(result.data, type);
                showNotification(`${type.charAt(0).toUpperCase() + type.slice(1)} sequence loaded successfully!`, 'success');
            } else {
                throw new Error(result.error || 'Failed to load sequence');
            }
            
        } catch (error) {
            console.error('Error:', error);
            showNotification(`Error loading ${type} sequence: ${error.message}`, 'error');
        } finally {
            if (loadBtn) {
                loadBtn.disabled = false;
                loadBtn.textContent = '📁 LOAD';
            }
        }
    }

    // Display loaded sequence in results area
    function displayLoadedSequence(data, type) {
        const resultsContent = type === 'reference' ? referenceResultsContent : sampleResultsContent;
        
        if (!resultsContent) return;
        
        // Show the loaded sequence with formatting
        if (data.formatted_sequence) {
            resultsContent.innerHTML = formatSequenceWithColors(data.formatted_sequence, type);
            
            // If we have stored SNP positions, apply highlighting to both reference and sample
            const storedSNPs = resultsContent.getAttribute('data-snp-positions');
            if (storedSNPs) {
                try {
                    const snps = JSON.parse(storedSNPs);
                    const currentSeq = data.sequence || data.formatted_sequence.replace(/[^ATGCN]/gi, '');
                    
                    // Apply SNP highlighting to the newly loaded sequence
                    setTimeout(() => {
                        if (type === 'sample') {
                            highlightSNPsInOutputSequences(snps, null, currentSeq);
                        } else if (type === 'reference') {
                            highlightSNPsInOutputSequences(snps, currentSeq, null);
                        }
                    }, 100); // Small delay to ensure DOM is updated
                } catch (error) {
                    console.error('Error applying stored SNP highlighting:', error);
                }
            }
        }
    }

    // Format sequence with colors and strand indicators
    function formatSequenceWithColors(formattedSequence, type) {
        const lines = formattedSequence.split('\n');
        let html = '';
        
        // Color scheme based on type
        const colorScheme = type === 'reference' ? 'blue' : 'green';
        
        for (const line of lines) {
            if (line.includes('|')) {
                const [lineNum, sequence] = line.split('|');
                
                // Apply DNA base coloring
                let coloredSequence = sequence.replace(/([ATGCNRYSWKMBDHV])/g, (match) => {
                    switch(match) {
                        case 'A': return `<span class="text-red-400">${match}</span>`;
                        case 'T': return `<span class="text-blue-400">${match}</span>`;
                        case 'G': return `<span class="text-green-400">${match}</span>`;
                        case 'C': return `<span class="text-yellow-400">${match}</span>`;
                        case 'N': return `<span class="text-gray-400">${match}</span>`;
                        // Ambiguous bases
                        case 'R': return `<span class="text-purple-400">${match}</span>`; // A or G
                        case 'Y': return `<span class="text-orange-400">${match}</span>`; // C or T
                        case 'S': return `<span class="text-cyan-400">${match}</span>`;   // G or C
                        case 'W': return `<span class="text-pink-400">${match}</span>`;   // A or T
                        case 'K': return `<span class="text-indigo-400">${match}</span>`; // G or T
                        case 'M': return `<span class="text-rose-400">${match}</span>`;   // A or C
                        case 'B': return `<span class="text-teal-400">${match}</span>`;   // C, G, or T
                        case 'D': return `<span class="text-lime-400">${match}</span>`;   // A, G, or T
                        case 'H': return `<span class="text-amber-400">${match}</span>`;  // A, C, or T
                        case 'V': return `<span class="text-emerald-400">${match}</span>`; // A, C, or G
                        default: return match;
                    }
                });
                
                // Format with line numbers and strand indicators
                html += `<div class="mb-1 w-full break-all flex">
                    <span class="text-gray-500 text-xs mr-2 flex-shrink-0 select-none" contenteditable="false">5'</span>
                    <span class="text-gray-400 mr-3 flex-shrink-0 font-mono text-sm select-none" contenteditable="false">${lineNum}:</span>
                    <span class="font-mono editable-sequence">${coloredSequence}</span>
                    <span class="text-gray-500 text-xs ml-2 flex-shrink-0 select-none" contenteditable="false">3'</span>
                </div>`;
            } else {
                // Fallback for lines without line numbers
                html += `<div class="mb-1 w-full break-all">${line}</div>`;
            }
        }
        
        return html;
    }

    // Toggle edit lock mode
    function toggleEditLockMode(type) {
        const editLockBtn = type === 'reference' ? referenceEditLockBtn : sampleEditLockBtn;
        const resultsContent = type === 'reference' ? referenceResultsContent : sampleResultsContent;
        
        if (!editLockBtn || !resultsContent) return;
        
        const editableElements = resultsContent.querySelectorAll('.editable-sequence');
        
        if (editLockBtn.textContent.includes('🔒')) {
            // Currently locked, unlock for editing
            editableElements.forEach(el => {
                el.contentEditable = true;
                el.style.backgroundColor = 'rgba(59, 130, 246, 0.1)'; // Blue highlight
                el.style.border = '1px dashed rgba(59, 130, 246, 0.3)';
                el.style.borderRadius = '2px';
            });
            editLockBtn.textContent = '🔓 Edit Mode';
            editLockBtn.className = editLockBtn.className.replace('orange-600', 'blue-600').replace('orange-500', 'blue-500').replace('orange-300', 'blue-300');
            showNotification(`${type.charAt(0).toUpperCase() + type.slice(1)} sequence is now editable`, 'info');
        } else {
            // Currently unlocked, lock editing
            editableElements.forEach(el => {
                el.contentEditable = false;
                el.style.backgroundColor = '';
                el.style.border = '';
                el.style.borderRadius = '';
            });
            editLockBtn.textContent = '🔒 Edit Lock';
            editLockBtn.className = editLockBtn.className.replace('blue-600', 'orange-600').replace('blue-500', 'orange-500').replace('blue-300', 'orange-300');
            showNotification(`${type.charAt(0).toUpperCase() + type.slice(1)} sequence is now locked`, 'info');
        }
    }

    // Parse FASTA format
    function parseFasta(content) {
        const lines = content.split('\n');
        let sequence = '';
        for (let line of lines) {
            if (!line.startsWith('>')) {
                sequence += line.trim();
            }
        }
        return sequence;
    }

    // Parse GenBank format
    function parseGenBank(content) {
        const lines = content.split('\n');
        let sequence = '';
        let inSequence = false;
        
        for (let line of lines) {
            if (line.startsWith('ORIGIN')) {
                inSequence = true;
                continue;
            }
            if (line.startsWith('//')) {
                break;
            }
            if (inSequence) {
                // Remove line numbers and spaces
                const cleanLine = line.replace(/^\s*\d+\s*/, '').replace(/\s/g, '');
                sequence += cleanLine;
            }
        }
        return sequence;
    }

    // Call Python script
    async function callPythonScript(operation, sequence, type = '') {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Python script
            const scriptPath = path.join(__dirname, 'assets', 'snp-highlighter.py');
            const args = [scriptPath, operation, sequence];
            
            // Add type argument if provided
            if (type) {
                args.push(type);
            }

            console.log('🐍 Calling SNP Highlighter Python script:', scriptPath);
            console.log('🔧 Arguments:', args);

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
                console.log(`🐍 Python process exited with code: ${code}`);
                
                if (code !== 0) {
                    console.error('❌ Python stderr:', stderr);
                    reject(new Error(`Python script failed with code ${code}: ${stderr}`));
                    return;
                }

                try {
                    console.log('📤 Python stdout:', stdout);
                    const result = JSON.parse(stdout);
                    resolve(result);
                } catch (error) {
                    console.error('❌ Failed to parse Python output:', stdout);
                    reject(new Error(`Failed to parse Python output: ${error.message}`));
                }
            });

            pythonProcess.on('error', (error) => {
                console.error('❌ Python process error:', error);
                reject(new Error(`Failed to start Python process: ${error.message}`));
            });
        });
    }

    // Show notification function
    function showNotification(message, type = 'info') {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `fixed top-4 right-4 z-50 px-4 py-3 rounded-lg shadow-lg transition-all duration-300 transform translate-x-full max-w-md`;
        
        // Set colors based on type
        switch(type) {
            case 'success':
                notification.className += ' bg-green-600/90 border border-green-500/50 text-green-100';
                break;
            case 'error':
                notification.className += ' bg-red-600/90 border border-red-500/50 text-red-100';
                break;
            case 'warning':
                notification.className += ' bg-orange-600/90 border border-orange-500/50 text-orange-100';
                break;
            case 'info':
            default:
                notification.className += ' bg-blue-600/90 border border-blue-500/50 text-blue-100';
                break;
        }
        
        notification.textContent = message;
        document.body.appendChild(notification);
        
        // Animate in
        setTimeout(() => {
            notification.style.transform = 'translateX(0)';
        }, 100);
        
        // Animate out and remove
        setTimeout(() => {
            notification.style.transform = 'translateX(100%)';
            setTimeout(() => {
                if (notification.parentNode) {
                    notification.parentNode.removeChild(notification);
                }
            }, 300);
        }, 3000);
    }

    // Initialize page state - ensure all result sections are hidden on load
    function initializePageState() {
        clearAllResults();
        
        // Reset progress bar
        updateProgressBar(0);
        updateQuickStatus('Ready');
        
        // Reset basic results
        if (quickTotalSNPs) quickTotalSNPs.textContent = '--';
        if (quickSimilarity) quickSimilarity.textContent = '--%';
        if (quickRefLength) quickRefLength.textContent = '-- bp';
        if (quickSampleLength) quickSampleLength.textContent = '-- bp';
    }

    // Analysis option elements
    const highlightSNPsCheckbox = document.getElementById('highlight-snps');
    const stopCodonDetectionCheckbox = document.getElementById('stop-codon-detection');
    const snpFrequencyCounterCheckbox = document.getElementById('snp-frequency-counter');
    const codingImpactAnalysisCheckbox = document.getElementById('coding-impact-analysis');

    // Result section elements
    const snpHighlightResults = document.getElementById('snp-highlight-results');
    const stopCodonResults = document.getElementById('stop-codon-results');
    const snpFrequencyResults = document.getElementById('snp-frequency-results');
    const synonymousResults = document.getElementById('synonymous-results');

    // Checkboxes only control what gets analyzed - no real-time show/hide
    // Results only appear after clicking "Analyze" button

    // Checkboxes control what gets analyzed - results only appear after clicking "Analyze"

    // Quick results elements
    const quickTotalSNPs = document.getElementById('quick-total-snps');
    const quickSimilarity = document.getElementById('quick-similarity');
    const quickRefLength = document.getElementById('quick-ref-length');
    const quickSampleLength = document.getElementById('quick-sample-length');
    const quickStatus = document.getElementById('quick-status');
    const progressBar = document.getElementById('progress-bar');
    const progressPercent = document.getElementById('progress-percent');

    // Analyze button
    if (analyzeBtn) {
        analyzeBtn.addEventListener('click', function() {
            const referenceSeq = referenceInput ? referenceInput.value.trim() : '';
            const sampleSeq = sampleInput ? sampleInput.value.trim() : '';
            
            if (!referenceSeq || !sampleSeq) {
                showNotification('Please enter both reference and sample DNA sequences', 'error');
                return;
            }
            
            // Get selected analysis options
            const analysisOptions = {
                highlightSNPs: highlightSNPsCheckbox ? highlightSNPsCheckbox.checked : false,
                stopCodonDetection: stopCodonDetectionCheckbox ? stopCodonDetectionCheckbox.checked : false,
                snpFrequencyCounter: snpFrequencyCounterCheckbox ? snpFrequencyCounterCheckbox.checked : false,
                codingImpactAnalysis: codingImpactAnalysisCheckbox ? codingImpactAnalysisCheckbox.checked : false
            };
            
            analyzeSNPs(referenceSeq, sampleSeq, analysisOptions);
        });
    }

    // Clear all previous results
    function clearAllResults() {
        // Hide all result sections
        if (snpHighlightResults) snpHighlightResults.style.display = 'none';
        if (stopCodonResults) stopCodonResults.style.display = 'none';
        if (snpFrequencyResults) snpFrequencyResults.style.display = 'none';
        if (synonymousResults) synonymousResults.style.display = 'none';
        
        // Hide detailed results container
        const detailedResults = document.getElementById('detailed-results');
        if (detailedResults) detailedResults.style.display = 'none';
        
        // Hide SNP breakdown
        const snpBreakdown = document.getElementById('snp-breakdown');
        if (snpBreakdown) snpBreakdown.style.display = 'none';
        
        // Clear mutation highlighting
        clearMutationHighlighting();
    }

    // Clear mutation highlighting from textareas
    function clearMutationHighlighting() {
        const referenceInput = document.getElementById('snp-reference-input');
        const sampleInput = document.getElementById('snp-sample-input');
        
        if (referenceInput) {
            referenceInput.style.borderColor = '';
            referenceInput.classList.remove('has-mutations');
        }
        
        if (sampleInput) {
            sampleInput.style.borderColor = '';
            sampleInput.classList.remove('has-mutations');
        }
        
        // Clear SNP highlighting classes
        clearSNPHighlighting();
        
        // Clear stored results
        window.lastCodingResults = null;
    }

    // Clear SNP highlighting classes
    function clearSNPHighlighting() {
        const sampleContent = document.getElementById('sample-results-content');
        const referenceContent = document.getElementById('reference-results-content');
        
        if (sampleContent) {
            sampleContent.classList.remove('snp-highlighting-active');
            sampleContent.removeAttribute('data-snp-positions');
        }
        
        if (referenceContent) {
            referenceContent.classList.remove('snp-highlighting-active');
            referenceContent.removeAttribute('data-snp-positions');
        }
        
        console.log('SNP highlighting classes cleared');
    }



    // Analyze SNPs function
    async function analyzeSNPs(referenceSeq, sampleSeq, options) {
        try {
            // Clear any previous results first
            clearAllResults();
            
            // Update UI to show analysis in progress
            updateProgressBar(0);
            updateQuickStatus('Analyzing sequences...');
            analyzeBtn.disabled = true;
            analyzeBtn.textContent = '🔄 Analyzing...';

            // Step 1: Basic SNP analysis
            updateProgressBar(25);
            const analysisResult = await callPythonScript('analyze', referenceSeq, sampleSeq);
            
            if (!analysisResult.success) {
                throw new Error(analysisResult.error || 'SNP analysis failed');
            }

            // Step 2: Update basic results
            updateProgressBar(50);
            updateBasicResults(analysisResult.data);

            // Step 3: Perform additional analyses based on selected options
            updateProgressBar(75);
            await performAdditionalAnalysesUpdated(referenceSeq, sampleSeq, analysisResult.data, options);
            
            // Step 4: Apply combined highlighting if both SNP and coding analysis are selected
            if (options.highlightSNPs && options.codingImpactAnalysis && window.lastCodingResults) {
                // Combine SNP and coding highlighting
                combinedHighlighting(analysisResult.data.snps, window.lastCodingResults, referenceSeq, sampleSeq);
            }

            // Step 4: Update UI with results
            updateProgressBar(100);
            updateQuickStatus('Analysis complete');
            showNotification('SNP analysis completed successfully!', 'success');

        } catch (error) {
            console.error('Analysis error:', error);
            showNotification(`Analysis failed: ${error.message}`, 'error');
            updateQuickStatus('Analysis failed');
        } finally {
            analyzeBtn.disabled = false;
            analyzeBtn.textContent = '🔍 Analyze SNPs';
        }
    }

    // Update progress bar
    function updateProgressBar(percent) {
        if (progressBar && progressPercent) {
            progressBar.style.width = `${percent}%`;
            progressPercent.textContent = `${percent}%`;
        }
    }

    // Update quick status
    function updateQuickStatus(status) {
        if (quickStatus) {
            quickStatus.textContent = status;
        }
    }

    // Update basic results
    function updateBasicResults(data) {
        if (quickTotalSNPs) quickTotalSNPs.textContent = data.total_snps || 0;
        if (quickSimilarity) quickSimilarity.textContent = `${(data.similarity || 0).toFixed(1)}%`;
        if (quickRefLength) quickRefLength.textContent = `${data.reference_length || 0} bp`;
        if (quickSampleLength) quickSampleLength.textContent = `${data.sample_length || 0} bp`;

        // Show SNP breakdown
        const snpBreakdown = document.getElementById('snp-breakdown');
        if (snpBreakdown && data.total_snps > 0) {
            snpBreakdown.style.display = 'block';
            
            const transitionsCount = document.getElementById('transitions-count');
            const transversionsCount = document.getElementById('transversions-count');
            const titvRatio = document.getElementById('titv-ratio');
            
            if (transitionsCount) transitionsCount.textContent = data.transitions || 0;
            if (transversionsCount) transversionsCount.textContent = data.transversions || 0;
            if (titvRatio) titvRatio.textContent = (data.ti_tv_ratio || 0).toFixed(2);
        }
    }

    // Perform additional analyses based on selected options
    async function performAdditionalAnalyses(referenceSeq, sampleSeq, basicData, options) {
        const detailedResults = document.getElementById('detailed-results');
        if (detailedResults) {
            detailedResults.style.display = 'block';
        }

        // SNP Highlighting
        if (options.highlightSNPs) {
            const snpHighlightResults = document.getElementById('snp-highlight-results');
            if (snpHighlightResults) {
                snpHighlightResults.style.display = 'block';
                const highlightedPositions = document.getElementById('highlighted-positions');
                if (highlightedPositions) {
                    highlightedPositions.textContent = `${basicData.total_snps} positions marked in sequences`;
                }
                
                // Highlight SNPs in the loaded sequences
                highlightSNPsInSequences(basicData.snps, referenceSeq, sampleSeq);
            }
        }

        // Stop Codon Detection
        if (options.stopCodonDetection) {
            try {
                const stopCodonResult = await callPythonScript('stop_codons', referenceSeq, sampleSeq);
                if (stopCodonResult.success) {
                    const stopCodonResults = document.getElementById('stop-codon-results');
                    if (stopCodonResults) {
                        stopCodonResults.style.display = 'block';
                        const refStopCodons = document.getElementById('ref-stop-codons');
                        const sampleStopCodons = document.getElementById('sample-stop-codons');
                        
                        if (refStopCodons) refStopCodons.textContent = `${stopCodonResult.data.reference_stop_codons || 0} stop codons`;
                        if (sampleStopCodons) sampleStopCodons.textContent = `${stopCodonResult.data.sample_stop_codons || 0} stop codons`;
                    }
                }
            } catch (error) {
                console.error('Stop codon analysis failed:', error);
            }
        }

        // SNP Frequency Analysis
        if (options.snpFrequencyCounter) {
            const snpFrequencyResults = document.getElementById('snp-frequency-results');
            if (snpFrequencyResults && basicData.snps) {
                snpFrequencyResults.style.display = 'block';
                
                // Calculate most common SNP type
                const snpTypes = basicData.snp_types || {};
                let mostCommon = 'None';
                let maxCount = 0;
                
                for (const [type, count] of Object.entries(snpTypes)) {
                    if (count > maxCount) {
                        maxCount = count;
                        mostCommon = type;
                    }
                }
                
                const mostCommonSNP = document.getElementById('most-common-snp');
                const snpDensity = document.getElementById('snp-density');
                
                if (mostCommonSNP) mostCommonSNP.textContent = maxCount > 0 ? `${mostCommon} (${maxCount})` : 'None';
                if (snpDensity) {
                    const density = basicData.reference_length > 0 ? ((basicData.total_snps / basicData.reference_length) * 1000).toFixed(2) : '0';
                    snpDensity.textContent = `${density} per 1000bp`;
                }
            }
        }

        // Synonymous/Non-synonymous Analysis
        if (options.synonymous || options.nonSynonymous || options.stopGained || options.startLost) {
            try {
                const codingResult = await callPythonScript('coding_impact', referenceSeq, sampleSeq);
                if (codingResult.success) {
                    const synonymousResults = document.getElementById('synonymous-results');
                    if (synonymousResults) {
                        synonymousResults.style.display = 'block';
                        
                        const synonymousCount = document.getElementById('synonymous-count');
                        const nonSynonymousCount = document.getElementById('non-synonymous-count');
                        const stopGainedCount = document.getElementById('stop-gained-count');
                        const startLostCount = document.getElementById('start-lost-count');
                        
                        if (synonymousCount) synonymousCount.textContent = codingResult.data.synonymous || 0;
                        if (nonSynonymousCount) nonSynonymousCount.textContent = codingResult.data.non_synonymous || 0;
                        if (stopGainedCount) stopGainedCount.textContent = codingResult.data.stop_gained || 0;
                        if (startLostCount) startLostCount.textContent = codingResult.data.start_lost || 0;
                    }
                }
            } catch (error) {
                console.error('Coding impact analysis failed:', error);
            }
        }
    }

    // Toggle result section visibility with smooth animation
    function toggleResultSection(element, show) {
        if (!element) return;
        
        if (show) {
            element.style.display = 'block';
            // Trigger reflow to ensure display change takes effect
            element.offsetHeight;
            element.style.opacity = '0';
            element.style.transform = 'translateY(-10px)';
            
            // Animate in
            setTimeout(() => {
                element.style.transition = 'all 0.3s ease-out';
                element.style.opacity = '1';
                element.style.transform = 'translateY(0)';
            }, 10);
        } else {
            element.style.transition = 'all 0.3s ease-out';
            element.style.opacity = '0';
            element.style.transform = 'translateY(-10px)';
            
            // Hide after animation
            setTimeout(() => {
                element.style.display = 'none';
            }, 300);
        }
    }

    // Perform additional analyses based on options selected when "Analyze" was clicked
    async function performAdditionalAnalysesUpdated(referenceSeq, sampleSeq, basicData, options) {
        const detailedResults = document.getElementById('detailed-results');
        
        // Check if any analysis options were selected when analyze was clicked
        const hasAnyAnalysis = options.highlightSNPs || options.stopCodonDetection || 
                              options.snpFrequencyCounter || options.codingImpactAnalysis;
        
        if (detailedResults) {
            detailedResults.style.display = hasAnyAnalysis ? 'block' : 'none';
        }

        // SNP Highlighting - only analyze and show if option was selected
        if (options.highlightSNPs && snpHighlightResults) {
            snpHighlightResults.style.display = 'block';
            const highlightedPositions = document.getElementById('highlighted-positions');
            if (highlightedPositions) {
                highlightedPositions.textContent = `${basicData.total_snps} positions marked in sequences`;
            }
            
            // Highlight SNPs in the output sequences
            highlightSNPsInOutputSequences(basicData.snps, referenceSeq, sampleSeq);
        } else if (snpHighlightResults) {
            snpHighlightResults.style.display = 'none';
        }

        // Stop Codon Detection - only analyze and show if option was selected
        if (options.stopCodonDetection && stopCodonResults) {
            try {
                const stopCodonResult = await callPythonScript('stop_codons', referenceSeq, sampleSeq);
                if (stopCodonResult.success) {
                    stopCodonResults.style.display = 'block';
                    const refStopCodons = document.getElementById('ref-stop-codons');
                    const sampleStopCodons = document.getElementById('sample-stop-codons');
                    
                    if (refStopCodons) refStopCodons.textContent = `${stopCodonResult.data.reference_stop_codons || 0} stop codons`;
                    if (sampleStopCodons) sampleStopCodons.textContent = `${stopCodonResult.data.sample_stop_codons || 0} stop codons`;
                }
            } catch (error) {
                console.error('Stop codon analysis failed:', error);
            }
        } else if (stopCodonResults) {
            stopCodonResults.style.display = 'none';
        }

        // SNP Frequency Analysis - only analyze and show if option was selected
        if (options.snpFrequencyCounter && snpFrequencyResults && basicData.snps) {
            snpFrequencyResults.style.display = 'block';
            
            // Calculate most common SNP type
            const snpTypes = basicData.snp_types || {};
            let mostCommon = 'None';
            let maxCount = 0;
            
            for (const [type, count] of Object.entries(snpTypes)) {
                if (count > maxCount) {
                    maxCount = count;
                    mostCommon = type;
                }
            }
            
            const mostCommonSNP = document.getElementById('most-common-snp');
            const snpDensity = document.getElementById('snp-density');
            
            if (mostCommonSNP) mostCommonSNP.textContent = maxCount > 0 ? `${mostCommon} (${maxCount})` : 'None';
            if (snpDensity) {
                const density = basicData.reference_length > 0 ? ((basicData.total_snps / basicData.reference_length) * 1000).toFixed(2) : '0';
                snpDensity.textContent = `${density} per 1000bp`;
            }
        } else if (snpFrequencyResults) {
            snpFrequencyResults.style.display = 'none';
        }

        // Coding Impact Analysis - only analyze and show if option was selected
        if (options.codingImpactAnalysis && synonymousResults) {
            try {
                const codingResult = await callPythonScript('coding_impact', referenceSeq, sampleSeq);
                if (codingResult.success) {
                    synonymousResults.style.display = 'block';
                    
                    const synonymousCount = document.getElementById('synonymous-count');
                    const nonSynonymousCount = document.getElementById('non-synonymous-count');
                    const stopGainedCount = document.getElementById('stop-gained-count');
                    const startLostCount = document.getElementById('start-lost-count');
                    
                    if (synonymousCount) synonymousCount.textContent = codingResult.data.synonymous || 0;
                    if (nonSynonymousCount) nonSynonymousCount.textContent = codingResult.data.non_synonymous || 0;
                    if (stopGainedCount) stopGainedCount.textContent = codingResult.data.stop_gained || 0;
                    if (startLostCount) startLostCount.textContent = codingResult.data.start_lost || 0;
                    
                    // Store results for later use
                    window.lastCodingResults = codingResult.data;
                    
                    // Display detailed mutations
                    displayDetailedMutations(codingResult.data);
                    
                    // Highlight coding mutations in the output sequences
                    highlightCodingMutationsInOutput(codingResult.data, referenceSeq, sampleSeq);
                }
            } catch (error) {
                console.error('Coding impact analysis failed:', error);
            }
        } else if (synonymousResults) {
            synonymousResults.style.display = 'none';
        }
    }

    // Display detailed mutations in the results
    function displayDetailedMutations(codingData) {
        const mutationsDetails = document.getElementById('mutations-details');
        if (!mutationsDetails || !codingData.mutations) return;
        
        mutationsDetails.innerHTML = '';
        
        // Group mutations by type for better display
        const mutationTypes = [
            { key: 'start_lost_mutations', label: 'Start Lost', color: 'text-yellow-400', bgColor: 'bg-yellow-500/10' },
            { key: 'stop_gained_mutations', label: 'Stop Gained', color: 'text-red-400', bgColor: 'bg-red-500/10' },
            { key: 'non_synonymous_mutations', label: 'Non-synonymous', color: 'text-orange-400', bgColor: 'bg-orange-500/10' },
            { key: 'synonymous_mutations', label: 'Synonymous', color: 'text-green-400', bgColor: 'bg-green-500/10' }
        ];
        
        mutationTypes.forEach(type => {
            const mutations = codingData[type.key] || [];
            if (mutations.length > 0) {
                mutations.forEach(mutation => {
                    const mutationDiv = document.createElement('div');
                    mutationDiv.className = `${type.bgColor} rounded p-2 text-xs`;
                    mutationDiv.innerHTML = `
                        <div class="flex items-center justify-between">
                            <span class="${type.color} font-medium">${type.label}</span>
                            <span class="text-gray-400">Pos ${mutation.position}</span>
                        </div>
                        <div class="mt-1 font-mono text-xs">
                            <span class="text-blue-300">${mutation.ref_codon}</span>
                            <span class="text-gray-500 mx-1">→</span>
                            <span class="text-green-300">${mutation.sample_codon}</span>
                            <span class="text-gray-500 ml-2">(${mutation.ref_aa} → ${mutation.sample_aa})</span>
                        </div>
                    `;
                    mutationsDetails.appendChild(mutationDiv);
                });
            }
        });
        
        if (mutationsDetails.children.length === 0) {
            mutationsDetails.innerHTML = '<div class="text-gray-500 text-xs text-center py-2">No mutations detected</div>';
        }
    }

    // Simple coding mutation highlighting - just add highlights, keep everything else the same
    function highlightCodingMutationsInOutput(codingData, referenceSeq, sampleSeq) {
        // Just show a simple notification for now
        if (codingData.mutations && codingData.mutations.length > 0) {
            console.log('Coding mutation highlighting enabled for', codingData.mutations.length, 'mutations');
        }
    }

    // Format sequence with line breaks every 60 characters
    function formatSequenceWithLineBreaks(sequence) {
        const lineLength = 60;
        let formatted = '';
        let i = 0;
        
        while (i < sequence.length) {
            // Find the next 60 characters, but be careful with HTML tags
            let line = '';
            let charCount = 0;
            
            while (charCount < lineLength && i < sequence.length) {
                if (sequence[i] === '<') {
                    // Skip HTML tags
                    let tagEnd = sequence.indexOf('>', i);
                    if (tagEnd !== -1) {
                        line += sequence.substring(i, tagEnd + 1);
                        i = tagEnd + 1;
                    } else {
                        line += sequence[i];
                        i++;
                        charCount++;
                    }
                } else {
                    line += sequence[i];
                    i++;
                    charCount++;
                }
            }
            
            formatted += line + '\n';
        }
        
        return formatted;
    }

    // Advanced SNP highlighting - highlight ONLY the mismatch positions on BOTH sides
    function highlightSNPsInOutputSequences(snps, referenceSeq, sampleSeq) {
        if (!snps || snps.length === 0) return;
        
        console.log('SNP highlighting enabled for', snps.length, 'positions');
        
        // Find both reference and sample content areas
        const sampleContent = document.getElementById('sample-results-content');
        const referenceContent = document.getElementById('reference-results-content');
        
        // Create a map of SNP positions for quick lookup (convert to 0-based)
        const snpMap = {};
        snps.forEach(snp => {
            snpMap[snp.position - 1] = snp; // Convert to 0-based indexing
        });
        
        // Highlight SAMPLE sequence (red/green for mutations)
        if (sampleContent) {
            sampleContent.classList.remove('snp-highlighting-active');
            
            let currentHTML = sampleContent.innerHTML;
            
            // If there's no sequence yet, just store the SNPs for later
            if (!currentHTML || currentHTML.includes('Sample sequence results will appear here')) {
                sampleContent.setAttribute('data-snp-positions', JSON.stringify(snps));
            } else {
                // FIRST: Remove any existing SNP highlights to prevent double highlighting
                const cleanHTML = removeExistingSNPHighlights(currentHTML);
                
                // THEN: Process the sequence and highlight only mismatch positions
                const highlightedHTML = highlightSpecificPositions(cleanHTML, snpMap, sampleSeq, 'sample');
                sampleContent.innerHTML = highlightedHTML;
            }
        }
        
        // Highlight REFERENCE sequence (blue for original bases)
        if (referenceContent) {
            referenceContent.classList.remove('snp-highlighting-active');
            
            let currentHTML = referenceContent.innerHTML;
            
            // If there's no sequence yet, just store the SNPs for later
            if (!currentHTML || currentHTML.includes('Reference sequence results will appear here')) {
                referenceContent.setAttribute('data-snp-positions', JSON.stringify(snps));
            } else {
                // FIRST: Remove any existing SNP highlights to prevent double highlighting
                const cleanHTML = removeExistingSNPHighlights(currentHTML);
                
                // THEN: Process the sequence and highlight only mismatch positions
                const highlightedHTML = highlightSpecificPositions(cleanHTML, snpMap, referenceSeq, 'reference');
                referenceContent.innerHTML = highlightedHTML;
            }
        }
        
        console.log('SNP highlighting applied to', snps.length, 'positions on both sequences');
    }
    
    // Helper function to remove existing SNP highlights before applying new ones
    function removeExistingSNPHighlights(html) {
        // Remove all existing SNP highlight spans but keep the nucleotide content
        let cleanHTML = html;
        
        // Remove sample SNP highlights: <span class="snp-highlight" title="...">X</span> → X
        cleanHTML = cleanHTML.replace(/<span class="snp-highlight"[^>]*>([ATGCN])<\/span>/gi, '$1');
        
        // Remove reference SNP highlights: <span class="snp-highlight-reference" title="...">X</span> → X
        cleanHTML = cleanHTML.replace(/<span class="snp-highlight-reference"[^>]*>([ATGCN])<\/span>/gi, '$1');
        
        console.log('Removed existing SNP highlights');
        return cleanHTML;
    }
    
    // Helper function to highlight specific nucleotide positions
    function highlightSpecificPositions(html, snpMap, sequence, sequenceType = 'sample') {
        if (!sequence) return html;
        
        console.log(`Processing ${sequenceType} HTML:`, html.substring(0, 200) + '...');
        
        // Clean sequence for position mapping
        const cleanSeq = sequence.replace(/[^ATGCN]/gi, '');
        
        // Different CSS classes for reference vs sample
        const cssClass = sequenceType === 'reference' ? 'snp-highlight-reference' : 'snp-highlight';
        
        // More robust HTML parsing - handle existing HTML tags properly
        let sequencePosition = 0;
        let result = '';
        let i = 0;
        
        while (i < html.length) {
            const char = html[i];
            
            // If we encounter an HTML tag, copy it entirely without modification
            if (char === '<') {
                let tagEnd = html.indexOf('>', i);
                if (tagEnd !== -1) {
                    result += html.substring(i, tagEnd + 1);
                    i = tagEnd + 1;
                    continue;
                }
            }
            
            // If it's a DNA nucleotide (not inside a tag)
            if (/[ATGCN]/i.test(char)) {
                // Check if this position has a SNP
                if (snpMap[sequencePosition]) {
                    const snp = snpMap[sequencePosition];
                    // Different tooltip text for reference vs sample
                    const tooltipText = sequenceType === 'reference' 
                        ? `Original: ${snp.reference} (will become ${snp.sample})`
                        : `SNP: ${snp.reference}→${snp.sample} (${snp.type})`;
                    
                    // Highlight this specific nucleotide
                    result += `<span class="${cssClass}" title="${tooltipText}">${char}</span>`;
                } else {
                    result += char;
                }
                sequencePosition++;
            } else {
                // Not a DNA nucleotide (space, number, newline, etc.)
                result += char;
            }
            
            i++;
        }
        
        console.log(`Processed ${sequenceType} HTML:`, result.substring(0, 200) + '...');
        return result;
    }

    // Main function called from analysis results
    function highlightSNPsInSequences(snps, referenceSeq, sampleSeq) {
        console.log('highlightSNPsInSequences called with', snps?.length || 0, 'SNPs');
        
        // Call our detailed highlighting function
        highlightSNPsInOutputSequences(snps, referenceSeq, sampleSeq);
    }



    // Simple combined highlighting
    function combinedHighlighting(snps, codingData, referenceSeq, sampleSeq) {
        console.log('Combined highlighting enabled');
    }

    // Store the last coding analysis results for highlighting
    window.lastCodingResults = null;

    // Copy and download button event listeners
    document.addEventListener('click', function(e) {
        if (e.target.tagName === 'BUTTON') {
            console.log('SNP Button clicked:', e.target.id, e.target.textContent);
            
            // Handle copy button - copy SAMPLE sequence (not reference)
            if (e.target.id === 'copy-snp-results') {
                console.log('Handling SNP copy button!');
                const sampleResultsEl = document.getElementById('sample-results-content');
                console.log('Sample results element:', sampleResultsEl);
                console.log('Sample results content:', sampleResultsEl?.textContent);
                
                if (sampleResultsEl && sampleResultsEl.textContent && sampleResultsEl.textContent.trim()) {
                    const rawText = sampleResultsEl.textContent;
                    console.log('Raw SNP sample text:', rawText);
                    
                    // Extract only the DNA sequence (remove line numbers, labels, etc.)
                    const cleanSequence = rawText
                        .replace(/^\d+\|\s*/gm, '') // Remove line numbers like "00001| "
                        .replace(/5'\s*-\s*3'/g, '') // Remove strand indicators
                        .replace(/3'\s*-\s*5'/g, '') // Remove reverse strand indicators
                        .replace(/\s+/g, '') // Remove all whitespace
                        .replace(/[^ATGCNRYSWKMBDHV]/gi, ''); // Keep only valid DNA bases
                    
                    console.log('Cleaned SNP sample sequence:', cleanSequence);
                    
                    if (cleanSequence) {
                        navigator.clipboard.writeText(cleanSequence).then(() => {
                            console.log('SNP copy successful!');
                            e.target.textContent = '✅ Copied!';
                            setTimeout(() => {
                                e.target.textContent = '📋 Copy Results';
                            }, 2000);
                        }).catch(err => {
                            console.error('SNP copy failed:', err);
                            showNotification('Failed to copy sample sequence', 'error');
                        });
                    } else {
                        console.log('No valid DNA sequence found in SNP sample');
                        showNotification('No valid DNA sequence found', 'error');
                    }
                } else {
                    console.log('No sample sequence to copy');
                    showNotification('No sample sequence to copy. Please run analysis first.', 'error');
                }
            }
            
            // Handle download button - download SAMPLE sequence (not reference)
            if (e.target.id === 'download-snp-results') {
                console.log('Handling SNP download button!');
                const sampleResultsEl = document.getElementById('sample-results-content');
                
                if (sampleResultsEl && sampleResultsEl.textContent && sampleResultsEl.textContent.trim()) {
                    const rawText = sampleResultsEl.textContent;
                    console.log('Raw SNP sample download text:', rawText);
                    
                    // Extract only the DNA sequence (remove line numbers, labels, etc.)
                    const cleanSequence = rawText
                        .replace(/^\d+\|\s*/gm, '') // Remove line numbers like "00001| "
                        .replace(/5'\s*-\s*3'/g, '') // Remove strand indicators
                        .replace(/3'\s*-\s*5'/g, '') // Remove reverse strand indicators
                        .replace(/\s+/g, '') // Remove all whitespace
                        .replace(/[^ATGCNRYSWKMBDHV]/gi, ''); // Keep only valid DNA bases
                    
                    console.log('Clean SNP sample sequence for download:', cleanSequence);
                    
                    if (cleanSequence) {
                        const blob = new Blob([`>Sample_Sequence_with_SNPs\n${cleanSequence}`], { type: 'text/plain' });
                        const url = URL.createObjectURL(blob);
                        const a = document.createElement('a');
                        a.href = url;
                        a.download = 'sample_sequence_snps.fasta';
                        document.body.appendChild(a);
                        a.click();
                        document.body.removeChild(a);
                        URL.revokeObjectURL(url);
                        console.log('SNP download initiated');
                    } else {
                        console.log('No valid DNA sequence found for SNP download');
                        showNotification('No valid DNA sequence found', 'error');
                    }
                } else {
                    console.log('No sample sequence to download');
                    showNotification('No sample sequence to download. Please run analysis first.', 'error');
                }
            }
        }
    });

    // Initialize the page state when everything is loaded
    initializePageState();
});