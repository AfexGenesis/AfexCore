// ORF Cleaner JavaScript
document.addEventListener('DOMContentLoaded', function() {
    console.log('ORF Cleaner script loaded!');
    
    // Elements
    const fileUpload = document.getElementById('orf-file-upload');
    const dnaInput = document.getElementById('orf-dna-input');
    const sequenceLength = document.getElementById('orf-sequence-length');
    const gcContent = document.getElementById('orf-gc-content');
    const clearBtn = document.getElementById('orf-clear-btn');
    const toggleAdvanced = document.getElementById('orf-toggle-advanced');
    const advancedOptions = document.getElementById('orf-advanced-options');
    const analyzeBtn = document.getElementById('orf-analyze-btn');
    const optimizeBtn = document.getElementById('orf-optimize-btn');
    
    // Output elements
    const cleanedOutput = document.getElementById('orf-cleaned-output');
    const statsOutput = document.getElementById('orf-stats-output');
    const emptyOutput = document.getElementById('orf-empty-output');
    const cleanedResult = document.getElementById('orf-cleaned-result');
    const improvementsCount = document.getElementById('orf-improvements-count');
    
    // Copy and download buttons
    const copyCleanedBtn = document.getElementById('copy-orf-cleaned');
    const downloadCleanedBtn = document.getElementById('download-orf-cleaned');
    
    console.log('Copy button:', copyCleanedBtn);
    console.log('Download button:', downloadCleanedBtn);
    console.log('Cleaned result element:', cleanedResult);
    
    // Check if buttons are actually found
    console.log('Copy button found:', !!copyCleanedBtn);
    console.log('Download button found:', !!downloadCleanedBtn);
    
    // Test if buttons are clickable at all
    if (copyCleanedBtn) {
        console.log('Copy button styles:', window.getComputedStyle(copyCleanedBtn));
        console.log('Copy button parent:', copyCleanedBtn.parentElement);
        console.log('Copy button visible:', copyCleanedBtn.offsetWidth > 0 && copyCleanedBtn.offsetHeight > 0);
    }
    
    // Add a simple test click handler that also handles copy/download
    document.addEventListener('click', function(e) {
        if (e.target.tagName === 'BUTTON') {
            console.log('Button clicked:', e.target.id, e.target.textContent);
            
            // Handle copy button
            if (e.target.id === 'copy-orf-cleaned') {
                console.log('Handling copy button in general handler!');
                const cleanedResultEl = document.getElementById('orf-cleaned-result');
                console.log('cleanedResult element:', cleanedResultEl);
                console.log('cleanedResult.textContent:', cleanedResultEl?.textContent);
                
                if (cleanedResultEl && cleanedResultEl.textContent && cleanedResultEl.textContent.trim()) {
                    const rawText = cleanedResultEl.textContent;
                    console.log('Raw text:', rawText);
                    
                    // Extract only the DNA sequence (remove line numbers, labels, etc.)
                    const cleanSequence = rawText
                        .replace(/^\d+\|\s*/gm, '') // Remove line numbers like "00001| "
                        .replace(/5'\s*-\s*3'/g, '') // Remove strand indicators
                        .replace(/3'\s*-\s*5'/g, '') // Remove reverse strand indicators
                        .replace(/\s+/g, '') // Remove all whitespace
                        .replace(/[^ATGCNRYSWKMBDHV]/gi, ''); // Keep only valid DNA bases
                    
                    console.log('Cleaned sequence:', cleanSequence);
                    
                    if (cleanSequence) {
                        navigator.clipboard.writeText(cleanSequence).then(() => {
                            console.log('Copy successful!');
                            e.target.textContent = '✅ Copied!';
                            setTimeout(() => {
                                e.target.textContent = '📋 Copy';
                            }, 2000);
                        }).catch(err => {
                            console.error('Copy failed:', err);
                            showNotification('Failed to copy sequence', 'error');
                        });
                    } else {
                        console.log('No valid DNA sequence found');
                        showNotification('No valid DNA sequence found', 'error');
                    }
                } else {
                    console.log('No text to copy');
                    showNotification('No sequence to copy', 'error');
                }
            }
            
            // Handle download button
            if (e.target.id === 'download-orf-cleaned') {
                console.log('Handling download button in general handler!');
                const cleanedResultEl = document.getElementById('orf-cleaned-result');
                
                if (cleanedResultEl && cleanedResultEl.textContent && cleanedResultEl.textContent.trim()) {
                    const rawText = cleanedResultEl.textContent;
                    console.log('Raw download text:', rawText);
                    
                    // Extract only the DNA sequence (remove line numbers, labels, etc.)
                    const cleanSequence = rawText
                        .replace(/^\d+\|\s*/gm, '') // Remove line numbers like "00001| "
                        .replace(/5'\s*-\s*3'/g, '') // Remove strand indicators
                        .replace(/3'\s*-\s*5'/g, '') // Remove reverse strand indicators
                        .replace(/\s+/g, '') // Remove all whitespace
                        .replace(/[^ATGCNRYSWKMBDHV]/gi, ''); // Keep only valid DNA bases
                    
                    console.log('Clean sequence for download:', cleanSequence);
                    
                    if (cleanSequence) {
                        const blob = new Blob([`>Optimized_ORF_Sequence\n${cleanSequence}`], { type: 'text/plain' });
                        const url = URL.createObjectURL(blob);
                        const a = document.createElement('a');
                        a.href = url;
                        a.download = 'optimized_orf_sequence.fasta';
                        document.body.appendChild(a);
                        a.click();
                        document.body.removeChild(a);
                        URL.revokeObjectURL(url);
                        console.log('Download initiated');
                    } else {
                        console.log('No valid DNA sequence found for download');
                        showNotification('No valid DNA sequence found', 'error');
                    }
                } else {
                    console.log('No sequence to download');
                    showNotification('No sequence to download', 'error');
                }
            }
        }
    });
    
    // Results indicators
    const outputGC = document.getElementById('orf-output-gc');
    const outputTm = document.getElementById('orf-output-tm');
    const loadBtn = document.getElementById('orf-load-btn');
    const editLockBtn = document.getElementById('orf-edit-lock-mode');

    // Statistics elements
    const statsOrfsFound = document.getElementById('stats-orfs-found');
    const statsOrfsOptimized = document.getElementById('stats-orfs-optimized');
    const statsCodonsImproved = document.getElementById('stats-codons-improved');
    const statsFinalOrfLength = document.getElementById('stats-final-orf-length');
    
    // Detailed issues elements
    const detailedIssuesSection = document.getElementById('orf-detailed-issues');
    const issuesList = document.getElementById('orf-issues-list');

    // File upload handler
    fileUpload.addEventListener('change', function(e) {
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
                
                dnaInput.value = content.toUpperCase();
                updateSequenceStats();
            };
            reader.readAsText(file);
        }
    });

    // DNA input handler
    dnaInput.addEventListener('input', function() {
        // Only allow valid DNA bases
        this.value = this.value.toUpperCase().replace(/[^ATGCNRYSWKMBDHV]/g, '');
        updateSequenceStats();
    });

    // Clear button
    clearBtn.addEventListener('click', function() {
        dnaInput.value = '';
        updateSequenceStats();
        hideResults();
    });

    // Toggle advanced options
    if (toggleAdvanced) {
        toggleAdvanced.addEventListener('click', function() {
            if (advancedOptions.style.display === 'none') {
                advancedOptions.style.display = 'block';
                this.textContent = '⚙️ Hide Advanced Options';
            } else {
                advancedOptions.style.display = 'none';
                this.textContent = '⚙️ Show Advanced Options';
            }
        });
    }

    // Analyze ORFs button (Preview Mode)
    analyzeBtn.addEventListener('click', function() {
        // Get sequence to analyze - prefer output sequence if available, otherwise use input
        let sequenceToAnalyze = '';
        
        // Check if there's a cleaned sequence in the output
        if (cleanedResult && cleanedResult.textContent && cleanedResult.textContent.trim()) {
            // Extract sequence from the formatted output (remove line numbers, strand indicators, and formatting)
            sequenceToAnalyze = cleanedResult.textContent
                .replace(/5'/g, '')  // Remove 5' indicators
                .replace(/3'/g, '')  // Remove 3' indicators
                .replace(/\d+:/g, '') // Remove line numbers with colons
                .replace(/\s+/g, '') // Remove all whitespace
                .replace(/[^ATGCN]/gi, ''); // Keep only valid DNA bases
            
            if (sequenceToAnalyze) {
                showNotification('Analyzing the current output sequence for ORF issues (preview mode)', 'info');
            }
        }
        
        // If no output sequence, use input sequence
        if (!sequenceToAnalyze) {
            const inputSequence = dnaInput.value.trim();
            if (!inputSequence) {
                showNotification('Please enter a DNA sequence', 'error');
                return;
            }
            sequenceToAnalyze = inputSequence;
        }
        
        previewOrfIssues(sequenceToAnalyze);
    });

    // Optimize ORFs button (Apply Fixes)
    optimizeBtn.addEventListener('click', function() {
        // Get sequence to optimize - prefer output sequence if available, otherwise use input
        let sequenceToOptimize = '';
        
        // Check if there's a cleaned sequence in the output
        if (cleanedResult && cleanedResult.textContent && cleanedResult.textContent.trim()) {
            // Extract sequence from the formatted output (remove line numbers, strand indicators, and formatting)
            sequenceToOptimize = cleanedResult.textContent
                .replace(/5'/g, '')  // Remove 5' indicators
                .replace(/3'/g, '')  // Remove 3' indicators
                .replace(/\d+:/g, '') // Remove line numbers with colons
                .replace(/\s+/g, '') // Remove all whitespace
                .replace(/[^ATGCN]/gi, ''); // Keep only valid DNA bases
            
            if (sequenceToOptimize) {
                showNotification('Optimizing the current output sequence (applying fixes)', 'info');
            }
        }
        
        // If no output sequence, use input sequence
        if (!sequenceToOptimize) {
            const inputSequence = dnaInput.value.trim();
            if (!inputSequence) {
                showNotification('Please enter a DNA sequence', 'error');
                return;
            }
            sequenceToOptimize = inputSequence;
        }
        
        optimizeOrfs(sequenceToOptimize);
    });

    // Load button functionality
    loadBtn.addEventListener('click', function() {
        const sequence = dnaInput.value.trim();
        
        if (!sequence) {
            showNotification('Please enter a DNA sequence to load', 'error');
            return;
        }
        
        loadSequence(sequence);
    });

    // Direct event listeners for copy and download buttons
    if (copyCleanedBtn) {
        copyCleanedBtn.addEventListener('click', function(e) {
            console.log('Copy button clicked directly!');
            const cleanedResultEl = document.getElementById('orf-cleaned-result');
            console.log('cleanedResult element:', cleanedResultEl);
            console.log('cleanedResult.textContent:', cleanedResultEl?.textContent);
            
            if (cleanedResultEl && cleanedResultEl.textContent && cleanedResultEl.textContent.trim()) {
                const text = cleanedResultEl.textContent;
                console.log('Copying text:', text);
                navigator.clipboard.writeText(text).then(() => {
                    console.log('Copy successful!');
                    this.textContent = '✅ Copied!';
                    setTimeout(() => {
                        this.textContent = '📋 Copy';
                    }, 2000);
                }).catch(err => {
                    console.error('Copy failed:', err);
                    showNotification('Failed to copy sequence', 'error');
                });
            } else {
                console.log('No text to copy');
                showNotification('No sequence to copy', 'error');
            }
        });
    }

    if (downloadCleanedBtn) {
        downloadCleanedBtn.addEventListener('click', function(e) {
            console.log('Download button clicked directly!');
            const cleanedResultEl = document.getElementById('orf-cleaned-result');
            
            if (cleanedResultEl && cleanedResultEl.textContent && cleanedResultEl.textContent.trim()) {
                const sequence = cleanedResultEl.textContent;
                console.log('Downloading sequence:', sequence);
                const blob = new Blob([`>Optimized_ORF_Sequence\n${sequence}`], { type: 'text/plain' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = 'optimized_orf_sequence.fasta';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(url);
                console.log('Download initiated');
            } else {
                console.log('No sequence to download');
                showNotification('No sequence to download', 'error');
            }
        });
    }

    // Update sequence statistics
    function updateSequenceStats() {
        const sequence = dnaInput.value.replace(/\s/g, '');
        const length = sequence.length;
        
        sequenceLength.textContent = `Length: ${length} bp`;
        
        if (length > 0) {
            const gcCount = (sequence.match(/[GC]/g) || []).length;
            const gcPercent = ((gcCount / length) * 100).toFixed(1);
            gcContent.textContent = `GC: ${gcPercent}%`;
        } else {
            gcContent.textContent = 'GC: 0%';
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
            } else if (line.startsWith('//')) {
                break;
            } else if (inSequence) {
                // Remove line numbers and spaces
                const seqLine = line.replace(/^\s*\d+/, '').replace(/\s+/g, '');
                sequence += seqLine;
            }
        }
        
        return sequence;
    }

    // Load sequence function
    async function loadSequence(sequence) {
        try {
            loadBtn.disabled = true;
            loadBtn.textContent = '🔄 Loading...';
            
            const result = await callPythonScript('load', sequence);
            
            if (result.success) {
                displayLoadedSequence(result.data);
                showNotification('Sequence loaded successfully!', 'success');
            } else {
                throw new Error(result.error || 'Failed to load sequence');
            }
            
        } catch (error) {
            console.error('Error:', error);
            showNotification('Error loading sequence: ' + error.message, 'error');
        } finally {
            loadBtn.disabled = false;
            loadBtn.textContent = '📁 LOAD';
        }
    }

    // Display loaded sequence in results area
    function displayLoadedSequence(data) {
        // Hide empty state
        emptyOutput.style.display = 'none';
        
        // Update GC and Tm indicators
        if (data.gc_content !== undefined) {
            outputGC.textContent = `${data.gc_content}%`;
        }
        if (data.melting_temp !== undefined) {
            outputTm.textContent = `${data.melting_temp}°C`;
        }
        
        // Show the loaded sequence in the cleaned output area
        if (data.formatted_sequence) {
            // Format sequence with proper line breaks and styling (like SNP Highlighter)
            const formattedHTML = formatSequenceWithColors(data.formatted_sequence);
            cleanedResult.innerHTML = formattedHTML;
            cleanedOutput.style.display = 'block';
            
            // Update the header to show it's a loaded sequence
            const cleanedHeader = document.querySelector('#orf-cleaned-output h4');
            if (cleanedHeader) {
                cleanedHeader.textContent = 'Loaded DNA Sequence';
                cleanedHeader.className = 'text-green-400 font-medium';
            }
            
            // Update improvements count to show sequence info
            improvementsCount.textContent = `${data.length} bp loaded`;
            improvementsCount.className = 'orf-badge orf-improvements-badge';
        }
    }

    // Preview ORF issues function (analyze without fixing)
    async function previewOrfIssues(sequence) {
        try {
            analyzeBtn.textContent = '🔍 Analyzing...';
            analyzeBtn.disabled = true;
            
            // Get analysis options
            const options = {
                min_orf_length: 30, // Fixed default value
                optimize_codons: document.getElementById('orf-optimize-codons')?.checked || true,
                target_organism: document.getElementById('orf-organism-selector')?.value || 'e-coli',
                fix_internal_starts: document.getElementById('orf-fix-internal-starts')?.checked || true,
                remove_internal_stops: document.getElementById('orf-fix-premature-stops')?.checked || true
            };
            
            const result = await callPythonScript('analyze', sequence, JSON.stringify(options));
            
            if (result.success) {
                displayPreviewResults(result);
                const totalIssues = result.total_issues || 0;
                if (totalIssues > 0) {
                    showNotification(`Found ${totalIssues} issues in ${result.orfs_with_issues} ORFs. Click OPTIMIZE to apply fixes.`, 'warning');
                } else {
                    showNotification('No ORF issues found! Your sequence is already optimized.', 'success');
                }
            } else {
                showNotification(result.error || 'Failed to analyze ORFs', 'error');
                return;
            }
            
        } catch (error) {
            console.error('Preview error:', error);
            showNotification(`Preview failed: ${error.message}`, 'error');
        } finally {
            analyzeBtn.textContent = '🔍 ANALYZE (Preview Issues)';
            analyzeBtn.disabled = false;
        }
    }

    // Optimize ORFs function (actually apply fixes)
    async function optimizeOrfs(sequence) {
        try {
            optimizeBtn.textContent = '🔄 Optimizing...';
            optimizeBtn.disabled = true;
            
            // Get optimization options
            const options = {
                min_orf_length: 30, // Fixed default value
                optimize_codons: document.getElementById('orf-optimize-codons')?.checked || true,
                target_organism: document.getElementById('orf-organism-selector')?.value || 'e-coli',
                fix_internal_starts: document.getElementById('orf-fix-internal-starts')?.checked || true,
                remove_internal_stops: document.getElementById('orf-fix-premature-stops')?.checked || true
            };
            
            const result = await callPythonScript('optimize', sequence, JSON.stringify(options));
            
            if (result.success) {
                displayOrfResults(result);
                const fixes = result.statistics.codons_fixed || 0;
                const improvements = result.statistics.codons_improved || 0;
                let message = `Successfully optimized ${result.statistics.orfs_optimized} ORFs`;
                if (fixes > 0 && improvements > 0) {
                    message += ` with ${fixes} problematic codons fixed and ${improvements} codon improvements!`;
                } else if (fixes > 0) {
                    message += ` with ${fixes} problematic codons fixed!`;
                } else if (improvements > 0) {
                    message += ` with ${improvements} codon improvements!`;
                } else {
                    message += `!`;
                }
                showNotification(message, 'success');
            } else {
                showNotification(result.error || 'Failed to optimize ORFs', 'error');
                return;
            }
            
        } catch (error) {
            console.error('Optimization error:', error);
            showNotification(`Optimization failed: ${error.message}`, 'error');
        } finally {
            optimizeBtn.textContent = '🧬 OPTIMIZE (Apply Fixes)';
            optimizeBtn.disabled = false;
        }
    }

    // Display preview results (issues found but not fixed)
    function displayPreviewResults(data) {
        // Hide empty state
        emptyOutput.style.display = 'none';
        
        // Update GC and Tm indicators
        if (data.sequence_analysis) {
            outputGC.textContent = `${data.sequence_analysis.gc_content}%`;
            outputTm.textContent = `${data.sequence_analysis.melting_temp}°C`;
        }
        
        // Show original sequence
        if (data.formatted_sequence) {
            const formattedHTML = formatSequenceWithColors(data.formatted_sequence);
            cleanedResult.innerHTML = formattedHTML;
            
            cleanedOutput.style.display = 'block';
            
            // Update header to show it's a preview
            const cleanedHeader = document.querySelector('#orf-cleaned-output h4');
            if (cleanedHeader) {
                cleanedHeader.textContent = 'ORF Analysis Preview';
                cleanedHeader.className = 'text-orange-400 font-medium';
            }
            
            // Update improvements count to show issues found
            const totalIssues = data.total_issues || 0;
            improvementsCount.textContent = `${totalIssues} issues found`;
            improvementsCount.className = totalIssues > 0 ? 
                'px-2 py-1 bg-orange-600/20 border border-orange-500/30 rounded text-orange-300 text-xs' : 
                'px-2 py-1 bg-green-600/20 border border-green-500/30 rounded text-green-300 text-xs';
        }
        
        // Show statistics
        if (data.statistics) {
            statsOrfsFound.textContent = data.statistics.orfs_found || 0;
            statsOrfsOptimized.textContent = data.statistics.orfs_with_issues || 0;
            statsCodonsImproved.textContent = data.statistics.issues_by_type?.suboptimal_codons || 0;
            statsFinalOrfLength.textContent = data.sequence_analysis?.length || 0;
            
            statsOutput.style.display = 'block';
        }
        
        // Show detailed issues
        console.log('Preview data received:', data);
        console.log('Statistics:', data.statistics);
        if (data.statistics && data.statistics.detailed_issues && data.statistics.detailed_issues.length > 0) {
            console.log('Detailed issues found:', data.statistics.detailed_issues);
            displayDetailedIssues(data.statistics.detailed_issues);
            detailedIssuesSection.classList.remove('hidden');
            detailedIssuesSection.style.display = 'block';
        } else {
            console.log('No detailed issues found or missing data');
            console.log('detailed_issues:', data.statistics?.detailed_issues);
            detailedIssuesSection.classList.add('hidden');
            detailedIssuesSection.style.display = 'none';
        }
    }

    // Display detailed issues
    function displayDetailedIssues(detailedIssues) {
        console.log('displayDetailedIssues called with:', detailedIssues);
        if (!issuesList) {
            console.log('issuesList element not found!');
            return;
        }
        
        issuesList.innerHTML = '';
        
        detailedIssues.forEach(issueType => {
            if (issueType.count > 0) {
                // Create issue type header
                const issueHeader = document.createElement('div');
                issueHeader.className = 'border-l-4 border-red-500 pl-3 mb-2';
                
                let issueIcon = '';
                let issueColor = '';
                let issueTitle = '';
                
                switch(issueType.type) {
                    case 'premature_stop_codon':
                        issueIcon = '⚠️';
                        issueColor = 'text-red-400';
                        issueTitle = 'Premature Stop Codons';
                        break;
                    case 'internal_stop_codon':
                        issueIcon = '🛑';
                        issueColor = 'text-orange-400';
                        issueTitle = 'Internal Stop Codons';
                        break;
                    case 'internal_start_codon':
                        issueIcon = '🔄';
                        issueColor = 'text-yellow-400';
                        issueTitle = 'Internal Start Codons';
                        break;
                    case 'suboptimal_codon':
                        issueIcon = '📊';
                        issueColor = 'text-blue-400';
                        issueTitle = 'Suboptimal Codons';
                        break;
                    default:
                        issueIcon = '❓';
                        issueColor = 'text-gray-400';
                        issueTitle = 'Other Issues';
                }
                
                issueHeader.innerHTML = `
                    <div class="flex items-center justify-between">
                        <div class="flex items-center gap-2">
                            <span>${issueIcon}</span>
                            <span class="${issueColor} font-medium">${issueTitle}</span>
                        </div>
                        <span class="px-2 py-1 bg-red-600/20 border border-red-500/30 rounded text-red-300 text-xs">
                            ${issueType.count} found
                        </span>
                    </div>
                `;
                
                issuesList.appendChild(issueHeader);
                console.log(`Added issue header for ${issueType.type} with ${issueType.count} issues`);
                
                // Create locations list
                if (issueType.locations && issueType.locations.length > 0) {
                    console.log(`Processing ${issueType.locations.length} locations for ${issueType.type}`);
                    const locationsList = document.createElement('div');
                    locationsList.className = 'ml-6 mt-2 space-y-1';
                    
                    issueType.locations.forEach((location, index) => {
                        if (index < 5) { // Show max 5 locations per issue type
                            const locationItem = document.createElement('div');
                            locationItem.className = 'text-gray-300 text-xs flex items-center gap-2';
                            
                            let locationText = '';
                            if (location.orf_position) {
                                locationText = `ORF ${location.orf_id}: Position ${location.orf_position} bp`;
                            } else {
                                locationText = `Global position ${location.position}`;
                            }
                            
                            if (location.codon) {
                                locationText += ` (${location.codon})`;
                            }
                            
                            locationItem.innerHTML = `
                                <span class="w-2 h-2 bg-red-400 rounded-full flex-shrink-0"></span>
                                <span>${locationText}</span>
                            `;
                            
                            locationsList.appendChild(locationItem);
                        }
                    });
                    
                    // Show "and X more" if there are more locations
                    if (issueType.locations.length > 5) {
                        const moreItem = document.createElement('div');
                        moreItem.className = 'text-gray-400 text-xs ml-4';
                        moreItem.textContent = `... and ${issueType.locations.length - 5} more`;
                        locationsList.appendChild(moreItem);
                    }
                    
                    issuesList.appendChild(locationsList);
                }
            }
        });
    }

    // Display ORF analysis results
    function displayOrfResults(data) {
        // Hide empty state
        emptyOutput.style.display = 'none';
        
        // Update GC and Tm indicators
        if (data.sequence_analysis) {
            outputGC.textContent = `${data.sequence_analysis.gc_content}%`;
            outputTm.textContent = `${data.sequence_analysis.melting_temp}°C`;
        }
        
        // Show cleaned sequence
        if (data.cleaned_sequence) {
            const formattedSeq = formatSequenceForDisplay(data.cleaned_sequence);
            const formattedHTML = formatSequenceWithColors(formattedSeq);
            cleanedResult.innerHTML = formattedHTML;
            cleanedOutput.style.display = 'block';
            
            // Update header
            const cleanedHeader = document.querySelector('#orf-cleaned-output h4');
            if (cleanedHeader) {
                cleanedHeader.textContent = 'Optimized ORF Sequence';
                cleanedHeader.className = 'text-yellow-400 font-medium';
            }
            
            // Update improvements count
            const totalImprovements = data.statistics.total_improvements || 0;
            improvementsCount.textContent = `${totalImprovements} improvements made`;
            improvementsCount.className = totalImprovements > 0 ? 
                'orf-badge orf-improvements-badge' : 
                'orf-badge orf-no-improvements-badge';
        }
        
        // Show statistics
        if (data.statistics) {
            statsOrfsFound.textContent = data.statistics.orfs_found || 0;
            statsOrfsOptimized.textContent = data.statistics.orfs_optimized || 0;
            statsCodonsImproved.textContent = data.statistics.codons_improved || 0;
            statsFinalOrfLength.textContent = data.sequence_analysis?.length || 0;
            
            statsOutput.style.display = 'block';
        }
    }

    // Format sequence for display
    function formatSequenceForDisplay(sequence, lineLength = 130) {
        const lines = [];
        for (let i = 0; i < sequence.length; i += lineLength) {
            const lineNum = String(i + 1).padStart(5, '0');
            const lineSeq = sequence.substring(i, i + lineLength);
            lines.push(`${lineNum}| ${lineSeq}`);
        }
        return lines.join('\n');
    }

    // Format sequence with colors and strand indicators (matching SNP Highlighter)
    function formatSequenceWithColors(formattedSequence) {
        const lines = formattedSequence.split('\n');
        let html = '';
        
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
                
                // Format with line numbers and strand indicators (exactly like SNP Highlighter)
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



    // Hide results
    function hideResults() {
        cleanedOutput.style.display = 'none';
        statsOutput.style.display = 'none';
        emptyOutput.style.display = 'block';
        outputGC.textContent = '0%';
        outputTm.textContent = '--°C';
        
        // Also hide detailed issues section
        if (detailedIssuesSection) {
            detailedIssuesSection.style.display = 'none';
        }
    }


    // Call Python script
    async function callPythonScript(operation, sequence, options = '') {
        try {
            const { spawn } = require('child_process');
            const path = require('path');
            
            const scriptPath = path.join(__dirname, 'assets', 'orf-cleaner.py');
            const args = [operation, sequence];
            
            if (options) {
                args.push(options);
            }
            
            return new Promise((resolve, reject) => {
                const python = spawn('python', [scriptPath, ...args]);
                let output = '';
                let error = '';
                
                python.stdout.on('data', (data) => {
                    output += data.toString();
                });
                
                python.stderr.on('data', (data) => {
                    error += data.toString();
                });
                
                python.on('close', (code) => {
                    if (code === 0) {
                        try {
                            const result = JSON.parse(output);
                            resolve(result);
                        } catch (e) {
                            reject(new Error('Failed to parse Python output: ' + output));
                        }
                    } else {
                        reject(new Error('Python script failed: ' + error));
                    }
                });
            });
            
        } catch (error) {
            throw new Error('Failed to execute Python script: ' + error.message);
        }
    }

    // Show notification
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

    // Initialize
    updateSequenceStats();
    hideResults();
});