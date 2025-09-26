// Restriction Site Cleaner JavaScript
document.addEventListener('DOMContentLoaded', function() {
    // Elements
    const fileUpload = document.getElementById('restriction-file-upload');
    const dnaInput = document.getElementById('restriction-dna-input');
    const sequenceLength = document.getElementById('restriction-sequence-length');
    const gcContent = document.getElementById('restriction-gc-content');
    const clearBtn = document.getElementById('restriction-clear-btn');
    const enzymeSelector = document.getElementById('enzyme-selector');
    const toggleAdvanced = document.getElementById('toggle-advanced');
    const advancedOptions = document.getElementById('advanced-options');
    const analyzeBtn = document.getElementById('analyze-restriction-btn');
    const analyzeOnlyBtn = document.getElementById('analyze-only-btn');
    
    // Output elements
    const cleanedOutput = document.getElementById('cleaned-output');
    const statsOutput = document.getElementById('stats-output');
    const reportOutput = document.getElementById('report-output');
    const emptyOutput = document.getElementById('restriction-empty-output');
    const cleanedResult = document.getElementById('cleaned-result');
    // statsContent and reportContent removed - now using individual elements
    const sitesRemovedCount = document.getElementById('sites-removed-count');
    
    // Copy and download buttons
    const copyCleanedBtn = document.getElementById('copy-cleaned');
    const downloadCleanedBtn = document.getElementById('download-cleaned');
    
    // Results indicators
    const outputGC = document.getElementById('restriction-output-gc');
    const outputTm = document.getElementById('restriction-output-tm');
    const loadBtn = document.getElementById('restriction-load-btn');
    const editLockBtn = document.getElementById('restriction-edit-lock-mode');

    // Track the last enzyme selection used for cleaning
    let lastEnzymeSelection = null;
    
    // Get additional elements for codon optimization
    const optimizeCodonsCheckbox = document.getElementById('optimize-codons');
    const organismSelector = document.getElementById('organism-selector');


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
    toggleAdvanced.addEventListener('click', function() {
        if (advancedOptions.style.display === 'none') {
            advancedOptions.style.display = 'block';
            this.textContent = '⚙️ Hide Advanced Options';
        } else {
            advancedOptions.style.display = 'none';
            this.textContent = '⚙️ Show Advanced Options';
        }
    });

    // Handle "Optimize codon usage" checkbox - force show advanced options
    if (optimizeCodonsCheckbox) {
        optimizeCodonsCheckbox.addEventListener('change', function() {
            if (this.checked) {
                // Force show advanced options when codon optimization is enabled
                advancedOptions.style.display = 'block';
                toggleAdvanced.textContent = '⚙️ Hide Advanced Options';
                
                // Show notification that organism selection is required
                showNotification('⚠️ Codon optimization enabled - please select target organism below', 'warning');
                
                // Highlight the organism selector
                if (organismSelector) {
                    organismSelector.style.borderColor = '#f59e0b';
                    organismSelector.style.boxShadow = '0 0 0 2px rgba(245, 158, 11, 0.2)';
                    
                    // Remove highlight after organism is selected
                    organismSelector.addEventListener('change', function() {
                        if (this.value) {
                            this.style.borderColor = '';
                            this.style.boxShadow = '';
                        }
                    }, { once: false });
                }
            }
        });
    }

    // Clean restriction sites button (works on output sequence if available)
    analyzeBtn.addEventListener('click', function() {
        // Get sequence to clean - prefer output sequence if available, otherwise use input
        let sequenceToClean = '';
        
        // Check if there's a cleaned sequence in the output
        if (cleanedResult && cleanedResult.textContent && cleanedResult.textContent.trim()) {
            // Extract sequence from the formatted output (remove line numbers and formatting)
            sequenceToClean = cleanedResult.textContent
                .replace(/^\d+\|/gm, '') // Remove line numbers like "00001|"
                .replace(/\s+/g, '') // Remove all whitespace
                .replace(/[^ATGCN]/gi, ''); // Keep only valid DNA bases
            
            if (sequenceToClean) {
                showNotification('Cleaning the current output sequence for additional processing', 'info');
            }
        }
        
        // If no output sequence, use input sequence
        if (!sequenceToClean) {
            const inputSequence = dnaInput.value.trim();
            if (!inputSequence) {
                showNotification('Please enter a DNA sequence', 'error');
                return;
            }
            sequenceToClean = inputSequence;
        }
        
        const selectedEnzyme = enzymeSelector ? String(enzymeSelector.value).trim() : null;
        
        if (!selectedEnzyme || selectedEnzyme === 'null' || selectedEnzyme === 'undefined' || selectedEnzyme.trim() === '') {
            showNotification('Please select restriction enzymes to remove', 'error');
            return;
        }
        
        // Check if codon optimization is enabled but no organism selected
        if (optimizeCodonsCheckbox && optimizeCodonsCheckbox.checked) {
            const selectedOrganism = organismSelector ? organismSelector.value.trim() : '';
            if (!selectedOrganism) {
                showNotification('⚠️ Codon optimization enabled - please select target organism in Advanced Options', 'error');
                
                // Force show advanced options and highlight organism selector
                advancedOptions.style.display = 'block';
                toggleAdvanced.textContent = '⚙️ Hide Advanced Options';
                
                if (organismSelector) {
                    organismSelector.style.borderColor = '#ef4444';
                    organismSelector.style.boxShadow = '0 0 0 2px rgba(239, 68, 68, 0.2)';
                    organismSelector.focus();
                }
                return;
            }
        }
        
        // Store the enzyme selection for later analysis
        lastEnzymeSelection = selectedEnzyme;
        analyzeRestrictionSites(sequenceToClean, selectedEnzyme);
    });

    // Analyze Only button (analyzes the cleaned/output sequence)
    analyzeOnlyBtn.addEventListener('click', function() {
        // Get the cleaned sequence from the output, or fall back to input if no cleaning done yet
        let sequenceToAnalyze = '';
        
        // Check if there's a cleaned sequence in the output
        if (cleanedResult && cleanedResult.textContent && cleanedResult.textContent.trim()) {
            // Extract sequence from the formatted output (remove line numbers and formatting)
            sequenceToAnalyze = cleanedResult.textContent
                .replace(/^\d+\|/gm, '') // Remove line numbers like "00001|"
                .replace(/\s+/g, '') // Remove all whitespace
                .replace(/[^ATGCN]/gi, ''); // Keep only valid DNA bases
        }
        
        // If no cleaned sequence, check if there's input sequence
        if (!sequenceToAnalyze) {
            const inputSequence = dnaInput.value.trim();
            if (!inputSequence) {
                showNotification('Please clean a sequence first, or enter a DNA sequence to analyze', 'error');
                return;
            }
            sequenceToAnalyze = inputSequence;
            // If analyzing input sequence, get current enzyme selection
            const currentEnzymeSelection = enzymeSelector ? String(enzymeSelector.value).trim() : null;
            if (currentEnzymeSelection && currentEnzymeSelection !== 'null' && currentEnzymeSelection !== 'undefined' && currentEnzymeSelection.trim() !== '') {
                lastEnzymeSelection = currentEnzymeSelection;
            }
            showNotification('No cleaned sequence found - analyzing input sequence with current enzyme selection', 'info');
        }
        
        analyzeOutputSequence(sequenceToAnalyze);
    });

    // Analyze restriction sites function
    async function analyzeRestrictionSites(sequence, enzymeSelection) {
        try {
            analyzeBtn.textContent = '🔄 Analyzing...';
            analyzeBtn.disabled = true;
            
            // Ensure enzymeSelection is a string
            enzymeSelection = String(enzymeSelection).trim();
            

            
            // Get cleaning method options
            const options = {
                preserve_orfs: document.getElementById('preserve-orfs')?.checked || false,
                silent_mutations: document.getElementById('silent-mutations')?.checked || false,
                optimize_codons: document.getElementById('optimize-codons')?.checked || false,
                preserve_regulatory: document.getElementById('preserve-regulatory')?.checked || false,
                target_organism: document.getElementById('organism-selector')?.value || '',
                gc_target: document.getElementById('gc-optimizer')?.value || ''
            };
            

            
            const result = await callPythonScript('analyze', sequence, enzymeSelection, JSON.stringify(options));
            

            
            if (result.success) {
                displayCleaningResults(result);
                showNotification(`Successfully cleaned ${result.statistics.sites_removed} restriction sites!`, 'success');
            } else {
                // Show error using HTML notification, not Electron popup
                showNotification(result.error || 'Failed to analyze restriction sites', 'error');
                return; // Don't throw error to avoid Electron popup
            }
            
        } catch (error) {
            console.error('Analysis error:', error);
            showNotification(`Analysis failed: ${error.message}`, 'error');
        } finally {
            analyzeBtn.textContent = '🧬 ANALYZE';
            analyzeBtn.disabled = false;
        }
    }

    // Analyze Output sequence function (analyzes the cleaned/output sequence)
    async function analyzeOutputSequence(sequence) {
        try {
            analyzeOnlyBtn.textContent = '🔄 Analyzing...';
            analyzeOnlyBtn.disabled = true;
            
            // Call Python script with 'scan' operation to analyze the output sequence
            // Use the same enzyme selection that was used for cleaning
            const result = await callPythonScript('scan', sequence, lastEnzymeSelection);
            
            if (result.success) {
                displayOutputAnalysisResults(result);
                if (result.total_sites === 0) {
                    showNotification('🎉 Perfect! No restriction sites found in the cleaned sequence!', 'success');
                } else {
                    showNotification(`⚠️ Found ${result.total_sites} remaining restriction sites across ${result.enzymes_found} enzymes`, 'warning');
                }
            } else {
                showNotification(result.error || 'Failed to analyze output sequence', 'error');
                return;
            }
            
        } catch (error) {
            console.error('Analysis error:', error);
            showNotification(`Analysis failed: ${error.message}`, 'error');
        } finally {
            analyzeOnlyBtn.textContent = '🔍 Analyze Output';
            analyzeOnlyBtn.disabled = false;
        }
    }

    // Display output analysis results (analyzes cleaned/output sequence)
    function displayOutputAnalysisResults(data) {
        // Hide empty state
        emptyOutput.style.display = 'none';
        
        // Show original sequence (no cleaning)
        if (data.sequence) {
            const formattedSeq = formatSequenceForDisplay(data.sequence);
            cleanedResult.innerHTML = formatSequenceWithColors(formattedSeq);
            cleanedOutput.style.display = 'block';
            
            // Update header to show it's output analysis
            const cleanedHeader = document.querySelector('#cleaned-output h4');
            if (cleanedHeader) {
                cleanedHeader.textContent = '🔍 Output Sequence Analysis';
                cleanedHeader.className = 'text-blue-400 font-medium';
            }
        }
        
        // Update GC and Tm indicators for original sequence
        if (data.sequence_analysis) {
            updateResultsIndicators(data.sequence_analysis.gc_content, data.sequence_analysis.melting_temp);
            if (data.total_sites === 0) {
                sitesRemovedCount.textContent = `✅ No sites remaining`;
                sitesRemovedCount.className = 'px-2 py-1 bg-green-600/20 border border-green-500/30 rounded text-green-300 text-xs';
            } else {
                sitesRemovedCount.textContent = `⚠️ ${data.total_sites} sites remaining`;
                sitesRemovedCount.className = 'px-2 py-1 bg-yellow-600/20 border border-yellow-500/30 rounded text-yellow-300 text-xs';
            }
        }
        
        // Show comprehensive statistics
        if (data.sequence_analysis) {
            displayAnalysisStats(data);
            statsOutput.style.display = 'block';
        }
        
        // Always show comprehensive site report (even if no sites found)
        displayComprehensiveSiteReport(data);
        reportOutput.style.display = 'block';
    }

    // Display analysis statistics (for analysis-only mode)
    function displayAnalysisStats(data) {
        const analysis = data.sequence_analysis;
        
        // Update statistics cards with color coding
        const sitesFoundElement = document.getElementById('stats-sites-found');
        sitesFoundElement.textContent = data.total_sites;
        sitesFoundElement.className = data.total_sites === 0 ? 'text-2xl font-bold text-green-400' : 'text-2xl font-bold text-yellow-400';
        
        document.getElementById('stats-sites-removed').textContent = '0';
        document.getElementById('stats-enzymes-processed').textContent = data.enzymes_found;
        document.getElementById('stats-final-length').textContent = analysis.length;
        
        // Update sequence properties
        document.getElementById('stats-gc-content').textContent = `${analysis.gc_content}%`;
        document.getElementById('stats-melting-temp').textContent = `${analysis.melting_temp}°C`;
        document.getElementById('stats-molecular-weight').textContent = `${analysis.molecular_weight.toLocaleString()} Da`;
        
        // Update base composition
        document.getElementById('stats-base-a').textContent = analysis.base_composition.A;
        document.getElementById('stats-base-t').textContent = analysis.base_composition.T;
        document.getElementById('stats-base-g').textContent = analysis.base_composition.G;
        document.getElementById('stats-base-c').textContent = analysis.base_composition.C;
        
        // Update stats header for output analysis mode
        const statsHeader = document.querySelector('#stats-output h4');
        if (statsHeader) {
            statsHeader.textContent = 'Output Sequence Analysis';
            statsHeader.className = 'text-blue-400 font-medium';
        }
    }

    // Display comprehensive site report (for output analysis mode)
    function displayComprehensiveSiteReport(data) {
        // Update enzymes found
        const enzymesContainer = document.getElementById('report-enzymes-targeted');
        const enzymesList = Object.keys(data.sites_found || {});
        
        if (enzymesList.length === 0) {
            enzymesContainer.innerHTML = '<span class="inline-block bg-green-500/20 text-green-300 px-2 py-1 rounded text-xs">🎉 No restriction sites found!</span>';
        } else {
            enzymesContainer.innerHTML = enzymesList.map(enzyme => 
                `<span class="inline-block bg-yellow-500/20 text-yellow-300 px-2 py-1 rounded text-xs mr-1 mb-1">${enzyme}</span>`
            ).join('');
        }
        
        // Update analysis info for output analysis
        const optionsContainer = document.getElementById('report-cleaning-options');
        if (data.total_sites === 0) {
            optionsContainer.innerHTML = `
                <div class="text-green-300">🎉 Perfect Cleaning!</div>
                <div class="text-green-400">All target sites removed</div>
                <div class="text-green-400">No remaining restrictions</div>
                <div class="text-green-300">Ready for cloning</div>
            `;
        } else {
            optionsContainer.innerHTML = `
                <div class="text-yellow-300">⚠️ Sites Remaining</div>
                <div class="text-yellow-400">${data.total_sites} sites still present</div>
                <div class="text-yellow-400">May need additional cleaning</div>
                <div class="text-blue-300">Output sequence analysis</div>
            `;
        }
        
        // Update sites container with detailed information
        const sitesContainer = document.getElementById('report-sites-container');
        let sitesHtml = '';
        
        if (!data.sites_found || Object.keys(data.sites_found).length === 0) {
            // No sites found - show success message
            sitesHtml = `
                <div class="bg-green-800/30 border border-green-500/30 rounded-lg p-4 text-center">
                    <div class="text-green-400 text-lg mb-2">🎉 Perfect Clean Sequence!</div>
                    <div class="text-green-300 text-sm mb-2">No restriction sites detected in the output sequence</div>
                    <div class="text-gray-300 text-xs">
                        ✅ Ready for cloning<br>
                        ✅ Compatible with standard vectors<br>
                        ✅ No unwanted cuts during digestion
                    </div>
                </div>
            `;
        } else {
            // Sites found - show detailed information
            Object.entries(data.sites_found).forEach(([enzyme, sites]) => {
            sitesHtml += `
                <div class="bg-slate-800/30 rounded-lg p-3">
                    <div class="flex items-center justify-between mb-2">
                        <h6 class="text-white font-medium">${enzyme}</h6>
                        <div class="flex gap-2">
                            <span class="bg-blue-500/20 text-blue-300 px-2 py-1 rounded text-xs">${sites.length} sites found</span>
                            <span class="bg-gray-500/20 text-gray-300 px-2 py-1 rounded text-xs">Recognition: ${sites[0]?.sequence || 'N/A'}</span>
                        </div>
                    </div>
                    
                    ${sites.length > 0 ? `
                        <div class="space-y-1">
                            <div class="text-xs text-gray-300 font-medium">Site Locations:</div>
                            <div class="max-h-32 overflow-y-auto">
                                ${sites.slice(0, 10).map(site => `
                                    <div class="text-xs bg-slate-900/50 rounded p-2 font-mono">
                                        <span class="text-gray-400">Position ${site.position}:</span>
                                        <span class="text-blue-300">${site.sequence}</span>
                                        <span class="text-gray-400 ml-2">(${site.cut_type || 'Standard cut'})</span>
                                    </div>
                                `).join('')}
                                ${sites.length > 10 ? `<div class="text-xs text-gray-400 text-center">... and ${sites.length - 10} more sites</div>` : ''}
                            </div>
                        </div>
                    ` : '<div class="text-xs text-gray-400">No sites found</div>'}
                </div>
            `;
            });
        }
        
        sitesContainer.innerHTML = sitesHtml;
        
        // Update benefits text for output analysis
        if (data.total_sites === 0) {
            document.getElementById('protein-function-text').textContent = 'Cleaning successful - sequence ready for cloning';
            document.getElementById('expression-text').textContent = 'No restriction sites to interfere with vectors';
            document.getElementById('regulatory-text').textContent = 'Clean sequence compatible with standard protocols';
        } else {
            document.getElementById('protein-function-text').textContent = `${data.total_sites} sites remain - may need additional cleaning`;
            document.getElementById('expression-text').textContent = 'Remaining sites may interfere with cloning';
            document.getElementById('regulatory-text').textContent = 'Consider additional cleaning steps';
        }
    }

    // Display cleaning results with detailed changes
    function displayCleaningResults(data) {
        // Hide empty state
        emptyOutput.style.display = 'none';
        
        // Show cleaned sequence
        if (data.cleaned_sequence) {
            const formattedCleanedSeq = formatSequenceForDisplay(data.cleaned_sequence);
            cleanedResult.innerHTML = formatSequenceWithColors(formattedCleanedSeq);
            cleanedOutput.style.display = 'block';
            
            // Update header to show it's cleaned
            const cleanedHeader = document.querySelector('#cleaned-output h4');
            if (cleanedHeader) {
                cleanedHeader.textContent = '🧬 Cleaned DNA Sequence';
                cleanedHeader.className = 'text-green-400 font-medium';
            }
        }
        
        // Update GC and Tm indicators for cleaned sequence
        if (data.sequence_analysis) {
            updateResultsIndicators(data.sequence_analysis.gc_content, data.sequence_analysis.melting_temp);
            sitesRemovedCount.textContent = `${data.statistics.sites_removed} sites removed`;
        }
        
        // Show detailed statistics
        if (data.statistics) {
            displayDetailedStats(data);
            statsOutput.style.display = 'block';
        }
        
        // Show detailed cleaning report
        if (data.changes_made && data.changes_made.length > 0) {
            displayCleaningReport(data);
            reportOutput.style.display = 'block';
        }
    }

    // Display detailed statistics
    function displayDetailedStats(data) {
        const stats = data.statistics;
        const analysis = data.sequence_analysis;
        
        // Update statistics cards
        document.getElementById('stats-sites-found').textContent = stats.total_sites_found;
        document.getElementById('stats-sites-removed').textContent = stats.sites_removed;
        document.getElementById('stats-enzymes-processed').textContent = stats.enzymes_processed;
        document.getElementById('stats-final-length').textContent = analysis.length;
        
        // Update sequence properties
        document.getElementById('stats-gc-content').textContent = `${analysis.gc_content}%`;
        document.getElementById('stats-melting-temp').textContent = `${analysis.melting_temp}°C`;
        document.getElementById('stats-molecular-weight').textContent = `${analysis.molecular_weight.toLocaleString()} Da`;
        
        // Update base composition
        document.getElementById('stats-base-a').textContent = analysis.base_composition.A;
        document.getElementById('stats-base-t').textContent = analysis.base_composition.T;
        document.getElementById('stats-base-g').textContent = analysis.base_composition.G;
        document.getElementById('stats-base-c').textContent = analysis.base_composition.C;
    }

    // Display detailed cleaning report
    function displayCleaningReport(data) {
        // Update enzymes targeted
        const enzymesContainer = document.getElementById('report-enzymes-targeted');
        enzymesContainer.innerHTML = data.enzymes_processed.map(enzyme => 
            `<span class="inline-block bg-red-500/20 text-red-300 px-2 py-1 rounded text-xs mr-1 mb-1">${enzyme}</span>`
        ).join('');
        
        // Update cleaning options
        const optionsContainer = document.getElementById('report-cleaning-options');
        optionsContainer.innerHTML = `
            <div class="${data.options_used.preserve_orfs ? 'text-green-300' : 'text-gray-500'}">${data.options_used.preserve_orfs ? '✅' : '❌'} Preserve ORFs</div>
            <div class="${data.options_used.silent_mutations ? 'text-green-300' : 'text-gray-500'}">${data.options_used.silent_mutations ? '✅' : '❌'} Silent Mutations</div>
            <div class="${data.options_used.optimize_codons ? 'text-green-300' : 'text-gray-500'}">${data.options_used.optimize_codons ? '✅' : '❌'} Codon Optimization</div>
            <div class="${data.options_used.preserve_regulatory ? 'text-green-300' : 'text-gray-500'}">${data.options_used.preserve_regulatory ? '✅' : '❌'} Preserve Regulatory</div>
        `;
        
        // Update sites container
        const sitesContainer = document.getElementById('report-sites-container');
        let sitesHtml = '';
        
        Object.entries(data.sites_found).forEach(([enzyme, sites]) => {
            const removedSites = data.changes_made.filter(change => change.enzyme === enzyme);
            
            sitesHtml += `
                <div class="bg-slate-800/30 rounded-lg p-3">
                    <div class="flex items-center justify-between mb-2">
                        <h6 class="text-white font-medium">${enzyme}</h6>
                        <div class="flex gap-2">
                            <span class="bg-red-500/20 text-red-300 px-2 py-1 rounded text-xs">${sites.length} found</span>
                            <span class="bg-green-500/20 text-green-300 px-2 py-1 rounded text-xs">${removedSites.length} removed</span>
                        </div>
                    </div>
                    <div class="text-xs text-gray-400 mb-2">Recognition site: <span class="font-mono text-white">${sites[0]?.sequence || 'N/A'}</span></div>
                    
                    ${removedSites.length > 0 ? `
                        <div class="space-y-1">
                            <div class="text-xs text-gray-300 font-medium">Changes Made:</div>
                            ${removedSites.slice(0, 5).map(change => `
                                <div class="text-xs bg-slate-900/50 rounded p-2 font-mono">
                                    <span class="text-gray-400">Position ${change.position}:</span>
                                    <span class="text-red-300">${change.original}</span>
                                    <span class="text-gray-400">→</span>
                                    <span class="text-green-300">${change.mutated}</span>
                                    <span class="text-blue-300 ml-2">(${change.change_type})</span>
                                </div>
                            `).join('')}
                            ${removedSites.length > 5 ? `<div class="text-xs text-gray-400 text-center">... and ${removedSites.length - 5} more changes</div>` : ''}
                        </div>
                    ` : '<div class="text-xs text-gray-400">No changes needed</div>'}
                </div>
            `;
        });
        
        sitesContainer.innerHTML = sitesHtml;
        
        // Update benefits text
        document.getElementById('protein-function-text').textContent = 
            data.options_used.silent_mutations ? 'Silent mutations preserve protein sequence' : 'Mutations may affect protein sequence';
        document.getElementById('expression-text').textContent = 
            data.options_used.optimize_codons ? 'Codon usage optimized for expression' : 'Original codon usage maintained';
        document.getElementById('regulatory-text').textContent = 
            data.options_used.preserve_regulatory ? 'Important regulatory sequences preserved' : 'Regulatory sequences may be affected';
    }

    // Update results indicators (GC content and Tm)
    function updateResultsIndicators(gcContent, meltingTemp) {
        if (outputGC) {
            outputGC.textContent = `${gcContent}%`;
        }
        if (outputTm) {
            outputTm.textContent = `${meltingTemp}°C`;
        }
    }

    // Copy cleaned sequence
    copyCleanedBtn.addEventListener('click', function() {
        const rawText = cleanedResult.textContent;
        console.log('Raw restriction text:', rawText);
        
        // Extract only the DNA sequence (remove line numbers, labels, etc.)
        const cleanSequence = rawText
            .replace(/^\d+\|\s*/gm, '') // Remove line numbers like "00001| "
            .replace(/5'\s*-\s*3'/g, '') // Remove strand indicators
            .replace(/3'\s*-\s*5'/g, '') // Remove reverse strand indicators
            .replace(/\s+/g, '') // Remove all whitespace
            .replace(/[^ATGCNRYSWKMBDHV]/gi, ''); // Keep only valid DNA bases
        
        console.log('Cleaned restriction sequence:', cleanSequence);
        
        if (cleanSequence) {
            navigator.clipboard.writeText(cleanSequence).then(() => {
                console.log('Restriction copy successful!');
                this.textContent = '✅ Copied!';
                setTimeout(() => {
                    this.textContent = '📋 Copy';
                }, 2000);
            }).catch(err => {
                console.error('Restriction copy failed:', err);
                showNotification('Failed to copy sequence', 'error');
            });
        } else {
            console.log('No valid DNA sequence found in restriction cleaner');
            showNotification('No valid DNA sequence found', 'error');
        }
    });

    // Download cleaned sequence
    downloadCleanedBtn.addEventListener('click', function() {
        const rawText = cleanedResult.textContent;
        console.log('Raw restriction download text:', rawText);
        
        // Extract only the DNA sequence (remove line numbers, labels, etc.)
        const cleanSequence = rawText
            .replace(/^\d+\|\s*/gm, '') // Remove line numbers like "00001| "
            .replace(/5'\s*-\s*3'/g, '') // Remove strand indicators
            .replace(/3'\s*-\s*5'/g, '') // Remove reverse strand indicators
            .replace(/\s+/g, '') // Remove all whitespace
            .replace(/[^ATGCNRYSWKMBDHV]/gi, ''); // Keep only valid DNA bases
        
        console.log('Clean restriction sequence for download:', cleanSequence);
        
        if (cleanSequence) {
            const blob = new Blob([`>Cleaned_Sequence\n${cleanSequence}`], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'cleaned_sequence.fasta';
            a.click();
            URL.revokeObjectURL(url);
            console.log('Restriction download initiated');
        } else {
            console.log('No valid DNA sequence found for restriction download');
            showNotification('No valid DNA sequence found', 'error');
        }
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

    // Edit Lock button functionality
    editLockBtn.addEventListener('click', function() {
        toggleEditLockMode();
    });

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
            }
            if (line.startsWith('//')) {
                break;
            }
            if (inSequence) {
                // Remove numbers and spaces from sequence lines
                const cleanLine = line.replace(/[^a-zA-Z]/g, '');
                sequence += cleanLine;
            }
        }
        return sequence;
    }

    // Load sequence using Python backend
    async function loadSequence(sequence) {
        try {
            loadBtn.disabled = true;
            loadBtn.textContent = '🔄 Loading...';
            
            const result = await callPythonScript('load', sequence);
            
            if (result.success) {
                displayLoadedSequence(result.data);
            } else {
                throw new Error(result.error || 'Failed to load sequence');
            }
            
        } catch (error) {
            console.error('Error:', error);
            alert('Error loading sequence: ' + error.message);
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
            cleanedResult.innerHTML = formatSequenceWithColors(data.formatted_sequence);
            cleanedOutput.style.display = 'block';
            
            // Update the header to show it's a loaded sequence
            const cleanedHeader = document.querySelector('#cleaned-output h4');
            if (cleanedHeader) {
                cleanedHeader.textContent = 'Loaded DNA Sequence';
                cleanedHeader.className = 'text-green-400 font-medium';
            }
            
            // Update sites removed count to show sequence info
            sitesRemovedCount.textContent = `${data.length} bp loaded`;
            sitesRemovedCount.className = 'px-2 py-1 bg-green-600/20 border border-green-500/30 rounded text-green-300 text-xs';
        }
        
        // Show basic statistics using new HTML structure
        if (data.base_composition) {
            // Update statistics cards (set to 0 for loaded sequence, will be updated after analysis)
            document.getElementById('stats-sites-found').textContent = '0';
            document.getElementById('stats-sites-removed').textContent = '0';
            document.getElementById('stats-enzymes-processed').textContent = '0';
            document.getElementById('stats-final-length').textContent = data.length;
            
            // Update sequence properties
            document.getElementById('stats-gc-content').textContent = `${data.gc_content}%`;
            document.getElementById('stats-melting-temp').textContent = `${data.melting_temp}°C`;
            document.getElementById('stats-molecular-weight').textContent = `${data.molecular_weight.toLocaleString()} Da`;
            
            // Update base composition
            document.getElementById('stats-base-a').textContent = data.base_composition.A;
            document.getElementById('stats-base-t').textContent = data.base_composition.T;
            document.getElementById('stats-base-g').textContent = data.base_composition.G;
            document.getElementById('stats-base-c').textContent = data.base_composition.C;
            
            statsOutput.style.display = 'block';
            
            // Update stats header
            const statsHeader = document.querySelector('#stats-output h4');
            if (statsHeader) {
                statsHeader.textContent = 'Sequence Analysis';
                statsHeader.className = 'text-cyan-400 font-medium';
            }
        }
    }



    // Format sequence with line breaks
    function formatSequence(sequence) {
        const lineLength = 80;
        let formatted = '';
        for (let i = 0; i < sequence.length; i += lineLength) {
            formatted += sequence.slice(i, i + lineLength) + '\n';
        }
        return formatted.trim();
    }

    // Format sequence with line numbers (JavaScript version of Python function)
    function formatSequenceForDisplay(sequence, lineLength = 140) {
        let formatted = '';
        for (let i = 0; i < sequence.length; i += lineLength) {
            const lineNum = i + 1;
            const lineSeq = sequence.substring(i, i + lineLength);
            // Zero-pad line numbers to 5 digits for consistent alignment
            const paddedLineNum = lineNum.toString().padStart(5, '0');
            formatted += `${paddedLineNum}|${lineSeq}\n`;
        }
        return formatted.trim();
    }

    // Format sequence with colors and line numbers like DNA Sequencer
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
                
                // Format with line numbers and strand indicators like DNA Sequencer
                // Separate non-editable elements from editable DNA sequence
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
        reportOutput.style.display = 'none';
        emptyOutput.style.display = 'block';
    }

    // Toggle Edit Lock Mode (exact copy from DNA Sequencer)
    function toggleEditLockMode() {
        const restrictionOutput = document.getElementById('cleaned-result');
        
        if (!editLockBtn || !restrictionOutput) {
            showNotification('Edit mode not available', 'error');
            return;
        }

        // Check if we're currently in edit mode
        const isEditMode = restrictionOutput.contentEditable === 'true';
        
        if (isEditMode) {
            // Exit edit mode - clean removal
            restrictionOutput.contentEditable = 'false';
            restrictionOutput.style.border = '';
            editLockBtn.innerHTML = '🔒 Edit Lock';
            editLockBtn.className = 'px-3 py-1 bg-orange-600/20 border border-orange-500/30 rounded text-orange-300 hover:bg-orange-600/30 transition-all text-sm';
            
            // Update the current sequence with edited content
            updateSequenceFromEdit();
            
            showNotification('Edit mode disabled - sequence locked', 'info');
        } else {
            // Enter edit mode - exact copy from DNA Sequencer
            if (!restrictionOutput.innerHTML || restrictionOutput.innerHTML.trim() === '') {
                showNotification('No sequence to edit - load a sequence first', 'error');
                return;
            }
            
            restrictionOutput.contentEditable = 'true';
            restrictionOutput.style.border = '2px solid #f59e0b';
            editLockBtn.innerHTML = '🔓 Exit Edit';
            editLockBtn.className = 'px-3 py-1 bg-green-600/20 border border-green-500/30 rounded text-green-300 hover:bg-green-600/30 transition-all text-sm';
            
            showNotification('Edit mode enabled - you can now freely edit the DNA sequence!', 'success');
            
            // Focus on the editable area
            restrictionOutput.focus();
        }
    }

    // Update sequence from edit (only extract DNA from editable spans)
    function updateSequenceFromEdit() {
        const restrictionOutput = document.getElementById('cleaned-result');
        if (!restrictionOutput) return;
        
        // Extract text only from editable sequence spans, not line numbers or strand indicators
        const editableSpans = restrictionOutput.querySelectorAll('.editable-sequence');
        let editedText = '';
        
        editableSpans.forEach(span => {
            const spanText = span.innerText || span.textContent || '';
            editedText += spanText;
        });
        
        // Clean up the text - remove non-DNA characters and normalize
        editedText = editedText.replace(/[^ATGCNRYSWKMBDHV\s]/gi, '').replace(/\s+/g, '').toUpperCase();
        
        if (editedText.length > 0) {
            // Update the GC content and Tm displays
            updateResultsIndicatorsFromSequence(editedText);
            
            // Reformat and redisplay the edited sequence with colors
            const formattedEditedSeq = formatSequenceForDisplay(editedText);
            restrictionOutput.innerHTML = formatSequenceWithColors(formattedEditedSeq);
            
            showNotification(`Sequence updated! Length: ${editedText.length} bp`, 'success');
        }
    }

    // Update results indicators from sequence (for edit mode)
    function updateResultsIndicatorsFromSequence(sequence) {
        const cleanSequence = sequence.replace(/\s/g, '').toUpperCase();
        
        if (cleanSequence.length > 0) {
            // Calculate GC content
            const gcCount = (cleanSequence.match(/[GC]/g) || []).length;
            const gcPercent = ((gcCount / cleanSequence.length) * 100).toFixed(1);
            outputGC.textContent = `${gcPercent}%`;
            
            // Calculate approximate melting temperature (basic formula)
            const gcContent = gcCount / cleanSequence.length;
            let meltingTemp;
            
            if (cleanSequence.length < 14) {
                // Short sequence formula
                const atCount = (cleanSequence.match(/[AT]/g) || []).length;
                meltingTemp = 4 * gcCount + 2 * atCount;
            } else {
                // Longer sequence formula
                meltingTemp = 64.9 + 41 * (gcContent - 0.164);
            }
            
            outputTm.textContent = `${Math.round(meltingTemp)}°C`;
        } else {
            outputGC.textContent = '0%';
            outputTm.textContent = '--°C';
        }
    }

    // Show notification (copied from DNA Sequencer)
    function showNotification(message, type = 'info') {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `fixed top-4 right-4 z-50 px-4 py-3 rounded-lg shadow-lg transition-all duration-300 transform translate-x-full`;
        
        // Set colors based on type
        switch(type) {
            case 'success':
                notification.className += ' bg-green-600/90 border border-green-500/50 text-green-100';
                break;
            case 'error':
                notification.className += ' bg-red-600/90 border border-red-500/50 text-red-100';
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

    // Call Python script using child_process.spawn
    async function callPythonScript(operation, sequence, enzymeSelection = null, options = null) {
        return new Promise((resolve, reject) => {
            const { spawn } = require('child_process');
            const path = require('path');

            // Prepare arguments for Python script
            const scriptPath = path.join(__dirname, 'assets', 'restriction-clearner.py');
            const args = [scriptPath, operation, sequence];
            
            // Add additional arguments for analyze operation
            if (operation === 'analyze' && enzymeSelection) {
                // Ensure enzymeSelection is a string
                const enzymeSelectionStr = String(enzymeSelection).trim();
                args.push(enzymeSelectionStr);
                
                if (options) {
                    const optionsStr = String(options);
                    args.push(optionsStr);
                }
            }

            console.log('🐍 Calling Restriction Cleaner Python script:', scriptPath);
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
                } catch (parseError) {
                    console.error('❌ Failed to parse Python output:', parseError);
                    console.error('📤 Raw output:', stdout);
                    reject(new Error(`Failed to parse Python output: ${parseError.message}`));
                }
            });

            pythonProcess.on('error', (error) => {
                console.error('❌ Failed to start Python process:', error);
                reject(new Error(`Failed to start Python process: ${error.message}`));
            });
        });
    }

    // Initialize
    updateSequenceStats();
});