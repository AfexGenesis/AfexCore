// ===== AFEXGENESIS™ RENDERER SCRIPT =====
// Main renderer process script for genetic laboratory interface
// Handles UI interactions, Python script execution, and data visualization

class AfexGenesisRenderer {
    constructor() {
        this.isElectron = typeof window.electronAPI !== 'undefined';
        this.currentTool = null;
        this.sequenceData = new Map();
        this.analysisResults = new Map();
        this.charts = new Map();
        this.init();
    }

    init() {
        console.log('AfexGenesis™ Renderer Initializing...');
        
        // Initialize core systems
        this.initializeEventListeners();
        this.initializeGeneticTools();
        this.initializeFileHandlers();
        this.initializeNotificationSystem();
        
        // Setup UI components
        this.setupProgressIndicators();
        this.setupTooltips();
        this.setupKeyboardShortcuts();
        
        console.log('AfexGenesis™ Renderer Ready');
    }

    // ===== EVENT LISTENERS =====
    initializeEventListeners() {
        // Tool launch buttons
        document.addEventListener('click', (e) => {
            if (e.target.classList.contains('genetic-launch-btn')) {
                const toolName = e.target.getAttribute('data-page');
                this.launchGeneticTool(toolName);
            }
        });

        // File input handlers
        document.addEventListener('change', (e) => {
            if (e.target.type === 'file') {
                this.handleFileInput(e.target);
            }
        });

        // Form submissions
        document.addEventListener('submit', (e) => {
            if (e.target.classList.contains('genetic-form')) {
                e.preventDefault();
                this.handleGeneticFormSubmission(e.target);
            }
        });

        // Real-time sequence validation
        document.addEventListener('input', (e) => {
            if (e.target.classList.contains('sequence-input')) {
                this.validateSequenceInput(e.target);
            }
        });
    }

    // ===== GENETIC TOOLS INITIALIZATION =====
    initializeGeneticTools() {
        this.geneticTools = {
            'dna-sequencer': {
                name: 'DNA Sequencer',
                script: 'dna-sequencer.py',
                icon: '🧬',
                description: 'Advanced genetic sequence analysis and reconstruction'
            },
            'dna-translator': {
                name: 'DNA Translator',
                script: 'dna-translator-direct.py',
                icon: '🔄',
                description: 'Translate DNA into RNA and proteins'
            },
            'gene-extractor': {
                name: 'Gene Extractor',
                script: 'gene-extractor.py',
                icon: '🧪',
                description: 'Extract genes from genome files'
            },
            'codon-finder': {
                name: 'Codon Finder',
                script: 'codon-finder.py',
                icon: '🔍',
                description: 'Find and analyze codons in sequences'
            },
            'genome-comparator': {
                name: 'Genome Comparator',
                script: 'gene-comparator.py',
                icon: '⚖️',
                description: 'Compare genetic sequences'
            },
            'reverse-transcriber': {
                name: 'Reverse Transcriber',
                script: 'reverse-transcriber.py',
                icon: '↩️',
                description: 'Convert RNA back to DNA'
            },
            'codon-optimizer': {
                name: 'Codon Optimizer',
                script: 'codon-optimizer.py',
                icon: '⚡',
                description: 'Optimize codon usage for target organisms'
            },
            'chromosome-splitter': {
                name: 'Chromosome Splitter',
                script: 'chromosome-splitter.py',
                icon: '✂️',
                description: 'Split multi-chromosome genome files'
            },
            'base-counter': {
                name: 'Base Counter',
                script: 'base-counter.py',
                icon: '📊',
                description: 'Count nucleotide bases and calculate GC content'
            },
            'restriction-site-cleaner': {
                name: 'Restriction Site Cleaner',
                script: 'restriction-cleaner.py',
                icon: '🧹',
                description: 'Remove unwanted restriction sites'
            },
            'snp-highlighter': {
                name: 'SNP Highlighter',
                script: 'snp-highlighter.py',
                icon: '🎯',
                description: 'Highlight single nucleotide polymorphisms'
            },
            'orf-cleaner': {
                name: 'ORF Cleaner',
                script: 'orf-cleaner.py',
                icon: '🔧',
                description: 'Clean and optimize open reading frames'
            }
        };
    }

    // ===== TOOL LAUNCHING =====
    async launchGeneticTool(toolName) {
        const tool = this.geneticTools[toolName];
        if (!tool) {
            this.showError(`Unknown tool: ${toolName}`);
            return;
        }

        this.currentTool = toolName;
        this.showNotification(`Launching ${tool.name}`, `${tool.icon} ${tool.description}`);
        
        // Show tool interface
        this.showToolInterface(toolName);
    }

    showToolInterface(toolName) {
        // Create dynamic tool interface
        const toolContainer = document.getElementById('tool-container') || this.createToolContainer();
        const tool = this.geneticTools[toolName];
        
        toolContainer.innerHTML = `
            <div class="genetic-tool-interface" data-tool="${toolName}">
                <div class="tool-header">
                    <div class="flex items-center gap-3">
                        <span class="text-2xl">${tool.icon}</span>
                        <h2 class="text-2xl font-bold text-white">${tool.name}</h2>
                    </div>
                    <button class="close-tool-btn" onclick="afexRenderer.closeTool()">
                        <svg class="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M6 18L18 6M6 6l12 12"></path>
                        </svg>
                    </button>
                </div>
                
                <div class="tool-content">
                    ${this.generateToolContent(toolName)}
                </div>
                
                <div class="tool-results" id="tool-results-${toolName}" style="display: none;">
                    <!-- Results will be displayed here -->
                </div>
            </div>
        `;
        
        toolContainer.style.display = 'block';
    }

    generateToolContent(toolName) {
        const commonInputs = `
            <div class="input-group">
                <label class="input-label">Input Sequence</label>
                <textarea class="sequence-input" id="sequence-input-${toolName}" 
                         placeholder="Enter your genetic sequence here (DNA, RNA, or protein)..."
                         rows="6"></textarea>
                <div class="sequence-info" id="sequence-info-${toolName}"></div>
            </div>
        `;

        const fileInput = `
            <div class="input-group">
                <label class="input-label">Or Upload File</label>
                <input type="file" class="file-input" id="file-input-${toolName}" 
                       accept=".fasta,.fa,.gb,.gbk,.txt" />
                <div class="file-info" id="file-info-${toolName}"></div>
            </div>
        `;

        const runButton = `
            <div class="action-group">
                <button class="run-analysis-btn" onclick="afexRenderer.runAnalysis('${toolName}')">
                    <span class="btn-icon">🚀</span>
                    Run Analysis
                </button>
            </div>
        `;

        // Tool-specific options
        let specificOptions = '';
        switch (toolName) {
            case 'dna-translator':
                specificOptions = `
                    <div class="input-group">
                        <label class="input-label">Translation Type</label>
                        <select id="translation-type-${toolName}" class="select-input">
                            <option value="dna-to-rna">DNA to RNA</option>
                            <option value="dna-to-protein">DNA to Protein</option>
                            <option value="rna-to-protein">RNA to Protein</option>
                        </select>
                    </div>
                `;
                break;
            case 'codon-optimizer':
                specificOptions = `
                    <div class="input-group">
                        <label class="input-label">Target Organism</label>
                        <select id="target-organism-${toolName}" class="select-input">
                            <option value="human">Human</option>
                            <option value="mouse">Mouse</option>
                            <option value="ecoli">E. coli</option>
                            <option value="yeast">Yeast</option>
                            <option value="drosophila">Drosophila</option>
                        </select>
                    </div>
                `;
                break;
            case 'genome-comparator':
                specificOptions = `
                    <div class="input-group">
                        <label class="input-label">Reference Sequence</label>
                        <textarea class="sequence-input" id="reference-sequence-${toolName}" 
                                 placeholder="Enter reference sequence for comparison..."
                                 rows="4"></textarea>
                    </div>
                `;
                break;
        }

        return `
            <form class="genetic-form" data-tool="${toolName}">
                ${commonInputs}
                ${fileInput}
                ${specificOptions}
                ${runButton}
            </form>
        `;
    }

    createToolContainer() {
        const container = document.createElement('div');
        container.id = 'tool-container';
        container.className = 'tool-container';
        container.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0, 0, 0, 0.9);
            backdrop-filter: blur(10px);
            z-index: 1000;
            display: none;
            overflow-y: auto;
            padding: 2rem;
        `;
        document.body.appendChild(container);
        return container;
    }

    closeTool() {
        const toolContainer = document.getElementById('tool-container');
        if (toolContainer) {
            toolContainer.style.display = 'none';
        }
        this.currentTool = null;
    }

    // ===== ANALYSIS EXECUTION =====
    async runAnalysis(toolName) {
        if (!this.isElectron) {
            this.showError('This feature requires the Electron desktop application');
            return;
        }

        const tool = this.geneticTools[toolName];
        const sequenceInput = document.getElementById(`sequence-input-${toolName}`);
        const sequence = sequenceInput.value.trim();

        if (!sequence) {
            this.showError('Please enter a sequence or upload a file');
            return;
        }

        // Validate sequence
        if (!this.validateSequence(sequence)) {
            this.showError('Invalid sequence format. Please check your input.');
            return;
        }

        this.showProgress(`Running ${tool.name}...`);

        try {
            let result;
            const options = this.getToolOptions(toolName);

            // Execute appropriate Python script
            switch (toolName) {
                case 'dna-sequencer':
                    result = await window.electronAPI.runDNASequencer(sequence, options);
                    break;
                case 'dna-translator':
                    const translationType = document.getElementById(`translation-type-${toolName}`).value;
                    result = await window.electronAPI.runDNATranslator(sequence, translationType);
                    break;
                case 'gene-extractor':
                    result = await window.electronAPI.runGeneExtractor(sequence, options);
                    break;
                case 'codon-finder':
                    result = await window.electronAPI.runCodonFinder(sequence);
                    break;
                case 'genome-comparator':
                    const refSequence = document.getElementById(`reference-sequence-${toolName}`).value;
                    result = await window.electronAPI.runGeneComparator(sequence, refSequence, options);
                    break;
                case 'reverse-transcriber':
                    result = await window.electronAPI.runReverseTranscriber(sequence, options);
                    break;
                case 'codon-optimizer':
                    const targetOrganism = document.getElementById(`target-organism-${toolName}`).value;
                    result = await window.electronAPI.runCodonOptimizer(sequence, targetOrganism);
                    break;
                case 'chromosome-splitter':
                    result = await window.electronAPI.runChromosomeSplitter(sequence, options);
                    break;
                case 'base-counter':
                    result = await window.electronAPI.runBaseCounter(sequence);
                    break;
                case 'restriction-site-cleaner':
                    result = await window.electronAPI.runRestrictionCleaner(sequence, options.restrictionSites);
                    break;
                case 'snp-highlighter':
                    const refSeq = document.getElementById(`reference-sequence-${toolName}`)?.value || '';
                    result = await window.electronAPI.runSNPHighlighter(sequence, refSeq, options);
                    break;
                case 'orf-cleaner':
                    result = await window.electronAPI.runORFCleaner(sequence, options);
                    break;
                default:
                    throw new Error(`Unknown tool: ${toolName}`);
            }

            this.hideProgress();
            this.displayResults(toolName, result);
            this.showNotification('Analysis Complete', `${tool.icon} Results are ready!`);

        } catch (error) {
            this.hideProgress();
            this.showError(`Analysis failed: ${error.message}`);
            console.error('Analysis error:', error);
        }
    }

    getToolOptions(toolName) {
        // Extract tool-specific options from the form
        const options = {};
        const form = document.querySelector(`[data-tool="${toolName}"] form`);
        
        if (form) {
            const inputs = form.querySelectorAll('input, select, textarea');
            inputs.forEach(input => {
                if (input.id && input.value && !input.id.includes('sequence-input')) {
                    const key = input.id.replace(`-${toolName}`, '').replace(/-/g, '_');
                    options[key] = input.value;
                }
            });
        }
        
        return options;
    }

    // ===== RESULTS DISPLAY =====
    displayResults(toolName, result) {
        const resultsContainer = document.getElementById(`tool-results-${toolName}`);
        if (!resultsContainer) return;

        try {
            const parsedResult = typeof result === 'string' ? JSON.parse(result) : result;
            
            resultsContainer.innerHTML = `
                <div class="results-header">
                    <h3 class="text-xl font-bold text-white mb-4">Analysis Results</h3>
                </div>
                <div class="results-content">
                    ${this.formatResults(toolName, parsedResult)}
                </div>
                <div class="results-actions">
                    <button onclick="afexRenderer.exportResults('${toolName}')" class="export-btn">
                        📄 Export Results
                    </button>
                    <button onclick="afexRenderer.saveResults('${toolName}')" class="save-btn">
                        💾 Save Analysis
                    </button>
                </div>
            `;
            
            resultsContainer.style.display = 'block';
            
            // Create visualizations if applicable
            this.createVisualization(toolName, parsedResult);
            
        } catch (error) {
            resultsContainer.innerHTML = `
                <div class="results-error">
                    <h3 class="text-xl font-bold text-red-400 mb-4">Raw Output</h3>
                    <pre class="bg-gray-800 p-4 rounded text-green-400 overflow-auto">${result}</pre>
                </div>
            `;
            resultsContainer.style.display = 'block';
        }
    }

    formatResults(toolName, results) {
        // Format results based on tool type
        switch (toolName) {
            case 'base-counter':
                return this.formatBaseCounterResults(results);
            case 'dna-translator':
                return this.formatTranslationResults(results);
            case 'codon-finder':
                return this.formatCodonResults(results);
            case 'genome-comparator':
                return this.formatComparisonResults(results);
            default:
                return `<pre class="bg-gray-800 p-4 rounded text-green-400 overflow-auto">${JSON.stringify(results, null, 2)}</pre>`;
        }
    }

    formatBaseCounterResults(results) {
        if (results.base_counts) {
            const { A, T, C, G } = results.base_counts;
            const total = A + T + C + G;
            const gcContent = results.gc_content || ((G + C) / total * 100);
            
            return `
                <div class="base-counter-results">
                    <div class="stats-grid">
                        <div class="stat-card">
                            <div class="stat-value">${A}</div>
                            <div class="stat-label">Adenine (A)</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">${T}</div>
                            <div class="stat-label">Thymine (T)</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">${C}</div>
                            <div class="stat-label">Cytosine (C)</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">${G}</div>
                            <div class="stat-label">Guanine (G)</div>
                        </div>
                    </div>
                    <div class="gc-content">
                        <h4 class="text-lg font-semibold text-white mb-2">GC Content</h4>
                        <div class="gc-bar">
                            <div class="gc-fill" style="width: ${gcContent}%"></div>
                        </div>
                        <div class="gc-percentage">${gcContent.toFixed(2)}%</div>
                    </div>
                </div>
            `;
        }
        return `<pre>${JSON.stringify(results, null, 2)}</pre>`;
    }

    // ===== FILE HANDLING =====
    initializeFileHandlers() {
        // Handle drag and drop
        document.addEventListener('dragover', (e) => {
            e.preventDefault();
            e.stopPropagation();
        });

        document.addEventListener('drop', (e) => {
            e.preventDefault();
            e.stopPropagation();
            
            const files = Array.from(e.dataTransfer.files);
            this.handleDroppedFiles(files);
        });
    }

    async handleDroppedFiles(files) {
        for (const file of files) {
            if (this.isGeneticFile(file)) {
                const content = await this.readFileContent(file);
                this.processGeneticFile(file.name, content);
            }
        }
    }

    isGeneticFile(file) {
        const validExtensions = ['.fasta', '.fa', '.gb', '.gbk', '.txt', '.seq'];
        return validExtensions.some(ext => file.name.toLowerCase().endsWith(ext));
    }

    async readFileContent(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = (e) => resolve(e.target.result);
            reader.onerror = (e) => reject(e);
            reader.readAsText(file);
        });
    }

    processGeneticFile(filename, content) {
        // Process FASTA, GenBank, or plain text files
        const sequences = this.parseGeneticFile(content);
        this.sequenceData.set(filename, sequences);
        
        this.showNotification('File Loaded', `${filename} loaded with ${sequences.length} sequence(s)`);
    }

    parseGeneticFile(content) {
        const sequences = [];
        
        if (content.startsWith('>')) {
            // FASTA format
            const entries = content.split('>').filter(entry => entry.trim());
            entries.forEach(entry => {
                const lines = entry.trim().split('\n');
                const header = lines[0];
                const sequence = lines.slice(1).join('').replace(/\s/g, '');
                sequences.push({ header, sequence });
            });
        } else {
            // Plain text or other format
            const cleanSequence = content.replace(/[^ATCGURYKMSWBDHVN]/gi, '');
            if (cleanSequence.length > 0) {
                sequences.push({ header: 'Imported Sequence', sequence: cleanSequence });
            }
        }
        
        return sequences;
    }

    // ===== VALIDATION =====
    validateSequence(sequence) {
        const cleanSequence = sequence.replace(/\s/g, '').toUpperCase();
        
        // Check for valid nucleotide characters
        const validDNA = /^[ATCG]+$/;
        const validRNA = /^[AUCG]+$/;
        const validProtein = /^[ACDEFGHIKLMNPQRSTVWY]+$/;
        const validAmbiguous = /^[ATCGURYKMSWBDHVN]+$/;
        
        return validDNA.test(cleanSequence) || 
               validRNA.test(cleanSequence) || 
               validProtein.test(cleanSequence) ||
               validAmbiguous.test(cleanSequence);
    }

    validateSequenceInput(input) {
        const sequence = input.value.trim();
        const infoElement = document.getElementById(input.id.replace('input', 'info'));
        
        if (!infoElement) return;
        
        if (sequence.length === 0) {
            infoElement.innerHTML = '';
            return;
        }
        
        const isValid = this.validateSequence(sequence);
        const cleanSequence = sequence.replace(/\s/g, '');
        const length = cleanSequence.length;
        
        infoElement.innerHTML = `
            <div class="sequence-validation ${isValid ? 'valid' : 'invalid'}">
                <span class="validation-icon">${isValid ? '✅' : '❌'}</span>
                <span class="sequence-length">Length: ${length} bp</span>
                <span class="validation-status">${isValid ? 'Valid sequence' : 'Invalid characters detected'}</span>
            </div>
        `;
    }

    // ===== UI HELPERS =====
    setupProgressIndicators() {
        const progressContainer = document.createElement('div');
        progressContainer.id = 'progress-container';
        progressContainer.className = 'progress-container hidden';
        progressContainer.innerHTML = `
            <div class="progress-modal">
                <div class="progress-content">
                    <div class="progress-spinner"></div>
                    <div class="progress-text" id="progress-text">Processing...</div>
                </div>
            </div>
        `;
        document.body.appendChild(progressContainer);
    }

    showProgress(message) {
        const container = document.getElementById('progress-container');
        const text = document.getElementById('progress-text');
        if (container && text) {
            text.textContent = message;
            container.classList.remove('hidden');
        }
    }

    hideProgress() {
        const container = document.getElementById('progress-container');
        if (container) {
            container.classList.add('hidden');
        }
    }

    // ===== NOTIFICATIONS =====
    initializeNotificationSystem() {
        const notificationContainer = document.createElement('div');
        notificationContainer.id = 'notification-container';
        notificationContainer.className = 'notification-container';
        document.body.appendChild(notificationContainer);
    }

    showNotification(title, message, type = 'info') {
        const container = document.getElementById('notification-container');
        if (!container) return;

        const notification = document.createElement('div');
        notification.className = `notification notification-${type}`;
        notification.innerHTML = `
            <div class="notification-content">
                <div class="notification-title">${title}</div>
                <div class="notification-message">${message}</div>
            </div>
            <button class="notification-close" onclick="this.parentElement.remove()">×</button>
        `;

        container.appendChild(notification);

        // Auto-remove after 5 seconds
        setTimeout(() => {
            if (notification.parentElement) {
                notification.remove();
            }
        }, 5000);
    }

    showError(message) {
        this.showNotification('Error', message, 'error');
    }

    // ===== KEYBOARD SHORTCUTS =====
    setupKeyboardShortcuts() {
        document.addEventListener('keydown', (e) => {
            if (e.ctrlKey || e.metaKey) {
                switch (e.key) {
                    case 'o':
                        e.preventDefault();
                        this.openFile();
                        break;
                    case 's':
                        e.preventDefault();
                        this.saveCurrentAnalysis();
                        break;
                    case 'Escape':
                        this.closeTool();
                        break;
                }
            }
        });
    }

    // ===== EXPORT/SAVE FUNCTIONS =====
    async exportResults(toolName) {
        if (!this.isElectron) {
            this.showError('Export feature requires the desktop application');
            return;
        }

        try {
            const results = this.analysisResults.get(toolName);
            if (!results) {
                this.showError('No results to export');
                return;
            }

            const options = {
                title: 'Export Analysis Results',
                defaultPath: `${toolName}-results.json`,
                filters: [
                    { name: 'JSON Files', extensions: ['json'] },
                    { name: 'Text Files', extensions: ['txt'] },
                    { name: 'All Files', extensions: ['*'] }
                ]
            };

            const result = await window.electronAPI.showSaveDialog(options);
            if (!result.canceled) {
                const content = JSON.stringify(results, null, 2);
                await window.electronAPI.writeFile(result.filePath, content);
                this.showNotification('Export Complete', `Results saved to ${result.filePath}`);
            }
        } catch (error) {
            this.showError(`Export failed: ${error.message}`);
        }
    }

    // ===== VISUALIZATION =====
    createVisualization(toolName, results) {
        // Create charts and visualizations based on results
        if (toolName === 'base-counter' && results.base_counts) {
            this.createBaseCompositionChart(toolName, results.base_counts);
        }
    }

    createBaseCompositionChart(toolName, baseCounts) {
        const chartContainer = document.createElement('div');
        chartContainer.innerHTML = `
            <div class="chart-container">
                <h4 class="chart-title">Base Composition</h4>
                <canvas id="chart-${toolName}" width="400" height="200"></canvas>
            </div>
        `;
        
        const resultsContainer = document.getElementById(`tool-results-${toolName}`);
        if (resultsContainer) {
            resultsContainer.appendChild(chartContainer);
            
            // Create Chart.js chart if available
            if (typeof Chart !== 'undefined') {
                const ctx = document.getElementById(`chart-${toolName}`).getContext('2d');
                new Chart(ctx, {
                    type: 'doughnut',
                    data: {
                        labels: ['Adenine (A)', 'Thymine (T)', 'Cytosine (C)', 'Guanine (G)'],
                        datasets: [{
                            data: [baseCounts.A, baseCounts.T, baseCounts.C, baseCounts.G],
                            backgroundColor: ['#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0']
                        }]
                    },
                    options: {
                        responsive: true,
                        plugins: {
                            legend: {
                                position: 'bottom'
                            }
                        }
                    }
                });
            }
        }
    }

    // ===== TOOLTIPS =====
    setupTooltips() {
        document.addEventListener('mouseover', (e) => {
            if (e.target.hasAttribute('data-tooltip')) {
                this.showTooltip(e.target, e.target.getAttribute('data-tooltip'));
            }
        });

        document.addEventListener('mouseout', (e) => {
            if (e.target.hasAttribute('data-tooltip')) {
                this.hideTooltip();
            }
        });
    }

    showTooltip(element, text) {
        const tooltip = document.createElement('div');
        tooltip.className = 'tooltip';
        tooltip.textContent = text;
        document.body.appendChild(tooltip);

        const rect = element.getBoundingClientRect();
        tooltip.style.left = rect.left + (rect.width / 2) - (tooltip.offsetWidth / 2) + 'px';
        tooltip.style.top = rect.top - tooltip.offsetHeight - 10 + 'px';
    }

    hideTooltip() {
        const tooltip = document.querySelector('.tooltip');
        if (tooltip) {
            tooltip.remove();
        }
    }
}

// ===== INITIALIZATION =====
// Initialize the renderer when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    window.afexRenderer = new AfexGenesisRenderer();
});

// ===== CONSOLE LOGGING =====
console.log(`
🧬 ===============================================
🔬     AFEXGENESIS™ RENDERER SCRIPT LOADED
🧪       Genetic Laboratory Interface Ready
🦕         Advanced Analysis Tools Available
===============================================
`);