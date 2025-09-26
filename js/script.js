// ===== AFEXGENESIS™ - GENETIC LABORATORY SOFTWARE =====
// Advanced Genetic Engineering Interface with WebView Integration

class AfexGenesisApp {
    constructor() {
        this.currentPage = 'home';
        this.activeOverlay = null;
        this.isResizing = false; // Track global resize state
        this.apiBaseUrl = 'http://localhost:5000/api';
        this.projectToDelete = null; // Store project data for deletion
        this.init();
    }

    init() {
        // Initialize file system handler
        this.bindSearchEvents();
        this.bindGeneticSearchEvents();
        this.bindResearchSearchEvents();
        this.bindChatBoxEvents();
        this.bindFinancialGridEvents();
        this.bindGeneticToolEvents();
        
        // Set initial page
        this.showPage('home');
        
        
        // Global click handler to close status dropdowns
        this.bindGlobalClickHandler();
    }


    // ===== GLOBAL EVENT HANDLERS =====
    bindGlobalClickHandler() {
        document.addEventListener('click', (e) => {
            // Close status dropdowns when clicking outside
            if (!e.target.closest('.status-circle-container')) {
                document.querySelectorAll('.status-dropdown').forEach(dropdown => {
                    dropdown.classList.add('hidden');
                });
            }
        });
    }



    showPage(pageId) {
        // Remove active class from all buttons
        document.querySelectorAll('.sidebar-btn[data-page]').forEach(btn => {
            btn.classList.remove('active');
        });
        
        // Add active class to clicked button (if it exists)
        const targetButton = document.querySelector(`[data-page="${pageId}"]`);
        if (targetButton) {
            targetButton.classList.add('active');
        }
        
        // Hide all pages (both main pages and financial pages)
        this.hideAllPages();
        
        // Show selected page
        const selectedPage = document.getElementById(`${pageId}-page`);
        console.log(`🔧 Looking for page: ${pageId}-page`);
        console.log(`🔧 Page found:`, !!selectedPage);
        if (selectedPage) {
            selectedPage.classList.add('active');
            selectedPage.style.display = 'block';
            console.log(`✅ Showing page: ${pageId}-page`);
        } else {
            console.error(`❌ Page not found: ${pageId}-page`);
        }
        
        // Control visibility of homepage-only elements
        this.toggleHomepageElements(pageId === 'home');
        
        // Load projects when switching to research page
        if (pageId === 'research') {
            console.log('📊 Switching to research page, loading projects...');
            this.loadProjects();
        }
        
        // Load purchase orders when switching to invoice generator page
        if (pageId === 'invoice-generator') {
            console.log('📄 Switching to invoice generator page, loading purchase orders...');
            // Call the function from invoice.js
            if (typeof loadPurchaseOrdersFromBackend === 'function') {
                loadPurchaseOrdersFromBackend();
            } else {
                console.error('❌ loadPurchaseOrdersFromBackend function not found');
            }
        }
        
        // Initialize ExpenseManager when switching to expense tracker page
        if (pageId === 'expense-tracker') {
            console.log('💰 Switching to expense tracker page, initializing ExpenseManager...');
            // Call the function from expense-manager.js
            if (typeof initializeExpenseManagerForPage === 'function') {
                initializeExpenseManagerForPage();
            } else {
                console.error('❌ initializeExpenseManagerForPage function not found');
            }
        }
        
        this.currentPage = pageId;
        console.log(`📄 Navigated to: ${pageId.toUpperCase()}`);
    }
    
    toggleHomepageElements(show) {
        const searchContainer = document.getElementById('search-bar-container');
        const chatboxContainer = document.getElementById('chatbox-container');
        const financialGridContainer = document.getElementById('financial-grid-container');
        
        const display = show ? 'block' : 'none';
        
        if (searchContainer) searchContainer.style.display = display;
        if (chatboxContainer) chatboxContainer.style.display = display;
        if (financialGridContainer) financialGridContainer.style.display = display;
        
        console.log(`🏠 Homepage elements ${show ? 'shown' : 'hidden'}`);
    }


    // ===== UTILITY METHODS =====
    getCurrentPage() {
        return this.currentPage;
    }

    getActiveOverlay() {
        return this.activeOverlay;
    }

    // ===== SEARCH FUNCTIONALITY =====
    bindSearchEvents() {
        const searchInput = document.getElementById('search-input');
        const searchClear = document.getElementById('search-clear');
        const searchResults = document.getElementById('search-results');
        
        if (!searchInput) return;
        
        let searchTimeout;
        
        // Search input event
        searchInput.addEventListener('input', (e) => {
            const query = e.target.value.trim();
            
            // Show/hide clear button
            searchClear.style.display = query ? 'flex' : 'none';
            
            // Clear previous timeout
            clearTimeout(searchTimeout);
            
            if (query.length === 0) {
                this.hideSearchResults();
                return;
            }
            
            // Debounce search
            searchTimeout = setTimeout(() => {
                this.performSearch(query);
            }, 300);
        });
        
        // Clear button event
        searchClear.addEventListener('click', () => {
            searchInput.value = '';
            searchClear.style.display = 'none';
            this.hideSearchResults();
            searchInput.focus();
        });
        
        // Hide results when clicking outside
        document.addEventListener('click', (e) => {
            if (!e.target.closest('.search-bar')) {
                this.hideSearchResults();
            }
        });
        
        // Keyboard navigation
        searchInput.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                this.hideSearchResults();
                searchInput.blur();
            }
        });
    }
    
    performSearch(query) {
        const searchResults = document.getElementById('search-results');
        if (!searchResults) return;
        
        // Sample search data - you can expand this with real data
        const searchData = [
            {
                title: 'DNA Sequencing Analysis',
                description: 'Advanced genetic sequencing tools and analysis methods',
                category: 'Genetic Lab',
                action: () => this.showPage('genetic')
            },
            {
                title: 'Research Documentation',
                description: 'Laboratory research papers and experimental data',
                category: 'Research',
                action: () => this.showPage('research')
            },
            {
                title: 'ChatGPT AI Assistant',
                description: 'AI-powered research assistant for genetic analysis',
                category: 'AI Tools',
                action: () => this.openOverlay('chatgpt')
            },
            {
                title: 'Gemini AI Analysis',
                description: 'Advanced AI for complex genetic pattern recognition',
                category: 'AI Tools',
                action: () => this.openOverlay('gemini')
            },
            {
                title: 'Custom AI Laboratory',
                description: 'Specialized AI tools for genetic engineering',
                category: 'AI Tools',
                action: () => this.openOverlay('ai')
            },
            {
                title: 'Web Browser Research',
                description: 'Access external genetic databases and research',
                category: 'Tools',
                action: () => this.openOverlay('browser')
            },
            {
                title: 'Genetic Engineering Protocols',
                description: 'Step-by-step genetic modification procedures',
                category: 'Genetic Lab',
                action: () => this.showPage('genetic')
            },
            {
                title: 'Laboratory Equipment Status',
                description: 'Monitor and control laboratory instruments',
                category: 'Research',
                action: () => this.showPage('research')
            }
        ];
        
        // Filter results based on query
        const filteredResults = searchData.filter(item => 
            item.title.toLowerCase().includes(query.toLowerCase()) ||
            item.description.toLowerCase().includes(query.toLowerCase()) ||
            item.category.toLowerCase().includes(query.toLowerCase())
        );
        
        this.displaySearchResults(filteredResults);
    }
    
    displaySearchResults(results) {
        const searchResults = document.getElementById('search-results');
        if (!searchResults) return;
        
        if (results.length === 0) {
            searchResults.innerHTML = `
                <div class="search-result-item">
                    <div class="search-result-title">No results found</div>
                    <div class="search-result-description">Try searching for genetic, research, AI, or laboratory terms</div>
                </div>
            `;
        } else {
            searchResults.innerHTML = results.map(result => `
                <div class="search-result-item" data-action="${results.indexOf(result)}">
                    <div class="search-result-title">${result.title}</div>
                    <div class="search-result-description">${result.description}</div>
                    <div class="search-result-category">${result.category}</div>
                </div>
            `).join('');
            
            // Add click handlers
            searchResults.querySelectorAll('.search-result-item').forEach((item, index) => {
                if (results[index] && results[index].action) {
                    item.addEventListener('click', () => {
                        results[index].action();
                        this.hideSearchResults();
                        document.getElementById('search-input').value = '';
                        document.getElementById('search-clear').style.display = 'none';
                    });
                }
            });
        }
        
        searchResults.style.display = 'block';
    }
    
    hideSearchResults() {
        const searchResults = document.getElementById('search-results');
        if (searchResults) {
            searchResults.style.display = 'none';
        }
    }

    // ===== GENETIC PAGE SEARCH FUNCTIONALITY =====
    bindGeneticSearchEvents() {
        const geneticSearchInput = document.getElementById('genetic-search-input');
        const geneticSearchClear = document.getElementById('genetic-search-clear');
        
        if (!geneticSearchInput) return;
        
        let searchTimeout;
        
        // Search input event
        geneticSearchInput.addEventListener('input', (e) => {
            const query = e.target.value.trim();
            
            // Show/hide clear button
            geneticSearchClear.style.display = query ? 'flex' : 'none';
            
            // Clear previous timeout
            clearTimeout(searchTimeout);
            
            // Debounce search
            searchTimeout = setTimeout(() => {
                this.filterGeneticCards(query);
            }, 200);
        });
        
        // Clear button event
        geneticSearchClear.addEventListener('click', () => {
            geneticSearchInput.value = '';
            geneticSearchClear.style.display = 'none';
            this.filterGeneticCards('');
            geneticSearchInput.focus();
        });
        
        // Keyboard navigation
        geneticSearchInput.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                geneticSearchInput.value = '';
                geneticSearchClear.style.display = 'none';
                this.filterGeneticCards('');
                geneticSearchInput.blur();
            }
        });
    }
    
    filterGeneticCards(query) {
        const geneticCards = document.querySelectorAll('#home-page .grid > div');
        
        if (!query) {
            // Show all cards when no search query
            geneticCards.forEach(card => {
                card.style.display = 'block';
            });
            return;
        }
        
        const searchQuery = query.toLowerCase();
        
        geneticCards.forEach(card => {
            const title = card.querySelector('h3')?.textContent.toLowerCase() || '';
            const description = card.querySelector('p')?.textContent.toLowerCase() || '';
            
            // Check if the card matches the search query
            const matches = title.includes(searchQuery) || description.includes(searchQuery);
            
            // Show or hide the card based on match
            card.style.display = matches ? 'block' : 'none';
        });
        
        console.log(`🔍 Genetic search: "${query}" - ${document.querySelectorAll('#home-page .grid > div[style*="block"]').length} results found`);
    }

    // ===== RESEARCH PAGE SEARCH FUNCTIONALITY =====
    bindResearchSearchEvents() {
        const researchSearchInput = document.getElementById('research-search-input');
        const researchSearchClear = document.getElementById('research-search-clear');
        
        if (!researchSearchInput) return;
        
        let searchTimeout;
        
        // Search input event
        researchSearchInput.addEventListener('input', (e) => {
            const query = e.target.value.trim();
            
            // Show/hide clear button
            researchSearchClear.style.display = query ? 'flex' : 'none';
            
            // Clear previous timeout
            clearTimeout(searchTimeout);
            
            // Debounce search
            searchTimeout = setTimeout(() => {
                this.filterResearchProjects(query);
            }, 200);
        });
        
        // Clear button event
        researchSearchClear.addEventListener('click', () => {
            researchSearchInput.value = '';
            researchSearchClear.style.display = 'none';
            this.filterResearchProjects('');
            researchSearchInput.focus();
        });
        
        // Keyboard navigation
        researchSearchInput.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                researchSearchInput.value = '';
                researchSearchClear.style.display = 'none';
                this.filterResearchProjects('');
                researchSearchInput.blur();
            }
        });
    }
    
    filterResearchProjects(query) {
        const projectCards = document.querySelectorAll('#research-projects-grid > div');
        
        if (!query) {
            // Show all cards when no search query
            projectCards.forEach(card => {
                card.style.display = 'block';
            });
            return;
        }
        
        const searchQuery = query.toLowerCase();
        
        projectCards.forEach(card => {
            const title = card.querySelector('h3')?.textContent.toLowerCase() || '';
            const description = card.querySelector('p')?.textContent.toLowerCase() || '';
            const owner = card.querySelector('.text-sm .text-gray-400')?.textContent.toLowerCase() || '';
            
            // Check if the card matches the search query
            const matches = title.includes(searchQuery) || 
                          description.includes(searchQuery) || 
                          owner.includes(searchQuery);
            
            // Show or hide the card based on match
            card.style.display = matches ? 'block' : 'none';
        });
        
        console.log(`🔍 Research search: "${query}" - ${document.querySelectorAll('#research-projects-grid > div[style*="block"]').length} results found`);
    }

    // ===== PROJECT CREATION FUNCTIONALITY =====
    bindProjectCreationEvents() {
        const createProjectCard = document.getElementById('create-project-card');
        const projectModal = document.getElementById('project-modal');
        const closeModalBtn = document.getElementById('close-project-modal');
        const cancelBtn = document.getElementById('cancel-project');
        const projectForm = document.getElementById('project-form');
        const datetimeInput = document.getElementById('project-datetime');
        
        if (!createProjectCard || !projectModal) return;
        
        // Open modal when clicking create project card
        createProjectCard.addEventListener('click', () => {
            this.openProjectModal();
        });
        
        // Close modal events
        closeModalBtn.addEventListener('click', () => {
            this.closeProjectModal();
        });
        
        cancelBtn.addEventListener('click', () => {
            this.closeProjectModal();
        });
        
        // Close modal when clicking outside
        projectModal.addEventListener('click', (e) => {
            if (e.target === projectModal) {
                this.closeProjectModal();
            }
        });
        
        // Handle form submission
        projectForm.addEventListener('submit', (e) => {
            e.preventDefault();
            this.createNewProject();
        });
        
        // Update datetime when modal opens
        projectModal.addEventListener('transitionend', () => {
            if (!projectModal.classList.contains('hidden')) {
                this.updateDateTime();
            }
        });
    }
    
    openProjectModal() {
        const projectModal = document.getElementById('project-modal');
        const projectForm = document.getElementById('project-form');
        
        // Reset form
        projectForm.reset();
        
        // Update current date/time
        this.updateDateTime();
        
        // Show modal
        projectModal.classList.remove('hidden');
        
        // Focus on first input
        setTimeout(() => {
            document.getElementById('project-name').focus();
        }, 100);
        
        console.log('📝 Project creation modal opened');
    }
    
    closeProjectModal() {
        const projectModal = document.getElementById('project-modal');
        projectModal.classList.add('hidden');
        console.log('❌ Project creation modal closed');
    }
    
    updateDateTime() {
        const datetimeInput = document.getElementById('project-datetime');
        const now = new Date();
        const formattedDateTime = now.toLocaleString('en-US', {
            year: 'numeric',
            month: '2-digit',
            day: '2-digit',
            hour: '2-digit',
            minute: '2-digit',
            hour12: false
        }).replace(/(\d+)\/(\d+)\/(\d+),\s*(\d+):(\d+)/, '$3-$1-$2 $4:$5');
        
        datetimeInput.value = formattedDateTime;
    }
    
    async createNewProject() {
        console.log('🆕 Creating new project...');
        
        const projectName = document.getElementById('project-name').value.trim();
        const projectOwner = document.getElementById('project-owner').value.trim();
        const projectDescription = document.getElementById('project-description').value.trim();
        const projectDateTime = document.getElementById('project-datetime').value;
        
        console.log('🆕 Project details:', { projectName, projectOwner, projectDescription, projectDateTime });
        
        if (!projectName || !projectOwner) {
            console.error('❌ Missing required fields');
            alert('Please fill in all required fields.');
            return;
        }
        
        try {
            // Send project data to API
            const response = await fetch(`${this.apiBaseUrl}/projects`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    name: projectName,
                    owner: projectOwner,
                    description: projectDescription || 'No description provided',
                    creation_date: projectDateTime
                })
            });
            
            const result = await response.json();
            
            if (result.success) {
                // Create new project card with the returned data
                const newProjectCard = this.generateProjectCard({
                    id: result.id,
                    name: result.name,
                    owner: result.owner,
                    description: result.description,
                    creation_date: result.creation_date,
                    status: result.status
                });
                
                // Add to projects grid (after the create project card)
                const projectsGrid = document.getElementById('research-projects-grid');
                const createCard = document.getElementById('create-project-card');
                
                // Insert after the create card (create card should always be first)
                if (createCard.nextSibling) {
                    projectsGrid.insertBefore(newProjectCard, createCard.nextSibling);
                } else {
                    projectsGrid.appendChild(newProjectCard);
                }
                
                // Close modal
                this.closeProjectModal();
                
                // Show success message
                console.log(`✅ New project created: "${projectName}" by ${projectOwner}`);
                this.showNotification(`Project "${projectName}" created successfully!`, 'success');
            } else {
                console.error('❌ Failed to create project:', result.error);
                this.showNotification(`Failed to create project: ${result.error}`, 'error');
            }
        } catch (error) {
            console.error('❌ Error creating project:', error);
            this.showNotification('Failed to create project. Please make sure the API server is running.', 'error');
        }
    }
    
    generateProjectCard(project) {
        const colors = ['cyan', 'blue', 'purple', 'green', 'yellow', 'red', 'orange'];
        const randomColor = colors[Math.floor(Math.random() * colors.length)];
        
        const cardDiv = document.createElement('div');
        cardDiv.className = `bg-slate-800/30 backdrop-blur-sm border border-${randomColor}-500/20 rounded-xl p-6 hover-lift relative`;
        cardDiv.dataset.projectId = project.id;
        
        // Get status color and info
        const statusInfo = this.getStatusInfo(project.status);
        
        cardDiv.innerHTML = `
            <!-- Status Circle (Top Right) -->
            <div class="absolute top-4 right-4">
                <div class="status-circle-container relative">
                    <div class="status-circle w-6 h-6 rounded-full cursor-pointer transition-all hover:scale-110 ${statusInfo.bgClass} ${statusInfo.borderClass} border-2 flex items-center justify-center" 
                         data-project-id="${project.id}" data-current-status="${project.status}">
                        <div class="w-2 h-2 rounded-full ${statusInfo.dotClass}"></div>
                    </div>
                    <!-- Status Dropdown -->
                    <div class="status-dropdown absolute top-8 right-0 bg-slate-800/95 backdrop-blur-sm border border-slate-600/50 rounded-lg shadow-xl z-50 min-w-[140px] hidden">
                        <div class="p-2">
                            <div class="text-xs text-gray-400 mb-2 px-2">Change Status:</div>
                            ${this.generateStatusOptions(project.status)}
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="flex items-center gap-3 mb-4 pr-8">
                <img src="emojis/test_tube.svg" class="w-8 h-8">
                <h3 class="text-xl font-semibold text-white">${project.name}</h3>
            </div>
            <p class="text-gray-300 mb-2 min-h-[3rem] overflow-y-auto max-h-32">${project.description}</p>
            <div class="text-sm text-gray-400 mb-4">
                <div>Owner: ${project.owner}</div>
                <div>Created: ${project.creation_date}</div>
                <div class="flex items-center gap-2">
                    <span>Status:</span>
                    <span class="status-text ${statusInfo.textClass} font-medium">${project.status}</span>
                </div>
            </div>
            <div class="flex gap-2">
                <button class="flex-1 py-2 px-4 bg-${randomColor}-600/20 border border-${randomColor}-500/30 rounded-lg text-${randomColor}-100 hover:bg-${randomColor}-600/30 transition-all project-open-btn">
                    Open
                </button>
                <button class="flex-1 py-2 px-4 bg-red-600/20 border border-red-500/30 rounded-lg text-red-100 hover:bg-red-600/30 transition-all project-delete-btn">
                    Delete
                </button>
            </div>
        `;
        
        // Add event listeners for the buttons
        const openBtn = cardDiv.querySelector('.project-open-btn');
        const deleteBtn = cardDiv.querySelector('.project-delete-btn');
        const statusCircle = cardDiv.querySelector('.status-circle');
        const statusDropdown = cardDiv.querySelector('.status-dropdown');
        
        openBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            this.openProject(project);
        });
        
        deleteBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            this.showDeleteModal(project);
        });
        
        // Status circle click handler
        statusCircle.addEventListener('click', (e) => {
            e.stopPropagation();
            this.toggleStatusDropdown(statusDropdown);
        });
        
        // Status option click handlers
        const statusOptions = cardDiv.querySelectorAll('.status-option');
        statusOptions.forEach(option => {
            option.addEventListener('click', (e) => {
                e.stopPropagation();
                const newStatus = option.dataset.status;
                this.updateProjectStatus(project.id, newStatus, cardDiv);
                this.hideStatusDropdown(statusDropdown);
            });
        });
        
        return cardDiv;
    }
    
    // ===== STATUS MANAGEMENT =====
    getStatusInfo(status) {
        const statusMap = {
            'Planning': {
                bgClass: 'bg-blue-600/20',
                borderClass: 'border-blue-500/50',
                dotClass: 'bg-blue-400',
                textClass: 'text-blue-400'
            },
            'Active': {
                bgClass: 'bg-green-600/20',
                borderClass: 'border-green-500/50',
                dotClass: 'bg-green-400',
                textClass: 'text-green-400'
            },
            'In Progress': {
                bgClass: 'bg-yellow-600/20',
                borderClass: 'border-yellow-500/50',
                dotClass: 'bg-yellow-400',
                textClass: 'text-yellow-400'
            },
            'Failed': {
                bgClass: 'bg-red-600/20',
                borderClass: 'border-red-500/50',
                dotClass: 'bg-red-400',
                textClass: 'text-red-400'
            },
            'Success': {
                bgClass: 'bg-emerald-600/20',
                borderClass: 'border-emerald-500/50',
                dotClass: 'bg-emerald-400',
                textClass: 'text-emerald-400'
            },
            'Abandoned': {
                bgClass: 'bg-gray-600/20',
                borderClass: 'border-gray-500/50',
                dotClass: 'bg-gray-400',
                textClass: 'text-gray-400'
            }
        };
        
        return statusMap[status] || statusMap['Planning'];
    }
    
    generateStatusOptions(currentStatus) {
        const statuses = ['Planning', 'Active', 'In Progress', 'Failed', 'Success', 'Abandoned'];
        
        return statuses.map(status => {
            const statusInfo = this.getStatusInfo(status);
            const isCurrentStatus = status === currentStatus;
            
            return `
                <div class="status-option flex items-center gap-2 px-2 py-1.5 rounded cursor-pointer hover:bg-slate-700/50 transition-colors ${isCurrentStatus ? 'bg-slate-700/30' : ''}" 
                     data-status="${status}">
                    <div class="w-3 h-3 rounded-full ${statusInfo.bgClass} ${statusInfo.borderClass} border flex items-center justify-center">
                        <div class="w-1.5 h-1.5 rounded-full ${statusInfo.dotClass}"></div>
                    </div>
                    <span class="text-sm ${statusInfo.textClass} ${isCurrentStatus ? 'font-medium' : ''}">${status}</span>
                    ${isCurrentStatus ? '<span class="text-xs text-gray-500 ml-auto">✓</span>' : ''}
                </div>
            `;
        }).join('');
    }
    
    toggleStatusDropdown(dropdown) {
        // Hide all other dropdowns first
        document.querySelectorAll('.status-dropdown').forEach(dd => {
            if (dd !== dropdown) {
                dd.classList.add('hidden');
            }
        });
        
        // Toggle current dropdown
        dropdown.classList.toggle('hidden');
    }
    
    hideStatusDropdown(dropdown) {
        dropdown.classList.add('hidden');
    }
    
    async updateProjectStatus(projectId, newStatus, cardElement) {
        try {
            const response = await fetch(`${this.apiBaseUrl}/projects/${projectId}/status`, {
                method: 'PUT',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    status: newStatus
                })
            });
            
            if (!response.ok) {
                throw new Error(`HTTP ${response.status}: ${response.statusText}`);
            }
            
            const result = await response.json();
            
            if (result.success) {
                // Update the UI elements
                const statusInfo = this.getStatusInfo(newStatus);
                const statusCircle = cardElement.querySelector('.status-circle');
                const statusText = cardElement.querySelector('.status-text');
                const statusDropdown = cardElement.querySelector('.status-dropdown');
                
                // Update circle appearance
                statusCircle.className = `status-circle w-6 h-6 rounded-full cursor-pointer transition-all hover:scale-110 ${statusInfo.bgClass} ${statusInfo.borderClass} border-2 flex items-center justify-center`;
                statusCircle.dataset.currentStatus = newStatus;
                
                // Update circle dot
                const dot = statusCircle.querySelector('div');
                dot.className = `w-2 h-2 rounded-full ${statusInfo.dotClass}`;
                
                // Update status text
                statusText.className = `status-text ${statusInfo.textClass} font-medium`;
                statusText.textContent = newStatus;
                
                // Update dropdown options
                statusDropdown.querySelector('.p-2').innerHTML = `
                    <div class="text-xs text-gray-400 mb-2 px-2">Change Status:</div>
                    ${this.generateStatusOptions(newStatus)}
                `;
                
                // Re-bind click events for new options
                const newStatusOptions = statusDropdown.querySelectorAll('.status-option');
                newStatusOptions.forEach(option => {
                    option.addEventListener('click', (e) => {
                        e.stopPropagation();
                        const status = option.dataset.status;
                        this.updateProjectStatus(projectId, status, cardElement);
                        this.hideStatusDropdown(statusDropdown);
                    });
                });
                
                console.log(`✅ Project status updated to: ${newStatus}`);
                this.showNotification(`Project status updated to "${newStatus}"`, 'success');
            } else {
                console.error('❌ Failed to update project status:', result.error);
                this.showNotification(`Failed to update status: ${result.error}`, 'error');
            }
        } catch (error) {
            console.error('❌ Error updating project status:', error);
            this.showNotification('Failed to update status. Please make sure the API server is running.', 'error');
        }
    }
    

    
    openProject(project) {
        console.log(`🔓 Opening project: ${project.name}`);
        this.currentProject = project;
        this.showResearchEditor(project);
    }
    
    showDeleteModal(project) {
        this.projectToDelete = project;
        
        const deleteModal = document.getElementById('delete-modal');
        const projectNameDisplay = document.getElementById('project-name-display');
        const confirmationInput = document.getElementById('delete-confirmation-input');
        const confirmButton = document.getElementById('confirm-delete');
        const errorMessage = document.getElementById('delete-error-message');
        
        // Set project name in display
        projectNameDisplay.textContent = project.name;
        
        // Reset input and button state
        confirmationInput.value = '';
        confirmButton.disabled = true;
        errorMessage.classList.add('hidden');
        
        // Show modal
        deleteModal.classList.remove('hidden');
        
        // Focus on input
        setTimeout(() => {
            confirmationInput.focus();
        }, 100);
        
        console.log(`🗑️ Delete modal opened for project: ${project.name}`);
    }
    
    closeDeleteModal() {
        const deleteModal = document.getElementById('delete-modal');
        deleteModal.classList.add('hidden');
        this.projectToDelete = null;
        console.log('❌ Delete modal closed');
    }
    
    async deleteProject() {
        if (!this.projectToDelete) return;
        
        const projectName = this.projectToDelete.name;
        
        try {
            const response = await fetch(`${this.apiBaseUrl}/projects/${this.projectToDelete.id}`, {
                method: 'DELETE'
            });
            
            if (!response.ok) {
                throw new Error(`HTTP ${response.status}: ${response.statusText}`);
            }
            
            const result = await response.json();
            
            if (result.success) {
                // Remove the project card from the DOM
                const projectCard = document.querySelector(`[data-project-id="${this.projectToDelete.id}"]`);
                if (projectCard) {
                    projectCard.remove();
                }
                
                // Close modal
                this.closeDeleteModal();
                
                // Show success message
                console.log(`✅ Project deleted: ${projectName}`);
                this.showNotification(`Project "${projectName}" deleted successfully!`, 'success');
            } else {
                console.error('❌ Failed to delete project:', result.error);
                this.showNotification(`Failed to delete project: ${result.error}`, 'error');
            }
        } catch (error) {
            console.error('❌ Error deleting project:', error);
            
            // Check if the project was actually deleted by trying to remove it from DOM
            const projectCard = document.querySelector(`[data-project-id="${this.projectToDelete.id}"]`);
            if (projectCard) {
                // Project card still exists, so deletion likely failed
                this.showNotification('Failed to delete project. Please make sure the API server is running.', 'error');
            } else {
                // Project card is gone, so deletion might have succeeded despite the error
                this.closeDeleteModal();
                console.log(`⚠️ Project "${projectName}" may have been deleted despite network error`);
                this.showNotification(`Project "${projectName}" deleted (with network error)`, 'success');
            }
        }
    }
    
    // ===== DELETE MODAL EVENTS =====
    bindDeleteModalEvents() {
        const deleteModal = document.getElementById('delete-modal');
        const closeModalBtn = document.getElementById('close-delete-modal');
        const cancelBtn = document.getElementById('cancel-delete');
        const confirmBtn = document.getElementById('confirm-delete');
        const confirmationInput = document.getElementById('delete-confirmation-input');
        const errorMessage = document.getElementById('delete-error-message');
        
        if (!deleteModal) return;
        
        // Close modal events
        closeModalBtn.addEventListener('click', () => {
            this.closeDeleteModal();
        });
        
        cancelBtn.addEventListener('click', () => {
            this.closeDeleteModal();
        });
        
        // Close modal when clicking outside
        deleteModal.addEventListener('click', (e) => {
            if (e.target === deleteModal) {
                this.closeDeleteModal();
            }
        });
        
        // Validate input and enable/disable confirm button
        confirmationInput.addEventListener('input', () => {
            const inputValue = confirmationInput.value.trim();
            const projectName = this.projectToDelete ? this.projectToDelete.name : '';
            
            if (inputValue === projectName) {
                confirmBtn.disabled = false;
                errorMessage.classList.add('hidden');
            } else {
                confirmBtn.disabled = true;
                if (inputValue.length > 0) {
                    errorMessage.classList.remove('hidden');
                } else {
                    errorMessage.classList.add('hidden');
                }
            }
        });
        
        // Handle confirm delete
        confirmBtn.addEventListener('click', () => {
            this.deleteProject();
        });
        
        // Handle Enter key in input
        confirmationInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter' && !confirmBtn.disabled) {
                this.deleteProject();
            }
        });
    }
    
    showNotification(message, type = 'info') {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `fixed top-4 right-4 z-50 px-4 py-3 rounded-lg border backdrop-blur-sm transition-all duration-300 transform translate-x-full`;
        
        // Set colors based on type
        const colors = {
            success: 'bg-green-600/20 border-green-500/30 text-green-100',
            error: 'bg-red-600/20 border-red-500/30 text-red-100',
            info: 'bg-blue-600/20 border-blue-500/30 text-blue-100'
        };
        
        notification.className += ` ${colors[type] || colors.info}`;
        notification.textContent = message;
        
        // Add to document
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

    // ===== CHATBOX FUNCTIONALITY =====
    bindChatBoxEvents() {
        const chatboxInput = document.getElementById('chatbox-input');
        const chatboxSend = document.getElementById('chatbox-send');
        const chatboxMinimize = document.getElementById('chatbox-minimize');
        const chatboxMessages = document.getElementById('chatbox-messages');
        
        if (!chatboxInput || !chatboxSend) return;
        
        // Send message on button click
        chatboxSend.addEventListener('click', () => {
            this.sendChatMessage();
        });
        
        // Send message on Enter key
        chatboxInput.addEventListener('keydown', (e) => {
            if (e.key === 'Enter' && !e.shiftKey) {
                e.preventDefault();
                this.sendChatMessage();
            }
        });
        
        // Minimize/maximize chatbox
        if (chatboxMinimize) {
            chatboxMinimize.addEventListener('click', () => {
                this.toggleChatBox();
            });
        }
        
        // Auto-scroll to bottom when new messages are added
        if (chatboxMessages) {
            const observer = new MutationObserver(() => {
                chatboxMessages.scrollTop = chatboxMessages.scrollHeight;
            });
            observer.observe(chatboxMessages, { childList: true });
        }
    }
    
    sendChatMessage() {
        const chatboxInput = document.getElementById('chatbox-input');
        const message = chatboxInput.value.trim();
        
        if (!message) return;
        
        // Add user message
        this.addChatMessage(message, 'user');
        
        // Clear input
        chatboxInput.value = '';
        
        // Simulate bot response (you can replace this with actual API calls later)
        setTimeout(() => {
            this.addBotResponse(message);
        }, 1000);
    }
    
    addChatMessage(message, sender = 'user') {
        const chatboxMessages = document.getElementById('chatbox-messages');
        if (!chatboxMessages) return;
        
        const messageElement = document.createElement('div');
        messageElement.className = `chat-message ${sender}-message`;
        
        const currentTime = new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' });
        
        const avatarIcon = sender === 'user' 
            ? '<circle cx="12" cy="8" r="5"></circle><path d="m20 21-16 0 2-4 12 0 2 4z"></path>'
            : '<circle cx="12" cy="8" r="5"></circle><path d="m20 21-16 0 2-4 12 0 2 4z"></path>';
        
        messageElement.innerHTML = `
            <div class="message-avatar">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    ${avatarIcon}
                </svg>
            </div>
            <div class="message-content">
                <div class="message-text">${message}</div>
                <div class="message-time">${currentTime}</div>
            </div>
        `;
        
        chatboxMessages.appendChild(messageElement);
        chatboxMessages.scrollTop = chatboxMessages.scrollHeight;
    }
    
    addBotResponse(userMessage) {
        // Simple response system - you can expand this with real AI integration
        const responses = {
            'hello': 'Hello! Welcome to the AFEX Genesis Lab. How can I assist you with your genetic research today?',
            'hi': 'Hi there! I\'m here to help you with any questions about genetic engineering, DNA analysis, or laboratory procedures.',
            'help': 'I can help you with:\n• DNA sequencing analysis\n• Genetic engineering protocols\n• Laboratory equipment guidance\n• Research documentation\n• AI tool recommendations',
            'dna': 'DNA analysis is one of our core capabilities! I can guide you through sequencing procedures, help interpret results, or recommend the best tools for your specific genetic research needs.',
            'research': 'For research assistance, I can help you organize your data, suggest methodologies, or connect you with relevant AI tools like ChatGPT or Gemini for advanced analysis.',
            'ai': 'Our lab has several AI assistants available:\n• ChatGPT for general research help\n• Gemini for complex pattern recognition\n• Custom AI tools for specialized genetic analysis',
            'default': 'That\'s an interesting question! While I\'m still learning about that specific topic, I can help you explore our lab resources or connect you with our specialized AI assistants for more detailed information.'
        };
        
        // Find the best response based on keywords
        let response = responses.default;
        const lowerMessage = userMessage.toLowerCase();
        
        for (const [keyword, reply] of Object.entries(responses)) {
            if (keyword !== 'default' && lowerMessage.includes(keyword)) {
                response = reply;
                break;
            }
        }
        
        this.addChatMessage(response, 'bot');
    }
    
    toggleChatBox() {
        const chatbox = document.querySelector('.chatbox');
        const minimizeBtn = document.getElementById('chatbox-minimize');
        
        if (!chatbox || !minimizeBtn) return;
        
        chatbox.classList.toggle('minimized');
        
        // Update minimize button icon
        const isMinimized = chatbox.classList.contains('minimized');
        minimizeBtn.innerHTML = isMinimized 
            ? '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="18,15 12,9 6,15"></polyline></svg>'
            : '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line></svg>';
    }

    // ===== FINANCIAL GRID FUNCTIONALITY =====
    bindFinancialGridEvents() {
        const gridItems = document.querySelectorAll('.grid-item');
        
        gridItems.forEach(item => {
            item.addEventListener('click', () => {
                const page = item.getAttribute('data-page');
                if (page) {
                    this.navigateToFinancialPage(page);
                }
            });
            
            // Add hover sound effect (optional)
            item.addEventListener('mouseenter', () => {
                // You can add a subtle hover sound here if desired
                item.style.transform = 'translateY(-2px)';
            });
            
            item.addEventListener('mouseleave', () => {
                item.style.transform = 'translateY(0)';
            });
        });
    }
    
    navigateToFinancialPage(page) {
        // Add a click animation
        const clickedItem = document.querySelector(`[data-page="${page}"]`);
        if (clickedItem) {
            clickedItem.style.transform = 'scale(0.95)';
            setTimeout(() => {
                clickedItem.style.transform = 'translateY(-2px)';
            }, 150);
        }
        
        // Hide current content and show the financial page
        this.hideAllPages();
        this.showFinancialPage(page);
        
        // Update URL without page reload
        window.history.pushState({page: page}, '', `#${page}`);
        
        console.log(`Navigating to: ${page}`);
    }
    
    showFinancialPage(page) {
        // Hide homepage elements when showing financial pages
        this.toggleHomepageElements(false);
        
        // Create or show the financial page content
        let pageContent = document.getElementById(`${page}-page`);
        
        if (!pageContent) {
            console.error(`Page content not found for: ${page}-page`);
            return;
        }
        
        pageContent.style.display = 'block';
        
        // Initialize specific page functionality
        if (page === 'expense-tracker') {
            // Initialize expense manager after a short delay to ensure DOM is ready
            setTimeout(() => {
                console.log('💰 Initializing ExpenseManager for expense tracker page...');
                if (typeof initializeExpenseManagerForPage === 'function') {
                    initializeExpenseManagerForPage();
                } else {
                    console.error('❌ initializeExpenseManagerForPage function not found');
                }
                
                // Also try the old function name for backward compatibility
                if (typeof initializeExpenseTracker === 'function') {
                    initializeExpenseTracker();
                }
            }, 100);
        } else if (page === 'revenue') {
            // Initialize revenue tracker charts after a short delay to ensure DOM is ready
            setTimeout(() => {
                if (window.revenueTracker) {
                    window.revenueTracker.refreshData();
                } else if (typeof RevenueTracker !== 'undefined') {
                    window.revenueTracker = new RevenueTracker();
                }
            }, 100);
        } else if (page === 'database-management') {
            // Initialize database management functionality
            setTimeout(() => {
                if (typeof initializeDatabaseManager === 'function') {
                    initializeDatabaseManager();
                }
            }, 100);
        }
        
        // Add fade-in animation
        setTimeout(() => {
            pageContent.style.opacity = '1';
        }, 50);
    }
    
    getFinancialPageInfo(page) {
        const pageData = {
            'expense-tracker': { icon: '💰', title: 'Expense Tracker' },
            'budget-planning': { icon: '📊', title: 'Budget Planning' },
            'revenue': { icon: '📈', title: 'Financial Company Revenue' },
            'invoice-generator': { icon: '📄', title: 'Invoice Generator' },
            'inventory-management': { icon: '📦', title: 'Inventory Management' },
            'asset-management': { icon: '🏢', title: 'Asset Management' },
            'tax-management': { icon: '🧾', title: 'Tax Management' },
            'accounts-payable': { icon: '💳', title: 'Accounts Payable' },
            'accounts-receivable': { icon: '💵', title: 'Accounts Receivable' },
            'purchase-order': { icon: '🛒', title: 'Purchase Order Management' },
            'database-management': { icon: '📂', title: 'Database Management' }
        };
        
        return pageData[page] || { icon: '❓', title: 'Unknown Page' };
    }
    
    goBack() {
        // Hide all financial pages
        const financialPages = document.querySelectorAll('.financial-page');
        financialPages.forEach(page => {
            page.style.opacity = '0';
            setTimeout(() => {
                page.style.display = 'none';
            }, 300);
        });
        
        // Show main content
        this.showPage('home');
        
        // Update URL
        window.history.pushState({page: 'home'}, '', '#home');
    }
    
    hideAllPages() {
        // Cleanup expense tracker if it's currently active
        if (typeof cleanupExpenseTracker === 'function') {
            cleanupExpenseTracker();
        }
        
        // Hide main pages (using both .page and .page-content classes)
        const pages = document.querySelectorAll('.page, .page-content');
        pages.forEach(page => {
            page.style.display = 'none';
            page.classList.remove('active');
        });
        
        // Hide financial pages
        const financialPages = document.querySelectorAll('.financial-page');
        financialPages.forEach(page => {
            page.style.display = 'none';
            page.style.opacity = '0';
        });
    }




    // ===== GENETIC TOOL FUNCTIONALITY =====
    bindGeneticToolEvents() {
        const geneticToolButtons = document.querySelectorAll('.genetic-launch-btn');
        
        geneticToolButtons.forEach(button => {
            button.addEventListener('click', () => {
                const page = button.getAttribute('data-page');
                if (page) {
                    this.navigateToGeneticTool(page);
                }
            });
            
            // Add hover effect
            button.addEventListener('mouseenter', () => {
                button.style.transform = 'translateY(-2px)';
            });
            
            button.addEventListener('mouseleave', () => {
                button.style.transform = 'translateY(0)';
            });
        });
    }
    
    navigateToGeneticTool(toolName) {
        // Add click animation
        const clickedButton = document.querySelector(`[data-page="${toolName}"]`);
        if (clickedButton) {
            clickedButton.style.transform = 'scale(0.95)';
            setTimeout(() => {
                clickedButton.style.transform = 'translateY(-2px)';
            }, 150);
        }
        
        // Hide all pages and show the genetic tool page
        this.hideAllPages();
        this.showGeneticToolPage(toolName);
        
        // Update URL without page reload
        window.history.pushState({page: toolName}, '', `#${toolName}`);
        
        console.log(`🧬 Navigating to genetic tool: ${toolName}`);
    }
    
    showGeneticToolPage(toolName) {
        // Hide homepage elements when showing genetic tool pages
        this.toggleHomepageElements(false);
        
        // Show the genetic tool page content
        const pageContent = document.getElementById(`${toolName}-page`);
        
        if (!pageContent) {
            console.error(`❌ Genetic tool page not found: ${toolName}-page`);
            return;
        }
        
        pageContent.style.display = 'block';
        pageContent.classList.add('active');
        
        console.log(`✅ Showing genetic tool page: ${toolName}-page`);
    }


}

// ===== SYSTEM INITIALIZATION =====
document.addEventListener('DOMContentLoaded', () => {
    // Initialize the main application
    window.afexGenesis = new AfexGenesisApp();
    
    // Make app globally available for genetic tools
    window.app = window.afexGenesis;
    

    // Performance monitoring
    console.log('⚡ Performance Metrics:', {
        loadTime: performance.now(),
        memory: performance.memory ? performance.memory.usedJSHeapSize : 'N/A',
        timestamp: new Date().toISOString()
    });
    
    // Welcome message
    console.log(`
    🧬 ===============================================
    🦕     AFEXGENESIS™ LABORATORY SYSTEM
    🔬         Genetic Engineering Interface
    🧪              System Status: ONLINE
    ===============================================
    🏠 Press Ctrl+1 for Home
    🧬 Press Ctrl+2 for Genetic Lab  
    🔬 Press Ctrl+3 for Research
    🌐 Press Ctrl+4 for Browser
    🤖 Press Ctrl+5 for ChatGPT
    💎 Press Ctrl+6 for Gemini
    🧠 Press Ctrl+7 for Custom AI
    ⌨️  Press ESC to close overlays
    ===============================================
    `);
});