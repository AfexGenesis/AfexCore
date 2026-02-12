const { ipcRenderer } = require('electron');

// Since contextIsolation is disabled, we can directly expose APIs to window
// This works with your current webPreferences settings
window.electronAPI = {
    // File operations (secure)
    openFile: () => ipcRenderer.invoke('dialog:openFile'),
    saveFile: (content) => ipcRenderer.invoke('dialog:saveFile', content),
    
    // Python execution (if needed - implement in main process)
    executePython: (script) => ipcRenderer.invoke('python:execute', script),
    
    // System info
    getVersion: () => ipcRenderer.invoke('app:getVersion'),
    
    // Security logging
    logSecurityEvent: (event) => ipcRenderer.invoke('security:log', event),
    
    // Safe navigation
    openExternal: (url) => ipcRenderer.invoke('shell:openExternal', url)
};

// Security: Freeze the API to prevent tampering
Object.freeze(window.electronAPI);

console.log('Preload script loaded - Secure API bridge established');