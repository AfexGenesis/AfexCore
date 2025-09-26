const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
const path = require('path');
const fs = require('fs');

function createWindow() {
    // Create the browser window with WebView enabled
    const mainWindow = new BrowserWindow({
        width: 1200,
        height: 800,
        //frame: false, // remove default OS frame
        //titleBarStyle: 'hidden', // macOS cleaner look
        autoHideMenuBar: true, // Hide the menu bar (File, Edit, View, etc.)
        webPreferences: {
            nodeIntegration: true, // Disable for security - use preload script instead NVM it's on back again, hey if you are seeing this code that means u are either wondering how this Software works or Editing it
            contextIsolation: false, // Enable context isolation for security NVM it's disabled again. wait u again?  well i do gotta thank you for Using our software :) soo yeah while it might contain bugs we are working on a Faster & More efficient Genetic Engineering Tool which keeping it Casual Friendly GUI. soo yeah :)
            webviewTag: true, // Enable webview tag
            allowRunningInsecureContent: false,
            webSecurity: true,
            //enableRemoteModule: false, // Disable remote module for security
            //sandbox: false, // Keep false for your Python integration needs
            preload: path.join(__dirname, 'preload.js'), // Add preload script for secure API exposure
            experimentalFeatures: false,
            plugins: false,
            //webgl: false, // Disable WebGL if not needed
            //webaudio: false // Disable Web Audio if not needed
        },
        icon: path.join(__dirname, 'emojis/t-rex.png'),
        title: 'AfexCore™ Powered by AfexGenesis™ | The Future of De-Extinction',
        show: false
    });

    // Load the app
    mainWindow.loadFile('index.html');

    // Show window when ready
    mainWindow.once('ready-to-show', () => {
        mainWindow.show();
        console.log('🧬 AfexGenesis™ Laboratory System Launched');
    });

    // Open DevTools in development
    if (process.env.NODE_ENV === 'development') {
        mainWindow.webContents.openDevTools();
    }

    // Handle window closed
    mainWindow.on('closed', () => {
        console.log('🔬 Laboratory System Closed');
    });

    return mainWindow;
}

// App event handlers
app.whenReady().then(createWindow);

app.on('window-all-closed', () => {
    if (process.platform !== 'darwin') {
        app.quit();
    }
});

app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) {
        createWindow();
    }
});

// Enhanced Security: Comprehensive web contents protection
app.on('web-contents-created', (event, contents) => {
    // Prevent new window creation
    contents.on('new-window', (navigationEvent, navigationURL, frameName, disposition, options) => {
        navigationEvent.preventDefault();
        console.log('🚫 New window blocked:', navigationURL);
    });

    // Prevent navigation to external URLs
    contents.on('will-navigate', (navigationEvent, navigationURL) => {
        const allowedDomains = [
            'https://cdn.tailwindcss.com',
            'https://cdnjs.cloudflare.com',
            'https://cdn.jsdelivr.net',
            'https://fonts.googleapis.com',
            'https://fonts.gstatic.com'
        ];
        
        const isAllowed = allowedDomains.some(domain => navigationURL.startsWith(domain));
        
        if (!navigationURL.startsWith('file://') && !isAllowed) {
            navigationEvent.preventDefault();
            console.log('🚫 Navigation blocked:', navigationURL);
        }
    });

    // Prevent opening external links
    contents.setWindowOpenHandler(({ url }) => {
        console.log('🚫 External link blocked:', url);
        return { action: 'deny' };
    });

    // Security headers and CSP enforcement
    contents.session.webRequest.onHeadersReceived((details, callback) => {
        callback({
            responseHeaders: {
                ...details.responseHeaders,
                'X-Frame-Options': ['DENY'],
                'X-Content-Type-Options': ['nosniff'],
                'Referrer-Policy': ['no-referrer'],
                'Permissions-Policy': ['geolocation=(), microphone=(), camera=()']
            }
        });
    });
});

// Secure IPC Handlers
ipcMain.handle('dialog:openFile', async () => {
    const result = await dialog.showOpenDialog({
        properties: ['openFile'],
        filters: [
            { name: 'Genetic Files', extensions: ['fasta', 'fa', 'gb', 'genbank', 'txt'] },
            { name: 'All Files', extensions: ['*'] }
        ]
    });
    
    if (!result.canceled && result.filePaths.length > 0) {
        try {
            const content = fs.readFileSync(result.filePaths[0], 'utf8');
            return { success: true, content, filePath: result.filePaths[0] };
        } catch (error) {
            console.error('🚫 File read error:', error);
            return { success: false, error: error.message };
        }
    }
    return { success: false, error: 'No file selected' };
});

ipcMain.handle('dialog:saveFile', async (event, content) => {
    const result = await dialog.showSaveDialog({
        filters: [
            { name: 'Text Files', extensions: ['txt'] },
            { name: 'FASTA Files', extensions: ['fasta', 'fa'] },
            { name: 'All Files', extensions: ['*'] }
        ]
    });
    
    if (!result.canceled) {
        try {
            fs.writeFileSync(result.filePath, content);
            return { success: true, filePath: result.filePath };
        } catch (error) {
            console.error('🚫 File save error:', error);
            return { success: false, error: error.message };
        }
    }
    return { success: false, error: 'Save cancelled' };
});

ipcMain.handle('app:getVersion', () => {
    return app.getVersion();
});

ipcMain.handle('security:log', (event, securityEvent) => {
    console.log('🔒 Security Event:', securityEvent);
    // You could log to file or send to monitoring service here
    return true;
});

ipcMain.handle('shell:openExternal', async (event, url) => {
    // Only allow HTTPS URLs for security
    if (url.startsWith('https://')) {
        await shell.openExternal(url);
        return { success: true };
    }
    console.log('🚫 Blocked non-HTTPS external URL:', url);
    return { success: false, error: 'Only HTTPS URLs allowed' };
});

// Python execution handler (implement based on your needs)
ipcMain.handle('python:execute', async (event, script) => {
    // Implement secure Python execution here
    // This is a placeholder - you'll need to implement based on your Python integration
    console.log('🐍 Python execution requested:', script.substring(0, 100) + '...');
    return { success: true, result: 'Python execution placeholder' };
});

console.log(`
🧬 ===============================================
🦕     AFEXGENESIS™ ELECTRON APPLICATION
🔬         Genetic Engineering Interface
🧪           Starting Electron Process...
===============================================
`);