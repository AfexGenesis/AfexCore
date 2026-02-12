// ===== PYTHON UTILITIES =====
// Helper functions for executing bundled Python scripts

class PythonUtils {
    static pythonPath = null;

    /**
     * Get the path to the bundled Python executable
     * Caches the result for performance
     */
    static async getPythonPath() {
        if (this.pythonPath) {
            return this.pythonPath;
        }

        try {
            const { ipcRenderer } = require('electron');
            this.pythonPath = await ipcRenderer.invoke('app:getPythonPath');
            console.log('Bundled Python path:', this.pythonPath);
            return this.pythonPath;
        } catch (error) {
            console.error('Failed to get Python path:', error);
            throw new Error('Unable to locate bundled Python executable');
        }
    }

    /**
     * Execute a Python script with the bundled Python
     * @param {string} scriptPath - Path to the Python script (relative to assets folder)
     * @param {string[]} args - Command line arguments to pass to the script
     * @param {Object} options - Additional spawn options
     * @returns {Promise<{success: boolean, output: string, error?: string}>}
     */
    static async executePython(scriptPath, args = [], options = {}) {
        return new Promise(async (resolve, reject) => {
            try {
                const { spawn } = require('child_process');
                const path = require('path');

                // Get the bundled Python executable path
                const pythonExe = await this.getPythonPath();

                // Resolve script path relative to assets folder
                const fullScriptPath = path.isAbsolute(scriptPath) 
                    ? scriptPath 
                    : path.join(__dirname, '..', 'assets', scriptPath);

                const pythonArgs = [fullScriptPath, ...args];

                console.log('Executing Python script:', fullScriptPath);
                console.log('Python executable:', pythonExe);
                console.log('Arguments:', pythonArgs);

                // Spawn the bundled Python process
                const pythonProcess = spawn(pythonExe, pythonArgs, {
                    cwd: path.join(__dirname, '..'),
                    stdio: ['pipe', 'pipe', 'pipe'],
                    ...options
                });

                let stdout = '';
                let stderr = '';

                pythonProcess.stdout.on('data', (data) => {
                    stdout += data.toString();
                    console.log('Python stdout:', data.toString());
                });

                pythonProcess.stderr.on('data', (data) => {
                    stderr += data.toString();
                    console.error('Python stderr:', data.toString());
                });

                pythonProcess.on('close', (code) => {
                    console.log(`Python process exited with code: ${code}`);

                    if (code === 0) {
                        try {
                            // Try to parse JSON output
                            const result = JSON.parse(stdout);
                            resolve({
                                success: true,
                                output: result
                            });
                        } catch (parseError) {
                            // If not JSON, return raw output
                            resolve({
                                success: true,
                                output: stdout
                            });
                        }
                    } else {
                        resolve({
                            success: false,
                            error: stderr || `Python process exited with code ${code}`
                        });
                    }
                });

                pythonProcess.on('error', (error) => {
                    console.error('Failed to start Python process:', error);
                    resolve({
                        success: false,
                        error: `Failed to execute Python: ${error.message}`
                    });
                });

            } catch (error) {
                console.error('Error in executePython:', error);
                resolve({
                    success: false,
                    error: error.message
                });
            }
        });
    }

    /**
     * Execute Python code directly (not from a file)
     * @param {string} code - Python code to execute
     * @param {Object} options - Additional spawn options
     * @returns {Promise<{success: boolean, output: string, error?: string}>}
     */
    static async executeCode(code, options = {}) {
        return new Promise(async (resolve, reject) => {
            try {
                const { spawn } = require('child_process');

                // Get the bundled Python executable path
                const pythonExe = await this.getPythonPath();

                console.log('Executing Python code');
                console.log('Python executable:', pythonExe);

                // Spawn the bundled Python process with -c flag
                const pythonProcess = spawn(pythonExe, ['-c', code], {
                    stdio: ['pipe', 'pipe', 'pipe'],
                    ...options
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
                    if (code === 0) {
                        try {
                            const result = JSON.parse(stdout);
                            resolve({
                                success: true,
                                output: result
                            });
                        } catch (parseError) {
                            resolve({
                                success: true,
                                output: stdout
                            });
                        }
                    } else {
                        resolve({
                            success: false,
                            error: stderr || `Python process exited with code ${code}`
                        });
                    }
                });

                pythonProcess.on('error', (error) => {
                    resolve({
                        success: false,
                        error: `Failed to execute Python: ${error.message}`
                    });
                });

            } catch (error) {
                resolve({
                    success: false,
                    error: error.message
                });
            }
        });
    }
}
