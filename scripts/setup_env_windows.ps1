# PowerShell setup script for FEM Contact Mechanics on Windows

Write-Host "=== FEM Contact Mechanics Environment Setup (Windows) ===" -ForegroundColor Green

# Check Python version
function Check-PythonVersion {
    $pythonCmd = Get-Command python -ErrorAction SilentlyContinue
    if (-not $pythonCmd) {
        Write-Host "Error: Python not found. Please install Python 3.10 or later." -ForegroundColor Red
        exit 1
    }
    
    $version = python --version 2>&1
    if ($version -match "Python (\d+)\.(\d+)") {
        $major = [int]$matches[1]
        $minor = [int]$matches[2]
        
        if ($major -ne 3 -or $minor -lt 10) {
            Write-Host "Error: Python 3.10 or later is required. Found: $version" -ForegroundColor Red
            exit 1
        }
        
        Write-Host "Found $version" -ForegroundColor Green
    }
}

# Parse command line arguments
$environment = "base"
if ($args.Count -gt 0) {
    switch ($args[0].ToLower()) {
        "ml" { $environment = "ml" }
        "viz" { $environment = "viz" }
        "dev" { $environment = "dev" }
        "all" { $environment = "all" }
        "base" { $environment = "base" }
        default {
            Write-Host "Usage: .\setup_env_windows.ps1 [base|ml|viz|dev|all]"
            Write-Host "  base: Core dependencies only (default)"
            Write-Host "  ml:   Base + Machine Learning dependencies"
            Write-Host "  viz:  Base + Visualization dependencies"
            Write-Host "  dev:  All dependencies + development tools"
            Write-Host "  all:  All dependencies"
            exit 1
        }
    }
}

# Check prerequisites
Check-PythonVersion

# Create virtual environment
if (Test-Path "venv") {
    Write-Host "Virtual environment already exists. Removing old environment..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force venv
}

Write-Host "Creating virtual environment..." -ForegroundColor Green
python -m venv venv

# Activate virtual environment
Write-Host "Activating virtual environment..." -ForegroundColor Green
& ".\venv\Scripts\Activate.ps1"

# Upgrade pip, setuptools, and wheel
Write-Host "Upgrading pip, setuptools, and wheel..." -ForegroundColor Green
python -m pip install --upgrade pip setuptools wheel

# Install build dependencies
Write-Host "Installing build dependencies..." -ForegroundColor Green
pip install "pybind11>=2.13.0,<2.14.0" "numpy>=1.26.0,<1.27.0"

# Install requirements based on selected environment
Write-Host "Installing $environment dependencies..." -ForegroundColor Green

switch ($environment) {
    "base" { pip install -r requirements\base.txt }
    "ml" { pip install -r requirements\ml.txt }
    "viz" { pip install -r requirements\viz.txt }
    "dev" { pip install -r requirements\dev.txt }
    "all" {
        pip install -r requirements\ml.txt
        pip install -r requirements\viz.txt
        pip install -r requirements\dev.txt
    }
}

# Build and install the package
Write-Host "Building C++ extensions..." -ForegroundColor Green
pip install -e .

# Verify installation
Write-Host "Verifying installation..." -ForegroundColor Green
try {
    python -c "import PyClasses.eigen_backend; print('C++ extensions loaded successfully!')"
} catch {
    Write-Host "Warning: C++ extensions may not have built correctly." -ForegroundColor Yellow
    Write-Host "Please ensure you have Visual Studio Build Tools and Eigen3 installed."
}

# Create activation script
@'
# Quick activation script for FEM environment
& ".\venv\Scripts\Activate.ps1"
Write-Host "FEM Contact Mechanics environment activated!" -ForegroundColor Green
Write-Host "Python: $(Get-Command python | Select-Object -ExpandProperty Source)"
Write-Host "Version: $(python --version)"
'@ | Out-File -FilePath "activate_env.ps1"

Write-Host "=== Setup Complete! ===" -ForegroundColor Green
Write-Host ""
Write-Host "To activate the environment, run:"
Write-Host "  .\venv\Scripts\Activate.ps1"
Write-Host "Or use the convenience script:"
Write-Host "  .\activate_env.ps1"
Write-Host ""
Write-Host "To deactivate, run:"
Write-Host "  deactivate"
Write-Host ""
Write-Host "Environment type: $environment" -ForegroundColor Green

# Run a simple test
Write-Host "`nRunning basic import test..." -ForegroundColor Yellow
python -c @"
import numpy as np
import scipy
print(f'NumPy version: {np.__version__}')
print(f'SciPy version: {scipy.__version__}')
try:
    import tensorflow as tf
    print(f'TensorFlow version: {tf.__version__}')
except ImportError:
    print('TensorFlow not installed (not in $environment environment)')
try:
    import matplotlib
    print(f'Matplotlib version: {matplotlib.__version__}')
except ImportError:
    print('Matplotlib not installed (not in $environment environment)')
"@

Write-Host "`nSetup completed successfully!" -ForegroundColor Green