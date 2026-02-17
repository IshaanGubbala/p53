@echo off
REM Setup CUDA OpenMM - Fixed version with TOS acceptance

echo ========================================
echo OpenMM CUDA Setup (Fixed)
echo ========================================
echo.

REM Initialize conda
set CONDA_ROOT=C:\ProgramData\miniconda3
set PATH=%CONDA_ROOT%;%CONDA_ROOT%\Scripts;%CONDA_ROOT%\Library\bin;%PATH%

REM Initialize conda for this shell
call "%CONDA_ROOT%\Scripts\activate.bat" "%CONDA_ROOT%"

echo [1/6] Accepting conda Terms of Service...
call conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
call conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
call conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2
echo.

echo [2/6] Creating environment 'openmm-cuda' with Python 3.12...
call conda create -n openmm-cuda python=3.12 -y
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Failed to create conda environment
    pause
    exit /b 1
)
echo.

echo [3/6] Activating openmm-cuda environment...
call conda activate openmm-cuda
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Failed to activate environment
    pause
    exit /b 1
)
echo.

echo [4/6] Installing CUDA-enabled OpenMM from conda-forge...
echo This may take a few minutes...
call conda install -c conda-forge openmm cudatoolkit=11.8 pdbfixer mdtraj -y
echo.

echo [5/6] Installing Python dependencies...
pip install numpy pandas torch transformers
echo.

echo [6/6] Verifying OpenMM CUDA platform...
python -c "import openmm; platforms = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]; print('\\nAvailable platforms:', platforms); cuda_ok = 'CUDA' in platforms; print('CUDA available:', cuda_ok); print('\\n' + ('SUCCESS: CUDA is ready!' if cuda_ok else 'WARNING: CUDA not found'))"
echo.

echo ========================================
echo Setup Complete!
echo ========================================
echo.
echo Environment: openmm-cuda
echo.
echo To use CUDA OpenMM:
echo   1. conda activate openmm-cuda
echo   2. cd C:\Users\rajee\OneDrive\Documents\p53
echo   3. python resume_campaign.py --run-id [your-run-id]
echo.
echo Current campaign continues with CPU (don't interrupt it)
echo CUDA will be used automatically in future runs!
echo.

pause
