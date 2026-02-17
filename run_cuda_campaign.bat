@echo off
REM Run fresh campaign with CUDA OpenMM
REM Open this in Command Prompt (not PowerShell)

set CONDA_ROOT=C:\ProgramData\miniconda3
set PROJECT=C:\Users\rajee\OneDrive\Documents\p53

call "%CONDA_ROOT%\Scripts\activate.bat" "%CONDA_ROOT%"
call conda activate openmm-cuda

cd /d "%PROJECT%"

echo.
echo [1/3] Installing p53cad project dependencies...
pip install -e ".[dev,drug]" --quiet
pip install psutil pynvml --quiet
echo Done.
echo.

echo [2/3] Verifying CUDA OpenMM...
python -c "import openmm; p=[openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]; print('Platforms:',p); print('CUDA ready:', 'CUDA' in p)"
echo.

echo [3/3] Starting fresh campaign with CUDA...
echo.
python scripts\run_full_campaign.py --budget medium

pause
