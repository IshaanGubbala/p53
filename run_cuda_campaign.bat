@echo off
REM Run fresh campaign with CUDA OpenMM
REM Open this in Command Prompt (not PowerShell)

set CONDA_ROOT=C:\ProgramData\miniconda3
set PROJECT=C:\Users\rajee\OneDrive\Documents\p53

call "%CONDA_ROOT%\Scripts\activate.bat" "%CONDA_ROOT%"
call conda activate openmm-cuda

cd /d "%PROJECT%"

echo.
echo [1/4] Installing p53cad project dependencies...
pip install -e ".[dev,drug]" --quiet
pip install psutil pynvml --quiet
echo Done.
echo.

echo [2/4] Installing CUDA torch (overrides CPU version from setup.py)...
pip install torch --index-url https://download.pytorch.org/whl/cu121 --force-reinstall --quiet
echo Done.
echo.

echo [3/4] Verifying CUDA torch + OpenMM...
python -c "import torch; print('PyTorch:', torch.__version__); print('CUDA torch:', torch.cuda.is_available()); import openmm; p=[openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]; print('OpenMM CUDA:', 'CUDA' in p)"
echo.

echo [4/4] Resuming campaign with CUDA...
echo.
REM Resumes existing run — skips already-completed scenarios automatically
python scripts\run_full_campaign.py --budget medium --run-id campaign_20260217_060021

pause
