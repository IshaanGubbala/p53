@echo off
REM Re-run physics validation (ESMFold + OpenMM energy) on completed campaign
REM *** Must be run by double-clicking in Windows Explorer OR from native cmd.exe ***
REM *** Do NOT run via PowerShell ./rerun_physics.bat — DLL paths break ***

set CONDA_ROOT=C:\ProgramData\miniconda3
set PROJECT=C:\Users\rajee\OneDrive\Documents\p53

REM Fix OMP conflict between Intel (OpenMM) and LLVM (PyTorch) OpenMP runtimes
set KMP_DUPLICATE_LIB_OK=TRUE
set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%
set MKL_NUM_THREADS=%NUMBER_OF_PROCESSORS%

call "%CONDA_ROOT%\Scripts\activate.bat" "%CONDA_ROOT%"
call conda activate openmm-cuda

cd /d "%PROJECT%"

echo.
echo [1/3] Verifying CUDA torch install...
python -c "import torch; v=torch.__version__; assert 'cu' in v, f'CPU-only: {v}'; print('  torch', v, '- OK')" 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo   CUDA torch not found, installing...
    pip install torch --index-url https://download.pytorch.org/whl/cu121 --quiet
)
echo Done.
echo.

echo [2/3] Quick sanity check...
python -c "import os; os.environ['KMP_DUPLICATE_LIB_OK']='TRUE'; import torch; print('  PyTorch:', torch.__version__); print('  CUDA:', torch.cuda.is_available())"
if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ERROR: PyTorch failed to load. Try running this file by double-clicking
    echo        it in Windows Explorer instead of from PowerShell/cmd.
    pause
    exit /b 1
)
echo.

echo [3/3] Running physics validation...
echo.
python scripts\rerun_physics.py

pause
