@echo off
REM Re-run physics validation (ESMFold + OpenMM energy) on completed campaign
REM *** Run from native cmd.exe or double-click — NOT from PowerShell ***

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
echo [1/3] Ensuring CUDA torch is installed...
python -c "import torch; v=torch.__version__; print('  torch', v, '- OK')" 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo   torch not importable — killing any Python processes holding DLLs...
    taskkill /F /IM python.exe /T 2>nul
    taskkill /F /IM pythonw.exe /T 2>nul
    timeout /t 3 /nobreak >nul
    echo   Installing torch CUDA 12.4...
    pip install torch --index-url https://download.pytorch.org/whl/cu124 --force-reinstall --quiet
)
echo Done.
echo.

echo [2/3] Quick sanity check...
python -c "import os; os.environ['KMP_DUPLICATE_LIB_OK']='TRUE'; import torch; print('  PyTorch:', torch.__version__); print('  CUDA:', torch.cuda.is_available())"
if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ERROR: PyTorch failed to load.
    echo   Try closing all Python/Jupyter windows and running this again.
    pause
    exit /b 1
)
echo.

echo [3/3] Running physics validation...
echo.
python scripts\rerun_physics.py

pause
