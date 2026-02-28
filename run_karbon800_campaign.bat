@echo off
REM Run campaign on OnLogic Karbon 800 (Intel 12th Gen / Xe GPU via IPEX + OpenCL OpenMM)
REM Open in Command Prompt (not PowerShell)
REM
REM Hardware: Karbon 800 with Intel Core i5/i7/i9 (Alder Lake) + Intel Xe integrated GPU
REM PyTorch: Intel Extension for PyTorch (IPEX) with XPU backend for Xe GPU
REM OpenMM:  OpenCL backend for energy minimization / MD on Intel Xe GPU
REM
REM If you later add a discrete NVIDIA GPU via PCIe, switch the conda env to
REM "openmm-cuda" and run run_cuda_campaign.bat instead.

set CONDA_ROOT=C:\ProgramData\miniconda3
set PROJECT=C:\Users\rajee\OneDrive\Documents\p53

call "%CONDA_ROOT%\Scripts\activate.bat" "%CONDA_ROOT%"
call conda activate base

cd /d "%PROJECT%"

echo.
echo [1/5] Installing p53cad project dependencies...
pip install -e ".[dev,drug]" --quiet
pip install psutil --quiet
echo Done.
echo.

echo [2/5] Installing Intel Extension for PyTorch (IPEX) for Xe GPU acceleration...
REM IPEX XPU wheel targets Intel Arc / Xe graphics via oneAPI
REM See: https://intel.github.io/intel-extension-for-pytorch/
pip install torch==2.1.0+cpu --index-url https://download.pytorch.org/whl/cpu --quiet
pip install intel_extension_for_pytorch --extra-index-url https://pytorch-extension.intel.com/release-whl/stable/xpu/us/ --quiet
echo Done.
echo.

echo [3/5] Installing OpenMM with OpenCL support for Xe GPU physics...
REM conda-forge OpenMM includes OpenCL plugin which works with Intel Xe integrated graphics
conda install -c conda-forge openmm --yes --quiet
echo Done.
echo.

echo [4/5] Verifying IPEX XPU + OpenMM OpenCL...
python -c "import torch; print('PyTorch:', torch.__version__)"
python -c "import intel_extension_for_pytorch as ipex; import torch; print('IPEX XPU available:', torch.xpu.is_available())" 2>nul || echo "IPEX XPU not detected - will run on CPU"
python -c "import openmm; platforms=[openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]; print('OpenMM platforms:', platforms); print('OpenCL available:', 'OpenCL' in platforms)"
echo.

echo [5/5] Starting campaign (auto-selects XPU if available, else CPU)...
echo.
REM --device auto will pick: CUDA > MPS > XPU > CPU
python scripts\run_full_campaign.py --budget medium

pause
