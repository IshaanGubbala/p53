import torch
import sys
import platform

def verify_environment():
    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")
    print(f"PyTorch: {torch.__version__}")
    
    # Check for MPS (Metal Performance Shaders) on Mac
    if torch.backends.mps.is_available():
        print("✅ MPS (Metal) Acceleration: AVAILABLE")
        device = torch.device("mps")
    elif torch.cuda.is_available():
        print("✅ CUDA Acceleration: AVAILABLE")
        device = torch.device("cuda")
    else:
        print("⚠️  No GPU Acceleration found. Using CPU (slow).")
        device = torch.device("cpu")
        
    print(f"Selected Device: {device}")
    
    # Simple tensor test
    try:
        x = torch.ones(5, device=device)
        y = x * 2
        print(f"Tensor Test: {y.tolist()}")
        print("✅ Environment is ready for Latent Manifold Rescue.")
    except Exception as e:
        print(f"❌ Tensor Test Failed: {e}")
        return False
        
    return True

if __name__ == "__main__":
    if not verify_environment():
        sys.exit(1)
