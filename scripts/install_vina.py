#!/usr/bin/env python3
"""
Download and install AutoDock Vina CLI for Windows.

This script downloads the pre-compiled AutoDock Vina binary and places it
in a location on your PATH.
"""

import os
import sys
import urllib.request
import zipfile
from pathlib import Path

VINA_VERSION = "1.2.7"
VINA_URL = f"https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v{VINA_VERSION}/vina_{VINA_VERSION}_windows_x86_64.zip"


def main():
    print(f"AutoDock Vina Installer v{VINA_VERSION}")
    print("=" * 60)

    # Get Python Scripts directory (on PATH)
    scripts_dir = Path(sys.executable).parent / "Scripts"
    if not scripts_dir.exists():
        print(f"ERROR: Scripts directory not found: {scripts_dir}")
        sys.exit(1)

    print(f"Installation directory: {scripts_dir}")
    print()

    # Download the zip file
    zip_path = Path.cwd() / "vina_temp.zip"
    print(f"Downloading Vina from {VINA_URL}...")

    try:
        urllib.request.urlretrieve(VINA_URL, zip_path)
        print(f"Downloaded to {zip_path}")
    except Exception as e:
        print(f"ERROR: Failed to download Vina: {e}")
        print()
        print("Manual installation:")
        print(f"1. Download from: {VINA_URL}")
        print(f"2. Extract vina.exe")
        print(f"3. Copy vina.exe to: {scripts_dir}")
        sys.exit(1)

    # Extract vina.exe
    print("Extracting vina.exe...")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Find vina.exe in the zip
            vina_exe = None
            for name in zip_ref.namelist():
                if name.endswith('vina.exe'):
                    vina_exe = name
                    break

            if not vina_exe:
                print("ERROR: vina.exe not found in the downloaded archive")
                sys.exit(1)

            # Extract to Scripts directory
            extracted_data = zip_ref.read(vina_exe)
            target_path = scripts_dir / "vina.exe"

            with open(target_path, 'wb') as f:
                f.write(extracted_data)

            print(f"Installed vina.exe to {target_path}")

    except Exception as e:
        print(f"ERROR: Failed to extract: {e}")
        sys.exit(1)

    finally:
        # Clean up
        if zip_path.exists():
            zip_path.unlink()
            print("Cleaned up temporary files")

    print()
    print("=" * 60)
    print("AutoDock Vina installation complete!")
    print()
    print("Verify installation by running:")
    print("  vina --help")
    print()


if __name__ == "__main__":
    main()
