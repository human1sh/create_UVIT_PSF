import os
import sys
import warnings

# 1. Suppress Gaia archive notice & Astropy deprecation warnings BEFORE importing any packages
warnings.filterwarnings("ignore", message=".*current archive release.*")
warnings.filterwarnings("ignore", message=".*TestRunner.*")
warnings.filterwarnings("ignore", message=".*astropy.tests.runner.*")

# 2. Add project root to sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import argparse
from src.config import setup_environment
from src.psf_pipeline import process_uvit_image

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    default_data_dir = os.path.join(base_dir, "data")
    default_output_dir = os.path.join(base_dir, "output")
    parser = argparse.ArgumentParser(description="Automated UVIT PSF Extraction & Interactive Selection Pipeline.")
    parser.add_argument("--filename", type=str, required=True, 
                        help="Path or relative path of the UVIT FITS file (e.g., data/F148W_RD=.../F148W_RD=....fits).")
    parser.add_argument("--data_dir", type=str, default=default_data_dir, help="Directory containing the input FITS image.")
    parser.add_argument("--output_dir", type=str, default=default_output_dir, help="Directory to save selected PSFs.")
    parser.add_argument("--size", type=int, default=81, help="Cutout size in pixels for candidate stars.")
    args = parser.parse_args()

    # Resolution logic for nested subdirectories and absolute/relative paths
    if os.path.exists(args.filename):
        fits_path = os.path.abspath(args.filename)
    else:
        fits_path = os.path.abspath(os.path.join(args.data_dir, args.filename))
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"Could not locate FITS image file at '{fits_path}'")
    setup_environment()
    process_uvit_image(fits_path=fits_path, output_dir=args.output_dir, cutout_size=args.size)


if __name__ == "__main__":
    main()