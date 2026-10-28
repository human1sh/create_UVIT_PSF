# UVIT PSF Extractor (`create_UVIT_PSF`)

An automated Python pipeline built to extract, analyze, and select Point Spread Functions (PSFs) from **UVIT (Ultra Violet Imaging Telescope)** FUV and NUV FITS images. 

The tool queries the **Gaia DR3** catalog for field star coordinates, extracts 2D cutouts, performs 1D isophotal fitting using a combined **Moffat + Exponential profile**, and provides an interactive visual environment for reviewing and selecting high-quality PSF candidates.

---

## Key Features

- **Dynamic UVIT Filter Detection**: Automatically identifies whether an image is **FUV** (e.g., `F148W`, `F154W`, `F169M`) or **NUV** (e.g., `N242W`, `N263M`) from FITS headers, filenames, or parent folder structures.
- **Gaia DR3 Cross-Matching**: Automatically retrieves stellar coordinates within the UVIT Field of View using `astroquery.gaia`.
- **Background & Star Detection**: Uses `IRAFStarFinder` combined with local background estimation to locate non-saturated candidate stars.
- **1D Isophotal Profile Analysis**: Fits candidate profiles using a two-stage model:
  $$\text{Model}(r) = I_{m0} \left(1 + \left(\frac{r}{\alpha}\right)^2\right)^{-\beta} + I_{e0} \exp\left(-\frac{r}{h_0}\right)$$
- **Interactive Candidate Selection UI**: Displays 2-panel diagnostic plots for each candidate star (Flux fit with embedded image thumbnail + residual plot) allowing direct user approval before saving outputs.
- **Native WCSAxes Overlay**: Renders field coverage maps without external legacy dependencies like `aplpy`.

---

## Repository Structure

```text
create_UVIT_PSF/
│
├── data/                       # Input UVIT FITS images (ignored by git)
│   └── F148W_RD=.../           # Subdirectory containing target FITS image
│
├── output/                     # Saved output FITS cutouts & PNG plots (ignored by git)
│
├── src/
│   ├── __init__.py
│   ├── config.py               # Environment setups & Matplotlib styling
│   └── uvit_psf_pipeline.py    # Main PSF extraction & fitting pipeline
│
├── main.py                     # Command-line entry point
├── .gitignore                  # Git ignore definitions
├── requirements.txt            # Python dependencies
└── README.md                   # Project documentation
```

---

## Installation & Requirements

### 1. Clone Repository
```bash
git clone https://github.com/YOUR_USERNAME/PSF.git
cd create_UVIT_PSF
```

### 2. Install Dependencies
It is recommended to use a virtual environment or Conda environment with Python $\ge$ 3.9.

```bash
pip install -r requirements.txt
```

#### `requirements.txt`
```text
numpy
pandas
matplotlib
astropy
astroquery
lmfit
photutils
```

---

## Quick Start / Usage

Place your input UVIT FITS files into the `data/` directory.

### Running the Pipeline

You can pass relative nested paths, direct file paths, or file names inside `data/`:

```bash
# Example 1: Passing nested subdirectory file path
python main.py --filename data/F148W_RD=217.532,63.855/F148W_RD=217.532,63.855.fits

# Example 2: Customizing cutout size (default: 81 pixels)
python main.py --filename data/F148W_RD=217.532,63.855/F148W_RD=217.532,63.855.fits --size 101
```

### Options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--filename` | `str` | *(Required)* | Path or filename of the target UVIT FITS file. |
| `--data_dir` | `str` | `data/` | Directory containing input data folders. |
| `--output_dir` | `str` | `output/` | Directory where accepted PSFs and plots are saved. |
| `--size` | `int` | `81` | Cutout side length in pixels for candidate stars. |

---

## Pipeline Workflow

1. **Gaia Cone Search**: Queries Gaia DR3 catalog over a 14 arcmin radius centered on the image WCS.
2. **Detection**: Runs `IRAFStarFinder` to isolate candidate stars within 2000 pixels from the center.
3. **Overview Map**: Plots the field with detected candidates on a native `WCSAxes` grid.
4. **Isophotal Profile Fitting**:
   - Calculates radial surface brightness profiles via `photutils.isophote.Ellipse`.
   - Fits outer tail ($r > 6''$) with an Exponential background profile.
   - Fits full profile with a Moffat + Exponential model to extract FWHM and shape parameters ($\alpha, \beta$).
5. **Interactive Review**:
   - Displays candidate cutout and 2-panel diagnostic plot.
   - Prompts in terminal: `Keep Candidate #i? [y/N]:`.
6. **Output Generation**:
   - `Star_i-PSF-FIT.fits`: FITS image cutout of the accepted star.
   - `Star_i-PSF-PAR.fits`: FITS table containing fitted parameters ($I_{m0}, \alpha, \beta, I_{e0}, h_0, \text{FWHM}$).
   - `Star_i-PSF-1DP.png`: Publication-ready 2-panel 1D profile figure.

---

## Output Plot Format

For every accepted candidate, the code saves a high-resolution 2-panel figure formatted as follows:
- **Top Panel**: Radial flux profile (log scale) showing the raw PSF data, total best-fit model, Moffat component, Exponential tail component, sky background RMS, and an inset thumbnail of the star image.
- **Bottom Panel**: Fractional residual profile: $\frac{\text{Flux}_{\text{data}} - \text{Flux}_{\text{model}}}{\text{Flux}_{\text{data}}}$

---

## License

This project is open-source and available under the [MIT License](LICENSE).
