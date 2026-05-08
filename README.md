# CardioSeg

CardioSeg is an interactive Python application for cell and nuclei segmentation, manual correction, spatial annotation, and quantitative analysis of H&E-stained cardiac histology images.

The application combines a napari-based desktop GUI with Cellpose segmentation, scikit-image measurements, and optional spatial transcriptomics annotation through Scanpy-compatible datasets.

> **Project status:** Alpha  
> CardioSeg is ready for source distribution and internal evaluation. Segmentation quality should still be validated against your own annotated cardiomyocyte datasets before clinical, diagnostic, or production use.

---

## Features

- Load common microscopy image formats including TIFF, PNG, JPEG, and CZI
- Segment cells and nuclei using configurable Cellpose diameter settings
- Inspect, recolor, and filter label layers in napari
- Manually correct segmentation masks and annotations
- Compute morphology and intensity measurements for labeled cells
- Export measurements to CSV for downstream workflows in Python, R, Fiji/ImageJ, or QuPath
- Attach optional spatial transcriptomics annotations and gene-expression overlays when compatible datasets are available

---

## Installation

CardioSeg is best installed with conda or mamba because of platform-specific dependencies such as Qt, PyTorch, napari, and microscopy image I/O libraries.

### Quick Start

```bash
git clone https://github.com/<username>/cardioseg.git
cd cardioseg

mamba env create -f environment.yml
mamba activate cardioseg

pip install -e .

python -m cardioseg.app
```

If you prefer conda:

```bash
conda env create -f environment.yml
conda activate cardioseg
```

### Windows CUDA Environment

Use the Windows CUDA environment only on systems with a CUDA 11.8-compatible NVIDIA GPU and driver.

```powershell
mamba env create -f environment_win.yml
mamba activate cardioseg-win

pip install -e .

python -m cardioseg.app
```

For Windows systems without CUDA support, use `environment.yml` instead.

---

## Running the Application

Launch CardioSeg with:

```bash
python -m cardioseg.app
```

---

## Basic Workflow

1. Open CardioSeg and load a microscopy image
2. Configure segmentation settings such as approximate object diameter
3. Run automated Cellpose segmentation
4. Inspect cells and nuclei in napari
5. Correct masks or annotations manually if needed
6. Compute measurements and export CSV outputs
7. Optionally load compatible spatial datasets for cell-type or gene-expression overlays

---

## Repository Layout

```text
cardioseg/
├── app.py                         # Application entry point
├── assignment/                    # Cell/nucleus relationship utilities
├── gui/                           # napari widgets, panels, and Qt styling
├── io/                            # Image and measurement I/O helpers
├── measurements/                  # Morphology and texture measurements
├── models/                        # Domain-specific data classes
├── segmentation/                  # Cellpose wrappers and mask utilities
├── spatial/                       # Spatial annotation and ROI helpers
└── utils/                         # Shared utility code

environment.yml                    # CPU/cross-platform conda environment
environment_win.yml                # Windows CUDA environment
pyproject.toml                     # Packaging metadata and console scripts
README.md                          # Project overview and usage
LICENSE                            # License information
```

---

## Outputs

CardioSeg currently focuses on CSV-based measurement export. The project also includes helper modules for image and annotation interoperability with microscopy and spatial-analysis ecosystems.

---

## Publication Readiness Checklist

Before publishing a release, validate the project from a clean checkout:

```bash
python -m pip install -e ".[dev]"
python -m compileall cardioseg
python -m pytest
```

Also validate installation in at least one fresh conda environment created from `environment.yml`. On CUDA-enabled Windows systems, validate `environment_win.yml` as well.

---

## System Requirements

CardioSeg has been tested on:

- macOS (Intel and Apple Silicon)
- Windows 10/11
- Python 3.10
- CUDA 11.8

Performance depends on image size and available GPU memory.

---

## License

This work is licensed under a Creative Commons Attribution 4.0 International License.

[![CC BY 4.0][cc-by-shield]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

---

## Acknowledgements

CardioSeg builds on the scientific Python imaging ecosystem, especially:

- napari
- Cellpose
- scikit-image
- Scanpy
- PyTorch
- tifffile

CardioSeg has been tested on:

- macOS (Apple Silicon)
- Windows 10
- Python 3.10
- CUDA 11.8 

Performance depends on image size and available GPU memory.
