# CardioSeg

CardioSeg is an interactive Python application for cell and nuclei segmentation,
manual correction, spatial annotation, and quantitative analysis of H&E-stained cardiac histology images.

The application combines a napari-based desktop GUI with Cellpose segmentation,
scikit-image measurements, and optional spatial transcriptomics annotation via
Scanpy-loaded datasets.

> **Project status:** alpha. The package is ready for source distribution and
> internal evaluation, but segmentation quality should be validated against your
> own annotated cardiomyocyte datasets before clinical, diagnostic, or production
> use.

## Features

- Load common microscopy image formats, including TIFF, PNG, JPEG, and CZI.
- Segment cells and nuclei with configurable Cellpose diameter settings.
- Inspect, recolor, and filter label layers in napari.
- Manually correct segmentation results with napari layers.
- Compute morphology and intensity measurements for labeled cells.
- Export measurements to CSV for downstream Python, R, Fiji/ImageJ, or QuPath
  workflows.
- Attach optional spatial transcriptomics annotations and gene-expression based
  coloring when compatible spatial data are available.

## Installation

CardioSeg is recommended to be installed using conda or mamba due to platform-specific dependencies including Qt, PyTorch, napari, and image I/O libraries.

### Quick start

```bash
git clone https://github.com/<username>/cardioseg.git
cd cardioseg

mamba env create -f environment.yml
mamba activate cardioseg

pip install -e .

python -m cardioseg.app
```

If you use conda instead of mamba:

```bash
conda env create -f environment.yml
conda activate cardioseg
```

### Windows CUDA environment

Use the Windows environment only on machines with a CUDA 11.8-compatible NVIDIA GPU and driver:

```powershell
mamba env create -f environment_win.yml
mamba activate cardioseg-win

pip install -e .

python -m cardioseg.app
```

For Windows machines without CUDA, use `environment.yml` instead.

````md
## Running the application

After installation, CardioSeg can be launched using:

```bash
python -m cardioseg.app

## Basic workflow

1. Open CardioSeg and load a microscopy image.
2. Choose segmentation settings, including approximate object diameter.
3. Run automated Cellpose segmentation.
4. Inspect cells and nuclei in napari.
5. Correct masks or annotations as needed.
6. Compute measurements and export CSV outputs.
7. Optionally load compatible spatial data for cell-type or gene overlays.

## Repository layout

```text
cardioseg/
├── app.py                         # Application entry point
├── assignment/                    # Cell/nucleus relationship utilities
├── gui/                           # napari widgets, panels, and Qt stylesheet
├── io/                            # Measurement and image I/O helpers
├── measurements/                  # Morphology and texture measurements
├── models/                        # Domain data classes
├── segmentation/                  # Cellpose wrappers and mask utilities
├── spatial/                       # Spatial annotation and ROI helpers
└── utils/                         # Shared utility code

environment.yml                    # CPU/cross-platform conda environment
environment_win.yml                # Windows CUDA conda environment
pyproject.toml                     # Packaging metadata and console script
README.md                          # Project overview and usage
LICENSE                            # CC-BY license
```

## Outputs

CardioSeg currently focuses on CSV measurement export. The package also contains
image and annotation helper modules intended to support interoperable workflows
with microscopy and spatial-analysis tools as the project evolves.

## Publication readiness checklist

Before publishing a release, verify the following from a clean checkout:

```bash
python -m pip install -e ".[dev]"
python -m compileall cardioseg
python -m pytest
```

For environment validation, create at least one fresh conda environment from
`environment.yml`; on CUDA-enabled Windows machines, also validate
`environment_win.yml`.

````md
## System requirements

CardioSeg has been tested on:

- macOS (Intel and Apple Silicon)
- Windows 10/11
- Python 3.10
- CUDA 11.8 

Performance depends on image size and available GPU memory.

## License

```md
[![CC BY 4.0][cc-by-shield]][cc-by]
```

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Acknowledgements

CardioSeg builds on the scientific Python imaging ecosystem, especially napari,
Cellpose, scikit-image, Scanpy, PyTorch, and tifffile.
