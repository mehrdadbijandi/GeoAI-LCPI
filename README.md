# GeoAI-LCPI

**Land Consolidation Potential Index (LCPI): a configurable GIS toolbox for parcel-level land consolidation prioritization, with GeoAI-based diagnostic profiling**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- After the first release, add the Zenodo DOI badge here -->

## Overview

This repository contains the source code and data accompanying the article:

> Bijandi M.(2026). *A configurable GIS tool for parcel-level land consolidation prioritization with GeoAI-based diagnostic profiling.* .

The LCPI ranks agricultural parcels by their spatial-structural need for land consolidation. It combines four parcel-level indicators:

| Indicator | What it measures |
|---|---|
| **SI** – Shape Suitability | overlap of the parcel with an ideal rectangle of equal area and a user-defined aspect ratio (new) |
| **BC** – Boundary Complexity | perimeter and number of vertices relative to the same reference rectangle (new) |
| **TF** – Tenure Fragmentation | spatial concentration of the parcels of a holding around its largest parcel |
| **AC** – Accessibility | road-adjacent frontage weighted by road width, or distance to the nearest road |

The indicators are normalized within each study area (inverted min–max, so that 1 = least favorable) and combined as a weighted sum. Five weighting schemes are available: Equal Weights, CRITIC, Standard Deviation, Shannon Entropy, and PCA.

The LCPI shows **where** consolidation is most needed. A Self-Organizing Map (SOM) of the four normalized indicators adds a diagnostic profile, showing **which deficiencies** drive that need.

If you use this code or the data, please cite the article above.

## Repository structure

```
GeoAI-LCPI/
├── toolbox/          ArcGIS Pro geoprocessing toolbox (indicators and LCPI)
│   ├── spatial_indices_module.py    script tool "IndicesCalculation"
│   ├── weighting_LCPI_module.py     script tool "LCPICalculation"
│   └── README.md
├── som/              GeoAI analysis with a Self-Organizing Map
│   ├── som_lcpi.py               SOM training, Ward clustering, figures
│   └── README.md
├── data/             normalized indicators and LCPI of the 363 study parcels
│   └── README.md
├── requirements.txt  Python packages for the SOM analysis
├── LICENSE
```

## Workflow

1. **Indicators.** In ArcGIS Pro, run *IndicesCalculation* on the parcel layer (parcel boundaries with an owner-ID field) and the road layer (centerlines with a width field). It adds the fields `SI_Index`, `BC_Index`, `TF_Index`, and `AC_Index`.
2. **LCPI.** Run *LCPICalculation* to normalize the indicators and compute the LCPI with the selected weighting scheme.
3. **Diagnostic profiling.** Export the normalized indicators to CSV (format in `data/README.md`) and run the scripts in `som/`.

See the README of each folder for parameters, inputs, and outputs.

## Requirements

| Component | Requirements |
|---|---|
| Toolbox | ArcGIS Pro 3.x (commercial licence) with its bundled Python 3 environment (ArcPy, NumPy, pandas) |
| SOM analysis | Python ≥ 3.9 (tested with 3.11), NumPy, pandas, SciPy, scikit-learn, Matplotlib; MiniSom is included (`som/minisom.py`) |

## Data availability

The cadastral parcel boundaries and ownership data of the study areas are subject to institutional restrictions and are not included. The folder `data/` contains the normalized indicator values and LCPI scores of the 363 parcels, without geometry or ownership information. These are sufficient to reproduce the SOM analysis.

## License

The code in this repository is released under the MIT License (see `LICENSE`). `som/minisom.py` is third-party code distributed under the licence of the MiniSom project.

## Contact

Mehrdad Bijandi, Assistant Professor, Department of Geodesy and Surveying Engineering, Tafresh University, bijandi@tafreshu.ac.ir
