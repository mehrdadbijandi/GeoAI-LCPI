# GeoAI-LCPI

**Land Consolidation Potential Index — Custom Geoprocessing 
Toolbox and GeoAI-based Validation for Agricultural Parcel 
Prioritization**

## Overview

This repository contains the source code and datasets 
accompanying the following publication:

> Bijandi,M. (2026). [Title]. *[Journal Name]*, 
> [Volume], [Pages]. DOI: []
If you use this code or dataset in your research, please cite the above publication.

## Components

### 1. Custom Geoprocessing Toolbox (ArcGIS Pro)
- **spatial_indices_module.py**: Computes four spatial 
  indicators (Shape Index, Boundary Complexity, Tenure 
  Fragmentation, Accessibility) for agricultural parcels
- **weighting_LCPI_module.py**: Calculates LCPI using five 
  weighting methods (CRITIC, Shannon Entropy, STD, PCA, 
  Equal Weights)

### 2. GeoAI Analysis (SOM)
- Self-Organizing Map training and clustering scripts
- Cluster validation and ANOVA analysis

### 3. Sample Data
- Normalized spatial indicator values for 363 agricultural 
  parcels across three study areas
- Five CSV files corresponding to five weighting scenarios

## Requirements

|  Component   |                    Requirements                     |
|--------------|-----------------------------------------------------|
| Toolbox      | ArcGIS Pro 3.x, Python 3.x, arcpy                   |
| SOM Analysis | Python 3.x, MiniSom 2.3.1, NumPy, SciPy, Matplotlib |

## License

This project is licensed under the MIT License 
see the LICENSE file for details.

## Citation

If you use this toolbox or dataset in your research, please 
cite:

```bibtex
@article{Bijandi2026lcpi,
  title={[Title]},
  author={[Authors]},
  journal={[Journal]},
  year={2026},
  doi={[DOI]}
}
Contact
Mehrdad Bijandi — bijandi@tafreshu.ac.ir
