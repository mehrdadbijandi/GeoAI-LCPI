# ArcGIS Pro Custom Geoprocessing Toolbox

## Module 1: spatial_indices_module.py
Automated computation of four spatial indicators (SI, BC, TF, AC)
for agricultural parcels.
- Parcel Shape Suitability Index (SI)
- Boundary Complexity Index (BC)
- Tenure Fragmentation Index (TF)
- Accessibility Index (AC)
### Requied Inputs:  
- The polygon feature layer representing agricultural parcels to be analyzed.
- The line feature layer representing the road network for accessibility analysis.
- The field containing unique identifiers for parcel owners (used to calculate Tenure Fragmentation).
- The numeric field in the roads layer that specifies the width of each road segment.
### Output:
- Results are appended as new attribute fields to the input feature class.

## Module 2: weighting_LCPI_module.py  
LCPI calculation using five objective weighting methods 
(CRITIC, Shannon Entropy, STD, PCA, Equal Weights).
### Requied Input:  
- The processed layer from spatial_indices_module.py (containing SI, BC, TF, and AC fields).
### Output:
- The Land Consolidation Potential Index (LCPI) for any parcels.

## Requirements
- ArcGIS Pro 3.x
- Python 3.x (with arcpy)

## Usage
Parameters and configuration are described in the manuscript.
