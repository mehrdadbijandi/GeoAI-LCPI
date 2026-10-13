# LCPI Toolbox for ArcGIS Pro

A geoprocessing toolbox with two Python script tools that compute the four parcel-level spatial indicators and the Land Consolidation Potential Index (LCPI). All model parameters are exposed as tool parameters, so the indicators can be adapted to other farming systems and cadastral data.

| Script tool | Source file | Purpose |
|---|---|---|
| **IndicesCalculation** | `IndicesCalculation.py` | Computes SI, BC, TF, and AC for each parcel |
| **LCPICalculation** | `LCPICalculation.py` | Normalizes the indicators, derives the weights, and computes the LCPI |

## Requirements

- ArcGIS Pro 3.x (commercial licence) with its bundled Python 3 environment (ArcPy, NumPy, pandas).
- No additional packages.

## Installation

1. Download or clone this repository.
2. In ArcGIS Pro, create a new toolbox (*Catalog* pane → *Toolboxes* → *New Toolbox*) and add two script tools (*New* → *Script*).
3. Point the tools to `IndicesCalculation.py` and `LCPICalculation.py`. Define the parameters listed below in the same order; the scripts read them by position.

## Tool 1 – IndicesCalculation

### Inputs

| Parameter | Type | Description |
|---|---|---|
| Input Parcels | Polygon feature layer (required) | Agricultural parcels in a projected coordinate system (metres) |
| Output Feature Class | Feature class | Copy of the parcels with the indicator fields (default `Score_Parcels`) |
| Owner ID Field | Field (required) | Field identifying the landholder; parcels with the same value form a holding (used by TF) |
| Input Roads | Polyline feature layer (required) | Road centerlines in the same coordinate system |
| Road Width Field | Field (required) | Road width in metres (used by AC) |
| Ideal Ratio | String, default `1:2` | Width-to-length ratio of the reference ("ideal") rectangle for SI and BC |
| Output Ideal Rectangles | Feature class | Ideal rectangle of each parcel, for checking (default `Ideal_Rectangles`) |
| Tolerance | Double, default `2` | Collinearity tolerance in degrees; vertices whose interior angle differs from 180° by less than this value are removed before vertices are counted for BC |
| K Factor | Double, default `500` | Distance-scaling constant K in metres, used in the proximity term K/(K + d) of TF and AC |
| Search Distance (m) | Double, default `5` | Road search distance added to half the road width to form the road buffer for AC |

### Indicators

| Field | Indicator | Definition (paper, Section 2.2) |
|---|---|---|
| `SI_Index` | Shape Suitability | Area of overlap between the parcel and its ideal rectangle, divided by the parcel area. The ideal rectangle has the parcel's area, is centred on its centroid, and is aligned with its minimum bounding rectangle. Range 0–1. |
| `BC_Index` | Boundary Complexity | (P_rect / P_parcel) × 4 / N, where N is the number of vertices after simplification. 1 for a parcel matching the reference rectangle. |
| `TF_Index` | Tenure Fragmentation | (A_R / A_max) × mean of K/(K + d) over all parcels of the holding, including R. R is the largest parcel of the holding. d = 0 for R and for parcels connected to R through shared boundaries (directly or via other parcels of the same holding; corner contact does not count); otherwise d is the centroid distance to R. The same value is assigned to all parcels of a holding. |
| `Dist_To_Ref` | – | Distance d of each parcel to the reference parcel of its holding (−1 for parcels without an owner ID) |
| `AC_Index` | Accessibility | Parcels with direct road access: 2·Σ(Lᵢ·Bᵢ) / (P·B_max), where Lᵢ is the boundary length within the road buffer and Bᵢ the adjacent road width. Parcels without direct access: ε·K/(K + D_road), with ε the smallest AC value among parcels with direct access. |

Higher values indicate more favorable conditions for all four indicators. SI lies in [0, 1]. BC equals 1 for a parcel that matches the reference rectangle and can slightly exceed 1 for more compact parcels (≈ 1.06 for a square). AC lies in [0, 2]. B_max is the largest value of the road-width field. The tool also adds the field `NEAR_DIST` (distance to the nearest road).

## Tool 2 – LCPICalculation

### Inputs

| Parameter | Type | Description |
|---|---|---|
| Input Feature Class | Polygon feature class | Output of IndicesCalculation |
| Output Feature Class | Feature class | Copy of the input with the normalized indicators, the weights, and the LCPI |
| SI Field (Shape), BC Field (Boundary), TF Field (Tenure), AC Field (Access) | Field (required) | The four indicator fields |
| Weighting Method | String (value list) | `Equal Weights` (used in the paper), `CRITIC`, `Shannon Entropy`, `Standard Deviation`, or `PCA - 1st Component` |

### Output fields

| Field | Content |
|---|---|
| `SI`, `BC`, `TF`, `AC` | Normalized indicators SIʹ–ACʹ (inverted min–max) |
| `W_SI_Cal`, `W_BC_Cal`, `W_TF_Cal`, `W_AC_Cal` | Weights of the selected method |
| `LCPI_EQ`, `LCPI_CR`, `LCPI_SH`, `LCPI_STD`, or `LCPI_PCA` | LCPI of the selected method |

### Method

1. Each indicator is converted to an inefficiency score by inverted min–max normalization: X′ = (max − x) / (max − min). After this step, 1 is the least favorable parcel and 0 the most favorable.
2. Weights are derived from the normalized indicators with the selected method. Standard deviations are sample standard deviations; null values are replaced by the mean of the field.
   - Equal Weights: wⱼ = 0.25.
   - Standard Deviation: proportional to σⱼ.
   - CRITIC: σⱼ·Σₖ(1 − rⱼₖ), normalized.
   - Shannon Entropy: proportional to 1 − Eⱼ, with pᵢⱼ = X′ᵢⱼ / Σᵢ X′ᵢⱼ and Eⱼ = −(1/ln n) Σᵢ pᵢⱼ ln pᵢⱼ.
   - PCA: absolute loadings of the first principal component, normalized.
3. LCPI = Σⱼ wⱼ X′ⱼ. Higher values indicate a greater need for consolidation.

Normalization and weights are computed over the parcels of the input layer. Run the tool separately for each study area if priorities are to be relative to the area, as in the paper. To compare areas, normalize against a common reference range.

## Notes

- All lengths and areas are computed in the units of the input coordinate system. Use a projected coordinate system in metres.
- The parameter values used in the paper are the defaults listed above.
