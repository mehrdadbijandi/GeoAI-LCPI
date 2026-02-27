# Data Description

## Input Datasets for SOM Training

All CSV files contain the following columns:

| Column |                Description               |
|--------|------------------------------------------|
|   ID   | Unique parcel identifier                 |
|   SI   | Normalized Shape Suitability Index [0–1] |
|   BC   | Normalized Boundary Complexity [0–1]     |
|   TF   | Normalized Tenure Fragmentation [0–1]    |
|   AC   | Normalized Accessibility Index [0–1]     |
|  LCPI  | Land Consolidation Potential Index [0–1] |

**Note:** SI, BC, TF, and AC values are identical across all 
five files. The LCPI column differs according to the weighting 
method applied.

|        File       |       Weighting Method       |
|-------------------|------------------------------|
| input_SOM_CR.csv  | CRITIC                       |
| input_SOM_EQ.csv  | Equal Weights                |
| input_SOM_PCA.csv | Principal Component Analysis |
| input_SOM_SH.csv  | Shannon Entropy              |
| input_SOM_STD.csv | Standard Deviation           |


## Data Restrictions
Cadastral parcel data are subject to institutional 
data-sharing restrictions and cannot be publicly distributed. 
The normalized indicator values provided here are sufficient 
to reproduce all SOM training, clustering, and ANOVA results 
reported in the manuscript.
