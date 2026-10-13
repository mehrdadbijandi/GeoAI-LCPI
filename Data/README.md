# Data

Normalized indicators and LCPI scores of the 363 paddy parcels of the three study areas (Mazandaran Province, Iran), as used in the article. No parcel geometry or ownership information is included.

| File | Weighting scheme |
|---|---|
| `plusEQ.csv` | Equal Weights (reference) |
| `plusCR.csv` | CRITIC |
| `plusSTD.csv` | Standard Deviation |
| `plusSH.csv` | Shannon Entropy |
| `plusPCA.csv` | Principal Component Analysis |

## Columns

| Column | Description |
|---|---|
| `ID` | Parcel ID |
| `SI`, `BC`, `TF`, `AC` | Normalized indicators SIʹ, BCʹ, TFʹ, ACʹ |
| `LCPI` | Land Consolidation Potential Index of the scheme |

The normalized indicators are inverted min–max values computed within each study area: 0 = most favorable and 1 = least favorable parcel of the area. The columns `SI`–`AC` are the same in all five files; only `LCPI` differs.

To analyse your own data, provide a CSV file with the same columns.
