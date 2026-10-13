# GeoAI diagnostic profiling with a Self-Organizing Map (SOM)

These scripts group parcels by their profiles of the four normalized indicators (SIʹ, BCʹ, TFʹ, ACʹ) and describe which deficiencies drive their LCPI. The SOM uses no indicator weights. The resulting clusters are used:

- as diagnostic profiles of the parcels (paper, Section 3.4);
- as an internal benchmark for the weighting schemes (paper, Section 3.5).

## Scripts

| Script | Purpose | Paper |
|---|---|---|
| `som_lcpi.py` | Trains the reference SOM, clusters the prototypes with Ward's method into four clusters, and reports quantization and topographic errors and figures | Sections 2.5 and 3.4; Figs. 6–8 |
| `cluster_validity.py` | Retrains the same SOM and computes the within-cluster sum of squares, silhouette score, and Davies–Bouldin index for 2–17 clusters | Section 3.4.2; Fig. A.1 |
| `cluster_separation.py` | Separation of the SOM clusters by the LCPI of each weighting scheme: one-way ANOVA (F, R²), Welch's ANOVA, Kruskal–Wallis | Section 3.5; Table A.2 |
| `minisom.py` | Copy of the MiniSom library used in the study (third-party) | – |

## Reference configuration

| Setting | Value |
|---|---|
| Input | `../data/plusEQ.csv` (SIʹ, BCʹ, TFʹ, ACʹ of 363 parcels) |
| Pre-processing | Yeo–Johnson transformation and standardization |
| Lattice | 9 × 13 neurons, hexagonal topology |
| Initialization | first two principal components |
| Neighbourhood | Gaussian; radius 4.0 → 1.0 |
| Learning rate | 0.5 → 0.01, exponential decay |
| Iterations | 100,000 (random sample order) |
| Seeds | random seed 42; training order (shuffle seed) 71 |
| Clustering | Ward agglomerative clustering of the prototypes, 4 clusters |

## Usage

```bash
pip install -r ../requirements.txt
python som_lcpi.py [input_csv] [shuffle_seed]     # defaults: ../data/plusEQ.csv 71
python cluster_validity.py [input_csv] [shuffle_seed]
python cluster_separation.py                      # after som_lcpi.py
```

All outputs are written to `som/outputs/`:

| File | Content |
|---|---|
| `Enhanced_Dataset.csv` | input data with the column `SOM_Cluster` |
| `SOM_Complete_Analysis_Hexagonal.png` | component planes, U-matrix, correlations of the neuron weights |
| `SOM_Clustering_Results_Hexagonal.png` | LCPI plane, cluster map, LCPI by cluster, cluster profiles |
| `cluster_validity.csv`, `cluster_validity.png` | validity indices for 2–17 clusters |
| `cluster_separation.csv` | ANOVA, Welch and Kruskal–Wallis statistics per weighting scheme |

## Expected results with the data provided

| Quantity | Value |
|---|---|
| Quantization error | 0.6285 |
| Topographic error | 0.0193 |
| Cluster sizes (clusters 0–3) | 80, 152, 79, 52 |
| Mean LCPI by cluster | 0.573, 0.615, 0.679, 0.821 |
| Silhouette score (4 clusters) | 0.361, the highest for 2–17 clusters |
| Davies–Bouldin index (4 clusters) | 0.897, the lowest |
| Equal Weights | ANOVA F = 139.3, R² = 0.54, the highest of the five schemes |

Small numerical differences can occur with other library versions or hardware.
