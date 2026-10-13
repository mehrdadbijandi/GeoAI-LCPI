# GeoAI diagnostic profiling with a Self-Organizing Map (SOM)

These scripts group parcels by their profiles of the four normalized indicators (SIʹ, BCʹ, TFʹ, ACʹ) and describe which deficiencies drive their LCPI. The SOM uses no indicator weights. The resulting clusters are used:

- as diagnostic profiles of the parcels (paper, Section 3.4);
- as an internal benchmark for the weighting schemes (paper, Section 3.5).

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
