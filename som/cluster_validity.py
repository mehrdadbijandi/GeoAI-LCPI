"""
Cluster validity indices for Ward clustering of the trained SOM prototypes (paper: Section 3.4.2, Fig. A.1).

Retrains the reference SOM exactly as som_lcpi.py (same data, seeds and settings), then computes for k = 2..12:
  within-cluster sum of squares (elbow), silhouette score and Davies-Bouldin index of the prototype clustering.

Usage:  python cluster_validity.py [input_csv] [shuffle_seed]
Outputs (./outputs): cluster_validity.csv, cluster_validity.png
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]
sys.argv = ['som_lcpi.py'] + (args if args else [os.path.join(HERE, '..', 'data', 'plusEQ.csv'), '71'])
# run the training part of som_lcpi.py (everything before the cluster assignment)
code = open(os.path.join(HERE, 'som_lcpi.py'), encoding='utf-8').read().split('# Assign clusters')[0]
exec(compile(code, 'som_lcpi.py', 'exec'))

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, davies_bouldin_score

W = som.get_weights().reshape(-1, X_scaled.shape[1])
rows = []
for k in range(2, 13):
    lab = AgglomerativeClustering(n_clusters=k, linkage='ward').fit_predict(W)
    wcss = sum(((W[lab == c] - W[lab == c].mean(0)) ** 2).sum() for c in range(k))
    rows.append(dict(k=k, wcss=wcss, silhouette=silhouette_score(W, lab), davies_bouldin=davies_bouldin_score(W, lab)))
t = pd.DataFrame(rows)
t.to_csv(os.path.join(OUT, 'cluster_validity.csv'), index=False)
print(t.round(3).to_string(index=False))

fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for a, col, lab, best in [(ax[0], 'wcss', 'Within-cluster sum of squares', None), (ax[1], 'silhouette', 'Silhouette score', 'max'),
                          (ax[2], 'davies_bouldin', 'Davies-Bouldin index', 'min')]:
    a.plot(t.k, t[col], 'o-', lw=1.4, ms=4)
    a.axvline(4, color='grey', ls='--', lw=0.8)
    a.set_xlabel('Number of clusters'); a.set_ylabel(lab); a.grid(alpha=0.3); a.set_xticks(range(2, 13, 2))
for i, a in enumerate(ax):
    a.set_title('(%s)' % 'abc'[i], loc='left')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'cluster_validity.png'), dpi=300)
