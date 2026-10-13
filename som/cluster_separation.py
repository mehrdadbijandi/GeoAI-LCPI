"""
Internal consistency check: separation of the SOM clusters by the LCPI of the five weighting schemes
(paper: Section 3.5, Table A.2).

For each scheme: one-way ANOVA (F, p, R^2 = share of LCPI variance explained by the clusters),
Welch's ANOVA and the Kruskal-Wallis test (H, p, epsilon^2).

Run som_lcpi.py first (it writes outputs/Enhanced_Dataset.csv with the column SOM_Cluster).
Usage:  python cluster_separation.py
Output: outputs/cluster_separation.csv
"""
import os
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, 'outputs')


def welch_anova(groups):
    k = len(groups); n = np.array([len(g) for g in groups]); m = np.array([g.mean() for g in groups])
    v = np.array([g.var(ddof=1) for g in groups]); w = n / v; mw = (w * m).sum() / w.sum()
    A = (w * (m - mw) ** 2).sum() / (k - 1); t = ((1 - w / w.sum()) ** 2 / (n - 1)).sum()
    F = A / (1 + 2 * (k - 2) / (k ** 2 - 1) * t); df2 = (k ** 2 - 1) / (3 * t)
    return F, df2, stats.f.sf(F, k - 1, df2)


def evaluate(lcpi, cl):
    g = [lcpi[cl == c] for c in np.unique(cl)]
    F, p = stats.f_oneway(*g)
    r2 = sum(len(v) * (v.mean() - lcpi.mean()) ** 2 for v in g) / ((lcpi - lcpi.mean()) ** 2).sum()
    wF, wdf2, wp = welch_anova(g); H, hp = stats.kruskal(*g)
    return dict(ANOVA_F=F, ANOVA_p=p, R2=r2, Welch_F=wF, Welch_df2=wdf2, Welch_p=wp, KW_H=H, KW_p=hp,
                KW_eps2=(H - len(g) + 1) / (len(lcpi) - len(g)))


cl = pd.read_csv(os.path.join(OUT, 'Enhanced_Dataset.csv')).set_index('ID').SOM_Cluster
rows = []
for key, name in [('EQ', 'Equal Weights'), ('CR', 'CRITIC'), ('STD', 'Standard Deviation'), ('SH', 'Shannon Entropy'), ('PCA', 'PCA')]:
    L = pd.read_csv(os.path.join(DATA, f'plus{key}.csv')).set_index('ID').LCPI.loc[cl.index]
    rows.append(dict(Scheme=name, **evaluate(L.values, cl.values)))
t = pd.DataFrame(rows); t.to_csv(os.path.join(OUT, 'cluster_separation.csv'), index=False)
print(t[['Scheme', 'ANOVA_F', 'R2', 'Welch_F', 'KW_H']].round(3).to_string(index=False))
