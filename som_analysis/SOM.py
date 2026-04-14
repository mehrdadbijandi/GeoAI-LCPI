import sys
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import BoundaryNorm, Normalize
from sklearn.preprocessing import PowerTransformer
from sklearn.cluster import AgglomerativeClustering
import warnings
import random

warnings.filterwarnings('ignore')

# Set random seeds for reproducibility
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

current_dir = os.path.dirname(os.path.abspath(__file__))

def install_package(package_name):
    """Install Python package if missing."""
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package_name])
        return True
    except Exception as e:
        print(f"Failed to install {package_name}: {e}")
        return False

# Ensure minisom is available
try:
    from minisom import MiniSom
except ImportError:
    install_package('minisom')
    from minisom import MiniSom

print("Loading data...")
csv_path = os.path.join(current_dir, 'input_SOM_EQ.csv')

try:
    df = pd.read_csv(csv_path)
except FileNotFoundError:
    print(f"Error: '{csv_path}' not found.")
    sys.exit(1)

features = ['SI', 'BC', 'TF', 'AC']
X = df[features].values

# Fix skewness to ensure equal feature contribution to Euclidean distances
power_transformer = PowerTransformer(method='yeo-johnson', standardize=True)
X_scaled = power_transformer.fit_transform(X)

# Define grid shape based on data size
map_dim = (9, 13)

class AdaptiveSOMTrainer:
    """Trainer class for SOM with stable decay parameters."""
    
    def __init__(self, X_scaled, map_dim, seed=SEED):
        self.X_scaled = X_scaled
        self.map_dim = map_dim
        self.seed = seed
        self.convergence_history = {
            'iterations': [], 'qe_values': [], 'te_values': [],
            'learning_rates': [], 'radii': []
        }
        
    def calculate_qe_te(self, som, winner_coordinates):
        """Vectorized computation of Quantization Error (QE) and Topographic Error (TE)."""
        weights = som.get_weights()
        n_samples = self.X_scaled.shape[0]
        x_dim, y_dim, n_features = weights.shape
        
        # Calculate average QE
        bmu_indices = winner_coordinates[:, 0] * y_dim + winner_coordinates[:, 1]
        weights_flat = weights.reshape(-1, n_features)
        bmu_weights = weights_flat[bmu_indices]
        avg_qe = np.mean(np.linalg.norm(self.X_scaled - bmu_weights, axis=1))
        
        # Calculate TE by finding two closest neurons
        diff = self.X_scaled[:, np.newaxis, :] - weights_flat[np.newaxis, :, :]
        distances = np.sqrt(np.sum(diff**2, axis=2))
        
        two_closest_indices = np.argpartition(distances, 2, axis=1)[:, :2]
        sample_indices = np.arange(n_samples)[:, np.newaxis]
        two_closest_distances = distances[sample_indices, two_closest_indices]
        
        sort_order = np.argsort(two_closest_distances, axis=1)
        two_closest_indices_sorted = two_closest_indices[sample_indices, sort_order]
        
        bmu_coords_flat = two_closest_indices_sorted[:, 0]
        second_bmu_coords_flat = two_closest_indices_sorted[:, 1]
        
        bmu_coords_2d = np.column_stack((bmu_coords_flat // y_dim, bmu_coords_flat % y_dim))
        second_bmu_coords_2d = np.column_stack((second_bmu_coords_flat // y_dim, second_bmu_coords_flat % y_dim))
        
        # Check topological preservation (hexagonal adjacency)
        te_count = 0
        for idx in range(n_samples):
            x1, y1 = bmu_coords_2d[idx]
            x2, y2 = second_bmu_coords_2d[idx]
            dx, dy = abs(x1 - x2), abs(y1 - y2)
            
            is_adjacent = False
            if dx == 0 and dy == 1:
                is_adjacent = True
            elif dx == 1:
                if x1 % 2 == 0:
                    if y2 == y1 or y2 == y1 - 1: is_adjacent = True
                else:
                    if y2 == y1 or y2 == y1 + 1: is_adjacent = True
            
            if not is_adjacent:
                te_count += 1
                
        return avg_qe, te_count / n_samples
        
    def train(self, check_early_stopping=False, stop_at_iteration=None):
        """Execute SOM training with PCA initialization and defined decay schedule."""
        np.random.seed(self.seed)
        random.seed(self.seed)
        
        initial_radius, final_radius = 4.0, 0.1
        initial_learning_rate, final_learning_rate = 0.5, 0.01
        
        som = MiniSom(x=self.map_dim[0], y=self.map_dim[1],
                    input_len=self.X_scaled.shape[1], sigma=initial_radius,
                    learning_rate=initial_learning_rate, neighborhood_function='gaussian',
                    topology='hexagonal', random_seed=self.seed)
        
        som.random_weights_init(self.X_scaled)
        
        max_iterations = 1000000
        n_samples = len(self.X_scaled)
        indices = list(range(n_samples))
        loop_limit = stop_at_iteration if stop_at_iteration else max_iterations
        stopped_early = False
        
        for iteration in range(loop_limit):
            if iteration % n_samples == 0:
                random.shuffle(indices)
            
            progress = iteration / max_iterations
            current_radius = initial_radius * (final_radius / initial_radius) ** progress
            current_lr = initial_learning_rate * (final_learning_rate / initial_learning_rate) ** progress
            
            som._sigma = current_radius
            som._learning_rate = current_lr

            sample_idx = indices[iteration % n_samples]
            x = self.X_scaled[sample_idx]
            som.update(x, som.winner(x), 0, 1)
            
            if (iteration + 1) % 1000 == 0:
                winner_coords = np.array([som.winner(x) for x in self.X_scaled])
                current_qe, current_te = self.calculate_qe_te(som, winner_coords)
                
                self.convergence_history['iterations'].append(iteration + 1)
                self.convergence_history['qe_values'].append(current_qe)
                self.convergence_history['te_values'].append(current_te)
                self.convergence_history['learning_rates'].append(current_lr)
                self.convergence_history['radii'].append(current_radius)
                
                if check_early_stopping and current_te < 0.06 and current_qe < 1:
                    stopped_early = True
                    break

        winner_coords = np.array([som.winner(x) for x in self.X_scaled])
        final_qe, final_te = self.calculate_qe_te(som, winner_coords)
        return som, winner_coords, final_qe, final_te, stopped_early

def perform_smart_training(X_scaled, map_dim):
    """Manage training phases to find optimal stopping points."""
    trainer = AdaptiveSOMTrainer(X_scaled, map_dim)
    print("Training SOM and searching for optimal convergence point...")
    som, winner_coords, final_qe, final_te, stopped_early = trainer.train(check_early_stopping=True)

    if stopped_early:
        return trainer, som, winner_coords, final_qe, final_te

    history_df = pd.DataFrame({
        'iter': trainer.convergence_history['iterations'],
        'qe': trainer.convergence_history['qe_values'],
        'te': trainer.convergence_history['te_values']
    })

    min_te_idx = history_df['te'].idxmin()
    min_qe_idx = history_df['qe'].idxmin()
    
    if min_te_idx == min_qe_idx:
        best_iter = history_df.loc[min_te_idx, 'iter']
    else:
        candidates = history_df[history_df['te'] < 0.06]
        if not candidates.empty:
            best_iter = candidates.loc[candidates['qe'].idxmin(), 'iter']
        else:
            best_iter = history_df.loc[min_te_idx, 'iter']

    trainer_retrain = AdaptiveSOMTrainer(X_scaled, map_dim)
    som, winner_coords, final_qe, final_te, _ = trainer_retrain.train(check_early_stopping=False, stop_at_iteration=int(best_iter))
    
    return trainer_retrain, som, winner_coords, final_qe, final_te

# ----------------- Execution & Post-processing -----------------
trainer, som, winner_coordinates, final_qe, final_te = perform_smart_training(X_scaled, map_dim)

df['SOM_X'] = winner_coordinates[:, 0]
df['SOM_Y'] = winner_coordinates[:, 1]
df['SOM_Neuron_ID'] = df['SOM_X'] * map_dim[1] + df['SOM_Y']

# Assign Clusters using Agglomerative Clustering
n_clusters = 4
neuron_weights = som.get_weights().reshape(-1, X_scaled.shape[1])
ac = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
neuron_clusters = ac.fit_predict(neuron_weights)

neuron_to_cluster = {i * map_dim[1] + j: int(neuron_clusters[i * map_dim[1] + j]) 
                     for i in range(map_dim[0]) for j in range(map_dim[1])}
df['SOM_Cluster'] = df['SOM_Neuron_ID'].map(neuron_to_cluster)

# Inverse transform weights for visualization mapping
weights_2d = som.get_weights().reshape(-1, X_scaled.shape[1])
weights_original = power_transformer.inverse_transform(weights_2d)
som_weights_original = weights_original.reshape(map_dim[0], map_dim[1], X_scaled.shape[1])

# Aggregate data for mapping visualizations
lcp_plane = np.zeros((map_dim[0], map_dim[1]))
count_plane = np.zeros((map_dim[0], map_dim[1]))
for i, (x, y) in enumerate(winner_coordinates):
    lcp_plane[x, y] += df.loc[i, 'LCPI']
    count_plane[x, y] += 1
lcp_plane_avg = np.full(lcp_plane.shape, np.nan)
np.divide(lcp_plane, count_plane, out=lcp_plane_avg, where=count_plane!=0)

cluster_grid = np.zeros((map_dim[0], map_dim[1]), dtype=int)
for i in range(map_dim[0]):
    for j in range(map_dim[1]):
        cluster_grid[i, j] = neuron_to_cluster[i * map_dim[1] + j]

box_data = [df[df['SOM_Cluster'] == i]['LCPI'].values for i in range(n_clusters)]
overall_mean = df['LCPI'].mean()
cluster_colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow']

# ----------------- Hexagonal Visualization Helpers -----------------
def hex_center_coords(grid_shape, radius=0.5):
    """Compute center coordinates for each hex cell."""
    dx = np.sqrt(3) * radius
    dy = 1.5 * radius
    centers = {}
    for i in range(grid_shape[0]):
        for j in range(grid_shape[1]):
            x = j * dx
            y = i * dy
            if i % 2 == 1:
                x += dx / 2
            centers[(i, j)] = (x, y)
    return centers

def hex_vertices(center, radius=0.5, orientation=0.0):
    """Return the 6 vertices (Nx2 array) of a flat-top hexagon centered at `center`."""
    cx, cy = center
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False) + orientation
    verts = np.stack([cx + radius * np.cos(angles), cy + radius * np.sin(angles)], axis=1)
    return verts

def draw_cluster_boundaries(ax, cluster_grid, radius=0.5, linewidth=3.0, alpha=1.0, cmap_name='tab20c'):
    """Draw cluster boundaries along the hexagon edges."""
    cluster_grid_vis = cluster_grid
    rows, cols = cluster_grid_vis.shape
    centers = hex_center_coords((rows, cols), radius=radius)
    verts_map = {}
    flat_id_map = {}
    
    for i in range(rows):
        for j in range(cols):
            center = centers[(i, j)]
            verts = hex_vertices(center, radius=radius, orientation=0.0)
            verts_map[(i, j)] = verts
            flat_id_map[(i, j)] = i * cols + j

    center_items = list(centers.items())
    coords = np.array([c for _, c in center_items])
    keys = [k for k, _, in center_items]
    dx = np.sqrt(3) * radius
    thresh = dx + 1e-6
    boundary_color = 'black'

    for idx, (cell_key, center) in enumerate(center_items):
        i, j = cell_key
        cell_cluster = int(cluster_grid_vis[i, j])
        dists = np.hypot(coords[:, 0] - center[0], coords[:, 1] - center[1])
        neighbor_indices = np.where((dists > 0) & (dists <= thresh))[0]
        
        for nidx in neighbor_indices:
            neigh_key = keys[nidx]
            if flat_id_map[neigh_key] <= flat_id_map[cell_key]:
                continue
            ni, nj = neigh_key
            neigh_cluster = int(cluster_grid_vis[ni, nj])
            if neigh_cluster == cell_cluster:
                continue

            verts1 = verts_map[cell_key]
            verts2 = verts_map[neigh_key]
            shared = []
            for v1 in verts1:
                for v2 in verts2:
                    if np.allclose(v1, v2, atol=1e-6, rtol=1e-6):
                        shared.append(tuple(v1))
            if len(shared) >= 2:
                a = shared[0]
                b = shared[1]
                ax.plot([a[0], b[0]], [a[1], b[1]],
                        color=boundary_color, linewidth=linewidth, alpha=alpha, solid_capstyle='butt', zorder=20)
            else:
                cx, cy = center
                ncx, ncy = centers[neigh_key]
                mx, my = (cx + ncx) / 2.0, (cy + ncy) / 2.0
                dx_seg = (ncx - cx)
                dy_seg = (ncy - cy)
                seg_len = radius * 0.6
                norm = np.hypot(dx_seg, dy_seg) + 1e-12
                px = -dy_seg / norm * seg_len
                py = dx_seg / norm * seg_len
                ax.plot([mx - px, mx + px], [my - py, my + py],
                        color=boundary_color, linewidth=linewidth, alpha=alpha, solid_capstyle='butt', zorder=20)

def create_hexagonal_grid(ax, grid_shape, values, cmap='viridis', edgecolor='black', linewidth=0.5, radius=0.5):
    hexagons, vals_list = [], []
    dx, dy = np.sqrt(3) * radius, 1.5 * radius

    for i in range(grid_shape[0]):
        for j in range(grid_shape[1]):
            x = j * dx + (dx / 2 if i % 2 == 1 else 0)
            hexagons.append(RegularPolygon((x, i * dy), numVertices=6, radius=radius,
                                         orientation=0, edgecolor=edgecolor, linewidth=linewidth))
            if hasattr(values, 'shape') and len(values.shape) == 2:
                vals_list.append(values[i, j])
            else:
                idx = i * grid_shape[1] + j
                vals_list.append(values[idx] if values is not None and idx < len(values) else np.nan)

    vals_arr = np.array(vals_list)
    collection = PatchCollection(hexagons, match_original=True, linewidth=linewidth)
    nan_mask = np.isnan(vals_arr)
    non_nan_vals = vals_arr[~nan_mask]
    
    if non_nan_vals.size > 0:
        unique_vals = np.unique(non_nan_vals)
        is_discrete = np.allclose(unique_vals, np.round(unique_vals)) and unique_vals.size <= 20

        if is_discrete:
            unique_sorted = np.sort(unique_vals)
            mapping = {val: idx for idx, val in enumerate(unique_sorted)}
            mapped_array = np.full(vals_arr.shape, -1, dtype=float)
            for i, val in enumerate(vals_arr):
                if not np.isnan(val): mapped_array[i] = mapping[val]

            cmap_obj = plt.get_cmap(cmap, len(unique_sorted))
            norm = BoundaryNorm(boundaries=np.arange(len(unique_sorted) + 1) - 0.5, ncolors=len(unique_sorted))
            
            face_rgba = np.ones((len(vals_arr), 4))
            edge_rgba = np.ones((len(vals_arr), 4))
            
            for i, idx in enumerate(mapped_array):
                if idx != -1:
                    face_rgba[i] = edge_rgba[i] = cmap_obj(norm(idx))
                else:
                    face_rgba[i], edge_rgba[i] = [1.0, 1.0, 1.0, 1.0], [0.827, 0.827, 0.827, 1.0]
                    
            collection.set_facecolor(face_rgba)
            collection.set_edgecolor(edge_rgba)
            collection.set_array(mapped_array)
            collection.set_cmap(cmap_obj)
            collection.set_norm(norm)
            collection._is_discrete = True
            collection._category_values = unique_sorted
        else:
            vmin, vmax = np.nanmin(vals_arr), np.nanmax(vals_arr)
            norm = Normalize(vmin=vmin, vmax=vmax)
            cmap_obj = plt.get_cmap(cmap)
            face_rgba = cmap_obj(norm(vals_arr))
            edge_rgba = face_rgba.copy()
            
            for i, is_nan in enumerate(nan_mask):
                if is_nan:
                    face_rgba[i], edge_rgba[i] = [1.0, 1.0, 1.0, 1.0], [0.827, 0.827, 0.827, 1.0]
            
            collection.set_facecolor(face_rgba)
            collection.set_edgecolor(edge_rgba)
            collection.set_array(vals_arr)
            collection.set_cmap(cmap_obj)
            collection.set_norm(norm)
            collection._is_discrete = False

    ax.add_collection(collection)
    ax.set_xlim(-dx / 2, (grid_shape[1] - 1) * dx + dx)
    ax.set_ylim(-radius * 1.1, (grid_shape[0] - 1) * dy + radius * 1.1)
    ax.set_aspect('equal')
    ax.axis('off')
    return collection

def plot_hexagonal_matrix(ax, matrix, title, cmap='viridis', add_colorbar=True):
    collection = create_hexagonal_grid(ax, matrix.shape, matrix, cmap=cmap)
    ax.set_title(title, fontweight='bold')
    if add_colorbar and hasattr(collection, 'get_array') and collection.get_array() is not None:
        cbar = plt.colorbar(collection, ax=ax, shrink=0.8)
        if getattr(collection, '_is_discrete', False):
            cats = collection._category_values
            cbar.set_ticks(np.arange(len(cats)))
            cbar.set_ticklabels([str(int(v)) for v in cats])
    return collection

# ----------------- Output 1: Complete SOM Analysis -----------------
print("Generating complete SOM analysis visualization...")
fig1 = plt.figure(figsize=(15, 10))
fig1.suptitle('SOM Analysis (LCPI) - Hexagonal Topology', fontsize=16, fontweight='bold')

# 1. U-Matrix
ax1 = plt.subplot(2, 3, 1)
u_matrix = som.distance_map()
plot_hexagonal_matrix(ax1, u_matrix.T, 'U-Matrix (Distance Map)', cmap='bone_r')

# 2-5. Component Planes
feature_names = ['Shape Suitability', 'Boundary Complexity', 'Accessibility', 'Ownership Fragmentation']
for idx in range(4):
    ax = plt.subplot(2, 3, 2 + idx)
    weight_matrix = som_weights_original[:, :, idx]
    plot_hexagonal_matrix(ax, weight_matrix.T, f'{feature_names[idx]}', cmap='coolwarm', add_colorbar=True)

# 6. Correlation Matrix
ax6 = plt.subplot(2, 3, 6)
som_weights_flat = som_weights_original.reshape(-1, len(features))
corr_df = pd.DataFrame(som_weights_flat, columns=features).dropna()
corr_matrix = corr_df.corr()

im6 = ax6.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
plt.colorbar(im6, ax=ax6, fraction=0.046, pad=0.04)
ax6.set_xticks(np.arange(len(corr_matrix.columns)))
ax6.set_yticks(np.arange(len(corr_matrix.columns)))
ax6.set_xticklabels(corr_matrix.columns)
ax6.set_yticklabels(corr_matrix.columns)
plt.setp(ax6.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

for i in range(len(corr_matrix.columns)):
    for j in range(len(corr_matrix.columns)):
        ax6.text(j, i, f"{corr_matrix.iloc[i, j]:.2f}",
                 ha="center", va="center",
                 color="black" if abs(corr_matrix.iloc[i, j]) < 0.5 else "white")
ax6.set_title('SOM Component Correlation', fontweight='bold')

plt.tight_layout()
fig1.savefig(os.path.join(current_dir, 'SOM_Complete_Analysis_Hexagonal.png'), dpi=300, bbox_inches='tight')
plt.close(fig1)

# ----------------- Output 2: SOM Clustering Results -----------------
print("Generating SOM clustering results...")
fig2, axes = plt.subplots(2, 2, figsize=(14, 12))
fig2.suptitle('SOM Clustering Results - Hexagonal Topology', fontsize=16, fontweight='bold', y=1.02)

# Top-Left: LCPI Component Plane with Boundaries
ax1 = axes[0, 0]
plot_hexagonal_matrix(ax1, lcp_plane_avg.T, 'LCPI Component Plane', cmap='RdYlGn')
draw_cluster_boundaries(ax1, cluster_grid.T, radius=0.5, linewidth=2.5, alpha=0.95)

# Top-Right: SOM Clusters
ax2 = axes[0, 1]
col2 = create_hexagonal_grid(ax2, cluster_grid.T.shape, cluster_grid.T, cmap='tab20c', edgecolor='black', linewidth=1)
ax2.set_title(f'SOM Clusters (K={n_clusters})', fontweight='bold', fontsize=12)
if getattr(col2, '_is_discrete', False):
    cbar2 = plt.colorbar(col2, ax=ax2, shrink=0.8)
    cbar2.set_ticks(np.arange(len(col2._category_values)))
    cbar2.set_ticklabels([f'Cluster {int(v)}' for v in col2._category_values])

# Bottom-Left: Combined Box and Violin Plot 
ax3 = axes[1, 0]
violin_parts = ax3.violinplot(box_data, showmeans=False, showmedians=False)
for i, pc in enumerate(violin_parts['bodies']):
    pc.set_facecolor(cluster_colors[i])
    pc.set_edgecolor('black')
    pc.set_alpha(0.4)

for part in ('cbars', 'cmins', 'cmaxes', 'cmeans', 'cmedians'):
    if part in violin_parts: violin_parts[part].set_linewidth(0)

bp = ax3.boxplot(box_data, positions=np.arange(1, n_clusters + 1), widths=0.15, patch_artist=True)
for patch, color in zip(bp['boxes'], cluster_colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('black')
    patch.set_alpha(0.9)
    patch.set_linewidth(1.5)

for element in ['whiskers', 'caps', 'medians', 'fliers']:
    for line in bp[element]:
        line.set_color('black')
        line.set_linewidth(1.5)
for flier in bp['fliers']:
    flier.set(marker='o', color='red', alpha=0.7, markersize=5)

for i, data in enumerate(box_data):
    ax3.scatter(i + 1, np.mean(data), color='red', zorder=5, s=120, marker='D', edgecolor='black', linewidth=1.5, label='Mean' if i == 0 else "")
    q1, median, q3 = np.percentile(data, [25, 50, 75])
    ax3.hlines(y=[q1, median, q3], xmin=i+0.8, xmax=i+1.2, colors=['gray', 'black', 'gray'], linewidths=[1.5, 2, 1.5], linestyles=[':', '-', ':'], zorder=4)
    x_jitter = np.random.normal(i + 1, 0.08, size=len(data))
    ax3.scatter(x_jitter, data, alpha=0.3, color='black', s=25, edgecolor='none', zorder=2)

ax3.axhline(y=overall_mean, color='blue', linestyle='--', alpha=0.8, linewidth=2.5, label=f'Overall Mean: {overall_mean:.3f}')
ax3.set_xticks(np.arange(1, n_clusters + 1))
ax3.set_xticklabels([f'Cluster {i}\n(n={len(box_data[i])})' for i in range(n_clusters)], fontsize=11)
ax3.set_title('Combined Box & Violin Plot: LCPI Distribution', fontweight='bold', fontsize=12)
ax3.set_ylabel('LCPI Value', fontsize=10)
ax3.set_xlabel('SOM Cluster', fontsize=10)
ax3.legend(loc='upper right', fontsize=9)
ax3.grid(True, alpha=0.3, linestyle='--')

info_text = f'Samples: {len(df)}\nMean LCPI: {overall_mean:.3f}\nQE: {final_qe:.4f}\nTE: {final_te:.4f}'
ax3.text(0.02, 0.98, info_text, transform=ax3.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Bottom-Right: Radar Chart
ax4 = axes[1, 1]
cluster_means = [[df[df['SOM_Cluster'] == i][f].mean() for f in features] for i in range(n_clusters)]
cluster_means_array = np.array(cluster_means)
for j in range(len(features)):
    min_val, max_val = cluster_means_array[:, j].min(), cluster_means_array[:, j].max()
    if max_val > min_val:
        cluster_means_array[:, j] = (cluster_means_array[:, j] - min_val) / (max_val - min_val)

angles = np.linspace(0, 2 * np.pi, len(features), endpoint=False).tolist() + [0]
for i in range(n_clusters):
    values = cluster_means_array[i].tolist() + [cluster_means_array[i][0]]
    ax4.plot(angles, values, 'o-', linewidth=2, label=f'Cluster {i}')
    ax4.fill(angles, values, alpha=0.25)

ax4.set_xticks(angles[:-1])
ax4.set_xticklabels(features, fontsize=10)
ax4.set_ylim(0, 1)
ax4.set_title('Normalized Feature Means by Cluster', fontweight='bold', fontsize=12)
ax4.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
ax4.grid(True)

plt.tight_layout()
fig2.savefig(os.path.join(current_dir, 'SOM_Clustering_Results_Hexagonal.png'), dpi=300, bbox_inches='tight')
plt.close(fig2)

# ----------------- Output 3: Enhanced CSV Export -----------------
print("Exporting enhanced dataset...")
enhanced_df = df.copy()
enhanced_df['SOM_Final_TE'] = final_te
enhanced_df['SOM_Final_QE'] = final_qe

output_filename = os.path.join(current_dir, 'Enhanced_Dataset.csv')
enhanced_df.to_csv(output_filename, index=False)

print("Processing complete. All required files have been generated and saved.")