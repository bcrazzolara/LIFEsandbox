import matplotlib.colors as mcolors
import numpy as np

colors = np.zeros((1000, 4))
colors[:,-1] = 1
colors[:299, 0] = np.linspace(0, 1, 299)
colors[299:, 0] = 1
colors[299:886, 1] = np.linspace(0, 1, 587)
colors[886:, 1] = 1
colors[886:, 2] = np.linspace(0, 1, 114)


colors_inv = np.zeros((1000, 4))
colors_inv[:,-1] = 1
colors_inv[-456:, 2] = np.linspace(1, 0, 456)
colors_inv[:-456, 2] = 1

colors_inv[-800:, 1] = np.linspace(1, 0, 800)
colors_inv[:-800, 1] = 1

#colors_inv[-701:-114, 1] = np.linspace(1, 0, 587)
#colors_inv[:-701, 1] = 1
colors_inv[200:544, 0] = np.linspace(0.3, 0, 344)
colors_inv[0:200, 0] = np.linspace(1, 0.3, 200)

# combine them and build a new colormap
colors_comb = np.vstack((colors_inv, colors))

sat = 0.8
Y = 0.299*colors_comb[:,0] + 0.587*colors_comb[:,1] + 0.114*colors_comb[:,2]
colors_comb[:,:3] = sat * colors_comb[:,:3]  + (1-sat)*Y[:,None]


neg_map = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors_comb[200:1000])
pos_map = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors_comb[-1000:-200])
full_map = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors_comb[200:-200])


