import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches


def convert_to_ranking(complexity_array):
    """
    Converts a 2D complexity array into a 2D ranking array.

    Args:
        complexity_array (numpy.ndarray): A 2D array of complexities.

    Returns:
        numpy.ndarray: A 2D array of rankings.
    """
    complexity_rankings = np.empty(complexity_array.shape, dtype=int)
    for i in range(complexity_rankings.shape[1]):
        pos = np.argsort(complexity_array[:, i])
        for j, pos in enumerate(pos):
            complexity_rankings[pos, i] = j
    return complexity_rankings


def ranking_plot(ranking_lists, line_colors, marker_size=10, line_width=1):
    """
    plotting a scatter plot of cells (genes) and connecting them with lines based on their complexity ranking

    Args:
        ranking_lists (list of numpy.ndarray): List of ranking arrays.
        line_colors: color of cells (genes).
        marker_size (int, optional): Size of scatter plot markers. Defaults to 10.
        line_width (int, optional): Width of connecting lines. Defaults to 1.

    Returns:
        matplotlib.figure.Figure, matplotlib.axes._subplots.AxesSubplot:
        The created figure and axes.
    """
    fig, ax = plt.subplots()
    indent = 0.8
    for ranking_list, color in zip(ranking_lists, line_colors):
        # Scatter plot
        ax.scatter(np.arange(len(ranking_list)), ranking_list, marker='o', color=color, s=marker_size, zorder=3)
        # Create path for connecting lines
        verts = [(i + d, rank) for i, rank in enumerate(ranking_list) for d in (-indent, 0, indent)][1:-1]
        codes = [Path.MOVETO] + [Path.CURVE4] * (len(verts) - 1)
        path = Path(verts, codes)
        # Add patch for connecting lines
        patch = patches.PathPatch(path, facecolor='none', lw=line_width, edgecolor=color)
        ax.add_patch(patch)

    return fig, ax

