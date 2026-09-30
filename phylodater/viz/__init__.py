"""
PhyloDater Visualization Module

This module provides visualization capabilities for phylogenetic dating results,
including time trees with geological timescales.

Based on ChronoPhylo visualization framework.
"""

from phylodater.viz.backends import (
    MatplotlibBackend,
    PlotlyBackend,
    PlotterBackend,
    get_backend,
    list_backends,
)
from phylodater.viz.dating_figure import DatingFigure
from phylodater.viz.geo_plot import (
    GeoInterval,
    GeoPlotter,
    create_geo_axis,
    get_default_geo_colors,
    get_geo_abbreviations,
    load_ics_data,
)
from phylodater.viz.tree_plot import (
    NodeCoordinate,
    TreePlotter,
    compute_node_coordinates,
    get_tree_segments,
)

__all__ = [
    "TreePlotter",
    "NodeCoordinate",
    "compute_node_coordinates",
    "get_tree_segments",
    "GeoPlotter",
    "GeoInterval",
    "load_ics_data",
    "get_default_geo_colors",
    "get_geo_abbreviations",
    "create_geo_axis",
    "DatingFigure",
    "PlotterBackend",
    "MatplotlibBackend",
    "PlotlyBackend",
    "get_backend",
    "list_backends",
]
