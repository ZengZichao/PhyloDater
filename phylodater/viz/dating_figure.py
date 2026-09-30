"""
DatingFigure - Main class for creating phylogenetic time tree visualizations.

Based on ChronoPhylo, adapted for PhyloDater's DatingResult.

Supports multiple visualization backends via the backend parameter:
- 'matplotlib': Default matplotlib-based backend
- 'plotly': Interactive Plotly-based backend (future extension)
"""

import math
from pathlib import Path
from types import TracebackType
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

from phylodater.infrastructure.logging import get_logger
from phylodater.models import DatingResult, PhylogeneticTree
from phylodater.viz.backends import PlotterBackend, get_backend
from phylodater.viz.geo_plot import GeoPlotter
from phylodater.viz.tree_plot import (
    TreePlotter,
    as_ete_tree,
    is_ete4_node,
    resolve_root_age,
)

if TYPE_CHECKING:
    pass

logger = get_logger()


class DatingFigure:
    """
    Main class for creating phylogenetic time tree visualizations.

    This class provides a fluent interface for building publication-quality
    phylogenetic tree figures with integrated geological timescales.

    ``root_age`` 不需要手填：未提供时按树自身的分支长度推断（root-to-tip 最长
    路径），因此下面这个文档示例是完整可用的（审阅项 B-13 修复前的行为是把它
    当成 0，于是所有节点塌到 x=0 且没有任何提示）。

    Example
    -------
    >>> from phylodater.viz import DatingFigure
    >>> from phylodater.models import PhylogeneticTree, DatingResult
    >>> tree = PhylogeneticTree.from_file("tree.nwk")
    >>> fig = (DatingFigure()
    ...     .add_tree(tree)
    ...     .add_geo_scale()
    ...     .render()
    ...     .save("output.pdf"))

    显式定标用 ``add_tree(tree, root_age=500)``；带定年结果时用
    ``add_result(result)`` 把节点年龄写进树，之后所有坐标都以注解年龄为准。
    """

    def __init__(
        self,
        width: float = 10.0,
        height: float = 8.0,
        dpi: int = 300,
        backend: str = "matplotlib",
    ) -> None:
        """
        Initialize a DatingFigure instance.

        Parameters
        ----------
        width : float
            Figure width in inches.
        height : float
            Figure height in inches.
        dpi : int
            Resolution for raster formats.
        backend : str
            Visualization backend to use ('matplotlib' or 'plotly').
        """
        self._width = width
        self._height = height
        self._dpi = dpi
        self._backend_name = backend
        if backend != "matplotlib":
            raise NotImplementedError(
                f"backend='{backend}' is not wired into the plotting path yet; "
                "only 'matplotlib' is currently supported."
            )
        self._backend: PlotterBackend = get_backend(backend)

        self._tree: Optional[Any] = None
        self._tree_plotter: Optional[TreePlotter] = None
        self._geo_plotter: Optional[GeoPlotter] = None

        self._figure: Optional[Figure] = None
        self._ax_tree: Optional[Any] = None
        self._ax_geo: Optional[Any] = None
        self._figure_created: bool = False

        self._geo_data: Optional[pd.DataFrame] = None
        self._show_geo_scale: bool = False
        self._custom_geo_path: Optional[Union[str, Path]] = None

        self._root_age: Optional[float] = None
        self._relative_time_mode: bool = False

        self._node_coords: Optional[Dict[Any, Any]] = None
        self._tree_dataframe: Optional[pd.DataFrame] = None

        self._color_mapping: Dict[str, Any] = {}
        self._size_mapping: Dict[str, Any] = {}

        self._metadata: Optional[pd.DataFrame] = None
        self._metadata_path: Optional[Union[str, Path]] = None

    @property
    def ax_tree(self) -> Any:
        """Get the tree axes."""
        if not self._figure_created or self._ax_tree is None:
            raise RuntimeError(
                "Please call render() first to initialize coordinates and axes."
            )
        return self._ax_tree

    @property
    def ax_geo(self) -> Any:
        """Get the geological timescale axes."""
        if not self._figure_created or self._ax_geo is None:
            raise RuntimeError(
                "Please call render() first to initialize coordinates and axes."
            )
        return self._ax_geo

    @property
    def figure(self) -> Figure:
        """Get the matplotlib Figure."""
        if not self._figure_created or self._figure is None:
            raise RuntimeError("Please call render() first to create the figure.")
        return self._figure

    def add_tree(
        self,
        tree: Any,
        root_age: Optional[float] = None,
        parse_labels: bool = True,
        metadata: Optional[Union[str, Path, pd.DataFrame]] = None,
    ) -> "DatingFigure":
        """
        Add a phylogenetic tree to the figure.

        Parameters
        ----------
        tree : Any
            The tree object (ete4.Tree, PhylogeneticTree, or Newick string).
        root_age : float, optional
            The root age in Ma. If not provided, it is inferred from the tree:
            the root node's ``age`` annotation first, otherwise the longest
            root-to-tip branch path (see :func:`resolve_root_age`). When the
            tree carries no branch lengths and no age annotation at all, a
            ``ValueError`` is raised rather than silently treating the root as
            0 Ma (review item B-13).
        parse_labels : bool
            Whether to parse taxonomic information from labels.
            Currently not implemented; reserved for future use.
        metadata : Union[str, Path, pd.DataFrame], optional
            Metadata file path or DataFrame for additional node information.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        logger.debug("Loading tree...")

        if metadata is not None:
            if isinstance(metadata, (str, Path)):
                self._metadata_path = Path(metadata)
                self._metadata = pd.read_csv(metadata)
            elif isinstance(metadata, pd.DataFrame):
                self._metadata = metadata

        if isinstance(tree, PhylogeneticTree):
            tree = as_ete_tree(tree)
        elif isinstance(tree, str):
            # 可视化层**不自行导入 ete**（审阅项 A-5/B-10：ete3 自 Python 3.13 起
            # 因 stdlib ``cgi`` 被移除而无法导入，且它是可选的 GPLv3 后端）。
            # Newick 字符串一律先交给模型层解析，再经 ``as_ete_tree()`` 取树后端，
            # 因此 models/tree.py 换成 ete4/dendropy 后端时本层无需改动；
            # 后端不可用时由 ``as_ete_tree`` 给出指名道姓的 ImportError。
            tree = as_ete_tree(PhylogeneticTree.from_newick(tree))

        # ete3/ete4 兼容性处理
        # ete3的树本身就是根节点，ete4的树有root属性
        root_node = tree.root if hasattr(tree, "root") else tree
        _is_ete4 = is_ete4_node(root_node)

        def _get_node_attr(node: Any, attr: str, default: Any = None) -> Any:
            """获取节点属性（兼容ete3和ete4）"""
            if _is_ete4:
                return node.props.get(attr, default)
            else:
                return getattr(node, attr, default)

        def _set_node_attr(node: Any, attr: str, value: Any) -> None:
            """设置节点属性（兼容ete3和ete4）"""
            if _is_ete4:
                node.add_prop(attr, value)
            else:
                node.add_feature(attr, value)

        # B-13：root_age 缺省时真实推断，而不是当成 0。
        root_age = resolve_root_age(
            root_node, root_age, context="DatingFigure.add_tree"
        )
        self._root_age = root_age

        for node in tree.traverse():
            is_root_node = node.is_root() if callable(node.is_root) else node.is_root
            if is_root_node:
                _set_node_attr(node, "age", root_age)
                break

        self._tree = tree
        self._get_node_attr = _get_node_attr
        self._set_node_attr = _set_node_attr

        if _get_node_attr(root_node, "relative_time_mode"):
            self._relative_time_mode = True
            logger.warning(
                "Tree is in relative time mode. Geological timescale will be disabled."
            )

        logger.debug("Tree loaded successfully")
        return self

    def add_result(
        self,
        result: DatingResult,
        method_name: Optional[str] = None,
    ) -> "DatingFigure":
        """
        Add a DatingResult to the figure.

        This will annotate the tree with node ages from the dating result.

        Parameters
        ----------
        result : DatingResult
            The dating result to visualize.
        method_name : str, optional
            The method to use if result contains multiple methods.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._tree is None:
            raise RuntimeError("Please call add_tree() first before adding results.")

        if method_name:
            logger.warning(
                f"method_name parameter is not supported for single DatingResult, ignoring '{method_name}'"
            )

        if len(result.node_ages) == 0:
            raise ValueError("DatingResult has no node ages.")

        for node_name, age_estimate in result.node_ages.items():
            matched = False
            for node in self._tree.traverse():
                if node.name == node_name:
                    self._set_node_attr(node, "age", age_estimate.mean_age)
                    if age_estimate.ci_lower is not None:
                        self._set_node_attr(node, "ci_lower", age_estimate.ci_lower)
                    if age_estimate.ci_upper is not None:
                        self._set_node_attr(node, "ci_upper", age_estimate.ci_upper)
                    matched = True
                    break
            if not matched:
                logger.warning(
                    f"Node '{node_name}' from DatingResult not found in tree"
                )

        logger.debug("Added dating result to tree")
        return self

    def add_geo_scale(
        self,
        custom_data: Optional[Union[str, Path, pd.DataFrame]] = None,
        show_labels: bool = True,
        label_rotation: int = 90,
        show_age_ticks: bool = True,
        age_tick_interval: float = 10.0,
        height_ratio: float = 0.15,
        split_rows: bool = True,
    ) -> "DatingFigure":
        """
        Add a geological timescale to the figure.

        Parameters
        ----------
        custom_data : Union[str, Path, pd.DataFrame], optional
            Custom timescale data file path or DataFrame.
        show_labels : bool
            Whether to show interval labels.
        label_rotation : int
            Rotation angle for labels.
        show_age_ticks : bool
            Whether to show age tick marks.
        age_tick_interval : float
            Interval between age ticks in Ma.
        height_ratio : float
            Height ratio of the geological axis.
        split_rows : bool
            Whether to split into two rows (Era on top, Eon on bottom).

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._relative_time_mode:
            logger.warning(
                "Cannot add geological timescale in relative time mode. "
                "Use absolute ages or provide root_age parameter."
            )
            return self

        self._show_geo_scale = True

        if isinstance(custom_data, (str, Path)):
            self._custom_geo_path = Path(custom_data)
        elif isinstance(custom_data, pd.DataFrame):
            self._geo_data = custom_data

        self._geo_plotter = GeoPlotter(
            show_labels=show_labels,
            label_rotation=label_rotation,
            show_age_ticks=show_age_ticks,
            age_tick_interval=age_tick_interval,
            height_ratio=height_ratio,
            split_rows=split_rows,
        )

        if self._custom_geo_path is not None:
            self._geo_data = self._geo_plotter.load_data(self._custom_geo_path)
        else:
            self._geo_data = self._geo_plotter.load_data()

        logger.debug("Geological timescale configured")
        return self

    def map_color(
        self,
        column_name: str,
        cmap: str = "viridis",
        palette: Optional[List[str]] = None,
        target: str = "branch",
    ) -> "DatingFigure":
        """
        Map colors to tree elements based on a data column.

        Parameters
        ----------
        column_name : str
            The column name to use for color mapping.
        cmap : str
            Matplotlib colormap name for continuous variables.
        palette : List[str], optional
            Color palette for categorical variables.
        target : str
            Target to color: "branch", "label", or "both".

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if target not in ("branch", "label", "both"):
            raise ValueError(
                f"Invalid target: {target}. Must be 'branch', 'label', or 'both'"
            )

        self._color_mapping = {
            "column": column_name,
            "cmap": cmap,
            "palette": palette,
            "target": target,
        }

        if self._tree_plotter is not None:
            self._tree_plotter.map_color(column_name, cmap, palette, target)

        logger.debug(f"Color mapping configured for column: {column_name}")
        return self

    def map_size(
        self,
        column_name: str,
        size_range: Tuple[float, float] = (2.0, 8.0),
    ) -> "DatingFigure":
        """
        Map sizes to tree elements based on a data column.

        Parameters
        ----------
        column_name : str
            The column name to use for size mapping.
        size_range : Tuple[float, float]
            Range for size scaling.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        self._size_mapping = {
            "column": column_name,
            "range": size_range,
        }

        if self._tree_plotter is not None:
            self._tree_plotter.map_size(column_name, size_range)

        logger.debug(f"Size mapping configured for column: {column_name}")
        return self

    def render(
        self,
        geo_data: Optional[pd.DataFrame] = None,
        **kwargs: Any,
    ) -> "DatingFigure":
        """
        Render the figure.

        This method creates the matplotlib figure and axes, computes node
        coordinates, and draws the tree and optional geological timescale.

        Parameters
        ----------
        geo_data : pd.DataFrame, optional
            Geological timescale data. Overrides any previously loaded data.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._tree is None:
            raise RuntimeError("No tree has been added. Call add_tree() first.")

        if geo_data is not None:
            self._geo_data = geo_data
            self._show_geo_scale = True

        self._create_figure()

        self._plot()

        logger.debug("Figure rendered successfully")
        return self

    def _create_figure(self) -> None:
        """Create the matplotlib figure and axes."""
        if self._figure_created:
            logger.debug("Reusing existing figure")
            return

        self._figure = plt.figure(
            figsize=(self._width, self._height),
            dpi=self._dpi,
        )

        if (
            self._show_geo_scale
            and self._geo_plotter is not None
            and not self._relative_time_mode
        ):
            gs = GridSpec(
                2,
                1,
                height_ratios=[
                    1 - self._geo_plotter.height_ratio,
                    self._geo_plotter.height_ratio,
                ],
                hspace=0.0,
            )
            self._ax_tree = self._figure.add_subplot(gs[0])
            self._ax_geo = self._figure.add_subplot(gs[1])
        else:
            self._ax_tree = self._figure.add_subplot(111)
            self._ax_geo = None

        self._figure_created = True
        logger.debug("Figure and axes created")

    def _plot(self) -> None:
        """Execute the actual plotting logic."""
        self._tree_plotter = TreePlotter(
            branch_color="#2C3E50",
            branch_linewidth=2.0,
            leaf_color="#E74C3C",
            leaf_size=8.0,
            internal_node_size=4.0,
            show_internal_labels=True,
            show_leaf_labels=True,
            label_offset=1.5,
        )

        color_column = self._color_mapping.get("column")
        color_cmap = self._color_mapping.get("cmap", "viridis")
        color_palette = self._color_mapping.get("palette")
        size_column = self._size_mapping.get("column")
        size_range = self._size_mapping.get("range", (2.0, 8.0))

        self._node_coords = self._tree_plotter.plot(
            ax=self._ax_tree,
            tree=self._tree,
            root_age=self._root_age,
            color_column=color_column,
            size_column=size_column,
            cmap=color_cmap,
            palette=color_palette,
            size_range=size_range,
        )

        # 退化时间轴一律拒绝出图（B-13 的后果链：root_age=0 时所有节点落在 x=0，
        # 图看上去"正常"却没有任何信息）。即使没开地质时标也照样检查。
        age_range = self._get_age_range()

        if (
            self._show_geo_scale
            and self._ax_geo is not None
            and self._geo_plotter is not None
        ):
            self._geo_plotter.plot(
                ax=self._ax_geo,
                age_range=age_range,
                geo_data=self._geo_data,
            )

            self._synchronize_axes()

    def _get_age_range(self) -> Tuple[float, float]:
        """Get the age range from the tree (descending: old on the left).

        Raises
        ------
        ValueError
            没有节点坐标、坐标里没有有限年龄、或所有节点落在同一时间点。
            此前这些情形会返回 (100, 0) / (0, 0) 之类的兜底值（审阅项 C-14），
            于是图上会出现一根没有任何数据的 100 Ma 时标，看起来像有结果。
        """
        if not self._node_coords:
            raise ValueError(
                "DatingFigure: 尚未计算出任何节点坐标，无法确定时间轴范围。"
                "请先成功渲染（render()）后再绘制地质时标。"
            )

        ages: List[float] = []
        for coord in self._node_coords.values():
            age = getattr(coord, "age", None)
            if age is None:
                continue
            try:
                value = float(age)
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                ages.append(value)

        if not ages:
            raise ValueError(
                "DatingFigure: 节点坐标里没有任何有限年龄值，无法确定时间轴范围"
                f"（共 {len(self._node_coords)} 个节点）。"
            )

        min_age = min(ages)
        max_age = max(ages)
        if max_age - min_age <= 0:
            raise ValueError(
                "DatingFigure: 时间轴退化——所有节点的年龄都相同"
                f"（{max_age:.6g} Ma），无法画出有宽度的时标。"
                "这通常意味着根年龄为 0 或树上没有分支长度信息；"
                "请显式提供 root_age 或检查输入树是否为时间树。"
            )

        padding = (max_age - min_age) * 0.05
        return (max_age + padding, max(0, min_age - padding))

    def _synchronize_axes(self) -> None:
        """Synchronize the X-axis limits between tree and geo axes.

        Uses the padded age range from _get_age_range() directly,
        avoiding double padding.
        """
        if self._ax_tree is None or self._ax_geo is None:
            return

        xlim = self._get_age_range()
        if xlim[0] - xlim[1] <= 0:
            raise ValueError(
                "DatingFigure: 计算得到的时间轴范围退化 "
                f"{xlim}（上界应比下界老），拒绝出图。"
            )
        # _get_age_range 已包含 5% padding，此处不再额外放大
        self._ax_tree.set_xlim(xlim)
        self._ax_geo.set_xlim(xlim)

    def set_title(
        self,
        title: str,
        fontsize: int = 14,
        **kwargs: Any,
    ) -> "DatingFigure":
        """
        Set the figure title.

        Parameters
        ----------
        title : str
            The title text.
        fontsize : int
            Font size for the title.
        **kwargs
            Additional keyword arguments passed to suptitle.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._figure is None:
            raise RuntimeError("Please call render() first.")

        self._figure.suptitle(title, fontsize=fontsize, **kwargs)
        return self

    def set_xlabel(
        self,
        label: str,
        fontsize: int = 10,
        **kwargs: Any,
    ) -> "DatingFigure":
        """
        Set the X-axis label.

        Parameters
        ----------
        label : str
            The label text.
        fontsize : int
            Font size for the label.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._ax_tree is None:
            raise RuntimeError("Please call render() first.")

        self._ax_tree.set_xlabel(label, fontsize=fontsize, **kwargs)
        return self

    def customize_axes(
        self,
        tree_kwargs: Optional[Dict[str, Any]] = None,
        geo_kwargs: Optional[Dict[str, Any]] = None,
    ) -> "DatingFigure":
        """
        Customize axes properties.

        This method allows direct customization of the tree and geo axes
        after rendering.

        Parameters
        ----------
        tree_kwargs : Dict[str, Any], optional
            Keyword arguments for the tree axes.
        geo_kwargs : Dict[str, Any], optional
            Keyword arguments for the geo axes.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if tree_kwargs and self._ax_tree is not None:
            for key, value in tree_kwargs.items():
                if hasattr(self._ax_tree, key):
                    getattr(self._ax_tree, "set_" + key)(value)

        if geo_kwargs and self._ax_geo is not None:
            for key, value in geo_kwargs.items():
                if hasattr(self._ax_geo, key):
                    getattr(self._ax_geo, "set_" + key)(value)

        return self

    def save(
        self,
        output_path: Union[str, Path],
        dpi: Optional[int] = None,
        format: Optional[str] = None,
        **kwargs: Any,
    ) -> Path:
        """
        Save the figure to a file.

        Parameters
        ----------
        output_path : Union[str, Path]
            The output file path.
        dpi : int, optional
            Dots per inch for raster formats. Overrides config.
        format : str, optional
            Output format. If not provided, inferred from extension.
        **kwargs
            Additional keyword arguments passed to savefig.

        Returns
        -------
        Path
            The path to the saved file.
        """
        if self._figure is None:
            raise RuntimeError("Please call render() first.")

        if dpi is None:
            dpi = self._dpi

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format is None:
            format = output_path.suffix.lstrip(".").lower()

        self._figure.savefig(output_path, format=format, dpi=dpi, **kwargs)
        logger.info(f"Figure saved to {output_path}")

        return output_path

    def show(self) -> None:
        """Display the figure."""
        if self._figure is None:
            raise RuntimeError("Please call render() first.")

        plt.show()

    def close(self) -> None:
        """Close the figure to free memory."""
        if self._figure is not None:
            plt.close(self._figure)
            self._figure = None
            self._figure_created = False
            self._ax_tree = None
            self._ax_geo = None

    def get_node_coordinates(self) -> Dict[Any, Any]:
        """
        Get the computed node coordinates.

        Returns
        -------
        Dict[Any, Any]
            Dictionary mapping nodes to their coordinates.
        """
        if self._node_coords is None:
            raise RuntimeError("Please call render() first to compute coordinates.")
        return self._node_coords

    def highlight_clade(
        self,
        node_name: str,
        color: str = "red",
        alpha: float = 0.3,
        **kwargs: Any,
    ) -> "DatingFigure":
        """
        Highlight a clade on the tree.

        Parameters
        ----------
        node_name : str
            Name of the node to highlight.
        color : str
            Color for the highlight.
        alpha : float
            Transparency of the highlight.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._ax_tree is None or self._node_coords is None or self._tree is None:
            # 三个都是 render() 才建起来的状态；_tree 尤其不能缺，
            # 下面的 traverse() 是逐节点定位的唯一手段。
            raise RuntimeError("Please call render() first.")

        target_node = None
        for node in self._tree.traverse():
            if node.name == node_name:
                target_node = node
                break

        if target_node is None:
            logger.warning(f"Node '{node_name}' not found")
            return self

        descendant_coords = []
        for node in target_node.traverse():
            if node in self._node_coords:
                coord = self._node_coords[node]
                descendant_coords.append((coord.x, coord.y))

        if not descendant_coords:
            return self

        x_vals = [c[0] for c in descendant_coords]
        y_vals = [c[1] for c in descendant_coords]

        x_min, x_max = min(x_vals), max(x_vals)
        y_min, y_max = min(y_vals), max(y_vals)

        from matplotlib.patches import Rectangle

        rect = Rectangle(
            (x_min, y_min),
            x_max - x_min,
            y_max - y_min,
            facecolor=color,
            alpha=alpha,
            zorder=0,
            **kwargs,
        )
        self._ax_tree.add_patch(rect)

        return self

    def collapse_clade(
        self,
        node_name: str,
        collapse_to: str = "triangle",
    ) -> "DatingFigure":
        """
        Mark a clade for collapse during rendering.

        Note: This modifies the tree and requires re-rendering.

        Parameters
        ----------
        node_name : str
            Name of the node to collapse.
        collapse_to : str
            Style of collapse representation ("triangle", "wedge").

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._tree is None:
            raise RuntimeError("No tree has been added.")

        for node in self._tree.traverse():
            if node.name == node_name:
                self._set_node_attr(node, "collapsed", True)
                self._set_node_attr(node, "collapse_style", collapse_to)
                logger.debug(f"Marked node '{node_name}' for collapse")
                break
        else:
            logger.warning(f"Node '{node_name}' not found")

        return self

    def rotate_clade(
        self,
        node_name: str,
    ) -> "DatingFigure":
        """
        Rotate a clade (swap children order).

        Note: This modifies the tree and requires re-rendering.

        Parameters
        ----------
        node_name : str
            Name of the node whose children to rotate.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._tree is None:
            raise RuntimeError("No tree has been added.")

        for node in self._tree.traverse():
            if node.name == node_name and len(node.children) == 2:
                node.children = node.children[::-1]
                logger.debug(f"Rotated clade at node '{node_name}'")
                break
        else:
            logger.warning(f"Node '{node_name}' not found or has != 2 children")

        return self

    def prune(
        self,
        taxa: List[str],
        preserve_branch_lengths: bool = True,
    ) -> "DatingFigure":
        """
        Prune the tree to keep only specified taxa.

        Note: This modifies the tree and requires re-rendering.

        Parameters
        ----------
        taxa : List[str]
            List of taxon names to keep.
        preserve_branch_lengths : bool
            Whether to preserve branch lengths during pruning.

        Returns
        -------
        DatingFigure
            Self for method chaining.
        """
        if self._tree is None:
            raise RuntimeError("No tree has been added.")

        try:
            self._tree.prune(taxa, preserve_branch_length=preserve_branch_lengths)
            logger.debug(f"Pruned tree to {len(taxa)} taxa")
        except Exception as e:
            logger.error(f"Failed to prune tree: {e}")
            raise RuntimeError(f"Failed to prune tree: {e}") from e

        return self

    def get_tree_dataframe(self) -> pd.DataFrame:
        """
        Get the tree data as a pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing node information.
        """
        if self._tree_dataframe is None:
            if self._tree is None:
                raise RuntimeError("No tree has been added.")
            self._tree_dataframe = self._tree_to_dataframe(self._tree)
        return self._tree_dataframe

    def _tree_to_dataframe(self, tree: Any) -> pd.DataFrame:
        """
        Convert tree to a pandas DataFrame.

        Parameters
        ----------
        tree : Any
            The tree to convert.

        Returns
        -------
        pd.DataFrame
            DataFrame with node information.
        """
        records = []
        for node in tree.traverse():
            is_leaf_node = node.is_leaf() if callable(node.is_leaf) else node.is_leaf
            is_root_node = node.is_root() if callable(node.is_root) else node.is_root

            record = {
                "name": node.name or "",
                "is_leaf": is_leaf_node,
                "is_root": is_root_node,
                "age": self._get_node_attr(node, "age"),
                "ci_lower": self._get_node_attr(node, "ci_lower"),
                "ci_upper": self._get_node_attr(node, "ci_upper"),
            }
            # 获取额外属性（兼容 ete3/ete4）
            if hasattr(node, "props"):
                # ete4
                for key, value in node.props.items():
                    if key not in record:
                        record[key] = value
            elif hasattr(node, "features"):
                # ete3: features 是 set
                for key in node.features:
                    if key not in record:
                        record[key] = getattr(node, key, None)
            records.append(record)
        return pd.DataFrame(records)

    def __enter__(self) -> "DatingFigure":
        """Context manager entry."""
        return self

    def __exit__(
        self,
        exc_type: Optional[type],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        """Context manager exit - close the figure."""
        self.close()

    def __repr__(self) -> str:
        """String representation."""
        status = "rendered" if self._figure_created else "not rendered"
        tree_status = "loaded" if self._tree is not None else "no tree"
        geo_status = "enabled" if self._show_geo_scale else "disabled"
        return f"DatingFigure({tree_status}, geo={geo_status}, {status})"
