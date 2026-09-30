"""
Plotter Backends - 可视化后端抽象

提供不同的可视化后端实现，支持：
- MatplotlibBackend: 基于 matplotlib 的后端（默认）
- PlotlyBackend: 基于 Plotly 的后端（未来扩展）
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt

from .tree_plot import NodeCoordinate

if TYPE_CHECKING:
    pass


class PlotterBackend(ABC):
    """
    可视化后端抽象基类

    定义所有可视化后端必须实现的接口
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """后端名称"""
        pass

    @abstractmethod
    def initialize_figure(self, width: float, height: float, dpi: int) -> Any:
        """初始化图表对象"""
        pass

    @abstractmethod
    def create_axes(
        self, figure: Any, layout: str = "single", **kwargs: Any
    ) -> Union[Any, Tuple[Any, Any]]:
        """创建坐标轴

        Args:
            figure: 图表对象
            layout: 布局类型 ("single", "tree_with_geo")
            **kwargs: 布局参数

        Returns:
            单坐标轴或坐标轴元组 (tree_ax, geo_ax)
        """
        pass

    @abstractmethod
    def plot_tree(
        self, ax: Any, tree: Any, coords: Dict[Any, NodeCoordinate], **kwargs: Any
    ) -> None:
        """绘制系统发育树"""
        pass

    @abstractmethod
    def plot_geo_timescale(
        self, ax: Any, age_range: Tuple[float, float], geo_data: Any, **kwargs: Any
    ) -> None:
        """绘制地质时间尺度"""
        pass

    @abstractmethod
    def synchronize_axes(self, ax1: Any, ax2: Any, **kwargs: Any) -> None:
        """同步两个坐标轴的显示范围"""
        pass

    @abstractmethod
    def set_title(self, figure: Any, title: str, **kwargs: Any) -> None:
        """设置图表标题"""
        pass

    @abstractmethod
    def set_xlabel(self, ax: Any, label: str, **kwargs: Any) -> None:
        """设置 X 轴标签"""
        pass

    @abstractmethod
    def save_figure(
        self, figure: Any, output_path: Union[str, Path], **kwargs: Any
    ) -> Path:
        """保存图表到文件"""
        pass

    @abstractmethod
    def close_figure(self, figure: Any) -> None:
        """关闭图表释放内存"""
        pass


class MatplotlibBackend(PlotterBackend):
    """
    Matplotlib 可视化后端

    基于 matplotlib 的标准后端实现
    """

    @property
    def name(self) -> str:
        return "matplotlib"

    def initialize_figure(self, width: float, height: float, dpi: int) -> Any:
        import matplotlib.pyplot as plt

        return plt.figure(figsize=(width, height), dpi=dpi)

    def create_axes(
        self, figure: Any, layout: str = "single", **kwargs: Any
    ) -> Union[Any, Tuple[Any, Any]]:
        from matplotlib.gridspec import GridSpec

        if layout == "tree_with_geo":
            height_ratio = kwargs.get("height_ratio", 0.15)
            gs = GridSpec(
                2, 1, height_ratios=[1 - height_ratio, height_ratio], hspace=0.0
            )
            tree_ax = figure.add_subplot(gs[0])
            geo_ax = figure.add_subplot(gs[1])
            return tree_ax, geo_ax
        else:
            return figure.add_subplot(111), None

    def plot_tree(
        self, ax: Any, tree: Any, coords: Dict[Any, NodeCoordinate], **kwargs: Any
    ) -> None:
        from matplotlib.collections import LineCollection

        branch_color = kwargs.get("branch_color", "#333333")
        leaf_color = kwargs.get("leaf_color", "#333333")
        leaf_size = kwargs.get("leaf_size", 5.0)
        show_leaf_labels = kwargs.get("show_leaf_labels", True)
        label_offset = kwargs.get("label_offset", 1.0)
        label_fontsize = kwargs.get("label_fontsize", 9)

        segments = []
        properties = []

        for node in tree.traverse():
            # ete3/ete4 兼容性
            is_root_node = node.is_root() if callable(node.is_root) else node.is_root
            is_leaf_node = node.is_leaf() if callable(node.is_leaf) else node.is_leaf

            if is_root_node:
                continue

            node_coord = coords.get(node)
            parent_coord = coords.get(node.up)

            if node_coord is None or parent_coord is None:
                continue

            x_child = node_coord.x
            y_child = node_coord.y
            x_parent = parent_coord.x

            horizontal_segment = [(x_parent, y_child), (x_child, y_child)]
            segments.append(horizontal_segment)
            properties.append(
                {
                    "is_horizontal": True,
                    "node_name": node_coord.name,
                    "is_leaf": is_leaf_node,
                    **node_coord.props,
                }
            )

        for node in tree.traverse():
            is_leaf_node = node.is_leaf() if callable(node.is_leaf) else node.is_leaf
            if is_leaf_node:
                continue

            node_coord = coords.get(node)
            if node_coord is None:
                continue

            child_ys = []
            for child in node.children:
                child_coord = coords.get(child)
                if child_coord is not None:
                    child_ys.append(child_coord.y)

            if len(child_ys) >= 1:
                x_parent = node_coord.x
                vertical_segment = [
                    (x_parent, min(child_ys)),
                    (x_parent, max(child_ys)),
                ]
                segments.append(vertical_segment)
                properties.append(
                    {
                        "is_horizontal": False,
                        "node_name": node_coord.name,
                        **node_coord.props,
                    }
                )

        if segments:
            colors = self._compute_colors(segments, properties, kwargs)
            linewidths = self._compute_linewidths(segments, properties, kwargs)

            lc = LineCollection(
                segments, colors=colors, linewidths=linewidths, zorder=1
            )
            ax.add_collection(lc)

        leaf_x, leaf_y = [], []
        internal_x, internal_y = [], []
        root_x, root_y = [], []

        for node, coord in coords.items():
            if coord.is_root:
                root_x.append(coord.x)
                root_y.append(coord.y)
            elif coord.is_leaf:
                leaf_x.append(coord.x)
                leaf_y.append(coord.y)
            else:
                internal_x.append(coord.x)
                internal_y.append(coord.y)

        if leaf_x and leaf_size > 0:
            ax.scatter(leaf_x, leaf_y, s=leaf_size**2, c=leaf_color, zorder=2)

        if internal_x and kwargs.get("internal_node_size", 0) > 0:
            ax.scatter(
                internal_x,
                internal_y,
                s=kwargs.get("internal_node_size", 0) ** 2,
                c=branch_color,
                zorder=2,
            )

        if root_x:
            ax.scatter(
                root_x, root_y, s=(leaf_size * 1.5) ** 2, c=branch_color, zorder=2
            )

        if show_leaf_labels:
            for node, coord in coords.items():
                if coord.is_leaf and coord.name:
                    ax.text(
                        coord.x - label_offset,
                        coord.y,
                        coord.name,
                        fontsize=label_fontsize,
                        va="center",
                        ha="left",
                    )

    def _compute_colors(
        self, segments: List, properties: List[Dict], kwargs: Dict
    ) -> List:
        import matplotlib.colors as mcolors

        default_color = kwargs.get("branch_color", "#333333")
        color_column = kwargs.get("color_column")

        if color_column is None:
            return [default_color] * len(segments)

        values = [props.get(color_column) for props in properties]
        unique_values = set(v for v in values if v is not None)

        if all(isinstance(v, (int, float)) for v in unique_values if v is not None):
            numeric_values = [v if v is not None else 0 for v in values]
            min_val, max_val = min(numeric_values), max(numeric_values)

            if max_val == min_val:
                normalized = [0.5] * len(numeric_values)
            else:
                normalized = [
                    (v - min_val) / (max_val - min_val) for v in numeric_values
                ]

            cmap = kwargs.get("cmap", "viridis")
            colormap = plt.colormaps.get(cmap, plt.colormaps["viridis"])
            return [colormap(n) for n in normalized]
        else:
            palette = kwargs.get("palette") or list(mcolors.TABLEAU_COLORS.values())
            value_to_color = {
                val: palette[i % len(palette)]
                for i, val in enumerate(sorted(unique_values))
            }
            return [value_to_color.get(v, default_color) for v in values]

    def _compute_linewidths(
        self, segments: List, properties: List[Dict], kwargs: Dict
    ) -> List:
        default_linewidth = kwargs.get("branch_linewidth", 1.5)
        size_column = kwargs.get("size_column")

        if size_column is None:
            return [default_linewidth] * len(segments)

        # 逐位把可调列转成 float：缺位当成 0.0（与旧行为一致）。
        values: List[float] = []
        for props in properties:
            raw = props.get(size_column)
            values.append(0.0 if raw is None else float(raw))
        min_val: float = min(values) if values else 0.0
        max_val: float = max(values) if values else 1.0
        size_range: Tuple[float, float] = kwargs.get("size_range", (2.0, 8.0))

        if max_val == min_val:
            normalized = [0.5] * len(values)
        else:
            normalized = [(v - min_val) / (max_val - min_val) for v in values]

        return [size_range[0] + n * (size_range[1] - size_range[0]) for n in normalized]

    def plot_geo_timescale(
        self, ax: Any, age_range: Tuple[float, float], geo_data: Any, **kwargs: Any
    ) -> None:
        import pandas as pd
        from matplotlib.patches import Rectangle

        show_labels = kwargs.get("show_labels", True)
        label_rotation = kwargs.get("label_rotation", 90)
        show_age_ticks = kwargs.get("show_age_ticks", True)
        age_tick_interval = kwargs.get("age_tick_interval", 10.0)
        split_rows = kwargs.get("split_rows", True)

        if geo_data is None or (isinstance(geo_data, pd.DataFrame) and geo_data.empty):
            ax.set_visible(False)
            return

        if isinstance(geo_data, (str, Path)):
            geo_data = pd.read_csv(geo_data)

        min_age, max_age = age_range
        # 本项目的时标区间既可能是 (新界, 老界) 也可能是 DatingFigure 给出的
        # 降序 (老界, 新界)：统一成升序，否则区间筛选与刻度都会算错。
        if min_age > max_age:
            min_age, max_age = max_age, min_age

        ax.set_xlim(max_age + 5, max(0, min_age - 5))
        ax.set_ylim(0, 1)

        for _, row in geo_data.iterrows():
            # 兼容 ICS 风格 (start/end) 与 age 风格 (max_age/min_age)
            age_top: Optional[float] = None
            age_bottom: Optional[float] = None
            if "max_age" in row and "min_age" in row and pd.notna(row.get("max_age")):
                age_top = float(row["max_age"])
                age_bottom = float(row["min_age"])
            elif "start" in row and "end" in row:
                a, b = float(row["start"]), float(row["end"])
                age_top, age_bottom = max(a, b), min(a, b)  # 统一为 (old, young)
            elif "age" in row:
                age_top = age_bottom = float(row["age"])
            # 三个分支要么把两个都设上、要么都不设，所以任一为 None 就意味着
            # 这一行没有可用的年龄，整行跳过（旧写法只查 age_top，于是
            # age_bottom 在静态视角下还是 Optional，后面的算术与比较全部落空）。
            if age_top is None or age_bottom is None:
                continue

            name = row.get("name", row.get("interval", ""))

            if age_top < min_age or age_bottom > max_age:
                continue

            y_base = 0.3 if split_rows else 0.2
            height = 0.4 if split_rows else 0.6

            color = row.get("color")
            if pd.notna(color):
                rect = Rectangle(
                    (age_top, y_base),
                    age_bottom - age_top,
                    height,
                    facecolor=color,
                    edgecolor="none",
                    alpha=0.7,
                )
                ax.add_patch(rect)

            if show_labels and name:
                label_y = y_base + height / 2
                ax.text(
                    (age_top + age_bottom) / 2,
                    label_y,
                    name,
                    ha="center",
                    va="center",
                    fontsize=7,
                    rotation=label_rotation,
                )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(left=False, labelleft=False, bottom=False, top=False)

        if show_age_ticks:
            ticks = []
            tick_labels = []
            tick = max_age - (max_age % age_tick_interval)
            while tick > min_age:
                ticks.append(tick)
                tick_labels.append(str(int(tick)))
                tick -= age_tick_interval

            ax.set_xticks(ticks)
            ax.set_xticklabels(tick_labels)
            ax.tick_params(bottom=True, labelbottom=True)

    def synchronize_axes(self, ax1: Any, ax2: Any, **kwargs: Any) -> None:
        padding_top = kwargs.get("padding_top", 5.0)
        padding_bottom = kwargs.get("padding_bottom", 0.0)

        xlim1 = ax1.get_xlim()
        xlim2 = ax2.get_xlim()

        new_xlim = (
            max(xlim1[0], xlim2[0]) + padding_top,
            min(xlim1[1], xlim2[1]) - padding_bottom,
        )

        ax1.set_xlim(new_xlim)
        ax2.set_xlim(new_xlim)

    def set_title(self, figure: Any, title: str, **kwargs: Any) -> None:
        fontsize = kwargs.get("fontsize", 14)
        figure.suptitle(title, fontsize=fontsize)

    def set_xlabel(self, ax: Any, label: str, **kwargs: Any) -> None:
        fontsize = kwargs.get("fontsize", 10)
        ax.set_xlabel(label, fontsize=fontsize)

    def save_figure(
        self, figure: Any, output_path: Union[str, Path], **kwargs: Any
    ) -> Path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        dpi = kwargs.get("dpi")
        format = kwargs.get("format")

        if dpi is None:
            dpi = 300
        if format is None:
            format = output_path.suffix.lstrip(".").lower()

        figure.savefig(output_path, format=format, dpi=dpi)
        return output_path

    def close_figure(self, figure: Any) -> None:
        import matplotlib.pyplot as plt

        plt.close(figure)


class PlotlyBackend(PlotterBackend):
    """
    Plotly 可视化后端（按需交互式后端）

    由于 plotly 不是核心依赖，本后端仅在运行时安装 plotly 后才可用。
    未安装时初始化会给出清晰的安装提示。
    """

    @property
    def name(self) -> str:
        return "plotly"

    def _require_plotly(self) -> Any:
        try:
            import plotly.graph_objects as go

            return go
        except ImportError as e:
            raise ImportError(
                "Plotly backend requires the 'plotly' package. "
                "Install it with: pip install plotly"
            ) from e

    def initialize_figure(self, width: float, height: float, dpi: int) -> Any:
        go = self._require_plotly()
        # Plotly 使用像素尺寸，将英寸按默认 96 dpi 换算
        px_width = width * 96
        px_height = height * 96
        return go.Figure(layout=dict(width=px_width, height=px_height))

    def create_axes(
        self, figure: Any, layout: str = "single", **kwargs: Any
    ) -> Union[Any, Tuple[Any, Any]]:
        self._require_plotly()
        # Plotly Figure 本身即画布，用 secondary_y 或 subplots 分隔；
        # 这里返回 figure 与 geo 占用的 row 索引信息，由调用方通过 make_subplots 处理。
        if layout == "tree_with_geo":
            try:
                from plotly.subplots import make_subplots
            except ImportError as e:
                raise ImportError(
                    "Plotly backend requires the 'plotly' package. "
                    "Install it with: pip install plotly"
                ) from e
            height_ratio = kwargs.get("height_ratio", 0.15)
            fig = make_subplots(
                rows=2,
                cols=1,
                row_heights=[1 - height_ratio, height_ratio],
                vertical_spacing=0.0,
                shared_xaxes=True,
            )
            return fig, (1, 2)
        return figure, None

    def plot_tree(
        self, ax: Any, tree: Any, coords: Dict[Any, NodeCoordinate], **kwargs: Any
    ) -> None:
        go = self._require_plotly()
        if not isinstance(ax, go.Figure):
            raise TypeError("PlotlyBackend.plot_tree expects a plotly Figure")

        branch_color = kwargs.get("branch_color", "#333333")
        branch_linewidth = kwargs.get("branch_linewidth", 1.5)
        leaf_color = kwargs.get("leaf_color", "#333333")
        leaf_size = kwargs.get("leaf_size", 5.0)
        show_leaf_labels = kwargs.get("show_leaf_labels", True)
        label_fontsize = kwargs.get("label_fontsize", 9)
        color_column = kwargs.get("color_column")
        cmap = kwargs.get("cmap", "viridis")

        # 构建分支线段
        x_lines, y_lines = [], []
        for node in tree.traverse():
            is_root_node = node.is_root() if callable(node.is_root) else node.is_root
            is_leaf_node = node.is_leaf() if callable(node.is_leaf) else node.is_leaf
            if is_root_node:
                continue
            node_coord = coords.get(node)
            parent_coord = coords.get(node.up)
            if node_coord is None or parent_coord is None:
                continue
            x_lines.extend([parent_coord.x, node_coord.x, None])
            y_lines.extend([node_coord.y, node_coord.y, None])
            if not is_leaf_node:
                child_ys = [coords[c].y for c in node.children if c in coords]
                if child_ys:
                    x_lines.extend([node_coord.x, node_coord.x, None])
                    y_lines.extend([min(child_ys), max(child_ys), None])

        ax.add_trace(
            go.Scatter(
                x=x_lines,
                y=y_lines,
                mode="lines",
                line=dict(color=branch_color, width=branch_linewidth),
                hoverinfo="skip",
                showlegend=False,
            )
        )

        # 节点散点
        leaf_x, leaf_y, leaf_text = [], [], []
        internal_x, internal_y = [], []
        root_x, root_y = [], []
        for node, coord in coords.items():
            if coord.is_root:
                root_x.append(coord.x)
                root_y.append(coord.y)
            elif coord.is_leaf:
                leaf_x.append(coord.x)
                leaf_y.append(coord.y)
                leaf_text.append(coord.name or "")
            else:
                internal_x.append(coord.x)
                internal_y.append(coord.y)

        if leaf_x and leaf_size > 0:
            marker = dict(size=leaf_size, color=leaf_color)
            if color_column and any(
                props.get(color_column)
                for props in [getattr(c, "props", {}) for c in coords.values()]
            ):
                # 颜色映射逻辑较复杂，此处仅支持连续数值列
                values = [
                    getattr(coords[n], "props", {}).get(color_column) for n in coords
                ]
                marker["color"] = values
                marker["colorscale"] = cmap
                marker["showscale"] = True
            ax.add_trace(
                go.Scatter(
                    x=leaf_x,
                    y=leaf_y,
                    mode="markers+text" if show_leaf_labels else "markers",
                    text=leaf_text if show_leaf_labels else None,
                    textposition="middle left",
                    textfont=dict(size=label_fontsize),
                    marker=marker,
                    hovertemplate="%{text}<br>Age: %{x:.3f}<extra></extra>",
                    showlegend=False,
                )
            )

        if internal_x and kwargs.get("internal_node_size", 0) > 0:
            ax.add_trace(
                go.Scatter(
                    x=internal_x,
                    y=internal_y,
                    mode="markers",
                    marker=dict(
                        size=kwargs.get("internal_node_size", 0), color=branch_color
                    ),
                    hoverinfo="skip",
                    showlegend=False,
                )
            )

        if root_x:
            ax.add_trace(
                go.Scatter(
                    x=root_x,
                    y=root_y,
                    mode="markers",
                    marker=dict(size=leaf_size * 1.5, color=branch_color),
                    hoverinfo="skip",
                    showlegend=False,
                )
            )

    def plot_geo_timescale(
        self, ax: Any, age_range: Tuple[float, float], geo_data: Any, **kwargs: Any
    ) -> None:
        go = self._require_plotly()
        if not isinstance(ax, go.Figure):
            raise TypeError("PlotlyBackend.plot_geo_timescale expects a plotly Figure")
        if geo_data is None or (hasattr(geo_data, "empty") and geo_data.empty):
            return
        import pandas as pd

        if isinstance(geo_data, (str, Path)):
            geo_data = pd.read_csv(geo_data)
        min_age, max_age = age_range
        if min_age > max_age:  # 同 MatplotlibBackend：统一成 (新界, 老界)
            min_age, max_age = max_age, min_age
        for _, row in geo_data.iterrows():
            # 兼容 ICS 风格 (start/end) 与 age 风格 (max_age/min_age)
            age_top: Optional[float] = None
            age_bottom: Optional[float] = None
            if "max_age" in row and "min_age" in row and pd.notna(row.get("max_age")):
                age_top = float(row["max_age"])
                age_bottom = float(row["min_age"])
            elif "start" in row and "end" in row:
                a, b = float(row["start"]), float(row["end"])
                age_top, age_bottom = max(a, b), min(a, b)  # 统一为 (old, young)
            elif "age" in row:
                age_top = age_bottom = float(row["age"])
            # 同上：两个年边界必须一起检查，否则 age_bottom 仍是 Optional。
            if age_top is None or age_bottom is None:
                continue
            if age_top < min_age or age_bottom > max_age:
                continue
            color = row.get("color")
            fillcolor = color if pd.notna(color) else "#cccccc"
            ax.add_trace(
                go.Scatter(
                    x=[age_top, age_bottom, age_bottom, age_top, age_top],
                    y=[0, 0, 1, 1, 0],
                    fill="toself",
                    fillcolor=fillcolor,
                    line=dict(width=0),
                    hoverinfo="text",
                    text=row.get("name", row.get("interval", "")),
                    showlegend=False,
                )
            )

    def synchronize_axes(self, ax1: Any, ax2: Any, **kwargs: Any) -> None:
        go = self._require_plotly()
        if not (isinstance(ax1, go.Figure) and isinstance(ax2, go.Figure)):
            return
        # Plotly 通过 shared_xaxes 在 create_axes 阶段已同步

    def set_title(self, figure: Any, title: str, **kwargs: Any) -> None:
        self._require_plotly()
        figure.update_layout(title_text=title)

    def set_xlabel(self, ax: Any, label: str, **kwargs: Any) -> None:
        go = self._require_plotly()
        if isinstance(ax, go.Figure):
            ax.update_layout(xaxis_title=label)

    def save_figure(
        self, figure: Any, output_path: Union[str, Path], **kwargs: Any
    ) -> Path:
        go = self._require_plotly()
        if not isinstance(figure, go.Figure):
            raise TypeError("Expected plotly Figure")
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        # Plotly 不支持 dpi 参数，因此忽略
        figure.write_image(str(output_path), scale=kwargs.get("scale", 3))
        return output_path

    def close_figure(self, figure: Any) -> None:
        # Plotly Figure 无需显式关闭，手动清理引用即可
        del figure


_BACKENDS = {
    "matplotlib": MatplotlibBackend,
    "plotly": PlotlyBackend,
}


def get_backend(name: str = "matplotlib") -> PlotterBackend:
    """获取指定名称的可视化后端"""
    if name not in _BACKENDS:
        raise ValueError(
            f"Unknown backend: {name}. Available: {list(_BACKENDS.keys())}"
        )
    return _BACKENDS[name]()


def list_backends() -> List[str]:
    """列出所有可用的可视化后端"""
    return list(_BACKENDS.keys())
