"""
Geological timescale plotting module for PhyloDater visualization.

Based on ChronoPhylo, provides functions for plotting ICS chronostratigraphic timescales
alongside phylogenetic trees.
"""

import importlib.resources
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.axes as mpl_axes
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from phylodater.infrastructure.logging import get_logger
from phylodater.viz import data as data_package

logger = get_logger()


@dataclass
class GeoInterval:
    """Data class for a geological time interval."""

    name: str
    start: float
    end: float
    color: str
    rank: str
    parent: Optional[str] = None


def load_ics_data(custom_path: Optional[Union[str, Path]] = None) -> pd.DataFrame:
    """
    Load ICS chronostratigraphic timescale data.

    Uses importlib.resources to ensure the data file is found after packaging.

    Parameters
    ----------
    custom_path : Union[str, Path], optional
        Path to a custom timescale CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the timescale data.
    """
    if custom_path is not None:
        path = Path(custom_path)
        logger.info(f"Loading custom timescale from {path}")
        return pd.read_csv(path)

    try:
        data_path = importlib.resources.files(data_package).joinpath(
            "ics_timescale.csv"
        )
        with importlib.resources.as_file(data_path) as f:
            df = pd.read_csv(f)
        logger.debug(f"Loaded ICS timescale with {len(df)} intervals")
        return df
    except Exception as e:
        logger.error(f"Failed to load ICS timescale: {e}")
        raise


def get_default_geo_colors() -> Dict[str, str]:
    """
    Get default colors for geological periods.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping period names to hex colors.
    """
    return {
        "Phanerozoic": "#9AD9DD",
        "Proterozoic": "#F73563",
        "Archean": "#F99BC1",
        "Hadean": "#6A0DAD",
        "Cenozoic": "#F2F91D",
        "Mesozoic": "#67C5CA",
        "Paleozoic": "#99C08D",
        "Neoproterozoic": "#FEB342",
        "Mesoproterozoic": "#FDB462",
        "Paleoproterozoic": "#F74370",
        "Quaternary": "#F9F97F",
        "Neogene": "#FFE619",
        "Paleogene": "#FD9A52",
        "Cretaceous": "#7FC64E",
        "Jurassic": "#34B2C9",
        "Triassic": "#812B92",
        "Permian": "#F04028",
        "Carboniferous": "#67A599",
        "Devonian": "#CB8C37",
        "Silurian": "#B3E1B6",
        "Ordovician": "#009270",
        "Cambrian": "#7FA056",
        "Ediacaran": "#FED96A",
        "Cryogenian": "#FECC5C",
        "Tonian": "#FEBF4E",
        "Holocene": "#FEEBD2",
        "Pleistocene": "#FFEFAF",
        "Pliocene": "#FFFF99",
        "Miocene": "#FFFF00",
        "Oligocene": "#FEC07A",
        "Eocene": "#FDB46C",
        "Paleocene": "#FDA75F",
        "Neoarchean": "#F99BC1",
        "Mesoarchean": "#F768A9",
        "Paleoarchean": "#F4449F",
        "Eoarchean": "#DC148C",
        "Stenian": "#FED99A",
        "Ectasian": "#FDCC8A",
        "Calymmian": "#FDC07A",
        "Statherian": "#F875A7",
        "Orosirian": "#F76898",
        "Rhyacian": "#F75B89",
        "Siderian": "#F74F7C",
        "Upper Cretaceous": "#A6D84A",
        "Lower Cretaceous": "#8CCD57",
        "Upper Jurassic": "#B3E3EE",
        "Middle Jurassic": "#80CFD8",
        "Lower Jurassic": "#42AED0",
        "Upper Triassic": "#BD8CC3",
        "Middle Triassic": "#B168B1",
        "Lower Triassic": "#983999",
        "Lopingian": "#FCC0B2",
        "Guadalupian": "#E38776",
        "Cisuralian": "#F04028",
        "Pennsylvanian": "#7EBCC6",
        "Mississippian": "#678F66",
        "Upper Devonian": "#F2EDAD",
        "Middle Devonian": "#F1D576",
        "Lower Devonian": "#E5B75A",
        "Pridoli": "#E6F5E1",
        "Ludlow": "#D9F0DF",
        "Wenlock": "#CCEBD1",
        "Llandovery": "#A6DCB5",
        "Upper Ordovician": "#A6DBAB",
        "Middle Ordovician": "#74C69C",
        "Lower Ordovician": "#33A97E",
        "Furongian": "#D9F0BB",
        "Miaolingian": "#B3D492",
        "Terreneuvian": "#99B575",
    }


def get_geo_abbreviations() -> Dict[str, str]:
    """
    Get abbreviations for geological periods and epochs.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping period/epoch names to abbreviations.
    """
    return {
        "Quaternary": "Q",
        "Neogene": "Ng",
        "Paleogene": "Pg",
        "Cretaceous": "K",
        "Jurassic": "J",
        "Triassic": "Tr",
        "Permian": "P",
        "Carboniferous": "C",
        "Devonian": "D",
        "Silurian": "S",
        "Ordovician": "O",
        "Cambrian": "Cm",
        "Ediacaran": "Ed",
        "Cryogenian": "Cr",
        "Tonian": "To",
        "Stenian": "Sn",
        "Ectasian": "Ec",
        "Calymmian": "Ca",
        "Statherian": "St",
        "Orosirian": "Or",
        "Rhyacian": "Rh",
        "Siderian": "Si",
        "Neoarchean": "NA",
        "Mesoarchean": "MA",
        "Paleoarchean": "PA",
        "Eoarchean": "EA",
        "Hadean": "H",
        "Cenozoic": "Cz",
        "Mesozoic": "Mz",
        "Paleozoic": "Pz",
        "Neoproterozoic": "NP",
        "Mesoproterozoic": "MP",
        "Paleoproterozoic": "PP",
        "Phanerozoic": "Ph",
        "Proterozoic": "Pr",
        "Archean": "Ar",
        "Precambrian": "Pc",
        "Holocene": "Ho",
        "Pleistocene": "Pl",
        "Pliocene": "Plc",
        "Miocene": "Mi",
        "Oligocene": "Ol",
        "Eocene": "Eo",
        "Paleocene": "Pa",
        "Upper Cretaceous": "Ku",
        "Lower Cretaceous": "Kl",
        "Upper Jurassic": "Ju",
        "Middle Jurassic": "Jm",
        "Lower Jurassic": "Jl",
        "Upper Triassic": "Tru",
        "Middle Triassic": "Trm",
        "Lower Triassic": "Trl",
        "Lopingian": "Lo",
        "Guadalupian": "Gu",
        "Cisuralian": "Ci",
        "Pennsylvanian": "Pn",
        "Mississippian": "Ms",
        "Upper Devonian": "Du",
        "Middle Devonian": "Dm",
        "Lower Devonian": "Dl",
        "Pridoli": "Pri",
        "Ludlow": "Lu",
        "Wenlock": "We",
        "Llandovery": "Ll",
        "Upper Ordovician": "Ou",
        "Middle Ordovician": "Om",
        "Lower Ordovician": "Olo",
        "Furongian": "Fu",
        "Miaolingian": "Mia",
        "Terreneuvian": "Te",
    }


class GeoPlotter:
    """
    Class for plotting geological timescales on matplotlib axes.

    This class handles the visualization of ICS chronostratigraphic data
    as colored bars with labels and age markers.
    """

    def __init__(
        self,
        show_labels: bool = True,
        label_rotation: int = 90,
        show_age_ticks: bool = True,
        age_tick_interval: float = 10.0,
        padding_top: float = 5.0,
        padding_bottom: float = 0.0,
        height_ratio: float = 0.15,
        split_rows: bool = False,
    ) -> None:
        """
        Initialize the GeoPlotter.

        Parameters
        ----------
        show_labels : bool
            Whether to show interval labels.
        label_rotation : int
            Rotation angle for labels.
        show_age_ticks : bool
            Whether to show age tick marks.
        age_tick_interval : float
            Interval between age ticks in Ma.
        padding_top : float
            Padding above the timescale in Ma.
        padding_bottom : float
            Padding below the timescale in Ma.
        height_ratio : float
            Height ratio of the geological axis relative to tree axis.
        split_rows : bool
            Whether to split into two rows (Era on top, Eon on bottom).
        """
        self.show_labels = show_labels
        self.label_rotation = label_rotation
        self.show_age_ticks = show_age_ticks
        self.age_tick_interval = age_tick_interval
        self.padding_top = padding_top
        self.padding_bottom = padding_bottom
        self.height_ratio = height_ratio
        self.split_rows = split_rows

        self._geo_data: Optional[pd.DataFrame] = None
        self._colors: Dict[str, str] = get_default_geo_colors()
        self._abbreviations: Dict[str, str] = get_geo_abbreviations()

    def load_data(
        self,
        custom_path: Optional[Union[str, Path]] = None,
    ) -> pd.DataFrame:
        """
        Load geological timescale data.

        Parameters
        ----------
        custom_path : Union[str, Path], optional
            Path to a custom timescale CSV file.

        Returns
        -------
        pd.DataFrame
            The loaded timescale data.
        """
        self._geo_data = load_ics_data(custom_path)
        return self._geo_data

    def plot(
        self,
        ax: mpl_axes.Axes,
        age_range: Tuple[float, float],
        geo_data: Optional[pd.DataFrame] = None,
        rank_filter: Optional[List[str]] = None,
    ) -> None:
        """
        Plot the geological timescale on the given axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        age_range : Tuple[float, float]
            The (min_age, max_age) range to display. 两个方向的传入都被接受：
            ``DatingFigure._get_age_range()`` 给出的是**降序** ``(老界, 新界)``
            （轴反向，老在左），这里统一换算成 (新界, 老界) 再使用。
        geo_data : pd.DataFrame, optional
            Timescale data. If not provided, uses loaded data.
        rank_filter : List[str], optional
            List of ranks to display (e.g., ["Period", "Epoch"]).
        """
        if geo_data is None:
            if self._geo_data is None:
                self._geo_data = load_ics_data()
            geo_data = self._geo_data

        if geo_data is None or len(geo_data) == 0:
            logger.warning("No geological data available for plotting")
            return

        # 统一成 (新界 min_age, 老界 max_age)：_filter_by_age / _plot_age_ticks /
        # _set_axes_style 都按升序口径实现，收到降序区间时（DatingFigure 的默认
        # 传参方式）此前会筛出错误的区间、画不出任何色块、也拿不到刻度。
        edge_a, edge_b = age_range
        min_age, max_age = (edge_a, edge_b) if edge_a <= edge_b else (edge_b, edge_a)

        filtered_data = self._filter_by_age(geo_data, min_age, max_age)

        if rank_filter:
            filtered_data = filtered_data[filtered_data["rank"].isin(rank_filter)]

        if self.split_rows:
            self._plot_intervals_split(ax, filtered_data, min_age, max_age)
        else:
            self._plot_intervals(ax, filtered_data, min_age, max_age)

        if self.show_age_ticks:
            self._plot_age_ticks(ax, min_age, max_age)

        self._set_axes_style(ax, min_age, max_age)

    def _filter_by_age(
        self,
        geo_data: pd.DataFrame,
        min_age: float,
        max_age: float,
    ) -> pd.DataFrame:
        """Filter intervals that overlap with the age range."""
        start_col = self._find_column(
            geo_data, ["start", "age_top", "top", "start_age"]
        )
        end_col = self._find_column(
            geo_data, ["end", "age_bottom", "bottom", "end_age"]
        )

        if start_col is None or end_col is None:
            logger.warning("Could not find age columns in geological data")
            return geo_data

        # 与区间两端哪一列更老无关：先归一化为 (新界, 老界) 再做重叠判定
        young = geo_data[[start_col, end_col]].min(axis=1)
        old = geo_data[[start_col, end_col]].max(axis=1)
        mask = (young <= max_age) & (old >= min_age)

        return geo_data[mask].copy()

    def _find_column(
        self,
        df: pd.DataFrame,
        candidates: List[str],
    ) -> Optional[str]:
        """Find a column from a list of candidate names."""
        for candidate in candidates:
            if candidate in df.columns:
                return candidate
        return None

    def _plot_intervals(
        self,
        ax: mpl_axes.Axes,
        geo_data: pd.DataFrame,
        min_age: float,
        max_age: float,
    ) -> None:
        """Plot colored intervals."""
        start_col = self._find_column(
            geo_data, ["start", "age_top", "top", "start_age"]
        )
        end_col = self._find_column(
            geo_data, ["end", "age_bottom", "bottom", "end_age"]
        )
        name_col = self._find_column(geo_data, ["name", "interval", "interval_name"])
        color_col = self._find_column(geo_data, ["color", "colour", "hex"])
        rank_col = self._find_column(geo_data, ["rank", "level", "type"])

        if start_col is None or end_col is None:
            return

        y_min, y_max = 0, 0.25

        for _, row in geo_data.iterrows():
            a = float(row[start_col])
            b = float(row[end_col])
            young, old = min(a, b), max(a, b)  # 统一为 (young, old)

            # 把区间裁剪到视窗内：矩形从"较新的一端"画到"较老的一端"。
            # （此前写成 max(old, min_age) / min(young, max_age)，两端取反，
            # plot_start >= plot_end 恒成立，整条时标一个色块都画不出来。）
            plot_start = max(young, min_age)
            plot_end = min(old, max_age)

            if plot_start >= plot_end:
                continue

            color = self._get_color(row, color_col)

            rect = patches.Rectangle(
                (plot_start, y_min),
                plot_end - plot_start,
                y_max - y_min,
                linewidth=0.5,
                edgecolor="black",
                facecolor=color,
            )
            ax.add_patch(rect)

            if self.show_labels and name_col:
                name = str(row[name_col])
                mid_age = (plot_start + plot_end) / 2
                width = plot_end - plot_start
                rank = str(row[rank_col]) if rank_col else ""

                if width > 1 and rank in ["Period", "Epoch"]:
                    label = self._abbreviations.get(name, name)
                    ax.text(
                        mid_age,
                        (y_min + y_max) / 2,
                        label,
                        rotation=0,
                        ha="center",
                        va="center",
                        fontsize=10,
                        clip_on=True,
                    )

    def _plot_intervals_split(
        self,
        ax: mpl_axes.Axes,
        geo_data: pd.DataFrame,
        min_age: float,
        max_age: float,
    ) -> None:
        """Plot colored intervals in two rows with intelligent rank selection."""
        start_col = self._find_column(
            geo_data, ["start", "age_top", "top", "start_age"]
        )
        end_col = self._find_column(
            geo_data, ["end", "age_bottom", "bottom", "end_age"]
        )
        name_col = self._find_column(geo_data, ["name", "interval", "interval_name"])
        color_col = self._find_column(geo_data, ["color", "colour", "hex"])

        if start_col is None or end_col is None:
            return

        if max_age <= 100:
            upper_rank = "Epoch"
            lower_rank = "Period"
        elif max_age <= 2500:
            upper_rank = "Period"
            lower_rank = "Era"
        else:
            upper_rank = "Era"
            lower_rank = "Eon"

        upper_data = geo_data[geo_data["rank"] == upper_rank].copy()
        lower_data = geo_data[geo_data["rank"] == lower_rank].copy()

        y_upper_min, y_upper_max = 0.25, 0.5
        y_lower_min, y_lower_max = 0, 0.25

        for _, row in upper_data.iterrows():
            start = float(row[start_col])
            end = float(row[end_col])
            young, old = min(start, end), max(start, end)

            plot_start = max(young, min_age)
            plot_end = min(old, max_age)

            if plot_start >= plot_end:
                continue

            color = self._get_color(row, color_col)

            rect = patches.Rectangle(
                (plot_start, y_upper_min),
                plot_end - plot_start,
                y_upper_max - y_upper_min,
                linewidth=0.5,
                edgecolor="black",
                facecolor=color,
            )
            ax.add_patch(rect)

            if self.show_labels and name_col:
                name = str(row[name_col])
                mid_age = (plot_start + plot_end) / 2
                width = plot_end - plot_start

                label_threshold = 30 if upper_rank == "Era" else 10
                if width > label_threshold:
                    label = self._abbreviations.get(name, name)
                    ax.text(
                        mid_age,
                        (y_upper_min + y_upper_max) / 2,
                        label,
                        rotation=0,
                        ha="center",
                        va="center",
                        fontsize=10,
                        clip_on=True,
                    )

        for _, row in lower_data.iterrows():
            start = float(row[start_col])
            end = float(row[end_col])
            young, old = min(start, end), max(start, end)

            plot_start = max(young, min_age)
            plot_end = min(old, max_age)

            if plot_start >= plot_end:
                continue

            color = self._get_color(row, color_col)

            rect = patches.Rectangle(
                (plot_start, y_lower_min),
                plot_end - plot_start,
                y_lower_max - y_lower_min,
                linewidth=0.5,
                edgecolor="black",
                facecolor=color,
            )
            ax.add_patch(rect)

            if self.show_labels and name_col:
                name = str(row[name_col])
                mid_age = (plot_start + plot_end) / 2
                width = plot_end - plot_start

                label_threshold = 50 if lower_rank == "Eon" else 20
                if width > label_threshold:
                    label = self._abbreviations.get(name, name)
                    ax.text(
                        mid_age,
                        (y_lower_min + y_lower_max) / 2,
                        label,
                        rotation=0,
                        ha="center",
                        va="center",
                        fontsize=11,
                        clip_on=True,
                    )

    def _get_color(
        self,
        row: pd.Series,
        color_col: Optional[str],
    ) -> str:
        """Get the color for an interval."""
        if color_col and pd.notna(row.get(color_col)):
            return str(row[color_col])

        name_col = self._find_column(
            row.to_frame().T, ["name", "interval", "interval_name"]
        )
        if name_col:
            name = str(row[name_col])
            if name in self._colors:
                return self._colors[name]

        return "#CCCCCC"

    def _plot_age_ticks(
        self,
        ax: mpl_axes.Axes,
        min_age: float,
        max_age: float,
    ) -> None:
        """Plot age tick marks and labels."""
        tick_start = (
            int(np.ceil(min_age / self.age_tick_interval)) * self.age_tick_interval
        )
        tick_end = (
            int(np.floor(max_age / self.age_tick_interval)) * self.age_tick_interval
        )

        ticks = np.arange(
            tick_start, tick_end + self.age_tick_interval, self.age_tick_interval
        )

        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{int(t)}" for t in ticks], fontsize=8)

    def _set_axes_style(
        self,
        ax: mpl_axes.Axes,
        min_age: float,
        max_age: float,
    ) -> None:
        """Set axes style and limits."""
        ax.set_xlim(max_age + self.padding_top, max(0, min_age - self.padding_bottom))
        ax.set_ylim(0, 1)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(False)

        ax.tick_params(
            left=False, labelleft=False, bottom=True, top=False, labelbottom=True
        )
        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")

        ax.set_xlabel("Age (Ma)", fontsize=9)

    def get_intervals_for_age(
        self,
        age: float,
        geo_data: Optional[pd.DataFrame] = None,
    ) -> List[Dict[str, Any]]:
        """
        Get geological intervals that contain a given age.

        Parameters
        ----------
        age : float
            The age in Ma.
        geo_data : pd.DataFrame, optional
            Timescale data.

        Returns
        -------
        List[Dict[str, Any]]
            List of intervals containing the age.
        """
        if geo_data is None:
            if self._geo_data is None:
                self._geo_data = load_ics_data()
            geo_data = self._geo_data

        start_col = self._find_column(
            geo_data, ["start", "age_top", "top", "start_age"]
        )
        end_col = self._find_column(
            geo_data, ["end", "age_bottom", "bottom", "end_age"]
        )

        if start_col is None or end_col is None:
            return []

        # 区间包含 age 的判据是 young <= age <= old；此前写成
        # (start >= age) & (end <= age)，对任何 start<end 的数据恒为空。
        young = geo_data[[start_col, end_col]].min(axis=1)
        old = geo_data[[start_col, end_col]].max(axis=1)
        mask = (young <= age) & (old >= age)
        matching = geo_data[mask]

        # pandas 的 to_dict("records") 在标注上是 dict[str, Any] 的列表；
        # 包一层 list() 把 Any 收拢成显式列表（运行期 to_dict 本来就返回 list）。
        return list(matching.to_dict("records"))


def create_geo_axis(
    fig: Any, position: Tuple[float, float, float, float], **kwargs: Any
) -> mpl_axes.Axes:
    """
    Create a geological timescale axis.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The parent figure.
    position : Tuple[float, float, float, float]
        The [left, bottom, width, height] position.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    matplotlib.axes.Axes
        The created axes.
    """
    ax: mpl_axes.Axes = fig.add_axes(position)
    return ax
