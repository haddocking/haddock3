"""Plotting functionalities."""

import json
import shutil

import numpy as np
import pandas as pd
import plotly.colors as px_colors
import plotly.express as px
import plotly.graph_objects as go
from plotly.io._utils import plotly_cdn_url
from plotly.offline.offline import get_plotlyjs
from plotly.subplots import make_subplots

from pathlib import Path

from haddock import log
from haddock.core.typing import (
    Any,
    DataFrameGroupBy,
    Figure,
    FilePath,
    ImgFormat,
    NDFloat,
    Optional,
    Union,
    )
from haddock.libs.assets import haddock_ui_path


SCATTER_PAIRS = [
    ("irmsd", "score"),
    ("irmsd", "desolv"),
    ("irmsd", "vdw"),
    ("irmsd", "elec"),
    ("irmsd", "air"),
    ("dockq", "score"),
    ("dockq", "desolv"),
    ("dockq", "vdw"),
    ("dockq", "elec"),
    ("dockq", "air"),
    ("lrmsd", "score"),
    ("lrmsd", "desolv"),
    ("lrmsd", "vdw"),
    ("lrmsd", "elec"),
    ("lrmsd", "air"),
    ("ilrmsd", "score"),
    ("ilrmsd", "desolv"),
    ("ilrmsd", "vdw"),
    ("ilrmsd", "elec"),
    ("ilrmsd", "air"),
    ("fnat", "score"),
    ("fnat", "desolv"),
    ("fnat", "vdw"),
    ("fnat", "elec"),
    ("fnat", "air"),
    ]


# if SCATTER_PAIRS changes, SCATTER_MATRIX_SIZE should change too!
SCATTER_MATRIX_SIZE = (5, 5)  # (number of rows, number of columns)

TITLE_NAMES = {
    "score": "HADDOCK score",
    "irmsd": "i-RMSD",
    "lrmsd": "l-RMSD",
    "ilrmsd": "il-RMSD",
    "dockq": "DOCKQ",
    "desolv": "Edesolv",
    "vdw": "Evdw",
    "elec": "Eelec",
    "air": "Eair",
    "fnat": "FCC",
    "bsa": "BSA",
    }

AXIS_NAMES = {
    "score": "HADDOCK score [a.u.]",
    "vdw": "Van der Waals Energy",
    "elec": "Electrostatic Energy",
    "air": "Restraints Energy",
    "desolv": "Desolvation Energy",
    "irmsd": "interface RMSD [A]",
    "lrmsd": "ligand RMSD [A]",
    "ilrmsd": "interface-ligand RMSD [A]",
    "fnat": "Fraction of Common Contacts",
    "dockq": "DOCKQ",
    "bsa": "Buried Surface Area [A^2]",
    }

ClRank = dict[int, int]
"""
A dict representing clusters' rank.

key  (int): cluster's id

value(int): cluster's rank
"""

HEATMAP_DEFAULT_PATH = Path('contacts.html')
SUPPORTED_OUTPUT_FORMATS = ('png', 'jpeg', 'webp', 'svg', 'pdf', 'eps', )

def create_html(
        json_content: str,
        plot_id: int = 1,
        plotly_js_import: Optional[str] = None,
        figure_height: int = 800,
        figure_width: int = 1000,
        ) -> str:
    """Create html content given a plotly json.

    Parameters
    ----------
    json_content : str
        plotly json content
    
    plot_id : int
        plot id to be used in the html content
    
    figure_height : int
        figure height (in pixels)
    
    figure_width : int
        figure width (in pixels)
    
    Returns
    -------
    html_content : str
        html content
    """
    # Check if plotly javascript must be flushed in this file
    if not plotly_js_import:
        plotly_js_import = f'<script src="{plotly_cdn_url()}"></script>'

    # Write HTML content
    html_content = f"""
    <div>
    <script type="text/javascript">window.PlotlyConfig = {{ MathJaxConfig: 'local' }};</script>
    {plotly_js_import}
    <div id="plot{plot_id}" class="plotly-graph-div" style="height:{figure_height}px; width:{figure_width}px;">
    </div>
    <script id="data{plot_id}" type="application/json">
    {json_content}
    </script>
    <script type="text/javascript">
        const dat{plot_id} = JSON.parse(document.getElementById("data{plot_id}").text)
        window.PLOTLYENV = window.PLOTLYENV || {{}};
        if (document.getElementById("plot{plot_id}")) {{
            Plotly.newPlot(
                "plot{plot_id}",
                dat{plot_id}.data,
                dat{plot_id}.layout,
                {{ responsive: true }},
            );
        }}
    </script>
    </div>
    """  # noqa : E501
    return html_content


def read_capri_table(
        capri_filename: FilePath,
        comment: str = "#",
        ) -> pd.DataFrame:
    """Read capri table with pandas.

    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    comment : str
        the string used to denote a commented line in capri tables

    Returns
    -------
    capri_df : pandas DataFrame
        dataframe of capri values
    """
    capri_df = pd.read_csv(capri_filename, sep="\t", comment=comment)
    return capri_df


def in_capri(column: str, df_columns: pd.Index) -> bool:
    """
    Check if the selected column is in the set of available columns.

    Parameters
    ----------
    column : str
        column name
    df_columns : pandas.DataFrame.columns
        columns of a pandas.DataFrame

    Returns
    -------
    resp : bool
        if True, the column is present
    """
    resp = True
    if column not in df_columns:
        log.warning(f"quantity {column} not present in capri table")
        resp = False
    return resp


def update_layout_plotly(
        fig: Figure,
        x_label: str,
        y_label: str,
        title: Optional[str] = None,
        ) -> Figure:
    """
    Update layout of plotly plot.

    Parameters
    ----------
    fig : plotly Figure
        figure
    x_label : str
        x axis name
    y_label : str
        y axis name
    title : str or None
        plot title
    """
    px_dict = {
        "title": title,
        "xaxis": dict(
            title=dict(text=x_label, font=dict(size=40)),
            tickfont_size=14,
            ),
        "yaxis": dict(
            title=dict(text=y_label, font=dict(size=40)),
            tickfont_size=14,
            ),
        "legend": dict(x=1.01, y=1.0, font_family="Helvetica", font_size=16),
        "hoverlabel": dict(font_size=16, font_family="Helvetica"),
        }
    fig.update_layout(px_dict)
    return fig


def box_plot_plotly(
        gb_full: pd.DataFrame,
        y_ax: str,
        cl_rank: dict[int, int],
        format: Optional[ImgFormat],
        scale: Optional[float],
        offline: bool = False,
        ) -> Figure:
    """
    Create a scatter plot in plotly.

    Parameters
    ----------
    gb_full : pandas DataFrame
        data to box plot
    y_ax : str
        variable to plot
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary
    format : str
        Produce images in the selected format.
    scale : int
        scale of image

    Returns
    -------
    fig_list : list
        a list of figures
    """
    colors = px_colors.qualitative.Dark24
    color_map = {}
    for cl_id in sorted(cl_rank.keys()):
        color_idx = (cl_rank[cl_id] - 1) % len(colors)  # color index

        # Note: the rank format (float/int) in "cl_rank" is different from
        # gb_full["cluster_ranking"]
        rns = gb_full[gb_full["cluster_id"] == cl_id]["cluster_ranking"]
        rn = rns.unique()[0]
        color_map[f"{rn}"] = colors[color_idx]

        # Choose a different color for "Other" like in scatter plots
        color_map["Other"] = "#DDDBDA"

    # to use color_discrete_map, cluster_ranking column should be str not int
    gb_full_string = gb_full.astype({"cluster_ranking": "string"})

    # Rename for a better name in legend
    gb_full_string.rename(
        columns={"cluster_ranking": "Cluster Rank"},
        inplace=True,
        )

    # "Cluster Rank" is equivalent to "capri_rank"!
    fig = px.box(
        gb_full_string,
        x="capri_rank",
        y=f"{y_ax}",
        color="Cluster Rank",
        color_discrete_map=color_map,
        boxmode="overlay",
        points="outliers",
        width=1000,
        height=800,
        hover_data=["caprieval_rank"],
        )
    # layout
    update_layout_plotly(fig, "Cluster Rank", AXIS_NAMES[y_ax])
    # save figure
    px_fpath = Path(f"{y_ax}_clt.html")
    json_content = fig.to_json()
    html_content = create_html(
        json_content,
        plotly_js_import=offline_js_manager(px_fpath, offline),
        )
    # write html_content to px_fname
    px_fpath.write_text(html_content)
    # create format boxplot if necessary
    if format:
        fig.write_image(f"{y_ax}_clt.{format}", scale=scale)
    return fig


def box_plot_data(capri_df: pd.DataFrame, cl_rank: ClRank) -> pd.DataFrame:
    """
    Retrieve box plot data.

    Parameters
    ----------
    capri_df : pandas DataFrame
        capri table dataframe
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary

    Returns
    -------
    gb_full : pandas DataFrame
        DataFrame of all the clusters to be plotted
    """
    gb_cluster = capri_df.groupby("cluster_id")
    gb_other = pd.DataFrame([])
    gb_good = pd.DataFrame([])
    for cl_id, cl_df in gb_cluster:
        if cl_id not in cl_rank.keys():
            gb_other = pd.concat([gb_other, cl_df])
        else:
            cl_df["capri_rank"] = cl_rank[cl_id]  # type: ignore
            gb_good = pd.concat([gb_good, cl_df])

    gb_other["cluster_id"] = "Other"
    gb_other["capri_rank"] = len(cl_rank.keys()) + 1
    gb_other["cluster_ranking"] = "Other"
    gb_full = pd.concat([gb_good, gb_other])

    # Sort based on "capri_rank"
    gb_full.sort_values(by=["capri_rank"], inplace=True)
    return gb_full


def box_plot_handler(
        capri_filename: FilePath,
        cl_rank: ClRank,
        format: Optional[ImgFormat],
        scale: Optional[float],
        offline: bool = False,
        ) -> list[Figure]:
    """Create box plots.

    The idea is that for each of the top X-ranked clusters we create a box plot
    showing how the basic statistics are distributed within each model.

    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary
    format : str
        Produce images in the selected format.
    scale : int
        scale for images.
    """
    # generating the correct dataframe
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_full = box_plot_data(capri_df, cl_rank)
    
    # iterate over the variables
    fig_list: list[Figure] = []
    for y_ax in AXIS_NAMES.keys():
        if not in_capri(y_ax, capri_df.columns):
            continue
        fig = box_plot_plotly(
            gb_full, y_ax, cl_rank, format, scale, offline=offline,
            )
        fig_list.append(fig)
    return fig_list


def scatter_plot_plotly(
        gb_cluster: DataFrameGroupBy,
        gb_other: pd.DataFrame,
        cl_rank: ClRank,
        x_ax: str,
        y_ax: str,
        colors: list[str],
        format: Optional[ImgFormat],
        scale: Optional[float],
        offline: bool = False,
        ) -> Figure:
    """Create a scatter plot in plotly.

    Parameters
    ----------
    gb_cluster : pandas DataFrameGroupBy
        capri DataFrame grouped by cluster_id
    gb_other : pandas DataFrame
        DataFrame of clusters not in the top cluster ranking
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary
    x_ax : str
        name of the x column
    y_ax : str
        name of the y column
    colors : list
        list of colors to be used
    format : str
        Produce images in the selected format.
    scale : int
        scale for images.

    Returns
    -------
    fig :
        an instance of plotly.graph_objects.Figure
    """

    def _build_hover_text(df):
        """Build a nice text for hover text."""
        text_list = []
        for _, row in df.iterrows():
            model_text = f"Model: {row['model'].split('/')[-1]}"
            score_text = f"Score: {row['score']}"
            caprieval_rank_text = f"Caprieval rank: {row['caprieval_rank']}"
            text_list.append(
                f"{model_text}<br>{score_text}"
                f"<br>{caprieval_rank_text}"
                )
        return text_list

    fig = go.Figure(layout={"width": 1000, "height": 800})
    traces: list[go.Scatter] = []
    n_colors = len(colors)
    cl_rank_swap = {v: k for k, v in cl_rank.items()}

    for cl_rn in sorted(cl_rank_swap.keys()):
        cl_id = cl_rank_swap[cl_rn]
        cl_df = gb_cluster.get_group(cl_id)
        if cl_id in cl_rank.keys():
            if cl_id == "-":
                cl_name = "Unclustered"
            else:
                cl_name = f"Cluster {cl_rank[cl_id]}"  # use rank
            color_idx = (cl_rank[cl_id] - 1) % n_colors  # color index
            
            traces.append(
                go.Scatter(
                    x=cl_df[x_ax],
                    y=cl_df[y_ax],
                    name=cl_name,
                    mode="markers",
                    text=_build_hover_text(cl_df),
                    legendgroup=cl_name,
                    marker_color=colors[color_idx],
                    hoverlabel=dict(
                        bgcolor=colors[color_idx],
                        font_size=16,
                        font_family="Helvetica",
                        ),
                    )
                )
            clt_text = f"{cl_name}<br>"
            
            # mean and std deviations for the top 4 members
            x_mean = np.mean(cl_df[x_ax].iloc[:4])
            y_mean = np.mean(cl_df[y_ax].iloc[:4])
            x_std = np.std(cl_df[x_ax].iloc[:4])
            y_std = np.std(cl_df[y_ax].iloc[:4])

            if "score" not in [x_ax, y_ax]:
                clt_text += f"Score: {np.mean(cl_df['score'].iloc[:4]):.3f}<br>"
            clt_text += f"{x_ax}: {x_mean:.3f}<br>{y_ax}: {y_mean:.3f}"

            clt_text_list = [clt_text]
            traces.append(
                go.Scatter(
                    x=[x_mean],
                    y=[y_mean],
                    # error bars
                    error_x=dict(
                        type="data", array=[x_std], visible=True
                        ),
                    error_y=dict(
                        type="data", array=[y_std], visible=True
                        ),
                    # color and text
                    marker_color=colors[color_idx],
                    text=clt_text_list,
                    legendgroup=cl_name,
                    showlegend=False,
                    mode="markers",
                    marker=dict(size=10, symbol="square-dot"),
                    hovertemplate=f"<b>{clt_text}</b><extra></extra>",
                    hoverlabel=dict(
                        bgcolor=colors[color_idx],
                        font_size=16,
                        font_family="Helvetica",
                        ),
                    )
                )
    # append trace other
    if not gb_other.empty:
        traces.append(
            go.Scatter(
                x=gb_other[x_ax],
                y=gb_other[y_ax],
                name="Other",
                mode="markers",
                text=_build_hover_text(gb_other),
                legendgroup="Other",
                marker=dict(
                    color="white",
                    line=dict(width=2, color="DarkSlateGrey"),
                    ),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=16,
                    font_family="Helvetica",
                    ),
                )
            )
    for trace in traces:
        fig.add_trace(trace)
    px_fpath = Path(f"{x_ax}_{y_ax}.html")
    update_layout_plotly(
        fig,
        TITLE_NAMES[x_ax],
        TITLE_NAMES[y_ax],
        title=f"{TITLE_NAMES[x_ax]} vs {TITLE_NAMES[y_ax]}",
        )
    json_content = fig.to_json()
    html_content = create_html(
        json_content,
        plotly_js_import=offline_js_manager(px_fpath, offline),
        )
    # write html_content to px_fname
    Path(px_fpath).write_text(html_content)

    # create format boxplot if necessary
    if format:
        fig.write_image(f"{x_ax}_{y_ax}.{format}", scale=scale)
    return fig


def scatter_plot_data(
        capri_df: pd.DataFrame,
        cl_rank: ClRank,
        ) -> tuple[DataFrameGroupBy, pd.DataFrame]:
    """Retrieve scatter plot data.

    Parameters
    ----------
    capri_df : pandas DataFrame
        capri table dataframe
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary

    Returns
    -------
    gb_cluster : pandas DataFrameGroupBy
        capri DataFrame grouped by cluster_id
    gb_other : pandas DataFrame
        DataFrame of clusters not in the top cluster ranking
    """
    gb_cluster = capri_df.groupby("cluster_id")
    gb_other = pd.DataFrame([])
    for cl_id, cl_df in gb_cluster:
        if cl_id not in cl_rank.keys():
            gb_other = pd.concat([gb_other, cl_df])
    return gb_cluster, gb_other


def scatter_plot_handler(
        capri_filename: FilePath,
        cl_rank: ClRank,
        format: Optional[ImgFormat],
        scale: Optional[float],
        offline: bool = False,
        ) -> list[Figure]:
    """Create scatter plots.

    The idea is that for each pair of variables of interest (SCATTER_PAIRS,
     declared as global) we create a scatter plot.
    If available, each scatter plot containts cluster information.

    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary
    format : str
        Produce images in the selected format.
    scale : int
        scale for images.

    Returns
    -------
    fig_list : list
        a list of figures
    """
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_cluster, gb_other = scatter_plot_data(capri_df, cl_rank)

    # defining colors
    colors = px_colors.qualitative.Dark24
    fig_list = []
    for x_ax, y_ax in SCATTER_PAIRS:
        if not in_capri(x_ax, capri_df.columns):
            continue
        if not in_capri(y_ax, capri_df.columns):
            continue
        fig = scatter_plot_plotly(
            gb_cluster,
            gb_other,
            cl_rank,
            x_ax,
            y_ax,
            colors,
            format,
            scale,
            offline=offline,
            )
        fig_list.append(fig)
    return fig_list


def _report_grid_size(plot_list: list[Figure]) -> tuple[int, int, int, int]:
    """
    Calculate the size of the grid in the report.

    By size, it means the number of rows/columns, and the width/height
    of an individual plot. In the report, some of the axes are shared. The
    settings for sharing axes depends on the type (scatters or boxes). The
    number of columns is set to the number of columns in SCATTER_MATRIX_SIZE. If
    the number of clusters is more than 5, it increases the width which causes
    horizontal scrolling in the report.

    Parameters
    ----------
    plot_list : list
        list of plots generated by analyse command

    Returns
    -------
    number_of_rows : int
        number of rows in the grid
    number_of_cols : int
        number of columns in the grid
    width : int
        the width of an individual plot
    height: int
        the height of an individual plot
    """
    # Calculate grid size for subplots
    number_of_plots = len(plot_list)
    number_of_clusters = len({trace.legendgroup for trace in plot_list[0].data})
    # same number of cols for both boxes and scatters
    number_of_cols = SCATTER_MATRIX_SIZE[1]
    number_of_rows = int(np.ceil(number_of_plots / number_of_cols))
    # enable horizontal scroll
    width = 600 if number_of_clusters > 5 else 350
    height = 600 if number_of_clusters > 5 else 350
    return number_of_rows, number_of_cols, width, height


def report_plots_handler(plots, shared_xaxes=False, shared_yaxes=False):
    """
    Create a figure that holds subplots.

    The idea is that for each type (scatters or boxes), the individual plots are
    considered subplots. In the report, some of the axes are shared. The
    settings for sharing axes depends on the type (scatters or boxes).

    Parameters
    ----------
    plots : list
        list of plots generated by `analyse` command
    shared_xaxes: boolean or str (default False)
        a parameter of plotly.subplots.make_subplots
    shared_yaxes: boolean or str (default False)
        a parameter of plotly.subplots.make_subplots

    Returns
    -------
    fig :
        an instance of plotly.graph_objects.Figure
    """
    number_of_rows, number_of_cols, width, height = _report_grid_size(plots)
    fig = make_subplots(
        rows=number_of_rows,
        cols=number_of_cols,
        shared_xaxes=shared_xaxes,
        shared_yaxes=shared_yaxes,
        vertical_spacing=(0.4 / number_of_rows),
        horizontal_spacing=(0.3 / number_of_cols),
        )
    for i, sub_fig in enumerate(plots):
        col_index = int((i % number_of_cols) + 1)
        row_index = int(np.floor(i / number_of_cols) + 1)
        # hide legend of plots except one
        if i != 0:
            sub_fig.for_each_trace(lambda trace: trace.update(showlegend=False))
        fig.add_traces(sub_fig.data, rows=row_index, cols=col_index)
        fig.update_yaxes(
            title_text=sub_fig.layout.yaxis.title.text,
            row=row_index,
            col=col_index,
            title_standoff=5,
            automargin=True,
            )
        # x title only on the last row
        if shared_xaxes == "all":
            row_index = number_of_rows
        fig.update_xaxes(
            title_text=sub_fig.layout.xaxis.title.text,
            row=row_index,
            col=col_index,
            title_standoff=5,
            automargin=True,
            )
        legend_title_text = sub_fig.layout.legend.title.text
    fig.update_layout(
        legend_title_text=legend_title_text,
        height=height * number_of_rows,
        width=width * number_of_cols,
        )
    return fig


def find_best_struct(
        df: pd.DataFrame,
        max_best_structs: int = 4,
        ) -> pd.DataFrame:
    """Find best structures for each cluster.

    Parameters
    ----------
    df: pd.DataFrame
        The loaded capri_ss.tsv dataframe
    max_best_structs: int
        The maximum number of best structures to return.

    Returns
    -------
    best_df: pd.DataFrame
        DataFrame of best structures with
        `cluster_id` and `best<model-cluster_ranking>` columns
        and empty strings for missing values.
    """
    df = df[["cluster_id", "model-cluster_ranking", "model"]]
    df = df[df["model-cluster_ranking"] <= max_best_structs]

    best_df = df.pivot(
        index="cluster_id",
        columns="model-cluster_ranking",
        values="model",
        )

    best_df = best_df.fillna('').reset_index()
    best_df.columns = [
        f"best{col}" if col != "cluster_id" else col
        for col in best_df.columns
        ]
    # Remove empty columns
    best_df = best_df.loc[:, (best_df != '').any(axis=0)]
    return best_df


def clean_capri_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a tidy capri table for the report.

    It also combines mean and std values in one column.
    Also it drops the columns that are not needed in the report.

    Makes inplace changes to the dataframe.

    Parameters
    ----------
    df : pandas DataFrame
        dataframe of capri values

    Returns
    -------
    pandas DataFrame
        DataFrame of capri table with new column names
    """
    for col_name in AXIS_NAMES.keys():
        if not in_capri(col_name, df.columns):
            continue
        mean_value = df[col_name]
        std_value = df[f"{col_name}_std"]
        df[col_name] = [
            {'mean': mean_value, 'std': std_value}
            for mean_value, std_value in zip(mean_value, std_value)
            ]

    # Drop columns ending with '_std'
    df = df.drop(df.filter(regex='_std$').columns, axis=1)
    return df


def create_other_cluster(
        clusters_df: pd.DataFrame,
        structs_df: pd.DataFrame,
        max_clusters: int,
        ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Combine all clusters with rank >= max_clusters into an "Other" cluster.

    Parameters
    ----------
    clusters_df : pandas DataFrame
        DataFrame of clusters
    structs_df : pandas DataFrame
        DataFrame of structures
    max_clusters : int
        From which cluster rank to consider as "Other"

    Returns
    -------
        tuple with clusters_df and structs_df
    """
    if len(clusters_df) <= max_clusters:
        return clusters_df, structs_df
    # other clusters
    other_structs_df = structs_df[structs_df['cluster_ranking'] >= max_clusters].copy()
    # drop other clusters from structs_df
    structs_df = structs_df[structs_df['cluster_ranking'] < max_clusters].copy()
    other_structs_df.loc[:,'cluster_id'] = 'Other'
    other_structs_df.loc[:,'cluster_ranking'] = max_clusters
    inner_rank = other_structs_df['caprieval_rank'].rank(method='first').astype(int)  # noqa : E501
    other_structs_df['model-cluster_ranking'] = inner_rank
    structs_df = pd.concat([structs_df, other_structs_df])

    clusters_df = clusters_df[clusters_df['cluster_rank'] < max_clusters]
    other_cluster = {
        'cluster_id': 'Other',
        'cluster_rank': max_clusters,
        'n': len(other_structs_df),
        'caprieval_rank': max_clusters,
        }
    for col in AXIS_NAMES.keys():
        if any([
                col not in other_structs_df.columns,
                col not in clusters_df.columns,
                ]):
            continue
        other_cluster[col] = other_structs_df[col].mean()
        other_cluster[col + '_std'] = other_structs_df[col].std()
    other_cluster_df = pd.DataFrame([other_cluster])
    clusters_df = pd.concat([clusters_df, other_cluster_df], ignore_index=True).round(2)
    return clusters_df, structs_df


def clt_table_handler(
        clt_file: FilePath,
        ss_file: FilePath,
        is_cleaned: bool = False,
        topX_clusters: int = 10,
        clustered_topX: int = 4,
        unclustered_topX: int = 10,
        top_ranked_mapping: Optional[dict[Path, Path]] = None,
        ) -> pd.DataFrame:
    """
    Create a dataframe including data for tables.

    The idea is to create tidy tables that report statistics available in
    capri_clt.tsv and capri_ss.tsv files.

    Parameters
    ----------
    clt_file : str or Path
        path to capri_clt.tsv file
    ss_file: str or Path
        path to capri_ss.tsv file
    is_cleaned: bool
        is the run going to be cleaned?

    Returns
    -------
    df_merged : pandas DataFrame
        a data frame including data for tables
    """
    # table of statistics
    clusters_df = read_capri_table(clt_file)
    structs_df = read_capri_table(ss_file)

    # Round all numbers to 2 decimal places
    clusters_df = clusters_df.round(2)
    structs_df = structs_df.round(2)

    # if the run will be cleaned, the structures are going to be gzipped
    if not top_ranked_mapping:
        if is_cleaned:
            # substitute the values in the df by adding .gz at the end
            structs_df['model'] = structs_df['model'].replace(
                to_replace=r"(\.pdb)$", value=r".pdb.gz", regex=True,
            )
    else:
        # ss_file is in NN_caprieval/ while report is in
        # analysis/NN_caprieval_analysis/
        # need to correct model paths by prepending ../
        def correct_relative_paths(
                path: str,
                top_ranked_mapping: Optional[dict[Path, Path]],
                ) -> str:
            """Prepend model paths in capri_ss files get their relative paths.

            Parameters
            ----------
            path : str
                Original path to a model file.
            top_ranked_mapping : Optional[dict[Path, Path]]
                Optional filepath mapping of top ranked models.

            Returns
            -------
            new_path : str
                New path to the file
            """
            try:
                # If top ranked is provided, use that information
                new_path = top_ranked_mapping[path]
            except (KeyError, TypeError, ):
                # Otherwise just prepend by `../`
                new_path = f"../{path}"
            return new_path
        
        structs_df['model'] = structs_df['model'].apply(
            lambda x: correct_relative_paths(x, top_ranked_mapping)
            )

    is_unclustered = clusters_df["cluster_rank"].unique().tolist() == ["-"]
    # If unclustered, we only want to show the top 10 structures in a table.
    if is_unclustered:
        structs_df = structs_df[:unclustered_topX]
        cols2keep = ['caprieval_rank', 'model'] + list(AXIS_NAMES.keys())
        structs_df = structs_df[cols2keep]
        # model has ../../01_rigidbody/rigidbody_62.pdb.gz
        # add id column with 62 as value
        structs_df['id'] = structs_df['model'].str.extract(r'(\d+).pdb')
        return structs_df

    clusters_df, structs_df = create_other_cluster(
        clusters_df,
        structs_df,
        max_clusters=topX_clusters + 1,
        )

    clusters_df = clean_capri_table(clusters_df)
    structs_df = find_best_struct(structs_df, max_best_structs=clustered_topX)
    df_merged = pd.merge(clusters_df, structs_df, on="cluster_id")
    return df_merged


def _css_styles_for_report(offline: bool) -> str:
    """
    Generate custom CSS styles for an analysis report.

    Parameters
    ----------
    offline : bool
        If True, the HTML will be generated for offline use.

    Returns
    -------
    The CSS styles as a string.
    """
    custom_css = """
    .title {
        font-family: Arial, sans-serif;
        font-size: 32px;
        font-weight: bold;
    }
    body {
      margin-left: 1em;
    }
    table {
        border-collapse: collapse;
    }
    th {
        background-color: #f2f2f2;
        padding: 8px;
        border: 1px solid #ddd;
        text-align: left;
    }
    th[scope="row"] {
        position: sticky;
        min-width: 16rem;
        left: 0;
        z-index: 1
    }
    td {
        border: 1px solid #ddd;
        padding: 8px;
        text-align: left;
    }
    tr:nth-child(even) {
        background-color: #f2f2f2
    }
    .js-plotly-plot .plotly .modebar svg {
	    display: inline;
    }
    """
    css_link = "https://cdn.jsdelivr.net/npm/@i-vresse/haddock3-ui@~0.3.0/dist/index.css"
    if offline:
        # copy the css file to the report directory
        src = haddock_ui_path / 'index.css'
        shutil.copyfile(str(src), "../data/ui/index.css")
        css_link = "../../data/ui/index.css"
    table_css = f' <link href="{css_link}" rel="stylesheet" />'
    return f"{table_css}<style>{custom_css}</style>"


def _generate_html_report(
        step: str,
        figures: list[Union[Figure, pd.DataFrame]],
        report_path: FilePath,
        offline: bool = False,
        ) -> str:
    """
    Generate an HTML report for a specific step of analysis, including figures.

    Parameters
    ----------
    step : str
        The step number.
    figures : list
        A list of figures to include in the HTML body.
        Each figure can be either a string representing a table or a
        plotly.graph_objects.Figure object.
    offline : bool
        If True, the HTML will be generated for offline use.

    Returns
    -------
    html_report : str
        The generated HTML report as a string.
    """
    html_report = "<!DOCTYPE html><html lang='en'>"
    html_report += _generate_html_head(step, offline)
    html_report += _generate_html_body(figures, report_path, offline=offline)
    html_report += "</html>"
    return html_report


def _generate_html_head(step, offline):
    """
    Generate the HTML head section for an analysis report.

    Parameters
    ----------
    step : str
        The step number.
    offline : bool
        If True, the HTML will be generated for offline use.

    Returns
    -------
    head : str
        The HTML head section as a string.
    """
    head = "<head>"
    head += f"<title>Analysis report of step {step}</title>"
    head += f"<p class='title'>Analysis report of step {step}</p>"
    head += _css_styles_for_report(offline)
    head += "</head>"
    return head


def _generate_unclustered_table_html(
        table_id: str,
        df: pd.DataFrame,
        bundle_url: str,
        ) -> str:
    data = df.to_json(orient='records')
    headers = [
        {'key': "caprieval_rank", 'label': "Structure Rank", 'sorted': "asc"},
        {'key': "model", 'label': "Structure",
         'sortable': False, 'type': "structure"
         },
        ] + [
            {'key': k, 'label': v, 'type': 'stats'}
            for k, v in AXIS_NAMES.items()
            if k in df.columns
            ] + [
        {'key': "id", 'label': "Structure ID"},
        ]
    return f"""
            <div id="{table_id}"></div>
            <script id="data{table_id}" type="application/json">
            {{
                "structures": {data},
                "headers": {json.dumps(headers)}
            }}
            </script>
            <script type="module">
            import {{ renderStructureTable }} from "{bundle_url}";

            const props = JSON.parse(document.getElementById("data{table_id}").text)

            renderStructureTable(document.getElementById('{table_id}'), props.headers, props.structures)
            </script>"""  # noqa : E501


def _generate_clustered_table_html(
        table_id: str,
        df: pd.DataFrame,
        bundle_url: str,
        ) -> str:
    data = df.to_json(orient='records')
    nr_best_columns = df.filter(like="best").shape[1]
    headers = [
        {'key': "cluster_rank", 'label': "Cluster Rank", 'sorted': "asc"},
        {'key': "cluster_id", 'label': "Cluster ID"},
        {'key': "n", 'label': "Cluster size"},
        ] + [
        {'key': k, 'label': v, 'type': 'stats'}
        for k, v in AXIS_NAMES.items()
        if k in df.columns
        ] + [
        {'key': f"best{i}", 'label': f"Nr {i} best structure",
         'sortable': False, 'type': "structure"
         }
        for i in range(1, nr_best_columns + 1)
        ]

    caption = ''
    if df['cluster_id'].isin(['Other']).any():
        caption = (
            'The "Other" cluster is not a real cluster it contains'
            'all structures that are not in the top 10 clusters.'
            )

    return f"""
            <div id="{table_id}"></div>
            <div>{caption}</div>
            <script id="data{table_id}" type="application/json">
            {{
                "clusters": {data},
                "headers": {json.dumps(headers)}
            }}
            </script>
            <script type="module">
            import {{ renderClusterTable }} from "{bundle_url}";
            const props = JSON.parse(document.getElementById("data{table_id}").text)

            renderClusterTable(document.getElementById('{table_id}'), props.headers, props.clusters);
            </script>"""  # noqa : E501


def _generate_html_body(
        figures: list[Union[Figure, pd.DataFrame]],
        report_path: FilePath,
        offline: bool = False,
        ) -> str:
    """
    Generate an HTML body section containing figures for an analysis report.

    Parameters
    ----------
    figures : list
        A list of figures to include in the HTML body.
        Each figure can be either a string representing a table or a
        plotly.graph_objects.Figure object.
    offline : bool
        If True, the HTML will be generated for offline use.

    Returns
    -------
    body : str
        The generated HTML body as a string.
    """
    body = "<body>"
    table_index: int = 1
    fig_index: int = 1
    for figure in figures:
        if isinstance(figure, pd.DataFrame):  # tables
            table_index += 1
            table_id = f"table{table_index}"

            is_unclustered = 'cluster_rank' not in figure
            bundle_url = "https://cdn.jsdelivr.net/npm/@i-vresse/haddock3-ui@~0.3.0/dist/report.bundle.js"
            if offline:
                # copy the bundle to the run_dir folder
                src = haddock_ui_path / 'report.bundle.js'
                shutil.copyfile(str(src), "../data/ui/report.bundle.js")
                bundle_url = "../../data/ui/report.bundle.js"
            if is_unclustered:
                inner_html = _generate_unclustered_table_html(table_id, figure, bundle_url)
            else:
                inner_html = _generate_clustered_table_html(table_id, figure, bundle_url)
        else:  # plots
            inner_json = figure.to_json()
            inner_html = create_html(
                inner_json,
                fig_index,
                plotly_js_import=offline_js_manager(report_path, offline),
                figure_height=figure.layout.height,
                figure_width=figure.layout.width,
                )
            fig_index += 1  # type: ignore
        body += "<br>"  # add a break between tables and plots
        body += inner_html
    body += "</body>"
    return body


def report_generator(
        boxes: list[Figure],
        scatters: list[Figure],
        tables: list,
        step: str,
        directory: FilePath = ".",
        offline: bool = False
        ) -> None:
    """
    Create a figure include plots and tables.

    The idea is to create a report.html file that includes all the plots and
    tables generated by the command `analyse`.

    Parameters
    ----------
    boxes : list
        list of box plots generated by box_plot_handler
    scatters: list
        list of scatter plots generated by scatter_plot_handler
    table: list
        a list including tables generated by clt_table_handler
    directory : Path
        path to the output folder
    offline: bool
        If True, the HTML will be generated for offline use.
    """
    figures = [tables]
    # Combine scatters
    figures.append(
        report_plots_handler(
            scatters,
            shared_xaxes="rows",
            shared_yaxes="columns",
            )
        )
    # Combine boxes"
    figures.append(report_plots_handler(boxes))

    if offline:
        Path('../data/ui').mkdir(parents=True, exist_ok=True)
    # Write everything to a html file
    report_path = Path(directory, "report.html")
    html_report = _generate_html_report(step, figures, report_path, offline)
    with open(report_path, "w", encoding="utf-8") as report:
        report.write(html_report)


def heatmap_plotly(
        matrix: NDFloat,
        labels: Optional[dict] = None,
        xlabels: Optional[list] = None,
        ylabels: Optional[list] = None,
        color_scale: str = 'Greys_r',
        title: Optional[str] = None,
        output_fname: Path = HEATMAP_DEFAULT_PATH,
        offline: bool = False,
        hovertemplate: Optional[str] = None,
        customdata: Optional[list[list[Any]]] = None,
        delineation_traces: Optional[list[dict[str, float]]] = None,
        ) -> Path:
    """Generate a `plotly heatmap` based on matrix content.

    Parameters
    ----------
    matrix : NDFloat
        The 2D matrix containing data to be shown.
    labels : dict
        Labels of the horizontal (x), vertical (y) and colorscale (color) axis.
    xlabels : list
        List of columns names.
    ylabels : list
        List of row names.
    color_scale : str
        Color scale to use.
    title : str
        Title of the figure.
    output_fname : Path
        Path to the output filename to generate.
    hovertemplate: Optional[str]
        Custrom string used to format data for hover annotation in plotly.
    customdata: Optional[list[list[list[int]]]]
        A matrix of cluster ids, used for extra hover annotation in plotly.
    delineation_traces: Optional[list[dict[str, float]]]
        A list of dict enabling to draw lines separating cluster ids.

    Return
    ------
    output_fname : Path
        Path to the generated filename
    """
    # Generate heatmap trace
    fig = px.imshow(
        matrix,
        labels=labels,
        x=xlabels,
        y=ylabels,
        color_continuous_scale=color_scale,
        title=title,
        )
    # Place X axis on top
    fig.update_xaxes(side="top")
    fig.update_traces(
        hovertemplate=hovertemplate,
        customdata=customdata,
        )
    # Add delineation traces
    if delineation_traces:
        # Loop over lines
        for trace in delineation_traces:
            # Draw them
            fig.add_shape(
                type="line",
                line={"dash": "5px"},
                x0=trace["x0"],
                x1=trace["x1"],
                y0=trace["y0"],
                y1=trace["y1"],
            )

    # Compute pixels
    nb_entries = matrix.shape[0]
    scaled_log = int(np.log(nb_entries)) * 200
    lower_bound = max(scaled_log, 1000)
    uppder_bound = min(lower_bound, 2000)
    # Set hight and width
    height = uppder_bound
    # Increment width for legend space
    width = height + 70

    # Save figure as html file
    export_plotly_figure(
        fig,
        output_fname,
        offline=offline,
        figure_height=height,
        figure_width=width,
        )

    return output_fname


def export_plotly_figure(
        fig: Figure,
        output_fname: Union[str, Path],
        figure_height: int = 1000,
        figure_width: int = 1000,
        offline: bool = False,
        ) -> None:
    """Write a plotly figure.

    Parameters
    ----------
    fig : Figure
        The plotly Figure object
    output_fname : Union[str, Path]
        Where to write it
    figure_height : int, optional
        Height of the figure (in pixels), by default 1000
    figure_width : int, optional
        Width of the figure (in pixels), by default 1000
    offline : bool, optional
        If True add the plotly js library to the file, by default False
    """
    # Detect output file extension
    _suffix = Path(output_fname).suffix
    suffix = _suffix[1:]
    # Check corresponding function
    if suffix == "html":
        fig_to_html(
            fig,
            output_fname,
            figure_height=figure_height,
            figure_width=figure_width,
            offline=offline,
            )
    elif suffix in SUPPORTED_OUTPUT_FORMATS:
        fig.write_image(output_fname)
    
  
def make_alascan_plot(
        df: pd.DataFrame,
        clt_id: int,
        scan_res: str = "ALA",
        offline: bool = False,
        ) -> None:
    """
    Make a plotly interactive plot.

    Score components are here **weighted** by their respective
    contribution to the total score.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the results of the alanine scan.
    clt_id : int
        Cluster ID.
    scan_res : str, optional
        Residue name used for the scan, by default "ALA"
    """
    plot_name = f"scan_clt_{clt_id}"
    log.info(f"Generating {scan_res} scanning plot {plot_name}")

    # create figure
    width, height = 2000, 1000
    fig = go.Figure(layout={"width": width, "height": height})
    # add traces
    fig.add_trace(
        go.Bar(
            x=df["full_resname"],
            y=df["delta_score"],
            name="delta_score",
            )
        )
    
    fig.add_trace(
        go.Bar(
            x=df["full_resname"],
            y=df["delta_vdw"],
            name="delta_vdw",
            )
        )
    # delta_elec is given its weight in the emscoring module
    fig.add_trace(
        go.Bar(
            x=df["full_resname"],
            y=0.2 * df["delta_elec"],
            name="delta_elec",
            )
        )

    fig.add_trace(
        go.Bar(
            x=df["full_resname"],
            y=df["delta_desolv"],
            name="delta_desolv",
            )
        )
    # prettifying layout
    fig.update_layout(
        title=f"{scan_res} scanning cluster {clt_id}",
        xaxis=dict(
            title=dict(text="Residue Name", font=dict(size=16)),
            tickfont_size=14,
            tick0=df["full_resname"],
            # in case we want to show less residues
            # dtick=10,
            ),
        yaxis=dict(
            title=dict(text="Average Delta (WT - mutant)", font=dict(size=16)),
            tickfont_size=14,
            ),
        legend=dict(x=1.01, y=1.0, font_family="Helvetica", font_size=16),
        barmode="group",
        bargap=0.05,
        bargroupgap=0.05,
        hovermode="x unified",
        hoverlabel=dict(font_size=16, font_family="Helvetica"),
        )
    for n in range(df.shape[0] - 1):
        fig.add_vline(x=0.5 + n, line_color="gray", opacity=0.2)

    # save html
    html_output_filename = f"{plot_name}.html"
    export_plotly_figure(
        fig,
        html_output_filename,
        figure_height=height,
        figure_width=width,
        offline=offline,
        )


def fig_to_html(
        fig: Figure,
        fpath: Union[str, Path],
        plot_id: int = 1,
        figure_height: int = 800,
        figure_width: int = 1000,
        offline: bool = False,
        ) -> None:
    """Workaround plotly html file generation.

    Parameters
    ----------
    fig : Figure
        A Figure object created by Plotly
    fpath : Union[str, Path]
        Where to write the content
    json_content : str
        plotly json content
    plot_id : int
        plot id to be used in the html content
    figure_height : int
        figure height (in pixels)
    figure_width : int
        figure width (in pixels)
    offline : bool
        If set to False, use the cdn url to obtain the javascript content
        for the rendering.
    """
    # Convert to json
    json_content = fig.to_json()
    # Create custom html file
    html_content = create_html(
        json_content,
        plot_id=plot_id,
        plotly_js_import=offline_js_manager(fpath, offline),
        figure_height=figure_height,
        figure_width=figure_width,
        )
    # Write it
    Path(fpath).write_text(html_content)


def offline_js_manager(fpath: FilePath, offline: bool) -> str:
    """Build string to access plotly javascript content.

    Parameters
    ----------
    fpath : FilePath
        Path to the figure about to be written.
    offline : bool
        if True use the offline approach.

    Returns
    -------
    plotly_js_import : str
        HTML solution for the importation of the plotly javascript content.
    """
    # Case where offline isrequired
    if offline:
        # Obtain directory where the figure should be written
        fig_dir = Path(fpath).parent
        # Set plotly js filepath
        plotly_js_fpath = Path(fig_dir, "plotly_bundle.js")
        # Check that this file do not already exists
        if not plotly_js_fpath.exists():
            # Write the plotly java script content
            plotly_js_fpath.write_text(get_plotlyjs())
        # Build HTML string
        plotly_js_import = f'<script src="{plotly_js_fpath}"></script>'
    else:
        # Use CDN url to obtain the script
        plotly_js_import = f'<script src="{plotly_cdn_url()}"></script>'
    return plotly_js_import
        

def make_traceback_plot(tr_subset, plot_filename, offline=False):
    """
    Create a traceback barplot with the 40 best ranked models.

    Parameters
    ----------
    tr_subset : pandas.DataFrame
        DataFrame containing the top traceback results
    plot_filename : Path
        Path to the output filename to generate
    """
    rank_columns = tr_subset.columns[tr_subset.columns.str.endswith("rank")]
    # for each row, plot a bar with the values of the rank columns
    fig = px.bar(tr_subset, x="Model", y=rank_columns)
    # vertical legend on the right with legend name 'Modules'
    fig.update_layout(legend_orientation="v", legend_title="Modules")
    # legend should have bigger font size
    fig.update_layout(
        legend=dict(
            title_font_size=24,
            font_size=24,
            ),
        )
    # y axis title 'Sum of Ranks'
    fig.update_layout(yaxis_title="Sum of Ranks", xaxis_title="Models")
    # bigger axis labels
    fig.update_layout(
        yaxis=dict(title_font_size=30, tickfont_size=16),
        xaxis=dict(title_font_size=30, tickfont_size=16),
        )
    # bigger title
    fig.update_layout(
        title_text=f"Top ranked {tr_subset.shape[0]} Models",
        title_font_size=30
        )
    fig_to_html(
        fig,
        plot_filename,
        figure_height=1200,
        figure_width=2000,
        offline=offline,
        )
    return plot_filename
