"""Plotting functionalities."""

import json

import numpy as np
import pandas as pd
import plotly.colors as px_colors
import plotly.express as px
import plotly.graph_objects as go
from plotly.io._utils import plotly_cdn_url
from plotly.subplots import make_subplots

from pathlib import Path

from haddock import log
from haddock.core.typing import (
    DataFrameGroupBy,
    Figure,
    FilePath,
    ImgFormat,
    NDFloat,
    Optional,
    Union,
)
from typing import Any


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


def create_html(
        json_content: str,
        plot_id: int = 1,
        figure_height: int = 800,
        figure_width: int = 1000,
        ) -> str:
    """
    create html content given a plotly json

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
    html_content = f"""
    <div>
    <script type="text/javascript">window.PlotlyConfig = {{ MathJaxConfig: 'local' }};</script>
    <script src="{plotly_cdn_url()}"></script>
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
    """
    return html_content


def read_capri_table(capri_filename: FilePath, comment: str = "#") -> pd.DataFrame:
    """
    Read capri table with pandas.

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
    fig: Figure, x_label: str, y_label: str, title: Optional[str] = None
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
            title=x_label,
            tickfont_size=14,
            titlefont_size=40,
        ),
        "yaxis": dict(
            title=y_label,
            tickfont_size=14,
            titlefont_size=40,
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
        # gb_full["cluster-ranking"]
        rns = gb_full[gb_full["cluster-id"] == cl_id]["cluster-ranking"]
        rn = rns.unique()[0]
        color_map[f"{rn}"] = colors[color_idx]

        # Choose a different color for "Other" like in scatter plots
        color_map["Other"] = "#DDDBDA"

    # to use color_discrete_map, cluster-ranking column should be str not int
    gb_full_string = gb_full.astype({"cluster-ranking": "string"})

    # Rename for a better name in legend
    gb_full_string.rename(columns={"cluster-ranking": "Cluster Rank"}, inplace=True)

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
    px_fname = f"{y_ax}_clt.html"
    json_content = fig.to_json()
    html_content = create_html(json_content)
    # write html_content to px_fname
    Path(px_fname).write_text(html_content)
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
    gb_cluster = capri_df.groupby("cluster-id")
    gb_other = pd.DataFrame([])
    gb_good = pd.DataFrame([])
    for cl_id, cl_df in gb_cluster:
        if cl_id not in cl_rank.keys():
            gb_other = pd.concat([gb_other, cl_df])
        else:
            cl_df["capri_rank"] = cl_rank[cl_id]  # type: ignore
            gb_good = pd.concat([gb_good, cl_df])

    gb_other["cluster-id"] = "Other"
    gb_other["capri_rank"] = len(cl_rank.keys()) + 1
    gb_other["cluster-ranking"] = "Other"
    gb_full = pd.concat([gb_good, gb_other])

    # Sort based on "capri_rank"
    gb_full.sort_values(by=["capri_rank"], inplace=True)
    return gb_full


def box_plot_handler(
    capri_filename: FilePath,
    cl_rank: ClRank,
    format: Optional[ImgFormat],
    scale: Optional[float],
) -> list[Figure]:
    """
    Create box plots.

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
        fig = box_plot_plotly(gb_full, y_ax, cl_rank, format, scale)
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
) -> Figure:
    """
    Create a scatter plot in plotly.

    Parameters
    ----------
    gb_cluster : pandas DataFrameGroupBy
        capri DataFrame grouped by cluster-id
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
            text_list.append(f"{model_text}<br>{score_text}<br>{caprieval_rank_text}")
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
                marker=dict(color="white", line=dict(width=2, color="DarkSlateGrey")),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=16,
                    font_family="Helvetica",
                ),
            )
        )
    for trace in traces:
        fig.add_trace(trace)
    px_fname = f"{x_ax}_{y_ax}.html"
    update_layout_plotly(
        fig,
        TITLE_NAMES[x_ax],
        TITLE_NAMES[y_ax],
        title=f"{TITLE_NAMES[x_ax]} vs {TITLE_NAMES[y_ax]}",
    )
    json_content = fig.to_json()
    html_content = create_html(json_content)
    # write html_content to px_fname
    Path(px_fname).write_text(html_content)
    # create format boxplot if necessary
    if format:
        fig.write_image(f"{x_ax}_{y_ax}.{format}", scale=scale)
    return fig


def scatter_plot_data(
    capri_df: pd.DataFrame, cl_rank: ClRank
) -> tuple[DataFrameGroupBy, pd.DataFrame]:
    """
    Retrieve scatter plot data.

    Parameters
    ----------
    capri_df : pandas DataFrame
        capri table dataframe
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary

    Returns
    -------
    gb_cluster : pandas DataFrameGroupBy
        capri DataFrame grouped by cluster-id
    gb_other : pandas DataFrame
        DataFrame of clusters not in the top cluster ranking
    """
    gb_cluster = capri_df.groupby("cluster-id")
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
) -> list[Figure]:
    """
    Create scatter plots.

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
            gb_cluster, gb_other, cl_rank, x_ax, y_ax, colors, format, scale
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


def find_best_struct(df: pd.DataFrame, max_best_structs: int = 4) -> pd.DataFrame:
    """Find best structures for each cluster.

    Args:
        df: DataFrame of capri_ss.tsv
        max_best_structs: The maximum number of best structures to return.

    Returns:
        DataFrame of best structures with
        `cluster_id` and `best<model-cluster-ranking>` columns
        and empty strings for missing values.
    """
    # capri_ss.tsv has a column named "cluster-id" 
    # capri_clt.tsv has a column named "cluster_id"
    # rename here to make the merge easier
    df.rename(columns={"cluster-id": "cluster_id"}, inplace=True)

    df = df[["cluster_id", "model-cluster-ranking", "model"]]
    df = df[df["model-cluster-ranking"] <= max_best_structs]

    best_df = df.pivot(index="cluster_id", columns="model-cluster-ranking", values="model")

    best_df = best_df.fillna('').reset_index()
    best_df.columns = [f"best{col}" if col != "cluster_id" else col for col in best_df.columns]
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
        df[col_name] = [{'mean': mean_value, 'std': std_value} for mean_value, std_value in zip(mean_value, std_value)]

    # Drop columns ending with '_std'
    df = df.drop(df.filter(regex='_std$').columns, axis=1)
    return df


def clt_table_handler(clt_file, ss_file, is_cleaned=False):
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

    # if the run will be cleaned, the structures are going to be gzipped
    if is_cleaned:
        # substitute the values in the df by adding .gz at the end
        structs_df['model'] = structs_df['model'].replace(
            to_replace=r"(\.pdb)$", value=r".pdb.gz", regex=True
        )

    # ss_file is in NN_caprieval/ while report is in analysis/NN_caprieval_analysis/
    # need to correct model paths by prepending ../ 
    structs_df['model'] = structs_df['model'].apply(lambda x: f"../{x}")

    is_unclustered = clusters_df["cluster_rank"].unique().tolist() == ["-"]
    # If unclustered, we only want to show the top 10 structures in a table.
    if is_unclustered:
        max_unstructured_structures = 10
        structs_df = structs_df[:max_unstructured_structures]
        cols2keep = ['caprieval_rank','model'] + list(AXIS_NAMES.keys())
        structs_df = structs_df[cols2keep]
        # model has ../../01_rigidbody/rigidbody_62.pdb.gz
        # add id column with 62 as value
        structs_df['id'] = structs_df['model'].str.extract(r'(\d+).pdb')
        return structs_df

    clusters_df = clean_capri_table(clusters_df)
    structs_df = find_best_struct(structs_df, max_best_structs=4)
    df_merged = pd.merge(clusters_df, structs_df, on="cluster_id")
    return df_merged


def _css_styles_for_report():
    """
    Generate custom CSS styles for an analysis report.

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

    """
    css_link = "https://esm.sh/@i-vresse/haddock3-analysis-components@~0.4.0-next.0/dist/style.css"  # noqa:E501
    table_css = f' <link href={css_link} rel="stylesheet" />'
    return f"{table_css}<style>{custom_css}</style>"


def _generate_html_report(step, figures):
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

    Returns
    -------
    html_report : str
        The generated HTML report as a string.
    """
    html_report = "<!DOCTYPE html><html lang='en'>"
    html_report += _generate_html_head(step)
    html_report += _generate_html_body(figures)
    html_report += "</html>"
    return html_report


def _generate_html_head(step):
    """
    Generate the HTML head section for an analysis report.

    Parameters
    ----------
    step : str
        The step number.

    Returns
    -------
    head : str
        The HTML head section as a string.
    """
    head = "<head>"
    head += f"<title>Analysis report of step {step}</title>"
    head += f"<p class='title'>Analysis report of step {step}</p>"
    head += _css_styles_for_report()
    head += """
            <script type="importmap">
            {
                "imports": {
                    "react": "https://esm.sh/react@^18.2.0",
                    "react-dom": "https://esm.sh/react-dom@^18.2.0",
                    "@i-vresse/haddock3-analysis-components": "https://esm.sh/@i-vresse/haddock3-analysis-components@~0.4.0-next.0?bundle"
                }
            }
            </script>"""
    head += "</head>"
    return head

def _generate_unclustered_table_html(
    table_id: str, df: pd.DataFrame
) -> str:
    data = df.to_json(orient='records')
    headers = [
        { 'key': "caprieval_rank", 'label': "Structure Rank", 'sorted': "asc" },
        { 'key': "model", 'label': "Structure", 'sortable': False, 'type': "structure" },
    ] + [
        {'key': k, 'label': v, 'type': 'stats'} for k, v in AXIS_NAMES.items()
    ] + [
        { 'key': "id", 'label': "Structure ID" },
    ]
    return f"""
            <div id="{table_id}"></div>
            <script type="module">
            import {{createRoot}} from "react-dom"
            import {{createElement}} from "react"
            import {{StructureTable}} from "@i-vresse/haddock3-analysis-components"

            const props = {{
                structures: {data},
                headers: {json.dumps(headers)}
            }}

            createRoot(document.getElementById('{table_id}')).render(
                createElement(StructureTable, props)
            )
            </script>"""

def _generate_clustered_table_html(
    table_id: str, df: pd.DataFrame
) -> str:
    data = df.to_json(orient='records')
    nr_best_columns = df.filter(like="best").shape[1]
    headers = [
        { 'key': "cluster_rank", 'label': "Cluster Rank", 'sorted': "asc" },
        { 'key': "cluster_id", 'label': "Cluster ID" },
        { 'key': "n", 'label': "Cluster size" },
    ] + [
        {'key': k, 'label': v, 'type': 'stats'} for k, v in AXIS_NAMES.items()
    ] + [
        { 'key': f"best{i}", 'label': f"Nr {i} best structure", 'sortable': False, 'type': "structure" } for i in range(1, nr_best_columns + 1)
    ]
    return f"""
            <div id="{table_id}"></div>
            <script type="module">
            import {{createRoot}} from "react-dom"
            import {{createElement}} from "react"
            import {{ClusterTable}} from "@i-vresse/haddock3-analysis-components"

            const clusters = {data}
            const headers = {json.dumps(headers)}

            createRoot(document.getElementById('{table_id}')).render(
                createElement(ClusterTable, {{ clusters, headers }})
            )
            </script>"""

def _generate_html_body(figures):
    """
    Generate an HTML body section containing figures for an analysis report.

    Parameters
    ----------
    figures : list
        A list of figures to include in the HTML body.
        Each figure can be either a string representing a table or a
        plotly.graph_objects.Figure object.

    Returns
    -------
    body : str
        The generated HTML body as a string.
    """
    body = "<body>"
    table_index = 1
    fig_index = 1
    for figure in figures:
        if isinstance(figure, pd.DataFrame):  # tables
            table_index += 1
            table_id = f"table{table_index}"

            is_unclustered = 'cluster_rank' not in figure
            if is_unclustered:
                inner_html = _generate_unclustered_table_html(
                    table_id, figure
                )
            else:
                inner_html = _generate_clustered_table_html(
                    table_id, figure
                )
        else:  # plots
            inner_json = figure.to_json()
            inner_html = create_html(inner_json, fig_index, figure.layout.height, figure.layout.width)
            fig_index += 1  # type: ignore
        body += "<br>"  # add a break between tables and plots
        body += inner_html
    body += "</body>"
    return body


def report_generator(boxes, scatters, tables, step):
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
    """
    figures = [tables]
    # Combine scatters
    figures.append(
        report_plots_handler(scatters, shared_xaxes="rows", shared_yaxes="columns")
    )
    # Combine boxes"
    figures.append(report_plots_handler(boxes))

    # Write everything to a html file
    html_report = _generate_html_report(step, figures)
    with open("report.html", "w", encoding="utf-8") as report:
        report.write(html_report)



def heatmap_plotly(
        matrix: NDFloat,
        labels: Optional[dict] = None,
        xlabels: Optional[list] = None,
        ylabels: Optional[list] = None,
        color_scale: str = 'Greys_r',  # Greys_r, gray
        title: Optional[str] = None,
        output_fname: Path = Path('contacts.html'),
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

    Return
    ------
    output_fname : Path
        Path to the generated filename
    """
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
    # Save figure as html file
    fig.write_html(output_fname)
    return output_fname

  
def make_alascan_plot(df, clt_id, scan_res="ALA"):
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
    fig = go.Figure(layout={"width": 2000, "height": 1000})
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
            title="Residue Name",
            tickfont_size=14,
            titlefont_size=16,
            tick0=df["full_resname"],
            # in case we want to show less residues
            # dtick=10,
            ),
        yaxis=dict(
            title="Weigted delta",
            titlefont_size=16,
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
    fig.write_html(html_output_filename)

