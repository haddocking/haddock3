"""Plotting functionalities."""

import numpy as np
import pandas as pd
import plotly.colors as px_colors
import plotly.express as px
import plotly.graph_objects as go

from haddock import log


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
    ]

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
    "fnat": "FCC"
    }

AXIS_NAMES = {
    "score": "HADDOCK score [a.u.]",
    "irmsd": "interface RMSD [A]",
    "lrmsd": "ligand RMSD [A]",
    "ilrmsd": "interface-ligand RMSD [A]",
    "desolv": "Desolvation Energy",
    "vdw": "Van der Waals Energy",
    "elec": "Electrostatic Energy",
    "air": "Restraints Energy",
    "fnat": "Fraction of Common Contacts",
    "dockq": "DOCKQ",
    }


def read_capri_table(capri_filename, comment="#"):
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
    capri_df = pd.read_csv(
        capri_filename,
        sep="\t",
        comment=comment)
    return capri_df


def in_capri(column, df_columns):
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


def update_layout_plotly(fig, x_label, y_label, title=None):
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
        "hoverlabel": dict(font_size=16, font_family="Helvetica")
        }
    fig.update_layout(px_dict)
    return


def box_plot_plotly(gb_full, y_ax, format, scale):
    """
    Create a scatter plot in plotly.
    
    Parameters
    ----------
    gb_full : pandas DataFrame
        data to box plot
    y_ax : str
        variable to plot
    format : str
        Produce images in the selected format.
    scale : int
        scale of image
    """
    fig = px.box(gb_full,
                 x="capri_rank",
                 y=f"{y_ax}",
                 color="cluster-id",
                 boxmode="overlay",
                 points="outliers",
                 width=1000,
                 height=800,
                 )
    # layout
    update_layout_plotly(fig, "Cluster rank", AXIS_NAMES[y_ax])
    # save figure
    px_fname = f"{y_ax}_clt.html"
    fig.write_html(px_fname, full_html=False, include_plotlyjs='cdn')
    # create format boxplot if necessary
    if format:
        fig.write_image(f"{y_ax}_clt.{format}", scale=scale)
    return


def box_plot_data(capri_df, cl_rank):
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
            cl_df["capri_rank"] = cl_rank[cl_id]
            gb_good = pd.concat([gb_good, cl_df])
            
    gb_other["cluster-id"] = "Other"
    gb_other["capri_rank"] = len(cl_rank.keys()) + 1
    gb_full = pd.concat([gb_good, gb_other])
    return gb_full


def box_plot_handler(capri_filename, cl_rank, format, scale):
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
    for y_ax in AXIS_NAMES.keys():
        if not in_capri(y_ax, capri_df.columns):
            continue
        box_plot_plotly(gb_full, y_ax, format, scale)
    return


def scatter_plot_plotly(gb_cluster, gb_other, cl_rank, x_ax, y_ax, colors, format, scale):  # noqa:E501
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
    """
    fig = go.Figure(layout={"width": 1000, "height": 800})
    traces = []
    n_colors = len(colors)
    for cl_id, cl_df in gb_cluster:
        if cl_id in cl_rank.keys():
            if cl_id == "-":
                cl_name = "Unclustered"
            else:
                cl_name = f"Cluster {cl_id}"
            color_idx = (cl_rank[cl_id] - 1) % n_colors  # color index
            x_mean = np.mean(cl_df[x_ax])
            y_mean = np.mean(cl_df[y_ax])
            text_list = [f"Model: {cl_df['model'].iloc[n].split('/')[-1]}<br>Score: {cl_df['score'].iloc[n]}" for n in range(cl_df.shape[0])]  # noqa:E501
            traces.append(
                go.Scatter(
                    x=cl_df[x_ax],
                    y=cl_df[y_ax],
                    name=cl_name,
                    mode="markers",
                    text=text_list,
                    legendgroup=cl_name,
                    marker_color=colors[color_idx],
                    hoverlabel=dict(
                        bgcolor=colors[color_idx],
                        font_size=16,
                        font_family="Helvetica"
                        )
                    )
                )
            clt_text = f"{cl_name}<br>"
            if 'score' not in [x_ax, y_ax]:
                clt_text += f"Score: {np.mean(cl_df['score']):.3f}<br>"
            clt_text += f"{x_ax}: {x_mean:.3f}<br>{y_ax}: {y_mean:.3f}"
            clt_text_list = [clt_text]
            traces.append(
                go.Scatter(
                    x=[x_mean],
                    y=[y_mean],
                    # error bars
                    error_x=dict(
                        type='data',
                        array=[np.std(cl_df[x_ax])],
                        visible=True),
                    error_y=dict(
                        type='data',
                        array=[np.std(cl_df[y_ax])],
                        visible=True),
                    # color and text
                    marker_color=colors[color_idx],
                    text=clt_text_list,
                    legendgroup=cl_name,
                    showlegend=False,
                    mode="markers",
                    marker=dict(
                        size=10,
                        symbol="square-dot"),
                    hovertemplate=f'<b>{clt_text}</b><extra></extra>',
                    hoverlabel=dict(
                        bgcolor=colors[color_idx],
                        font_size=16,
                        font_family="Helvetica"
                        )
                    )
                )
    # append trace other
    if not gb_other.empty:
        text_list_other = [f"Model: {gb_other['model'].iloc[n].split('/')[-1]}<br>Score: {gb_other['score'].iloc[n]}" for n in range(gb_other.shape[0])]  # noqa:E501
        traces.append(
            go.Scatter(
                x=gb_other[x_ax],
                y=gb_other[y_ax],
                name="Other",
                mode="markers",
                text=text_list_other,
                legendgroup="Other",
                marker=dict(
                    color="white",
                    line=dict(
                        width=2,
                        color='DarkSlateGrey')
                    ),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=16,
                    font_family="Helvetica"
                    )
                )
            )
    for trace in traces:
        fig.add_trace(trace)
    px_fname = f"{x_ax}_{y_ax}.html"
    update_layout_plotly(fig,
                         TITLE_NAMES[x_ax],
                         TITLE_NAMES[y_ax],
                         title=f"{TITLE_NAMES[x_ax]} vs {TITLE_NAMES[y_ax]}")
    fig.write_html(px_fname, full_html=False, include_plotlyjs='cdn')
    # create format boxplot if necessary
    if format:
        fig.write_image(f"{x_ax}_{y_ax}.{format}", scale=scale)
    return


def scatter_plot_data(capri_df, cl_rank):
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


def scatter_plot_handler(capri_filename, cl_rank, format, scale):
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
    """
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_cluster, gb_other = scatter_plot_data(capri_df, cl_rank)
    
    # defining colors
    colors = px_colors.qualitative.Alphabet
    for x_ax, y_ax in SCATTER_PAIRS:
        if not in_capri(x_ax, capri_df.columns):
            continue
        if not in_capri(y_ax, capri_df.columns):
            continue
        scatter_plot_plotly(gb_cluster,
                            gb_other,
                            cl_rank,
                            x_ax,
                            y_ax,
                            colors,
                            format,
                            scale)
    return
