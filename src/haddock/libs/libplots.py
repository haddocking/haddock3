"""Plotting functionalities."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.colors as px_colors
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
    return fig


def update_layout_matplotlib(x_label, y_label):
    """
    Update layout of matplotlib plot.

    Parameters
    ----------
    x_label : str
        x axis name
    y_label : str
        y axis name
    """
    plt.title(f"{x_label} vs {y_label}", fontsize=40)
    plt.xlabel(x_label, fontsize=40)
    plt.ylabel(y_label, fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    legend = plt.legend(bbox_to_anchor=(1.0, 1.0), fontsize=30)
    for handle in legend.legendHandles:
        handle.set_sizes([150])
    plt.tight_layout()
    return


def box_plot_plotly(gb_full, y_ax, cl_rank):
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
    """
    colors = px_colors.qualitative.Alphabet
    color_map = {}
    for cl_id in sorted(cl_rank.keys()):
        color_idx = (cl_rank[cl_id] - 1) % len(colors)  # color index
        color_map[f"{cl_id}"] = colors[color_idx]
    # to use color_discrete_map, cluster-id column should be str not int
    gb_full_string= gb_full.astype({"cluster-id":"string"})
    fig = px.box(gb_full_string,
                 x="capri_rank",
                 y=f"{y_ax}",
                 color="cluster-id",
                 color_discrete_map=color_map,
                 boxmode="overlay",
                 points="outliers"
                 )
    # layout
    update_layout_plotly(fig, "Cluster rank", AXIS_NAMES[y_ax])
    # save figure
    px_fname = f"{y_ax}_clt.html"
    fig.write_html(px_fname, full_html=False, include_plotlyjs='cdn')
    return fig


def box_plot_matplotlib(gbcl, y_ax, legend_labels, dpi):
    """
    Create a scatter plot in plotly.

    Parameters
    ----------
    gbcl : pandas DataFrameGroupBy
        capri data grouped by cluster-id
    y_ax : str
        variable to plot
    legend_labels : list
        list of legend labels
    dpi : int
        DPI for png images.

    """
    plt.figure(figsize=(20, 16))
    fig, ax = plt.subplots()
    # box plot must take element 1 of groupby instance
    bplot = ax.boxplot([el[1][y_ax] for el in gbcl],
                       patch_artist=True,
                       showmeans=False,
                       )
    # assigning colors to boxes and medians
    colors = plt.cm.tab20.colors
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    for median, color in zip(bplot['medians'], colors):
        median.set_color(color)

    # plot details
    ax.legend(bplot["boxes"], legend_labels, bbox_to_anchor=(1.0, 1.0))
    plt.ylabel(AXIS_NAMES[y_ax])
    plt.xlabel("Cluster rank")
    plt.tight_layout()
    # save figure
    plt_fname = f"{y_ax}_clt.png"
    plt.savefig(fname=plt_fname, dpi=dpi)
    plt.close()
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
    gb_other : pandas DataFrame
        DataFrame of clusters not in the top cluster ranking
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
    return gb_full, gb_other


def box_plot_handler(capri_filename, cl_rank, png, dpi):
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
    png : bool
        Produce png images.
    dpi : int
        DPI for png images.
    """
    # generating the correct dataframe
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_full, gb_other = box_plot_data(capri_df, cl_rank)

    if png:
        # when png is true we have to manually sort according to capri_rank
        gb_sorted = gb_full.sort_values("capri_rank")
        gbcl = gb_sorted.groupby("cluster-id", sort=False)
        legend_labels = [f"Cluster {el[0]}" for el in gbcl]
        if not gb_other.empty:
            legend_labels[-1] = "Other"

    # iterate over the variables
    fig_list = []
    for y_ax in AXIS_NAMES.keys():
        if not in_capri(y_ax, capri_df.columns):
            continue
        fig = box_plot_plotly(gb_full, y_ax, cl_rank)
        fig_list.append(fig)
        # create png boxplot if necessary
        if png:
            box_plot_matplotlib(gbcl, y_ax, legend_labels, dpi)
    return fig_list


def scatter_plot_plotly(gb_cluster, gb_other, cl_rank, x_ax, y_ax, colors):
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
    return fig


def scatter_plot_matplotlib(gbcl, gb_other, cl_rank, x_ax, y_ax, colors, dpi):
    """
    Create a scatter plot in matplotlib.

    Parameters
    ----------
    gbcl : pandas DataFrameGroupBy
        capri DataFrame grouped by cluster-id.
    gb_other : pandas DataFrame
        DataFrame of clusters not in the top cluster ranking.
    cl_rank : dict
        {cluster_id : cluster_rank} dictionary.
    x_ax : str
        name of the x column.
    y_ax : str
        name of the y column.
    colors : list
        list of colors to be used.
    dpi : int
        DPI for png images.
    """
    plt.figure(figsize=(20, 16))
    for cl_id, cl_df in gbcl:
        if cl_id not in cl_rank.keys():
            gb_other = pd.concat([gb_other, cl_df])
        else:
            if cl_id == "-":
                cl_name = "Unclustered"
            else:
                cl_name = f"Cluster {cl_id}"
            color_idx = (cl_rank[cl_id] - 1) % len(colors)  # color index
            plt.scatter(x=cl_df[x_ax],
                        y=cl_df[y_ax],
                        color=colors[color_idx],
                        label=cl_name,
                        s=40)
    update_layout_matplotlib(TITLE_NAMES[x_ax], TITLE_NAMES[y_ax])
    plt_fname = f"{x_ax}_{y_ax}.png"
    plt.savefig(plt_fname, dpi=dpi)
    # creating full picture
    if not gb_other.empty:
        plt.scatter(
            x=gb_other[x_ax],
            y=gb_other[y_ax],
            color="gray",
            label="Other",
            s=10,
            )
        # details
        update_layout_matplotlib(TITLE_NAMES[x_ax], TITLE_NAMES[y_ax])
        plt_fname = f"{x_ax}_{y_ax}_full.png"
        plt.savefig(plt_fname, dpi=dpi)
    plt.close()
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


def scatter_plot_handler(capri_filename, cl_rank, png, dpi):
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
    png : bool
        Produce png images.
    dpi : int
        DPI for png images.
    """
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_cluster, gb_other = scatter_plot_data(capri_df, cl_rank)

    # defining colors
    colors = px_colors.qualitative.Alphabet
    fig_list = []
    for x_ax, y_ax in SCATTER_PAIRS:
        if not in_capri(x_ax, capri_df.columns):
            continue
        if not in_capri(y_ax, capri_df.columns):
            continue
        fig = scatter_plot_plotly(gb_cluster,
                                gb_other,
                                cl_rank,
                                x_ax,
                                y_ax,
                                colors)
        fig_list.append(fig)
        if png:
            scatter_plot_matplotlib(gb_cluster,
                                    gb_other,
                                    cl_rank,
                                    x_ax,
                                    y_ax,
                                    colors,
                                    dpi)
    return fig_list


def _report_grid(plot_list):
    # Calculate grid size for subplots
    number_of_plots = len(plot_list)
    number_of_clusters = len({trace.legendgroup for trace in plot_list[0].data})
    number_of_cols = 5
    number_of_rows = int(np.ceil(number_of_plots / number_of_cols))
    if number_of_clusters > 5:
        number_of_cols = 2
        number_of_rows = int(np.ceil(number_of_plots / number_of_cols))
    return number_of_rows, number_of_cols


def report_plots_handler(plots, plot_title, shared_xaxes=False):
    number_of_rows, number_of_cols = _report_grid(plots)
    fig = make_subplots(
        rows=number_of_rows,
        cols=number_of_cols,
        shared_xaxes=shared_xaxes,
        vertical_spacing=(0.4 / number_of_rows),
        horizontal_spacing=(0.3 / number_of_cols),
        )
    for i, sub_fig in enumerate(plots):
        col_index = int((i % number_of_cols) + 1)
        row_index = int(np.floor(i / number_of_cols) + 1)
        # hide legend of plots except one
        if i !=0:
            sub_fig.for_each_trace(lambda trace:trace.update(showlegend=False))
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
        title_text=plot_title,
        legend_title_text = legend_title_text,
        height=900)
    return fig


def find_best_struct(ss_file, number_of_struct=4):
    dfss = read_capri_table(ss_file)
    dfss = dfss.sort_values(by=["cluster-id", "model-cluster-ranking"])
    best_struct_df = dfss.groupby("cluster-id").head(number_of_struct).copy()
    number_of_cluster = len(best_struct_df["cluster-id"].unique())
    col_names = [
        f"Nr {number + 1} best structure" for number in range(number_of_struct)
        ] * number_of_cluster
    best_struct_df = best_struct_df.assign(Structure=col_names)
    best_struct_df = best_struct_df.pivot_table(
        index=["cluster-id"], columns=['Structure'], values="model", aggfunc=lambda x:x,
        )
    best_struct_df.reset_index(inplace=True)
    best_struct_df.rename(columns={"cluster-id":"Cluster ID"}, inplace=True)
    return best_struct_df


def clean_capri_table(dfcl):
    dfcl = dfcl.sort_values(by=["cluster_id"])
    # what metrics are in both dfcl and AXIS_NAMES
    col_list = dfcl.columns.intersection(list(AXIS_NAMES.keys())).tolist()
    # columns of the final table
    table_col = ["Cluster ID", "Cluster size"]
    for col_name in col_list:
        dfcl[AXIS_NAMES[col_name]] = dfcl[col_name].astype(str) + " ± " + dfcl[f"{col_name}_std"].astype(str)
        table_col.append(AXIS_NAMES[col_name])
    dfcl.drop(columns=col_list, inplace=True)
    dfcl.rename(columns={"cluster_id":"Cluster ID", "n":"Cluster size"}, inplace=True)
    return dfcl[table_col]


def clt_table_handler(clt_file, ss_file):
    dfcl = read_capri_table(clt_file)
    statistics_df = clean_capri_table(dfcl)
    structs_df = find_best_struct(ss_file, number_of_struct=4)
    fig = make_subplots(
        rows=2,
        cols=1,
        specs=[[{"type": "table"}], [{"type": "table"}]],
        vertical_spacing=0.1,
        )
    for i, df in enumerate([statistics_df, structs_df]):
        table = go.Table(
            columnwidth = 200,
            header=dict(values=list(df.columns),
            align='left'),
            cells=dict(values=df.transpose().values.tolist(),
            # fill_color='lavender',
            align='left',
            height=40),
            )
        fig.add_trace(table, row=i+1, col=1)
    fig.update_layout(
        title_text="Summary",
        height=700)
    fig.write_html("clt_table.html", full_html=False, include_plotlyjs='cdn')
    return fig


def report_generator(boxes, scatters, table, step):
    figures = [table]
    # Combine scatters
    # TODO fix number of rows and share axes for scatters with horizontal scrolling
    plot_title = f"Scatter plots of {step}"
    figures.append(report_plots_handler(scatters, plot_title))
    # Combine boxes
    plot_title = f"Box plots of {step}"
    figures.append(report_plots_handler(boxes, plot_title, shared_xaxes="all"))
    report_filename = "report.html"
    report = open(report_filename, 'w')
    report.write("<html><head></head><body>" + "\n")
    include_plotlyjs = 'cdn'
    for figure in figures:
        inner_html = figure.to_html(full_html=False, include_plotlyjs=include_plotlyjs)
        report.write(inner_html)
        include_plotlyjs = False
    report.write("</body></html>" + "\n")