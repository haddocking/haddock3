"""Plotting functionalities."""

from pathlib import Path
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

# if SCATTER_PAIRS changes, SCATTER_MATRIX_SIZE should change too!
SCATTER_MATRIX_SIZE = (4, 5)  # (number of rows, number of columns)

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
    capri_df = pd.read_csv(capri_filename, sep="\t", comment=comment)
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
        "hoverlabel": dict(font_size=16, font_family="Helvetica"),
        }
    fig.update_layout(px_dict)
    return fig


def box_plot_plotly(gb_full, y_ax, cl_rank, format, scale):
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
    colors = px_colors.qualitative.Alphabet
    color_map = {}
    for cl_id in sorted(cl_rank.keys()):
        color_idx = (cl_rank[cl_id] - 1) % len(colors)  # color index
        color_map[f"{cl_id}"] = colors[color_idx]
    # to use color_discrete_map, cluster-id column should be str not int
    gb_full_string = gb_full.astype({"cluster-id": "string"})

    fig = px.box(gb_full_string,
                 x="capri_rank",
                 y=f"{y_ax}",
                 color="cluster-id",
                 color_discrete_map=color_map,
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
    return fig


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
    fig_list = []
    for y_ax in AXIS_NAMES.keys():
        if not in_capri(y_ax, capri_df.columns):
            continue
        fig = box_plot_plotly(gb_full, y_ax, cl_rank, format, scale)
        fig_list.append(fig)
    return fig_list


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

    Returns
    -------
    fig :
        an instance of plotly.graph_objects.Figure
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
            # refactored due to E501 line too long
            text_list = []
            for n in range(cl_df.shape[0]):
                model_text = cl_df['model'].iloc[n].split('/')[-1]
                score_text = cl_df['score'].iloc[n]
                text_list.append(f"Model: {model_text}<br>Score: {score_text}")
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
                        font_family="Helvetica",
                        ),
                    )
                )
            clt_text = f"{cl_name}<br>"
            if "score" not in [x_ax, y_ax]:
                clt_text += f"Score: {np.mean(cl_df['score']):.3f}<br>"
            clt_text += f"{x_ax}: {x_mean:.3f}<br>{y_ax}: {y_mean:.3f}"
            clt_text_list = [clt_text]
            traces.append(
                go.Scatter(
                    x=[x_mean],
                    y=[y_mean],
                    # error bars
                    error_x=dict(
                        type="data", array=[np.std(cl_df[x_ax])], visible=True
                        ),
                    error_y=dict(
                        type="data", array=[np.std(cl_df[y_ax])], visible=True
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
        # refactored due to E501 line too long
        text_list_other = []
        for n in range(gb_other.shape[0]):
            model_text = gb_other['model'].iloc[n].split('/')[-1]
            score_text = gb_other['score'].iloc[n]
            text_list_other.append(
                f"Model: {model_text}<br>Score: {score_text}"
                )
        traces.append(
            go.Scatter(
                x=gb_other[x_ax],
                y=gb_other[y_ax],
                name="Other",
                mode="markers",
                text=text_list_other,
                legendgroup="Other",
                marker=dict(
                    color="white", line=dict(width=2, color="DarkSlateGrey")
                    ),
                hoverlabel=dict(
                    bgcolor="white", font_size=16, font_family="Helvetica",
                    ),
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
    return fig


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

    Returns
    -------
    fig_list : list
        a list of figures
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
                                  colors,
                                  format,
                                  scale)
        fig_list.append(fig)
    return fig_list


def _report_grid_size(plot_list):
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


def report_plots_handler(
        plots,
        shared_xaxes=False,
        shared_yaxes=False
        ):
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
        a paramtere of plotly.subplots.make_subplots
    shared_yaxes: boolean or str (default False)
        a paramtere of plotly.subplots.make_subplots

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


def find_best_struct(ss_file, number_of_struct=10):
    """
    Find best structures.

    It inspects model-cluster-ranking recorded in capri_ss.tsv file and finds
    the best models (models with lower ranks).
    By default, it selects the 10 best models.

    Parameters
    ----------
    ss_file : path
        path to capri_ss.tsv
    number_of_struct: int
        number of models with lower model-cluster-ranking

    Returns
    -------
    best_struct_df : pandas DataFrame
        DataFrame of best structures
    """
    dfss = read_capri_table(ss_file)
    dfss = dfss.sort_values(by=["cluster-id", "model-cluster-ranking"])
    # TODO need a check for "Unclustered"

    # count values within each cluster
    # and select the column model-cluster-ranking
    dfss_grouped = dfss.groupby("cluster-id").count()["model-cluster-ranking"]

    # number of structs can be different per each cluster,
    # so min value is picked here
    max_number_of_struct = dfss_grouped.min()

    # number_of_struct cannot be greater than max_number_of_struct
    number_of_struct = min(number_of_struct, max_number_of_struct)

    # select the best `number_of_struct` e.g. 4 structures for each cluster
    best_struct_df = dfss.groupby("cluster-id").head(number_of_struct).copy()

    # define names for best structures, e.g.
    # Nr 1 best structure, Nr 2 best structure, ...
    number_of_cluster = len(best_struct_df["cluster-id"].unique())
    # zero pad number so after pivot columns are sorted correctly
    col_names = [
        f"Nr {(number + 1):02d} best structure" for number in range(number_of_struct)
        ] * number_of_cluster

    # add a new column `Structure` to the dataframe
    best_struct_df = best_struct_df.assign(Structure=col_names)

    # reshape data frame where columns are
    # cluster-id, model,.., Nr 1 best structure, Nr 2 best structure, ...
    best_struct_df = best_struct_df.pivot_table(
        index=["cluster-id"],
        columns=["Structure"],
        values="model",
        aggfunc=lambda x: x,
        )
    best_struct_df.reset_index(inplace=True)
    best_struct_df.rename(columns={"cluster-id": "Cluster ID"}, inplace=True)
    return best_struct_df


def _ngl_viewer():
    # http://nglviewer.org/ngl/gallery/
    return f"""<script src="https://cdn.rawgit.com/arose/ngl/v2.1.0/dist/ngl.js"></script>
            <script>
                var dialog;
                var stage;

                document.addEventListener("DOMContentLoaded", function () {{
                    // Setup to load data from rawgit
                    NGL.DatasourceRegistry.add(
                        "data", new NGL.StaticDatasource( "//cdn.rawgit.com/arose/ngl/v2.1.0/data/" )
                    );

                    // Create NGL Stage object
                    stage = new NGL.Stage( "viewport" );

                    // Handle window resizing
                    window.addEventListener( "resize", function( event ){{
                        stage.handleResize();
                    }}, false );

                }});

                function showStructure(file_name) {{
                    dialog = document.getElementById("structureViewerDialog");
                    dialog.showModal();
                    stage.loadFile(file_name).then(function (o) {{
                        o.addRepresentation("cartoon");
                        o.autoView();
                        stage.handleResize();
                    }});
                }}
            </script>
            <body>
                <dialog id="structureViewerDialog">
                    <div id="viewport" style="width:800px; height:600px;"></div>
                    <form>
                        <button value="cancel" formmethod="dialog">X</button>
                    </form>
                </dialog>
            </body>"""


def _add_viewers(df):
    def format_cell(row):
        pdb_file = Path(row)
        correct_path = Path(f"{row}.gz") if not pdb_file.exists() else pdb_file
        return f'''\
            <p>&#8595;&nbsp;<a href="{str(correct_path)}">Download</a>&nbsp;\
            <a onClick="showStructure('{str(correct_path)}')"
            style="cursor:pointer;">&#x1F441; View</a></p>\
            '''

    table_df = df.copy()
    for col_name in table_df.columns[1:]:
        table_df[col_name] = df[col_name].apply(format_cell)
    return table_df


def clean_capri_table(dfcl):
    """
    Craete a tidy capri table for the report.

    It changes the clomuns names of capri table (a dataframe that is read by
    read_capri_table). It also combines mean and std values in one column.

    Parameters
    ----------
    dfcl : pandas DataFrame
        dataframe of capri values

    Returns
    -------
    dfcl : pandas DataFrame
        DataFrame of capri table with new column names
    """
    dfcl = dfcl.sort_values(by=["score"])
    # what metrics are in both dfcl and AXIS_NAMES
    col_list = dfcl.columns.intersection(list(AXIS_NAMES.keys())).tolist()
    # columns of the final table
    table_col = ["Cluster ID", "Cluster size"]
    for col_name in col_list:
        mean_value = dfcl[col_name].astype(str)
        std_value = dfcl[f"{col_name}_std"].astype(str)
        dfcl[AXIS_NAMES[col_name]] = (mean_value + " &#177; " + std_value)
        table_col.append(AXIS_NAMES[col_name])
    dfcl.drop(columns=col_list, inplace=True)
    dfcl.rename(
        columns={"cluster_id": "Cluster ID", "n": "Cluster size"},
        inplace=True,
        )
    return dfcl[table_col]


def _pandas_df_to_html_table(df, table_id):
    return df.to_html(
        table_id=table_id,
        index= False,
        index_names=False,
        classes="table",
        escape=False,
        )


def clt_table_handler(clt_file, ss_file):
    """
    Create a figure include tables.

    The idea is to create tidy tables that report statistics available in
    capri_clt.tsv and capri_ss.tsv files.

    Parameters
    ----------
    clt_file : str or Path
        path to capri_clt.tsv file
    ss_file: str or Path
        path to capri_ss.tsv file

    Returns
    -------
    fig :
        an instance of plotly.graph_objects.Figure
    """
    # TODO be able to sort inside the table
    # TODO add some of the ngl tools e.g. rotate, show water
    # table of statistics
    dfcl = read_capri_table(clt_file)
    statistics_df = clean_capri_table(dfcl)
    # Convert dataframe to html table
    statistics_table = _pandas_df_to_html_table(statistics_df, table_id="statistics")

    # table of structures
    structs_df = find_best_struct(ss_file, number_of_struct=10)

    # Order structs by best (lowest score) cluster on top
    structs_df = structs_df.set_index('Cluster ID')
    structs_df = structs_df.reindex(index=statistics_df['Cluster ID'])
    structs_df = structs_df.reset_index()

    # Add download links and viewer
    struct_df = _add_viewers(structs_df)
    # Convert dataframe to html table
    struct_table = _pandas_df_to_html_table(struct_df, table_id="structures")

    return [statistics_table, struct_table]


def _css_styles_for_report():
    custom_css = '''
    .table {
        font-family: Arial, sans-serif;
        border-collapse: collapse;
        width: 100%;
        font-size: 12px;
        }
    .table th {
        font-weight: bold;
        background-color: #f2f2f2;
        padding: 8px;
        border: 1px solid #ddd;
        text-align: left;
        }
    .table td {
        border: 1px solid #ddd;
        padding: 8px;
        text-align: left;
        }
    .table tr:nth-child(even) {
        background-color: #f2f2f2;
        }
    .title {
        font-family: Arial, sans-serif;
        font-size: 16px;
        }

    '''
    return f'<style>{custom_css}</style>'


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
    figures = tables
    # Combine scatters
    figures.append(
        report_plots_handler(
            scatters, shared_xaxes="rows", shared_yaxes="columns"
            )
        )
    # Combine boxes"
    figures.append(report_plots_handler(boxes, shared_xaxes="all"))
    report_filename = "report.html"
    with open(report_filename, "w", encoding="utf-8") as report:
        # TODO move html code to a function
        report.write("<!DOCTYPE html><html lang='en'><head>")
        report.write(f"<title>Analysis report of step {step}</title>")
        report.write(f"<p class='title'>Analysis report of step {step}</p>")
        report.write(_css_styles_for_report())
        report.write("</head><body>")
        include_plotlyjs = "cdn"
        for figure in figures:
            if isinstance(figure, str): # tables
                inner_html = figure
            else: # plots
                inner_html = figure.to_html(
                    full_html=False, include_plotlyjs=include_plotlyjs
                    )
                include_plotlyjs = False
            report.write("<br>")
            report.write(inner_html)
        report.write(_ngl_viewer())
        report.write("</body></html>")
