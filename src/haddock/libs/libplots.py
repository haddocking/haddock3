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


def box_plots(capri_filename, cl_ranking):
    """
    Create box plots.

    The idea is that for each of the top X-ranked clusters we create a box plot
    showing how the basic statistics are distributed within each model.
    
    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    cl_ranking : dict
        {cluster_id : cluster_rank} dictionary
    """
    # generating the correct dataframe
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_cluster = capri_df.groupby("cluster-id")
    gb_other = pd.DataFrame([])
    gb_good = pd.DataFrame([])
    for cl_id, cl_df in gb_cluster:
        if cl_id not in cl_ranking.keys():
            gb_other = pd.concat([gb_other, cl_df])
        else:
            cl_df["capri_rank"] = cl_ranking[cl_id]
            gb_good = pd.concat([gb_good, cl_df])
    gb_other["cluster-id"] = "Other"
    gb_other["capri_rank"] = len(cl_ranking.keys()) + 1
    gb_cluster = pd.concat([gb_good, gb_other])
    
    for x_ax in AXIS_NAMES.keys():
        if x_ax not in capri_df.columns:
            log.info(f"x axis quantity {x_ax} not present in capri table")
            continue
        fig = px.box(gb_cluster,
                     x="capri_rank",
                     y=f"{x_ax}",
                     color="cluster-id",
                     boxmode="overlay",
                     points="outliers"
                     )
        # layout
        px_fname = f"{x_ax}_clt.html"
        fig.update_layout(
            xaxis=dict(
                title="Cluster rank",
                tickfont_size=14,
                titlefont_size=40,
                ),
            yaxis=dict(
                title=AXIS_NAMES[x_ax],
                tickfont_size=14,
                titlefont_size=40,
                ),
            legend=dict(x=1.01, y=1.0, font_family="Helvetica", font_size=16),
            hoverlabel=dict(font_size=16, font_family="Helvetica"),
            )

        fig.write_html(px_fname, full_html=False, include_plotlyjs='cdn')
        

def scatter_plots(capri_filename, cl_ranking):
    """
    Create scatter plots.

    The idea is that for each pair of variables of interest (SCATTER_PAIRS,
     declared as global) we create a scatter plot.
    If available, each scatter plot containts cluster information.

    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    cl_ranking : dict
        {cluster_id : cluster_rank} dictionary
    """
    capri_df = read_capri_table(capri_filename, comment="#")
    for x_ax, y_ax in SCATTER_PAIRS:
        if x_ax not in capri_df.columns:
            log.warning(f"x axis quantity {x_ax} not present in capri table")
            continue
        if y_ax not in capri_df.columns:
            log.warning(f"y axis quantity {y_ax} not present in capri table")
            continue
        
        gb_cluster = capri_df.groupby("cluster-id")
        gb_other = pd.DataFrame([])
        fig = go.Figure(layout={"width": 1000, "height": 800})
        traces = []
        # defining colors
        colors = px_colors.qualitative.Alphabet
        n_colors = len(colors)
        for cl_id, cl_df in gb_cluster:
            if cl_id not in cl_ranking.keys():
                gb_other = pd.concat([gb_other, cl_df])
            else:
                if cl_id == "-":
                    cl_name = "Unclustered"
                else:
                    cl_name = f"Cluster {cl_id}"
                color_idx = (cl_ranking[cl_id] - 1) % n_colors  # color index
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

        fig.update_layout(
            title=f"{TITLE_NAMES[x_ax]} vs {TITLE_NAMES[y_ax]}",
            xaxis=dict(
                title=AXIS_NAMES[x_ax],
                tickfont_size=14,
                titlefont_size=40,
                ),
            yaxis=dict(
                title=AXIS_NAMES[y_ax],
                tickfont_size=14,
                titlefont_size=40,
                ),
            legend=dict(x=1.01, y=1.0, font_family="Helvetica", font_size=16),
            hoverlabel=dict(font_size=16, font_family="Helvetica"),
            )

        fig.write_html(px_fname, full_html=False, include_plotlyjs='cdn')
