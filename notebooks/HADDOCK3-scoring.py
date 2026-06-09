import marimo

__generated_with = "0.23.9"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    from pathlib import Path as _Path
    work_dir_input = mo.ui.text(
        value=str(_Path.cwd()),
        label="Working directory",
        full_width=True,
    )
    return (work_dir_input,)


@app.cell
def _(mo):
    mo.md(r"""
    # HADDOCK3 Scoring Workflow

    Upload one or more PDB complex files to score and cluster them with HADDOCK3.
    Results include per-model and per-cluster HADDOCK statistics displayed as
    interactive sortable tables, 3D visualisation of the best model from each
    cluster, and a contact-map chord chart for each cluster.

    **Workflow steps:**

    1. **`topoaa`** — builds CNS all-atom topology and parameters
    2. **`emscoring`** — energy minimization scoring via CNS
    3. **`clustfcc`** — FCC-based clustering of scored models
    4. **`caprieval`** — per-model and per-cluster HADDOCK statistics
    5. **`contactmap`** — contact-map chord chart and heatmap for the top models

    > **Note:** this notebook supports standard protein, DNA/RNA, and ion inputs only.
    > Small-molecule ligands requiring custom topology/parameter files are **not** supported.
    > Meaningful clustering requires at least a few models; a single uploaded file will
    > produce one unclustered model.

    **Requirements:** CNS must be installed and discoverable (`haddock3-cfg cns_exec`).
    """)
    return


@app.cell
def _(mo):
    pdb_upload = mo.ui.file(
        filetypes=[".pdb"],
        label="Drop or select PDB complex file(s) — each file is one model of the complex",
        kind="area",
        multiple=True,
    )
    return (pdb_upload,)


@app.cell
def _(mo):
    import os as _os
    _ncores_max = _os.cpu_count() or 1
    ncores_slider = mo.ui.slider(
        1, _ncores_max, value=_ncores_max, step=1,
        label="CPU cores", show_value=True,
    )
    clust_cutoff_slider = mo.ui.slider(
        0.50, 1.00, value=0.75, step=0.05,
        label="FCC clustering cutoff", show_value=True,
    )
    min_population_slider = mo.ui.slider(
        1, 20, value=4, step=1,
        label="Min cluster population", show_value=True,
    )
    topX_slider = mo.ui.slider(
        1, 20, value=5, step=1,
        label="Top models for contact map", show_value=True,
    )
    chordchart_toggle = mo.ui.switch(label="Chord chart", value=True)
    heatmap_toggle = mo.ui.switch(label="Heatmap", value=True)
    return (
        chordchart_toggle,
        clust_cutoff_slider,
        heatmap_toggle,
        min_population_slider,
        ncores_slider,
        topX_slider,
    )


@app.cell
def _(
    chordchart_toggle,
    clust_cutoff_slider,
    heatmap_toggle,
    min_population_slider,
    mo,
    ncores_slider,
    pdb_upload,
    topX_slider,
    work_dir_input,
):
    mo.vstack([
        mo.md("### Configuration"),
        mo.hstack([
            mo.vstack([mo.md("**Input models**"), pdb_upload], align="start"),
            mo.vstack([mo.md("**Working directory**"), work_dir_input], align="start"),
            mo.vstack([mo.md("**Execution**"), ncores_slider], align="start"),
            mo.vstack([
                mo.md("**Clustering (clustfcc)**"),
                clust_cutoff_slider,
                min_population_slider,
            ], align="start"),
            mo.vstack([
                mo.md("**Contact map**"),
                topX_slider,
                chordchart_toggle,
                heatmap_toggle,
            ], align="start"),
        ], gap=2, justify="start"),
    ])
    return


@app.cell
def _(mo):
    run_btn = mo.ui.run_button(label="▶  Run HADDOCK3 Scoring Workflow", kind="success")
    run_btn
    return (run_btn,)


@app.cell
def _(
    chordchart_toggle,
    clust_cutoff_slider,
    heatmap_toggle,
    min_population_slider,
    mo,
    ncores_slider,
    pdb_upload,
    run_btn,
    topX_slider,
    work_dir_input,
):
    import logging
    import os
    import traceback
    from contextlib import contextmanager
    from datetime import datetime
    from pathlib import Path

    mo.stop(not run_btn.value)
    mo.stop(
        not pdb_upload.value,
        mo.callout(
            mo.md("**No PDB file loaded.** Upload at least one PDB complex above, then click Run."),
            kind="warn",
        ),
    )

    # ── Run directory ─────────────────────────────────────────────────────────
    _files = pdb_upload.value
    _stem = Path(_files[0].name).stem if len(_files) == 1 else "scoring"
    _timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    _base_dir = Path(work_dir_input.value)
    _base_dir.mkdir(parents=True, exist_ok=True)
    _run_dir = _base_dir / f"{_stem}_{_timestamp}"
    _run_dir.mkdir(exist_ok=True)

    _pdb_paths = []
    for _uf in _files:
        _dest = _run_dir / _uf.name
        _dest.write_bytes(_uf.contents)
        _pdb_paths.append(Path(_uf.name))

    _ncores = ncores_slider.value

    @contextmanager
    def _cwd(path):
        _prev = Path.cwd()
        os.chdir(path)
        try:
            yield
        finally:
            os.chdir(_prev)

    # ── Live log ───────────────────────────────────────────────────────────────
    _log_lines = []

    def _render_log(done=False):
        _label = "✅ Completed" if done else "⏳ Running…"
        _body = "\n".join(_log_lines) if _log_lines else "(waiting for output…)"
        _panel = mo.Html(
            '<div style="height:400px;overflow-y:auto;background:#f0f0f0;'
            'color:#1a1a1a;font-family:monospace;padding:10px 14px;'
            'border-radius:6px;font-size:12px;white-space:pre-wrap;line-height:1.5;">'
            + _body.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            + "</div>"
        )
        mo.output.replace(mo.accordion({f"HADDOCK3 Log — {_label}": _panel}))

    class _LogHandler(logging.Handler):
        def emit(self, record):
            _log_lines.append(self.format(record))
            _render_log(done=False)

    _handler = _LogHandler()
    _handler.setFormatter(logging.Formatter("[%(asctime)s %(module)s %(levelname)s] %(message)s"))
    _render_log()

    # ── Workflow parameters ───────────────────────────────────────────────────
    _workflow_params = {
        "topoaa.1": {
            "molecules": _pdb_paths,
            "autohis": True,
            "mode": "local",
            "ncores": _ncores,
            "clean": False,
            "offline": False,
        },
        "emscoring.1": {
            "mode": "local",
            "ncores": _ncores,
            "clean": False,
            "offline": False,
        },
        "clustfcc.1": {
            "clust_cutoff": clust_cutoff_slider.value,
            "min_population": min_population_slider.value,
            "mode": "local",
            "ncores": _ncores,
            "clean": False,
            "offline": False,
        },
        "caprieval.1": {
            "mode": "local",
            "ncores": _ncores,
            "clean": False,
            "offline": False,
        },
        "contactmap.1": {
            "generate_chordchart": chordchart_toggle.value,
            "generate_heatmap": heatmap_toggle.value,
            "single_model_analysis": False,
            "topX": topX_slider.value,
            "mode": "local",
            "ncores": _ncores,
            "clean": False,
            "offline": False,
        },
    }

    # ── Run ───────────────────────────────────────────────────────────────────
    _success = False
    _error = None

    import haddock as _haddock_pkg
    _haddock_pkg.log.addHandler(_handler)
    try:
        from haddock.libs.libworkflow import WorkflowManager
        with _cwd(_run_dir):
            _wf = WorkflowManager(_workflow_params, start=0)
            _wf.run()
        _success = True
    except Exception:
        _error = traceback.format_exc()
        _log_lines.append(_error)
    finally:
        _haddock_pkg.log.removeHandler(_handler)

    _render_log(done=_success)

    run_result = {"success": _success, "error": _error, "run_dir": _run_dir}
    return (run_result,)


@app.cell
def _(mo, run_result):
    import html as _html
    import json as _json
    import pandas as _pd
    from pathlib import Path as _Path

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _iframe(fragment: str, height: int = 900) -> mo.Html:
        """Embed an HTML fragment via srcdoc so inline <script> tags execute."""
        _doc = (
            "<!DOCTYPE html><html><head><meta charset='utf-8'></head>"
            f"<body style='margin:0;padding:0'>{fragment}</body></html>"
        )
        _escaped = _html.escape(_doc, quote=True)
        return mo.Html(
            f'<iframe srcdoc="{_escaped}" style="width:100%;height:{height}px;border:none;"></iframe>'
        )

    def _mol_viewer(pdb_path: _Path, height: int = 500) -> mo.Html:
        """3Dmol.js structure viewer embedded in an srcdoc iframe."""
        _pdb_json = _json.dumps(pdb_path.read_text())
        _doc = (
            "<!DOCTYPE html><html>"
            "<head><meta charset='utf-8'>"
            "<script src='https://3Dmol.csb.pitt.edu/build/3Dmol-min.js'></script>"
            "</head>"
            f"<body style='margin:0;padding:0;overflow:hidden'>"
            f"<div id='v' style='width:100%;height:{height}px;position:relative;'></div>"
            "<script>"
            "let v=$3Dmol.createViewer(document.getElementById('v'),{backgroundColor:'white'});"
            f"v.addModel({_pdb_json},'pdb');"
            "v.setStyle({},{cartoon:{colorscheme:'chain'}});"
            "v.zoomTo();"
            "v.render();"
            "</script>"
            "</body></html>"
        )
        _escaped = _html.escape(_doc, quote=True)
        return mo.Html(
            f'<iframe srcdoc="{_escaped}" style="width:100%;height:{height}px;border:none;"></iframe>'
        )

    def _resolve_path(model_str: str, run_dir: _Path) -> _Path:
        """Resolve a model path that may be absolute or relative to the run dir."""
        p = _Path(model_str)
        if p.is_absolute():
            return p
        # capri_ss.tsv paths are often relative to the caprieval step dir
        for base in [run_dir, run_dir / "3_caprieval"]:
            candidate = (base / p).resolve()
            if candidate.exists():
                return candidate
        return p

    # ── Results ───────────────────────────────────────────────────────────────

    if not run_result["success"]:
        _out = mo.callout(
            mo.vstack([
                mo.md("**Workflow failed.**"),
                mo.code(run_result["error"], language="text"),
            ]),
            kind="danger",
        )
    else:
        _run_dir = run_result["run_dir"]
        _caprieval_dir = _run_dir / "3_caprieval"
        _contactmap_dir = _run_dir / "4_contactmap"

        _sections = [
            mo.callout(
                mo.md(f"Workflow completed. Results saved in `{_run_dir}`"),
                kind="success",
            )
        ]

        # ── Single-structure statistics ───────────────────────────────────────
        _ss_file = _caprieval_dir / "capri_ss.tsv"
        _df_ss = None

        if _ss_file.exists():
            _sections.append(mo.md("---\n## Single-structure statistics"))
            try:
                _df_ss = _pd.read_csv(_ss_file, sep="\t", comment="#")
                _fc = _df_ss.select_dtypes(include="float").columns
                _df_ss[_fc] = _df_ss[_fc].round(3)
                _sections.append(mo.ui.table(
                    _df_ss,
                    pagination=False,
                    selection=None,
                    show_column_summaries=True,
                ))
            except Exception as _e:
                _sections.append(mo.callout(
                    mo.md(f"Could not read `capri_ss.tsv`: {_e}"), kind="warn",
                ))
        else:
            _sections.append(mo.callout(
                mo.md(f"`capri_ss.tsv` not found in `{_caprieval_dir}`."), kind="warn",
            ))

        # ── Cluster statistics ────────────────────────────────────────────────
        _clt_file = _caprieval_dir / "capri_clt.tsv"

        if _clt_file.exists():
            _sections.append(mo.md("---\n## Cluster statistics"))
            try:
                _df_clt = _pd.read_csv(_clt_file, sep="\t", comment="#")
                _fc = _df_clt.select_dtypes(include="float").columns
                _df_clt[_fc] = _df_clt[_fc].round(3)
                _sections.append(mo.ui.table(
                    _df_clt,
                    pagination=False,
                    selection=None,
                    show_column_summaries=True,
                ))
            except Exception as _e:
                _sections.append(mo.callout(
                    mo.md(f"Could not read `capri_clt.tsv`: {_e}"), kind="warn",
                ))

        # ── Best model per cluster: 3D viewer + chord chart ───────────────────
        if _df_ss is not None and "cluster_id" in _df_ss.columns:
            _clustered = _df_ss[_df_ss["cluster_id"] != -1].copy()
            _unclustered = _df_ss[_df_ss["cluster_id"] == -1].copy()

            if not _clustered.empty:
                _sections.append(mo.md("---\n## Best model per cluster"))
                # Order clusters by mean score (best first)
                _clt_order = (
                    _clustered.groupby("cluster_id")["score"]
                    .mean().sort_values().index.tolist()
                )
                for _clt_id in _clt_order:
                    _clt_rows = _clustered[_clustered["cluster_id"] == _clt_id]
                    _best = _clt_rows.loc[_clt_rows["score"].idxmin()]
                    _model_path = _resolve_path(str(_best["model"]), _run_dir)
                    _stem = _model_path.stem

                    _sections.append(mo.md(
                        f"### Cluster {int(_clt_id)}"
                        f"&nbsp;|&nbsp; {len(_clt_rows)} models"
                        f"&nbsp;|&nbsp; best score: **{_best['score']:.3f}**"
                    ))

                    if _model_path.exists():
                        _sections.append(_mol_viewer(_model_path, height=500))
                    else:
                        _sections.append(mo.callout(
                            mo.md(f"Model file not found: `{_model_path}`"), kind="warn",
                        ))

                    # Chord chart for this model (contactmap uses the model stem)
                    if _contactmap_dir.exists():
                        _cf = _contactmap_dir / f"{_stem}_chordchart.html"
                        _hf = _contactmap_dir / f"{_stem}_heatmap.html"
                        if _cf.exists():
                            _m = __import__("re").search(r'style="height:(\d+)px', _cf.read_text())
                            _ch = int(_m.group(1)) + 30 if _m else 930
                            _sections.append(_iframe(_cf.read_text(), height=_ch))
                        if _hf.exists():
                            _sections.append(_iframe(_hf.read_text(), height=900))

            elif not _unclustered.empty:
                # No clusters formed — show top unclustered models
                _sections.append(mo.md("---\n## Top-ranked models (no clusters formed)"))
                for _, _row in _unclustered.nsmallest(5, "score").iterrows():
                    _model_path = _resolve_path(str(_row["model"]), _run_dir)
                    _stem = _model_path.stem
                    _sections.append(mo.md(
                        f"### `{_model_path.name}` — score: **{_row['score']:.3f}**"
                    ))
                    if _model_path.exists():
                        _sections.append(_mol_viewer(_model_path, height=500))
                    if _contactmap_dir.exists():
                        _cf = _contactmap_dir / f"{_stem}_chordchart.html"
                        if _cf.exists():
                            _sections.append(_iframe(_cf.read_text(), height=930))

        _out = mo.vstack(_sections)

    _out
    return


if __name__ == "__main__":
    app.run()
