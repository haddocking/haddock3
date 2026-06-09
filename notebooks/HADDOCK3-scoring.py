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

    Upload a PDB file of complex containing one or an ensemble of model to score and cluster them with HADDOCK3.
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
        0.40, 1.00, value=0.60, step=0.05,
        label="FCC clustering cutoff", show_value=True,
    )
    min_population_slider = mo.ui.slider(
        1, 20, value=4, step=1,
        label="Min cluster population", show_value=True,
    )
    topX_slider = mo.ui.slider(
        1, 20, value=4, step=1,
        label="Top models for contact map", show_value=True,
    )
    chordchart_toggle = mo.ui.switch(label="Chord chart", value=True)
    heatmap_toggle = mo.ui.switch(label="Heatmap", value=False)
    config_edit_toggle = mo.ui.switch(label="Edit config before running", value=False)
    return (
        chordchart_toggle,
        clust_cutoff_slider,
        config_edit_toggle,
        heatmap_toggle,
        min_population_slider,
        ncores_slider,
        topX_slider,
    )


@app.cell
def _(
    chordchart_toggle,
    clust_cutoff_slider,
    config_edit_toggle,
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
        mo.hstack([mo.md("**Advanced:**"), config_edit_toggle], justify="start"),
    ])
    return


@app.cell
def _(
    chordchart_toggle,
    clust_cutoff_slider,
    config_edit_toggle,
    heatmap_toggle,
    min_population_slider,
    mo,
    topX_slider,
):
    def _make_cfg_str():
        _v = lambda b: "true" if b else "false"
        return "\n".join([
            "# HADDOCK3 scoring workflow — module parameters",
            "# CPU cores, mode, input files and run directory are set by the panel above.",
            "# Add, remove or reorder [module] sections to change the workflow steps.",
            "",
            "[topoaa]",
            "autohis = true",
            "",
            "[emscoring]",
            "",
            "[clustfcc]",
            f"clust_cutoff = {clust_cutoff_slider.value:.2f}",  # default 0.60
            f"min_population = {min_population_slider.value}",
            "",
            "[caprieval]",
            "",
            "[contactmap]",
            f"generate_chordchart = {_v(chordchart_toggle.value)}",
            f"generate_heatmap = {_v(heatmap_toggle.value)}",
            "single_model_analysis = false",
            f"topX = {topX_slider.value}",
            "",
        ])

    _cfg_str = _make_cfg_str()

    if config_edit_toggle.value:
        config_textarea = mo.ui.text_area(
            value=_cfg_str,
            rows=22,
            full_width=True,
            label="Workflow configuration (TOML)",
        )
        _out = mo.vstack([
            mo.callout(
                mo.md(
                    "**Config editor active.** Edit parameters below — these override the panel "
                    "settings for module-specific options. Changing a slider above will regenerate "
                    "this editor and discard manual edits."
                ),
                kind="warn",
            ),
            config_textarea,
        ])
    else:
        config_textarea = mo.ui.text_area(value="")
        _out = mo.accordion(
            {"View workflow configuration": mo.md(f"```toml\n{_cfg_str}\n```")}
        )

    _out
    return (config_textarea,)


@app.cell
def _(mo):
    run_btn = mo.ui.run_button(label="▶  Run HADDOCK3 Scoring Workflow", kind="success")
    run_btn
    return (run_btn,)


@app.cell
def _(mo):
    import sys as _sys
    import threading as _t
    if "_h3nb_stop" not in _sys.modules:
        _m = type(_sys)("_h3nb_stop")
        _m.event = _t.Event()
        _sys.modules["_h3nb_stop"] = _m
    stop_btn = mo.ui.run_button(label="⏹  Stop workflow", kind="danger")
    if stop_btn.value:
        _sys.modules["_h3nb_stop"].event.set()
    stop_btn
    return


@app.cell
async def _(
    chordchart_toggle,
    clust_cutoff_slider,
    config_edit_toggle,
    config_textarea,
    heatmap_toggle,
    min_population_slider,
    mo,
    ncores_slider,
    pdb_upload,
    run_btn,
    topX_slider,
    work_dir_input,
):
    import asyncio as _asyncio
    import logging
    import os
    import sys as _sys
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

    def _render_log(done=False, stopped=False):
        if stopped:
            _label = "🛑 Stopped"
        elif done:
            _label = "✅ Completed"
        else:
            _label = "⏳ Running…"
        _body = "\n".join(_log_lines) if _log_lines else "(waiting for output…)"
        _panel = mo.Html(
            '<div style="height:400px;overflow-y:auto;background:#f0f0f0;'
            'color:#1a1a1a;font-family:monospace;padding:10px 14px;'
            'border-radius:6px;font-size:12px;white-space:pre-wrap;line-height:1.5;">'
            + _body.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            + "</div>"
        )
        mo.output.replace(mo.accordion({f"HADDOCK3 Log — {_label}": _panel}))

    # Log handler only appends; rendering happens after each step in the async loop.
    class _LogHandler(logging.Handler):
        def emit(self, record):
            _log_lines.append(self.format(record))

    _handler = _LogHandler()
    _handler.setFormatter(logging.Formatter("[%(asctime)s %(module)s %(levelname)s] %(message)s"))
    _render_log()

    # ── Workflow parameters ───────────────────────────────────────────────────
    _exec = {"mode": "local", "ncores": _ncores, "clean": False, "offline": False}

    if config_edit_toggle.value and config_textarea.value.strip():
        from haddock.gear.config import loads as _cfg_loads
        try:
            _parsed = _cfg_loads(config_textarea.value)["final_cfg"]
        except Exception as _cfg_err:
            mo.stop(True, mo.callout(
                mo.md(f"**Config parse error:** {_cfg_err}"),
                kind="danger",
            ))
        _workflow_params = {}
        for _step_key, _step_params in _parsed.items():
            _p = {**_exec, **_step_params}
            if _step_key.startswith("topoaa"):
                _p.setdefault("molecules", _pdb_paths)
                _p.setdefault("autohis", True)
            _workflow_params[_step_key] = _p
    else:
        _workflow_params = {
            "topoaa.1": {
                "molecules": _pdb_paths,
                "autohis": True,
                **_exec,
            },
            "emscoring.1": {**_exec},
            "clustfcc.1": {
                "clust_cutoff": clust_cutoff_slider.value,
                "min_population": min_population_slider.value,
                **_exec,
            },
            "caprieval.1": {**_exec},
            "contactmap.1": {
                "generate_chordchart": chordchart_toggle.value,
                "generate_heatmap": heatmap_toggle.value,
                "single_model_analysis": False,
                "topX": topX_slider.value,
                **_exec,
            },
        }

    # ── Run ───────────────────────────────────────────────────────────────────
    # Reset the stop flag from any previous run.
    _stop = _sys.modules.get("_h3nb_stop")
    if _stop:
        _stop.event.clear()

    _success = False
    _error = None
    _stopped = False

    import haddock as _haddock_pkg
    _haddock_pkg.log.addHandler(_handler)
    try:
        from haddock.libs.libworkflow import WorkflowManager
        with _cwd(_run_dir):
            _wf = WorkflowManager(_workflow_params, start=0)
            # Run steps one at a time so the event loop can process the stop button
            # between steps (asyncio.to_thread yields to the event loop while each
            # step executes in a thread-pool worker).
            for _i, _step in enumerate(_wf.recipe.steps):
                if _stop and _stop.event.is_set():
                    _log_lines.append("[notebook] Stop requested — workflow halted after current step.")
                    _stopped = True
                    break
                await _asyncio.to_thread(_step.execute)
                _render_log()  # update log panel after each completed step
                # Replicate WorkflowManager.run() output-param forwarding
                _op = getattr(_step.module, "_output_params", {})
                if _op:
                    for _fs in _wf.recipe.steps[_i + 1:]:
                        for _k, _v in _op.items():
                            if _fs.config.get(_k) in (None, ""):
                                _fs.config[_k] = _v
            else:
                _success = True
    except _asyncio.CancelledError:
        _log_lines.append("[notebook] Workflow execution was cancelled.")
        _render_log(done=False)
        raise
    except Exception:
        _error = traceback.format_exc()
        _log_lines.append(_error)
    finally:
        _haddock_pkg.log.removeHandler(_handler)

    _render_log(done=_success, stopped=_stopped)

    run_result = {"success": _success, "error": _error, "run_dir": _run_dir, "stopped": _stopped}
    return (run_result,)


@app.cell
def _(mo):
    clt_show_std = mo.ui.switch(label="Show standard deviations", value=False)
    return (clt_show_std,)


@app.cell
def _(clt_show_std, mo, run_result):
    import pandas as _pd
    from pathlib import Path as _Path

    # Default exports — overwritten when a successful run with clusters exists.
    clt_table = mo.ui.table(_pd.DataFrame(), selection="single", pagination=False)
    df_ss = None

    if run_result.get("stopped"):
        _out = mo.callout(
            mo.md(f"Workflow stopped by user. Partial results may be available in `{run_result['run_dir']}`."),
            kind="warn",
        )
    elif not run_result["success"]:
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

        _sections = [
            mo.callout(
                mo.md(f"Workflow completed. Results saved in `{_run_dir}`"),
                kind="success",
            )
        ]

        # ── Cluster statistics (row selection drives the visualisation below) ──
        _clt_file = _caprieval_dir / "capri_clt.tsv"

        if _clt_file.exists():
            _sections.append(mo.md("---\n## Cluster statistics"))
            _sections.append(mo.md(
                "*Click a row to display the best model and contact map for that cluster.*"
            ))
            try:
                _df_clt = _pd.read_csv(_clt_file, sep="\t", comment="#")
                _fc = _df_clt.select_dtypes(include="float").columns
                _df_clt[_fc] = _df_clt[_fc].round(3)
                _std_cols = [c for c in _df_clt.columns if c.endswith("_std")]
                clt_table = mo.ui.table(
                    _df_clt,
                    selection="single",
                    pagination=False,
                    show_column_summaries=True,
                    hidden_columns=[] if clt_show_std.value else _std_cols,
                )
                _sections.append(mo.hstack([clt_show_std], justify="end"))
                _sections.append(clt_table)
            except Exception as _e:
                _sections.append(mo.callout(
                    mo.md(f"Could not read `capri_clt.tsv`: {_e}"), kind="warn",
                ))

        # ── Single-structure statistics (collapsed by default) ────────────────
        _ss_file = _caprieval_dir / "capri_ss.tsv"

        if _ss_file.exists():
            try:
                df_ss = _pd.read_csv(_ss_file, sep="\t", comment="#")
                _fc = df_ss.select_dtypes(include="float").columns
                df_ss[_fc] = df_ss[_fc].round(3)
                _sections.append(mo.accordion({
                    f"Single-structure statistics ({len(df_ss)} models)": mo.ui.table(
                        df_ss,
                        pagination=False,
                        selection=None,
                        show_column_summaries=True,
                    )
                }))
            except Exception as _e:
                _sections.append(mo.callout(
                    mo.md(f"Could not read `capri_ss.tsv`: {_e}"), kind="warn",
                ))
        else:
            _sections.append(mo.callout(
                mo.md(f"`capri_ss.tsv` not found in `{_caprieval_dir}`."), kind="warn",
            ))

        _out = mo.vstack(_sections)

    _out
    return (clt_table, df_ss)


@app.cell
def _(clt_table, df_ss, mo, run_result):
    import html as _html
    import json as _json
    import re as _re
    from pathlib import Path as _Path2

    def _iframe(fragment: str, height: int = 900) -> mo.Html:
        _doc = (
            "<!DOCTYPE html><html><head><meta charset='utf-8'></head>"
            f"<body style='margin:0;padding:0'>{fragment}</body></html>"
        )
        _esc = _html.escape(_doc, quote=True)
        return mo.Html(f'<iframe srcdoc="{_esc}" style="width:100%;height:{height}px;border:none;"></iframe>')

    def _mol_viewer(pdb_path: _Path2, height: int = 500) -> mo.Html:
        _pdb_json = _json.dumps(pdb_path.read_text())
        _doc = (
            "<!DOCTYPE html><html><head><meta charset='utf-8'>"
            "<script src='https://3Dmol.csb.pitt.edu/build/3Dmol-min.js'></script></head>"
            f"<body style='margin:0;padding:0;overflow:hidden'>"
            f"<div id='v' style='width:100%;height:{height}px;position:relative;'></div>"
            "<script>"
            "let v=$3Dmol.createViewer(document.getElementById('v'),{backgroundColor:'white'});"
            f"v.addModel({_pdb_json},'pdb');"
            "v.setStyle({},{cartoon:{colorscheme:'chain'}});"
            "v.addStyle({elem:'C'},{stick:{colorscheme:'chain',radius:0.12}});"
            "v.addStyle({elem:'N'},{stick:{color:'#3050F8',radius:0.12}});"
            "v.addStyle({elem:'O'},{stick:{color:'#FF0D0D',radius:0.12}});"
            "v.addStyle({elem:'S'},{stick:{color:'#FFFF30',radius:0.12}});"
            "v.addStyle({elem:'P'},{stick:{color:'#FF8000',radius:0.12}});"
            "v.zoomTo();v.render();"
            "</script></body></html>"
        )
        _esc = _html.escape(_doc, quote=True)
        return mo.Html(f'<iframe srcdoc="{_esc}" style="width:100%;height:{height}px;border:none;"></iframe>')

    def _resolve_path(model_str: str, run_dir: _Path2) -> _Path2:
        p = _Path2(model_str)
        if p.is_absolute():
            return p
        for base in [run_dir, run_dir / "3_caprieval"]:
            candidate = (base / p).resolve()
            if candidate.exists():
                return candidate
        return p

    def _add_contactmap(base_name: str, contactmap_dir: _Path2, sections: list) -> None:
        _cf = contactmap_dir / f"{base_name}_chordchart.html"
        _hf = contactmap_dir / f"{base_name}_heatmap.html"
        if _cf.exists():
            _txt = _cf.read_text()
            _m = _re.search(r'style="height:(\d+)px', _txt)
            _ch = int(_m.group(1)) + 30 if _m else 930
            sections.append(_iframe(_txt, height=_ch))
        if _hf.exists():
            sections.append(_iframe(_hf.read_text(), height=900))

    # ── Visualisation — driven by cluster table selection ─────────────────────

    if not run_result["success"] or run_result.get("stopped"):
        _viz_out = mo.md("")
    elif clt_table.value.empty:
        _viz_out = mo.callout(
            mo.md("Select a row in the cluster table above to display the best model and contact map."),
            kind="info",
        )
    else:
        try:
            _sel = clt_table.value.iloc[0]
            # cluster_rank may be int or float depending on pandas dtype — normalise to int
            _clt_rank = int(float(_sel["cluster_rank"]))
            _clt_id   = int(float(_sel["cluster_id"]))
            _score    = float(_sel["score"])

            _run_dir       = run_result["run_dir"]
            _contactmap_dir = _run_dir / "4_contactmap"

            _vsections = [mo.md(
                f"---\n## Cluster {_clt_id}"
                f"&nbsp;|&nbsp; rank {_clt_rank}"
                f"&nbsp;|&nbsp; mean score: **{_score:.3f}**"
            )]

            # ── 3D viewer ──────────────────────────────────────────────────
            if df_ss is None:
                _vsections.append(mo.callout(
                    mo.md("Single-structure data (`capri_ss.tsv`) not available — cannot show 3D model."),
                    kind="warn",
                ))
            elif "cluster_ranking" not in df_ss.columns:
                _vsections.append(mo.callout(
                    mo.md(f"`cluster_ranking` column not found in capri_ss.tsv. Columns: {list(df_ss.columns)}"),
                    kind="warn",
                ))
            else:
                # cluster_ranking may also be float — compare as float
                _clt_models = df_ss[df_ss["cluster_ranking"].apply(
                    lambda x: int(float(x)) == _clt_rank if x == x else False  # NaN-safe
                )]
                if _clt_models.empty:
                    _vsections.append(mo.callout(
                        mo.md(
                            f"No models found for cluster rank **{_clt_rank}** in capri_ss.tsv. "
                            f"Available ranks: {sorted(df_ss['cluster_ranking'].dropna().unique().tolist())}"
                        ),
                        kind="warn",
                    ))
                else:
                    _best = _clt_models.loc[_clt_models["score"].idxmin()]
                    _model_path = _resolve_path(str(_best["model"]), _run_dir)
                    _vsections.append(mo.md(
                        f"Best model: `{_model_path.name}` &nbsp;|&nbsp; score: **{_best['score']:.3f}**"
                    ))
                    if _model_path.exists():
                        _vsections.append(_mol_viewer(_model_path, height=500))
                    else:
                        _vsections.append(mo.callout(
                            mo.md(f"Model file not found: `{_model_path}`"), kind="warn",
                        ))

            # ── Chord chart / heatmap ──────────────────────────────────────
            if not _contactmap_dir.exists():
                _vsections.append(mo.callout(
                    mo.md(f"Contact map directory not found: `{_contactmap_dir}`"), kind="warn",
                ))
            else:
                _cf = _contactmap_dir / f"cluster{_clt_rank}_chordchart.html"
                _hf = _contactmap_dir / f"cluster{_clt_rank}_heatmap.html"
                if not _cf.exists() and not _hf.exists():
                    _vsections.append(mo.callout(
                        mo.md(
                            f"No contact map files found for cluster rank **{_clt_rank}**. "
                            f"Expected `{_cf.name}` in `{_contactmap_dir}`."
                        ),
                        kind="warn",
                    ))
                else:
                    _add_contactmap(f"cluster{_clt_rank}", _contactmap_dir, _vsections)

            _viz_out = mo.vstack(_vsections)

        except Exception as _exc:
            import traceback as _tb
            _viz_out = mo.callout(
                mo.vstack([mo.md("**Error rendering cluster visualisation:**"),
                           mo.code(_tb.format_exc(), language="text")]),
                kind="danger",
            )

    _viz_out
    return


if __name__ == "__main__":
    app.run()
