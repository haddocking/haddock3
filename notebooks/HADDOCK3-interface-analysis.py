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
    # HADDOCK3 Interface Analysis

    Upload a multi-chain PDB complex to run an energy-minimization scoring
    workflow and analyse intermolecular contacts. Results are displayed as an
    interactive chord chart and heatmap. An optional alanine scanning step
    estimates the energetic contribution of each interface residue by
    systematically mutating it to alanine and reporting the change in HADDOCK
    score.

    **Workflow steps:**

    1. **`topoaa`** — builds CNS all-atom topology and parameters
    2. **`emscoring`** — energy minimization scoring via CNS
    3. **`contactmap`** — computes residue–residue contacts and generates figures
    4. **`alascan`** *(optional)* — alanine scanning of interface residues

    > **Note:** This notebook supports standard protein, DNA/RNA, and ion
    > inputs only. Small-molecule ligands requiring custom topology and
    > parameter files are **not** supported.

    > **Dislaimer:** This notebook was generated with help of AI.
    """)
    return


@app.cell
def _(mo):
    pdb_upload = mo.ui.file(
        filetypes=[".pdb"],
        label="Drop or select a PDB complex file",
        kind="area",
    )
    return (pdb_upload,)


@app.cell
def _(mo):
    import os as _os

    _std_aa = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]

    _ncores_max = _os.cpu_count() or 1
    ncores_slider = mo.ui.slider(
        1,
        _ncores_max,
        value=_ncores_max,
        step=1,
        label="CPU cores",
        show_value=True,
    )
    chordchart_toggle = mo.ui.switch(label="Generate chord chart", value=True)
    heatmap_toggle = mo.ui.switch(label="Generate heatmap", value=False)
    single_model_analysis_toggle = mo.ui.switch(label="Per-model analysis", value=False)
    alascan_toggle = mo.ui.switch(label="Run alanine scan", value=False)
    alascan_scan_residue = mo.ui.dropdown(
        options=_std_aa,
        value="ALA",
        label="Scan residue",
    )
    alascan_plot_toggle = mo.ui.switch(label="Generate plot", value=True)
    alascan_show_std = mo.ui.switch(label="Show std deviations", value=False)
    config_edit_toggle = mo.ui.switch(label="Edit config before running", value=False)
    return (
        alascan_plot_toggle,
        alascan_scan_residue,
        alascan_show_std,
        alascan_toggle,
        chordchart_toggle,
        config_edit_toggle,
        heatmap_toggle,
        ncores_slider,
        single_model_analysis_toggle,
    )


@app.cell
def _(mo):
    get_topX, set_topX = mo.state(1)
    return get_topX, set_topX


@app.cell
def _(get_topX, mo):
    topX_slider = mo.ui.slider(
        1,
        20,
        value=get_topX(),
        step=1,
        label="Models to analyse (topX)",
        show_value=True,
    )
    return (topX_slider,)


@app.cell
def _(set_topX, single_model_analysis_toggle, topX_slider):
    if single_model_analysis_toggle.value and topX_slider.value < 10:
        set_topX(10)
    return


@app.cell
def _(
    alascan_plot_toggle,
    alascan_scan_residue,
    alascan_show_std,
    alascan_toggle,
    chordchart_toggle,
    config_edit_toggle,
    heatmap_toggle,
    mo,
    ncores_slider,
    pdb_upload,
    single_model_analysis_toggle,
    topX_slider,
    work_dir_input,
):
    _alascan_sub = (
        mo.vstack(
            [alascan_scan_residue, alascan_plot_toggle, alascan_show_std], align="start"
        )
        if alascan_toggle.value
        else mo.md("")
    )
    mo.vstack(
        [
            mo.md("### Configuration"),
            mo.hstack(
                [
                    mo.vstack(
                        [
                            mo.md("**Input file**"),
                            pdb_upload,
                        ],
                        align="start",
                    ),
                    mo.vstack(
                        [
                            mo.md("**Working directory**"),
                            work_dir_input,
                        ],
                        align="start",
                    ),
                    mo.vstack(
                        [
                            mo.md("**Execution**"),
                            ncores_slider,
                            topX_slider,
                        ],
                        align="start",
                    ),
                    mo.vstack(
                        [
                            mo.md("**contactmap outputs**"),
                            chordchart_toggle,
                            heatmap_toggle,
                            single_model_analysis_toggle,
                        ],
                        align="start",
                    ),
                    mo.vstack(
                        [
                            mo.md("**alanine scan**"),
                            alascan_toggle,
                            _alascan_sub,
                        ],
                        align="start",
                    ),
                ],
                gap=2,
                justify="start",
            ),
            mo.hstack([mo.md("**Advanced:**"), config_edit_toggle], justify="start"),
        ]
    )
    return


@app.cell
def _(
    alascan_plot_toggle,
    alascan_scan_residue,
    alascan_toggle,
    chordchart_toggle,
    config_edit_toggle,
    heatmap_toggle,
    mo,
    single_model_analysis_toggle,
    topX_slider,
):
    def _make_cfg_str():
        _v = lambda b: "true" if b else "false"
        _topX = (
            max(topX_slider.value, 10) if single_model_analysis_toggle.value else topX_slider.value
        )
        _lines = [
            "# HADDOCK3 interface analysis workflow — module parameters",
            "# CPU cores, mode, input file and run directory are set by the panel above.",
            "# Add, remove or reorder [module] sections to change the workflow steps.",
            "",
            "[topoaa]",
            "autohis = true",
            "",
            "[emscoring]",
            "per_interface_scoring = true",
            "",
            "[caprieval]",
            "[contactmap]",
            f"generate_chordchart = {_v(chordchart_toggle.value)}",
            f"generate_heatmap = {_v(heatmap_toggle.value)}",
            f"single_model_analysis = {_v(single_model_analysis_toggle.value)}",
            f"topX = {_topX}",
            "",
        ]
        if alascan_toggle.value:
            _lines += [
                "[alascan]",
                f'scan_residue = "{alascan_scan_residue.value}"',
                f"plot = {_v(alascan_plot_toggle.value)}",
                "",
            ]
        return "\n".join(_lines)

    _cfg_str = _make_cfg_str()

    if config_edit_toggle.value:
        config_textarea = mo.ui.text_area(
            value=_cfg_str,
            rows=20,
            full_width=True,
            label="Workflow configuration (TOML)",
        )
        _out = mo.vstack(
            [
                mo.callout(
                    mo.md(
                        "**Config editor active.** Edit parameters below — these override the panel "
                        "settings for module-specific options. Changing a slider or toggle above will "
                        "regenerate this editor and discard manual edits."
                    ),
                    kind="warn",
                ),
                config_textarea,
            ]
        )
    else:
        config_textarea = mo.ui.text_area(value="")
        _out = mo.accordion(
            {"View workflow configuration": mo.md(f"```toml\n{_cfg_str}\n```")}
        )

    _out
    return (config_textarea,)


@app.cell
def _(mo):
    run_btn = mo.ui.run_button(label="▶  Run HADDOCK3 Workflow", kind="success")
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
    stop_btn
    return (stop_btn,)


@app.cell
def _(stop_btn):
    import sys as _sys2

    if stop_btn.value:
        _sys2.modules["_h3nb_stop"].event.set()
    return


@app.cell
async def _(
    alascan_plot_toggle,
    alascan_scan_residue,
    alascan_toggle,
    chordchart_toggle,
    config_edit_toggle,
    config_textarea,
    heatmap_toggle,
    mo,
    ncores_slider,
    pdb_upload,
    run_btn,
    single_model_analysis_toggle,
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
            mo.md(
                "**No PDB file loaded.** Upload a PDB complex above, then click Run."
            ),
            kind="warn",
        ),
    )

    # ── Run directory ─────────────────────────────────────────────────────────
    _pdb_name = pdb_upload.name()
    _pdb_stem = Path(_pdb_name).stem
    _timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    _base_dir = Path(work_dir_input.value)
    _base_dir.mkdir(parents=True, exist_ok=True)
    _run_dir = _base_dir / f"{_pdb_stem}_{_timestamp}"
    _run_dir.mkdir(exist_ok=True)
    (_run_dir / _pdb_name).write_bytes(pdb_upload.contents())

    _ncores = ncores_slider.value

    @contextmanager
    def _cwd(path):
        _prev = Path.cwd()
        os.chdir(path)
        try:
            yield
        finally:
            os.chdir(_prev)

    # ── Live log ──────────────────────────────────────────────────────────────
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
            "color:#1a1a1a;font-family:monospace;padding:10px 14px;"
            "border-radius:6px;font-size:12px;white-space:pre-wrap;"
            'line-height:1.5;">'
            + _body.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            + "</div>"
        )
        mo.output.replace(mo.accordion({f"HADDOCK3 Log — {_label}": _panel}))

    # Log handler only appends; rendering happens after each step in the async loop.
    class _LogHandler(logging.Handler):
        def emit(self, record):
            _log_lines.append(self.format(record))

    _handler = _LogHandler()
    _handler.setFormatter(
        logging.Formatter("[%(asctime)s %(module)s %(levelname)s] %(message)s")
    )
    _render_log()

    # ── Workflow parameters ───────────────────────────────────────────────────
    _exec = {"mode": "local", "ncores": _ncores, "clean": False, "offline": False}

    if config_edit_toggle.value and config_textarea.value.strip():
        from haddock.gear.config import loads as _cfg_loads

        try:
            _parsed = _cfg_loads(config_textarea.value)["final_cfg"]
        except Exception as _cfg_err:
            mo.stop(
                True,
                mo.callout(
                    mo.md(f"**Config parse error:** {_cfg_err}"),
                    kind="danger",
                ),
            )
        _workflow_params = {}
        for _step_key, _step_params in _parsed.items():
            _p = {**_exec, **_step_params}
            if _step_key.startswith("topoaa"):
                _p.setdefault("molecules", [Path(_pdb_name)])
                _p.setdefault("autohis", True)
            _workflow_params[_step_key] = _p
    else:
        _workflow_params = {
            "topoaa.1": {
                "molecules": [Path(_pdb_name)],
                "autohis": True,
                **_exec,
            },
            "emscoring.1": {"per_interface_scoring": True, **_exec},
            "contactmap.1": {
                "generate_chordchart": chordchart_toggle.value,
                "generate_heatmap": heatmap_toggle.value,
                "single_model_analysis": single_model_analysis_toggle.value,
                "topX": (
                    max(topX_slider.value, 10)
                    if single_model_analysis_toggle.value
                    else topX_slider.value
                ),
                **_exec,
            },
        }
        if alascan_toggle.value:
            _workflow_params["alascan.1"] = {
                "scan_residue": alascan_scan_residue.value,
                "plot": alascan_plot_toggle.value,
                **_exec,
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
                    _log_lines.append(
                        "[notebook] Stop requested — workflow halted after current step."
                    )
                    _stopped = True
                    break
                await _asyncio.to_thread(_step.execute)
                _render_log()  # update log panel after each completed step
                # Replicate WorkflowManager.run() output-param forwarding
                _op = getattr(_step.module, "_output_params", {})
                if _op:
                    for _fs in _wf.recipe.steps[_i + 1 :]:
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

    run_result = {
        "success": _success,
        "error": _error,
        "run_dir": _run_dir,
        "stopped": _stopped,
    }
    return (run_result,)


@app.cell
def _(mo, run_result):
    import pandas as _pd
    from pathlib import Path as _Path

    score_table = None
    model_map: dict = {}

    if not run_result["success"] or run_result.get("stopped"):
        _score_out = mo.md("")
    else:
        try:
            _run_dir = run_result["run_dir"]
            _em_dir = next(iter(sorted(_run_dir.glob("*_emscoring"))), None)
            _main_tsv = _em_dir / "emscoring.tsv" if _em_dir else None

            if not _main_tsv or not _main_tsv.exists():
                _score_out = mo.callout(
                    mo.md(f"`emscoring.tsv` not found in `{_em_dir}`."),
                    kind="warn",
                )
            else:
                # ── Main HADDOCK score ─────────────────────────────────────────
                _df = _pd.read_csv(_main_tsv, sep="\t")
                _df["model"] = _df["original_name"].apply(lambda x: _Path(x).name)
                _df = _df.rename(columns={"score": "HADDOCK score"})
                _df = _df.drop(columns=["original_name", "md5"], errors="ignore")

                # ── Energy terms from PDB REMARK headers ──────────────────────
                from haddock.gear.haddockmodel import HaddockModel as _HaddockModel

                _extra = []
                for _struct in _df["structure"]:
                    _pdb = _em_dir / _struct
                    _row: dict = {}
                    if _pdb.exists():
                        try:
                            _e = _HaddockModel(str(_pdb)).energies
                            _row = {
                                "Evdw": _e.get("vdw"),
                                "Eelec": _e.get("elec"),
                                "Edesolv": _e.get("desolv"),
                                "BSA": _e.get("bsa"),
                            }
                        except Exception:
                            pass
                    _extra.append(_row)
                _df = _pd.concat(
                    [_df.reset_index(drop=True), _pd.DataFrame(_extra)], axis=1
                )

                # ── Per-interface scores from emscoring_{interface}.tsv ────────
                for _iface_tsv in sorted(_em_dir.glob("emscoring_*.tsv")):
                    _iface = _iface_tsv.stem[len("emscoring_") :]  # e.g., "A_B"
                    _di = _pd.read_csv(_iface_tsv, sep="\t")[["structure", "score"]]
                    _di = _di.rename(columns={"score": f"score_{_iface}"})
                    _df = _df.merge(_di, on="structure", how="left")

                # ── Final column order ─────────────────────────────────────────
                # Keep "structure" (the emscoring PDB stem) so the visualisation
                # cell can match it against contactmap filenames; hide it in the UI.
                _base = ["structure", "model", "HADDOCK score", "Evdw", "Eelec", "Edesolv", "BSA"]
                _iface_cols = [c for c in _df.columns if c.startswith("score_")]
                _df = _df[[c for c in _base if c in _df.columns] + _iface_cols]
                _fc = _df.select_dtypes(include="float").columns
                _df[_fc] = _df[_fc].round(3)

                # Map structure stem → original model name for chart labelling.
                model_map = dict(
                    zip(
                        _df["structure"].apply(lambda s: s.rsplit(".", 1)[0]),
                        _df["model"],
                    )
                ) if "structure" in _df.columns and "model" in _df.columns else {}

                # Replace "model" column with HTML download links via base64 data URI.
                import base64 as _b64
                _model_links = []
                for _, _row in _df.iterrows():
                    _pdb_path = (_em_dir / str(_row["structure"])) if _em_dir else None
                    _display = str(_row["model"])
                    if _pdb_path and _pdb_path.exists():
                        _b64_data = _b64.b64encode(_pdb_path.read_bytes()).decode()
                        _dl_name = _display.rsplit(".", 1)[0] + "_EM.pdb"
                        _model_links.append(
                            mo.Html(
                                f'<a href="data:chemical/x-pdb;base64,{_b64_data}"'
                                f' download="{_dl_name}">{_display}</a>'
                            )
                        )
                    else:
                        _model_links.append(_display)
                _df = _df.copy()
                _df["model"] = _model_links

                score_table = mo.ui.table(
                    _df,
                    pagination=False,
                    selection="single",
                    show_column_summaries=False,
                    hidden_columns=["structure"],
                )
                _hint = (
                    mo.callout(
                        mo.md(
                            "Click a row to select a model and filter the "
                            "contact-map visualisations below.\n\n"
                            "Click on a model name to download the corresponding energy-minimized PDB file."
                        ),
                        kind="info",
                    )
                    if len(_df) > 1
                    else mo.md("")
                )
                _score_out = mo.vstack(
                    [
                        mo.md("---\n## HADDOCK Scoring"),
                        _hint,
                        score_table,
                    ]
                )

        except Exception:
            import traceback as _tb

            _score_out = mo.callout(
                mo.vstack(
                    [
                        mo.md("**Error loading scoring data:**"),
                        mo.code(_tb.format_exc(), language="text"),
                    ]
                ),
                kind="danger",
            )

    _score_out
    return (model_map, score_table,)


@app.cell
def _(mo, model_map, run_result, score_table):
    import html as _html_3d
    import json as _json_3d
    from pathlib import Path as _Path_3d

    def _mol_viewer(pdb_path: _Path_3d) -> mo.Html:
        _pdb_json = _json_3d.dumps(pdb_path.read_text())
        _doc = (
            "<!DOCTYPE html><html style='height:100%'><head><meta charset='utf-8'>"
            "<script src='https://3Dmol.csb.pitt.edu/build/3Dmol-min.js'></script></head>"
            "<body style='margin:0;padding:0;overflow:hidden;height:100%'>"
            "<div id='v' style='width:100%;height:100%;position:relative;'></div>"
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
        _esc = _html_3d.escape(_doc, quote=True)
        return mo.Html(
            '<iframe srcdoc="' + _esc + '" style="width:100%;aspect-ratio:3/2;border:none;"></iframe>'
        )

    def _resolve_em_path(structure: str, run_dir: _Path_3d) -> _Path_3d:
        _em_dir = next(iter(sorted(run_dir.glob("*_emscoring"))), None)
        if _em_dir:
            _candidate = _em_dir / structure
            if _candidate.exists():
                return _candidate
        # Fallback: search the whole run directory tree
        for _p in run_dir.rglob(structure):
            return _p
        return _Path_3d(structure)

    if not run_result["success"] or run_result.get("stopped"):
        _3d_out = mo.md("")
    elif score_table is None or (
        score_table.value.empty and len(model_map) != 1
    ):
        _3d_out = mo.callout(
            mo.md("Select a row in the scoring table above to display the model in 3D."),
            kind="info",
        )
    else:
        try:
            if not score_table.value.empty:
                _sel = score_table.value.iloc[0]
                _structure = str(_sel.get("structure", ""))
                _score = float(_sel.get("HADDOCK score", 0))
                # "model" column now contains mo.Html; use model_map for the display name.
                _model_name = model_map.get(_structure.rsplit(".", 1)[0], _structure)
            else:
                # Single model: auto-display without requiring a row click
                _struct_stem = next(iter(model_map))
                _structure = _struct_stem + ".pdb"
                _model_name = model_map[_struct_stem]
                _score = float("nan")
            _run_dir = run_result["run_dir"]
            _model_path = _resolve_em_path(_structure, _run_dir)

            import math as _math_3d
            _score_str = f"&nbsp;|&nbsp; HADDOCK score: **{_score:.3f}**" if not _math_3d.isnan(_score) else ""
            _vsections = [
                mo.md(f"---\n## 3D Visualisation — `{_model_name}`{_score_str}")
            ]
            if _model_path.exists():
                _vsections.append(_mol_viewer(_model_path))
            else:
                _vsections.append(
                    mo.callout(
                        mo.md(f"Model file not found: `{_model_path}`"),
                        kind="warn",
                    )
                )
            _3d_out = mo.vstack(_vsections)
        except Exception:
            import traceback as _tb_3d

            _3d_out = mo.callout(
                mo.vstack(
                    [
                        mo.md("**Error rendering 3D visualisation:**"),
                        mo.code(_tb_3d.format_exc(), language="text"),
                    ]
                ),
                kind="danger",
            )

    _3d_out
    return


@app.cell
def _(alascan_show_std, mo, model_map, run_result, score_table):
    import html as _html

    def _iframe(fragment: str, height: int = 900) -> mo.Html:
        """Wrap a plotly HTML fragment in a full document served via srcdoc.

        mo.Html() injects content via innerHTML, which browsers block from
        executing <script> tags. An <iframe srcdoc> runs in its own document
        context so scripts (including the Plotly CDN loader) execute normally.
        """
        _doc = (
            "<!DOCTYPE html><html><head><meta charset='utf-8'></head>"
            f"<body style='margin:0;padding:0'>{fragment}</body></html>"
        )
        _escaped = _html.escape(_doc, quote=True)
        return mo.Html(
            f'<iframe srcdoc="{_escaped}"'
            f' style="width:100%;height:{height}px;border:none;"></iframe>'
        )

    if run_result.get("stopped"):
        _out = mo.callout(
            mo.md(
                f"Workflow stopped by user. Partial results may be available in `{run_result['run_dir']}`."
            ),
            kind="warn",
        )
    elif not run_result["success"]:
        _out = mo.callout(
            mo.vstack(
                [
                    mo.md("**Workflow failed.**"),
                    mo.code(run_result["error"], language="text"),
                ]
            ),
            kind="danger",
        )
    else:
        _run_dir = run_result["run_dir"]

        # ── Resolve contactmap directory (step number varies with workflow) ───
        _contactmap_dir = next(iter(sorted(_run_dir.glob("*_contactmap"))), None)
        _chord_files = (
            sorted(_contactmap_dir.glob("*_chordchart.html")) if _contactmap_dir else []
        )
        _heatmap_files = (
            sorted(_contactmap_dir.glob("*_heatmap.html")) if _contactmap_dir else []
        )

        # ── Filter by selected model ──────────────────────────────────────────
        _has_selection = score_table is not None and not score_table.value.empty
        _multi_charts = len(_chord_files) > 1 or len(_heatmap_files) > 1

        if _has_selection:
            # "structure" is the emscoring PDB name (e.g. "emscoring_2.pdb").
            # Chart files embed that stem (e.g. "Unclustered_emscoring_2_chordchart").
            # Use "_{stem}_" to avoid "emscoring_1" matching "emscoring_10".
            _sel_structure = score_table.value.iloc[0].get("structure", "")
            _sel_stem = _sel_structure.rsplit(".", 1)[0] if _sel_structure else None
            if _sel_stem and _multi_charts:
                _needle = f"_{_sel_stem}_"
                _chord_files = [f for f in _chord_files if _needle in f.stem]
                _heatmap_files = [f for f in _heatmap_files if _needle in f.stem]
        elif _multi_charts:
            # Multiple per-model charts but no row selected yet — defer display.
            _chord_files = []
            _heatmap_files = []

        _sections = [
            mo.callout(
                mo.md(f"Workflow completed. Results saved in `{_run_dir}`"),
                kind="success",
            )
        ]

        def _chart_label(f) -> str:
            _key = next((k for k in model_map if f"_{k}_" in f.stem), None)
            if _key:
                return f"{_key} — {model_map[_key]}"
            return f.stem

        if _chord_files:
            _sections.append(mo.md("---\n## Chord Chart"))
            for _f in _chord_files:
                _sections.append(mo.md(f"**{_chart_label(_f)}**"))
                _sections.append(_iframe(_f.read_text(), height=900))

        if _heatmap_files:
            _sections.append(mo.md("---\n## Heatmap"))
            for _f in _heatmap_files:
                _sections.append(mo.md(f"**{_chart_label(_f)}**"))
                _sections.append(_iframe(_f.read_text(), height=900))

        if not _chord_files and not _heatmap_files:
            if _multi_charts and not _has_selection:
                _sections.append(
                    mo.callout(
                        mo.md("Select a row in the scoring table above to display the contact map."),
                        kind="info",
                    )
                )
            elif _contactmap_dir:
                _sections.append(
                    mo.callout(
                        mo.md(
                            "No chord chart or heatmap files found in "
                            f"`{_contactmap_dir}`. "
                            "Check that at least one of the output options is enabled."
                        ),
                        kind="warn",
                    )
                )

        # ── Alanine scan results (shown when the step directory exists) ───────
        import pandas as _pd

        _ala_dir = next(iter(sorted(_run_dir.glob("*_alascan"))), None)
        if _ala_dir is not None:
            _sections.append(mo.md("---\n## Alanine Scan"))
            _tsv_files = sorted(_ala_dir.glob("scan_clt_*.tsv"))
            _ala_html_files = sorted(_ala_dir.glob("scan_clt_*.html"))
            for _tsv in _tsv_files:
                _sections.append(mo.md(f"### Results — {_tsv.stem}"))
                try:
                    _df = _pd.read_csv(_tsv, sep="\t", comment="#")
                    _df = _df.drop(columns=["full_resname"], errors="ignore")
                    _float_cols = _df.select_dtypes(include="float").columns
                    _df[_float_cols] = _df[_float_cols].round(3)
                    _std_cols = [c for c in _df.columns if c.endswith("_std")]
                    _sections.append(
                        mo.ui.table(
                            _df,
                            pagination=False,
                            selection=None,
                            show_column_summaries=True,
                            hidden_columns=[] if alascan_show_std.value else _std_cols,
                        )
                    )
                except Exception as _e:
                    _sections.append(
                        mo.callout(
                            mo.md(f"Could not read `{_tsv.name}`: {_e}"), kind="warn"
                        )
                    )
            for _hf in _ala_html_files:
                import re as _re

                _hf_content = _hf.read_text()
                _m = _re.search(r'style="height:(\d+)px', _hf_content)
                _fig_height = int(_m.group(1)) + 30 if _m else 1030
                _sections.append(mo.md(f"### Plot — {_hf.stem}"))
                _sections.append(_iframe(_hf_content, height=_fig_height))
            if not _tsv_files and not _ala_html_files:
                _sections.append(
                    mo.callout(
                        mo.md(
                            f"No output files found in `{_ala_dir.name}/`. "
                            "Interface residues may not have been detected."
                        ),
                        kind="warn",
                    )
                )

        _out = mo.vstack(_sections)

    _out
    return


if __name__ == "__main__":
    app.run()
