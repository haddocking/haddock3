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
    # HADDOCK3 Contact Map Analysis

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

    > **Note:** this notebook supports standard protein, DNA/RNA, and ion
    > inputs only. Small-molecule ligands requiring custom topology and
    > parameter files are **not** supported.

    **Requirements:** CNS must be installed and discoverable (see `haddock3-cfg cns_exec`).
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
    _std_aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
               "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
               "THR", "TRP", "TYR", "VAL"]

    _ncores_max = _os.cpu_count() or 1
    ncores_slider = mo.ui.slider(
        1, _ncores_max, value=_ncores_max, step=1,
        label="CPU cores",
        show_value=True,
    )
    chordchart_toggle = mo.ui.switch(label="Generate chord chart", value=True)
    heatmap_toggle = mo.ui.switch(label="Generate heatmap", value=True)
    alascan_toggle = mo.ui.switch(label="Run alanine scan", value=False)
    alascan_scan_residue = mo.ui.dropdown(
        options=_std_aa, value="ALA", label="Scan residue",
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
    )


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
    work_dir_input,
):
    _alascan_sub = (
        mo.vstack([alascan_scan_residue, alascan_plot_toggle, alascan_show_std], align="start")
        if alascan_toggle.value
        else mo.md("")
    )
    mo.vstack([
        mo.md("### Configuration"),
        mo.hstack([
            mo.vstack([
                mo.md("**Input file**"),
                pdb_upload,
            ], align="start"),
            mo.vstack([
                mo.md("**Working directory**"),
                work_dir_input,
            ], align="start"),
            mo.vstack([
                mo.md("**Execution**"),
                ncores_slider,
            ], align="start"),
            mo.vstack([
                mo.md("**contactmap outputs**"),
                chordchart_toggle,
                heatmap_toggle,
            ], align="start"),
            mo.vstack([
                mo.md("**alanine scan**"),
                alascan_toggle,
                _alascan_sub,
            ], align="start"),
        ], gap=2, justify="start"),
        mo.hstack([mo.md("**Advanced:**"), config_edit_toggle], justify="start"),
    ])
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
):
    def _make_cfg_str():
        _v = lambda b: "true" if b else "false"
        _lines = [
            "# HADDOCK3 interface analysis workflow — module parameters",
            "# CPU cores, mode, input file and run directory are set by the panel above.",
            "# Add, remove or reorder [module] sections to change the workflow steps.",
            "",
            "[topoaa]",
            "autohis = true",
            "",
            "[emscoring]",
            "",
            "[contactmap]",
            f"generate_chordchart = {_v(chordchart_toggle.value)}",
            f"generate_heatmap = {_v(heatmap_toggle.value)}",
            "single_model_analysis = false",
            "topX = 1",
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
        _out = mo.vstack([
            mo.callout(
                mo.md(
                    "**Config editor active.** Edit parameters below — these override the panel "
                    "settings for module-specific options. Changing a slider or toggle above will "
                    "regenerate this editor and discard manual edits."
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
    run_btn = mo.ui.run_button(label="▶  Run HADDOCK3 Workflow", kind="success")
    return (run_btn,)


@app.cell
def _(mo, run_btn):
    import sys as _sys
    import threading as _t
    if "_h3nb_stop" not in _sys.modules:
        _m = type(_sys)("_h3nb_stop")
        _m.event = _t.Event()
        _sys.modules["_h3nb_stop"] = _m
    stop_btn = mo.ui.run_button(label="⏹  Stop", kind="danger")
    if stop_btn.value:
        _sys.modules["_h3nb_stop"].event.set()
    mo.hstack([run_btn, stop_btn], gap=2, justify="start")
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
            mo.md("**No PDB file loaded.** Upload a PDB complex above, then click Run."),
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
            'color:#1a1a1a;font-family:monospace;padding:10px 14px;'
            'border-radius:6px;font-size:12px;white-space:pre-wrap;'
            'line-height:1.5;">'
            + _body.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            + "</div>"
        )
        mo.output.replace(
            mo.accordion({f"HADDOCK3 Log — {_label}": _panel})
        )

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
            mo.stop(True, mo.callout(
                mo.md(f"**Config parse error:** {_cfg_err}"),
                kind="danger",
            ))
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
            "emscoring.1": {**_exec},
            "contactmap.1": {
                "generate_chordchart": chordchart_toggle.value,
                "generate_heatmap": heatmap_toggle.value,
                "single_model_analysis": False,
                "topX": 1,
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
def _(alascan_show_std, mo, run_result):
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
        _contactmap_dir = _run_dir / "2_contactmap"

        _chord_files = sorted(_contactmap_dir.glob("*_chordchart.html"))
        _heatmap_files = sorted(_contactmap_dir.glob("*_heatmap.html"))

        _sections = [
            mo.callout(
                mo.md(f"Workflow completed. Results saved in `{_run_dir}`"),
                kind="success",
            )
        ]

        if _chord_files:
            _sections.append(mo.md("---\n## Chord Chart"))
            for _f in _chord_files:
                _sections.append(mo.md(f"**{_f.stem}**"))
                _sections.append(_iframe(_f.read_text(), height=900))

        if _heatmap_files:
            _sections.append(mo.md("---\n## Heatmap"))
            for _f in _heatmap_files:
                _sections.append(mo.md(f"**{_f.stem}**"))
                _sections.append(_iframe(_f.read_text(), height=900))

        if not _chord_files and not _heatmap_files:
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
        _ala_dir = _run_dir / "3_alascan"
        if _ala_dir.exists():
            _sections.append(mo.md("---\n## Alanine Scan"))
            _tsv_files = sorted(_ala_dir.glob("scan_clt_*.tsv"))
            _ala_html_files = sorted(_ala_dir.glob("scan_clt_*.html"))
            for _tsv in _tsv_files:
                _sections.append(mo.md(f"### Results — {_tsv.stem}"))
                try:
                    _df = _pd.read_csv(_tsv, sep="\t", comment="#")
                    # drop the redundant combined label column
                    _df = _df.drop(columns=["full_resname"], errors="ignore")
                    _float_cols = _df.select_dtypes(include="float").columns
                    _df[_float_cols] = _df[_float_cols].round(3)
                    _std_cols = [c for c in _df.columns if c.endswith("_std")]
                    _sections.append(mo.ui.table(
                        _df,
                        pagination=False,
                        selection=None,
                        show_column_summaries=True,
                        hidden_columns=[] if alascan_show_std.value else _std_cols,
                    ))
                except Exception as _e:
                    _sections.append(mo.callout(mo.md(f"Could not read `{_tsv.name}`: {_e}"), kind="warn"))
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
                        mo.md("No output files found in `3_alascan/`. Interface residues may not have been detected."),
                        kind="warn",
                    )
                )

        _out = mo.vstack(_sections)

    _out
    return



if __name__ == "__main__":
    app.run()
