/* ============================================================
   Clonalyzer-2  –  app.js
   Manages Pyodide initialization, file handling, UI state,
   interactive Plotly chart rendering, and ZIP download.
   ============================================================ */

"use strict";

// ── State ──────────────────────────────────────────────────────────────────────
let pyodide = null;
let results = null;

// ── DOM refs (populated after DOMContentLoaded) ────────────────────────────────
let elInitBar, elInitMsg, elInitSection;
let elUploadSection, elDropZone, elFileInput;
let elProcessSection, elProcessBar, elProcessMsg;
let elResultsSection;
let elExpPhaseStartInput, elExpPhaseInput, elSkipFirstRow, elUseVolume;
let elFluorGFP, elFluorTMRM, elFluorBODIPY, elFluorCellROX;

// ── Plotly config used for all interactive charts ──────────────────────────────
const PLOTLY_CONFIG = {
    responsive:              true,
    displayModeBar:          true,
    modeBarButtonsToRemove:  ["select2d", "lasso2d", "autoScale2d"],
    displaylogo:             false,
    toImageButtonOptions:    { format: "png", width: 900, height: 500, scale: 2 },
};

// ── Init ───────────────────────────────────────────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
    elInitBar        = document.getElementById("init-bar");
    elInitMsg        = document.getElementById("init-msg");
    elInitSection    = document.getElementById("init-section");
    elUploadSection  = document.getElementById("upload-section");
    elDropZone       = document.getElementById("drop-zone");
    elFileInput      = document.getElementById("file-input");
    elProcessSection = document.getElementById("process-section");
    elProcessBar     = document.getElementById("process-bar");
    elProcessMsg     = document.getElementById("process-msg");
    elResultsSection     = document.getElementById("results-section");
    elExpPhaseStartInput = document.getElementById("exp-phase-start");
    elExpPhaseInput      = document.getElementById("exp-phase-end");
    elSkipFirstRow       = document.getElementById("skip-first-row");
    elUseVolume          = document.getElementById("use-volume");
    elFluorGFP           = document.getElementById("fluor-gfp");
    elFluorTMRM          = document.getElementById("fluor-tmrm");
    elFluorBODIPY        = document.getElementById("fluor-bodipy");
    elFluorCellROX       = document.getElementById("fluor-cellrox");

    setupDropZone();
    setupTabResizeListener();
    initPyodide();
});

// ── Pyodide initialization ─────────────────────────────────────────────────────
async function initPyodide() {
    try {
        setInitProgress("Loading Python runtime…", 15);
        pyodide = await loadPyodide();

        setInitProgress("Installing scientific packages…", 40);
        await pyodide.loadPackage(["pandas", "numpy", "scipy"]);

        setInitProgress("Loading Clonalyzer…", 80);
        const resp = await fetch("clonalyzer.py");
        if (!resp.ok) throw new Error("Could not load clonalyzer.py");
        const code = await resp.text();
        await pyodide.runPythonAsync(code);

        setInitProgress("Ready!", 100);
        await sleep(400);
        show(elUploadSection);
        hide(elInitSection);
    } catch (err) {
        elInitMsg.textContent = "Error loading Python environment: " + err.message;
        elInitMsg.style.color = "#ef4444";
        console.error(err);
    }
}

function setInitProgress(msg, pct) {
    elInitMsg.textContent       = msg;
    elInitBar.style.width       = pct + "%";
    elInitBar.setAttribute("aria-valuenow", pct);
}

// ── Drop zone ──────────────────────────────────────────────────────────────────
function setupDropZone() {
    elDropZone.addEventListener("dragover", e => {
        e.preventDefault();
        elDropZone.classList.add("drag-over");
    });
    elDropZone.addEventListener("dragleave", () => {
        elDropZone.classList.remove("drag-over");
    });
    elDropZone.addEventListener("drop", e => {
        e.preventDefault();
        elDropZone.classList.remove("drag-over");
        const file = e.dataTransfer.files[0];
        if (file) handleFile(file);
    });
    elDropZone.addEventListener("click", () => elFileInput.click());
    elFileInput.addEventListener("change", e => {
        if (e.target.files[0]) handleFile(e.target.files[0]);
    });
}

// ── Resize Plotly charts when a tab becomes visible ───────────────────────────
function setupTabResizeListener() {
    document.getElementById("plot-tabs")?.addEventListener("shown.bs.tab", e => {
        const paneId = e.target.getAttribute("data-bs-target");
        const pane = document.querySelector(paneId);
        pane?.querySelectorAll(".js-plotly-plot").forEach(div => Plotly.Plots.resize(div));
    });
}

// ── File processing ────────────────────────────────────────────────────────────
async function handleFile(file) {
    if (!file.name.toLowerCase().endsWith(".csv")) {
        alert("Please select a CSV file.");
        return;
    }

    hide(elUploadSection);
    hide(elResultsSection);
    show(elProcessSection);
    setProcessProgress("Reading file…", 2);

    try {
        const csvText       = await file.text();
        const expPhaseStart = parseFloat(elExpPhaseStartInput.value) || 0.0;
        const expPhaseEnd   = parseFloat(elExpPhaseInput.value) || 96.0;
        const skipFirstRow  = elSkipFirstRow.checked;
        const useVolume     = elUseVolume.checked;

        // Build comma-separated list of enabled fluorescence channels
        const fluorParts = [];
        if (elFluorGFP?.checked)     fluorParts.push("GFP");
        if (elFluorTMRM?.checked)    fluorParts.push("TMRM");
        if (elFluorBODIPY?.checked)  fluorParts.push("BODIPY");
        if (elFluorCellROX?.checked) fluorParts.push("CellROX");
        const fluorChannels = fluorParts.join(",");  // "" means none selected → Python includes all

        pyodide.globals.set("_progress_cb",    (msg, pct) => setProcessProgress(msg, pct));
        pyodide.globals.set("_csv_text",        csvText);
        pyodide.globals.set("_exp_start",       expPhaseStart);
        pyodide.globals.set("_exp_end",         expPhaseEnd);
        pyodide.globals.set("_skip_first_row",  skipFirstRow);
        pyodide.globals.set("_use_volume",      useVolume);
        pyodide.globals.set("_fluor_channels",  fluorChannels);

        const pyResult = await pyodide.runPythonAsync(
            `run_analysis(_csv_text, _exp_start, _exp_end, _skip_first_row, _progress_cb, _use_volume, _fluor_channels)`
        );

        results = deepConvert(pyResult);

        displayResults(results);
        hide(elProcessSection);
        show(elResultsSection);
        elResultsSection.scrollIntoView({ behavior: "smooth" });

    } catch (err) {
        setProcessProgress("Error: " + err.message, 0);
        elProcessBar.classList.add("bg-danger");
        console.error(err);
        setTimeout(() => {
            hide(elProcessSection);
            show(elUploadSection);
            elProcessBar.classList.remove("bg-danger");
        }, 4000);
    }
}

function setProcessProgress(msg, pct) {
    elProcessMsg.textContent     = msg;
    elProcessBar.style.width     = pct + "%";
    elProcessBar.setAttribute("aria-valuenow", pct);
}

// ── Results rendering ──────────────────────────────────────────────────────────
function displayResults(r) {
    renderInfoCards(r.info);
    renderColorPickers(r.info.clones, r.info.palette || {});
    renderPlotTab("tab-scatter",      r.plots.scatter);
    renderPlotTab("tab-lines",        r.plots.lines);
    renderPlotTab("tab-bars",         r.plots.bars);
    renderPlotTab("tab-correlations", r.plots.correlations);
    document.querySelector("#plot-tabs .nav-link")?.click();
    populateCustomCorrSelects(r.avail_cols || []);
    populateMultiAxisTool(r.info.clones || [], r.avail_cols || []);
    document.getElementById("custom-corr-gallery").innerHTML = "";
    document.getElementById("multi-axis-gallery").innerHTML = "";
    // Resize first-tab charts after the browser has finished layout
    requestAnimationFrame(() => requestAnimationFrame(() => {
        document.querySelectorAll(".tab-pane.active .js-plotly-plot")
            .forEach(div => Plotly.Plots.resize(div));
    }));
}

// ── Clone color pickers ────────────────────────────────────────────────────────
function renderColorPickers(clones, palette) {
    const container = document.getElementById("clone-color-pickers");
    container.innerHTML = clones.map((c, i) => `
        <div class="d-flex align-items-center gap-2">
            <input type="color"
                   id="clone-color-${i}"
                   value="${palette[c] || "#000000"}"
                   class="clone-color-input"
                   title="Color for ${c}" />
            <span class="small fw-medium">${c}</span>
        </div>
    `).join("");
    show(document.getElementById("clone-colors-section"));
}

async function applyCloneColors() {
    const btn = document.getElementById("btn-apply-colors");
    btn.disabled = true;
    btn.textContent = "Applying…";
    try {
        const clones = results.info.clones;
        const palette = {};
        clones.forEach((c, i) => {
            const el = document.getElementById(`clone-color-${i}`);
            if (el) palette[c] = el.value;
        });

        pyodide.globals.set("_palette_json", JSON.stringify(palette));
        const pyResult = await pyodide.runPythonAsync("regenerate_plots(_palette_json)");
        const newPlots = deepConvert(pyResult);

        results.plots = newPlots;
        renderPlotTab("tab-scatter",      newPlots.scatter);
        renderPlotTab("tab-lines",        newPlots.lines);
        renderPlotTab("tab-bars",         newPlots.bars);
        renderPlotTab("tab-correlations", newPlots.correlations);

        // Restore active tab resize
        const activeBtn = document.querySelector("#plot-tabs .nav-link.active");
        if (activeBtn) {
            const paneId = activeBtn.getAttribute("data-bs-target");
            document.querySelector(paneId)
                ?.querySelectorAll(".js-plotly-plot")
                .forEach(div => Plotly.Plots.resize(div));
        }
        // Clear custom correlations (they'll be regenerated with new colors on demand)
        document.getElementById("custom-corr-gallery").innerHTML = "";
        document.getElementById("multi-axis-gallery").innerHTML = "";
    } catch (err) {
        alert("Error applying colors: " + err.message);
        console.error(err);
    } finally {
        btn.disabled = false;
        btn.textContent = "Apply to all charts";
    }
}

function renderInfoCards(info) {
    const el = document.getElementById("info-cards");
    const activeSet    = new Set(info.active_fluor  || []);
    const enabledFluor = info.enabled_fluor || [];
    let cytoTag = "";
    if (enabledFluor.length === 0) {
        cytoTag = `<span class="badge bg-secondary ms-2">No fluorescence</span>`;
    } else {
        cytoTag = enabledFluor.map(ch =>
            activeSet.has(ch)
                ? `<span class="badge bg-success ms-1">${ch} ✓</span>`
                : `<span class="badge bg-light border ms-1 text-muted">${ch} –</span>`
        ).join("");
    }
    const scenarioTag = info.scenario === "variable_volume"
        ? `<span class="badge bg-info text-dark ms-2">Variable volume</span>`
        : `<span class="badge bg-warning text-dark ms-2">Constant volume</span>`;

    el.innerHTML = `
        <div class="col-6 col-md-3">
            <div class="card text-center shadow-sm">
                <div class="card-body">
                    <div class="fs-1 fw-bold text-primary">${info.n_clones}</div>
                    <div class="text-muted small">Clones</div>
                    <div class="mt-1">${info.clones.join(", ")}</div>
                </div>
            </div>
        </div>
        <div class="col-6 col-md-3">
            <div class="card text-center shadow-sm">
                <div class="card-body">
                    <div class="fs-1 fw-bold text-primary">${info.n_reps}</div>
                    <div class="text-muted small">Replicates</div>
                </div>
            </div>
        </div>
        <div class="col-6 col-md-3">
            <div class="card text-center shadow-sm">
                <div class="card-body">
                    <div class="fs-1 fw-bold text-primary">${info.n_timepoints}</div>
                    <div class="text-muted small">Timepoints</div>
                </div>
            </div>
        </div>
        <div class="col-6 col-md-3">
            <div class="card text-center shadow-sm">
                <div class="card-body">
                    <div class="fs-1 fw-bold text-primary">${info.n_rows}</div>
                    <div class="text-muted small">Total rows ${cytoTag}</div>
                    <div class="mt-1">${scenarioTag}</div>
                </div>
            </div>
        </div>
    `;
}

/**
 * Render a tab of Plotly charts.
 * Keys starting with "__section__" are treated as section headers, not plots.
 * @param {string} tabId  - id of the inner <div> inside the tab pane
 * @param {Object} plotsData - { fname: { traces, layout } | null }
 */
function renderPlotTab(tabId, plotsData) {
    const container = document.getElementById(tabId);
    if (!container) return;

    const entries = Object.entries(plotsData || {}).filter(([k]) => !k.startsWith("__section__"));
    if (entries.length === 0) {
        container.innerHTML = `<p class="text-muted text-center py-4">No plots available for this category.</p>`;
        return;
    }

    // Build card grid, inserting section headers for __section__ sentinels
    const allEntries = Object.entries(plotsData || {});
    container.innerHTML = allEntries.map(([fname, spec]) => {
        if (fname.startsWith("__section__")) {
            const title = fname.slice("__section__".length);
            return `<div class="col-12"><div class="plots-section-header">${title}</div></div>`;
        }
        return `
        <div class="col-12 col-lg-6 mb-4">
            <div class="card shadow-sm h-100">
                <div class="card-header py-2 px-3 small text-muted">
                    ${fname.replace(/_/g, " ")}
                </div>
                <div class="card-body p-1">
                    <div id="plot-${tabId}-${fname}" class="plotly-chart"></div>
                </div>
            </div>
        </div>`;
    }).join("");

    // Render each Plotly chart (skip sentinels)
    for (const [fname, spec] of allEntries) {
        if (!fname.startsWith("__section__") && spec) {
            Plotly.newPlot(`plot-${tabId}-${fname}`, spec.traces, spec.layout, PLOTLY_CONFIG);
        }
    }
}

// ── ZIP download ───────────────────────────────────────────────────────────────
async function downloadZip() {
    if (!results) return;

    const btn = document.getElementById("btn-download");
    btn.disabled = true;

    try {
        const zip = new JSZip();
        zip.file("data_kinetics_processed.csv", results.processed_csv);
        zip.file("data_exp_phase_summary.csv",  results.summary_csv);

        const folders = {
            "plots/01_scatter":        results.plots.scatter,
            "plots/02_lines":          results.plots.lines,
            "plots/03_bars_exp_phase": results.plots.bars,
            "plots/04_correlations":   results.plots.correlations,
        };

        // Count total plots for progress display (exclude section sentinel keys)
        let total = 0;
        for (const plots of Object.values(folders))
            total += Object.keys(plots || {}).filter(k => !k.startsWith("__section__")).length;
        let done = 0;

        for (const [folder, plots] of Object.entries(folders)) {
            for (const [name, spec] of Object.entries(plots || {})) {
                if (name.startsWith("__section__")) continue;
                btn.textContent = `Exporting ${++done}/${total}…`;
                const b64 = await specToPng(spec);
                if (b64) zip.file(`${folder}/${name}.png`, b64, { base64: true });
            }
        }

        btn.textContent = "Compressing…";
        const blob = await zip.generateAsync({ type: "blob", compression: "DEFLATE" });
        saveAs(blob, "clonalyzer_results.zip");

    } catch (err) {
        alert("Error generating ZIP: " + err.message);
        console.error(err);
    } finally {
        btn.disabled    = false;
        btn.textContent = "⬇ Download results (.zip)";
    }
}

/**
 * Render a Plotly spec in an off-screen temporary div and return a base64 PNG string.
 */
async function specToPng(spec, width = 900, height = 500) {
    const tempDiv = document.createElement("div");
    tempDiv.style.cssText = `position:fixed;left:-9999px;top:0;width:${width}px;height:${height}px;`;
    document.body.appendChild(tempDiv);
    try {
        const exportLayout = { ...spec.layout, width, height, autosize: false };
        await Plotly.newPlot(tempDiv, spec.traces, exportLayout, { staticPlot: true });
        const dataUrl = await Plotly.toImage(tempDiv, { format: "png", width, height });
        return dataUrl.replace("data:image/png;base64,", "");
    } finally {
        Plotly.purge(tempDiv);
        document.body.removeChild(tempDiv);
    }
}

// ── Analyze another file ───────────────────────────────────────────────────────
function analyzeAnother() {
    // Destroy existing Plotly charts to free memory
    document.querySelectorAll(".js-plotly-plot").forEach(div => Plotly.purge(div));
    results = null;
    elFileInput.value = "";
    hide(elResultsSection);
    hide(document.getElementById("clone-colors-section"));
    document.getElementById("custom-corr-gallery").innerHTML = "";
    document.getElementById("multi-axis-gallery").innerHTML = "";
    show(elUploadSection);
    elUploadSection.scrollIntoView({ behavior: "smooth" });
}

// ── Custom correlations ────────────────────────────────────────────────────────
function populateCustomCorrSelects(availCols) {
    const xSel = document.getElementById("custom-x");
    const ySel = document.getElementById("custom-y");
    xSel.innerHTML = "";
    ySel.innerHTML = "";
    for (const [label, col] of availCols) {
        const makeOpt = sel => {
            const o = document.createElement("option");
            o.value = col;
            o.textContent = label;
            sel.appendChild(o);
        };
        makeOpt(xSel);
        makeOpt(ySel);
    }
    if (xSel.options.length > 0) xSel.selectedIndex = 0;
    if (ySel.options.length > 1) ySel.selectedIndex = 1;
}

async function generateCustomCorr() {
    const xSel   = document.getElementById("custom-x");
    const ySel   = document.getElementById("custom-y");
    const btn    = document.getElementById("btn-custom-corr");
    const xcol   = xSel.value;
    const ycol   = ySel.value;
    const xlabel = xSel.options[xSel.selectedIndex].text;
    const ylabel = ySel.options[ySel.selectedIndex].text;

    if (!xcol || !ycol) return;

    btn.disabled    = true;
    btn.textContent = "Generating…";

    try {
        pyodide.globals.set("_xcol", xcol);
        pyodide.globals.set("_ycol", ycol);
        const pyResult = await pyodide.runPythonAsync(
            `make_custom_correlation(_xcol, _ycol)`
        );
        const spec = deepConvert(pyResult);
        if (spec) addCustomCorrPlot(xlabel, ylabel, spec);
        else alert("Not enough data for this combination.");
    } catch (err) {
        alert("Error: " + err.message);
        console.error(err);
    } finally {
        btn.disabled    = false;
        btn.textContent = "Generate plot";
    }
}

function addCustomCorrPlot(xlabel, ylabel, spec) {
    const gallery   = document.getElementById("custom-corr-gallery");
    const id        = `cc-${Date.now()}`;
    const plotDivId = `plot-${id}`;

    const div = document.createElement("div");
    div.className = "col-12 col-lg-6 mb-4";
    div.id = id;
    div.innerHTML = `
        <div class="card shadow-sm h-100">
            <div class="card-header py-2 px-3 d-flex justify-content-between align-items-center">
                <span class="small text-muted fw-semibold">${xlabel} vs ${ylabel}</span>
                <button class="btn btn-sm btn-outline-danger py-0 px-2 lh-1"
                        onclick="Plotly.purge('${plotDivId}'); document.getElementById('${id}').remove()"
                        title="Remove">✕</button>
            </div>
            <div class="card-body p-1">
                <div id="${plotDivId}" class="plotly-chart"></div>
            </div>
        </div>`;
    gallery.prepend(div);

    Plotly.newPlot(plotDivId, spec.traces, spec.layout, PLOTLY_CONFIG);
}

// ── Multi-axis time comparison ────────────────────────────────────────────────
function populateMultiAxisTool(clones, availCols) {
    const cloneSel = document.getElementById("multi-clone");
    const sels = [
        document.getElementById("multi-y1"),
        document.getElementById("multi-y2"),
        document.getElementById("multi-y3"),
    ];

    cloneSel.innerHTML = "";
    {
        const opt = document.createElement("option");
        opt.value = "__all__";
        opt.textContent = "All clones";
        cloneSel.appendChild(opt);
    }
    for (const clone of clones) {
        const opt = document.createElement("option");
        opt.value = clone;
        opt.textContent = clone;
        cloneSel.appendChild(opt);
    }

    for (const sel of sels) sel.innerHTML = "";
    for (const [label, col] of availCols) {
        for (const sel of sels) {
            const opt = document.createElement("option");
            opt.value = col;
            opt.textContent = label;
            sel.appendChild(opt);
        }
    }

    const blankOpt = document.createElement("option");
    blankOpt.value = "";
    blankOpt.textContent = "None";
    sels[2].prepend(blankOpt);

    if (sels[0].options.length > 0) sels[0].selectedIndex = 0;
    if (sels[1].options.length > 1) sels[1].selectedIndex = 1;
    sels[2].selectedIndex = 0;

    ["multi-y1", "multi-y2", "multi-y3"].forEach(id => {
        document.getElementById(id)?.addEventListener("change", updateMultiAxisLabels);
    });
    updateMultiAxisLabels();
}

function updateMultiAxisLabels() {
    const pairs = [
        ["multi-y1", "multi-y1-label"],
        ["multi-y2", "multi-y2-label"],
        ["multi-y3", "multi-y3-label"],
    ];
    for (const [selId, labelId] of pairs) {
        const sel = document.getElementById(selId);
        const label = document.getElementById(labelId);
        if (!sel || !label) continue;
        label.textContent = sel.value
            ? (sel.options[sel.selectedIndex]?.text || "Axis")
            : "Unused axis";
    }
}

function readAxisRange(prefix) {
    const minEl = document.getElementById(`${prefix}-min`);
    const maxEl = document.getElementById(`${prefix}-max`);
    const min = minEl?.value?.trim() ? parseFloat(minEl.value) : null;
    const max = maxEl?.value?.trim() ? parseFloat(maxEl.value) : null;
    return [Number.isFinite(min) ? min : null, Number.isFinite(max) ? max : null];
}

async function generateMultiAxisPlot() {
    const btn = document.getElementById("btn-multi-axis");
    const clone = document.getElementById("multi-clone")?.value;
    const cols = [
        document.getElementById("multi-y1")?.value,
        document.getElementById("multi-y2")?.value,
        document.getElementById("multi-y3")?.value,
    ].filter(Boolean);

    const uniqueCols = [...new Set(cols)];
    if (!clone || uniqueCols.length < 2) {
        alert("Choose one clone and at least two different variables.");
        return;
    }

    const ranges = [
        readAxisRange("multi-y1"),
        readAxisRange("multi-y2"),
        readAxisRange("multi-y3"),
    ];

    btn.disabled = true;
    btn.textContent = "Generating…";

    try {
        pyodide.globals.set("_multi_clone", clone);
        pyodide.globals.set("_multi_cols_json", JSON.stringify(uniqueCols));
        pyodide.globals.set("_multi_ranges_json", JSON.stringify(ranges));
        const pyResult = await pyodide.runPythonAsync(
            `make_multi_axis_timeseries(_multi_clone, _multi_cols_json, _multi_ranges_json)`
        );
        const spec = deepConvert(pyResult);
        if (spec) {
            addMultiAxisPlot(clone, uniqueCols, spec);
        } else {
            alert("Not enough data for this multi-axis combination.");
        }
    } catch (err) {
        alert("Error: " + err.message);
        console.error(err);
    } finally {
        btn.disabled = false;
        btn.textContent = "Generate plot";
    }
}

function addMultiAxisPlot(clone, cols, spec) {
    const gallery = document.getElementById("multi-axis-gallery");
    const id = `ma-${Date.now()}`;
    const plotDivId = `plot-${id}`;
    const cloneLabel = clone === "__all__" ? "All clones" : clone;

    const labels = (spec.traces || []).map(t => t.name).join(" · ");
    const div = document.createElement("div");
    div.className = "col-12 mb-4";
    div.id = id;
    div.innerHTML = `
        <div class="card shadow-sm h-100">
            <div class="card-header py-2 px-3 d-flex justify-content-between align-items-center">
                <span class="small text-muted fw-semibold">${cloneLabel} · ${labels}</span>
                <button class="btn btn-sm btn-outline-danger py-0 px-2 lh-1"
                        onclick="Plotly.purge('${plotDivId}'); document.getElementById('${id}').remove()"
                        title="Remove">✕</button>
            </div>
            <div class="card-body p-1">
                <div id="${plotDivId}" class="plotly-chart"></div>
            </div>
        </div>`;
    gallery.prepend(div);

    Plotly.newPlot(plotDivId, spec.traces, spec.layout, PLOTLY_CONFIG);
}

// ── Utilities ──────────────────────────────────────────────────────────────────
function show(el) { el?.classList.remove("d-none"); }
function hide(el) { el?.classList.add("d-none"); }
function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }

// ── MathJax re-render when Methods modal opens ─────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
    document.getElementById("methodsModal")?.addEventListener("shown.bs.modal", () => {
        if (window.MathJax) MathJax.typesetPromise();
    });
});

/** Recursively convert Pyodide Map/PyProxy objects to plain JS objects/arrays. */
function deepConvert(obj) {
    if (obj && typeof obj.toJs === "function") {
        return deepConvert(obj.toJs({ dict_converter: Object.fromEntries }));
    }
    if (obj instanceof Map) {
        const out = {};
        for (const [k, v] of obj) out[k] = deepConvert(v);
        return out;
    }
    if (Array.isArray(obj)) return obj.map(deepConvert);
    if (obj !== null && typeof obj === "object") {
        const out = {};
        for (const [k, v] of Object.entries(obj)) out[k] = deepConvert(v);
        return out;
    }
    return obj;
}
