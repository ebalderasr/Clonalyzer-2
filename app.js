/* ============================================================
   Clonalyzer-2  –  app.js
   Manages Pyodide initialization, file handling, UI state,
   plot rendering, and ZIP download.
   ============================================================ */

"use strict";

// ── State ──────────────────────────────────────────────────────────────────────
let pyodide   = null;
let results   = null;

// ── DOM refs (populated after DOMContentLoaded) ────────────────────────────────
let elInitBar, elInitMsg, elInitSection;
let elUploadSection, elDropZone, elFileInput;
let elProcessSection, elProcessBar, elProcessMsg;
let elResultsSection;
let elExpPhaseStartInput, elExpPhaseInput, elSkipFirstRow;

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

    setupDropZone();
    initPyodide();
});

// ── Pyodide initialization ─────────────────────────────────────────────────────
async function initPyodide() {
    try {
        setInitProgress("Loading Python runtime…", 15);
        pyodide = await loadPyodide();

        setInitProgress("Installing scientific packages…", 40);
        await pyodide.loadPackage(["pandas", "numpy", "matplotlib", "scipy"]);

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

// ── File processing ────────────────────────────────────────────────────────────
async function handleFile(file) {
    if (!file.name.endsWith(".csv")) {
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

        // Pass progress callback from JS into Python
        pyodide.globals.set("_progress_cb",   (msg, pct) => setProcessProgress(msg, pct));
        pyodide.globals.set("_csv_text",       csvText);
        pyodide.globals.set("_exp_start",      expPhaseStart);
        pyodide.globals.set("_exp_end",        expPhaseEnd);
        pyodide.globals.set("_skip_first_row", skipFirstRow);

        const pyResult = await pyodide.runPythonAsync(
            `run_analysis(_csv_text, _exp_start, _exp_end, _skip_first_row, _progress_cb)`
        );

        // Convert Python proxy → JS object (nested)
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
    renderPlotTab("tab-scatter",      r.plots.scatter,      "Scatter");
    renderPlotTab("tab-lines",        r.plots.lines,        "Lines");
    renderPlotTab("tab-bars",         r.plots.bars,         "Bars");
    renderPlotTab("tab-correlations", r.plots.correlations, "Correlations");
    // Activate first tab
    document.querySelector("#plot-tabs .nav-link")?.click();
}

function renderInfoCards(info) {
    const el = document.getElementById("info-cards");
    const cytoTag = info.has_cyto
        ? `<span class="badge bg-success ms-2">GFP / TMRM ✓</span>`
        : `<span class="badge bg-secondary ms-2">No cytometry</span>`;

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
                </div>
            </div>
        </div>
    `;
}

function renderPlotTab(tabId, plots, label) {
    const container = document.getElementById(tabId);
    if (!container) return;

    const entries = Object.entries(plots || {});
    if (entries.length === 0) {
        container.innerHTML = `<p class="text-muted text-center py-4">No plots available for this category.</p>`;
        return;
    }

    container.innerHTML = entries.map(([fname, b64]) => `
        <div class="col-12 col-lg-6 mb-4">
            <div class="card shadow-sm h-100">
                <div class="card-header py-2 px-3 small text-muted">${fname.replace(/_/g," ")}</div>
                <div class="card-body p-2 text-center">
                    <img src="data:image/png;base64,${b64}"
                         class="img-fluid plot-img"
                         alt="${fname}"
                         loading="lazy">
                </div>
            </div>
        </div>
    `).join("");
}

// ── ZIP download ───────────────────────────────────────────────────────────────
async function downloadZip() {
    if (!results) return;

    const btn = document.getElementById("btn-download");
    btn.disabled    = true;
    btn.textContent = "Preparing ZIP…";

    try {
        const zip = new JSZip();
        zip.file("data_kinetics_processed.csv", results.processed_csv);
        zip.file("data_exp_phase_summary.csv",  results.summary_csv);

        const folders = {
            "plots/01_scatter":         results.plots.scatter,
            "plots/02_lines":           results.plots.lines,
            "plots/03_bars_exp_phase":  results.plots.bars,
            "plots/04_correlations":    results.plots.correlations,
        };
        for (const [folder, plots] of Object.entries(folders)) {
            for (const [name, b64] of Object.entries(plots || {})) {
                zip.file(`${folder}/${name}.png`, b64, { base64: true });
            }
        }

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

// ── Analyze another file ───────────────────────────────────────────────────────
function analyzeAnother() {
    results = null;
    elFileInput.value = "";
    hide(elResultsSection);
    show(elUploadSection);
    elUploadSection.scrollIntoView({ behavior: "smooth" });
}

// ── Utilities ──────────────────────────────────────────────────────────────────
function show(el) { el?.classList.remove("d-none"); }
function hide(el) { el?.classList.add("d-none"); }
function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }

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
