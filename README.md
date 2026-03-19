# Clonalyzer 2

A browser-based tool for kinetics analysis of CHO cell fed-batch cultures. Upload a CSV, configure the exponential phase window, and instantly get publication-ready plots and downloadable results — no installation, no server, runs entirely in your browser.

**Live app → [ebalderasr.github.io/Clonalyzer-2](https://ebalderasr.github.io/Clonalyzer-2/)**

---

## Features

- **Zero installation** — runs fully client-side via [Pyodide](https://pyodide.org) (Python in WebAssembly)
- **Batch and fed-batch** support through the `is_post_feed` flag
- **Configurable exponential phase** — set start and end time independently
- **Four plot families** — scatter, mean ± SD lines, exponential-phase bar charts, and correlation plots
- **Custom correlations** — choose any two variables and generate on demand
- **Download ZIP** — all processed CSVs and PNG plots in a single file

---

## How to use

1. Open [ebalderasr.github.io/Clonalyzer-2](https://ebalderasr.github.io/Clonalyzer-2/)
2. Wait ~30 s for the Python environment to load (one-time per session)
3. Set the exponential phase window (default: 0 – 96 h)
4. Check or uncheck **Skip first row** depending on your CSV format
5. Drag & drop your CSV file or click to browse
6. Explore results in the tabbed plot viewer
7. Click **Download results (.zip)** to save everything

---

## Input CSV format

The tool accepts two CSV layouts:

### With metadata row (default — "Skip first row" checked)

```
Culture time (h), Clone ID, Biological replicate, ...   ← row 1: human labels (skipped)
t_hr, Clone, Rep, is_post_feed, VCD, ...                ← row 2: column names
0, Control, 1, FALSE, 1.80E+05, ...                     ← data
```

### Without metadata row ("Skip first row" unchecked)

```
t_hr, Clone, Rep, is_post_feed, VCD, ...                ← row 1: column names
0, Control, 1, FALSE, 1.80E+05, ...                     ← data
```

### Required columns

| Column | Description | Units |
|---|---|---|
| `t_hr` | Culture time | h |
| `Clone` | Clone identifier | text |
| `Rep` | Biological replicate | integer |
| `is_post_feed` | `TRUE` after a feed event, `FALSE` otherwise | boolean |
| `VCD` | Viable cell density | cells/mL |
| `DCD` | Dead cell density | cells/mL |
| `Viab_pct` | Viability | % |
| `rP_mg_L` | Recombinant protein titer | mg/L |
| `Glc_g_L` | Glucose concentration | g/L |
| `Lac_g_L` | Lactate concentration | g/L |
| `Gln_mM` | Glutamine concentration | mM |
| `Glu_mM` | Glutamate concentration | mM |
| `Vol_mL` | Culture volume | mL |

### Optional columns

If present, fluorescence kinetics and additional correlation plots are generated automatically.

| Column | Description | Units |
|---|---|---|
| `GFP_mean` | Mean GFP fluorescence | A.U. |
| `GFP_std` | GFP standard deviation | A.U. |
| `TMRM_mean` | Mean TMRM fluorescence | A.U. |
| `TMRM_std` | TMRM standard deviation | A.U. |

> European decimal commas (`1,5`) are accepted and converted automatically.

---

## Calculations

All specific rates use a **trapezoid approximation** for the average cell count between consecutive timepoints:

$$\bar{N} = \frac{1}{2}(VCD_1 \cdot V_1 + VCD_2 \cdot V_2) \times 10^3$$

### Specific rates computed

| Parameter | Formula | Units |
|---|---|---|
| μ | ln(VCD₂/VCD₁) / Δt | 1/h |
| qGlc | Δ(Glc·V) / Δt / N̄ | pmol/cell/day |
| qLac | Δ(Lac·V) / Δt / N̄ | pmol/cell/day |
| qGln | Δ(Gln·V) / Δt / N̄ | pmol/cell/h |
| qGlu | Δ(Glu·V) / Δt / N̄ | pmol/cell/h |
| qP | Δ(rP·V) / Δt / N̄ | pg/cell/day |
| Y Lac/Glc | ΔLac·V / ΔGlc·V | g/g |
| Y Glu/Gln | ΔGlu·V / ΔGln·V | mol/mol |
| IVCD | ∫ VCD dt (cumulative) | cells·h/mL |
| dGFP/dt | (GFP₂ − GFP₁) / Δt | A.U./h |
| dTMRM/dt | (TMRM₂ − TMRM₁) / Δt | A.U./h |

### Fed-batch handling

When `is_post_feed = TRUE`, the interval from the previous pre-feed row to that post-feed row is **skipped** from all rate calculations. The apparent change in metabolite concentrations at a feed event is caused by medium dilution, not cellular activity. All other intervals — including post-feed → next pre-feed — are used normally, since volume-corrected mass balances (S × V) account for any remaining dilution.

### Phase classification

Each row is assigned to a phase based on the configurable time window:

- **Exponential phase** — t_start ≤ t_hr ≤ t_end (default: 0 – 96 h)
- **Stationary phase** — all other timepoints

Exponential-phase summaries (bar charts) are computed as the mean of interval rates across all replicates within the defined window.

---

## Output

### ZIP contents

```
clonalyzer_results.zip
├── data_kinetics_processed.csv     ← full dataset with all computed parameters
├── data_exp_phase_summary.csv      ← per clone × replicate summary
└── plots/
    ├── 01_scatter/                 ← individual data points coloured by clone
    ├── 02_lines/                   ← clone mean ± SD across time
    ├── 03_bars_exp_phase/          ← exponential phase bar charts
    └── 04_correlations/            ← default correlation plots
```

Custom correlation plots generated in the browser are not included in the ZIP. Use your browser's right-click → "Save image" to export individual plots.

---

## Plot conventions

| Element | Encoding |
|---|---|
| **Color** | Clone |
| **Marker** | Phase — ● circle = Exponential, ▲ triangle = Stationary |
| **Regression line** | Solid = Exponential, Dashed = Stationary |
| **Error bands** | ± 1 SD across biological replicates |

---

## Project structure

```
Clonalyzer-2/
├── index.html       ← single-page app (Bootstrap 5)
├── app.js           ← Pyodide initialization, UI logic, ZIP generation
└── clonalyzer.py    ← Python analysis module (runs in browser via Pyodide)
```

The entire analysis pipeline runs client-side. No data is sent to any server.

---

## Related

[Clonalyzer](https://github.com/ebalderasr/clonalyzer) — command-line version of the same pipeline, for scripted or batch use.

---

## Author

**Emiliano Balderas Ramírez**
