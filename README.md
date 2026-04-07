<div align="center">

# Clonalyzer 2

### Fed-batch kinetics analysis for CHO cell cultures

<br>

**[→ Open the live app](https://ebalderasr.github.io/Clonalyzer-2/)**

<br>

[![Stack](https://img.shields.io/badge/Stack-Pyodide_·_Pandas_·_Matplotlib-4A90D9?style=for-the-badge)]()
[![Focus](https://img.shields.io/badge/Focus-Fed--batch_Kinetics_·_CHO-34C759?style=for-the-badge)]()
[![License](https://img.shields.io/badge/License-MIT-blue?style=for-the-badge)](./LICENSE)
[![Part of](https://img.shields.io/badge/Part_of-HostCell_Lab_Suite-5856D6?style=for-the-badge)](https://github.com/ebalderasr)

</div>

---

## What is Clonalyzer 2?

Clonalyzer 2 is a **browser-based kinetics analysis tool** for CHO fed-batch cultures. Upload a raw measurement CSV and it computes over 20 kinetic parameters — specific rates, yields, integral cell counts — then generates publication-ready plots and packages everything into a downloadable ZIP.

It runs entirely in the browser through [Pyodide](https://pyodide.org), the full scientific Python stack compiled to WebAssembly. No Python environment to install, no scripts to run, no data leaves your machine.

---

## Why it matters

After each fed-batch run, extracting kinetic parameters from raw viability counter and metabolite analyzer data typically means hours of manual spreadsheet work — repeated for every clone and replicate. Without a dedicated tool:

- Specific rates must be calculated interval by interval from raw concentration and volume data
- Feed events require manual correction to avoid confounding dilution with cellular activity
- Clone comparisons and replicate statistics are assembled by hand from separate files

Clonalyzer 2 automates the entire pipeline in a single CSV drop.

---

## How it works

### 1. Upload

Drag and drop your CSV file or click to browse. Two layouts are accepted:

**With metadata row** *(default — toggle "Skip first row")*
```
Culture time (h),  Clone ID,  Biological replicate, ...   ← ignored
t_hr,              Clone,     Rep,                  ...   ← column names
0,                 Control,   1,                    ...   ← data
```

**Without metadata row**
```
t_hr,  Clone,  Rep,  ...   ← column names
0,     Control,  1,  ...   ← data
```

> European decimal commas (`1,5`) are accepted and converted automatically.

### 2. Set the exponential phase window

Define the start and end time (h) of the exponential phase. All specific rates are summarized separately for the exponential and stationary phases. Default: 0–96 h.

### 3. Explore results

Results are organized across four tabbed plot families:

| Tab | Contents |
|---|---|
| Scatter | Individual data points colored by clone |
| Mean ± SD | Time-resolved profiles with error bars across biological replicates |
| Exponential phase bars | Clone-level comparison of exponential-phase parameter means ± SD |
| Correlations | Custom pairwise plots segmented by clone and phase |

### 4. Download

Click **Download results (.zip)** to export all CSVs and plots:

```
clonalyzer_results.zip
├── data_kinetics_processed.csv      ← all computed parameters, one row per sample
├── data_exp_phase_summary.csv       ← per clone × replicate exponential-phase means
└── plots/
    ├── 01_scatter/
    ├── 02_lines/
    ├── 03_bars_exp_phase/
    └── 04_correlations/
```

---

## Methods

Clonalyzer 2 computes specific rates interval by interval between consecutive timepoints using a trapezoidal approximation. Two scenarios are supported depending on whether culture volume data are available.

### Specific growth rate

μ is always calculated from viable cell density, regardless of the volume scenario:

$$\mu = \frac{\ln(VCD_2 / VCD_1)}{\Delta t} \quad \text{[h}^{-1}\text{]}$$

> Using total cell counts (TC = VCD × V) for μ would produce artefactual negative values whenever volume drops during sampling. VCD is the correct basis.

---

### Scenario A — Constant volume (concentration-based)

Used when `Vol_mL` is absent or the *Variable volume* checkbox is disabled. Rates are normalized by the **Integral Viable Cell Density (IVCD)**:

$$\Delta IVCD = \frac{VCD_1 + VCD_2}{2} \cdot \Delta t \quad \left[\frac{\text{cells} \cdot \text{h}}{\text{mL}}\right]$$

The specific rate of metabolite $i$ is the concentration change per unit of biomass exposure:

$$q_i = \frac{\Delta C_i}{\Delta IVCD} \quad \text{[pg or pmol / cell / day]}$$

| Parameter | Numerator | Sign convention |
|---|---|---|
| qGlc | $Glc_1 - Glc_2$ | + consumed, reported as pmol/cell/day |
| qLac | $Lac_2 - Lac_1$ | + produced, reported as pmol/cell/day |
| qP | $rP_2 - rP_1$ | + produced |
| qGln | $Gln_1 - Gln_2$ | + consumed, reported as pmol/cell/day |
| qGlu | $Glu_2 - Glu_1$ | + produced, reported as pmol/cell/day |

Metabolic yields are the ratio of concentration changes:

$$Y_{Lac/Glc} = \frac{Lac_2 - Lac_1}{Glc_1 - Glc_2} \quad Y_{Glu/Gln} = \frac{Glu_2 - Glu_1}{Gln_1 - Gln_2}$$

---

### Scenario B — Variable volume (mass-balance)

Used when `Vol_mL` is present and the *Variable volume* checkbox is enabled. Volume changes between timepoints — due to sampling, evaporation, or feeding — are accounted for explicitly via a mass balance.

**Total viable cells** in the reactor at each timepoint:

$$TC = VCD \times V \quad \text{[cells]}$$

**Integral Total Viable Cells (ITVC):**

$$\Delta ITVC = \frac{TC_1 + TC_2}{2} \cdot \Delta t \quad \text{[cells} \cdot \text{h]}$$

**Total mass** of each metabolite in the reactor:

$$M_i = C_i \times \frac{V}{1000} \quad \text{[g or mmol]}$$

The specific rate is the change in total mass normalized by ITVC:

$$q_i = \frac{\Delta M_i}{\Delta ITVC} \quad \text{[pg or pmol / cell / day]}$$

| Parameter | Numerator | Sign convention |
|---|---|---|
| qGlc | $M_{Glc,1} - M_{Glc,2}$ | + consumed, reported as pmol/cell/day |
| qLac | $M_{Lac,2} - M_{Lac,1}$ | + produced, reported as pmol/cell/day |
| qP | $M_{rP,2} - M_{rP,1}$ | + produced |
| qGln | $M_{Gln,1} - M_{Gln,2}$ | + consumed, reported as pmol/cell/day |
| qGlu | $M_{Glu,2} - M_{Glu,1}$ | + produced, reported as pmol/cell/day |

Metabolic yields are the ratio of mass changes:

$$Y_{Lac/Glc} = \frac{M_{Lac,2} - M_{Lac,1}}{M_{Glc,1} - M_{Glc,2}} \quad Y_{Glu/Gln} = \frac{M_{Glu,2} - M_{Glu,1}}{M_{Gln,1} - M_{Gln,2}}$$

---

### Fed-batch correction

Intervals where `is_post_feed` transitions from `FALSE` to `TRUE` are excluded from rate calculations. The apparent concentration change at that point reflects medium addition, not cellular activity. All other intervals are used normally; the mass-balance approach ($M = C \times V$) corrects for any remaining dilution in variable-volume runs.

---

## Features

| | |
|---|---|
| **Zero installation** | Runs fully client-side via Pyodide — no Python, no pip, no server |
| **20+ kinetic parameters** | μ, qGlc, qLac, qP, qGln, qGlu, yields, IVCD/ITVC, fluorescence kinetics |
| **Two calculation scenarios** | Constant-volume (ΔC/ΔIVCD) or variable-volume (ΔM/ΔITVC) — auto-detected from the input data |
| **4 plot families** | Scatter · Mean ± SD · Exponential phase bars · Correlations |
| **Fed-batch aware** | Feed events excluded automatically from rate calculations |
| **Configurable phase window** | Set exponential phase start and end independently |
| **Custom correlations** | Build any pairwise plot on demand, segmented by clone and phase |
| **Download ZIP** | All CSVs and PNGs packaged in one click |
| **No data upload** | Everything runs locally in the browser — no data leaves your machine |

---

## Input format

### Required columns

| Column | Description | Units |
|---|---|---|
| `t_hr` | Culture time | h |
| `Clone` | Clone identifier | — |
| `Rep` | Biological replicate number | — |
| `is_post_feed` | `TRUE` if sampled after a feed event | boolean |
| `VCD` | Viable cell density | cells/mL |
| `DCD` | Dead cell density | cells/mL |
| `Viab_pct` | Viability | % |
| `rP_mg_L` | Recombinant protein titer | mg/L |
| `Glc_g_L` | Glucose | g/L |
| `Lac_g_L` | Lactate | g/L |
| `Gln_mM` | Glutamine | mM |
| `Glu_mM` | Glutamate | mM |

### Optional columns

If present, the corresponding features are enabled automatically.

| Column | Description |
|---|---|
| `Vol_mL` | Culture volume (mL) — enables variable-volume mass-balance calculations |
| `GFP_mean` / `GFP_std` | GFP fluorescence intensity (A.U.) |
| `TMRM_mean` / `TMRM_std` | TMRM fluorescence intensity (A.U.) |

---

## Tech stack

**Analysis (in-browser via WebAssembly)**

![Pyodide](https://img.shields.io/badge/Pyodide-3776AB?style=flat-square&logo=python&logoColor=white)
![Pandas](https://img.shields.io/badge/Pandas-150458?style=flat-square&logo=pandas&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-013243?style=flat-square&logo=numpy&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557C?style=flat-square&logo=python&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?style=flat-square&logo=scipy&logoColor=white)

**Frontend**

![HTML5](https://img.shields.io/badge/HTML5-E34F26?style=flat-square&logo=html5&logoColor=white)
![CSS3](https://img.shields.io/badge/CSS3-1572B6?style=flat-square&logo=css3&logoColor=white)
![JavaScript](https://img.shields.io/badge/JavaScript-F7DF1E?style=flat-square&logo=javascript&logoColor=black)

**Packaging**

![JSZip](https://img.shields.io/badge/JSZip-555555?style=flat-square)
![FileSaver](https://img.shields.io/badge/FileSaver.js-555555?style=flat-square)

---

## Project structure

```
Clonalyzer-2/
├── index.html          ← markup only
├── style.css           ← all custom styles
├── app.js              ← Pyodide init, UI logic, ZIP generation
├── clonalyzer.py       ← analysis pipeline (runs in-browser via Pyodide)
└── Fig/                ← screenshot assets for this README
```

---

## Author

**Emiliano Balderas Ramírez**
Bioengineer · PhD Candidate in Biochemical Sciences
Instituto de Biotecnología (IBt), UNAM

[![LinkedIn](https://img.shields.io/badge/LinkedIn-emilianobalderas-0A66C2?style=flat-square&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/emilianobalderas/)
[![Email](https://img.shields.io/badge/Email-ebalderas%40live.com.mx-D14836?style=flat-square&logo=gmail&logoColor=white)](mailto:ebalderas@live.com.mx)

---

## Related

[**CellSplit**](https://github.com/ebalderasr/CellSplit) — Neubauer cell counting and passage planning for CHO cultures.

[**CellBlock**](https://github.com/ebalderasr/CellBlock) — shared biosafety cabinet scheduling for cell culture research groups.

---

<div align="center"><i>Clonalyzer 2 — drop a CSV, get your kinetics.</i></div>
