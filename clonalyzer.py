"""
Clonalyzer processing module – runs inside Pyodide (browser WebAssembly).

Returns all results as a Python dict that JavaScript converts to a JS object.
Plot data is returned as Plotly-compatible {traces, layout} dicts for
interactive rendering in the browser via Plotly.js.
"""

import io
import math
import numpy as np
import pandas as pd
from scipy import stats

# ── Constants ──────────────────────────────────────────────────────────────────
H_PER_DAY    = 24.0
PG_PER_G     = 1e12
MMOL_TO_PMOL = 1e9
X_RANGE      = (-10, 270)

# Molar masses for pmol conversions (g/mol)
MW_GLC = 180.16   # glucose
MW_LAC = 90.08    # lactate / lactic acid

PALETTE = ["#000000","#FF0066","#107F80","#F0E442",
           "#0072B2","#D55E00","#CC79A7","#999999"]

# Columns that must be present in every dataset
REQUIRED_COLS = ["t_hr","Clone","Rep",
                 "VCD","DCD","Viab_pct","rP_mg_L",
                 "Glc_g_L","Lac_g_L","Gln_mM","Glu_mM"]

# Vol_mL is OPTIONAL — its presence/absence determines the calculation scenario:
#   present  → variable-volume: mass balances with ITVC normalization
#   absent   → constant-volume: concentration-based with IVCD normalization
VOL_COL = "Vol_mL"

# is_post_feed is optional; defaults to False when absent (batch/constant-volume datasets)
FEED_COL = "is_post_feed"

OPTIONAL_COLS = [
    "GFP_mean","GFP_std","TMRM_mean","TMRM_std",
    "Bodipy_mean","Bodipy_std","CellRox_mean","CellRox_std",
]

SCENARIO_VAR   = "variable_volume"   # Vol_mL present
SCENARIO_CONST = "constant_volume"   # Vol_mL absent

PHASE_EXP  = "Exponential"
PHASE_STAT = "Stationary"

# output column names
MU        = "mu_per_h"
QGLC      = "qGlc_pg_cell_day"
QLAC      = "qLac_pg_cell_day"
QGLC_PMOL = "qGlc_pmol_cell_day"
QLAC_PMOL = "qLac_pmol_cell_day"
QP        = "qP_pg_cell_day"
QGLN_H    = "qGln_pmol_cell_h"
QGLN_D    = "qGln_pmol_cell_day"
QGLU_H    = "qGlu_pmol_cell_h"
QGLU_D    = "qGlu_pmol_cell_day"
Y_LG      = "Y_Lac_per_Glc_g_per_g"
Y_GQ      = "Y_Glu_per_Gln_mol_per_mol"
IVCD_INT  = "IVCD_interval_cells_h_per_mL"   # VCD trapezoid — constant-volume normalizer
IVCD_CUM  = "IVCD_cum_cells_h_per_mL"
ITVC_INT  = "ITVC_interval_cells_h"          # TC trapezoid  — variable-volume normalizer
ITVC_CUM  = "ITVC_cum_cells_h"
DGFP      = "dGFP_dt"
DTMRM     = "dTMRM_dt"
DBODIPY   = "dBodipy_dt"
DCELLROX  = "dCellRox_dt"
PHASE     = "Fase_Cultivo"

# Fluorescence channels: (intensity column, rate column, display label)
FLUOR_CHANNELS = [
    ("GFP_mean",     DGFP,     "GFP"),
    ("TMRM_mean",    DTMRM,    "TMRM"),
    ("Bodipy_mean",  DBODIPY,  "BODIPY"),
    ("CellRox_mean", DCELLROX, "CellROX"),
]

# plot specs: (column, y-label, filename-stem, ylim or None)
TS_SPECS = [
    ("VCD",          "VCD (cells/mL)",           "01_VCD",         None),
    ("Viab_pct",     "Viability (%)",             "02_Viability",   [0, 100]),
    (MU,             "μ (1/h)",                   "03_Mu",          None),
    ("Glc_g_L",      "Glucose (g/L)",             "04_Glc",         None),
    ("Lac_g_L",      "Lactate (g/L)",             "05_Lac",         None),
    ("Gln_mM",       "Glutamine (mM)",            "06_Gln",         None),
    ("Glu_mM",       "Glutamate (mM)",            "07_Glu",         None),
    ("rP_mg_L",      "rP Titer (mg/L)",           "08_Product",     None),
    (QGLC_PMOL,      "qGlc (pmol/cell/day)",      "09_qGlc",        None),
    (QLAC_PMOL,      "qLac (pmol/cell/day)",      "10_qLac",        None),
    (QGLN_D,         "qGln (pmol/cell/day)",      "11_qGln",        None),
    (QGLU_D,         "qGlu (pmol/cell/day)",      "12_qGlu",        None),
    (QP,             "qP (pg/cell/day)",          "13_qP",          None),
    (Y_LG,           "Y Lac/Glc (g/g)",           "14_YLacGlc",     None),
    (Y_GQ,           "Y Glu/Gln (mol/mol)",       "15_YGluGln",     None),
    (IVCD_CUM,       "IVCD (cells·h/mL)",          "16_IVCD",        None),
    (ITVC_CUM,       "ITVC (cells·h)",             "16b_ITVC",       None),
    ("GFP_mean",     "GFP intensity (A.U.)",      "17_GFP",         None),
    ("TMRM_mean",    "TMRM intensity (A.U.)",     "18_TMRM",        None),
    (DGFP,           "dGFP/dt (A.U./h)",          "19_dGFP_dt",     None),
    (DTMRM,          "dTMRM/dt (A.U./h)",         "20_dTMRM_dt",    None),
    ("Bodipy_mean",  "BODIPY intensity (A.U.)",   "21_Bodipy",      None),
    ("CellRox_mean", "CellROX intensity (A.U.)",  "22_CellRox",     None),
    (DBODIPY,        "dBODIPY/dt (A.U./h)",       "23_dBodipy_dt",  None),
    (DCELLROX,       "dCellROX/dt (A.U./h)",      "24_dCellRox_dt", None),
]

BAR_SPECS = [
    ("mu_exp",           "μ exp (1/h)",               "03_Mu_Exp"),
    ("qGlc_pmol_exp",    "qGlc exp (pmol/cell/day)",  "09_qGlc_Exp"),
    ("qLac_pmol_exp",    "qLac exp (pmol/cell/day)",  "10_qLac_Exp"),
    ("qP_exp",           "qP exp (pg/cell/day)",      "13_qP_Exp"),
    ("qGln_exp",         "qGln exp (pmol/cell/day)",  "11_qGln_Exp"),
    ("qGlu_exp",         "qGlu exp (pmol/cell/day)",  "12_qGlu_Exp"),
    ("Y_Lac_Glc_exp",    "Y Lac/Glc exp (g/g)",       "14_YLacGlc_Exp"),
    ("Y_Glu_Gln_exp",    "Y Glu/Gln exp (mol/mol)",   "15_YGluGln_Exp"),
    ("dGFP_dt_exp",      "dGFP/dt exp (A.U./h)",      "19_dGFP_dt_Exp"),
    ("dTMRM_dt_exp",     "dTMRM/dt exp (A.U./h)",     "20_dTMRM_dt_Exp"),
    ("dBODIPY_dt_exp",   "dBODIPY/dt exp (A.U./h)",   "21_dBodipy_dt_Exp"),
    ("dCellROX_dt_exp",  "dCellROX/dt exp (A.U./h)",  "22_dCellRox_dt_Exp"),
]

CORR_SPECS = [
    # GFP / TMRM
    (QP,           "GFP_mean",   "qP (pg/cell/day)",      "GFP (A.U.)",           "40_qP_vs_GFP"),
    (QP,           "TMRM_mean",  "qP (pg/cell/day)",      "TMRM (A.U.)",          "41_qP_vs_TMRM"),
    ("GFP_mean",   "TMRM_mean",  "GFP (A.U.)",            "TMRM (A.U.)",          "42_GFP_vs_TMRM"),
    (QGLC_PMOL,    "TMRM_mean",  "qGlc (pmol/cell/day)",  "TMRM (A.U.)",          "50_qGlc_vs_TMRM"),
    (QP,           QGLC_PMOL,    "qP (pg/cell/day)",      "qGlc (pmol/cell/day)", "52_qP_vs_qGlc"),
    (QLAC_PMOL,    QGLC_PMOL,    "qLac (pmol/cell/day)",  "qGlc (pmol/cell/day)", "55_qLac_vs_qGlc"),
    (QP,           DGFP,         "qP (pg/cell/day)",      "dGFP/dt (A.U./h)",     "60_qP_vs_dGFPdt"),
    # BODIPY
    (QP,           "Bodipy_mean","qP (pg/cell/day)",      "BODIPY (A.U.)",        "61_qP_vs_Bodipy"),
    (QGLC_PMOL,    "Bodipy_mean","qGlc (pmol/cell/day)",  "BODIPY (A.U.)",        "62_qGlc_vs_Bodipy"),
    ("GFP_mean",   "Bodipy_mean","GFP (A.U.)",            "BODIPY (A.U.)",        "63_GFP_vs_Bodipy"),
    (QP,           DBODIPY,      "qP (pg/cell/day)",      "dBODIPY/dt (A.U./h)",  "64_qP_vs_dBodipydt"),
    # CellROX
    (QP,           "CellRox_mean","qP (pg/cell/day)",     "CellROX (A.U.)",       "65_qP_vs_CellRox"),
    (QGLC_PMOL,    "CellRox_mean","qGlc (pmol/cell/day)", "CellROX (A.U.)",       "66_qGlc_vs_CellRox"),
    ("Bodipy_mean","CellRox_mean","BODIPY (A.U.)",        "CellROX (A.U.)",       "67_Bodipy_vs_CellRox"),
    (QP,           DCELLROX,     "qP (pg/cell/day)",      "dCellROX/dt (A.U./h)", "68_qP_vs_dCellRoxdt"),
]

# phase → (Plotly marker symbol, Plotly line dash, short label)
_PHASE_STYLE = {
    PHASE_EXP:  ("circle",      "solid", "Exp"),
    PHASE_STAT: ("triangle-up", "dash",  "Stat"),
}

# Column catalogue for custom correlations
CUSTOM_CORR_COLS = [
    ("VCD (cells/mL)",           "VCD"),
    ("Viability (%)",            "Viab_pct"),
    ("μ (1/h)",                  MU),
    ("Glucose (g/L)",            "Glc_g_L"),
    ("Lactate (g/L)",            "Lac_g_L"),
    ("Glutamine (mM)",           "Gln_mM"),
    ("Glutamate (mM)",           "Glu_mM"),
    ("rP Titer (mg/L)",          "rP_mg_L"),
    ("qGlc (pmol/cell/day)",     QGLC_PMOL),
    ("qLac (pmol/cell/day)",     QLAC_PMOL),
    ("qGln (pmol/cell/day)",     QGLN_D),
    ("qGlu (pmol/cell/day)",     QGLU_D),
    ("qP (pg/cell/day)",         QP),
    ("Y Lac/Glc (g/g)",          Y_LG),
    ("Y Glu/Gln (mol/mol)",      Y_GQ),
    ("IVCD (cells·h/mL)",        IVCD_CUM),
    ("ITVC (cells·h)",           ITVC_CUM),
    ("GFP intensity (A.U.)",     "GFP_mean"),
    ("TMRM intensity (A.U.)",    "TMRM_mean"),
    ("dGFP/dt (A.U./h)",         DGFP),
    ("dTMRM/dt (A.U./h)",        DTMRM),
    ("BODIPY intensity (A.U.)",  "Bodipy_mean"),
    ("CellROX intensity (A.U.)", "CellRox_mean"),
    ("dBODIPY/dt (A.U./h)",      DBODIPY),
    ("dCellROX/dt (A.U./h)",     DCELLROX),
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def _palette(clones):
    return {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(clones)}

def _active_fluor(df, fluor_set=None):
    """Return list of (mean_col, rate_col, label) for channels that are both
    present in the dataframe AND enabled by the user selection (fluor_set).
    fluor_set=None means all channels that have data are included."""
    return [
        (mc, rc, lbl) for mc, rc, lbl in FLUOR_CHANNELS
        if mc in df.columns and (fluor_set is None or lbl in fluor_set)
    ]

def _v(x):
    """Scalar → Python float; NaN/inf/None → None (JSON-safe)."""
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(f) or math.isinf(f)) else f

def _to_list(series):
    """Pandas Series → plain Python list with NaN/inf replaced by None."""
    return [_v(x) for x in series.tolist()]

def _std0(x):
    """Std value: NaN → 0.0 (for error bar arrays)."""
    v = _v(x)
    return 0.0 if v is None else v


# ── Loader ─────────────────────────────────────────────────────────────────────

def _load(csv_text: str, skip_first_row: bool = True) -> pd.DataFrame:
    header_row = 1 if skip_first_row else 0
    df = pd.read_csv(io.StringIO(csv_text), header=header_row)

    # ── Required columns ──────────────────────────────────────────────────────
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # ── Optional: Vol_mL (scenario detection) ─────────────────────────────────
    # Presence triggers variable-volume path; absence → constant-volume path.
    # Do NOT synthesize a default here — _compute() handles the vol=1 fallback
    # so that the scenario flag is unambiguous.

    # ── Optional: is_post_feed (feed-event marker) ────────────────────────────
    # Irrelevant for constant-volume / batch datasets; defaults to False.
    if FEED_COL not in df.columns:
        df[FEED_COL] = False
    elif df[FEED_COL].dtype == object:
        df[FEED_COL] = (
            df[FEED_COL].str.strip().str.upper()
            .map({"TRUE": True, "FALSE": False, "1": True, "0": False})
            .fillna(False)
        )
    df[FEED_COL] = df[FEED_COL].astype(bool)

    # ── Numeric coercion ──────────────────────────────────────────────────────
    numeric = ["t_hr","Rep","VCD","DCD","Viab_pct","rP_mg_L",
               "Glc_g_L","Lac_g_L","Gln_mM","Glu_mM"] + \
              ([VOL_COL] if VOL_COL in df.columns else []) + \
              [c for c in OPTIONAL_COLS if c in df.columns]
    for col in numeric:
        if df[col].dtype == object:
            df[col] = df[col].str.replace(",", ".", regex=False)
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # ── Cytometry forward-fill ────────────────────────────────────────────────
    cyto = [c for c in OPTIONAL_COLS if c in df.columns]
    if cyto:
        df = df.sort_values(["Clone","Rep","t_hr"]).reset_index(drop=True)
        df[cyto] = df.groupby(["Clone","Rep"])[cyto].transform(lambda s: s.ffill())

    return df.sort_values(["Clone","Rep","t_hr"]).reset_index(drop=True)


# ── Kinetics ───────────────────────────────────────────────────────────────────

def _compute(df: pd.DataFrame, exp_start: float, exp_end: float,
             use_volume: bool = True, fluor_set=None) -> pd.DataFrame:
    """
    Compute kinetic and metabolic parameters for every consecutive interval
    within each Clone × Rep group.

    Two independent calculation paths (PDF §2 and §3):
    ────────────────────────────────────────────────────
    • variable_volume  (Vol_mL present AND use_volume=True):  PDF §3
        - TC    = VCD × V  (total viable cells)
        - M_i   = C_i × V  (total mass / moles of metabolite i)
        - ΔITVC = (TC_1 + TC_2) / 2 × Δt
        - q_i   = ΔM_i / ΔITVC
        - Y     = ΔM_product / ΔM_substrate

    • constant_volume  (Vol_mL absent OR use_volume=False):   PDF §2
        - No volume data used at all
        - ΔIVCD = (VCD_1 + VCD_2) / 2 × Δt
        - q_i   = ΔC_i / ΔIVCD
        - Y     = ΔC_product / ΔC_substrate

    Critical exception (both scenarios):
    ─────────────────────────────────────
    μ ALWAYS uses VCD (concentration), never TC (total cells).
    Using TC for μ when volume decreases (sampling) yields spurious
    negative values that do not reflect true cellular growth rate.
    """
    df = df.copy()

    # ── Scenario detection ────────────────────────────────────────────────────
    # has_vol is True only when the column is present AND the user has opted in
    has_vol = (VOL_COL in df.columns) and bool(use_volume)

    out_cols = [MU, QGLC, QLAC, QGLC_PMOL, QLAC_PMOL, QP,
                QGLN_H, QGLN_D, QGLU_H, QGLU_D,
                Y_LG, Y_GQ, IVCD_INT, IVCD_CUM]
    if has_vol:
        out_cols += [ITVC_INT, ITVC_CUM]
    out_cols += [rc for _, rc, _ in _active_fluor(df, fluor_set)]
    for col in out_cols:
        df[col] = np.nan

    for (_, _rep), grp in df.groupby(["Clone","Rep"], sort=False):
        idx = grp.index.tolist()
        ivcd_cum = itvc_cum = 0.0

        for k in range(1, len(idx)):
            ip, ic = idx[k-1], idx[k]
            rp, rc = df.loc[ip], df.loc[ic]

            # Skip the pre-feed → post-feed dilution interval (not cellular activity)
            if (not rp[FEED_COL]) and rc[FEED_COL]:
                continue
            dt = rc["t_hr"] - rp["t_hr"]
            if dt <= 0:
                continue

            v1, v2 = rp["VCD"], rc["VCD"]

            # ── μ: ALWAYS VCD-based  (PDF §2.1 / §3.4) ───────────────────────
            # Using TC for μ produces spurious negative values whenever volume
            # drops (e.g. sampling removes medium), confusing dilution with death.
            if v1 > 0 and v2 > 0:
                df.at[ic, MU] = np.log(v2 / v1) / dt

            if has_vol:
                # ══════════════════════════════════════════════════════════════
                # VARIABLE VOLUME — mass-balance approach  (PDF §3)
                #
                #   TC    = VCD × V          (total viable cells in reactor)
                #   M_i   = C_i × V          (total mass/moles of metabolite i)
                #   ΔITVC = (TC_1 + TC_2)/2 × Δt
                #   q_i   = ΔM_i / ΔITVC
                # ══════════════════════════════════════════════════════════════
                vol1, vol2 = rp[VOL_COL], rc[VOL_COL]

                # TC = VCD [cells/mL] × V [mL]  →  [cells]
                TC1 = v1 * vol1
                TC2 = v2 * vol2

                # ΔITVC = (TC_1 + TC_2) / 2 × Δt  [cells·h]  (PDF eq. 3.3)
                delta_ITVC = 0.5 * (TC1 + TC2) * dt
                if delta_ITVC <= 0:
                    continue

                # M_i = C_i [g/L or mM] × V [mL] / 1000  →  [g] or [mmol]
                # Dividing by 1000 converts V from mL to L to get SI mass units.
                M_Glc1 = rp["Glc_g_L"] * vol1 / 1000   # [g]
                M_Glc2 = rc["Glc_g_L"] * vol2 / 1000
                M_Lac1 = rp["Lac_g_L"] * vol1 / 1000   # [g]
                M_Lac2 = rc["Lac_g_L"] * vol2 / 1000
                M_Gln1 = rp["Gln_mM"]  * vol1 / 1000   # [mmol]
                M_Gln2 = rc["Gln_mM"]  * vol2 / 1000
                M_Glu1 = rp["Glu_mM"]  * vol1 / 1000   # [mmol]
                M_Glu2 = rc["Glu_mM"]  * vol2 / 1000
                rP1 = rp["rP_mg_L"] if pd.notna(rp["rP_mg_L"]) else 0.0
                rP2 = rc["rP_mg_L"] if pd.notna(rc["rP_mg_L"]) else 0.0
                M_rP1 = rP1 * vol1 / 1000               # [mg]
                M_rP2 = rP2 * vol2 / 1000

                # q_i = ΔM_i / ΔITVC  with unit conversions to pg/cell/day
                # [g]    / [cells·h] × 10^12 pg/g    × 24 h/day  → pg/cell/day
                # [mg]   / [cells·h] × 10^9  pg/mg   × 24 h/day  → pg/cell/day
                # [mmol] / [cells·h] × 10^9  pmol/mmol           → pmol/cell/h
                qglc_pg = (M_Glc1 - M_Glc2) / delta_ITVC * PG_PER_G  * H_PER_DAY
                qlac_pg = (M_Lac2 - M_Lac1) / delta_ITVC * PG_PER_G  * H_PER_DAY
                qP_pg   = (M_rP2  - M_rP1 ) / delta_ITVC * 1e9       * H_PER_DAY
                qgh     = (M_Gln1 - M_Gln2) / delta_ITVC * MMOL_TO_PMOL
                qguh    = (M_Glu2 - M_Glu1) / delta_ITVC * MMOL_TO_PMOL

                # Y = ΔM_product / ΔM_substrate  (PDF §2.3 extended to masses)
                # The /1000 factors cancel in the ratio.
                glc_consumed = M_Glc1 - M_Glc2
                lac_produced = M_Lac2 - M_Lac1
                gln_consumed = M_Gln1 - M_Gln2
                glu_produced = M_Glu2 - M_Glu1

                # ITVC cumulative integral  (PDF eq. 3.3)
                itvc_cum += delta_ITVC
                df.at[ic, ITVC_INT] = delta_ITVC
                df.at[ic, ITVC_CUM] = itvc_cum

            else:
                # ══════════════════════════════════════════════════════════════
                # CONSTANT VOLUME — concentration-based approach  (PDF §2)
                #
                #   ΔIVCD = (VCD_1 + VCD_2)/2 × Δt
                #   q_i   = ΔC_i / ΔIVCD
                #
                # No volume data needed — concentration changes alone are used.
                # ══════════════════════════════════════════════════════════════

                # ΔIVCD = (VCD_1 + VCD_2) / 2 × Δt  [cells/mL · h]  (PDF eq. 2.2)
                delta_IVCD = 0.5 * (v1 + v2) * dt
                if delta_IVCD <= 0:
                    continue

                # Concentration changes — sign convention: positive = consumed/produced
                dGlc = rp["Glc_g_L"] - rc["Glc_g_L"]   # [g/L]   positive = consumed
                dLac = rc["Lac_g_L"] - rp["Lac_g_L"]   # [g/L]   positive = produced
                dGln = rp["Gln_mM"]  - rc["Gln_mM"]    # [mM]    positive = consumed
                dGlu = rc["Glu_mM"]  - rp["Glu_mM"]    # [mM]    positive = produced
                rP1  = rp["rP_mg_L"] if pd.notna(rp["rP_mg_L"]) else 0.0
                rP2  = rc["rP_mg_L"] if pd.notna(rc["rP_mg_L"]) else 0.0
                drP  = rP2 - rP1                         # [mg/L]  positive = produced

                # q_i = ΔC_i / ΔIVCD  with unit conversions to pg/cell/day
                # VCD is in cells/mL and concentrations are in /L, so the
                # mixed unit [g/L] / [cells/mL · h] carries a mL/L = 1e-3 factor.
                # [g/L]  / [cells/mL · h] × 1e-3 × 10^12 pg/g  × 24 h/day → pg/cell/day
                # [mg/L] / [cells/mL · h] × 1e-3 × 10^9  pg/mg × 24 h/day → pg/cell/day
                # [mM]   / [cells/mL · h] × 1e-3 × 10^9  pmol/mmol        → pmol/cell/h
                qglc_pg = dGlc / delta_IVCD * 1e-3 * PG_PER_G  * H_PER_DAY
                qlac_pg = dLac / delta_IVCD * 1e-3 * PG_PER_G  * H_PER_DAY
                qP_pg   = drP  / delta_IVCD * 1e-3 * 1e9       * H_PER_DAY
                qgh     = dGln / delta_IVCD * 1e-3 * MMOL_TO_PMOL
                qguh    = dGlu / delta_IVCD * 1e-3 * MMOL_TO_PMOL

                # Y = ΔC_product / ΔC_substrate  (PDF §2.3)
                glc_consumed = dGlc
                lac_produced = dLac
                gln_consumed = dGln
                glu_produced = dGlu

            # ── Write rates to DataFrame (both scenarios) ─────────────────────
            df.at[ic, QGLC]      = qglc_pg
            df.at[ic, QGLC_PMOL] = qglc_pg / MW_GLC
            df.at[ic, QLAC]      = qlac_pg
            df.at[ic, QLAC_PMOL] = qlac_pg / MW_LAC
            df.at[ic, QP]        = qP_pg
            df.at[ic, QGLN_H]    = qgh
            df.at[ic, QGLN_D]    = qgh * H_PER_DAY
            df.at[ic, QGLU_H]    = qguh
            df.at[ic, QGLU_D]    = qguh * H_PER_DAY

            # ── Yields (both scenarios) ───────────────────────────────────────
            if glc_consumed > 0:
                df.at[ic, Y_LG] = lac_produced / glc_consumed
            if gln_consumed > 0:
                df.at[ic, Y_GQ] = glu_produced / gln_consumed

            # ── IVCD: always computed  (PDF §2.2) ─────────────────────────────
            # IVCD_2 = IVCD_1 + (VCD_1 + VCD_2)/2 × Δt
            ivcd_int  = 0.5 * (v1 + v2) * dt
            ivcd_cum += ivcd_int
            df.at[ic, IVCD_INT] = ivcd_int
            df.at[ic, IVCD_CUM] = ivcd_cum

            # ── Fluorescence rates of change ──────────────────────────────────
            for mean_col, rate_col, _ in _active_fluor(df, fluor_set):
                f1, f2 = rp[mean_col], rc[mean_col]
                if pd.notna(f1) and pd.notna(f2):
                    df.at[ic, rate_col] = (f2 - f1) / dt

    df[PHASE] = np.where(
        (df["t_hr"] >= exp_start) & (df["t_hr"] <= exp_end),
        PHASE_EXP, PHASE_STAT
    )
    return df


def _summarise(df: pd.DataFrame, fluor_set=None) -> pd.DataFrame:
    exp = df[df[PHASE] == PHASE_EXP]
    rate_map = [(QGLC_PMOL,"qGlc_pmol_exp"),(QLAC_PMOL,"qLac_pmol_exp"),(QP,"qP_exp"),
                (QGLN_D,"qGln_exp"),(QGLU_D,"qGlu_exp"),
                (Y_LG,"Y_Lac_Glc_exp"),(Y_GQ,"Y_Glu_Gln_exp")]
    rate_map += [(rc, f"d{lbl}_dt_exp") for _, rc, lbl in _active_fluor(df, fluor_set)]

    records = []
    for (clone, rep), grp in exp.groupby(["Clone","Rep"], sort=False):
        row = {"Clone": clone, "Rep": rep}
        vcd_v = grp["VCD"].dropna()
        if len(vcd_v) >= 2:
            t0 = grp.loc[vcd_v.index[0],  "t_hr"]
            tf = grp.loc[vcd_v.index[-1], "t_hr"]
            dt = tf - t0
            if dt > 0 and vcd_v.iloc[0] > 0 and vcd_v.iloc[-1] > 0:
                row["mu_exp"] = np.log(vcd_v.iloc[-1] / vcd_v.iloc[0]) / dt
        for col, key in rate_map:
            if col in grp.columns:
                row[key] = grp[col].mean(skipna=True)
        records.append(row)
    return pd.DataFrame(records)


# ── Plotly chart data builders ──────────────────────────────────────────────────

def _ts_layout(ylabel, ylim=None, hovermode="closest"):
    """Shared layout for time-series charts (scatter & lines)."""
    layout = {
        "plot_bgcolor":  "white",
        "paper_bgcolor": "white",
        "font":   {"size": 12},
        "height": 380,
        "margin": {"t": 20, "r": 20, "b": 55, "l": 65},
        "xaxis": {
            "title":     "Culture time (h)",
            "range":     list(X_RANGE),
            "showgrid":  True,
            "gridcolor": "#e5e5e5",
            "zeroline":  False,
        },
        "yaxis": {
            "title":     ylabel,
            "showgrid":  True,
            "gridcolor": "#e5e5e5",
            "zeroline":  False,
        },
        "legend":    {"title": {"text": "Clone"}},
        "hovermode": hovermode,
    }
    if ylim is not None:
        layout["yaxis"]["range"] = list(ylim)
    return layout


def _scatter_data(df, clones, pal):
    """One scatter trace per clone, one chart per TS_SPEC."""
    out = {}
    for col, ylabel, fname, ylim in TS_SPECS:
        if col not in df.columns or df[col].isna().all():
            continue
        traces = []
        for c in clones:
            sub = df[df["Clone"] == c].dropna(subset=[col])
            if sub.empty:
                continue
            traces.append({
                "x":    _to_list(sub["t_hr"]),
                "y":    _to_list(sub[col]),
                "name": str(c),
                "mode": "markers",
                "type": "scatter",
                "marker": {"color": pal[c], "size": 8, "opacity": 0.8},
                "hovertemplate": (
                    f"<b>{c}</b><br>"
                    f"t: %{{x}} h<br>"
                    f"{ylabel}: %{{y:.4g}}"
                    "<extra></extra>"
                ),
            })
        if traces:
            out[fname] = {"traces": traces, "layout": _ts_layout(ylabel, ylim)}
    return out


def _lines_data(df, clones, pal):
    """Mean ± SD line per clone, pre-feed rows only."""
    out = {}
    pre = df[~df["is_post_feed"]]
    for col, ylabel, fname, ylim in TS_SPECS:
        if col not in df.columns or df[col].isna().all():
            continue
        s = (pre.groupby(["Clone","t_hr"])[col]
             .agg(mean="mean", std="std").reset_index())
        traces = []
        for c in clones:
            sub = s[s["Clone"] == c].sort_values("t_hr")
            if sub.empty:
                continue
            traces.append({
                "x":    _to_list(sub["t_hr"]),
                "y":    _to_list(sub["mean"]),
                "error_y": {
                    "type":      "data",
                    "array":     [_std0(v) for v in sub["std"].tolist()],
                    "visible":   True,
                    "thickness": 1.5,
                    "width":     4,
                },
                "name": str(c),
                "mode": "lines+markers",
                "type": "scatter",
                "line":   {"color": pal[c], "width": 1.8},
                "marker": {"color": pal[c], "size": 6},
                "hovertemplate": (
                    f"<b>{c}</b><br>"
                    f"t: %{{x}} h<br>"
                    f"Mean: %{{y:.4g}}"
                    "<extra></extra>"
                ),
            })
        if traces:
            out[fname] = {"traces": traces, "layout": _ts_layout(ylabel, ylim, "x unified")}
    return out


def _bars_data(summary, clones, pal):
    """Bar (mean ± SD) + replicate scatter overlay, one chart per BAR_SPEC."""
    out = {}
    for col, ylabel, fname in BAR_SPECS:
        if col not in summary.columns or summary[col].isna().all():
            continue
        agg = (summary.groupby("Clone")[col]
               .agg(mean="mean", std="std").reindex(clones))
        traces = []
        for c in clones:
            if c not in agg.index:
                continue
            mean_val = _v(float(agg.loc[c, "mean"]))
            std_val  = _std0(agg.loc[c, "std"])
            rep_vals = [_v(v) for v in
                        summary.loc[summary["Clone"] == c, col].dropna().tolist()]

            traces.append({
                "type": "bar",
                "x":    [str(c)],
                "y":    [mean_val],
                "error_y": {
                    "type":      "data",
                    "array":     [std_val],
                    "visible":   True,
                    "thickness": 1.5,
                    "width":     6,
                },
                "name":        str(c),
                "marker":      {"color": pal[c], "opacity": 0.85},
                "showlegend":  True,
                "legendgroup": str(c),
                "hovertemplate": (
                    f"<b>{c}</b><br>Mean: %{{y:.4g}}<extra></extra>"
                ),
            })
            if rep_vals:
                traces.append({
                    "type": "scatter",
                    "x":    [str(c)] * len(rep_vals),
                    "y":    rep_vals,
                    "mode": "markers",
                    "name": str(c),
                    "marker": {
                        "color": pal[c],
                        "size":  9,
                        "line":  {"color": "white", "width": 1.5},
                    },
                    "showlegend":  False,
                    "legendgroup": str(c),
                    "hovertemplate": (
                        f"<b>{c}</b><br>Replicate: %{{y:.4g}}<extra></extra>"
                    ),
                })
        if not traces:
            continue
        layout = {
            "plot_bgcolor":  "white",
            "paper_bgcolor": "white",
            "font":   {"size": 12},
            "height": 380,
            "margin": {"t": 20, "r": 20, "b": 55, "l": 65},
            "xaxis": {
                "title":         "Clone",
                "type":          "category",
                "categoryorder": "array",
                "categoryarray": [str(c) for c in clones],
                "showgrid":      False,
            },
            "yaxis": {
                "title":     ylabel,
                "showgrid":  True,
                "gridcolor": "#e5e5e5",
                "zeroline":  False,
            },
            "barmode": "group",
            "legend":  {"title": {"text": "Clone"}},
        }
        out[fname] = {"traces": traces, "layout": layout}
    return out


def _make_corr_data(df, xcol, ycol, xlabel, ylabel, clones, pal):
    """Scatter + regression line per clone × phase. Returns {traces, layout} or None."""
    if PHASE not in df.columns:
        return None
    sub = df[[xcol, ycol, "Clone", PHASE]].dropna()
    if sub.empty:
        return None

    traces = []
    for clone in clones:
        for phase, (symbol, dash, short) in _PHASE_STYLE.items():
            cd = sub[(sub["Clone"] == clone) & (sub[PHASE] == phase)]
            if cd.empty:
                continue
            grp_key = f"{clone}_{phase}"

            traces.append({
                "x":    _to_list(cd[xcol]),
                "y":    _to_list(cd[ycol]),
                "name": f"{clone} – {short}",
                "mode": "markers",
                "type": "scatter",
                "marker": {
                    "color":  pal[clone],
                    "symbol": symbol,
                    "size":   8,
                    "opacity": 0.8,
                },
                "legendgroup": grp_key,
                "showlegend":  True,
                "hovertemplate": (
                    f"<b>{clone} ({short})</b><br>"
                    f"{xlabel}: %{{x:.4g}}<br>"
                    f"{ylabel}: %{{y:.4g}}"
                    "<extra></extra>"
                ),
            })
            if len(cd) >= 3 and cd[xcol].nunique() > 1:
                sl, ic_r, r, _, _ = stats.linregress(cd[xcol], cd[ycol])
                x_fit = [float(cd[xcol].min()), float(cd[xcol].max())]
                y_fit = [float(sl * x + ic_r) for x in x_fit]
                traces.append({
                    "x":    x_fit,
                    "y":    y_fit,
                    "name": f"R²={r**2:.2f} ({clone}/{short})",
                    "mode": "lines",
                    "type": "scatter",
                    "line": {"color": pal[clone], "dash": dash, "width": 1.8},
                    "legendgroup": grp_key,
                    "showlegend":  True,
                    "hoverinfo":   "skip",
                })

    if not traces:
        return None

    layout = {
        "plot_bgcolor":  "white",
        "paper_bgcolor": "white",
        "font":   {"size": 12},
        "height": 400,
        "margin": {"t": 20, "r": 20, "b": 55, "l": 65},
        "xaxis": {
            "title":     xlabel,
            "showgrid":  True,
            "gridcolor": "#e5e5e5",
            "zeroline":  False,
        },
        "yaxis": {
            "title":     ylabel,
            "showgrid":  True,
            "gridcolor": "#e5e5e5",
            "zeroline":  False,
        },
        "legend": {
            "title": {"text": "Clone – Phase"},
            "font":  {"size": 10},
        },
    }
    return {"traces": traces, "layout": layout}


def _correlations_data(df, clones, pal):
    out = {}
    for xcol, ycol, xlabel, ylabel, fname in CORR_SPECS:
        if xcol not in df.columns or ycol not in df.columns:
            continue
        result = _make_corr_data(df, xcol, ycol, xlabel, ylabel, clones, pal)
        if result:
            out[fname] = result
    return out


# ── Analysis state (persists for custom correlations) ──────────────────────────
_state = {}


def make_custom_correlation(xcol, ycol):
    """
    Generate a single custom correlation plot using the last loaded dataset.
    Called from JavaScript after run_analysis has been executed.
    Returns a {traces, layout} dict or None if insufficient data.
    """
    if not _state:
        raise RuntimeError("No analysis loaded. Run run_analysis() first.")
    label_map = {col: label for label, col in CUSTOM_CORR_COLS}
    xlabel = label_map.get(xcol, xcol)
    ylabel = label_map.get(ycol, ycol)
    return _make_corr_data(
        _state["df"], xcol, ycol, xlabel, ylabel,
        _state["clones"], _state["pal"],
    )


# ── Public entry point ─────────────────────────────────────────────────────────

def run_analysis(csv_text, exp_phase_start=0.0, exp_phase_end=96.0,
                 skip_first_row=True, progress_cb=None, use_volume=True,
                 fluor_channels=""):
    """
    Main function called from JavaScript via Pyodide.

    Parameters
    ----------
    use_volume : bool
        When True (default) and Vol_mL is present in the CSV, use the
        variable-volume (mass-balance) scenario. When False, force the
        constant-volume (concentration-based) scenario regardless of
        whether Vol_mL is in the CSV.
    fluor_channels : str
        Comma-separated list of fluorescence channel labels to include,
        e.g. "GFP,TMRM". Empty string or None means all channels with
        data are included.

    Returns a dict with keys:
        info           – dataset metadata (includes 'scenario', 'active_fluor')
        processed_csv  – full kinetics CSV as string
        summary_csv    – exp-phase summary CSV as string
        avail_cols     – list of [label, colname] for custom correlation selectors
        plots          – {scatter, lines, bars, correlations}
                         each value is {fname: {traces, layout}} for Plotly.js
    """
    def cb(msg, pct):
        if progress_cb is not None:
            progress_cb(msg, pct)

    # Parse user-selected fluorescence channels ("" or None → include all)
    if fluor_channels:
        fluor_set = {lbl.strip() for lbl in str(fluor_channels).split(",") if lbl.strip()}
    else:
        fluor_set = None

    cb("Loading and cleaning data…", 5)
    df = _load(csv_text, bool(skip_first_row))

    cb("Computing kinetics…", 20)
    df_kin  = _compute(df, float(exp_phase_start), float(exp_phase_end),
                       use_volume=bool(use_volume), fluor_set=fluor_set)
    summary = _summarise(df_kin, fluor_set=fluor_set)

    clones = df_kin["Clone"].unique().tolist()
    pal    = _palette(clones)

    _state.clear()
    _state.update({"df": df_kin, "clones": clones, "pal": pal})

    avail_cols = [
        [label, col] for label, col in CUSTOM_CORR_COLS
        if col in df_kin.columns and not df_kin[col].isna().all()
    ]

    cb("Building scatter charts…", 35)
    scatter = _scatter_data(df_kin, clones, pal)

    cb("Building line charts…", 52)
    lines   = _lines_data(df_kin, clones, pal)

    cb("Building bar charts…", 68)
    bars    = _bars_data(summary, clones, pal)

    cb("Building correlation charts…", 84)
    corr    = _correlations_data(df_kin, clones, pal)

    cb("Done!", 100)

    scenario = SCENARIO_VAR if ((VOL_COL in df_kin.columns) and bool(use_volume)) else SCENARIO_CONST
    active_fluor = [lbl for _, _, lbl in _active_fluor(df_kin, fluor_set)]

    return {
        "info": {
            "n_clones":     int(df_kin["Clone"].nunique()),
            "n_reps":       int(df_kin["Rep"].nunique()),
            "n_timepoints": int(df_kin["t_hr"].nunique()),
            "n_rows":       int(len(df_kin)),
            "clones":       clones,
            "active_fluor": active_fluor,
            "scenario":     scenario,
        },
        "avail_cols":    avail_cols,
        "processed_csv": df_kin.to_csv(index=False),
        "summary_csv":   summary.to_csv(index=False),
        "plots": {
            "scatter":      scatter,
            "lines":        lines,
            "bars":         bars,
            "correlations": corr,
        },
    }
