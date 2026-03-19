"""
Clonalyzer processing module – runs inside Pyodide (browser WebAssembly).

Returns all results as a Python dict that JavaScript converts to a JS object.
"""

import io, base64
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# ── Constants ──────────────────────────────────────────────────────────────────
H_PER_DAY    = 24.0
PG_PER_G     = 1e12
MMOL_TO_PMOL = 1e9
FIG_SIZE     = (9, 6)
DPI          = 100
X_RANGE      = (-10, 270)

# Molar masses for pmol conversions (g/mol)
MW_GLC = 180.16   # glucose
MW_LAC = 90.08    # lactate / lactic acid

PALETTE = ["#000000","#FF0066","#107F80","#F0E442",
           "#0072B2","#D55E00","#CC79A7","#999999"]

REQUIRED_COLS = ["t_hr","Clone","Rep","is_post_feed",
                 "VCD","DCD","Viab_pct","rP_mg_L",
                 "Glc_g_L","Lac_g_L","Gln_mM","Glu_mM","Vol_mL"]
OPTIONAL_COLS = ["GFP_mean","GFP_std","TMRM_mean","TMRM_std"]

PHASE_EXP  = "Exponential"
PHASE_STAT = "Stationary"
PHASE_CLR  = {PHASE_EXP: "#0072B2", PHASE_STAT: "#D55E00"}

# output column names
MU        = "mu_per_h"
QGLC      = "qGlc_pg_cell_day"
QLAC      = "qLac_pg_cell_day"
QGLC_PMOL = "qGlc_pmol_cell_day"
QLAC_PMOL = "qLac_pmol_cell_day"
QP       = "qP_pg_cell_day"
QGLN_H   = "qGln_pmol_cell_h"
QGLN_D   = "qGln_pmol_cell_day"
QGLU_H   = "qGlu_pmol_cell_h"
QGLU_D   = "qGlu_pmol_cell_day"
Y_LG     = "Y_Lac_per_Glc_g_per_g"
Y_GQ     = "Y_Glu_per_Gln_mol_per_mol"
IVCD_INT = "IVCD_interval_cells_h_per_mL"
IVCD_CUM = "IVCD_cum_cells_h_per_mL"
IVC_INT  = "IVC_interval_cells_h"
IVC_CUM  = "IVC_cum_cells_h"
DGFP     = "dGFP_dt"
DTMRM    = "dTMRM_dt"
PHASE    = "Fase_Cultivo"

# plot specs: (column, y-label, filename-stem)
TS_SPECS = [
    ("VCD",      "VCD (cells/mL)",          "01_VCD"),
    ("Viab_pct", "Viability (%)",            "02_Viability"),
    (MU,         "μ (1/h)",                  "03_Mu"),
    ("Glc_g_L",  "Glucose (g/L)",            "04_Glc"),
    ("Lac_g_L",  "Lactate (g/L)",            "05_Lac"),
    ("Gln_mM",   "Glutamine (mM)",           "06_Gln"),
    ("Glu_mM",   "Glutamate (mM)",           "07_Glu"),
    ("rP_mg_L",  "rP Titer (mg/L)",          "08_Product"),
    (QGLC_PMOL,  "qGlc (pmol/cell/day)",     "09_qGlc"),
    (QLAC_PMOL,  "qLac (pmol/cell/day)",     "10_qLac"),
    (QGLN_D,     "qGln (pmol/cell/day)",     "11_qGln"),
    (QGLU_D,     "qGlu (pmol/cell/day)",     "12_qGlu"),
    (QP,         "qP (pg/cell/day)",         "13_qP"),
    (Y_LG,       "Y Lac/Glc (g/g)",         "14_YLacGlc"),
    (Y_GQ,       "Y Glu/Gln (mol/mol)",     "15_YGluGln"),
    (IVCD_CUM,   "IVCD (cells·h/mL)",        "16_IVCD"),
    ("GFP_mean", "GFP intensity (A.U.)",     "17_GFP"),
    ("TMRM_mean","TMRM intensity (A.U.)",    "18_TMRM"),
    (DGFP,       "dGFP/dt (A.U./h)",         "19_dGFP_dt"),
    (DTMRM,      "dTMRM/dt (A.U./h)",        "20_dTMRM_dt"),
]

BAR_SPECS = [
    ("mu_exp",           "μ exp (1/h)",               "03_Mu_Exp"),
    ("qGlc_pmol_exp",    "qGlc exp (pmol/cell/day)",  "09_qGlc_Exp"),
    ("qLac_pmol_exp",    "qLac exp (pmol/cell/day)",  "10_qLac_Exp"),
    ("qP_exp",        "qP exp (pg/cell/day)",     "13_qP_Exp"),
    ("qGln_exp",      "qGln exp (pmol/cell/day)", "11_qGln_Exp"),
    ("qGlu_exp",      "qGlu exp (pmol/cell/day)", "12_qGlu_Exp"),
    ("Y_Lac_Glc_exp", "Y Lac/Glc exp (g/g)",     "14_YLacGlc_Exp"),
    ("Y_Glu_Gln_exp", "Y Glu/Gln exp (mol/mol)", "15_YGluGln_Exp"),
    ("dGFP_dt_exp",   "dGFP/dt exp (A.U./h)",    "19_dGFP_dt_Exp"),
    ("dTMRM_dt_exp",  "dTMRM/dt exp (A.U./h)",   "20_dTMRM_dt_Exp"),
]

CORR_SPECS = [
    (QP,        "GFP_mean",  "qP (pg/cell/day)",      "GFP (A.U.)",          "40_qP_vs_GFP"),
    (QP,        "TMRM_mean", "qP (pg/cell/day)",      "TMRM (A.U.)",         "41_qP_vs_TMRM"),
    ("GFP_mean","TMRM_mean", "GFP (A.U.)",            "TMRM (A.U.)",         "42_GFP_vs_TMRM"),
    (QGLC_PMOL, "TMRM_mean", "qGlc (pmol/cell/day)",  "TMRM (A.U.)",         "50_qGlc_vs_TMRM"),
    (QP,        QGLC_PMOL,   "qP (pg/cell/day)",      "qGlc (pmol/cell/day)","52_qP_vs_qGlc"),
    (QLAC_PMOL, QGLC_PMOL,   "qLac (pmol/cell/day)",  "qGlc (pmol/cell/day)","55_qLac_vs_qGlc"),
    (QP,        DGFP,        "qP (pg/cell/day)",      "dGFP/dt (A.U./h)",    "60_qP_vs_dGFPdt"),
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def _setup_mpl():
    plt.rcParams.update({
        "axes.facecolor":  "white",
        "figure.facecolor":"white",
        "axes.grid":       True,
        "grid.color":      "#dddddd",
        "grid.linewidth":  0.8,
        "font.size":       13,
        "axes.spines.top":   False,
        "axes.spines.right": False,
    })

def _palette(clones):
    return {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(clones)}

def _to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64

def _has_cyto(df):
    return all(c in df.columns for c in ("GFP_mean", "TMRM_mean"))


# ── Loader ─────────────────────────────────────────────────────────────────────

def _load(csv_text: str) -> pd.DataFrame:
    df = pd.read_csv(io.StringIO(csv_text), header=1)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    numeric = ["t_hr","Rep","VCD","DCD","Viab_pct","rP_mg_L",
               "Glc_g_L","Lac_g_L","Gln_mM","Glu_mM","Vol_mL"] + \
              [c for c in OPTIONAL_COLS if c in df.columns]
    for col in numeric:
        if df[col].dtype == object:
            df[col] = df[col].str.replace(",", ".", regex=False)
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if df["is_post_feed"].dtype == object:
        df["is_post_feed"] = (
            df["is_post_feed"].str.strip().str.upper()
            .map({"TRUE": True, "FALSE": False, "1": True, "0": False})
            .fillna(False)
        )
    df["is_post_feed"] = df["is_post_feed"].astype(bool)

    cyto = [c for c in OPTIONAL_COLS if c in df.columns]
    if cyto:
        df = df.sort_values(["Clone","Rep","t_hr"]).reset_index(drop=True)
        df[cyto] = df.groupby(["Clone","Rep"])[cyto].transform(lambda s: s.ffill())

    return df.sort_values(["Clone","Rep","t_hr"]).reset_index(drop=True)


# ── Kinetics ───────────────────────────────────────────────────────────────────

def _compute(df: pd.DataFrame, exp_start: float, exp_end: float) -> pd.DataFrame:
    df = df.copy()
    out_cols = [MU,QGLC,QLAC,QGLC_PMOL,QLAC_PMOL,QP,QGLN_H,QGLN_D,QGLU_H,QGLU_D,
                Y_LG,Y_GQ,IVCD_INT,IVCD_CUM,IVC_INT,IVC_CUM]
    if _has_cyto(df):
        out_cols += [DGFP, DTMRM]
    for col in out_cols:
        df[col] = np.nan

    for (_, rep), grp in df.groupby(["Clone","Rep"], sort=False):
        idx = grp.index.tolist()
        ivcd_cum = ivc_cum = 0.0

        for k in range(1, len(idx)):
            ip, ic = idx[k-1], idx[k]
            rp, rc  = df.loc[ip], df.loc[ic]

            if (not rp["is_post_feed"]) and rc["is_post_feed"]:
                continue
            dt = rc["t_hr"] - rp["t_hr"]
            if dt <= 0:
                continue

            v1, v2   = rp["VCD"],   rc["VCD"]
            vol1, vol2 = rp["Vol_mL"], rc["Vol_mL"]

            if v1 > 0 and v2 > 0:
                df.at[ic, MU] = np.log(v2/v1) / dt

            ntrap = 0.5 * (v1*1e3*vol1 + v2*1e3*vol2)
            if ntrap <= 0:
                continue

            dglc = rp["Glc_g_L"]*vol1 - rc["Glc_g_L"]*vol2
            qglc_pg = dglc / dt / ntrap * PG_PER_G * H_PER_DAY
            df.at[ic, QGLC]      = qglc_pg
            df.at[ic, QGLC_PMOL] = qglc_pg / MW_GLC

            dlac = rc["Lac_g_L"]*vol2 - rp["Lac_g_L"]*vol1
            qlac_pg = dlac / dt / ntrap * PG_PER_G * H_PER_DAY
            df.at[ic, QLAC]      = qlac_pg
            df.at[ic, QLAC_PMOL] = qlac_pg / MW_LAC

            p1 = rp["rP_mg_L"] if pd.notna(rp["rP_mg_L"]) else 0.0
            p2 = rc["rP_mg_L"] if pd.notna(rc["rP_mg_L"]) else 0.0
            df.at[ic, QP] = (p2*vol2 - p1*vol1) / dt / ntrap * 1e9 * H_PER_DAY

            dgln = rp["Gln_mM"]*vol1 - rc["Gln_mM"]*vol2
            qgh  = dgln / dt / ntrap * MMOL_TO_PMOL
            df.at[ic, QGLN_H]  = qgh
            df.at[ic, QGLN_D]  = qgh * H_PER_DAY

            dglu = rc["Glu_mM"]*vol2 - rp["Glu_mM"]*vol1
            qguh = dglu / dt / ntrap * MMOL_TO_PMOL
            df.at[ic, QGLU_H]  = qguh
            df.at[ic, QGLU_D]  = qguh * H_PER_DAY

            glc_c = rp["Glc_g_L"]*vol1 - rc["Glc_g_L"]*vol2
            lac_p = rc["Lac_g_L"]*vol2 - rp["Lac_g_L"]*vol1
            if glc_c > 0:
                df.at[ic, Y_LG] = lac_p / glc_c

            gln_c = rp["Gln_mM"]*vol1 - rc["Gln_mM"]*vol2
            glu_p = rc["Glu_mM"]*vol2 - rp["Glu_mM"]*vol1
            if gln_c > 0:
                df.at[ic, Y_GQ] = glu_p / gln_c

            ivcd_int  = 0.5*(v1+v2)*dt
            ivcd_cum += ivcd_int
            df.at[ic, IVCD_INT] = ivcd_int
            df.at[ic, IVCD_CUM] = ivcd_cum

            ivc_int  = ntrap * dt
            ivc_cum += ivc_int
            df.at[ic, IVC_INT] = ivc_int
            df.at[ic, IVC_CUM] = ivc_cum

            if _has_cyto(df):
                g1, g2 = rp["GFP_mean"],  rc["GFP_mean"]
                if pd.notna(g1) and pd.notna(g2):
                    df.at[ic, DGFP] = (g2-g1) / dt
                t1, t2 = rp["TMRM_mean"], rc["TMRM_mean"]
                if pd.notna(t1) and pd.notna(t2):
                    df.at[ic, DTMRM] = (t2-t1) / dt

    df[PHASE] = np.where(
        (df["t_hr"] >= exp_start) & (df["t_hr"] <= exp_end),
        PHASE_EXP, PHASE_STAT
    )
    return df


def _summarise(df: pd.DataFrame) -> pd.DataFrame:
    exp = df[df[PHASE] == PHASE_EXP]
    rate_map = [(QGLC_PMOL,"qGlc_pmol_exp"),(QLAC_PMOL,"qLac_pmol_exp"),(QP,"qP_exp"),
                (QGLN_D,"qGln_exp"),(QGLU_D,"qGlu_exp"),
                (Y_LG,"Y_Lac_Glc_exp"),(Y_GQ,"Y_Glu_Gln_exp")]
    if _has_cyto(df):
        rate_map += [(DGFP,"dGFP_dt_exp"),(DTMRM,"dTMRM_dt_exp")]

    records = []
    for (clone, rep), grp in exp.groupby(["Clone","Rep"], sort=False):
        row = {"Clone": clone, "Rep": rep}
        # μ_exp: first-to-last VCD within the exponential window
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


# ── Plot generators ────────────────────────────────────────────────────────────

def _scatter(df, clones, pal):
    out = {}
    for col, ylabel, fname in TS_SPECS:
        if col not in df.columns or df[col].isna().all():
            continue
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        for c in clones:
            sub = df[df["Clone"] == c]
            ax.scatter(sub["t_hr"], sub[col], color=pal[c], label=c, s=60, alpha=0.8)
        ax.set_xlabel("Culture time (h)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(X_RANGE)
        ax.legend(title="Clone", bbox_to_anchor=(1.02,1), loc="upper left")
        fig.tight_layout()
        out[fname] = _to_b64(fig)
    return out


def _lines(df, clones, pal):
    out = {}
    pre = df[~df["is_post_feed"]]
    for col, ylabel, fname in TS_SPECS:
        if col not in df.columns or df[col].isna().all():
            continue
        s = (pre.groupby(["Clone","t_hr"])[col]
             .agg(mean="mean", std="std").reset_index())
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        for c in clones:
            sub = s[s["Clone"] == c].sort_values("t_hr")
            if sub.empty:
                continue
            color = pal[c]
            ax.plot(sub["t_hr"], sub["mean"], marker="o", color=color, label=c)
            ax.fill_between(sub["t_hr"],
                            sub["mean"] - sub["std"].fillna(0),
                            sub["mean"] + sub["std"].fillna(0),
                            color=color, alpha=0.15)
        ax.set_xlabel("Culture time (h)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(X_RANGE)
        ax.legend(title="Clone", bbox_to_anchor=(1.02,1), loc="upper left")
        fig.tight_layout()
        out[fname] = _to_b64(fig)
    return out


def _bars(summary, clones, pal):
    out = {}
    for col, ylabel, fname in BAR_SPECS:
        if col not in summary.columns or summary[col].isna().all():
            continue
        agg = (summary.groupby("Clone")[col]
               .agg(mean="mean", std="std").reindex(clones).reset_index())
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        x = np.arange(len(clones))
        ax.bar(x, agg["mean"], yerr=agg["std"].fillna(0),
               color=[pal[c] for c in clones], capsize=5,
               width=0.6, alpha=0.85, error_kw={"elinewidth":1.5})
        for xi, c in enumerate(clones):
            vals = summary.loc[summary["Clone"] == c, col].dropna()
            ax.scatter(np.full(len(vals), xi), vals,
                       color=pal[c], edgecolors="white", s=70, zorder=5)
        ax.set_xticks(x)
        ax.set_xticklabels(clones, rotation=30, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Clone")
        fig.tight_layout()
        out[fname] = _to_b64(fig)
    return out


def _correlations(df, clones, pal):
    """
    One scatter + regression trace per clone × phase combination.
    Color  = clone (from palette)
    Marker = phase  (circle = Exponential, triangle = Stationary)
    """
    # phase → (marker, linestyle, short label)
    PHASE_STYLE = {
        PHASE_EXP:  ("o", "-",  "Exp"),
        PHASE_STAT: ("^", "--", "Stat"),
    }

    out = {}
    if PHASE not in df.columns:
        return out

    for xcol, ycol, xlabel, ylabel, fname in CORR_SPECS:
        if xcol not in df.columns or ycol not in df.columns:
            continue
        sub = df[[xcol, ycol, "Clone", PHASE]].dropna()
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=FIG_SIZE)

        for clone in clones:
            for phase, (marker, ls, short) in PHASE_STYLE.items():
                cd = sub[(sub["Clone"] == clone) & (sub[PHASE] == phase)]
                if cd.empty:
                    continue
                label = f"{clone} – {short}"
                ax.scatter(cd[xcol], cd[ycol],
                           color=pal[clone], marker=marker,
                           alpha=0.80, s=60, label=label, zorder=3)
                if len(cd) >= 3:
                    sl, ic, r, _, _ = stats.linregress(cd[xcol], cd[ycol])
                    xf = np.linspace(cd[xcol].min(), cd[xcol].max(), 100)
                    ax.plot(xf, sl*xf+ic,
                            color=pal[clone], linestyle=ls, linewidth=1.8,
                            label=f"R²={r**2:.2f}")

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left",
                  fontsize=9, title="Clone – Phase")
        fig.tight_layout()
        out[fname] = _to_b64(fig)
    return out


# ── Public entry point ─────────────────────────────────────────────────────────

def run_analysis(csv_text, exp_phase_start=0.0, exp_phase_end=96.0, progress_cb=None):
    """
    Main function called from JavaScript via Pyodide.

    Returns a dict with keys:
        info           – dataset metadata
        processed_csv  – full kinetics CSV as string
        summary_csv    – exp-phase summary CSV as string
        plots          – {scatter, lines, bars, correlations} → {fname: base64_png}
    """
    def cb(msg, pct):
        if progress_cb is not None:
            progress_cb(msg, pct)

    _setup_mpl()

    cb("Loading and cleaning data…", 5)
    df = _load(csv_text)

    cb("Computing kinetics…", 20)
    df_kin  = _compute(df, float(exp_phase_start), float(exp_phase_end))
    summary = _summarise(df_kin)

    clones = df_kin["Clone"].unique().tolist()
    pal    = _palette(clones)

    cb("Generating scatter plots…", 35)
    scatter = _scatter(df_kin, clones, pal)

    cb("Generating line plots…", 55)
    lines   = _lines(df_kin, clones, pal)

    cb("Generating bar charts…", 72)
    bars    = _bars(summary, clones, pal)

    cb("Generating correlation plots…", 88)
    corr    = _correlations(df_kin, clones, pal)

    cb("Done!", 100)

    return {
        "info": {
            "n_clones":     int(df_kin["Clone"].nunique()),
            "n_reps":       int(df_kin["Rep"].nunique()),
            "n_timepoints": int(df_kin["t_hr"].nunique()),
            "n_rows":       int(len(df_kin)),
            "clones":       clones,
            "has_cyto":     bool(_has_cyto(df_kin)),
        },
        "processed_csv": df_kin.to_csv(index=False),
        "summary_csv":   summary.to_csv(index=False),
        "plots": {
            "scatter":      scatter,
            "lines":        lines,
            "bars":         bars,
            "correlations": corr,
        },
    }
