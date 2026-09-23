#!/usr/bin/env python3
"""Calibrate effective fire-protection conductivity from a manufacturer table.

Manufacturers publish, for a given critical steel temperature, the required
insulation thickness for combinations of fire-resistance class (R30, R60, ...)
and section factor (Ap/V). This tool inverts such a table into a single
**effective thermal conductivity** lambda_p for use in the simplified
EN 1993-1-2 protected-steel model of ``steel_temperature.py`` - so the simple
model reproduces the manufacturer's certified performance.

The density and specific heat are taken from the datasheet (or a reasonable
assumption); only lambda_p is fitted.

Method
------
For a table cell (section factor A_p/V, R-class -> ISO 834 duration t_R,
thickness d_p) the protected-steel model reaches the critical temperature
theta_cr at t_R for exactly one lambda_p. Because published thicknesses are
rounded UP to the nearest available board, and low cells sit at the minimum
board (over-protected), the fit:

* de-quantizes each governing cell to the midpoint of its rounding interval
  (d_p - step/2) as an unbiased estimate of the true required thickness;
* excludes cells at the minimum-thickness floor (they only bound lambda_p);
* fits one lambda_p minimising the steel-temperature residual at those
  de-quantized thicknesses;
* validates by predicting each cell's thickness, rounding it back onto the
  board ladder, and comparing with the table (match rate + mm error).

The physics (ISO 834 curve, protected-steel solver, temperature-dependent
steel specific heat) is imported from ``steel_temperature.py`` so the
calibrated lambda_p is guaranteed consistent with the post-processor.

Usage
-----
    python calibrate_protection.py --table taulukko.xlsx --critical-temp 450 \\
        --rho 164 --cp 1030 --declared-lambda 0.039 --product "Stone wool"
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import steel_temperature as st

# Time step for the calibration integrations (s). EN 1993-1-2 permits <= 30 s
# for protected members; the solver sub-steps to at most this.
CAL_DT = 30.0
# Search brackets.
LAMBDA_BRACKET = (1e-3, 1.0)      # W/mK
THICKNESS_BRACKET = (1e-3, 0.5)   # m


# ============================================================
# TABLE READING
# ============================================================

def parse_rclass_minutes(label):
    """'R 30' / 'R30' / '30' / 'R 90 min' -> 30/90 (minutes) or None."""
    m = re.search(r"(\d+)", str(label))
    return int(m.group(1)) if m else None


def read_table(path):
    """Read the manufacturer table.

    First column = section factor Ap/V [1/m]; remaining columns = R-classes;
    cells = required thickness [mm]. Returns a tidy DataFrame with columns
    section_factor, r_minutes, thickness_mm (missing cells dropped).
    """
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"Table not found: {path}")
    if path.suffix.lower() in (".xlsx", ".xls"):
        raw = pd.read_excel(path)
    else:
        raw = pd.read_csv(path)

    sf_col = raw.columns[0]
    records = []
    for col in raw.columns[1:]:
        minutes = parse_rclass_minutes(col)
        if minutes is None:
            continue
        for _, row in raw.iterrows():
            sf = row[sf_col]
            th = row[col]
            if pd.isna(sf) or pd.isna(th):
                continue
            records.append({
                "section_factor": float(sf),
                "r_minutes": minutes,
                "thickness_mm": float(th),
            })
    df = pd.DataFrame.from_records(records)
    if df.empty:
        raise SystemExit("No usable (section factor, R-class, thickness) cells found.")
    return df


def thickness_ladder(df):
    """Sorted unique available board thicknesses [mm]."""
    return np.array(sorted(df["thickness_mm"].unique()))


def dequantize(dp_mm, ladder):
    """Unbiased true-thickness estimate: midpoint of the round-up interval."""
    below = ladder[ladder < dp_mm - 1e-9]
    prev = below.max() if below.size else 0.0
    return dp_mm - (dp_mm - prev) / 2.0


def ceil_to_ladder(d_mm, ladder):
    """Round a required thickness up to the next available board."""
    above = ladder[ladder >= d_mm - 1e-9]
    return float(above.min()) if above.size else float(ladder.max())


# ============================================================
# MODEL EVALUATIONS (reuse steel_temperature physics)
# ============================================================

def steel_temp_at(t_r_s, section_factor, dp_m, lam, rho_p, cp):
    """Protected-steel temperature at t_R under ISO 834 [deg C]."""
    time = np.arange(0.0, t_r_s + CAL_DT, CAL_DT)
    iso = st.iso834_temperature(time)
    prot = {"thickness": dp_m, "conductivity": lam, "rho": rho_p, "cp": cp}
    return float(st.solve_protected_steel(time, iso, section_factor, prot)[-1])


def _bisect(f, lo, hi, ftol=0.05, xtol=1e-4, itmax=40):
    """Bisect f (a temperature residual, deg C) to ftol, or bracket to xtol."""
    flo, fhi = f(lo), f(hi)
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if flo * fhi > 0:
        return None  # not bracketed
    mid = 0.5 * (lo + hi)
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) < ftol or (hi - lo) < xtol:
            return mid
        if flo * fmid < 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return mid


def invert_lambda(t_r_s, section_factor, dp_m, theta_cr, rho_p, cp):
    """Effective lambda_p that makes Ts(t_R) == theta_cr at thickness dp."""
    f = lambda lam: steel_temp_at(t_r_s, section_factor, dp_m, lam, rho_p, cp) - theta_cr
    return _bisect(f, *LAMBDA_BRACKET)


def invert_thickness(t_r_s, section_factor, lam, theta_cr, rho_p, cp):
    """Required thickness [m] so Ts(t_R) == theta_cr for conductivity lam."""
    f = lambda d: steel_temp_at(t_r_s, section_factor, d, lam, rho_p, cp) - theta_cr
    return _bisect(f, *THICKNESS_BRACKET)


def fit_lambda(gov_cells, theta_cr, rho_p, cp):
    """Single lambda_p minimising the Ts residual over governing cells."""
    def sse(lam):
        return sum(
            (steel_temp_at(c["r_minutes"] * 60.0, c["section_factor"],
                           c["d_target_mm"] / 1000.0, lam, rho_p, cp) - theta_cr) ** 2
            for c in gov_cells
        )
    # Golden-section minimisation on a unimodal residual.
    lo, hi = LAMBDA_BRACKET
    gr = (np.sqrt(5) - 1) / 2
    a, b = lo, hi
    c1 = b - gr * (b - a)
    c2 = a + gr * (b - a)
    f1, f2 = sse(c1), sse(c2)
    for _ in range(60):
        if abs(b - a) < 1e-5:
            break
        if f1 < f2:
            b, c2, f2 = c2, c1, f1
            c1 = b - gr * (b - a)
            f1 = sse(c1)
        else:
            a, c1, f1 = c1, c2, f2
            c2 = a + gr * (b - a)
            f2 = sse(c2)
    return 0.5 * (a + b)


# ============================================================
# MAIN
# ============================================================

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Calibrate effective protection conductivity from a "
                    "manufacturer thickness table (EN 1993-1-2).")
    p.add_argument("-t", "--table", required=True,
                   help="Manufacturer table (.xlsx/.csv): col1=Ap/V, cols=R-classes, cells=mm.")
    p.add_argument("--critical-temp", type=float, required=True,
                   help="Critical steel temperature the table is tabulated for [deg C].")
    p.add_argument("--rho", type=float, required=True,
                   help="Protection density rho_p [kg/m3].")
    p.add_argument("--cp", type=float, required=True,
                   help="Protection specific heat c_p [J/kgK].")
    p.add_argument("--declared-lambda", type=float, default=None,
                   help="Datasheet (ambient) conductivity, for comparison [W/mK].")
    p.add_argument("--product", default="protection",
                   help="Product name (used in output filenames and the report).")
    p.add_argument("-o", "--output-prefix", default=None,
                   help="Output prefix (default: derived from --product).")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    theta_cr, rho_p, cp = args.critical_temp, args.rho, args.cp
    prefix = args.output_prefix or re.sub(r"\W+", "_", args.product.strip().lower())

    # Use the EN-permissible 30 s protected-member step for a fast calibration.
    st.MAX_TIME_STEP = CAL_DT

    df = read_table(args.table)
    ladder = thickness_ladder(df)
    floor = float(ladder.min())
    print(f"Read {len(df)} cells; board ladder (mm): "
          f"{', '.join(f'{v:g}' for v in ladder)}", file=sys.stderr)
    print(f"Critical temp {theta_cr:g} C; rho_p={rho_p:g} kg/m3; cp={cp:g} J/kgK",
          file=sys.stderr)

    # Per-cell inverted lambda + de-quantized target thickness.
    rows = []
    gov_cells = []
    for r in df.to_dict("records"):
        dp = r["thickness_mm"]
        t_r = r["r_minutes"] * 60.0
        floored = dp <= floor + 1e-9
        # Per-cell lambda is only meaningful for governing (non-floored) cells.
        lam_cell = None
        if not floored:
            lam_cell = invert_lambda(t_r, r["section_factor"], dp / 1000.0,
                                     theta_cr, rho_p, cp)
            gov_cells.append({**r, "d_target_mm": dequantize(dp, ladder)})
        rows.append({
            "Section factor (1/m)": r["section_factor"],
            "R (min)": r["r_minutes"],
            "Table thickness (mm)": dp,
            "Governing": not floored,
            "Per-cell lambda (W/mK)": round(lam_cell, 4) if lam_cell else None,
        })

    if not gov_cells:
        raise SystemExit("All cells sit at the minimum board; nothing to fit.")

    lam_fit = fit_lambda(gov_cells, theta_cr, rho_p, cp)

    # Per-R-class fit: stone-wool conductivity rises with temperature, so one
    # lambda per fire rating fits far better than a single global value.
    per_class = {}
    for r in sorted({c["r_minutes"] for c in gov_cells}):
        subset = [c for c in gov_cells if c["r_minutes"] == r]
        per_class[r] = fit_lambda(subset, theta_cr, rho_p, cp)

    # Validation: predict each cell, round onto the ladder, compare to the table,
    # for both the single global lambda and the per-R-class lambda.
    n_match = n_match_c = n_total = 0
    for rec, r in zip(rows, df.to_dict("records")):
        t_r = r["r_minutes"] * 60.0
        d1 = invert_thickness(t_r, r["section_factor"], lam_fit, theta_cr, rho_p, cp)
        d1_mm = ceil_to_ladder(d1 * 1000.0 if d1 else float("nan"), ladder)
        rec["Model on ladder (mm)"] = d1_mm
        rec["Error (mm)"] = round(d1_mm - r["thickness_mm"], 1)
        rec["Match"] = abs(d1_mm - r["thickness_mm"]) < 1e-6

        lam_c = per_class.get(r["r_minutes"], lam_fit)
        d2 = invert_thickness(t_r, r["section_factor"], lam_c, theta_cr, rho_p, cp)
        d2_mm = ceil_to_ladder(d2 * 1000.0 if d2 else float("nan"), ladder)
        rec["Model per-class (mm)"] = d2_mm
        rec["Match per-class"] = abs(d2_mm - r["thickness_mm"]) < 1e-6

        n_total += 1
        n_match += int(rec["Match"])
        n_match_c += int(rec["Match per-class"])

    report = pd.DataFrame(rows)
    gov = report[report["Governing"]]
    gov_lams = gov["Per-cell lambda (W/mK)"].dropna()

    # --- outputs --------------------------------------------------------
    csv_path = f"{prefix}_calibration.csv"
    report.to_csv(csv_path, index=False)

    plot_path = f"{prefix}_calibration.png"
    _plot(report, lam_fit, args.declared_lambda, args.product, plot_path)

    # --- console report -------------------------------------------------
    print("\n" + "=" * 64)
    print(f"PROTECTION CALIBRATION - {args.product}")
    print("=" * 64)
    print(f"  Fitted effective lambda_p : {lam_fit:.4f} W/mK")
    if args.declared_lambda:
        print(f"  Declared (ambient) lambda : {args.declared_lambda:.4f} W/mK "
              f"(ratio {lam_fit / args.declared_lambda:.2f})")
    print(f"  Per-cell lambda spread    : {gov_lams.min():.4f} - {gov_lams.max():.4f} "
          f"(median {gov_lams.median():.4f}) over {len(gov_lams)} governing cells")
    print(f"  Ladder-match (single lam) : {n_match}/{n_total} cells "
          f"({100.0 * n_match / n_total:.0f}%)")
    print(f"  Ladder-match (per R-class): {n_match_c}/{n_total} cells "
          f"({100.0 * n_match_c / n_total:.0f}%)")
    err = report["Error (mm)"].abs()
    print(f"  Thickness error (single)  : mean {err.mean():.1f} mm, max {err.max():.0f} mm")
    print("-" * 64)
    print("  Effective lambda_p by fire rating (recommended - captures the")
    print("  temperature dependence of stone-wool conductivity):")
    for r in sorted(per_class):
        print(f"      R{r:<4d}  lambda_p = {per_class[r]:.4f} W/mK")
    print("=" * 64)
    print("\nRecommended: use the per-R-class lambda_p matching each member's target")
    print("rating. Example protection block (R90) for steel_config:")
    r_ex = sorted(per_class)[len(per_class) // 2]
    print("    protection:")
    print(f"      thickness: 0.020         # m - set per member from the table")
    print(f"      conductivity: {per_class[r_ex]:.4f}    # W/mK (calibrated for R{r_ex})")
    print(f"      rho: {rho_p:g}")
    print(f"      cp: {cp:g}")
    print(f"\nSingle-value fallback: conductivity = {lam_fit:.4f} W/mK", file=sys.stderr)
    print(f"Saved {csv_path} and {plot_path}", file=sys.stderr)


def _plot(report, lam_fit, declared, product, path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Per-cell lambda vs section factor, coloured by R-class.
    for r in sorted(report["R (min)"].unique()):
        sub = report[(report["R (min)"] == r) & report["Governing"]]
        if not sub.empty:
            ax1.plot(sub["Section factor (1/m)"], sub["Per-cell lambda (W/mK)"],
                     "o-", label=f"R{r}")
    ax1.axhline(lam_fit, color="black", lw=2, label=f"fitted {lam_fit:.4f}")
    if declared:
        ax1.axhline(declared, color="grey", ls="--", label=f"declared {declared:.3f}")
    ax1.set_xlabel("Section factor Ap/V (1/m)")
    ax1.set_ylabel("Effective lambda_p (W/mK)")
    ax1.set_title(f"{product}: per-cell effective conductivity")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    # Table vs model-on-ladder thickness.
    ax2.plot(report["Table thickness (mm)"], report["Model on ladder (mm)"],
             "o", alpha=0.6)
    lim = [0, max(report["Table thickness (mm)"].max(),
                  report["Model on ladder (mm)"].max()) + 10]
    ax2.plot(lim, lim, "k--", lw=1)
    ax2.set_xlim(lim)
    ax2.set_ylim(lim)
    ax2.set_xlabel("Manufacturer table thickness (mm)")
    ax2.set_ylabel("Calibrated model thickness on ladder (mm)")
    ax2.set_title("Reproduction of the table (fitted lambda_p)")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
