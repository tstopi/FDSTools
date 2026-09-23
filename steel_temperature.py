#!/usr/bin/env python3
"""FDS AST -> steel temperature post-processor (EN 1993-1-2).

Reads adiabatic-surface-temperature (AST) gas devices from an FDS
``<CHID>_devc.csv`` file, groups them by structural member and longitudinal
location, and integrates the EN 1993-1-2 steel heating equations for each
location (protected or unprotected, configurable per member type).

For every member the hottest location is identified and assessed against a
per-member critical temperature, giving the maximum steel temperature, the
utilisation ratio and the interpolated time at which the critical temperature
is first exceeded.

Outputs (each prefixed with the scenario CHID, e.g. ``Case_A_...``)
-------
* ``<CHID>_steel_temperature_locations.csv`` - steel temperatures, all locations
* ``<CHID>_steel_fire_results.xlsx``         - workbook (summary/locations/peaks/config)
* ``<CHID>_member_plots/*.png``              - one temperature plot per member
* ``<CHID>_failure_map_tcrit_<view>.png`` / ``_cre_<view>.png`` - failure-
  location maps, one file per view, plus ``_tcrit_heatmap.png``,
  ``_tcrit_contourf.png`` and ``_tcrit_contour.png`` (plan heat map and
  contours); only with ``--fds``, see below

The CHID is taken from ``--chid`` / the config, or derived from the input
filename (``<CHID>_devc.csv``), so several scenarios can be post-processed in
the same directory without overwriting each other.

Device naming convention (``<Member>_<Face>_<Location>``)::

    AP_F1_001      YP_F4_012      D_01_F2_007      AP-2_F3_005

Members are matched to a configuration entry by longest family-prefix match:
``AP``, ``AP-1`` and ``AP-15`` all inherit the ``AP`` entry, and ``D_01``,
``D_02`` ... inherit ``D_``. Matching stops at a letter boundary, so ``AP``
does not swallow an unrelated ``APRON``.

For each member a time-equivalence of fire exposure is also reported using the
Cumulative Radiant Energy (CRE) method: the ISO 834 standard-fire duration whose
cumulative radiant energy equals that of the natural (CFD/AST) exposure.

Physics reference: EN 1993-1-2:2005, sections 3.4.1.2 (specific heat),
4.2.5.1 (unprotected members) and 4.2.5.2 (protected members).

Usage
-----
Run with the built-in defaults below::

    python steel_temperature.py

or drive everything from an external YAML/JSON config (see
``steel_config.example.yaml``)::

    python steel_temperature.py --config my_project.yaml
    python steel_temperature.py --config my_project.yaml --input CHID_devc.csv

Failure-location maps
---------------------
Given the FDS input file (``--fds model.fds``), the &DEVC lines - including
those in &CATF include files - are matched to the AST devices, and each member
location is placed at the centroid of its face devices. Two maps are written -
the time to reach the critical temperature (``_tcrit``) and the CRE equivalent
time (``_cre``) - each as separate plan, elevation and isometric files
(``_plan``, ``_elevation_xz``, ``_elevation_yz``, ``_isometric``). The plan and
elevations are drawn to true scale with axes covering only the data range; a
view that collapses to a line (e.g. the plan of a planar truss) is skipped.
In plan, the failure time is also shown as a heat map (``_tcrit_heatmap``),
filled contours (``_tcrit_contourf``) and labelled isochrones
(``_tcrit_contour``). All three use one field: each location takes the
earliest failure of any steel within the plan resolution (the spacing of
locations along the members), so stacked chords and diagonals collapse to the
governing one; square cells of ``--heatmap-cell`` metres (default 0.5) take
the earliest of their locations; empty cells are interpolated linearly
between them; and locations that never fail count as the end of the
simulation. Times are rounded down to whole
minutes and grouped into classes starting at the first failure, either every
``--failure-interval`` minutes or automatically into about
``--failure-classes`` (default 4) classes of a readable width. The CRE map
always uses automatic classes::

    python steel_temperature.py --input CHID_devc.csv --fds CHID.fds
    python steel_temperature.py --input CHID_devc.csv --fds CHID.fds --failure-interval 5
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import textwrap
import types
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: no display required
import matplotlib.patches
import matplotlib.patheffects
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.tri
import numpy as np
import pandas as pd

# ============================================================
# USER SETTINGS  (built-in defaults; override with --config)
# ============================================================

INPUT_CSV = "fds_devc.csv"

# FDS job identifier. Left as None it is derived from the input filename
# (``<CHID>_devc.csv`` -> ``<CHID>``) and prepended to every output name, so
# several scenarios can be post-processed side by side without overwriting.
CHID = None

# How the AST devices on the faces of one location are combined into a single
# thermal boundary condition. "max" is conservative and recommended.
GROUPING_METHOD = "max"  # "max" or "mean"

# FDS input file whose &DEVC lines give the device coordinates. When set (or
# given with --fds) the failure-location maps are produced; None disables them.
FDS_INPUT = None

# Failure-time classes for the maps: times are rounded down to whole minutes and
# grouped from the first failure in steps of FAILURE_INTERVAL minutes, or - when
# None - in a readable step chosen to give about FAILURE_CLASSES classes.
FAILURE_INTERVAL = None
FAILURE_CLASSES = 4

# Cell size (m) of the plan heat map and contours of failure time.
HEATMAP_CELL = 0.5

INITIAL_STEEL_TEMP = 20.0  # deg C

# Maximum integration time step (s). The (possibly coarse) FDS output interval
# is subdivided into steps no larger than this, with the AST linearly
# interpolated in between. EN 1993-1-2 recommends <= 30 s; smaller is safer.
MAX_TIME_STEP = 5.0

# ------------------------------------------------------------
# Fire / boundary conditions (unprotected members)
# ------------------------------------------------------------
# The AST already embeds the radiative + convective environment, so it is used
# directly as the driving temperature in the standard EN net-heat-flux
# expression (convection + radiation to a surface at the steel temperature).
ALPHA_C = 25.0          # convective coefficient, W/m2K (25 standard, 35 hydrocarbon)
EMISSIVITY = 0.7        # resultant emissivity eps_m * eps_f (0.7 carbon steel)
CONFIG_FACTOR = 1.0     # view / configuration factor Phi
SHADOW_FACTOR = 1.0     # correction factor k_sh for the shadow effect
SIGMA = 5.67e-8         # Stefan-Boltzmann constant, W/m2K4

# ------------------------------------------------------------
# Steel properties
# ------------------------------------------------------------
RHO_STEEL = 7850.0      # kg/m3
# Specific heat is temperature-dependent per EN 1993-1-2 (see steel_specific_heat).

# ============================================================
# MEMBER PROPERTIES  (longest-prefix match)
# ============================================================
# Each entry may set:
#   AmV           - section factor Am/V (unprotected) or Ap/V (protected), 1/m
#   critical_temp - critical steel temperature, deg C (None to skip assessment)
#   protected     - True to use the protected-steel formulation
#   protection    - dict with thickness/conductivity/rho/cp (protected only)

MEMBER_PROPERTIES = {
    "AP": {
        "AmV": 115.0,
        "critical_temp": 550.0,
        "protected": False,
    },
    "YP": {
        "AmV": 115.0,
        "critical_temp": 550.0,
        "protected": False,
    },
    "D_": {
        "AmV": 165.0,
        "critical_temp": 620.0,
        "protected": True,
        "protection": {
            "thickness": 0.025,     # dp, m
            "conductivity": 0.12,   # lambda_p, W/mK
            "rho": 300.0,           # rho_p, kg/m3
            "cp": 1700.0,           # c_p, J/kgK
        },
    },
    "__default__": {
        "AmV": 150.0,
        "critical_temp": None,
        "protected": False,
    },
}

# Fallback protection used when a member is flagged protected but supplies no
# "protection" dict of its own.
DEFAULT_PROTECTION = {
    "thickness": 0.025,
    "conductivity": 0.12,
    "rho": 300.0,
    "cp": 1700.0,
}

# ============================================================
# OUTPUT FILES  (base names; the CHID is prepended at run time
# unless an explicit path is given in the config's "output" block)
# ============================================================

LOCATION_OUTPUT = "steel_temperature_locations.csv"
EXCEL_OUTPUT = "steel_fire_results.xlsx"
PLOT_DIR = Path("member_plots")
# Stem of the failure maps; "_tcrit_<view>.png" and "_cre_<view>.png" are
# appended (views: plan, elevation_xz, elevation_yz, isometric).
FAILURE_MAP_OUTPUT = "failure_map"


# ============================================================
# EN 1993-1-2 MATERIAL MODEL
# ============================================================

def steel_specific_heat(temp_c):
    """Temperature-dependent specific heat of carbon steel, J/kgK.

    EN 1993-1-2:2005 eq. (3.2a)-(3.2d). Accepts a scalar or ndarray.
    """
    t = np.asarray(temp_c, dtype=float)
    ca = np.empty_like(t)

    r1 = t < 600.0
    r2 = (t >= 600.0) & (t < 735.0)
    r3 = (t >= 735.0) & (t < 900.0)
    r4 = t >= 900.0

    ca[r1] = (
        425.0
        + 0.773 * t[r1]
        - 1.69e-3 * t[r1] ** 2
        + 2.22e-6 * t[r1] ** 3
    )
    # Guard the singularities at 735 C and 731 C in the EN expressions.
    ca[r2] = 666.0 + 13002.0 / np.clip(738.0 - t[r2], 1e-6, None)
    ca[r3] = 545.0 + 17820.0 / np.clip(t[r3] - 731.0, 1e-6, None)
    ca[r4] = 650.0

    return float(ca) if ca.ndim == 0 else ca


# ============================================================
# CONFIGURATION HELPERS
# ============================================================

@dataclass
class MemberConfig:
    """Resolved configuration for one member."""

    AmV: float
    critical_temp: float | None
    protected: bool
    protection: dict | None = field(default=None)


def _prefix_matches(member, prefix):
    """True if *member* belongs to the *prefix* family.

    A prefix matches when the member equals it, or continues with a non-letter
    boundary. So ``"AP"`` matches ``AP``, ``AP-1``, ``AP_3`` and ``AP15`` but
    not ``APRON``; ``"D_"`` matches ``D_01``.
    """
    if not member.startswith(prefix):
        return False
    if len(member) == len(prefix):
        return True
    return not member[len(prefix)].isalpha()


def get_member_config(member):
    """Return the :class:`MemberConfig` for *member* by longest-prefix match."""
    best_prefix = None
    for prefix in MEMBER_PROPERTIES:
        if prefix == "__default__":
            continue
        if _prefix_matches(member, prefix):
            if best_prefix is None or len(prefix) > len(best_prefix):
                best_prefix = prefix

    props = MEMBER_PROPERTIES[best_prefix] if best_prefix else MEMBER_PROPERTIES["__default__"]

    protected = bool(props.get("protected", False))
    protection = None
    if protected:
        protection = dict(DEFAULT_PROTECTION)
        protection.update(props.get("protection", {}))

    return MemberConfig(
        AmV=float(props["AmV"]),
        critical_temp=props.get("critical_temp"),
        protected=protected,
        protection=protection,
    )


# ============================================================
# DEVICE PARSING / GROUPING
# ============================================================

_DEVICE_RE = re.compile(r"(?P<member>.+)_F\d+_(?P<location>\d+)$")


def parse_device_name(name):
    """Split ``<Member>_F<face>_<location>`` into (member, location) or None."""
    m = _DEVICE_RE.match(str(name).strip())
    if not m:
        return None
    return m.group("member"), m.group("location")


def build_location_groups(columns):
    """Map ``(member, location)`` -> list of device column names."""
    groups: dict[tuple[str, str], list[str]] = {}
    skipped = []
    for col in columns:
        parsed = parse_device_name(col)
        if parsed is None:
            skipped.append(col)
            continue
        groups.setdefault(parsed, []).append(col)
    if skipped:
        print(f"  Ignored {len(skipped)} column(s) not matching the AST naming "
              f"convention (e.g. {skipped[0]!r}).")
    return groups


def aggregate_ast(df, cols):
    """Combine the face AST columns of one location into a single series."""
    method = GROUPING_METHOD.lower()
    if method == "max":
        return df[cols].max(axis=1)
    if method == "mean":
        return df[cols].mean(axis=1)
    raise ValueError(f"Unknown GROUPING_METHOD: {GROUPING_METHOD!r} (use 'max' or 'mean')")


# ============================================================
# STEEL TEMPERATURE SOLVERS  (EN 1993-1-2)
# ============================================================

def _substep_times(t0, t1):
    """Return sub-step boundaries covering [t0, t1] with dt <= MAX_TIME_STEP."""
    span = t1 - t0
    if span <= MAX_TIME_STEP:
        return np.array([t0, t1])
    n = int(np.ceil(span / MAX_TIME_STEP))
    return np.linspace(t0, t1, n + 1)


def solve_unprotected_steel(time, ast, AmV):
    """Integrate the unprotected-member heating equation (EN 1993-1-2 4.2.5.1)."""
    time = np.asarray(time, dtype=float)
    ast = np.asarray(ast, dtype=float)
    Ts = np.empty(len(time))
    Ts[0] = INITIAL_STEEL_TEMP
    ts = Ts[0]

    for i in range(1, len(time)):
        edges = _substep_times(time[i - 1], time[i])
        for k in range(1, len(edges)):
            dt = edges[k] - edges[k - 1]
            # AST linearly interpolated at the start of the sub-step.
            frac = (edges[k - 1] - time[i - 1]) / (time[i] - time[i - 1] or 1.0)
            tg = ast[i - 1] + frac * (ast[i] - ast[i - 1])

            h_conv = ALPHA_C * (tg - ts)
            h_rad = (CONFIG_FACTOR * EMISSIVITY * SIGMA
                     * ((tg + 273.15) ** 4 - (ts + 273.15) ** 4))
            h_net = h_conv + h_rad

            ca = steel_specific_heat(ts)
            ts += SHADOW_FACTOR * AmV / (ca * RHO_STEEL) * h_net * dt
        Ts[i] = ts

    return Ts


def solve_protected_steel(time, ast, AmV, protection):
    """Integrate the protected-member heating equation (EN 1993-1-2 4.2.5.2)."""
    time = np.asarray(time, dtype=float)
    ast = np.asarray(ast, dtype=float)
    dp = protection["thickness"]
    lam = protection["conductivity"]
    rho_p = protection["rho"]
    cp = protection["cp"]

    Ts = np.empty(len(time))
    Ts[0] = INITIAL_STEEL_TEMP
    ts = Ts[0]

    for i in range(1, len(time)):
        edges = _substep_times(time[i - 1], time[i])
        for k in range(1, len(edges)):
            dt = edges[k] - edges[k - 1]
            span = time[i] - time[i - 1] or 1.0
            f0 = (edges[k - 1] - time[i - 1]) / span
            f1 = (edges[k] - time[i - 1]) / span
            tg0 = ast[i - 1] + f0 * (ast[i] - ast[i - 1])
            tg1 = ast[i - 1] + f1 * (ast[i] - ast[i - 1])
            dtg = tg1 - tg0  # gas-temperature increment over the sub-step

            ca = steel_specific_heat(ts)
            # phi: ratio of protection heat capacity to steel heat capacity.
            phi = (cp * rho_p / (ca * RHO_STEEL)) * dp * AmV

            d_ts = (
                (lam * AmV) / (dp * ca * RHO_STEEL)
                * (tg0 - ts) / (1.0 + phi / 3.0) * dt
                - (np.exp(phi / 10.0) - 1.0) * dtg
            )
            # EN 1993-1-2: no negative increment while the gas is heating.
            if d_ts < 0.0 and dtg > 0.0:
                d_ts = 0.0
            ts += d_ts
        Ts[i] = ts

    return Ts


# ============================================================
# CRITICAL-TEMPERATURE ASSESSMENT
# ============================================================

def find_critical_time(time, temperature, critical_temp):
    """First time the temperature reaches *critical_temp* (linear interp)."""
    if critical_temp is None:
        return None
    time = np.asarray(time, dtype=float)
    temperature = np.asarray(temperature, dtype=float)

    if temperature[0] >= critical_temp:
        return float(time[0])

    for i in range(1, len(time)):
        t1, t2 = temperature[i - 1], temperature[i]
        if t1 < critical_temp <= t2:
            if abs(t2 - t1) < 1e-12:
                return float(time[i])
            frac = (critical_temp - t1) / (t2 - t1)
            return float(time[i - 1] + frac * (time[i] - time[i - 1]))
    return None


# ============================================================
# TIME EQUIVALENCE - CUMULATIVE RADIANT ENERGY (CRE)
# ============================================================
# The severity of a fire exposure is characterised by the cumulative radiant
# energy delivered to a surface. Radiation dominates at fire temperatures, so
# the energy is proportional to the time integral of (T + 273.15)^4 measured
# above the ambient baseline; the emissivity and Stefan-Boltzmann constant are
# identical for the natural and standard exposures and therefore cancel.
#
# The CRE-equivalent time t_e is the ISO 834 standard-fire duration whose
# cumulative radiant energy equals that of the natural (CFD/AST) exposure:
#
#     integral_0^t_fire  [ (T_ast + 273.15)^4 - (T0 + 273.15)^4 ]_+ dt
#   = integral_0^t_e     [ (T_iso + 273.15)^4 - (T0 + 273.15)^4 ] dt
#
# The exposure temperature is the AST, so t_e depends only on the fire, not on
# the section factor or protection.

# Time-equivalence settings (overridable via the config "time_equivalence" block).
TIME_EQUIVALENCE = True     # compute the CRE-equivalent time
CRE_AMBIENT_TEMP = 20.0     # T0 baseline for the excess-radiant-energy integral
CRE_MAX_EQUIV_TIME = 6 * 3600.0  # cap on the ISO search, s


def iso834_temperature(t_seconds):
    """ISO 834 standard fire temperature, deg C (t in seconds)."""
    t = np.asarray(t_seconds, dtype=float)
    return 20.0 + 345.0 * np.log10(8.0 * (t / 60.0) + 1.0)


def _excess_radiant(temp_c, t0=None):
    """Excess radiant emissive term (T+273.15)^4 - (T0+273.15)^4, floored at 0."""
    if t0 is None:
        t0 = CRE_AMBIENT_TEMP
    base = (t0 + 273.15) ** 4
    val = (np.asarray(temp_c, dtype=float) + 273.15) ** 4 - base
    return np.clip(val, 0.0, None)


# np.trapz was renamed to np.trapezoid in NumPy 2.x (and later removed).
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def cumulative_radiant_energy(time, exposure_temp):
    """Total cumulative radiant energy of an exposure (trapezoidal integral)."""
    time = np.asarray(time, dtype=float)
    return float(_trapz(_excess_radiant(exposure_temp), time))


def equivalent_time_cre(time, exposure_temp):
    """CRE-equivalent ISO 834 exposure time, s (None if the exposure is nil)."""
    e_nat = cumulative_radiant_energy(time, exposure_temp)
    if e_nat <= 0.0:
        return 0.0

    # Cumulative radiant energy of the ISO 834 curve on a 1 s grid.
    dt = 1.0
    grid = np.arange(0.0, CRE_MAX_EQUIV_TIME + dt, dt)
    iso_excess = _excess_radiant(iso834_temperature(grid))
    e_iso = np.concatenate(([0.0], np.cumsum((iso_excess[1:] + iso_excess[:-1]) * dt / 2.0)))

    if e_nat >= e_iso[-1]:
        return None  # exceeds the search cap; exposure hotter/longer than cap

    i = int(np.searchsorted(e_iso, e_nat))
    e1, e2 = e_iso[i - 1], e_iso[i]
    if e2 - e1 < 1e-30:
        return float(grid[i])
    frac = (e_nat - e1) / (e2 - e1)
    return float(grid[i - 1] + frac * dt)


# ============================================================
# CSV READER (handles the two-row FDS devc.csv header)
# ============================================================

def read_devc_csv(path):
    """Read an FDS ``_devc.csv`` (units row + name row) or a plain CSV."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")

    with path.open("r") as fh:
        first = fh.readline()
        second = fh.readline()

    # FDS writes a units row first, then the device-ID row containing "Time".
    header_row = 1 if "time" in second.lower() and "time" not in first.lower() else 0
    df = pd.read_csv(path, header=header_row)
    df.columns = [str(c).strip() for c in df.columns]

    time_col = next((c for c in df.columns if c.upper() == "TIME"), None)
    if time_col is None:
        raise RuntimeError("Could not find a 'Time' column in the input CSV.")

    # Coerce everything to numeric; drop rows without a valid time.
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=[time_col]).reset_index(drop=True)
    return df, time_col


# ============================================================
# EXCEL WORKBOOK
# ============================================================

def write_excel(path, summary_df, location_df, peaks_df, config_df, map_df=None):
    """Write the results workbook with conditional formatting.

    *map_df*, when given, is written as a "Location Map" sheet after the
    member summary.
    """
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter

    green = PatternFill("solid", fgColor="C6EFCE")
    red = PatternFill("solid", fgColor="FFC7CE")
    header_fill = PatternFill("solid", fgColor="305496")
    header_font = Font(bold=True, color="FFFFFF")

    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        summary_df.to_excel(xl, sheet_name="Member Summary", index=False)
        if map_df is not None:
            map_df.to_excel(xl, sheet_name="Location Map", index=False)
        location_df.to_excel(xl, sheet_name="Location Temperatures", index=False)
        peaks_df.to_excel(xl, sheet_name="Member Peaks", index=False)
        config_df.to_excel(xl, sheet_name="Configuration", index=False)

        for ws in xl.book.worksheets:
            # Style the header row and auto-size columns.
            for cell in ws[1]:
                cell.fill = header_fill
                cell.font = header_font
            for col_cells in ws.columns:
                width = max((len(str(c.value)) for c in col_cells if c.value is not None),
                            default=8)
                ws.column_dimensions[get_column_letter(col_cells[0].column)].width = min(width + 3, 40)
            ws.freeze_panes = "A2"

        # Colour the Member Summary rows by critical-temperature exceedance.
        ws = xl.book["Member Summary"]
        cols = {c.value: c.column for c in ws[1]}
        util_c = cols.get("Utilisation")
        tmax_c = cols.get("Maximum Temperature (C)")
        tcrit_c = cols.get("Critical Temperature (C)")
        for row in range(2, ws.max_row + 1):
            tcrit = ws.cell(row=row, column=tcrit_c).value if tcrit_c else None
            tmax = ws.cell(row=row, column=tmax_c).value if tmax_c else None
            if tcrit in (None, "") or tmax in (None, ""):
                continue
            fill = red if float(tmax) >= float(tcrit) else green
            if util_c:
                ws.cell(row=row, column=util_c).fill = fill
            ws.cell(row=row, column=tmax_c).fill = fill


# ============================================================
# PLOTTING
# ============================================================

def plot_member(member, time, series, hottest_location, config):
    plt.figure(figsize=(10, 6))
    plt.plot(time, series, linewidth=2, color="#1f4e79", label=f"{member} steel temp")

    crit = config.critical_temp
    if crit is not None:
        plt.axhline(crit, color="red", linestyle="--", linewidth=1.5,
                    label=f"Critical = {crit:.0f} °C")
        t_crit = find_critical_time(time, np.asarray(series), crit)
        if t_crit is not None:
            plt.plot(t_crit, crit, "o", color="red", markersize=8,
                     label=f"Reached @ {t_crit:.1f} s")
            plt.annotate(f"{t_crit:.1f} s", xy=(t_crit, crit), xytext=(10, 10),
                         textcoords="offset points", color="red",
                         arrowprops={"arrowstyle": "->", "color": "red"})
        else:
            plt.text(0.98, 0.05, "Critical temperature not reached",
                     transform=plt.gca().transAxes, ha="right", va="bottom",
                     color="green",
                     bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "green"})

    tmax = float(np.max(series))
    i_max = int(np.argmax(series))
    plt.plot(time[i_max], tmax, "^", color="black", markersize=7,
             label=f"Max = {tmax:.0f} °C")

    plt.xlabel("Time (s)")
    plt.ylabel("Steel Temperature (°C)")
    plt.title(f"{member}   (hottest location: {hottest_location})")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    PLOT_DIR.mkdir(exist_ok=True)
    plt.savefig(PLOT_DIR / f"{member}.png", dpi=200)
    plt.close()


# ============================================================
# FDS INPUT - DEVICE LOCATIONS
# ============================================================
# Each namelist record starts with "&NAME" as the first non-blank text on a
# line and ends at the first "/" outside a quoted string, so a record may span
# several lines and its strings may contain slashes (e.g. file paths).

_NAMELIST_START = re.compile(r"(?m)^[ \t]*&([A-Za-z_]+)")
_NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"


def _iter_namelists(text):
    """Yield ``(NAME, body)`` for every namelist record in an FDS input."""
    pos = 0
    while True:
        m = _NAMELIST_START.search(text, pos)
        if not m:
            return
        quote = None
        j = m.end()
        while j < len(text):
            ch = text[j]
            if quote:
                if ch == quote:
                    quote = None
            elif ch in "'\"":
                quote = ch
            elif ch == "/":
                break
            j += 1
        yield m.group(1).upper(), text[m.end():j]
        pos = j + 1


def _string_param(body, key):
    """Value of a quoted parameter, e.g. ID='AP_F1_001' (None if absent)."""
    m = re.search(rf"(?<![A-Za-z0-9_]){key}\s*=\s*(['\"])(.*?)\1", body,
                  re.IGNORECASE | re.DOTALL)
    return m.group(2).strip() if m else None


def _numeric_param(body, key, count):
    """First *count* numbers of a parameter, e.g. XYZ=1,2,3 (None if absent)."""
    m = re.search(rf"(?<![A-Za-z0-9_]){key}\s*=\s*((?:{_NUMBER}[\s,]*)+)", body,
                  re.IGNORECASE)
    if not m:
        return None
    values = re.findall(_NUMBER, m.group(1))
    if len(values) < count:
        return None
    return np.array([float(v.replace("d", "e").replace("D", "e"))
                     for v in values[:count]])


def parse_devc_locations(fds_file, _seen=None):
    """Map device ID -> {"xyz", "orientation", "quantity", "source"}.

    Reads every &DEVC record of *fds_file*, following &CATF OTHER_FILES
    includes (relative to the including file). A device given by XB is placed
    at the centre of its box.
    """
    path = Path(fds_file)
    seen = set() if _seen is None else _seen
    resolved = path.resolve()
    if resolved in seen:
        return {}
    seen.add(resolved)
    if not path.exists():
        print(f"  WARNING: FDS file not found: {path}")
        return {}

    devices = {}
    text = path.read_text(errors="replace")
    for name, body in _iter_namelists(text):
        if name == "CATF":
            m = re.search(r"OTHER_FILES\s*=(.*)", body, re.IGNORECASE | re.DOTALL)
            for _, other in re.findall(r"(['\"])(.*?)\1", m.group(1) if m else ""):
                included = parse_devc_locations(path.parent / other.strip(), seen)
                for dev_id, info in included.items():
                    if dev_id in devices:
                        print(f"  WARNING: duplicate DEVC ID {dev_id!r} in "
                              f"{info['source']}; keeping the first.")
                        continue
                    devices[dev_id] = info
            continue
        if name != "DEVC":
            continue

        dev_id = _string_param(body, "ID")
        if not dev_id:
            continue
        xyz = _numeric_param(body, "XYZ", 3)
        if xyz is None:
            xb = _numeric_param(body, "XB", 6)
            if xb is None:
                continue
            xyz = np.array([(xb[0] + xb[1]) / 2, (xb[2] + xb[3]) / 2,
                            (xb[4] + xb[5]) / 2])
        if dev_id in devices:
            print(f"  WARNING: duplicate DEVC ID {dev_id!r} in {path.name}; "
                  "keeping the first.")
            continue
        devices[dev_id] = {
            "xyz": xyz,
            "orientation": _numeric_param(body, "ORIENTATION", 3),
            "quantity": _string_param(body, "QUANTITY"),
            "source": path.name,
        }
    return devices


def location_coordinates(groups, devices):
    """Centroid of each location's face devices.

    Returns ``({location_key: (xyz, n_devices)}, [device IDs without a &DEVC])``.
    """
    coords = {}
    missing = []
    for (member, location), devs in groups.items():
        points = [devices[d]["xyz"] for d in devs if d in devices]
        missing.extend(d for d in devs if d not in devices)
        if points:
            coords[f"{member}_{location}"] = (np.mean(points, axis=0), len(points))
    return coords, missing


# ============================================================
# FAILURE-TIME CLASSES
# ============================================================
# Times are rounded DOWN to whole minutes (conservative, no sub-minute values)
# and grouped into classes that start at the first failure t0:
#
#     class(t) = t0 + floor((t - t0) / step) * step        [minutes]
#
# The step is FAILURE_INTERVAL when given, otherwise the smallest readable step
# that yields at most FAILURE_CLASSES classes.

READABLE_STEPS_MIN = (1, 2, 3, 5, 10, 15, 20, 30, 60)


def floor_minutes(t_seconds):
    """Whole minutes, rounded down; None for a missing time."""
    if t_seconds is None or pd.isna(t_seconds):
        return None
    return int(np.floor(t_seconds / 60.0 + 1e-9))


def choose_interval(t_first, t_last, n_classes):
    """Smallest readable step (min) giving at most *n_classes* classes."""
    need = max(1, int(np.ceil((t_last - t_first + 1) / max(1, int(n_classes)))))
    for step in READABLE_STEPS_MIN:
        if step >= need:
            return step
    return int(np.ceil(need / 60.0)) * 60


def failure_classes(minutes, interval=None, n_classes=4):
    """Group whole-minute times into classes.

    *minutes* maps a key to whole minutes (or None). Returns None when no key
    has a time, else a dict with ``t0`` (first time), ``step``, ``assign``
    (key -> class start or None) and ``starts`` (non-empty class starts).
    """
    values = [m for m in minutes.values() if m is not None]
    if not values:
        return None
    t0, t_last = min(values), max(values)
    step = choose_interval(t0, t_last, n_classes) if not interval else int(interval)
    step = max(1, step)
    assign = {k: (None if m is None else t0 + ((m - t0) // step) * step)
              for k, m in minutes.items()}
    starts = sorted({s for s in assign.values() if s is not None})
    return {"t0": t0, "step": step, "assign": assign, "starts": starts}


def class_label(start, step):
    return f"{start} min" if step == 1 else f"{start}–{start + step - 1} min"


# ============================================================
# FAILURE-LOCATION MAPS
# ============================================================

# One-hue ordinal ramp (light -> dark), validated for up to five classes on the
# light surface; more classes are interpolated and also get distinct markers.
_ORDINAL_RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
                 "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
_CLASS_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]
_SURFACE = "#fcfcfb"
_INK = "#0b0b0b"
_INK_2 = "#52514e"
_MUTED = "#c3c2b7"
_MEMBER_LINE = "#dcdbd5"
_MEMBER_LINE_DARK = "#b9b8b2"   # members on a plain background (contour lines)


def _class_styles(n):
    """(colour, marker) per class, ordered light -> dark."""
    if n <= 5:
        idx = [6] if n == 1 else np.round(np.linspace(0, 9, n)).astype(int)
        colours = [_ORDINAL_RAMP[i] for i in idx]
        markers = ["o"] * n
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("ord", _ORDINAL_RAMP)
        colours = [matplotlib.colors.to_hex(cmap(v)) for v in np.linspace(0, 1, n)]
        markers = [_CLASS_MARKERS[i % len(_CLASS_MARKERS)] for i in range(n)]
    return list(zip(colours, markers))


# View layout (inches): axes are sized so both axes share one scale (true
# proportions) and the limits cover only the data range plus a pad that keeps
# edge markers whole.
_MAX_AXES_W = 12.0
_MAX_AXES_H = 9.0
_MARKER_PAD_IN = 0.14
_VIEWS_2D = [("plan", "x", "y", "Plan (x–y)"),
             ("elevation_xz", "x", "z", "Elevation (x–z)"),
             ("elevation_yz", "y", "z", "Elevation (y–z)")]


def _row_major(items, ncol):
    """Reorder *items* so a column-filling legend reads row by row."""
    nrow = int(np.ceil(len(items) / ncol))
    grid = [items[r * ncol:(r + 1) * ncol] for r in range(nrow)]
    return [row[c] for c in range(ncol) for row in grid if c < len(row)]


def _save_trimmed(fig, path, top_in, bottom_in, pad_in=0.2):
    """Save *fig*, dropping empty rows between the title and legend bands.

    A 3D box never fills its axes, so the blank margin above and below it is
    removed from the rendered image (the title and legend bands are kept).
    """
    from PIL import Image  # bundled with matplotlib

    dpi = fig.dpi
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
    top = int(round(top_in * dpi))
    bottom = img.shape[0] - int(round(bottom_in * dpi))
    band = img[top:bottom]
    background = np.array(matplotlib.colors.to_rgb(_SURFACE)) * 255
    used = np.where((np.abs(band.astype(int) - background).max(axis=2) > 8).any(axis=1))[0]
    if used.size:
        pad = int(round(pad_in * dpi))
        band = band[max(0, used[0] - pad):min(len(band), used[-1] + 1 + pad)]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(np.concatenate([img[:top], band, img[bottom:]])).save(
        path, dpi=(dpi, dpi))


def _legend_layout(labels, fig_w):
    """Columns and height (inches) of a legend wrapped to the figure width."""
    entry_w = max(len(s) for s in labels) * 0.075 + 0.55
    ncol = max(1, min(len(labels), int((fig_w - 0.4) // entry_w)))
    return ncol, 0.28 * int(np.ceil(len(labels) / ncol)) + 0.15


def _title_and_legend(fig, fig_w, fig_h, heading, subtitle, handles, labels,
                      legend_top_in):
    """Title block at the top and a row-major legend below the axes."""
    fig.text(0.35 / fig_w, 1 - 0.2 / fig_h, heading,
             color=_INK, fontsize=13, ha="left", va="top")
    fig.text(0.35 / fig_w, 1 - 0.55 / fig_h, subtitle, color=_INK_2,
             fontsize=9, ha="left", va="top")
    ncol, _ = _legend_layout(labels, fig_w)
    fig.legend(_row_major(handles, ncol), _row_major(labels, ncol),
               loc="upper center", bbox_to_anchor=(0.5, legend_top_in / fig_h),
               ncol=ncol, frameon=False, fontsize=9, labelcolor=_INK)


def plot_failure_map(stem, points, classes, *, title, early_is_severe,
                     highlight, highlight_label, unassigned_label):
    """Write one PNG per view (plan, elevations, isometric); return the paths.

    *points* has columns key, member, location, x, y, z. *classes* is the
    result of :func:`failure_classes`. When *early_is_severe* the earliest
    class is drawn darkest (failure time); otherwise the latest (CRE).
    *highlight* is a list of keys ringed and listed as *highlight_label*.
    Files are named ``<stem>_plan.png``, ``<stem>_elevation_xz.png``,
    ``<stem>_elevation_yz.png`` and ``<stem>_isometric.png``; a view that
    collapses to a line (an axis spanning < 2 % of the model, e.g. the plan
    of a planar truss) is skipped.
    """
    starts = classes["starts"]
    step = classes["step"]
    styles = _class_styles(len(starts))
    if early_is_severe:
        styles = styles[::-1]
    style_of = dict(zip(starts, styles))
    # Most severe class drawn last so it sits on top.
    draw_order = starts[::-1] if early_is_severe else starts

    pts = points.copy()
    pts["cls"] = pts["key"].map(classes["assign"])
    unassigned = pts[pts["cls"].isna()]
    ring = pts[pts["key"].isin(highlight)]

    lo = {c: float(pts[c].min()) for c in ("x", "y", "z")}
    spans = {c: float(pts[c].max()) - lo[c] for c in ("x", "y", "z")}
    biggest = max(spans.values()) or 1.0
    flat = {c for c, v in spans.items() if v < 0.02 * biggest}

    def draw(ax, cols):
        # Members as faint lines through their locations, in location order.
        for _, grp in pts.groupby("member"):
            grp = grp.sort_values("location")
            ax.plot(*(grp[c] for c in cols), color=_MEMBER_LINE, lw=1.0, zorder=1)
        if len(unassigned):
            ax.scatter(*(unassigned[c] for c in cols), s=10, color=_MUTED,
                       linewidths=0, zorder=2)
        for z, start in enumerate(draw_order, start=3):
            sel = pts[pts["cls"] == start]
            colour, marker = style_of[start]
            ax.scatter(*(sel[c] for c in cols), s=42, color=colour, marker=marker,
                       edgecolors=_SURFACE, linewidths=0.8, zorder=z)
        # Rings keep full contrast in 3D (no depth shading), like the legend.
        no_fade = {"depthshade": False} if len(cols) == 3 else {}
        ax.scatter(*(ring[c] for c in cols), s=170, facecolors="none",
                   edgecolors=_INK, linewidths=1.4, zorder=len(starts) + 4, **no_fade)

    # Legend entries: classes from most to least severe, then neutral entries.
    handles, labels = [], []
    for start in (starts if early_is_severe else starts[::-1]):
        colour, marker = style_of[start]
        n = int((pts["cls"] == start).sum())
        handles.append(plt.Line2D([], [], ls="", marker=marker, ms=8, color=colour,
                                  markeredgecolor=_SURFACE))
        labels.append(f"{class_label(start, step)}  ({n})")
    if len(unassigned):
        handles.append(plt.Line2D([], [], ls="", marker="o", ms=5, color=_MUTED))
        labels.append(f"{unassigned_label}  ({len(unassigned)})")
    handles.append(plt.Line2D([], [], ls="", marker="o", ms=11, mfc="none",
                              mec=_INK, mew=1.4))
    labels.append(highlight_label)
    handles.append(plt.Line2D([], [], color=_MEMBER_LINE, lw=1.5))
    labels.append("Member")
    subtitle = (f"Times rounded down to whole minutes; classes of {step} min "
                f"starting at {classes['t0']} min.")

    def finish(fig, fig_w, fig_h, view_label, legend_top_in):
        _title_and_legend(fig, fig_w, fig_h, f"{title} – {view_label}",
                          subtitle, handles, labels, legend_top_in)

    written = []
    for suffix, a, b, view_label in _VIEWS_2D:
        if a in flat or b in flat:
            continue
        # One scale (in/m) for both axes, limits = data range + marker pad.
        scale = min((_MAX_AXES_W - 2 * _MARKER_PAD_IN) / max(spans[a], 1e-9),
                    (_MAX_AXES_H - 2 * _MARKER_PAD_IN) / max(spans[b], 1e-9))
        ax_w = spans[a] * scale + 2 * _MARKER_PAD_IN
        ax_h = spans[b] * scale + 2 * _MARKER_PAD_IN
        pad_a = pad_b = _MARKER_PAD_IN / scale

        fig_w = max(ax_w + 1.3, 8.0)
        _, legend_h = _legend_layout(labels, fig_w)
        top, xlab = 0.95, 0.6
        fig_h = top + ax_h + xlab + legend_h
        fig = plt.figure(figsize=(fig_w, fig_h), facecolor=_SURFACE)
        left = max(0.9, (fig_w - ax_w) / 2)
        ax = fig.add_axes([left / fig_w, (xlab + legend_h) / fig_h,
                           ax_w / fig_w, ax_h / fig_h])
        draw(ax, (a, b))
        ax.set_xlim(lo[a] - pad_a, lo[a] + spans[a] + pad_a)
        ax.set_ylim(lo[b] - pad_b, lo[b] + spans[b] + pad_b)
        ax.set_xlabel(f"{a} (m)", color=_INK_2)
        ax.set_ylabel(f"{b} (m)", color=_INK_2)
        ax.grid(True, color="#e6e5e0", lw=0.6)
        ax.set_facecolor(_SURFACE)
        ax.tick_params(colors=_INK_2, labelsize=8)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(
            nbins=max(2, int(ax_h / 0.45)), steps=[1, 2, 2.5, 5, 10]))
        for spine in ax.spines.values():
            spine.set_color("#d4d3cd")
        finish(fig, fig_w, fig_h, view_label, legend_h)

        path = f"{stem}_{suffix}.png"
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, facecolor=_SURFACE)
        plt.close(fig)
        written.append(path)

    if not flat:
        fig_w, fig_h = 11.0, 9.0
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=200, facecolor=_SURFACE)
        _, legend_h = _legend_layout(labels, fig_w)
        top_band, bottom_band = 0.9, legend_h + 0.1
        ax3d = fig.add_axes([0.02, bottom_band / fig_h, 0.96,
                             1 - (top_band + bottom_band) / fig_h], projection="3d")
        draw(ax3d, ("x", "y", "z"))
        for c in ("x", "y", "z"):
            getattr(ax3d, f"set_{c}lim")(lo[c], lo[c] + spans[c])
            getattr(ax3d, f"set_{c}label")(f"{c} (m)", color=_INK_2)
        # Keep thin directions at least a quarter of the longest so the view
        # stays readable (the plan and elevations carry the true proportions).
        ax3d.set_box_aspect([max(spans[c], 0.25 * biggest) for c in ("x", "y", "z")],
                           zoom=1.15)
        for axis in (ax3d.xaxis, ax3d.yaxis, ax3d.zaxis):
            axis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax3d.view_init(elev=25, azim=-60)
        ax3d.set_facecolor(_SURFACE)
        ax3d.tick_params(colors=_INK_2, labelsize=7)
        finish(fig, fig_w, fig_h, "Isometric", legend_h)
        path = f"{stem}_isometric.png"
        _save_trimmed(fig, path, top_band, bottom_band)
        plt.close(fig)
        written.append(path)

    return written


def plan_resolution(map_df):
    """Typical spacing (m) between neighbouring locations along a member."""
    gaps = []
    for _, grp in map_df.groupby("Member"):
        xyz = grp.sort_values("Location")[["X (m)", "Y (m)", "Z (m)"]].to_numpy(float)
        if len(xyz) > 1:
            gaps.extend(np.linalg.norm(np.diff(xyz, axis=0), axis=1))
    return float(np.median(gaps)) if gaps else 0.0


def plan_envelope(x, y, values, h):
    """Earliest value of any location within plan distance *h* of each location.

    Chords and diagonals stacked above one another interleave in plan; taking
    the governing (earliest) value within the plan resolution keeps the field
    from zig-zagging between them, whatever the cell size.
    """
    if h <= 0:
        return values.copy()
    out = np.empty_like(values)
    for s in range(0, len(x), 2000):                     # chunked to bound memory
        d2 = (x[s:s + 2000, None] - x) ** 2 + (y[s:s + 2000, None] - y) ** 2
        out[s:s + 2000] = np.where(d2 <= (h * 1.001) ** 2, values, np.inf).min(axis=1)
    return out


def plan_failure_field(map_df, classes, *, t_end_min, cell):
    """Failure-time field in plan, shared by the heat map and contour plots.

    Each location first takes the earliest failure of any steel within the plan
    resolution (the spacing of locations along the members), so stacked chords
    and diagonals collapse to the governing one; a location that never fails
    counts as the end of the simulation, *t_end_min*. Each square cell (side
    *cell* m) containing steel then takes the earliest of its locations. Empty cells are interpolated linearly between the
    cells with steel (Delaunay triangulation of their centres), without
    extrapolation. Returns None when the plan collapses to a line (e.g. a
    single planar truss), else a namespace with the cell grid (xe, ye), the
    triangulation (tri, or None if too few cells) and its values (tri_values),
    the rounded-down field and its class index per cell, and the classes.
    """
    x = map_df["X (m)"].to_numpy(float)
    y = map_df["Y (m)"].to_numpy(float)
    minutes = map_df["Failure Time (min)"].to_numpy(float)   # NaN = not reached
    span_x, span_y = float(np.ptp(x)), float(np.ptp(y))
    if min(span_x, span_y) < 0.02 * max(span_x, span_y, 1e-9):
        return None

    # Classes from the first failure to the simulation end (last one may be
    # partial). The field spans more time than the failures themselves, so an
    # automatic class width is widened to the next readable step until at most
    # max(FAILURE_CLASSES, 5) classes remain (five is the most the colour ramp
    # separates); an explicit FAILURE_INTERVAL is always honoured.
    t0, step = classes["t0"], classes["step"]
    last_minute = int(np.ceil(t_end_min)) - 1
    if not FAILURE_INTERVAL:
        cap = max(int(FAILURE_CLASSES), 5)
        while (last_minute - t0) // step + 1 > cap:
            larger = [s_ for s_ in READABLE_STEPS_MIN if s_ > step]
            step = larger[0] if larger else step + 60
    edges = [float(t0 + k * step) for k in range(int(np.ceil((t_end_min - t0) / step)))]
    edges = [e for e in edges if e < t_end_min] + [float(t_end_min)]
    n_cls = len(edges) - 1
    labels = [class_label(int(a), step) if b - a == step else
              (f"{int(a)} min" if last_minute <= a else f"{int(a)}–{last_minute} min")
              for a, b in zip(edges[:-1], edges[1:])]

    # Earliest failure per cell (not reached -> simulation end).
    nx = max(1, int(np.ceil(span_x / cell)))
    ny = max(1, int(np.ceil(span_y / cell)))
    xe = x.min() + cell * np.arange(nx + 1)
    ye = y.min() + cell * np.arange(ny + 1)
    ix = np.clip(((x - x.min()) // cell).astype(int), 0, nx - 1)
    iy = np.clip(((y - y.min()) // cell).astype(int), 0, ny - 1)
    resolution = plan_resolution(map_df)
    governing = plan_envelope(x, y, np.where(np.isnan(minutes), t_end_min, minutes),
                              resolution)
    measured = np.full((ny, nx), np.inf)
    np.minimum.at(measured, (iy, ix), governing)
    has = np.isfinite(measured)
    measured[~has] = np.nan
    # Centroid of each cell's locations: the contours are triangulated there,
    # so they reach the members themselves rather than stopping at cell centres.
    sums = np.zeros((3, ny, nx))
    for k, v in enumerate((x, y, np.ones_like(x))):
        np.add.at(sums[k], (iy, ix), v)
    px_c, py_c = sums[0][has] / sums[2][has], sums[1][has] / sums[2][has]
    try:
        contour_tri = matplotlib.tri.Triangulation(px_c, py_c)
    except (RuntimeError, ValueError):
        contour_tri = None

    # Linear interpolation between the cells with steel.
    xc, yc = np.meshgrid((xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2)
    field = measured.copy()
    tri = None
    try:
        tri = matplotlib.tri.Triangulation(xc[has], yc[has])
        interp = matplotlib.tri.LinearTriInterpolator(tri, measured[has])
        field = np.ma.filled(interp(xc, yc), np.nan)
        field[has] = measured[has]
    except (RuntimeError, ValueError):
        tri = None
        print("  Too few non-collinear cells to interpolate - the plan heat map "
              "shows cells with steel only, and no contours are drawn.")
    # Round down to whole minutes. The tolerance keeps an interpolated
    # 59.99999 (roundoff between values of 60) out of the class below; the
    # contours use the same shifted values, so all plots agree.
    field = np.floor(field + 1e-6)
    cls_idx = np.searchsorted(edges, field, side="right") - 1.0
    cls_idx = np.where(np.isnan(field), np.nan, np.clip(cls_idx, 0, n_cls))

    return types.SimpleNamespace(
        xe=xe, ye=ye, cell=cell, resolution=resolution, tri=contour_tri,
        tri_values=measured[has] + 1e-6,
        field=field, cls_idx=cls_idx, t0=t0, step=step, edges=edges,
        labels=labels, n_cls=n_cls, t_end_min=t_end_min,
        colours=[c for c, _ in _class_styles(n_cls)][::-1] + [_MUTED])


def _plan_figure(field, map_df, *, heading, subtitle, handles, labels,
                 member_colour, extent=None):
    """True-scale plan figure; returns (fig, ax).

    The axes cover *extent* (x0, x1, y0, y1), by default the field's cell grid.
    Members and first-failure rings are drawn above the data (zorder 2-3), so
    the caller adds the field itself at zorder 1.
    """
    x0, x1, y0, y1 = extent or (field.xe[0], field.xe[-1], field.ye[0], field.ye[-1])
    gx, gy = x1 - x0, y1 - y0
    scale = min(_MAX_AXES_W / gx, _MAX_AXES_H / gy)
    ax_w, ax_h = gx * scale, gy * scale
    fig_w = max(ax_w + 1.3, 8.0)
    _, legend_h = _legend_layout(labels, fig_w)
    # Wrap a long subtitle to the figure width (about 0.068 in per character
    # at 9 pt) and make room for the extra lines.
    lines = textwrap.wrap(subtitle, width=max(40, int((fig_w - 0.6) / 0.068)))
    subtitle = "\n".join(lines)
    top, xlab = 0.95 + 0.16 * (len(lines) - 1), 0.6
    fig_h = top + ax_h + xlab + legend_h
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=_SURFACE)
    left = max(0.9, (fig_w - ax_w) / 2)
    ax = fig.add_axes([left / fig_w, (xlab + legend_h) / fig_h, ax_w / fig_w, ax_h / fig_h])
    for _, grp in map_df.groupby("Member"):
        grp = grp.sort_values("Location")
        ax.plot(grp["X (m)"], grp["Y (m)"], color=member_colour, lw=1.0, zorder=2)
    first = map_df[map_df["Failure Time (min)"] == field.t0]
    ax.scatter(first["X (m)"], first["Y (m)"], s=170, facecolors="none",
               edgecolors=_INK, linewidths=1.4, zorder=3, clip_on=False)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("x (m)", color=_INK_2)
    ax.set_ylabel("y (m)", color=_INK_2)
    ax.set_facecolor(_SURFACE)
    ax.tick_params(colors=_INK_2, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#d4d3cd")
    _title_and_legend(fig, fig_w, fig_h, heading, subtitle, handles, labels, legend_h)
    return fig, ax


def _class_legend(field):
    """Legend entries for the filled plots: classes, not reached, members, rings."""
    handles = [matplotlib.patches.Patch(color=c) for c in field.colours]
    labels = list(field.labels) + [f"Not reached in {field.t_end_min:g} min"]
    outline = [matplotlib.patheffects.withStroke(linewidth=3, foreground=_MUTED)]
    handles.append(plt.Line2D([], [], color="white", lw=1.5, path_effects=outline))
    labels.append("Member (plan)")
    handles.append(plt.Line2D([], [], ls="", marker="o", ms=11, mfc="none",
                              mec=_INK, mew=1.4))
    labels.append(f"First failure ({field.t0} min)")
    return handles, labels


def _save(fig, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, facecolor=_SURFACE)
    plt.close(fig)
    return path


def plot_failure_heatmap(path, map_df, field, *, title):
    """Plan heat map: each cell coloured by its (interpolated) failure class."""
    handles, labels = _class_legend(field)
    fig, ax = _plan_figure(
        field, map_df, heading=f"{title} – Heat map (plan)",
        subtitle=(f"Earliest failure within {field.resolution:g} m in plan, per "
                  f"{field.cell:g} m cell; linear between cells; rounded down to whole "
                  f"minutes; classes of {field.step} min."),
        handles=handles, labels=labels, member_colour="white")
    cmap = matplotlib.colors.ListedColormap(field.colours)
    norm = matplotlib.colors.BoundaryNorm(np.arange(-0.5, field.n_cls + 1.5), cmap.N)
    ax.pcolormesh(field.xe, field.ye, np.ma.masked_invalid(field.cls_idx),
                  cmap=cmap, norm=norm, zorder=1)
    return _save(fig, path)


def plot_failure_contours(stem, map_df, field, *, title):
    """Filled contours and labelled isochrones of the failure-time field.

    Both are drawn on the triangulation of the cells with steel (the field
    the heat map interpolates), at the class boundaries. Contouring the
    continuous field at whole minutes equals classing the rounded-down times,
    since floor(t) >= n exactly when t >= n. Returns the written paths.
    """
    if field.tri is None:
        return []
    written = []
    vals = field.tri_values
    top = max(float(vals.max()), field.edges[-1]) + 1.0
    # The contoured field spans the triangulated cell centroids: fit the axes.
    tri = field.tri
    extent = (tri.x.min(), tri.x.max(), tri.y.min(), tri.y.max())

    handles, labels = _class_legend(field)
    fig, ax = _plan_figure(
        field, map_df, heading=f"{title} – Filled contours (plan)",
        subtitle=(f"Bands of the earliest failure within {field.resolution:g} m in plan, "
                  f"linear between {field.cell:g} m cells; rounded down to whole minutes; "
                  f"classes of {field.step} min."),
        handles=handles, labels=labels, member_colour="white", extent=extent)
    ax.tricontourf(field.tri, vals, levels=field.edges + [top],
                   colors=field.colours, zorder=1)
    written.append(_save(fig, f"{stem}_contourf.png"))

    handles = [plt.Line2D([], [], color=_INK, lw=1.0),
               plt.Line2D([], [], color=_MEMBER_LINE_DARK, lw=1.5),
               plt.Line2D([], [], ls="", marker="o", ms=11, mfc="none", mec=_INK, mew=1.4)]
    labels = [f"Isochrone (labelled; {field.t_end_min:g} min = not reached)",
              "Member (plan)", f"First failure ({field.t0} min)"]
    fig, ax = _plan_figure(
        field, map_df, heading=f"{title} – Contours (plan)",
        subtitle=(f"Isochrones every {field.step} min of the earliest failure within "
                  f"{field.resolution:g} m in plan, linear between {field.cell:g} m cells; "
                  "rounded down to whole minutes."),
        handles=handles, labels=labels, member_colour=_MEMBER_LINE_DARK,
        extent=extent)
    lines = ax.tricontour(field.tri, vals, levels=field.edges[1:], colors=_INK,
                          linewidths=0.9, zorder=2.5)
    ax.clabel(lines, fmt=lambda v: f"{int(round(v))} min", fontsize=8, inline=True)
    written.append(_save(fig, f"{stem}_contour.png"))
    return written


# ============================================================
# EXTERNAL CONFIGURATION (YAML / JSON)
# ============================================================
# Everything under USER SETTINGS / MEMBER PROPERTIES above acts as the built-in
# default. An optional --config file (YAML or JSON) overrides any of it, so the
# same script can be reused across projects without editing the source. See
# steel_config.example.yaml for the full schema.

# Maps a config key -> the module-level global it overrides.
_SCALAR_KEYS = {
    "input_csv": "INPUT_CSV",
    "chid": "CHID",
    "grouping_method": "GROUPING_METHOD",
    "initial_steel_temp": "INITIAL_STEEL_TEMP",
    "max_time_step": "MAX_TIME_STEP",
    "fds_input": "FDS_INPUT",
}
_FIRE_KEYS = {
    "alpha_c": "ALPHA_C",
    "emissivity": "EMISSIVITY",
    "config_factor": "CONFIG_FACTOR",
    "shadow_factor": "SHADOW_FACTOR",
    "sigma": "SIGMA",
}
_OUTPUT_KEYS = {
    "locations_csv": "LOCATION_OUTPUT",
    "excel": "EXCEL_OUTPUT",
    "plot_dir": "PLOT_DIR",
    "failure_map": "FAILURE_MAP_OUTPUT",
}

# Output globals whose path the user set explicitly (config/CLI) and which must
# therefore NOT receive the automatic CHID prefix.
_OUTPUT_OVERRIDDEN = set()


def derive_chid(input_csv):
    """Return the CHID from an FDS input filename (``<CHID>_devc.csv``)."""
    name = Path(input_csv).name
    if name.lower().endswith("_devc.csv"):
        return name[: -len("_devc.csv")]
    return Path(name).stem


def prefix_with_chid(path, chid):
    """Prepend ``<chid>_`` to the file/dir name, keeping any parent directory."""
    p = Path(path)
    return p.parent / f"{chid}_{p.name}"


def load_config(path):
    """Load a YAML or JSON config file into a dict (auto-detected by suffix)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    text = path.read_text()
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError(
                "PyYAML is required for YAML config files (pip install pyyaml), "
                "or use a .json config instead."
            ) from exc
        data = yaml.safe_load(text)
    elif path.suffix.lower() == ".json":
        data = json.loads(text)
    else:
        raise ValueError(f"Unsupported config extension: {path.suffix!r} "
                         "(use .yaml, .yml or .json)")
    if not isinstance(data, dict):
        raise ValueError("Config file must contain a mapping at the top level.")
    return data


def apply_config(cfg):
    """Override the module-level settings from a loaded config dict."""
    g = globals()

    for key, name in _SCALAR_KEYS.items():
        if key in cfg:
            g[name] = cfg[key]

    for key, name in _FIRE_KEYS.items():
        if key in cfg.get("fire", {}):
            g[name] = cfg["fire"][key]

    if "rho" in cfg.get("steel", {}):
        g["RHO_STEEL"] = cfg["steel"]["rho"]

    te = cfg.get("time_equivalence", {})
    if "enabled" in te:
        g["TIME_EQUIVALENCE"] = bool(te["enabled"])
    if "ambient_temp" in te:
        g["CRE_AMBIENT_TEMP"] = te["ambient_temp"]
    if "max_equiv_time" in te:
        g["CRE_MAX_EQUIV_TIME"] = te["max_equiv_time"]

    fm = cfg.get("failure_map", {})
    if "interval_min" in fm:
        g["FAILURE_INTERVAL"] = fm["interval_min"]
    if "n_classes" in fm:
        g["FAILURE_CLASSES"] = fm["n_classes"]
    if "heatmap_cell" in fm:
        g["HEATMAP_CELL"] = float(fm["heatmap_cell"])

    for key, name in _OUTPUT_KEYS.items():
        if key in cfg.get("output", {}):
            value = cfg["output"][key]
            g[name] = Path(value) if name == "PLOT_DIR" else value
            _OUTPUT_OVERRIDDEN.add(name)

    if "default_protection" in cfg:
        g["DEFAULT_PROTECTION"] = dict(cfg["default_protection"])

    # A `members` block, if present, replaces MEMBER_PROPERTIES wholesale so the
    # Configuration audit sheet reflects exactly what was supplied.
    if "members" in cfg:
        members = dict(cfg["members"])
        members.setdefault("__default__", MEMBER_PROPERTIES["__default__"])
        g["MEMBER_PROPERTIES"] = members


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="FDS AST -> steel temperature post-processor (EN 1993-1-2).")
    parser.add_argument(
        "-c", "--config", metavar="FILE",
        help="YAML or JSON configuration file (overrides the built-in defaults).")
    parser.add_argument(
        "-i", "--input", metavar="CSV",
        help="FDS <CHID>_devc.csv file (overrides input_csv from config/defaults).")
    parser.add_argument(
        "--chid", metavar="CHID",
        help="Scenario identifier prepended to output names "
             "(default: derived from the input filename).")
    parser.add_argument(
        "--fds", metavar="FDS",
        help="FDS input file; its &DEVC lines locate the devices and enable "
             "the failure-location maps.")
    parser.add_argument(
        "--failure-interval", metavar="MIN", type=int,
        help="Failure-time class width in whole minutes (default: automatic).")
    parser.add_argument(
        "--failure-classes", metavar="N", type=int,
        help="Target number of failure-time classes when the interval is "
             "automatic (default 4).")
    parser.add_argument(
        "--heatmap-cell", metavar="M", type=float,
        help="Cell size of the plan heat map and contours in metres "
             "(default 0.5).")
    return parser.parse_args(argv)


# ============================================================
# FAILURE-MAP ASSEMBLY
# ============================================================

def _whole_number_setting(value, name, minimum=1):
    """Coerce a setting to a whole number >= minimum, rounding down (None passes)."""
    if value is None:
        return None
    whole = max(minimum, int(np.floor(float(value))))
    if whole != value:
        print(f"  NOTE: {name}={value} used as {whole} (whole number, >= {minimum}).")
    return whole


def locate_devices(fds_file, groups, csv_devices):
    """Parse the FDS input and return the location centroids (with warnings)."""
    if not Path(fds_file).exists():
        raise SystemExit(f"FDS input file not found: {fds_file}")
    devices = parse_devc_locations(fds_file)
    coords, missing = location_coordinates(groups, devices)
    print(f"FDS:     {fds_file} - {len(devices)} &DEVC records, "
          f"{len(coords)}/{len(groups)} locations placed.")
    if missing:
        print(f"  WARNING: {len(missing)} AST device(s) in the CSV have no &DEVC "
              f"line (e.g. {missing[0]!r}); their locations use the other faces "
              "or are left off the map.")
    unused = [d for d in devices
              if parse_device_name(d) is not None and d not in csv_devices]
    if unused:
        print(f"  NOTE: {len(unused)} AST-style &DEVC ID(s) are not in the CSV "
              f"(e.g. {unused[0]!r}).")
    return coords


def build_location_map(time, member_locations, location_results, location_ast,
                       coords):
    """One row per placed location: coordinates, peak, failure and CRE times."""
    rows = []
    for member in sorted(member_locations):
        cfg = get_member_config(member)
        for key in member_locations[member]:
            if key not in coords:
                continue
            xyz, n_dev = coords[key]
            steel = location_results[key].to_numpy(dtype=float)
            t_crit = find_critical_time(time, steel, cfg.critical_temp)
            row = {
                "Location": key,
                "Member": member,
                "X (m)": round(float(xyz[0]), 4),
                "Y (m)": round(float(xyz[1]), 4),
                "Z (m)": round(float(xyz[2]), 4),
                "Devices": n_dev,
                "Maximum Temperature (C)": round(float(steel.max()), 1),
                "Critical Temperature (C)": cfg.critical_temp,
                "Critical Time (s)": round(t_crit, 1) if t_crit is not None else None,
                "Failure Time (min)": floor_minutes(t_crit),
            }
            if TIME_EQUIVALENCE:
                t_eq = equivalent_time_cre(time, location_ast[key])
                row["CRE Equivalent Time (min)"] = floor_minutes(t_eq)
            rows.append(row)
    return pd.DataFrame(rows)


def classify_location_map(map_df):
    """Add class columns to *map_df*; return {"tcrit": classes, "cre": classes}.

    FAILURE_INTERVAL applies to the failure times; the CRE equivalent times
    always use automatic classes (about FAILURE_CLASSES of them).
    """
    result = {}
    columns = [("tcrit", "Failure Time (min)", "Failure Class", FAILURE_INTERVAL)]
    if "CRE Equivalent Time (min)" in map_df:
        columns.append(("cre", "CRE Equivalent Time (min)", "CRE Class", None))
    for name, source, target, interval in columns:
        minutes = {k: (None if pd.isna(v) else int(v))
                   for k, v in zip(map_df["Location"], map_df[source])}
        classes = failure_classes(minutes, interval, FAILURE_CLASSES)
        result[name] = classes
        if classes is None:
            map_df[target] = None
            continue
        map_df[target] = [
            None if classes["assign"][k] is None
            else class_label(classes["assign"][k], classes["step"])
            for k in map_df["Location"]]
    return result


def write_failure_maps(map_df, classes, chid, t_end_s):
    """Plot the failure-time and CRE maps; print a short class summary.

    *t_end_s* is the simulation end time; the plan heat map treats locations
    that never fail as failing then.
    """
    points = pd.DataFrame({
        "key": map_df["Location"],
        "member": map_df["Member"],
        "location": [k[len(m) + 1:] for k, m in zip(map_df["Location"], map_df["Member"])],
        "x": map_df["X (m)"], "y": map_df["Y (m)"], "z": map_df["Z (m)"],
    })
    stem = str(FAILURE_MAP_OUTPUT)
    written = []

    tc = classes.get("tcrit")
    print("\nFAILURE MAP")
    if tc is None:
        print("  Critical temperature not reached at any placed location - "
              "no failure-time map.")
    else:
        first = map_df.loc[map_df["Failure Time (min)"] == tc["t0"], "Location"].tolist()
        written += plot_failure_map(
            f"{stem}_tcrit", points, tc,
            title=f"{chid}: time to reach the critical steel temperature",
            early_is_severe=True, highlight=first,
            highlight_label=f"First failure ({tc['t0']} min)",
            unassigned_label="Not reached / no criterion")
        field = plan_failure_field(map_df, tc, t_end_min=t_end_s / 60.0,
                                   cell=HEATMAP_CELL)
        if field is None:
            print("  Plan collapses to a line - no plan heat map or contours.")
        else:
            title = f"{chid}: time to reach the critical steel temperature"
            written.append(plot_failure_heatmap(
                f"{stem}_tcrit_heatmap.png", map_df, field, title=title))
            written += plot_failure_contours(f"{stem}_tcrit", map_df, field, title=title)
        extra = f" (+{len(first) - 1} more)" if len(first) > 1 else ""
        print(f"  First failure: {tc['t0']} min at {first[0]}{extra}")
        counts = map_df["Failure Class"].value_counts()
        parts = [f"{class_label(s, tc['step'])}: {counts.get(class_label(s, tc['step']), 0)}"
                 for s in tc["starts"]]
        n_none = int(map_df["Failure Class"].isna().sum())
        print(f"  Classes of {tc['step']} min: " + " | ".join(parts)
              + f" | not reached: {n_none}")

    cre = classes.get("cre")
    if cre is not None:
        top = int(map_df["CRE Equivalent Time (min)"].max())
        worst = map_df.loc[map_df["CRE Equivalent Time (min)"] == top, "Location"].tolist()
        written += plot_failure_map(
            f"{stem}_cre", points, cre,
            title=f"{chid}: CRE equivalent time of fire exposure",
            early_is_severe=False, highlight=worst,
            highlight_label=f"Most severe exposure ({top} min)",
            unassigned_label="Beyond CRE search cap")
        print(f"  CRE equivalent time: {cre['t0']}–{top} min "
              f"(classes of {cre['step']} min)")
    for path in written:
        print(f"  Saved {path}")


# ============================================================
# MAIN
# ============================================================

def main(argv=None):
    args = parse_args(argv)
    if args.config:
        print(f"Config:  {args.config}")
        apply_config(load_config(args.config))
    if args.input:
        globals()["INPUT_CSV"] = args.input
    if args.chid:
        globals()["CHID"] = args.chid
    if args.fds:
        globals()["FDS_INPUT"] = args.fds
    if args.failure_interval is not None:
        globals()["FAILURE_INTERVAL"] = args.failure_interval
    if args.failure_classes is not None:
        globals()["FAILURE_CLASSES"] = args.failure_classes
    if args.heatmap_cell is not None:
        globals()["HEATMAP_CELL"] = args.heatmap_cell
    if not HEATMAP_CELL > 0:
        raise SystemExit(f"Heat-map cell size must be positive (got {HEATMAP_CELL}).")
    globals()["FAILURE_INTERVAL"] = _whole_number_setting(
        FAILURE_INTERVAL, "failure interval (min)")
    globals()["FAILURE_CLASSES"] = _whole_number_setting(
        FAILURE_CLASSES, "failure classes")

    # Resolve the scenario CHID and prefix any output not set explicitly.
    chid = CHID or derive_chid(INPUT_CSV)
    for name in _OUTPUT_KEYS.values():
        if name not in _OUTPUT_OVERRIDDEN:
            globals()[name] = prefix_with_chid(globals()[name], chid)
    print(f"CHID:    {chid}")

    print(f"Reading: {INPUT_CSV}")
    df, time_col = read_devc_csv(INPUT_CSV)
    time = df[time_col].to_numpy(dtype=float)

    device_columns = [c for c in df.columns if c != time_col]
    groups = build_location_groups(device_columns)
    if not groups:
        raise RuntimeError("No AST devices matched the naming convention "
                           "'<Member>_F<face>_<location>'.")
    print(f"Found {len(groups)} member locations "
          f"across {len({m for m, _ in groups})} members.")

    # --- integrate every location ---------------------------------------
    location_data: dict[str, np.ndarray] = {"Time": time}
    location_ast: dict[str, np.ndarray] = {}
    member_locations: dict[str, list[str]] = {}

    for (member, location), devs in sorted(groups.items()):
        cfg = get_member_config(member)
        ast = aggregate_ast(df, devs).to_numpy(dtype=float)

        if cfg.protected:
            steel = solve_protected_steel(time, ast, cfg.AmV, cfg.protection)
        else:
            steel = solve_unprotected_steel(time, ast, cfg.AmV)

        key = f"{member}_{location}"
        location_data[key] = steel
        location_ast[key] = ast
        member_locations.setdefault(member, []).append(key)

    # Build in one shot to avoid DataFrame fragmentation.
    location_results = pd.DataFrame(location_data)
    location_results.to_csv(LOCATION_OUTPUT, index=False)
    print(f"Saved {LOCATION_OUTPUT}")

    # --- device coordinates (failure maps) ------------------------------
    coords = locate_devices(FDS_INPUT, groups, set(device_columns)) if FDS_INPUT else {}

    def hottest_xyz(key):
        if key not in coords:
            return {}
        xyz = coords[key][0]
        return {"Hottest X (m)": round(float(xyz[0]), 4),
                "Hottest Y (m)": round(float(xyz[1]), 4),
                "Hottest Z (m)": round(float(xyz[2]), 4)}

    # --- hottest location + assessment per member -----------------------
    summary_rows = []
    peaks_data: dict[str, np.ndarray] = {"Time": time}
    member_hottest = {}

    for member in sorted(member_locations):
        cfg = get_member_config(member)
        locs = member_locations[member]

        hottest_key = max(locs, key=lambda k: location_results[k].max())
        series = location_results[hottest_key].to_numpy(dtype=float)
        tmax = float(series.max())

        crit = cfg.critical_temp
        t_crit = find_critical_time(time, series, crit)
        util = (tmax / crit) if crit else None

        member_hottest[member] = (hottest_key, series, cfg)
        peaks_data[member] = series

        row = {
            "Member": member,
            "Hottest Location": hottest_key,
            **hottest_xyz(hottest_key),
            "Maximum Temperature (C)": round(tmax, 1),
            "Critical Temperature (C)": crit,
            "Utilisation": round(util, 3) if util is not None else None,
            "Critical Time (s)": round(t_crit, 1) if t_crit is not None else None,
            "Protected": cfg.protected,
            "Protection Thickness (m)": cfg.protection["thickness"] if cfg.protected else None,
        }

        if TIME_EQUIVALENCE:
            # CRE-equivalent time from the most severe exposure (max cumulative
            # radiant energy) among the member's locations.
            worst_loc = max(locs, key=lambda k: cumulative_radiant_energy(time, location_ast[k]))
            t_eq = equivalent_time_cre(time, location_ast[worst_loc])
            row["CRE Equivalent Time (s)"] = round(t_eq, 1) if t_eq is not None else None
            row["CRE Equivalent Time (min)"] = round(t_eq / 60.0, 1) if t_eq is not None else None
            row["CRE Exposure Location"] = worst_loc

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    peaks = pd.DataFrame(peaks_data)

    # --- configuration audit sheet --------------------------------------
    config_rows = []
    for prefix, props in MEMBER_PROPERTIES.items():
        prot = props.get("protection", {}) if props.get("protected") else {}
        config_rows.append({
            "Prefix": prefix,
            "Am/V (1/m)": props.get("AmV"),
            "Critical Temperature (C)": props.get("critical_temp"),
            "Protected": bool(props.get("protected", False)),
            "Protection Thickness (m)": prot.get("thickness"),
            "Protection Conductivity (W/mK)": prot.get("conductivity"),
            "Protection Density (kg/m3)": prot.get("rho"),
            "Protection Specific Heat (J/kgK)": prot.get("cp"),
        })
    config_df = pd.DataFrame(config_rows)

    map_df = None
    map_classes = {}
    if coords:
        map_df = build_location_map(time, member_locations, location_results,
                                    location_ast, coords)
        map_classes = classify_location_map(map_df)

    write_excel(EXCEL_OUTPUT, summary_df, location_results, peaks, config_df,
                map_df)
    print(f"Saved {EXCEL_OUTPUT}")

    # --- plots ----------------------------------------------------------
    for member, (hottest_key, series, cfg) in member_hottest.items():
        plot_member(member, time, series, hottest_key, cfg)
    print(f"Saved {len(member_hottest)} plot(s) to {PLOT_DIR}/")

    # --- console report -------------------------------------------------
    print("\n" + "=" * 64)
    print("STEEL FIRE ASSESSMENT SUMMARY")
    print("=" * 64)
    for r in summary_df.sort_values("Member").to_dict("records"):
        util_val = r["Utilisation"]
        util = f"{util_val:.2f}" if not pd.isna(util_val) else "n/a"
        t_crit = r["Critical Time (s)"]
        crit_t = f"{t_crit:.0f} s" if not pd.isna(t_crit) else "not reached"
        status = ""
        if not pd.isna(r["Critical Temperature (C)"]):
            status = "EXCEEDED" if not pd.isna(util_val) and util_val >= 1.0 else "OK"
        teq = ""
        if TIME_EQUIVALENCE:
            t_eq = r.get("CRE Equivalent Time (min)")
            teq = f"  t_eq(CRE)={t_eq:>5.1f} min" if not pd.isna(t_eq) else "  t_eq(CRE)=  >cap"
        print(f"  {r['Member']:<8} hottest {r['Hottest Location']:<14} "
              f"Tmax={r['Maximum Temperature (C)']:>6.1f} C  "
              f"util={util:>5}  t_crit={crit_t:<12} {status}{teq}")
    print("=" * 64)

    if map_df is not None and not map_df.empty:
        write_failure_maps(map_df, map_classes, chid, float(time[-1]))


if __name__ == "__main__":
    main()
