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

The CHID is taken from ``--chid`` / the config, or derived from the input
filename (``<CHID>_devc.csv``), so several scenarios can be post-processed in
the same directory without overwriting each other.

Device naming convention (``<Member>_<Face>_<Location>``)::

    AP_F1_001      YP_F4_012      D_01_F2_007      AP-2_F3_005

Members are matched to a configuration entry by longest family-prefix match:
``AP``, ``AP-1`` and ``AP-15`` all inherit the ``AP`` entry, and ``D_01``,
``D_02`` ... inherit ``D_``. Matching stops at a letter boundary, so ``AP``
does not swallow an unrelated ``APRON``.

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
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: no display required
import matplotlib.pyplot as plt
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

def write_excel(path, summary_df, location_df, peaks_df, config_df):
    """Write the four-sheet results workbook with conditional formatting."""
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter

    green = PatternFill("solid", fgColor="C6EFCE")
    red = PatternFill("solid", fgColor="FFC7CE")
    header_fill = PatternFill("solid", fgColor="305496")
    header_font = Font(bold=True, color="FFFFFF")

    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        summary_df.to_excel(xl, sheet_name="Member Summary", index=False)
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
    return parser.parse_args(argv)


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
        member_locations.setdefault(member, []).append(key)

    # Build in one shot to avoid DataFrame fragmentation.
    location_results = pd.DataFrame(location_data)
    location_results.to_csv(LOCATION_OUTPUT, index=False)
    print(f"Saved {LOCATION_OUTPUT}")

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

        summary_rows.append({
            "Member": member,
            "Hottest Location": hottest_key,
            "Maximum Temperature (C)": round(tmax, 1),
            "Critical Temperature (C)": crit,
            "Utilisation": round(util, 3) if util is not None else None,
            "Critical Time (s)": round(t_crit, 1) if t_crit is not None else None,
            "Protected": cfg.protected,
            "Protection Thickness (m)": cfg.protection["thickness"] if cfg.protected else None,
        })

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

    write_excel(EXCEL_OUTPUT, summary_df, location_results, peaks, config_df)
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
        print(f"  {r['Member']:<8} hottest {r['Hottest Location']:<14} "
              f"Tmax={r['Maximum Temperature (C)']:>6.1f} C  "
              f"util={util:>5}  t_crit={crit_t:<12} {status}")
    print("=" * 64)


if __name__ == "__main__":
    main()
