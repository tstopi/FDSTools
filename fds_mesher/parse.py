"""FDS namelist reader and extraction of the entities the mesher needs.

Reads &MESH, &MULT, &OBST, &HOLE and &VENT. Everything else is ignored.
"""

import bisect
import re
from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np


class Namelist(NamedTuple):
    """One parsed namelist: raw string values keyed by upper-case name."""
    name: str
    params: dict
    line: int        # 1-based first line
    line_end: int    # 1-based last line
    span: tuple      # (start, end) character offsets in the source text


@dataclass
class Box:
    """An XB box (OBST, HOLE, MESH) or a VENT (xb or an MB face name)."""
    xb: object = None
    surf: str = None
    mb: str = None
    line: int = 0


@dataclass
class Model:
    text: str = ""
    namelists: list = field(default_factory=list)
    meshes: list = field(default_factory=list)
    obsts: list = field(default_factory=list)
    holes: list = field(default_factory=list)
    vents: list = field(default_factory=list)
    geoms: list = field(default_factory=list)
    warnings: list = field(default_factory=list)


# ============================================================
# NAMELIST SCANNER
# ============================================================

_COMMENT_LINE = re.compile(r"^[ \t]*!.*$", re.M)
_START = re.compile(r"&([A-Za-z_]\w*)")
_TOKEN = re.compile(r"'[^']*'|\"[^\"]*\"|!|/")
_QUOTED_OR_COMMENT = re.compile(r"'[^']*'|\"[^\"]*\"|!.*")
_QUOTED = re.compile(r"'[^']*'|\"[^\"]*\"")
_KEY = re.compile(r"([A-Za-z_]\w*)\s*(?:\([^)]*\))?\s*=")


def parse_namelists(text):
    """Return a list of Namelist. A '/' inside a quoted string does not end one."""
    clean = _COMMENT_LINE.sub(lambda m: " " * len(m.group()), text)
    newlines = [m.start() for m in re.finditer("\n", clean)]

    def line_of(pos):
        return bisect.bisect_left(newlines, pos) + 1

    out, pos = [], 0
    while True:
        m = _START.search(clean, pos)
        if not m:
            break
        p = m.end()
        while True:
            t = _TOKEN.search(clean, p)
            if t is None:
                raise ValueError(f"unterminated &{m.group(1)} namelist at "
                                 f"line {line_of(m.start())}")
            if t.group() == "/":
                end, p = t.start(), t.end()
                break
            if t.group() == "!":
                nxt = clean.find("\n", t.end())
                p = len(clean) if nxt < 0 else nxt
            else:
                p = t.end()
        body = _QUOTED_OR_COMMENT.sub(
            lambda g: g.group() if g.group()[0] in "'\"" else "",
            clean[m.end():end])
        out.append(Namelist(m.group(1).upper(), _split_body(body),
                            line_of(m.start()), line_of(end),
                            (m.start(), p)))
        pos = p
    return out


def _split_body(body):
    masked = _QUOTED.sub(lambda g: "x" * len(g.group()), body)
    keys = list(_KEY.finditer(masked))
    params = {}
    for i, k in enumerate(keys):
        e = keys[i + 1].start() if i + 1 < len(keys) else len(body)
        params[k.group(1).upper()] = body[k.end():e].strip().rstrip(",").strip()
    return params


# ============================================================
# VALUE HELPERS
# ============================================================

def get_floats(raw):
    """Parse '1, 2.5D0, 3*0.0' into a list of floats."""
    out = []
    raw = re.sub(r"(?<=[\d.])[dD](?=[+-]?\d)", "E", raw or "")
    for tok in re.split(r"[,\s]+", raw.strip()):
        if not tok:
            continue
        if "*" in tok:
            n, v = tok.split("*", 1)
            out.extend([float(v)] * int(n))
        else:
            out.append(float(tok))
    return out


def get_str(raw):
    """First quoted string of a value (or the bare value)."""
    if raw is None:
        return None
    m = re.search(r"'([^']*)'|\"([^\"]*)\"", raw)
    if m:
        return m.group(1) if m.group(1) is not None else m.group(2)
    return raw.strip()


def get_int(raw, default=None):
    return int(round(get_floats(raw)[0])) if raw else default


def get_flag(raw):
    return bool(raw) and raw.strip().upper().lstrip(".")[:1] == "T"


# ============================================================
# MODEL EXTRACTION
# ============================================================

def parse_text(text):
    """Parse FDS input text into a Model."""
    model = Model(text=text, namelists=parse_namelists(text))
    n_removable = n_hole_ctrl = 0
    for nl in model.namelists:
        p = nl.params
        if nl.name == "MESH":
            model.meshes.append(_box(nl, "XB"))
        elif nl.name == "OBST":
            model.obsts.append(_box(nl, "XB", surf=True))
            if any(k in p for k in ("CTRL_ID", "DEVC_ID")) or \
                    get_flag(p.get("REMOVABLE")):
                n_removable += 1
        elif nl.name == "HOLE":
            model.holes.append(_box(nl, "XB"))
            if any(k in p for k in ("CTRL_ID", "DEVC_ID")):
                n_hole_ctrl += 1
        elif nl.name == "VENT":
            model.vents.append(_vent(nl))
    if n_removable:
        model.warnings.append(f"{n_removable} OBST are removable / controlled; "
                              "treated as present (solid at t=0)")
    if n_hole_ctrl:
        model.warnings.append(f"{n_hole_ctrl} HOLE are controlled; "
                              "treated as open")
    return model


def _box(nl, key, surf=False):
    raw = nl.params.get(key)
    if raw is None:
        raise ValueError(f"&{nl.name} at line {nl.line} has no {key}")
    xb = np.array(get_floats(raw), float)
    if xb.size != 6:
        raise ValueError(f"&{nl.name} at line {nl.line}: {key} needs 6 values")
    xb = xb.reshape(3, 2)
    xb.sort(axis=1)
    return Box(xb.ravel(), get_str(nl.params.get("SURF_ID")) if surf else None,
               None, nl.line)


def _vent(nl):
    p = nl.params
    surf = get_str(p.get("SURF_ID"))
    if "XB" in p:
        b = _box(nl, "XB")
        b.surf = surf
        return b
    if "MB" in p:
        return Box(None, surf, get_str(p["MB"]).upper(), nl.line)
    raise ValueError(f"&VENT at line {nl.line} has neither XB nor MB")
