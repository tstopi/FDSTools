"""FDS namelist reader and extraction of the entities the mesher needs.

Reads &MESH, &MULT, &OBST, &HOLE, &VENT and triangle &GEOM. Everything else
(&MOVE, &INIT, &ZONE, ...) is ignored.
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
class Geom:
    id: str
    verts: np.ndarray    # (nv, 3)
    faces: np.ndarray    # (nt, 3) zero-based vertex indices


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
    raw = re.sub(r"(?<=[\d.])[dD](?=[+-]?\d)", "E", raw or "")
    if "*" not in raw:      # fast path (large VERTS / FACES arrays)
        return np.array(raw.replace(",", " ").split(), float).tolist()
    out = []
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
    """Parse FDS input text into a Model (MULT_ID already expanded)."""
    model = Model(text=text, namelists=parse_namelists(text))
    mults = {}
    for nl in model.namelists:
        if nl.name == "MULT":
            mults[get_str(nl.params.get("ID"))] = nl
    n_removable = n_hole_ctrl = 0
    for nl in model.namelists:
        p = nl.params
        if nl.name == "MESH":
            model.meshes += _expand(_box(nl, "XB"), nl, mults)
        elif nl.name == "OBST":
            model.obsts += _expand(_box(nl, "XB", surf=True), nl, mults)
            if any(k in p for k in ("CTRL_ID", "DEVC_ID")) or \
                    get_flag(p.get("REMOVABLE")):
                n_removable += 1
        elif nl.name == "HOLE":
            model.holes += _expand(_box(nl, "XB"), nl, mults)
            if any(k in p for k in ("CTRL_ID", "DEVC_ID")):
                n_hole_ctrl += 1
        elif nl.name == "VENT":
            model.vents += _expand(_vent(nl), nl, mults)
        elif nl.name == "GEOM":
            model.geoms.append(_geom(nl))
    if n_removable:
        model.warnings.append(f"{n_removable} OBST are removable / controlled; "
                              "treated as present (solid at t=0)")
    if n_hole_ctrl:
        model.warnings.append(f"{n_hole_ctrl} HOLE are controlled; "
                              "treated as open")
    return model


# ============================================================
# MULT
# ============================================================

_MULT_UNSUPPORTED = ("DXB", "DX0", "DY0", "DZ0")


def mult_offsets(nl):
    """(n, 3) array of translation offsets of a &MULT namelist.

    Supports DX/DY/DZ with I/J/K_LOWER/UPPER, or N_LOWER/N_UPPER (offset
    n*(DX,DY,DZ)). DXB, DX0-style origins and *_SKIP / *_NEW are rejected.
    """
    p = nl.params
    for k in p:
        if k in _MULT_UNSUPPORTED or k.endswith("_SKIP") or k.endswith("_NEW"):
            raise NotImplementedError(
                f"&MULT ID={get_str(p.get('ID'))!r} (line {nl.line}): {k} is "
                "not supported")
    d = np.array([get_floats(p.get(k, "0"))[0] for k in ("DX", "DY", "DZ")])
    if "N_LOWER" in p or "N_UPPER" in p:
        n = np.arange(get_int(p.get("N_LOWER"), 0),
                      get_int(p.get("N_UPPER"), 0) + 1)
        return n[:, None] * d[None, :]
    rng = [np.arange(get_int(p.get(a + "_LOWER"), 0),
                     get_int(p.get(a + "_UPPER"), 0) + 1) for a in "IJK"]
    i, j, k = np.meshgrid(*rng, indexing="ij")
    idx = np.stack([i.ravel(), j.ravel(), k.ravel()], axis=1)
    return idx * d[None, :]


def _expand(box, nl, mults):
    """Copies of a Box for its MULT_ID (or [box] if it has none)."""
    mid = get_str(nl.params.get("MULT_ID"))
    if mid is None:
        return [box]
    if mid not in mults:
        raise ValueError(f"&{nl.name} line {nl.line}: unknown MULT_ID '{mid}'")
    if box.xb is None:
        raise NotImplementedError(f"&VENT line {nl.line}: MULT_ID with MB")
    out = []
    for off in mult_offsets(mults[mid]):
        out.append(Box(box.xb + np.repeat(off, 2), box.surf, box.mb, box.line))
    return out


# ============================================================
# GEOM
# ============================================================

_GEOM_UNSUPPORTED = ("BINARY_FILE", "ZVALS", "XB", "POLY", "IJK")


def resolve_stride(n_values):
    """FACES stride: 4 (3 vertices + surface index, BlenderFDS) or 3.

    Same rule as geom_to_ast_devc.py: prefer 4 when the length allows.
    """
    if n_values % 4 == 0:
        return 4
    if n_values % 3 == 0:
        return 3
    raise ValueError(f"FACES length {n_values} is divisible by neither 3 nor 4")


def _geom(nl):
    p = nl.params
    gid = get_str(p.get("ID")) or f"line {nl.line}"
    bad = [k for k in p if k in _GEOM_UNSUPPORTED or k.startswith(("SPHERE_",
                                                                  "CYLINDER_"))]
    if bad:
        raise NotImplementedError(f"GEOM '{gid}': {bad[0]} is not supported "
                                  "(only VERTS/FACES triangle meshes)")
    if "VERTS" not in p or "FACES" not in p:
        raise NotImplementedError(f"GEOM '{gid}': no VERTS/FACES")
    vals = get_floats(p["VERTS"])
    if len(vals) % 3 or not vals:
        raise ValueError(f"GEOM '{gid}': VERTS length {len(vals)} is not a "
                         "multiple of 3")
    verts = np.array(vals).reshape(-1, 3)
    ints = np.rint(np.array(get_floats(p["FACES"]))).astype(np.int64)
    stride = resolve_stride(len(ints))
    faces = ints.reshape(-1, stride)[:, :3] - 1
    if faces.size and (faces.min() < 0 or faces.max() >= len(verts)):
        raise ValueError(f"GEOM '{gid}': face index out of range "
                         f"(stride {stride})")
    return Geom(gid, verts, faces)


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
