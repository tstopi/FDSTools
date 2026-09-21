#!/usr/bin/env python3
"""Generate FDS AST devices on the faces of GEOM members.

Reads the ``&GEOM`` namelists from an FDS input file, reconstructs the
prismatic (usually rectangular) members, and writes ``&DEVC`` lines carrying
adiabatic-surface-temperature (AST) gas devices spread along every exposed
longitudinal face. The device IDs follow the convention consumed by
``steel_temperature.py``::

    <Member>_F<Face>_<Location>      e.g.  AP_F1_003,  D_01_F4_012

The ``<Location>`` index is a **shared axial station**: for a given member,
``_003`` denotes the same longitudinal position on every face, so the steel
post-processor can group all faces at one cross-section correctly - even for
mitred or unequal-length faces.

For each member the workflow is:

1. Reconstruct quad faces by merging coplanar triangle pairs.
2. Find the member's long axis by PCA of its vertices.
3. Drop the end-cap faces (normal roughly parallel to the long axis).
4. Along each remaining face, place devices at the shared axial stations that
   fall within that face, offset slightly into the gas and oriented outward.
5. Skip any device that falls inside another member's bounding box.

Usage
-----
    python geom_to_ast_devc.py --input model-GEOM.fds
    python geom_to_ast_devc.py -i model-GEOM.fds -o model_ast_devc.fds --delta 0.25
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# ============================================================
# DEFAULT SETTINGS  (override on the command line)
# ============================================================

INPUT_FILE = "Uusi_puoli_ver0p2-GEOM.fds"

DELTA = 0.50            # spacing of stations along the member length [m]
OFFSET = 0.005          # offset of the device outside the face [m]

# AST gas-phase quantity; this is what the device measures.
QUANTITY = "ADIABATIC SURFACE TEMPERATURE GAS"

# Faces whose |normal . long-axis| exceeds this are treated as end caps.
ENDCAP_ALIGNMENT = 0.90

# Number of integers per face in the FACES array. BlenderFDS writes 4
# (3 vertex indices + 1 surface index); a bare triangle list is 3. "auto"
# infers it from the array length, but the length is ambiguous when it is a
# multiple of 12, so the BlenderFDS default of 4 is used here.
FACES_STRIDE = 4


# ============================================================
# PARSE GEOMS
# ============================================================

def _strip_comments(text):
    """Remove full-line FDS comments (lines whose first token starts with !)."""
    out = []
    for line in text.splitlines():
        if line.lstrip().startswith("!"):
            continue
        out.append(line)
    return "\n".join(out)


def _floats(blob):
    return [float(x) for x in blob.replace("\n", " ").split(",") if x.strip()]


def _resolve_stride(n_faces_values, n_verts):
    """Determine the FACES stride (3 or 4)."""
    if FACES_STRIDE in (3, 4):
        stride = FACES_STRIDE
        if n_faces_values % stride != 0:
            raise ValueError(
                f"FACES length {n_faces_values} is not divisible by "
                f"the configured stride {stride}.")
        return stride
    # auto: prefer 4 (BlenderFDS) when possible, else 3.
    if n_faces_values % 4 == 0:
        if n_faces_values % 3 == 0:
            print("  FACES length is a multiple of both 3 and 4; assuming "
                  "stride 4 (BlenderFDS). Set FACES_STRIDE=3 if wrong.",
                  file=sys.stderr)
        return 4
    if n_faces_values % 3 == 0:
        return 3
    raise ValueError(
        f"FACES length {n_faces_values} is divisible by neither 3 nor 4; "
        "set FACES_STRIDE explicitly.")


def parse_geoms(filename):
    """Return a list of {id, verts (N,3), faces (list of 0-based triplets)}."""
    path = Path(filename)
    if not path.exists():
        raise SystemExit(f"Input file not found: {filename}")
    text = _strip_comments(path.read_text())
    blocks = re.findall(r"&GEOM\b(.*?)/", text, re.DOTALL)

    geoms = []
    for block in blocks:
        id_match = re.search(r"ID\s*=\s*'([^']+)'", block)
        verts_match = re.search(r"VERTS\s*=(.*?)FACES\s*=", block, re.DOTALL)
        faces_match = re.search(r"FACES\s*=(.*)", block, re.DOTALL)
        if not (id_match and verts_match and faces_match):
            continue

        verts_data = _floats(verts_match.group(1))
        if len(verts_data) % 3 != 0 or not verts_data:
            print(f"  Skipping GEOM {id_match.group(1)!r}: VERTS length "
                  f"{len(verts_data)} is not a multiple of 3.", file=sys.stderr)
            continue
        verts = np.array(verts_data).reshape((-1, 3))

        face_ints = [int(round(v)) for v in _floats(faces_match.group(1))]
        if not face_ints:
            continue
        stride = _resolve_stride(len(face_ints), len(verts))

        faces = []
        for i in range(0, len(face_ints), stride):
            tri = face_ints[i:i + 3]
            if len(tri) == 3:
                # store as 0-based indices
                faces.append([t - 1 for t in tri])

        # validate indices
        max_idx = max((max(t) for t in faces), default=-1)
        if max_idx >= len(verts) or min((min(t) for t in faces), default=0) < 0:
            print(f"  Skipping GEOM {id_match.group(1)!r}: face index out of "
                  f"range (stride={stride}). Try setting FACES_STRIDE.",
                  file=sys.stderr)
            continue

        geoms.append({"id": id_match.group(1), "verts": verts, "faces": faces})

    return geoms


# ============================================================
# ORIENTED BOUNDING BOX (for occlusion tests)
# ============================================================

def build_obb(vertices):
    center = vertices.mean(axis=0)
    x = vertices - center
    _, eigvecs = np.linalg.eigh(np.cov(x.T))
    local = x @ eigvecs
    return {
        "center": center,
        "R": eigvecs,
        "mins": local.min(axis=0),
        "maxs": local.max(axis=0),
    }


def point_inside_obb(pt, box, tol=1e-6):
    p = (pt - box["center"]) @ box["R"]
    return np.all(p >= box["mins"] - tol) and np.all(p <= box["maxs"] + tol)


def inside_any_member(pt, boxes, own_index):
    for i, box in enumerate(boxes):
        if i != own_index and point_inside_obb(pt, box):
            return True
    return False


# ============================================================
# TRIANGLE / GEOMETRY UTILITIES
# ============================================================

def tri_normal(tri):
    n = np.cross(tri[1] - tri[0], tri[2] - tri[0])
    mag = np.linalg.norm(n)
    return None if mag < 1e-12 else n / mag


def coplanar(n1, n2, tol=1e-3):
    if n1 is None or n2 is None:
        return False
    return abs(abs(np.dot(n1, n2)) - 1.0) < tol


def member_axis(vertices):
    """Return the unit long-axis of the member by PCA of its vertices."""
    x = vertices - vertices.mean(axis=0)
    eigvals, eigvecs = np.linalg.eigh(np.cov(x.T))
    axis = eigvecs[:, np.argmax(eigvals)]
    return axis / np.linalg.norm(axis)


def build_quads(geom):
    """Merge coplanar triangle pairs that share an edge into quad faces."""
    verts, faces = geom["verts"], geom["faces"]
    quads = []
    used = set()

    for i, f1 in enumerate(faces):
        if i in used:
            continue
        nodes1 = set(f1)
        n1 = tri_normal(verts[np.array(f1)])

        for j in range(i + 1, len(faces)):
            if j in used:
                continue
            f2 = faces[j]
            nodes2 = set(f2)
            if len(nodes1 & nodes2) != 2:
                continue
            n2 = tri_normal(verts[np.array(f2)])
            if not coplanar(n1, n2):
                continue
            all_nodes = list(nodes1 | nodes2)
            if len(all_nodes) != 4:
                continue
            quads.append(verts[np.array(all_nodes)])
            used.add(i)
            used.add(j)
            break

    return quads


# ============================================================
# DEVICE GENERATION
# ============================================================

def axial_stations(gmin, gmax, delta):
    """Centred station positions along [gmin, gmax] at spacing *delta*."""
    stations = []
    s = gmin + delta / 2.0
    while s < gmax:
        stations.append(s)
        s += delta
    return stations


def generate_asts(geom, boxes, own_index):
    """Return the &DEVC lines for one member."""
    lines = []
    verts = geom["verts"]
    member_center = verts.mean(axis=0)
    u = member_axis(verts)

    # Shared axial stations from the member's full extent, so a given station
    # index maps to the same cross-section on every face.
    proj_all = verts @ u
    stations = axial_stations(proj_all.min(), proj_all.max(), DELTA)
    tol = 1e-9

    quads = build_quads(geom)
    face_id = 0

    for face in quads:
        n = tri_normal(np.array([face[0], face[1], face[2]]))
        if n is None:
            continue
        face_center = face.mean(axis=0)

        # Ensure the normal points outward from the member centre.
        if np.dot(n, face_center - member_center) < 0:
            n = -n

        # Skip end caps (normal roughly parallel to the long axis).
        if abs(np.dot(n, u)) > ENDCAP_ALIGNMENT:
            continue

        face_id += 1
        proj = face @ u
        smin, smax = proj.min(), proj.max()
        center_proj = float(np.dot(face_center, u))

        for k, s in enumerate(stations, start=1):
            if s < smin - tol or s > smax + tol:
                continue  # station not covered by this face
            p = face_center + (s - center_proj) * u
            ast = p + OFFSET * n
            if inside_any_member(ast, boxes, own_index):
                continue
            x, y, z = ast
            ox, oy, oz = n
            lines.append(
                f"&DEVC ID='{geom['id']}_F{face_id}_{k:03d}', "
                f"QUANTITY='{QUANTITY}', "
                f"XYZ={x:.4f},{y:.4f},{z:.4f}, "
                f"ORIENTATION={ox:.6f},{oy:.6f},{oz:.6f}/"
            )

    return lines


# ============================================================
# CLI / MAIN
# ============================================================

def default_output_name(input_file):
    """Derive an output filename, e.g. model-GEOM.fds -> model_ast_devc.fds."""
    stem = Path(input_file).stem
    if stem.upper().endswith("-GEOM"):
        stem = stem[:-len("-GEOM")]
    return f"{stem}_ast_devc.fds"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate FDS AST &DEVC lines on GEOM member faces.")
    p.add_argument("-i", "--input", default=INPUT_FILE,
                   help="FDS input file containing the &GEOM namelists.")
    p.add_argument("-o", "--output",
                   help="Output .fds include file (default: <base>_ast_devc.fds).")
    p.add_argument("--delta", type=float, default=DELTA,
                   help="Station spacing along the member length [m].")
    p.add_argument("--offset", type=float, default=OFFSET,
                   help="Device offset outside the face [m].")
    p.add_argument("--faces-stride", choices=["auto", "3", "4"],
                   default=str(FACES_STRIDE),
                   help="Integers per face in FACES (BlenderFDS=4). Default 4.")
    p.add_argument("--stdout", action="store_true",
                   help="Also echo the generated DEVC lines to stdout.")
    return p.parse_args(argv)


def main(argv=None):
    global DELTA, OFFSET, FACES_STRIDE
    args = parse_args(argv)
    DELTA = args.delta
    OFFSET = args.offset
    FACES_STRIDE = args.faces_stride if args.faces_stride == "auto" else int(args.faces_stride)
    output = args.output or default_output_name(args.input)

    geoms = parse_geoms(args.input)
    print(f"Found {len(geoms)} GEOM object(s) in {args.input}", file=sys.stderr)
    if not geoms:
        raise SystemExit("No usable &GEOM namelists found.")

    boxes = [build_obb(g["verts"]) for g in geoms]

    all_lines = []
    for idx, geom in enumerate(geoms):
        devcs = generate_asts(geom, boxes, idx)
        all_lines.extend(devcs)
        print(f"  {geom['id']:<10} {len(devcs):>5} devices", file=sys.stderr)

    header = (
        f"! AST devices generated by geom_to_ast_devc.py\n"
        f"! source={args.input}  members={len(geoms)}  "
        f"delta={DELTA}  offset={OFFSET}\n"
    )
    Path(output).write_text(header + "\n".join(all_lines) + "\n")

    print(f"Wrote {len(all_lines)} AST devices to {output}", file=sys.stderr)
    if args.stdout:
        print(header + "\n".join(all_lines))


if __name__ == "__main__":
    main()
