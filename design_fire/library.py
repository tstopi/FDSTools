"""Design fire library: JSON entries shipped with the package plus the user's.

Each entry is one JSON file::

    {"name": "...", "description": "...",
     "curve": {"type": "power-law", ...} or {"type": "tabulated", ...},
     "fuels": [["Polyurethane foam (flexible)", 1.0]],
     "area": 4.0}

``fuels`` (names from ``fuels.FUELS`` with mass fractions) and ``area``
(m^2) are optional. Built-in entries live in ``design_fire/library/``; user
entries in ``$FDSTOOLS_FIRE_LIBRARY`` or ``~/.fdstools/design_fires/``.
"""

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from .curves import TabulatedCurve, curve_from_dict
from .fuels import FUELS, species_id

BUILTIN_DIR = Path(__file__).with_name("library")


def user_dir():
    env = os.environ.get("FDSTOOLS_FIRE_LIBRARY")
    if env:
        return Path(env).expanduser()
    return Path.home() / ".fdstools" / "design_fires"


@dataclass
class LibraryEntry:
    name: str
    curve: object                                 # a curves.Curve
    fuels: list = field(default_factory=list)     # [(fuel name, fraction)]
    area: float = None                            # m^2, optional
    description: str = ""
    builtin: bool = False
    path: Path = None

    def __post_init__(self):
        if not self.name.strip():
            raise ValueError("library entry needs a name")
        unknown = [n for n, _ in self.fuels if n not in FUELS]
        if unknown:
            raise ValueError(f"{self.name}: unknown fuels {unknown}")
        if any(w <= 0 for _, w in self.fuels):
            raise ValueError(f"{self.name}: mass fractions must be positive")
        if self.area is not None and self.area <= 0:
            raise ValueError(f"{self.name}: fire area must be positive")

    def components(self):
        """[(Fuel, mass fraction)] for ``writer.DesignFire``."""
        return [(FUELS[n], w) for n, w in self.fuels]

    def to_dict(self):
        d = {"name": self.name, "description": self.description,
             "curve": self.curve.to_dict(),
             "fuels": [[n, w] for n, w in self.fuels]}
        if self.area is not None:
            d["area"] = self.area
        return d

    @classmethod
    def from_dict(cls, data, builtin=False, path=None):
        return cls(name=data["name"], curve=curve_from_dict(data["curve"]),
                   fuels=[(n, float(w)) for n, w in data.get("fuels", [])],
                   area=data.get("area"),
                   description=data.get("description", ""),
                   builtin=builtin, path=path)


def _load_dir(directory, builtin):
    entries, errors = [], []
    if not directory.is_dir():
        return entries, errors
    for path in sorted(directory.glob("*.json")):
        try:
            with open(path, encoding="utf-8") as fh:
                entries.append(LibraryEntry.from_dict(json.load(fh), builtin,
                                                      path))
        except (OSError, ValueError, KeyError, TypeError) as e:
            errors.append(f"{path.name}: {e}")
    return entries, errors


def load_library(builtin_dir=None, user=None):
    """(entries, errors): built-in entries first, then the user's.

    Files that fail to load are skipped and reported in ``errors``.
    """
    entries, errors = _load_dir(Path(builtin_dir or BUILTIN_DIR), True)
    more, errs = _load_dir(Path(user or user_dir()), False)
    return entries + more, errors + errs


def entry_path(name, directory=None):
    stem = species_id(name).lower() or "design_fire"
    return Path(directory or user_dir()) / f"{stem}.json"


def save_entry(entry, directory=None):
    """Write ``entry`` to the user library and return its path."""
    path = entry_path(entry.name, directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(entry.to_dict(), fh, indent=2)
        fh.write("\n")
    entry.path = path
    entry.builtin = False
    return path


def delete_entry(entry):
    if entry.builtin:
        raise ValueError("built-in library entries cannot be deleted")
    if entry.path is not None:
        Path(entry.path).unlink()


def read_csv_curve(path):
    """Tabulated curve from a two-column time (s), HRR (kW) text file.

    Commas, tabs or spaces separate the columns, or semicolons with decimal
    commas; lines that do not start with two numbers (headers, comments)
    are skipped.
    """
    table = []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            if ";" in line:     # semicolon columns may use decimal commas
                cells = [c.replace(",", ".") for c in line.split(";")]
            else:
                cells = re.split(r"[,\s]+", line.strip())
            try:
                table.append((float(cells[0]), float(cells[1])))
            except (ValueError, IndexError):
                continue
    return TabulatedCurve(table)
