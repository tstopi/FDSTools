"""Fuel property library.

Values are well-ventilated effective (chemical) heats of combustion and
CO/soot yields from the SFPE Handbook (Tewarson's tables). Formulas are per
monomer or representative unit. Check them against the edition you cite
before using them in a design.
"""

from dataclasses import dataclass, field

ATOMIC_MASS = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007, "Cl": 35.45}
ELEMENTS = tuple(ATOMIC_MASS)


@dataclass
class Fuel:
    name: str
    formula: dict                # element -> atoms per molecule
    heat_of_combustion: float    # kJ/kg, effective
    co_yield: float              # kg CO / kg fuel
    soot_yield: float            # kg soot / kg fuel
    fds_id: str = ""             # FDS species ID; blank = derived from name
    predefined: bool = False     # True when FDS knows the species by ID
    note: str = field(default="", compare=False)

    def __post_init__(self):
        unknown = set(self.formula) - set(ELEMENTS)
        if unknown:
            raise ValueError(f"{self.name}: unsupported elements {sorted(unknown)}")
        if self.formula.get("C", 0) <= 0:
            raise ValueError(f"{self.name}: fuel must contain carbon")
        if self.heat_of_combustion <= 0:
            raise ValueError(f"{self.name}: heat of combustion must be positive")
        if self.co_yield < 0 or self.soot_yield < 0:
            raise ValueError(f"{self.name}: yields must be non-negative")
        if not self.fds_id:
            self.fds_id = species_id(self.name)

    @property
    def molar_mass(self):
        """g/mol"""
        return sum(ATOMIC_MASS[e] * n for e, n in self.formula.items())

    def atoms(self, element):
        return self.formula.get(element, 0.0)


def species_id(name):
    """FDS-safe species ID from a free-form name."""
    out = "".join(c if c.isalnum() else "_" for c in name.upper())
    return "_".join(p for p in out.split("_") if p)


def formula_string(formula, digits=4):
    """'C2H3Cl' style string; non-integer counts keep up to ``digits`` decimals."""
    parts = []
    for e in ELEMENTS:
        n = formula.get(e, 0.0)
        if n <= 0:
            continue
        if abs(n - 1) < 10 ** -digits:
            parts.append(e)
        elif abs(n - round(n)) < 10 ** -digits:
            parts.append(f"{e}{int(round(n))}")
        else:
            parts.append(f"{e}{n:.{digits}f}".rstrip("0").rstrip("."))
    return "".join(parts)


FUELS = {f.name: f for f in [
    Fuel("Methane", {"C": 1, "H": 4}, 49600, 0.0, 0.0, "METHANE", True),
    Fuel("Propane", {"C": 3, "H": 8}, 43700, 0.005, 0.024, "PROPANE", True),
    Fuel("Heptane", {"C": 7, "H": 16}, 41200, 0.010, 0.037, "N-HEPTANE", True),
    Fuel("Polyethylene", {"C": 2, "H": 4}, 38400, 0.024, 0.060),
    Fuel("Polypropylene", {"C": 3, "H": 6}, 38600, 0.024, 0.059),
    Fuel("Polystyrene", {"C": 8, "H": 8}, 27000, 0.060, 0.164),
    Fuel("PMMA", {"C": 5, "H": 8, "O": 2}, 24200, 0.010, 0.022),
    Fuel("Nylon", {"C": 6, "H": 11, "O": 1, "N": 1}, 27100, 0.038, 0.075),
    Fuel("PVC", {"C": 2, "H": 3, "Cl": 1}, 5700, 0.063, 0.172),
    Fuel("Polyurethane foam (flexible)",
         {"C": 1, "H": 1.75, "O": 0.32, "N": 0.07}, 17600, 0.031, 0.227),
    Fuel("Wood (red oak)", {"C": 1, "H": 1.7, "O": 0.72, "N": 0.001},
         12400, 0.004, 0.015),
]}
