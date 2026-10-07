"""FDS input text for a design fire: &RAMP, &SPEC, &REAC and &SURF."""

from dataclasses import dataclass

from .fuels import formula_string
from .reactions import balance, blend, normalise

# FDS default ambient air composition by mass
AMBIENT = {"OXYGEN": 0.232378, "CARBON DIOXIDE": 0.000595}
PRODUCTS = ["CARBON DIOXIDE", "WATER VAPOR", "CARBON MONOXIDE", "SOOT",
            "HYDROGEN CHLORIDE", "NITROGEN"]


def _num(v, digits=6):
    """Compact FDS real: always has a decimal point."""
    s = f"{v:.{digits}g}"
    if "e" in s:
        mant, exp = s.split("e")
        if "." not in mant:
            mant += "."
        return f"{mant}E{int(exp)}"
    return s if "." in s else s + "."


def _str_list(items):
    return ",".join(f"'{i}'" for i in items)


def _num_list(items, digits=6):
    return ",".join(_num(v, digits) for v in items)


@dataclass
class DesignFire:
    curve: object                 # a curves.Curve
    area: float                   # m^2 burning area
    components: list              # [(Fuel, mass fraction), ...]
    per_fuel_reactions: bool = False
    surf_id: str = "FIRE"
    ramp_id: str = "FIRE_RAMP"

    def __post_init__(self):
        if self.area <= 0:
            raise ValueError("fire area must be positive")
        self.components = normalise(self.components)
        ids = [f.fds_id for f, _ in self.components]
        if len(set(ids)) != len(ids):
            raise ValueError("each fuel may appear only once in the mixture")

    @property
    def hrrpua(self):
        """kW/m^2 at peak"""
        return self.curve.peak_hrr / self.area

    @property
    def mixture(self):
        return blend(self.components)

    def reactions(self):
        if self.per_fuel_reactions:
            return [balance(f) for f, _ in self.components]
        return [balance(self.mixture)]

    def to_fds(self):
        reacs = self.reactions()
        lines = [self._header(), ""]
        lines += self._ramp() + [""]
        lines += self._species(reacs) + [""]
        lines += [self._reac(r) for r in reacs] + [""]
        lines.append(self._surf())
        return "\n".join(lines) + "\n"

    def _header(self):
        mix = ", ".join(f"{f.name} {w:.0%}" for f, w in self.components)
        c = self.curve
        out = [f"! Design fire: {c.name}, peak {c.peak_hrr:g} kW, "
               f"area {self.area:g} m2, HRRPUA {self.hrrpua:.4g} kW/m2",
               f"! Fuel: {mix}"]
        out += [f"! {line}" for line in c.describe()]
        return "\n".join(out)

    def _ramp(self):
        return [f"&RAMP ID='{self.ramp_id}', T={_num(t, 5)}, F={_num(f, 5)} /"
                for t, f in self.curve.points()]

    def _species(self, reacs):
        used = set()
        for r in reacs:
            used.update(r.nu)
        lines = ["&SPEC ID='NITROGEN', BACKGROUND=.TRUE. /"]
        for sp, y0 in AMBIENT.items():
            lines.append(f"&SPEC ID='{sp}', MASS_FRACTION_0={_num(y0)} /")
        for sp in PRODUCTS:
            if sp in used and sp not in AMBIENT and sp != "NITROGEN":
                lines.append(f"&SPEC ID='{sp}' /")
        for r in reacs:
            f = r.fuel
            if f.predefined:
                lines.append(f"&SPEC ID='{f.fds_id}' /")
            else:
                lines.append(f"&SPEC ID='{f.fds_id}', "
                             f"FORMULA='{formula_string(f.formula)}' /")
        return lines

    def _reac(self, r):
        sp = list(r.nu)
        return (f"&REAC ID='{r.fuel.fds_id}', FUEL='{r.fuel.fds_id}', "
                f"HEAT_OF_COMBUSTION={_num(r.fuel.heat_of_combustion)},\n"
                f"      SPEC_ID_NU={_str_list(sp)},\n"
                f"      NU={_num_list([r.nu[s] for s in sp], 8)} /")

    def _surf(self):
        if not self.per_fuel_reactions or len(self.components) == 1:
            return (f"&SURF ID='{self.surf_id}', HRRPUA={_num(self.hrrpua)}, "
                    f"RAMP_Q='{self.ramp_id}', COLOR='RED' /")
        # Split the burning rate between fuels by mass fraction so the total
        # heat release still follows HRRPUA * ramp.
        burn = self.hrrpua / self.mixture.heat_of_combustion   # kg/m2/s
        k = len(self.components)
        ids = [f.fds_id for f, _ in self.components]
        flux = [w * burn for _, w in self.components]
        return (f"&SURF ID='{self.surf_id}', COLOR='RED',\n"
                f"      SPEC_ID(1:{k})={_str_list(ids)},\n"
                f"      MASS_FLUX(1:{k})={_num_list(flux, 8)},\n"
                f"      RAMP_MF(1:{k})={_str_list([self.ramp_id] * k)} /")
