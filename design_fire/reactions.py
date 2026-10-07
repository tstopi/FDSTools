"""Stoichiometry for fuels and mass-fraction fuel blends.

One mole of fuel C_x H_y O_z N_n Cl_c burns to

    nu_O2 O2 -> nu_CO2 CO2 + nu_H2O H2O + nu_CO CO + nu_s SOOT
                + c HCl + n/2 N2

CO and soot come from the yields, all chlorine leaves as HCl, the remaining
hydrogen as water and the remaining carbon as CO2. Soot is pure carbon (the
FDS default). Oxygen closes the balance.
"""

from dataclasses import dataclass

from .fuels import ATOMIC_MASS, ELEMENTS, Fuel

SPECIES_FORMULA = {
    "OXYGEN": {"O": 2},
    "CARBON DIOXIDE": {"C": 1, "O": 2},
    "WATER VAPOR": {"H": 2, "O": 1},
    "CARBON MONOXIDE": {"C": 1, "O": 1},
    "SOOT": {"C": 1},
    "HYDROGEN CHLORIDE": {"H": 1, "Cl": 1},
    "NITROGEN": {"N": 2},
}


def molar_mass(formula):
    return sum(ATOMIC_MASS[e] * n for e, n in formula.items())


@dataclass
class Reaction:
    fuel: Fuel
    nu: dict   # species ID -> moles per mole of fuel; fuel and O2 negative

    def element_imbalance(self):
        """element -> (products - reactants) atoms per mole of fuel."""
        out = {}
        for e in ELEMENTS:
            total = -self.fuel.atoms(e)
            for sp, nu in self.nu.items():
                if sp == self.fuel.fds_id:
                    continue
                total += nu * SPECIES_FORMULA[sp].get(e, 0.0)
            out[e] = total
        return out

    def mass_imbalance(self):
        """(products - reactants) g per mole of fuel."""
        m = -self.fuel.molar_mass
        for sp, nu in self.nu.items():
            if sp != self.fuel.fds_id:
                m += nu * molar_mass(SPECIES_FORMULA[sp])
        return m


def balance(fuel):
    """Balanced complete-combustion reaction with CO, soot and HCl."""
    x, y, z, n, c = (fuel.atoms(e) for e in ELEMENTS)
    m_f = fuel.molar_mass
    nu_co = fuel.co_yield * m_f / molar_mass(SPECIES_FORMULA["CARBON MONOXIDE"])
    nu_s = fuel.soot_yield * m_f / molar_mass(SPECIES_FORMULA["SOOT"])
    nu_co2 = x - nu_co - nu_s
    if nu_co2 < 0:
        raise ValueError(f"{fuel.name}: CO and soot yields need more carbon "
                         "than the fuel contains")
    if c > y:
        raise ValueError(f"{fuel.name}: not enough hydrogen to form HCl "
                         "from all the chlorine")
    nu_h2o = (y - c) / 2.0
    nu_o2 = (2 * nu_co2 + nu_h2o + nu_co - z) / 2.0
    if nu_o2 <= 0:
        raise ValueError(f"{fuel.name}: fuel needs no oxygen to burn")

    nu = {fuel.fds_id: -1.0, "OXYGEN": -nu_o2,
          "CARBON DIOXIDE": nu_co2, "WATER VAPOR": nu_h2o}
    if nu_co > 0:
        nu["CARBON MONOXIDE"] = nu_co
    if nu_s > 0:
        nu["SOOT"] = nu_s
    if c > 0:
        nu["HYDROGEN CHLORIDE"] = c
    if n > 0:
        nu["NITROGEN"] = n / 2.0
    return Reaction(fuel, nu)


def blend(components, name="Fuel mix", fds_id="FUEL_MIX"):
    """One effective fuel from [(Fuel, mass fraction), ...].

    The formula is per carbon atom. Heat of combustion and yields are mass
    weighted, so energy and CO/soot per kg of burnt mixture are preserved.
    """
    components = normalise(components)
    if len(components) == 1:
        return components[0][0]

    mol = {e: sum(w * f.atoms(e) / f.molar_mass for f, w in components)
           for e in ELEMENTS}
    formula = {e: v / mol["C"] for e, v in mol.items() if v > 0}
    return Fuel(
        name=name,
        formula=formula,
        heat_of_combustion=sum(w * f.heat_of_combustion for f, w in components),
        co_yield=sum(w * f.co_yield for f, w in components),
        soot_yield=sum(w * f.soot_yield for f, w in components),
        fds_id=fds_id,
    )


def normalise(components):
    """Mass fractions scaled to sum to one, zero entries dropped."""
    components = [(f, float(w)) for f, w in components if w > 0]
    total = sum(w for _, w in components)
    if total <= 0:
        raise ValueError("fuel mixture is empty")
    return [(f, w / total) for f, w in components]
