import math
import re
import unittest

import design_fire as d
from design_fire.fuels import formula_string
from design_fire.reactions import SPECIES_FORMULA, molar_mass


def parse_nu(text):
    """[(species list, nu list)] from each &REAC in the FDS text."""
    out = []
    for m in re.finditer(r"SPEC_ID_NU=(.*?),\s*NU=(.*?) /", text, re.S):
        sp = re.findall(r"'([^']*)'", m.group(1))
        nu = [float(v) for v in m.group(2).split(",")]
        out.append((sp, nu))
    return out


class Curves(unittest.TestCase):
    def test_growth_rates_reach_1055_kw(self):
        for name, tg in [("slow", 600), ("medium", 300), ("fast", 150),
                         ("ultrafast", 75)]:
            c = d.TSquaredCurve(d.GROWTH_RATES[name], 5000, 1000)
            self.assertAlmostEqual(c.hrr(tg), 1055.0, places=6)

    def test_plateau_and_ramp_points(self):
        c = d.TSquaredCurve(d.GROWTH_RATES["fast"], 2000, 600, n_growth=10)
        self.assertAlmostEqual(c.t_peak, 150 * math.sqrt(2000 / 1055))
        self.assertEqual(c.hrr(400), 2000)
        pts = c.points()
        self.assertEqual(pts[0], (0.0, 0.0))
        self.assertAlmostEqual(pts[-2][0], c.t_peak)
        self.assertAlmostEqual(pts[-2][1], 1.0)
        self.assertEqual(pts[-1], (600.0, 1.0))
        self.assertEqual(len(pts), 12)
        times = [t for t, _ in pts]
        self.assertEqual(times, sorted(times))

    def test_duration_shorter_than_growth(self):
        c = d.TSquaredCurve(d.GROWTH_RATES["slow"], 5000, 300)
        pts = c.points()
        self.assertAlmostEqual(pts[-1][0], 300)
        self.assertLess(pts[-1][1], 1.0)

    def test_registry(self):
        self.assertIs(d.CURVE_TYPES["t-squared"], d.TSquaredCurve)

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            d.TSquaredCurve(0, 1000, 600)
        with self.assertRaises(ValueError):
            d.TSquaredCurve(0.01, -1, 600)


class Reactions(unittest.TestCase):
    def assertBalanced(self, r):
        for e, v in r.element_imbalance().items():
            self.assertAlmostEqual(v, 0.0, places=9, msg=e)
        self.assertAlmostEqual(r.mass_imbalance(), 0.0, places=9)

    def test_all_library_fuels_balance(self):
        for f in d.FUELS.values():
            with self.subTest(fuel=f.name):
                self.assertBalanced(d.balance(f))

    def test_pvc_products(self):
        f = d.FUELS["PVC"]
        r = d.balance(f)
        self.assertEqual(r.nu["HYDROGEN CHLORIDE"], 1.0)
        self.assertAlmostEqual(r.nu["WATER VAPOR"], 1.0)
        # yields come back out of the stoichiometry
        m_f = f.molar_mass
        self.assertAlmostEqual(
            r.nu["SOOT"] * molar_mass(SPECIES_FORMULA["SOOT"]) / m_f,
            f.soot_yield)
        self.assertAlmostEqual(
            r.nu["CARBON MONOXIDE"]
            * molar_mass(SPECIES_FORMULA["CARBON MONOXIDE"]) / m_f,
            f.co_yield)
        # HCl mass yield = Cl mass fraction * M_HCl / M_Cl, about 0.58
        y_hcl = molar_mass(SPECIES_FORMULA["HYDROGEN CHLORIDE"]) / m_f
        self.assertAlmostEqual(y_hcl, 0.5834, places=3)

    def test_propane_complete(self):
        r = d.balance(d.Fuel("p", {"C": 3, "H": 8}, 46000, 0, 0))
        self.assertAlmostEqual(r.nu["OXYGEN"], -5)
        self.assertAlmostEqual(r.nu["CARBON DIOXIDE"], 3)
        self.assertAlmostEqual(r.nu["WATER VAPOR"], 4)
        self.assertNotIn("SOOT", r.nu)

    def test_impossible_yields(self):
        with self.assertRaises(ValueError):
            d.balance(d.Fuel("x", {"C": 1, "H": 4}, 50000, 0.0, 0.9))
        with self.assertRaises(ValueError):
            d.balance(d.Fuel("x", {"C": 1, "H": 1, "Cl": 3}, 5000, 0, 0))

    def test_blend_preserves_element_mass(self):
        pvc, pe = d.FUELS["PVC"], d.FUELS["Polyethylene"]
        mix = d.blend([(pvc, 3), (pe, 7)])   # weights are normalised
        # Cl mass fraction of the blend = 0.3 * Cl fraction of PVC
        cl = lambda f: f.atoms("Cl") * 35.45 / f.molar_mass
        self.assertAlmostEqual(cl(mix), 0.3 * cl(pvc))
        self.assertAlmostEqual(mix.heat_of_combustion,
                               0.3 * 5700 + 0.7 * 38400)
        self.assertAlmostEqual(mix.soot_yield, 0.3 * 0.172 + 0.7 * 0.060)
        self.assertEqual(mix.atoms("C"), 1.0)
        self.assertBalanced(d.balance(mix))

    def test_blend_single_returns_fuel(self):
        pvc = d.FUELS["PVC"]
        self.assertIs(d.blend([(pvc, 1), (d.FUELS["PMMA"], 0)]), pvc)

    def test_formula_string(self):
        self.assertEqual(formula_string({"C": 2, "H": 3, "Cl": 1}), "C2H3Cl")
        self.assertEqual(formula_string({"C": 1, "H": 1.75, "O": 0.32}),
                         "CH1.75O0.32")


class Writer(unittest.TestCase):
    def setUp(self):
        self.curve = d.TSquaredCurve(d.GROWTH_RATES["fast"], 2000, 600)
        self.mix = [(d.FUELS["PVC"], 0.3), (d.FUELS["Polyethylene"], 0.7)]

    def check_atoms(self, text):
        """Every &REAC in the text conserves atoms using the written values."""
        formulas = dict(SPECIES_FORMULA)
        for m in re.finditer(r"&SPEC ID='([^']*)', FORMULA='([^']*)'", text):
            formulas[m.group(1)] = {
                e: float(n) if n else 1.0
                for e, n in re.findall(r"([A-Z][a-z]?)([0-9.]*)", m.group(2))}
        for f in d.FUELS.values():
            formulas.setdefault(f.fds_id, f.formula)
        for sp, nu in parse_nu(text):
            for e in ("C", "H", "O", "N", "Cl"):
                bal = sum(n * formulas[s].get(e, 0) for s, n in zip(sp, nu))
                self.assertAlmostEqual(bal, 0, places=3, msg=(sp[0], e))

    def test_blended(self):
        text = d.DesignFire(self.curve, 4, self.mix).to_fds()
        self.assertIn("HRRPUA=500.", text)
        self.assertIn("RAMP_Q='FIRE_RAMP'", text)
        self.assertEqual(text.count("&REAC"), 1)
        self.assertIn("&SPEC ID='HYDROGEN CHLORIDE' /", text)
        self.assertIn("&SPEC ID='NITROGEN', BACKGROUND=.TRUE. /", text)
        self.assertEqual(text.count("&RAMP"), len(self.curve.points()))
        self.check_atoms(text)

    def test_per_fuel(self):
        fire = d.DesignFire(self.curve, 4, self.mix, per_fuel_reactions=True)
        text = fire.to_fds()
        self.assertEqual(text.count("&REAC"), 2)
        self.assertNotIn("HRRPUA=", text.split("&SURF")[1])
        self.check_atoms(text)
        flux = [float(v) for v in
                re.search(r"MASS_FLUX\(1:2\)=([^\n]*),\n", text).group(1)
                .split(",")]
        hrr = flux[0] * 5700 + flux[1] * 38400
        self.assertAlmostEqual(hrr, 500.0, places=3)
        self.assertAlmostEqual(flux[0] / flux[1], 0.3 / 0.7, places=4)

    def test_predefined_fuel_has_no_formula(self):
        text = d.DesignFire(self.curve, 1, [(d.FUELS["Propane"], 1)]).to_fds()
        self.assertIn("&SPEC ID='PROPANE' /", text)
        self.assertNotIn("HYDROGEN CHLORIDE", text)

    def test_duplicate_fuel_rejected(self):
        with self.assertRaises(ValueError):
            d.DesignFire(self.curve, 1, [(d.FUELS["PVC"], 1),
                                         (d.FUELS["PVC"], 1)])


if __name__ == "__main__":
    unittest.main()
