import json
import math
import re
import tempfile
import unittest
from pathlib import Path

import design_fire as d
from design_fire import library as lib
from design_fire.curves import curve_from_dict
from design_fire.envelope import fit_power_law, pointwise_max
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

    def test_linear_decay_after_plateau(self):
        c = d.TSquaredCurve(d.GROWTH_RATES["fast"], 2000, 1200,
                            decay_start=600, decay_time=400)
        self.assertEqual(c.hrr(600), 2000)
        self.assertAlmostEqual(c.hrr(800), 1000)
        self.assertEqual(c.hrr(1000), 0)
        self.assertEqual(c.hrr(1100), 0)
        pts = c.points()
        self.assertEqual(pts[-3:], [(600.0, 1.0), (1000.0, 0.0),
                                    (1200.0, 0.0)])
        # ramp reproduces the curve between its points
        for (ta, fa), (tb, fb) in zip(pts[-3:], pts[-2:]):
            tm = (ta + tb) / 2
            self.assertAlmostEqual(c.hrr(tm) / 2000, (fa + fb) / 2)

    def test_decay_before_peak(self):
        c = d.TSquaredCurve(d.GROWTH_RATES["slow"], 5000, 1000,
                            decay_start=300, decay_time=100)
        q0 = c.alpha * 300 ** 2
        self.assertAlmostEqual(c.hrr(350), q0 / 2)
        pts = c.points()
        self.assertAlmostEqual(pts[-3][0], 300)
        self.assertEqual(pts[-2], (400.0, 0.0))
        times = [t for t, _ in pts]
        self.assertEqual(times, sorted(set(times)))

    def test_decay_past_duration(self):
        c = d.TSquaredCurve(d.GROWTH_RATES["fast"], 2000, 700,
                            decay_start=600, decay_time=400)
        pts = c.points()
        self.assertEqual(pts[-1][0], 700)
        self.assertAlmostEqual(pts[-1][1], 0.75)

    def test_registry(self):
        self.assertIs(d.CURVE_TYPES["power-law"], d.PowerLawCurve)
        self.assertIs(d.CURVE_TYPES["tabulated"], d.TabulatedCurve)
        c = d.TSquaredCurve(d.GROWTH_RATES["fast"], 2000, 600)
        self.assertIsInstance(c, d.PowerLawCurve)
        self.assertAlmostEqual(c.growth_time, 150)

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            d.TSquaredCurve(0, 1000, 600)
        with self.assertRaises(ValueError):
            d.TSquaredCurve(0.01, -1, 600)
        with self.assertRaises(ValueError):
            d.TSquaredCurve(0.01, 1000, 600, decay_start=300)


class PowerLaw(unittest.TestCase):
    def test_free_exponent_and_growth_time(self):
        c = d.PowerLawCurve(3, 200, 4000, 900)
        self.assertAlmostEqual(c.hrr(200), 1055)
        self.assertAlmostEqual(c.hrr(100), 1055 / 8)
        self.assertAlmostEqual(c.alpha, 1055 / 200 ** 3)
        self.assertAlmostEqual(c.hrr(c.t_peak), 4000)
        c = d.PowerLawCurve(1.5, 100, 500, 900, q_ref=250)
        self.assertAlmostEqual(c.hrr(100), 250)

    def test_start_time(self):
        c = d.PowerLawCurve(2, 150, 2000, 900, t_start=60, n_growth=10)
        self.assertEqual(c.hrr(60), 0)
        self.assertAlmostEqual(c.hrr(210), 1055)
        pts = c.points()
        self.assertEqual(pts[:2], [(0.0, 0.0), (60.0, 0.0)])
        self.assertAlmostEqual(pts[-2][0], c.t_peak)
        self.assertEqual(len(pts), 13)

    def test_decay_exponent(self):
        c = d.PowerLawCurve(2, 150, 2000, 1200, decay_start=600,
                            decay_time=400, decay_exponent=2)
        self.assertAlmostEqual(c.hrr(800), 500)
        self.assertEqual(c.hrr(1000), 0)
        times = [t for t, _ in c.points()]
        self.assertEqual(len([t for t in times if 600 <= t <= 1000]), 11)

    def test_round_trip(self):
        for c in (d.PowerLawCurve(2.5, 120, 3000, 900, q_ref=1000,
                                  t_start=30, decay_start=500,
                                  decay_time=200, decay_exponent=1.5),
                  d.TabulatedCurve([(0, 0), (60, 100), (120, 0)])):
            c2 = curve_from_dict(json.loads(json.dumps(c.to_dict())))
            self.assertEqual(type(c2), type(c))
            self.assertEqual(c2.points(), c.points())

    def test_bad_input(self):
        for kw in (dict(exponent=0), dict(growth_time=-1), dict(q_ref=0),
                   dict(t_start=-5), dict(t_start=100, decay_start=50,
                                          decay_time=10),
                   dict(decay_start=50, decay_time=10, decay_exponent=0)):
            args = dict(exponent=2, growth_time=150, peak_hrr=1000,
                        duration=600)
            args.update(kw)
            with self.subTest(**kw), self.assertRaises(ValueError):
                d.PowerLawCurve(**args)


class Tabulated(unittest.TestCase):
    def test_interpolation_and_holds(self):
        c = d.TabulatedCurve([(10, 100), (20, 300), (40, 0)])
        self.assertEqual(c.hrr(0), 100)
        self.assertEqual(c.hrr(15), 200)
        self.assertEqual(c.hrr(30), 150)
        self.assertEqual(c.hrr(100), 0)
        self.assertEqual(c.peak_hrr, 300)
        self.assertEqual(c.duration, 40)
        self.assertEqual(c.points()[1], (20.0, 1.0))

    def test_bad_input(self):
        for table in ([(0, 1)], [(0, 0), (0, 1)], [(0, 0), (10, -1)],
                      [(0, 0), (10, 0)], [(-1, 0), (10, 5)]):
            with self.subTest(table=table), self.assertRaises(ValueError):
                d.TabulatedCurve(table)

    def test_writer(self):
        c = d.TabulatedCurve([(0, 0), (60, 500), (300, 1000), (600, 0)])
        text = d.DesignFire(c, 2, [(d.FUELS["Propane"], 1)]).to_fds()
        self.assertIn("HRRPUA=500.", text)
        self.assertIn("&RAMP ID='FIRE_RAMP', T=60., F=0.5 /", text)
        self.assertIn("! tabulated HRR, 4 points to 600 s", text)


def grid(curves, n=2000):
    t_end = max(c.duration for c in curves) * 1.1
    return [t_end * i / n for i in range(n + 1)]


def ramp(curve, t):
    """HRR at ``t`` as FDS sees it: the &RAMP points joined linearly."""
    from design_fire.curves import interpolate
    return interpolate([(a, f * curve.peak_hrr) for a, f in curve.points()],
                       t)


class Envelope(unittest.TestCase):
    def test_pointwise_max_with_crossing(self):
        a = d.TabulatedCurve([(0, 0), (100, 1000), (200, 0)])
        b = d.TabulatedCurve([(0, 0), (50, 200), (300, 700)])
        env = pointwise_max([a, b])
        times = [t for t, _ in env.table]
        # crossing of a's fall and b's rise
        self.assertTrue(any(150 < t < 200 for t in times))
        self.assertNotIn(250.0, times)    # collinear point dropped
        for t in grid([a, b]):
            self.assertAlmostEqual(env.hrr(t), max(a.hrr(t), b.hrr(t)),
                                   places=6)

    def test_pointwise_max_of_power_laws(self):
        cs = [d.PowerLawCurve(2, 300, 3000, 1200, decay_start=700,
                              decay_time=300),
              d.PowerLawCurve(3, 100, 1500, 900, t_start=40)]
        env = pointwise_max(cs)
        self.assertEqual(env.duration, 1200)
        for t in grid(cs):
            self.assertAlmostEqual(env.hrr(t), max(ramp(c, t) for c in cs),
                                   places=6)

    def test_fit_same_exponent(self):
        cs = [d.PowerLawCurve(2, g, p, 1200) for g, p in
              ((150, 2000), (300, 5000), (600, 3000))]
        fit = fit_power_law(cs)
        self.assertAlmostEqual(fit.growth_time, 150)
        self.assertEqual(fit.peak_hrr, 5000)
        self.assertIsNone(fit.decay_start)
        self.assertEqual(fit.t_start, 0)

    def test_fit_covers_with_decay(self):
        cs = [d.TabulatedCurve([(0, 0), (60, 50), (180, 800), (360, 2400),
                                (480, 2000), (720, 600), (1080, 0)]),
              d.PowerLawCurve(2, 300, 1500, 1000, decay_start=600,
                              decay_time=300, decay_exponent=2),
              d.PowerLawCurve(3, 200, 1200, 800, t_start=30,
                              decay_start=500, decay_time=200)]
        for m in (1, 2, 0.5):
            with self.subTest(decay_exponent=m):
                fit = fit_power_law(cs, exponent=2, decay_exponent=m)
                self.assertEqual(fit.peak_hrr, 2400)
                self.assertEqual(fit.decay_start, 360)
                self.assertEqual(fit.duration, 1080)
                for t in grid(cs):
                    if t < 60:      # straight rise from zero, see docstring
                        continue
                    for c in cs:
                        self.assertGreaterEqual(fit.hrr(t) + 1e-6, c.hrr(t),
                                                msg=(t, c.name))

    def test_fit_tight(self):
        # the decay touches the tabulated curve at 480 s
        c = d.TabulatedCurve([(0, 0), (360, 2400), (480, 2000), (1080, 0)])
        fit = fit_power_law([c])
        self.assertAlmostEqual(fit.hrr(480), 2000)
        self.assertAlmostEqual(fit.t_end, 1080)

    def test_fit_delayed_start(self):
        cs = [d.PowerLawCurve(2, 150, 2000, 900, t_start=120),
              d.PowerLawCurve(2, 300, 2000, 900, t_start=200)]
        fit = fit_power_law(cs)
        self.assertEqual(fit.t_start, 120)
        self.assertAlmostEqual(fit.growth_time, 150)

    def test_fit_errors(self):
        with self.assertRaises(ValueError):
            fit_power_law([])
        with self.assertRaises(ValueError):
            fit_power_law([d.TabulatedCurve([(0, 100), (60, 200)])])


class Library(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_builtins_load(self):
        entries, errors = lib.load_library(user=self.dir)
        self.assertEqual(errors, [])
        self.assertGreaterEqual(len(entries), 5)
        self.assertTrue(all(e.builtin for e in entries))
        for e in entries:
            d.DesignFire(e.curve, e.area or 1, e.components()).to_fds()

    def test_save_load_delete(self):
        curve = d.PowerLawCurve(2.5, 200, 3000, 900, t_start=10)
        e = lib.LibraryEntry("Sofa / test", curve, [("PMMA", 0.3),
                                                    ("Nylon", 0.7)], 2.5,
                             "note")
        path = lib.save_entry(e, self.dir)
        self.assertEqual(path.name, "sofa_test.json")
        entries, errors = lib.load_library(self.dir / "none", self.dir)
        self.assertEqual(errors, [])
        (e2,) = entries
        self.assertEqual((e2.name, e2.fuels, e2.area, e2.description),
                         (e.name, e.fuels, e.area, e.description))
        self.assertEqual(e2.curve.points(), curve.points())
        self.assertFalse(e2.builtin)
        lib.delete_entry(e2)
        self.assertFalse(path.exists())

    def test_bad_files_reported(self):
        (self.dir / "broken.json").write_text("{")
        (self.dir / "fuel.json").write_text(json.dumps(
            {"name": "x", "curve": {"type": "tabulated",
                                    "table": [[0, 0], [1, 1]]},
             "fuels": [["Unobtainium", 1]]}))
        entries, errors = lib.load_library(self.dir / "none", self.dir)
        self.assertEqual(entries, [])
        self.assertEqual(len(errors), 2)

    def test_builtin_cannot_be_deleted(self):
        entries, _ = lib.load_library(user=self.dir)
        with self.assertRaises(ValueError):
            lib.delete_entry(entries[0])

    def test_read_csv(self):
        for text in ("time (s),HRR (kW)\n0,0\n60,500\n120,1000\n",
                     "# test\nt;Q\n0;0\n60;500\n120;1000\n",
                     "0\t0\n60\t500\n120\t1000\n",
                     "0 0\n60  500\n120 1000\n"):
            with self.subTest(text=text):
                p = self.dir / "hrr.csv"
                p.write_text(text)
                c = lib.read_csv_curve(p)
                self.assertEqual(c.table, [(0, 0), (60, 500), (120, 1000)])


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
