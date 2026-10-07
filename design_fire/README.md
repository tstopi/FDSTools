# design_fire

A small tkinter app that turns a few choices into FDS input lines for a design
fire: the heat release rate curve (`&RAMP` + `&SURF`) and the combustion
chemistry (`&SPEC` + `&REAC`). Standard library only (needs Python's tkinter).

## Run

```
python -m design_fire          # from the repository root
fds-design-fire                # after `pip install -e .`
```

Pick the growth rate, peak HRR, fire area and duration, build the fuel
mixture, then copy the text or save it to a `.fds` file to paste into your
model. Put `SURF_ID='FIRE'` on the burner `&VENT` or `&OBST` yourself.

## What it writes

- **Curve.** t² growth (`Q = αt²`) to the peak HRR, then a plateau until the
  duration. Growth rates are the NFPA/SFPE slow, medium, fast and ultrafast
  values (1055 kW at 600, 300, 150, 75 s) or a custom α. `HRRPUA` is the peak
  HRR divided by the fire area and `RAMP_Q` scales it in time.
- **Fuels.** Formula, effective heat of combustion and CO/soot yields from the
  SFPE Handbook (Tewarson's tables), in `fuels.py`. Check them against the
  edition you cite before relying on them.
- **Reaction.** Each fuel `CxHyOzNnClc` is balanced explicitly: CO and soot
  (pure carbon) from the yields, all chlorine to HCl, remaining hydrogen to
  water, remaining carbon to CO2, nitrogen to N2. This needs tracked species
  rather than FDS simple chemistry, so the output defines `&SPEC` lines with
  nitrogen as the background species and air O2/CO2 as initial mass
  fractions.
- **Mixtures.** By default the fuels are blended by mass fraction into one
  effective fuel (`FUEL_MIX`, formula per carbon atom, mass-weighted heat of
  combustion and yields) with a single `&REAC`. With *One &REAC per fuel*
  each fuel gets its own species and reaction, and the `&SURF` uses
  `MASS_FLUX` per fuel, split by mass fraction and sized so the total heat
  release still follows the ramp.

## Extending

New curve types subclass `curves.Curve` (implement `hrr`, `peak_hrr`,
`points`) and register in `CURVE_TYPES`. New fuels go in `fuels.FUELS`.

## Tests

```
python -m unittest discover -s design_fire/tests -t .
```
