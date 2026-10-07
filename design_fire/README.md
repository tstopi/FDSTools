# design_fire

A small tkinter app that turns a few choices into FDS input lines for a design
fire: the heat release rate curve (`&RAMP` + `&SURF`) and the combustion
chemistry (`&SPEC` + `&REAC`). Standard library only (needs Python's tkinter).

## Run

```
python -m design_fire          # from the repository root
fds-design-fire                # after `pip install -e .`
```

On the **Design fire** tab, pick the curve, fire area and duration, build the
fuel mixture, then copy the text or save it to a `.fds` file to paste into
your model. Put `SURF_ID='FIRE'` on the burner `&VENT` or `&OBST` yourself.

The **Library** tab lists saved design fires. Tick several to overlay them
and get a design fire that envelopes them all, then send it to the Design
fire tab or save it back to the library.

## What it writes

- **Curve.** Two types. `HRRPUA` is the peak HRR divided by the fire area
  and `RAMP_Q` scales it in time.
  - *Power law:* `Q = Q_ref ((t − t_start) / t_g)^n` to the peak HRR, then a
    plateau until the duration. The exponent `n`, growth time `t_g` (time to
    reach `Q_ref`), `Q_ref` (default 1055 kW) and start time are free. The
    slow, medium, fast and ultrafast presets are the NFPA/SFPE t² fires
    (n = 2, 1055 kW at 600, 300, 150, 75 s). An optional decay starts at a
    set time and takes the HRR from its value there to zero over the decay
    duration as `(1 − τ)^m`; m = 1 is linear.
  - *Tabulated:* time (s) and HRR (kW) pairs loaded from a CSV or text
    file, written to the `&RAMP` as they are.
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

## Library

Each library fire is one JSON file: name, description, curve, and
optionally the fuel mixture and fire area.

```json
{"name": "Sofa", "description": "...",
 "curve": {"type": "tabulated", "table": [[0, 0], [60, 120], [300, 1800]]},
 "fuels": [["Polyurethane foam (flexible)", 1.0]], "area": 2.0}
```

Built-in fires ship in `design_fire/library/`: the four t² growth rates to
5 MW and an illustrative tabulated example (a made-up shape that shows the
format, not test data). Your own fires are saved in
`~/.fdstools/design_fires/`, or in the folder named by the
`FDSTOOLS_FIRE_LIBRARY` environment variable. *Save design…* stores the
current Design fire tab, *Import CSV…* stores a tabulated curve, and *Load
into design* opens an entry for editing (save it again under the same name
to update it). Built-in entries cannot be deleted.

## Enveloping design fire

Tick library fires to stack them. Two methods:

- **Pointwise maximum.** At every time, the highest of the ticked fires.
  It is exact for the curves as FDS runs them (their `&RAMP` points joined
  by straight lines, last value held), with a point at every input point
  and every crossing. The result is a tabulated curve.
- **Fitted power law.** For the chosen exponent `n` and `Q_ref`, the
  slowest power-law curve (longest `t_g`) that stays at or above every
  ticked fire, with a plateau at the largest peak. It starts when the first
  fire starts. When every fire has burnt out by the end, it adds the
  shortest decay with the chosen decay exponent that still covers the
  tails; otherwise the plateau runs to the end. A straight rise from zero
  in a tabulated curve cannot be covered by an n > 1 power law starting at
  the same moment, so growth is checked at the curves' own points. Add
  points to the table if that first segment matters.

*Reaction from* picks the fuel mixture sent with the envelope: one of the
ticked fires' mixtures, or whatever the Design fire tab already has. The
fire area is left as it is on the Design fire tab.

## Extending

New curve types subclass `curves.Curve` (implement `hrr`, `peak_hrr`,
`points`, `describe`, `to_dict` and `from_dict`) and register in
`CURVE_TYPES`. New fuels go in `fuels.FUELS`.

## Tests

```
python -m unittest discover -s design_fire/tests -t .
```
