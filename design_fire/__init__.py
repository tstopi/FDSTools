"""Design fire generator for FDS: HRR curve plus combustion reaction."""

from .curves import (CURVE_TYPES, GROWTH_RATES, GROWTH_TIMES, Curve,
                     PowerLawCurve, TabulatedCurve, TSquaredCurve)
from .envelope import fit_power_law, pointwise_max
from .fuels import FUELS, Fuel
from .library import LibraryEntry, load_library
from .reactions import Reaction, balance, blend
from .writer import DesignFire

__all__ = ["CURVE_TYPES", "GROWTH_RATES", "GROWTH_TIMES", "Curve",
           "PowerLawCurve", "TabulatedCurve", "TSquaredCurve",
           "fit_power_law", "pointwise_max", "FUELS", "Fuel", "LibraryEntry",
           "load_library", "Reaction", "balance", "blend", "DesignFire"]
