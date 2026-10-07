"""Heat release rate curves.

Every curve type subclasses ``Curve`` and is registered in ``CURVE_TYPES`` so
the GUI and writer can offer new types without changes elsewhere.
"""

import math

# t^2 growth coefficients (kW/s^2): time to reach 1055 kW is 600, 300, 150
# and 75 s (NFPA 72 / SFPE Handbook).
GROWTH_RATES = {
    "slow": 1055.0 / 600.0 ** 2,
    "medium": 1055.0 / 300.0 ** 2,
    "fast": 1055.0 / 150.0 ** 2,
    "ultrafast": 1055.0 / 75.0 ** 2,
}


class Curve:
    """HRR as a function of time; ``points`` gives the &RAMP table."""

    name = ""

    def hrr(self, t):
        raise NotImplementedError

    @property
    def peak_hrr(self):
        raise NotImplementedError

    def points(self):
        """List of (t [s], fraction of peak HRR) for an FDS &RAMP."""
        raise NotImplementedError


class TSquaredCurve(Curve):
    """Q = alpha t^2 up to the peak HRR, then a constant plateau."""

    name = "t-squared"

    def __init__(self, alpha, peak_hrr, duration, n_growth=20):
        if alpha <= 0:
            raise ValueError("growth coefficient alpha must be positive")
        if peak_hrr <= 0:
            raise ValueError("peak HRR must be positive")
        if duration <= 0:
            raise ValueError("duration must be positive")
        if n_growth < 2:
            raise ValueError("n_growth must be at least 2")
        self.alpha = float(alpha)
        self._peak = float(peak_hrr)
        self.duration = float(duration)
        self.n_growth = int(n_growth)

    @property
    def peak_hrr(self):
        return self._peak

    @property
    def t_peak(self):
        """Time (s) at which the growth phase reaches the peak HRR."""
        return math.sqrt(self._peak / self.alpha)

    def hrr(self, t):
        if t <= 0:
            return 0.0
        return min(self.alpha * t * t, self._peak)

    def points(self):
        t_grow = min(self.t_peak, self.duration)
        pts = []
        for i in range(self.n_growth + 1):
            t = t_grow * i / self.n_growth
            pts.append((t, self.hrr(t) / self._peak))
        if self.duration > t_grow:
            pts.append((self.duration, self.hrr(self.duration) / self._peak))
        return pts


CURVE_TYPES = {TSquaredCurve.name: TSquaredCurve}
