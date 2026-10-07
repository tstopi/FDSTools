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
    """Q = alpha t^2 up to the peak HRR, then a constant plateau.

    With ``decay_start`` set, the HRR falls linearly from its value at that
    time to zero over ``decay_time`` seconds.
    """

    name = "t-squared"

    def __init__(self, alpha, peak_hrr, duration, n_growth=20,
                 decay_start=None, decay_time=None):
        if alpha <= 0:
            raise ValueError("growth coefficient alpha must be positive")
        if peak_hrr <= 0:
            raise ValueError("peak HRR must be positive")
        if duration <= 0:
            raise ValueError("duration must be positive")
        if n_growth < 2:
            raise ValueError("n_growth must be at least 2")
        if decay_start is not None:
            if decay_start <= 0:
                raise ValueError("decay start must be positive")
            if decay_time is None or decay_time <= 0:
                raise ValueError("decay duration must be positive")
        self.alpha = float(alpha)
        self._peak = float(peak_hrr)
        self.duration = float(duration)
        self.n_growth = int(n_growth)
        self.decay_start = None if decay_start is None else float(decay_start)
        self.decay_time = None if decay_time is None else float(decay_time)

    @property
    def peak_hrr(self):
        return self._peak

    @property
    def t_peak(self):
        """Time (s) at which the growth phase reaches the peak HRR."""
        return math.sqrt(self._peak / self.alpha)

    @property
    def t_end(self):
        """Time (s) the fire is out, or None without a decay phase."""
        if self.decay_start is None:
            return None
        return self.decay_start + self.decay_time

    def _growth(self, t):
        return min(self.alpha * t * t, self._peak)

    def hrr(self, t):
        if t <= 0:
            return 0.0
        if self.decay_start is None or t <= self.decay_start:
            return self._growth(t)
        q0 = self._growth(self.decay_start)
        return max(q0 * (1.0 - (t - self.decay_start) / self.decay_time), 0.0)

    def points(self):
        t_grow = min(self.t_peak, self.duration)
        if self.decay_start is not None:
            t_grow = min(t_grow, self.decay_start)
        pts = []
        for i in range(self.n_growth + 1):
            t = t_grow * i / self.n_growth
            pts.append((t, self.hrr(t) / self._peak))
        knees = [self.duration]
        if self.decay_start is not None:
            knees += [self.decay_start, self.t_end]
        for t in sorted(k for k in set(knees)
                        if t_grow < k <= self.duration):
            pts.append((t, self.hrr(t) / self._peak))
        return pts


CURVE_TYPES = {TSquaredCurve.name: TSquaredCurve}
