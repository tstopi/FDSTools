"""Heat release rate curves.

Every curve type subclasses ``Curve`` and is registered in ``CURVE_TYPES`` so
the GUI, writer and library can offer new types without changes elsewhere.
"""

import math

Q_REF = 1055.0   # kW, reference HRR for growth times (NFPA 72 / SFPE)

# Time (s) for a t^2 fire to reach 1055 kW (NFPA 72 / SFPE Handbook).
GROWTH_TIMES = {"slow": 600.0, "medium": 300.0, "fast": 150.0,
                "ultrafast": 75.0}
# Matching t^2 growth coefficients (kW/s^2).
GROWTH_RATES = {k: Q_REF / tg ** 2 for k, tg in GROWTH_TIMES.items()}


class Curve:
    """HRR as a function of time; ``points`` gives the &RAMP table."""

    name = ""
    duration = 0.0

    def hrr(self, t):
        raise NotImplementedError

    @property
    def peak_hrr(self):
        raise NotImplementedError

    def points(self):
        """List of (t [s], fraction of peak HRR) for an FDS &RAMP."""
        raise NotImplementedError

    def describe(self):
        """Comment lines for the FDS header, without the leading '!'."""
        return []

    def to_dict(self):
        raise NotImplementedError

    @classmethod
    def from_dict(cls, data):
        raise NotImplementedError


class PowerLawCurve(Curve):
    """Q = Q_ref ((t - t_start) / t_g)^n up to the peak HRR, then a plateau.

    ``t_g`` is the time the growth takes to reach ``q_ref`` after
    ``t_start``; with n = 2 and Q_ref = 1055 kW this is the usual t^2 fire.
    With ``decay_start`` set, the HRR falls from its value q0 at that time as
    q0 (1 - (t - decay_start) / decay_time)^m, reaching zero after
    ``decay_time`` seconds; m = 1 is a linear decay.
    """

    name = "power-law"

    def __init__(self, exponent, growth_time, peak_hrr, duration,
                 q_ref=Q_REF, t_start=0.0, decay_start=None, decay_time=None,
                 decay_exponent=1.0, n_growth=20, n_decay=10):
        if exponent <= 0:
            raise ValueError("growth exponent must be positive")
        if growth_time <= 0:
            raise ValueError("growth time must be positive")
        if q_ref <= 0:
            raise ValueError("reference HRR must be positive")
        if peak_hrr <= 0:
            raise ValueError("peak HRR must be positive")
        if duration <= 0:
            raise ValueError("duration must be positive")
        if t_start < 0:
            raise ValueError("start time must not be negative")
        if n_growth < 2:
            raise ValueError("n_growth must be at least 2")
        if decay_start is not None:
            if decay_start <= t_start:
                raise ValueError("decay must start after the fire starts")
            if decay_time is None or decay_time <= 0:
                raise ValueError("decay duration must be positive")
            if decay_exponent <= 0:
                raise ValueError("decay exponent must be positive")
        self.exponent = float(exponent)
        self.growth_time = float(growth_time)
        self._peak = float(peak_hrr)
        self.duration = float(duration)
        self.q_ref = float(q_ref)
        self.t_start = float(t_start)
        self.decay_start = None if decay_start is None else float(decay_start)
        self.decay_time = None if decay_time is None else float(decay_time)
        self.decay_exponent = float(decay_exponent)
        self.n_growth = int(n_growth)
        self.n_decay = int(n_decay)

    @property
    def peak_hrr(self):
        return self._peak

    @property
    def alpha(self):
        """Growth coefficient Q_ref / t_g^n (kW/s^n)."""
        return self.q_ref / self.growth_time ** self.exponent

    @property
    def t_peak(self):
        """Time (s) at which the growth phase reaches the peak HRR."""
        return self.t_start + self.growth_time * (
            self._peak / self.q_ref) ** (1.0 / self.exponent)

    @property
    def t_end(self):
        """Time (s) the fire is out, or None without a decay phase."""
        if self.decay_start is None:
            return None
        return self.decay_start + self.decay_time

    def _growth(self, t):
        if t <= self.t_start:
            return 0.0
        s = (t - self.t_start) / self.growth_time
        return min(self.q_ref * s ** self.exponent, self._peak)

    def hrr(self, t):
        if self.decay_start is None or t <= self.decay_start:
            return self._growth(t)
        q0 = self._growth(self.decay_start)
        f = max(1.0 - (t - self.decay_start) / self.decay_time, 0.0)
        return q0 * f ** self.decay_exponent

    def points(self):
        t_grow = min(self.t_peak, self.duration)
        if self.decay_start is not None:
            t_grow = min(t_grow, self.decay_start)
        ts = min(self.t_start, t_grow)
        times = [0.0, ts, t_grow]
        n = self.n_growth
        times += [ts + (t_grow - ts) * i / n for i in range(1, n)]
        if self.decay_start is not None:
            if self.decay_exponent == 1.0:
                times += [self.decay_start, self.t_end]
            else:
                k = self.n_decay
                times += [self.decay_start + self.decay_time * i / k
                          for i in range(k + 1)]
        times.append(self.duration)
        times = sorted(set(t for t in times if 0.0 <= t <= self.duration))
        return [(t, self.hrr(t) / self._peak) for t in times]

    def describe(self):
        n = self.exponent
        out = [f"Q = {self.q_ref:g} kW * ((t - {self.t_start:g} s) / "
               f"{self.growth_time:.4g} s)^{n:g}, alpha {self.alpha:.5g} "
               f"kW/s^{n:g}, peak reached at {self.t_peak:.0f} s"]
        if self.decay_start is not None:
            kind = ("linear" if self.decay_exponent == 1.0
                    else f"exponent {self.decay_exponent:g}")
            out.append(f"{kind} decay from {self.decay_start:g} s to zero at "
                       f"{self.t_end:g} s")
        return out

    def to_dict(self):
        d = {"type": PowerLawCurve.name, "exponent": self.exponent,
             "growth_time": self.growth_time, "peak_hrr": self._peak,
             "duration": self.duration, "q_ref": self.q_ref,
             "t_start": self.t_start}
        if self.decay_start is not None:
            d.update(decay_start=self.decay_start, decay_time=self.decay_time,
                     decay_exponent=self.decay_exponent)
        return d

    @classmethod
    def from_dict(cls, data):
        keys = ("exponent", "growth_time", "peak_hrr", "duration", "q_ref",
                "t_start", "decay_start", "decay_time", "decay_exponent")
        return PowerLawCurve(**{k: data[k] for k in keys if k in data})


class TSquaredCurve(PowerLawCurve):
    """Q = alpha t^2 up to the peak HRR, then a constant plateau.

    A power-law curve with n = 2 given by its growth coefficient. With
    ``decay_start`` set, the HRR falls linearly from its value at that time
    to zero over ``decay_time`` seconds.
    """

    name = "t-squared"

    def __init__(self, alpha, peak_hrr, duration, n_growth=20,
                 decay_start=None, decay_time=None):
        if alpha <= 0:
            raise ValueError("growth coefficient alpha must be positive")
        super().__init__(2.0, math.sqrt(Q_REF / alpha), peak_hrr, duration,
                         decay_start=decay_start, decay_time=decay_time,
                         n_growth=n_growth)


class TabulatedCurve(Curve):
    """HRR given as (t [s], Q [kW]) pairs, linearly interpolated.

    Like an FDS &RAMP, the first value holds before the first time and the
    last value after the last time. The duration is the last time.
    """

    name = "tabulated"

    def __init__(self, table):
        table = [(float(t), float(q)) for t, q in table]
        if len(table) < 2:
            raise ValueError("a tabulated curve needs at least two points")
        times = [t for t, _ in table]
        if times[0] < 0:
            raise ValueError("times must not be negative")
        if any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("times must be strictly increasing")
        if any(q < 0 for _, q in table):
            raise ValueError("HRR values must not be negative")
        if max(q for _, q in table) <= 0:
            raise ValueError("peak HRR must be positive")
        self.table = table
        self.duration = times[-1]

    @property
    def peak_hrr(self):
        return max(q for _, q in self.table)

    def hrr(self, t):
        return interpolate(self.table, t)

    def points(self):
        p = self.peak_hrr
        return [(t, q / p) for t, q in self.table]

    def describe(self):
        return [f"tabulated HRR, {len(self.table)} points to "
                f"{self.duration:g} s"]

    def to_dict(self):
        return {"type": TabulatedCurve.name,
                "table": [[t, q] for t, q in self.table]}

    @classmethod
    def from_dict(cls, data):
        return TabulatedCurve(data["table"])


def interpolate(table, t):
    """Linear interpolation in sorted (t, value) pairs, holding the ends."""
    if t <= table[0][0]:
        return table[0][1]
    if t >= table[-1][0]:
        return table[-1][1]
    lo, hi = 0, len(table) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if table[mid][0] <= t:
            lo = mid
        else:
            hi = mid
    (ta, qa), (tb, qb) = table[lo], table[hi]
    return qa + (qb - qa) * (t - ta) / (tb - ta)


def curve_from_dict(data):
    try:
        cls = CURVE_TYPES[data["type"]]
    except KeyError:
        raise ValueError(f"unknown curve type {data.get('type')!r}") from None
    return cls.from_dict(data)


CURVE_TYPES = {c.name: c for c in (PowerLawCurve, TabulatedCurve)}
