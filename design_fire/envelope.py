"""Design fires that envelope several others.

``pointwise_max`` is exact for the curves as FDS sees them: each curve's
&RAMP points joined by straight lines and held past the last point. The
result is a tabulated curve with a point at every knot of every input curve
and at every crossing between them.

``fit_power_law`` finds the slowest power-law curve with a chosen exponent
that stays at or above all the curves: the longest growth time t_g, a
plateau at the largest peak, and (when every curve has burnt out) the
shortest decay that still covers the tails. See its docstring for where
coverage is checked.
"""

from .curves import Q_REF, PowerLawCurve, TabulatedCurve, interpolate

REL_TOL = 1e-9
DECAY_SUBSTEPS = 32


def _knots(curve):
    p = curve.peak_hrr
    return [(t, f * p) for t, f in curve.points()]


def _drop_collinear(table):
    out = [table[0]]
    for i in range(1, len(table) - 1):
        (ta, qa), (tb, qb), (tc, qc) = out[-1], table[i], table[i + 1]
        q_line = qa + (qc - qa) * (tb - ta) / (tc - ta)
        scale = max(abs(qa), abs(qb), abs(qc), 1.0)
        if abs(q_line - qb) > REL_TOL * scale:
            out.append(table[i])
    out.append(table[-1])
    return out


def pointwise_max(curves):
    """Tabulated curve equal to the highest of the curves at every time."""
    curves = list(curves)
    if not curves:
        raise ValueError("select at least one design fire")
    knots = [_knots(c) for c in curves]
    times = sorted({0.0} | {t for k in knots for t, _ in k})
    extra = []
    for a, b in zip(times, times[1:]):
        va = [interpolate(k, a) for k in knots]
        vb = [interpolate(k, b) for k in knots]
        eps = REL_TOL * max(b, 1.0)
        for i in range(len(knots)):
            for j in range(i + 1, len(knots)):
                da, db = va[i] - va[j], vb[i] - vb[j]
                if da * db < 0:
                    tc = a + (b - a) * da / (da - db)
                    if a + eps < tc < b - eps:
                        extra.append(tc)
    times = sorted(set(times) | set(extra))
    table = [(t, max(interpolate(k, t) for k in knots)) for t in times]
    return TabulatedCurve(_drop_collinear(table))


def _exact_max(curves, t):
    return max(c.hrr(t) for c in curves)


def _segment_critical(table, t0, n):
    """On each straight segment Q = c + b s (s = t - t0), Q / s^n peaks at
    s = -n c / (b (n - 1)); return those times that fall inside a segment."""
    out = []
    if n == 1.0:
        return out
    for (ta, qa), (tb, qb) in zip(table, table[1:]):
        sa, sb = ta - t0, tb - t0
        if sb <= 0 or tb == ta:
            continue
        b = (qb - qa) / (tb - ta)
        if b == 0:
            continue
        c = qa - b * sa
        s = -n * c / (b * (n - 1))
        if max(sa, 0.0) < s < sb:
            out.append(t0 + s)
    return out


def fit_power_law(curves, exponent=2.0, q_ref=Q_REF, decay_exponent=1.0):
    """Slowest ``PowerLawCurve`` with the given exponent covering all curves.

    The growth is checked against each curve at its own &RAMP points and,
    on straight tabulated segments, where a power law is hardest pressed.
    A straight rise from zero cannot be covered by any n > 1 power law
    starting at the same time, so between a zero and the next point of a
    tabulated curve the fit can fall below the straight line joining them;
    add points there if that matters. The plateau and decay are checked
    against the highest curve at every point and between them.
    """
    curves = list(curves)
    if not curves:
        raise ValueError("select at least one design fire")
    if exponent <= 0:
        raise ValueError("growth exponent must be positive")
    if decay_exponent <= 0:
        raise ValueError("decay exponent must be positive")
    knot_times = sorted({0.0} | {t for c in curves for t, _ in c.points()})
    duration = max(c.duration for c in curves)
    peak = max(_exact_max(curves, t) for t in knot_times)

    # start: the last knot before which every curve is still at zero
    if _exact_max(curves, 0.0) > 0:
        raise ValueError("a curve already burns at t = 0; a power-law curve "
                         "starting from zero cannot cover it")
    t0 = 0.0
    for t in knot_times:
        if _exact_max(curves, t) > 0:
            break
        t0 = t

    # growth: largest Q/s^n over the check points sets the growth time
    worst = 0.0
    for c in curves:
        checks = [t for t, _ in c.points()]
        if isinstance(c, TabulatedCurve):
            checks += _segment_critical(c.table, t0, exponent)
        for t in checks:
            q = c.hrr(t)
            if t > t0 and q > 0:
                worst = max(worst, q / (t - t0) ** exponent)
    growth_time = (q_ref / worst) ** (1.0 / exponent)

    # decay: only when every curve has burnt out by the end
    decay = {}
    end_q = _exact_max(curves, duration)
    t_last_peak = max(t for t in knot_times
                      if _exact_max(curves, t) >= peak * (1 - REL_TOL))
    if end_q <= 0 and t_last_peak < duration:
        ds = t_last_peak
        tail = [t for t in knot_times if t > ds]
        samples = set(tail)
        for a, b in zip([ds] + tail, tail):
            samples.update(a + (b - a) * i / DECAY_SUBSTEPS
                           for i in range(1, DECAY_SUBSTEPS))
        # the fit must be out no earlier than the curves
        t_out = next(t for t in tail
                     if all(_exact_max(curves, u) <= 0 for u in tail
                            if u >= t))
        need = t_out - ds
        # just after ds the decay must fall no faster than the curves:
        # P (1 - s/D)^m ~ P - m P s / D against a slope of -k
        h = 1e-6 * max(duration, 1.0)
        k = (peak - _exact_max(curves, ds + h)) / h
        if k > 0:
            need = max(need, decay_exponent * peak / k)
        for t in samples:
            q = _exact_max(curves, t)
            if q <= 0:
                continue
            frac = min((q / peak) ** (1.0 / decay_exponent), 1 - REL_TOL)
            need = max(need, (t - ds) / (1.0 - frac))
        if need > 0:
            decay = dict(decay_start=ds, decay_time=need,
                         decay_exponent=decay_exponent)

    return PowerLawCurve(exponent, growth_time, peak, duration, q_ref=q_ref,
                         t_start=t0, **decay)
