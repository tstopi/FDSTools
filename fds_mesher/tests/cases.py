"""Synthetic FDS input builders for the tests. Default dx is 0.2 m."""


def _f(*v):
    return ",".join(f"{x:g}" for x in v)


def head(xb=None, ijk=None):
    s = "&HEAD CHID='t', TITLE='synthetic / test' /\n"
    if xb:
        s += f"&MESH IJK={_f(*ijk)}, XB={_f(*xb)} /\n"
    return s


def obst(xb, extra=""):
    return f"&OBST XB={_f(*xb)}{extra} /\n"


def hole(xb):
    return f"&HOLE XB={_f(*xb)} /\n"


def open_vent(xb):
    return f"&VENT XB={_f(*xb)}, SURF_ID='OPEN' /\n"


def rock_tunnel(L=16.0, W=4.0, H=4.0, ty=(1, 3), tz=(1, 3), dx=0.2):
    """Solid box with a straight tunnel along x carved by a HOLE."""
    ijk = (round(L / dx), round(W / dx), round(H / dx))
    s = head((0, L, 0, W, 0, H), ijk)
    s += obst((0, L, 0, W, 0, H))
    s += hole((0, L, ty[0], ty[1], tz[0], tz[1]))
    s += open_vent((0, 0, ty[0], ty[1], tz[0], tz[1]))
    s += open_vent((L, L, ty[0], ty[1], tz[0], tz[1]))
    return s


def shell_tunnel(t=0.2, gap=False, L=16.0, W=8.0, H=8.0, dx=0.2):
    """Empty box with a thin-wall tunnel (interior y,z in 3..5) along x."""
    ijk = (round(L / dx), round(W / dx), round(H / dx))
    s = head((0, L, 0, W, 0, H), ijk)
    lo, hi = 3 - t, 5 + t
    if gap:   # floor with a one-cell hole at x=8..8.2
        s += obst((0, 8, lo, hi, lo, 3)) + obst((8.2, L, lo, hi, lo, 3))
    else:
        s += obst((0, L, lo, hi, lo, 3))
    s += obst((0, L, lo, hi, 5, hi))
    s += obst((0, L, lo, 3, 3, 5)) + obst((0, L, 5, hi, 3, 5))
    s += open_vent((0, 0, 3, 5, 3, 5)) + open_vent((L, L, 3, 5, 3, 5))
    return s
