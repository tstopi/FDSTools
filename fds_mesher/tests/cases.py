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


def l_shaft_tunnel(dx=0.2):
    """Solid rock 10 x 10 x 8 with an L-shaped tunnel and a vertical shaft.

    Leg 1 runs along x (portal at x=0), leg 2 along y, the shaft rises to the
    domain top where an OPEN vent sits.
    """
    s = head((0, 10, 0, 10, 0, 8), (round(10 / dx), round(10 / dx), round(8 / dx)))
    s += obst((0, 10, 0, 10, 0, 8))
    s += hole((0, 8, 1, 3, 1, 3))      # leg 1
    s += hole((6, 8, 1, 9, 1, 3))      # leg 2
    s += hole((6, 8, 7, 9, 1, 8))      # shaft
    s += open_vent((0, 0, 1, 3, 1, 3))
    s += open_vent((6, 8, 7, 9, 8, 8))
    return s


def geom_text(verts, faces, gid="tube", stride=4):
    v = ",".join(f"{x:.6f}" for p in verts for x in p)
    f = []
    for t in faces:
        f += [t[0] + 1, t[1] + 1, t[2] + 1] + ([1] if stride == 4 else [])
    return f"&GEOM ID='{gid}', SURF_ID='INERT',\n VERTS={v},\n FACES={','.join(map(str, f))} /\n"


def tube_geom(L=16.0, r=0.9, c=(4.0, 4.0), n=24, rings=None, stride=4):
    """Polygonal cylinder shell along x, end caps removed (open ends)."""
    import math
    rings = rings or max(1, round(L / 0.5))
    verts = []
    for i in range(rings + 1):
        x = L * i / rings
        for k in range(n):
            a = 2 * math.pi * k / n
            verts.append((x, c[0] + r * math.cos(a), c[1] + r * math.sin(a)))
    faces = []
    for i in range(rings):
        for k in range(n):
            a, b = i * n + k, i * n + (k + 1) % n
            c2, d = a + n, b + n
            faces += [(a, b, c2), (b, d, c2)]
    return verts, faces


def geom_tunnel(dx=0.2):
    """Empty 16 x 8 x 8 box, GEOM tube (open ends) with OPEN vents at its ends."""
    s = head((0, 16, 0, 8, 0, 8), (round(16 / dx), round(8 / dx), round(8 / dx)))
    s += geom_text(*tube_geom())
    s += open_vent((0, 0, 3.3, 4.7, 3.3, 4.7)) + open_vent((16, 16, 3.3, 4.7, 3.3, 4.7))
    return s


def odd_rock_tunnel(dx=0.2):
    """Rock tunnel in a 16 x 7.4 x 7 box: y has 37 cells; the tunnel spans y 0.2..7.2."""
    s = head((0, 16, 0, 7.4, 0, 7), (80, 37, 35))
    s += obst((0, 16, 0, 7.4, 0, 7)) + hole((0, 16, 0.2, 7.2, 3, 5))
    s += open_vent((0, 0, 0.2, 7.2, 3, 5)) + open_vent((16, 16, 0.2, 7.2, 3, 5))
    return s
