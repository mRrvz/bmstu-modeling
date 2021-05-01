from collections import namedtuple
from math import pow
from scipy.interpolate import InterpolatedUnivariateSpline

Params = namedtuple('Params', 'Np l T0 Tconst, sigma F0 alpha')

params = Params(
    1.4, 0.2, 300, 400, 5.668 * pow(10, -12), 100, 0.05
)

fst_table = (
    (300, 1.36 * pow(10, -2)),
    (500, 1.63 * pow(10, -2)),
    (800, 1.81 * pow(10, -2)),
    (1100, 1.98 * pow(10, -2)),
    (2000, 2.50 * pow(10, -2)),
    (2400, 2.74 * pow(10, -2)),
)

snd_table = (
    (293, 2.0 * pow(10, -2)),
    (1278, 5.0 * pow(10, -2)),
    (1528, 7.8 * pow(10, -2)),
    (1677, 1.0 * pow(10, -1)),
    (2000, 1.3 * pow(10, -1)),
    (2400, 2.0 * pow(10, -1)),
)

def T(x):
    return params.Tconst - (params.Tconst - params.T0) * x / params.l


def interpolate(x, x_pts, y_pts, order=1):
    return InterpolatedUnivariateSpline(x_pts, y_pts, k=order)(x)