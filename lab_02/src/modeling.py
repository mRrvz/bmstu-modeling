import math

from collections import namedtuple
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline

RungeCoeffs = namedtuple('RungeCoeffs', 'kn pn')
Params = namedtuple('Params', 'R Le Lk Ck Rk Uc0 I0 Tw')

params = Params(
    0.35, 12, 60e-6, 150e-6, 1, 1500, 0.5, 2000
)

fst_table = [
    [0.5, 6730, 0.50],
    [1.0, 6790, 0.55],
    [5.0, 7150, 1.7],
    [10.0, 7270, 3],
    [50.0, 8010, 11],
    [200.0, 9185, 32],
    [400.0, 10010, 40],
    [800.0, 11140, 41],
    [1200.0, 12010, 39]
]

snd_table = [
    [4000,  0.031],
    [5000, 0.27],
    [6000, 2.05],
    [7000, 6.06],
    [8000, 12],
    [9000, 19.9],
    [10000, 29.6],
    [11000, 41.1],
    [12000, 54.1],
    [13000, 67.7],
    [14000, 81.5]
]


def T(z, T0, m):
    return T0 + (params.Tw - T0) * z ** m


def f(x, y, z, Rp):
    return -((params.Rk + Rp) * y - z) / params.Lk


def phi(x, y, z):
    return -y / params.Ck


def sigma(T):
    return interpolate(T, list(map(lambda x: x[0], snd_table)), list(map(lambda x: x[1], snd_table)))


def get_T0(I):
    return interpolate(I, list(map(lambda x: x[0], fst_table)), list(map(lambda x: x[1], fst_table)))


def get_m(I):
    return interpolate(I, list(map(lambda x: x[0], fst_table)), list(map(lambda x: x[2], fst_table)))


def get_Rp(I):
    integral = integrate.quad(lambda z: sigma(T(z, get_T0(I), get_m(I))) * z, 0, 1)
    return params.Le / (2 * math.pi * params.R ** 2 * integral[0])


def interpolate(x, x_pts, y_pts, order=1):
    spline = InterpolatedUnivariateSpline(x_pts, y_pts, k=order)
    return spline(x)


def get_current_addition(h, coeffs, i, order):
    if i == 0:
        return 0, 0, 0
    elif i == order - 1:
        return h, coeffs.kn, coeffs.pn

    return h / 2, coeffs.kn / 2, coeffs.pn / 2


def get_next_members(current_y, current_z, coeffs):
    k_sum = 0
    p_sum = 0

    for i in range(len(coeffs)):
        if i > 0 and i < len(coeffs) - 1:
            k_sum += 2 * coeffs[i].kn
            p_sum += 2 * coeffs[i].pn
        else:
            k_sum += coeffs[i].kn
            p_sum += coeffs[i].pn

    divider = 2 * (len(coeffs) - 2) + 2
    return current_y + k_sum / divider, current_z + p_sum / divider


def get_runge_kutta(x, y, z, h, Rp, order=4):
    coeffs = [RungeCoeffs(0, 0) for x in range(order)]

    for i in range(order):
        curr_h, y_add, z_add = get_current_addition(h, coeffs[i - 1], i, order)
        coeffs[i] = RungeCoeffs(h * f(x + curr_h, y + y_add, z + z_add, Rp), h * phi(x + curr_h, y + y_add, z + z_add))


    return get_next_members(y, z, coeffs)