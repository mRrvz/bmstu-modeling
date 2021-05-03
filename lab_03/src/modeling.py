from collections import namedtuple
from math import pow
from scipy.interpolate import InterpolatedUnivariateSpline
from numpy import arange

import plot

Params = namedtuple('Params', 'Np l T0 Tconst, sigma F0 alpha h')

params = Params(
    1.4, 0.2, 300, 400, 5.668 * pow(10, -12), 100, 0.05, pow(10, -4)
)

fst_table = (
    (
        300, 500, 800, 1100, 2000, 2400
    ),
    (
        1.36 * pow(10, -2), 1.63 * pow(10, -2), 1.81 * pow(10, -2),
        1.98 * pow(10, -2), 2.50 * pow(10, -2), 2.74 * pow(10, -2)
    )
)

snd_table = (
    (
        293, 1278, 1528, 1677, 2000, 2400
    ),
    (
        2.0 * pow(10, -2), 5.0 * pow(10, -2), 7.8 * pow(10, -2),
        1.0 * pow(10, -1), 1.3 * pow(10, -1), 2.0 * pow(10, -1)
    ),
)

def interpolate(x_pts, y_pts, order=1):
    return InterpolatedUnivariateSpline(x_pts, y_pts, k=order)


def p(k_t, t, n):
    return 4 * params.Np * params.Np * params.sigma * k_t(t[n]) * pow(t[n], 3)


def f(k_t, t, n):
    return 4 * params.Np *params.Np + params.sigma * k_t(t[n]) * pow(params.T0, 4)


def x_right(l_t, t, n):
    return (l_t(t[n]) + l_t(t[n + 1])) / 2


def x_left(l_t, t, n):
    return (l_t(t[n]) + l_t(t[n - 1])) / 2


def p_right(k_t, t, n):
    return (p(k_t, t, n) + p(k_t, t, n + 1)) / 2


def p_left(k_t, t, n):
    return (p(k_t, t, n) + p(k_t, t, n - 1)) / 2


def f_right(k_t, t, n):
    return (f(k_t, t, n) + f(k_t, t, n + 1)) / 2


def f_left(k_t, t, n):
    return (f(k_t, t, n) + f(k_t, t, n - 1)) / 2


def A(l_t, t, n):
    return (l_t(t[n]) + l_t(t[n - 1])) / 2 / params.h


def B(l_t, k_t, t, n):
    return A(l_t, t, n) + C(l_t, t, n) + 4 * params.Np * params.Np * params.sigma * k_t(t[n]) * pow(t[n], 3) * params.h


def C(l_t, t, n):
    return (l_t(t[n]) + l_t(t[n + 1])) / 2 / params.h


def D(k_t, t, n):
    return 4 * params.Np * params.Np + params.sigma * k_t(t[n]) * pow(params.T0, 4) * params.h


def get_right_conditions(l_t, k_t, t):
    K0 = x_right(l_t, t, 0) + pow(params.h, 2) / 8 * p_right(k_t, t, 0) + pow(params.h, 2) / 4 * p(k_t, t, 0)
    M0 = pow(params.h, 2) / 8 * p_right(k_t, t, 0) - x_right(l_t, t, 0)
    P0 = params.h * params.F0 + pow(params.h, 2) / 4 * (f_right(k_t, t, 0) + f_left(k_t, t, 0))

    return K0, M0, P0


def get_left_conditions(k_t, l_t, t, n):
    Kn = x_left(l_t, t, n) / params.h - params.alpha - params.h * p(k_t, t, n) / 4 - params.h * p_left(k_t, t, n) / 8
    Mn = x_left(l_t, t, n) / params.h - params.h * p_left(k_t, t, n) / 8
    Pn = -(params.alpha * params.T0 + (f_right(k_t, t, n) + f_left(k_t, t, n)) / 4 * params.h)

    return Kn, Mn, Pn


def start():
    l_t = interpolate(fst_table[0], fst_table[1])
    k_t = interpolate(snd_table[0], snd_table[1])
    t = [0 for _ in range(int(1 / params.h) + 2)]

    K0, M0, P0 = get_right_conditions(l_t, k_t, t)

    xi_list = [0]
    eta_list = [0]
    x_list = list()

    x = 0
    n = 0
    while x + params.h < 1:
        x_list.append(x)

        xi_list.append(C(l_t, t, n) / (B(l_t, k_t, t, n) - A(l_t, t, n) * xi_list[n]))
        eta_list.append((D(k_t, t, n) + A(l_t, t, n) * xi_list[n]) / (B(l_t, k_t, t, n) - A(l_t, t, n) * xi_list[n]))

        n += 1
        x += params.h

    x_list.extend([x + params.h, x + params.h * 2])

    Kn, Mn, Pn = get_left_conditions(k_t, l_t, t, n)

    t[n] = (Pn - Mn * xi_list[n]) / (Kn + Mn * xi_list[n])
    for i in range(n - 1, -1, -1):
        t[i] = xi_list[i + 1] * t[i + 1] + eta_list[i + 1]

    plot.add_figure(figure_id=1, subplot=321, x=x_list, y=t, label='T (x)', x_label='x', y_label='T', grid=True)
    plot.show()
