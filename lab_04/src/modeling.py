from collections import namedtuple
from math import pow, fabs

Params = namedtuple('Params', 'a1 b1 c1 m1 a2 b2 c2 m2 alpha0 alphaN l T0 R F0 h t eps')

params = Params(
    0.0134, 1, 4.35e-4, 1, 2.049, 0.563e-3, 0.528e5, 1, 0.05, 0.01, 10, 300, 0.5, 50, 1e-3, 1, 1e-2
)


def approximation_plus(func, n, step):
    return (func(n) + func(n + step)) / 2


def approximation_minus(func, n, step):
    return (func(n) + func(n - step)) / 2


def k(T):
    return params.a1 * (params.b1 + params.c1 * pow(T, params.m1))


def c(T):
    return params.a2 + params.b2 * pow(T, params.m2) - (params.c2 / pow(T, 2))


def alpha(x):
    d = (params.alphaN * params.l) / (params.alphaN - params.alpha0)
    c = -params.alpha0 * d
    return c / (x - d)


def p(x) :
    return alpha(x) * 2 / params.R


def f(x):
    return alpha(x) * 2 * params.T0 / params.R


def A(T):
    return params.t / params.h * approximation_minus(k, T, params.t)


def D(T):
    return params.t / params.h * approximation_plus(k, T, params.t)


def B(x, T):
    return A(T) + D(T) + params.h * c(T) + params.h * params.t * p(x)


def F(x, T):
    return params.h * params.t * f(x) + T * params.h * c(T)


def get_left_conditions(T):
    c_plus = approximation_plus(c, T[0], params.t)
    k_plus = approximation_plus(k, T[0], params.t)

    K0 = params.h / 8 * c_plus + params.h / 4 * c(T[0]) + params.t / params.h * k_plus + \
            params.t * params.h / 8 * p(params.h / 2) + params.t * params.h / 4 * p(0)

    M0 = params.h / 8 * c_plus - params.t / params.h * k_plus + params.t * params.h / 8 * p(params.h / 2)

    P0 = params.h / 8 * c_plus * (T[0] + T[1]) + params.h / 4 * c(T[0]) * T[0] + \
            params.F0 * params.t + params.t * params.h / 8 * (3 * f(0) + f(params.h))

    return K0, M0, P0


def get_right_conditions(T):
    c_minus = approximation_minus(c, T[-1], params.t)
    k_minus = approximation_minus(k, T[-1], params.t)

    KN = params.h / 8 * c_minus + params.h / 4 * c(T[-1]) + params.t / params.h * k_minus + \
            params.t * params.alphaN + params.t * params.h/ 8 * p(params.l - params.h / 2) + \
            params.t * params.h / 4 * p(params.l)

    MN = params.h / 8 * c_minus - params.t / params.h * k_minus + \
            params.t * params.h / 8 * p(params.l - params.h / 2)

    PN = params.h / 8 * c_minus * (T[-1] + T[-2]) + params.h / 4 * c(T[-1]) * T[-1] + \
            params.t * params.alphaN * params.T0 + params.t * params.h / 4 * (f(params.l) + f(params.l - params.h / 2))

    return KN, MN, PN


def get_new_T(T):
    K0, M0, P0 = get_left_conditions(T)
    KN, MN, PN = get_right_conditions(T)

    xi = [0, -M0 / K0]
    eta = [0, P0 / K0]

    x = params.h
    n = 1

    while x + params.h < params.l:
        Tn = T[n]
        denominator = (B(x, Tn) - A(Tn) * xi[n])

        next_xi = D(Tn) / denominator
        next_eta = (F(x, Tn) + A(Tn) * eta[n]) / denominator

        xi.append(next_xi)
        eta.append(next_eta)

        n += 1
        x += params.h

    T_new = [0 for _ in range(n + 1)]
    T_new[n] = (PN - MN * eta[n]) / (KN + MN * xi[n])

    for i in range(n - 1, -1, -1):
        T_new[i] = xi[i + 1] * T_new[i + 1] + eta[i + 1]

    return T_new


def simple_iteration_method():
    T = [params.T0 for _ in range(int(params.l / params.h) + 1)]
    T_new = [0 for _ in range(int(params.l / params.h) + 1)]

    result = [T]
    ti = 0

    epsilon_condition = True
    while epsilon_condition:
        T_prev = T
        current_max = 1

        while current_max >= 1:
            T_new = get_new_T(T_prev)
            current_max = fabs((T[0] - T_new[0]) / T_new[0])

            for T_i, Tnew_i in zip(T, T_new):
                d = fabs(T_i - Tnew_i) / Tnew_i
                if d > current_max:
                    current_max = d

            T_prev = T_new

        result.append(T_new)
        ti += params.t

        epsilon_condition = False
        for T_i, Tnew_i in zip(T, T_new):
            if fabs((T_i - Tnew_i) / Tnew_i) > params.eps:
                epsilon_condition = True

        T = T_new

    return result, ti
