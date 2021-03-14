
from numpy import arange

import plot
from modeling import params, get_runge_kutta, get_Rp, get_T0, get_m

def main():
    I = params.I0
    Uc = params.Uc0
    h = 1e-6

    I_res = []
    Rp_res = []
    U_res = []
    T_res = []
    IRp_res = []
    t_res = []

    for t in arange(0, 0.0003, h):
        T0 = get_T0(I)
        Rp = get_Rp(I, T0, get_m(I))

        I, Uc = get_runge_kutta(t, I, Uc, h, Rp)

        if t > h:
            t_res.append(t)
            I_res.append(I)
            Rp_res.append(Rp)
            U_res.append(Uc)
            T_res.append(T0)
            IRp_res.append(I * Rp)

        print(f"DEBUG Rp: {Rp}, I: {I}, Uc: {Uc}, T0: {T0}, Rk: {params.Rk}")

    plot.add_figure(figure_id=1, subplot=321, x=t_res, y=I_res, label='I (t)', x_label='t', y_label='I', grid=True)
    plot.add_figure(subplot=322, x=t_res, y=Rp_res, label='Rp (t)', x_label='t', y_label='Rp', grid=True)
    plot.add_figure(subplot=323, x=t_res, y=U_res, label='Uc (t)', x_label='t', y_label='U', grid=True)
    plot.add_figure(subplot=324, x=t_res, y=IRp_res, label='I * Rp (t)', x_label='t', y_label='I * Rp', grid=True)
    plot.add_figure(subplot=325, x=t_res, y=T_res, label='T0 (t)', x_label='t', y_label='T0', grid=True)
    plot.show()


if __name__ == "__main__":
    main()
