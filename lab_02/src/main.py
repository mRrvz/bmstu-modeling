import plot

from numpy import arange
from modeling import params, get_runge_kutta, get_Rp

def main():
    I = params.I0
    Uc = params.Uc0
    h = 1e-6

    I_plot = []
    t_plot = []
    for t in arange(0, 0.0005, h):
        Rp = get_Rp(I)
        I, Uc = get_runge_kutta(t, I, Uc, h, Rp)

        t_plot.append(t)
        I_plot.append(I)
        print(f"DEBUG Rp: {Rp}, I: {I}, Uc: {Uc}")

    plot.add_figure(figure_id=1, subplot=321, x=t_plot, y=I_plot, label='I', x_label='t', y_label='I', grid=True)
    plot.show()

if __name__ == "__main__":
    main()