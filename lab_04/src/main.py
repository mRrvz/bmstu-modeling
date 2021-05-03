from numpy import arange

import plot
import modeling as m


def main():
    result, ti = m.simple_iteration_method()

    x = [i for i in arange(0, m.params.l, m.params.h)]
    points = []

    for i, res in enumerate(result):
        if i % 2 == 0:
            points.append({"x": x, "y": res[:-1]})

    points.append({"x": x, "y": result[-1][:-1]})
    plot.add_figures(figures=points, x_label = "x, cм", y_label = "T, K")
    plot.show()

    points = []
    te = [i for i in range(0, ti, m.params.t)]
    point = []

    for s in arange(0, m.params.l / 3, 0.1):
        point = [j[int(s / m.params.h)] for j in result]
        points.append({"x": te, "y": point[:-1]})

    plot.add_figures(figures=points, x_label="t, ceк", y_label="T, K")
    plot.show()


if __name__ == "__main__":
    main()