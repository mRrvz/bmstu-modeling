from matplotlib import pyplot as plt


def add_figures(*args, **kwargs):
    for figure in kwargs['figures']:
        plt.plot(figure['x'], figure['y'])

    plt.xlabel(kwargs['x_label'])
    plt.ylabel(kwargs['y_label'])
    plt.grid()


def show():
    plt.show()
