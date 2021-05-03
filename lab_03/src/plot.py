from matplotlib import pyplot as plt

def add_figure(*args, **kwargs):
    if 'figure_id' in kwargs:
        plt.figure(kwargs['figure_id'])

    plt.subplot(kwargs['subplot'])
    plt.plot(kwargs['x'], kwargs['y'], label=kwargs['label'])
    plt.xlabel(kwargs['x_label'])
    plt.ylabel(kwargs['y_label'])
    plt.grid(kwargs['grid'])


def show():
    plt.show()
