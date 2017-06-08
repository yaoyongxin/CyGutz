import matplotlib.pyplot as plt
from future.builtins.iterators import zip

'''Simple pdf plot for quick visual check.
'''


def xy_plot(x_list, y_list, xlabel='x', ylabel='y', fsave='test.pdf'):
    plt.figure()
    plt.plot(x_list, y_list)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(fsave)


def xy2_plot(x_list, y_list, pattern_list, label_list,
        xlabel='x', ylabel='y', fsave='test.pdf'):
    plt.figure()
    for x, y, pattern, label in zip(x_list, y_list, pattern_list, label_list):
        plt.plot(x, y, pattern, label=label)
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(fsave)
