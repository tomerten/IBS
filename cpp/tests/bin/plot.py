import glob

import matplotlib.gridspec as gridspec
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

if __name__ == "__main__":
    csvfiles = glob.glob("*.csv")

    dfs = [pd.read_csv(f) for f in csvfiles]

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(ncols=3, nrows=1)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    ax1.set_title(r"$\varepsilon_x$")
    ax2.set_title(r"$\varepsilon_y$")
    ax3.set_title(r"$\sigma_s$")

    for df in dfs:
        ax1.plot(df.t, df.ex)
        ax2.plot(df.t, df.ey)
        ax3.plot(df.t, df.sigs)

        # ax1.set_yscale("log")
        ax2.set_yscale("log")

    # plt.tick_params(axis="y", which="minor")
    # ax2.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    plt.show()
