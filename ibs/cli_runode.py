# -*- coding: utf-8 -*-

"""Command line interface runode."""

import sys
from json import load
from math import sqrt
from posixpath import abspath

import click
import IBSLib as ibslib
import matplotlib
import pandas as pd
from matplotlib import pyplot as plt

# Print options menu

_MODEL_MAP = {
    0: "All",
    1: "Piwinski Smooth",
    2: "Piwinski Lattice",
    3: "Piwinski Lattice Modified",
    4: "Nagaitsev",
    5: "Nagaitsev Tailcut",
    6: "Zimmerman Simpson Decade Scaling",
    7: "Zimmerman Simpson Decade Scaling with Tailcut",
    8: "Bjorken-Mtingwa Standard Simpson",
    9: "Bjorken-Mtingwa Simpson Decade",
    10: "Bjorken-Mtingwa Simpson Decade Tailcut",
    11: "Conte-Martini Simpson Decade",
    12: "Conte-Martini Simpson Decade Tailcut",
    13: "Zimmerman Simpson Decade",
}


def read_sim_input(infile: str) -> dict:
    """Read the simulation input file (JSON).

    Parameters
    ----------
    infile : str
        JSON input file with simulation settings.

    Returns
    -------
    dict
        Simulation settings as dict.
    """
    with open(infile, "r") as f:
        sim_input = load(f)

    return sim_input


def check_sim_input(sim_input: dict) -> str:
    """Check correctness of content of simulation input file.

    Parameters
    ----------
    sim_input : dict
        output for read_sim_input

    Returns
    -------
    str
        man | auto - sets the timesteps in the simulation to man or auto.
    """
    _REQUIRED_KEYS = [
        "twissfile",
        "harmon",
        "voltages",
        "ex",
        "ey",
        "sigs",
        "model",  # 0 is all
        "method",
        "pnumber",
        "outfile",
        "plotfile",
    ]
    _MAN_KEYS = ["nsteps", "stepsize"]

    assert all(k in sim_input for k in _REQUIRED_KEYS)
    assert isinstance(sim_input["twissfile"], str)
    assert isinstance(sim_input["harmon"], list)
    assert isinstance(sim_input["voltages"], list)
    assert sim_input["model"] in list(range(14))
    assert isinstance(sim_input["outfile"], str)
    assert isinstance(sim_input["plotfile"], str)
    assert sim_input["method"] in ["rlx", "der"]

    if all(k in sim_input for k in _MAN_KEYS):
        return "man"
    else:
        return "auto"


def run_single(sim_input: dict, sim_type: str) -> dict:
    """Run simulation for single IBS model.

    Parameters
    ----------
    sim_input : dict
        output of read_sim_input
    sim_type : str
        output of check_sim_input

    Returns
    -------
    dict
        Returns a dict with the simultion output (keys: t, ex, ey, sigs, model).
    """
    twissheader = ibslib.GetTwissHeader(sim_input["twissfile"])
    twisstable = ibslib.GetTwissTable(sim_input["twissfile"])
    twisstable = ibslib.updateTwiss(twisstable)

    t = [0.0]
    exa = [sim_input["ex"]]
    eya = [sim_input["ey"]]
    sigsa = [sim_input["sigs"]]
    sigea = []

    if sim_type == "auto":
        res = ibslib.runODE(
            twissheader,
            twisstable,
            sim_input["harmon"],
            sim_input["voltages"],
            t,
            exa,
            eya,
            sigsa,
            sigea,
            sim_input["model"],
            sim_input["pnumber"],
            sim_input["coupling"],
            sim_input["threshold"],
            sim_input["method"],
        )
        return res
    else:
        res = ibslib.runODE(
            twissheader,
            twisstable,
            sim_input["harmon"],
            sim_input["voltages"],
            t,
            exa,
            eya,
            sigsa,
            sigea,
            sim_input["model"],
            sim_input["pnumber"],
            sim_input["nsteps"],
            sim_input["stepsize"],
            sim_input["coupling"],
            sim_input["method"],
        )
        return res


def run_all(sim_input: dict, sim_type: str) -> dict:
    """Runs the simulation for all available IBS models.

    Parameters
    ----------
    sim_input : dict
        output of read_sim_input
    sim_type : str
        output of check_sim_input

    Returns
    -------
    dict
        Dict with all simulation outputs,
        model column is used to distinguish the results for the different models.
    """
    result = {}
    for m in range(1, 14):
        sim_input["model"] = m
        result[_MODEL_MAP[m]] = run_single(sim_input, sim_type)
    sim_input["model"] = 0
    return result


SPINE_COLOR = "gray"


def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert columns in [1, 2]

    if fig_width is None:
        fig_width = 3.39 if columns == 1 else 6.9  # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print(
            "WARNING: fig_height too large:"
            + fig_height
            + "so will reduce to"
            + MAX_HEIGHT_INCHES
            + "inches."
        )
        fig_height = MAX_HEIGHT_INCHES

    params = {
        "backend": "ps",
        "text.latex.preamble": "\\usepackage{gensymb}",
        "axes.labelsize": 20,  # fontsize for x and y labels (was 10)
        "axes.titlesize": 20,
        "font.size": 20,  # was 10
        "legend.fontsize": 20,  # was 10
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "text.usetex": False,
        "figure.figsize": [fig_width, fig_height],
        "font.family": "serif",
        "axes.formatter.limits": [-3, 3],
        "axes.formatter.use_mathtext": True,
    }

    matplotlib.rcParams.update(params)


def format_axes(ax):

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction="out", color=SPINE_COLOR)

    return ax


def plot(df: pd.DataFrame, sim_input: dict, save=True) -> None:
    """Method to plot the ODE simulation results quickly.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the simulation results.
    save: bool
        save plot to file
    sim_input : dict
        output of read_sim_input
    """
    # latexify(columns=2)
    if sim_input["model"] == 0:
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(ncols=2, nrows=2)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])

        ax1.set_title(r"$\epsilon_x$")
        ax2.set_title(r"$\epsilon_y$")
        ax3.set_title(r"$\sigma_s$")

        ax1.grid()
        ax2.grid()
        ax3.grid()

        for m, g in df.groupby("model"):
            ax1.plot(g.t, g.ex, lw=1)
            ax2.plot(g.t, g.ey, lw=1)
            ax3.plot(g.t, g.sigs, lw=1)

            ax1.set_yscale("log")
            ax2.set_yscale("log")

            ax1.set_xlabel("t[s]")
            ax2.set_xlabel("t[s]")
            ax3.set_xlabel("t[s]")
        fig.legend(
            list(_MODEL_MAP.values())[1:],
            loc="upper center",
            bbox_to_anchor=(0.5, 1.25),
            ncol=2,
            prop={"size": 8},
        )
        # plt.tick_params(axis="y", which="minor")
        # ax2.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
        if save:
            plt.savefig(sim_input["plotfile"], bbox_inches="tight")
    else:
        print("plotting")
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(ncols=2, nrows=2)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])

        ax1.grid()
        ax2.grid()
        ax3.grid()

        ax1.set_title(r"$\epsilon_x$")
        ax2.set_title(r"$\epsilon_y$")
        ax3.set_title(r"$\sigma_s$")

        ax1.plot(df.t, df.ex, lw=1)
        ax2.plot(df.t, df.ey, lw=1)
        ax3.plot(df.t, df.sigs, lw=1)

        ax1.set_yscale("log")
        ax2.set_yscale("log")

        ax1.set_xlabel("t[s]")
        ax2.set_xlabel("t[s]")
        ax3.set_xlabel("t[s]")

        if save:
            plt.savefig(sim_input["plotfile"])


@click.command()
@click.argument("infile")
def main(infile):
    sim_input = read_sim_input(infile)
    sim_type = check_sim_input(sim_input)

    model = sim_input["model"]
    print(model)
    if model == 0:
        result = run_all(sim_input, sim_type)
        frames = []

        for model, dc in result.items():
            frame = pd.DataFrame.from_dict(dc, orient="columns")
            frame["model"] = model
            frames.append(frame)

        df = pd.concat(frames)
        print(df)
        df.to_csv(sim_input["outfile"], index=False)

        plot(df, sim_input)
    else:
        print("in single")
        result = run_single(sim_input, sim_type)
        df = pd.DataFrame.from_dict(result, orient="columns")
        df.to_csv(sim_input["outfile"], index=False)
        print(df)
        plot(df, sim_input)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
# eodf
