# -*- coding: utf-8 -*-
"""Command line interface runode (no sub-commands)."""

import sys
from json import load
from posixpath import abspath

import click
import IBSLib as ibslib
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
    12: "Conte-MArtini Simpson Decade Tailcut",
    13: "Zimmerman Simpson Decade",
}


def read_sim_input(infile: str) -> dict:
    with open(infile, "r") as f:
        sim_input = load(f)

    return sim_input


def check_sim_input(sim_input: dict) -> str:
    _REQUIRED_KEYS = [
        "twissfile",
        "harmon",
        "voltages",
        "ex",
        "ey",
        "sigs",
        "model",  # 0 is all
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

    if all(k in sim_input for k in _MAN_KEYS):
        return "man"
    else:
        return "auto"


def run_single(sim_input: dict, sim_type: str) -> dict:
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
        )
        return res
    else:
        return ibslib.runODE(
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
        )


def run_all(sim_input: dict, sim_type: str) -> dict:
    result = {}
    for m in range(1, 14):
        sim_input["model"] = m
        result[_MODEL_MAP[m]] = run_single(sim_input, sim_type)
    sim_input["model"] = 0
    return result


def plot(df: pd.DataFrame, sim_input: dict) -> None:
    if sim_input["model"] == 0:
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(ncols=2, nrows=2)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])

        ax1.set_title(r"$\epsilon_x$")
        ax2.set_title(r"$\epsilon_y$")
        ax3.set_title(r"$\sigma_s$")

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
        plt.savefig(sim_input["plotfile"], bbox_inches="tight")
    else:
        fig = plt.figure(constrained_layout=True)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])

        ax1.set_title(r"$\varepsilon_x$")
        ax2.set_title(r"$\varepsilon_y$")
        ax3.set_title(r"$\sigma_s$")

        ax1.plot(df.t, df.ex)
        ax2.plot(df.t, df.ey)
        ax3.plot(df.t, df.sigs)
        plt.savefig(sim_input["plotfile"])


@click.command()
@click.argument("infile")
def main(infile):
    sim_input = read_sim_input(infile)
    sim_type = check_sim_input(sim_input)

    model = sim_input["model"]

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
        result = run_single(sim_input, sim_type)
        df = pd.DataFrame.from_dict(result, orient="columns")
        plot(df, sim_input)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
# eodf
