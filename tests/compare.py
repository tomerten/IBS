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
    12: "Conte-Martini Simpson Decade Tailcut",
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
            sim_input["coupling"],
            sim_input["threshold"],
        )
        return res
    else:
        res = ibslib.runODEBMAD(
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
        )
        return res


def run_all(sim_input: dict, sim_type: str) -> dict:
    result = {}
    for m in range(1, 14):
        sim_input["model"] = m
        result[_MODEL_MAP[m]] = run_single(sim_input, sim_type)
    sim_input["model"] = 0
    return result


def main(infile="sim_input_test_single.json"):
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

        # plot(df, sim_input)
    else:
        print("in single")
        result = run_single(sim_input, sim_type)
        df = pd.DataFrame.from_dict(result, orient="columns")
        df.to_csv(sim_input["outfile"], index=False)

        sim_input["coupling"] = 100
        sim_input["nsteps"] = 30
        sim_input["stepsize"] = 0.2
        sim_type = "man"
        result = run_single(sim_input, sim_type)
        df = pd.DataFrame.from_dict(result, orient="columns")
        df.to_csv(sim_input["outfile"] + "2", index=False)

        # plot(df, sim_input)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
# eodf
