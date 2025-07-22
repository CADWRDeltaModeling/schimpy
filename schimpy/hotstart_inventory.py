#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import re
import datetime as dtm
import os
import pandas as pd
from schimpy.param import *
from vtools import days, seconds


def expected_hotstarts(run_start, dt, nday, hot_freq):
    """
    Generate expected hotstarts based on run parameters.
    Parameters
    ----------
    run_start : datetime
        The start time of the simulation run.
    dt : int
        The time step in seconds.
    nday : int
        The number of days for the simulation run.
    hot_freq : pandas.tseries.offsets.DateOffset
        The frequency at which hotstarts are generated.
    Returns
    -------
    pandas.DataFrame
        A DataFrame with a DatetimeIndex representing the expected hotstart times
        and a column "iteration" indicating the corresponding iteration numbers.
    Notes
    -----
    This function calculates the expected hotstart times and their corresponding
    iteration numbers based on the provided simulation parameters. It does not
    perform any file inventory or check for existing hotstart files.
    """
    """Generate expected hotstarts based on run parameters."""

    end = run_start + days(nday)
    t = run_start
    iters = []
    times = []
    dt_freq = pd.tseries.offsets.Second(dt)
    iters_per_hot = hot_freq / dt_freq
    itr = 0
    print("Projected hotstarts (not inventory of files)")
    while t < end:
        times.append(t)
        iters.append(itr)
        t = t + hot_freq
        itr = itr + int(iters_per_hot)
        print(f"t={t}, itr={itr}")
    ndx = pd.DatetimeIndex(times)
    df = pd.DataFrame(index=ndx, data=iters)
    df.columns = ["iteration"]
    df.index.name = "datetime"
    return df


def hotstart_inventory(
    run_start=None,
    dt=None,
    nday=None,
    workdir=".",
    paramfile=None,
    hot_freq=None,
    expected=False,
):
    """
    Create an inventory of existing hotstarts or expected hotstarts.

    Whether the inventory is for existing or expected hotstarts depends on
    whether the workdir is an outputs directory or a study directory.

    Parameters
    ----------
    run_start : str or pandas.Timestamp, optional
        The start time of the simulation run in ISO-like format (e.g., '2013-12-03').
    dt : int, optional
        The time step in seconds of the model. Default is None.
    nday : int, optional
        The number of days in the simulation run or the maximum to catalog. Default is None.
    workdir : str, optional
        The working directory, which may be a launch or outputs directory. Default is '.'.
    paramfile : str, optional
        The name of the param.nml file if it is used to infer runtime. Default is None.
    hot_freq : str or pandas.tseries.offsets.DateOffset, optional
        The hotstart frequency in pandas frequency terms (e.g., '5D'). Default is None.
    expected : bool, optional
        Flag to generate expected hotstarts instead of an inventory of existing files. Default is False.

    Returns
    -------
    pandas.DataFrame or None
        A DataFrame containing the inventory of hotstarts with a DatetimeIndex and a column "iteration".
        If no hotstarts are found and `expected` is False, returns None.
    """

    run_start = pd.to_datetime(run_start)
    if type(hot_freq) == str:
        hot_freq = pd.tseries.frequencies.to_offset(hot_freq)

    hots = glob.glob(os.path.join(workdir, "hotstart_000000_*.nc"))
    is_existing = len(hots) > 0
    print(workdir)
    param_needed = (
        (dt is None)
        or (run_start is None)
        or (hot_freq is None and not is_existing)
        or pd.isnull(run_start)
    )
    print(
        "Existing hotstarts found:",
        is_existing,
        " params read from file: " + str(param_needed),
    )
    if param_needed:
        if paramfile is None or paramfile == "":
            paramfile = "param.nml"
        print("got here", paramfile, workdir)
        if os.path.exists(os.path.join(workdir, "..", paramfile)):
            paramfile = os.path.join(workdir, "..", paramfile)
        params = read_params(paramfile)
        dt_param = params["dt"]
        run_start_param = params.run_start
        run_len_param = params["rnday"]
        hot_freq_param = params.hotstart_freq
        print(
            f"dt={dt_param}, run_start={run_start_param}, run_len={run_len_param}, hot_freq={hot_freq_param}"
        )
        if run_start is None or pd.isnull(run_start):
            run_start = run_start_param
        if dt is None:
            dt = dt_param
        if nday is None or nday == 0:
            nday = run_len_param
        if hot_freq is None:
            hot_freq = hot_freq_param

    if expected:
        expected_df = expected_hotstarts(run_start, dt, nday, hot_freq)
        print(expected_df)
        return expected_df

    if is_existing:
        return hotstart_inventory_exist(run_start, dt, workdir)
    else:
        return None


def hotstart_inventory_exist(start, dt=90, workdir=".", do_print=True):
    """
    Check for the existence of hotstart inventory files and generate a DataFrame
    mapping iterations to corresponding datetime values.
    Parameters
    ----------
    start : str or pandas.Timestamp
        The start time as a string or pandas.Timestamp object.
    dt : int, optional
        Time step in seconds between iterations (default is 90).
    workdir : str, optional
        Directory to search for hotstart files (default is '.').
    do_print : bool, optional
        Whether to print the iterations and corresponding times (default is True).
    Returns
    -------
    pandas.DataFrame
        A DataFrame with datetime as the index and iteration numbers as the column.
    Notes
    -----
    - The function searches for hotstart files in the specified `workdir` directory.
    - If no files are found in `workdir`, it searches in the current directory.
    - The function assumes hotstart files follow the naming pattern
      "hotstart_000000_<iteration>.nc".
    """

    print("Hotstart Inventory")
    if isinstance(start, str):
        start = pd.Timestamp(start)
    hots = glob.glob(os.path.join(workdir, "hotstart_000000_*.nc"))
    if len(hots) == 0:
        hots = glob.glob(os.path.join("hotstart_0000_*.nc"))
    hots.sort()
    iters = [int(x.split("_")[2].replace(".nc", "")) for x in hots]

    iters.sort()
    times = [start + seconds(x * dt) for x in iters]
    df = pd.DataFrame(index=times, data=iters)
    df.columns = ["iteration"]
    df.index.name = "datetime"

    if do_print:
        for it, t in zip(iters, times):
            print("{}: {}".format(it, t))

    return df


def create_arg_parser():
    parser = argparse.ArgumentParser(
        "Lookup station metadata by partial string match on id or name"
    )
    parser.add_argument(
        "--dt", default=90, type=int, help="Time step in seconds of model"
    )
    parser.add_argument(
        "--run_start", default="", help="Start time in iso-like format, e.g. 2013-12-03"
    )
    parser.add_argument(
        "--nday",
        default=0,
        type=int,
        help="Number of days in simulation (rnday) or maximum to catalog",
    )
    parser.add_argument(
        "--workdir",
        default=".",
        type=str,
        help="Working directory, which is the outputs dir",
    )
    parser.add_argument(
        "--paramfile",
        default="",
        type=str,
        help="Name of param.nml file if file is used to infer runtime. If neither params nor paramfile provided, ./param.nmo or ../param.nml will be tried. ",
    )
    parser.add_argument(
        "--hot_freq",
        default=None,
        help="Hotstart frequency in pandas freq terms (e.g. '5D')",
    )
    parser.add_argument(
        "--expected",
        action="store_true",
        help="Flag to generate expected hotstarts instead of inventory of existing files",
    )
    return parser


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    hotstart_inventory(
        args.run_start,
        args.dt,
        args.nday,
        args.workdir,
        args.paramfile,
        args.hot_freq,
        args.expected,
    )


if __name__ == "__main__":
    main()
