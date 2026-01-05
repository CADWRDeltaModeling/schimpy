#!/usr/bin/env python
# -*- coding: utf-8 -*-

import click
import pandas as pd
from vtools import hours


def test_date_label_correct(fname_th):
    """Make sure that the date column in the file doens't have comment marker and raise error if it does"""
    label_correct = False
    with open(fname_th, "r") as testread:
        for line in testread:
            if (
                line.startswith("date ")
                or line.startswith("time ")
                or line.startswith("datetime ")
            ):
                good = True
            if line.startswith("#time") or line.startswith("#date"):
                raise ValueError(
                    f"Input th file ({fname_th}) has comment marker (#) on date column name. This is deprecated"
                )


def test_template_times(template, dt):
    rounded_index = template.index.round(dt)
    if (rounded_index != template.index).any(axis=None):
        raise ValueError("The template times and dt must be neatly aligned")


def interpolate_structure(
    template_th,
    output_th,
    dt=None,
    int_cols=["install", "ndup_weir", "ndup_culvert", "ndup_pipe"],
):
    """Interpolate a dated th "template" for a structure
    The interpolation adds time detail, maintain data types.
    Comments are removed and floating precision is hardwired to two digits

    The input template file must have a header and time stamps that are neat
    with respect to dt. This is checked.

    Parameters
    ----------
    template_th : str
        Path to the template file
    output_th : str
        Name of the output file
    dt : interval_type
        Number of hours between changes, default is 1 hour. Can be any time (string, delta, offset) that can be coerced to offset
    int_cols : list
        List of column names that should be integers. This can be a superset, and it shouldn't change much if you name your template columns the standard way.
    """
    if dt is None:
        dt = hours(1)
    else:
        dt = pd.tseries.frequencies.to_offset(dt)
    test_date_label_correct(template_th)
    th_orig = pd.read_csv(
        template_th,
        comment="#",
        sep=r"\s+",
        header=0,
        dtype={"install": int, "ndup": int},
        index_col=0,
        parse_dates=[0],
    )
    print(th_orig)
    test_template_times(th_orig, dt)
    th = th_orig.copy()
    cols = th_orig.columns
    agg_rules = {label: "first" for label in cols}
    # agg_rules = { "A": "mean", "B": "sum", "C": "first", "D": "last",}
    th = th_orig.resample(dt).mean()
    for c in cols:
        if c in int_cols:
            th[c] = th[c].ffill().astype(int)
        else:
            th[c] = th[c].interpolate()
    th = th.drop_duplicates(keep="first")
    th.to_csv(
        "test.th",
        sep=" ",
        float_format="%.2f",
        date_format="%Y-%m-%dT%H:%M",
        header=True,
    )


def ensure_offset(arg):
    if isinstance(arg, str):
        return pandas.tseries.frequencies.to_offset()
    else:
        return arg


@click.command()
@click.option(
    "--template_th",
    default=None,
    help="th file path containing skeleton of operations.",
)
@click.option("--output_th", default=None, help="Output file path")
@click.option(
    "--dt",
    default=None,
    help="Time step of output. Input template timestamps must be neat with respect to this",
)
@click.option(
    "--int_cols",
    multiple=True,
    default=None,
    help="List of column names to treat as integers. Default is install, ndup, ndup_culvert, ndup_pipe, ndup_weir.",
)
def interpolate_structure_cli(template_th, output_th, dt, int_cols):
    """Interpolate a dated th template for a structure"""
    if dt is None:
        dt = hours(1)
    if int_cols is None or len(int_cols) == 0:
        int_cols = ["install", "ndup_weir", "ndup_culvert", "ndup_pipe"]
    interpolate_structure(
        template_th=template_th, output_th=output_th, dt=dt, int_cols=int_cols
    )


if __name__ == "__main__":
    interpolate_structure_cli()
