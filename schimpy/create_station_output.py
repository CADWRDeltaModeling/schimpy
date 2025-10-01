#!/usr/bin/env python
import click
import os
import sys
import pandas as pd
import schimpy.station as station
import glob
import numpy as np


def read_staout_lite(fname, station_infile):
    """Read a SCHISM staout_* file into a pandas DataFrame
    A "light" version of `read_staout` from schimpy.
    The original `read_staout` requires `reftime`, which requires a prior knowledge of the start time.

    Parameters
    ----------
    fname : file path
        Path to input staout file or a variable name in ["elev", "air pressure", "wind_x", "wind_y",  "temp", "salt", "u", "v", "w"] whose
        1-index will be mapped to a name like staout_1 for elev

    station_infile : str or DataFrame
        Path to station.in file or DataFrame from read_station_in
    """
    if isinstance(station_infile, str):
        station_in = station.read_station_in(station_infile)
    else:
        station_in = station_infile

    staout = pd.read_csv(fname, index_col=0, sep=r"\s+", header=None)
    station_index = station_in.index.copy()
    staout.columns = station_index
    staout.columns = [f"{loc}_{subloc}" for loc, subloc in staout.columns]

    return staout


def write_staout(staout: pd.DataFrame, out_path):
    staout.index = staout.index.map("{:0.6E}".format)
    staout.to_csv(
        out_path,
        sep="\t",
        header=None,
        na_rep="-999",
        float_format="%.6e",
    )


def create_station_output(
    input_type, old_input, new_input, out_dir, overwrite_existing, target
):
    file_list = glob.glob(target)

    # check for overwriting
    for file in file_list:
        out_path = os.path.join(out_dir, os.path.basename(file))
        if os.path.exists(out_path):
            if overwrite_existing == True:
                print(f"Warning: {out_path} already exists and will be overwritten.")
            else:
                raise Exception(
                    f"Error: {out_path} already exists. Rename it or set 'overwrite_existing=True' to overwrite."
                )

    # station.in
    if input_type == "station":

        station_new = station.read_station_in(new_input)

        # Get the list of columns from station_new
        out_columns = [f"{loc}_{subloc}" for loc, subloc in station_new.index]

        for file in file_list:
            staout_old = read_staout_lite(file, old_input)

            # Identify columns in station_new that are missing from staout_old, fill with NaN
            missing_cols = [col for col in out_columns if col not in staout_old.columns]
            for col in missing_cols:
                staout_old[col] = np.nan

            # Reorder new staout to match the column order of B
            staout_new = staout_old[out_columns]

            write_staout(staout_new, os.path.join(out_dir, os.path.basename(file)))


@click.command()
@click.option(
    "--input_type",
    default=None,
    required=True,
    type=str,
    help="Type of input. Must be one of: `station`, `fluxflag`, `fluxlines`",
)
@click.option(
    "--old_input",
    default=None,
    required=True,
    type=click.Path(exists=True),
    help="Reference input file from a previous run",
)
@click.option(
    "--new_input",
    default=None,
    required=True,
    type=click.Path(exists=True),
    help="station.in file for the current simulation",
)
@click.option(
    "--out_dir",
    default=".",
    required=False,
    type=click.Path(exists=True),
    help="output location",
)
@click.option(
    "--overwrite_existing",
    default=None,
    required=False,
    type=bool,
    help="If true, existing output files will be overwritten. If false, warning given without generating file.",
)
@click.argument("target", metavar="target")
@click.help_option("-h", "--help")
def create_station_output_cli(
    input_type, old_input, new_input, out_dir, overwrite_existing, target
):
    """
    Prepares station-related SCHISM inputs
    target: SCHISM input associated with `input_type`. For example, use following input_type and target pair:\n
    station : staout*
    """

    create_station_output(
        input_type, old_input, new_input, out_dir, overwrite_existing, target
    )


if __name__ == "__main__":
    create_station_output_cli()
