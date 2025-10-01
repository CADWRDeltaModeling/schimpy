#!/usr/bin/env python
import click
import os
import sys
import pandas as pd
import schimpy.station as station
from schimpy.util.yaml_load import yaml_from_file
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


def write_station_output(
    staout: pd.DataFrame,
    out_path,
    index_map,
    float_format,
    na_rep,
):
    staout.index = staout.index.map(index_map.format)
    staout.to_csv(
        out_path,
        sep="\t",
        header=None,
        na_rep=na_rep,
        float_format=float_format,
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

            write_station_output(
                staout_new,
                os.path.join(out_dir, os.path.basename(file)),
                index_map="{:0.6E}",
                float_format="%.6e",
                na_rep="-999",
            )

    # fluxflag.prop
    elif input_type == "fluxflag":

        fluxflag_loc_new = station.station_names_from_file(new_input)
        flux_out_old = station.read_flux_out(target, old_input, reftime=None)
        fluxflag_loc_new = [loc.lower() for loc in fluxflag_loc_new]

        # Identify missing columns and fill with NaN
        out_columns = fluxflag_loc_new
        missing_cols = [col for col in out_columns if col not in flux_out_old.columns]
        for col in missing_cols:
            flux_out_old[col] = np.nan

        # Reorder to match column order
        flux_out_new = flux_out_old[out_columns]

        write_station_output(
            flux_out_new,
            os.path.join(out_dir, os.path.basename(file)),
            index_map="{:0.6f}",
            float_format="%.4e",
            na_rep="-999",
        )

    # fluxflag.prop
    elif input_type == "fluxline":

        fluxline_old = yaml_from_file(old_input)
        fluxline_new = yaml_from_file(new_input)

        # extract station names
        fluxline_loc_old = [item["station_id"] for item in fluxline_old["linestrings"]]
        fluxline_loc_old = [loc.lower() for loc in fluxline_loc_old]
        fluxline_loc_new = [item["station_id"] for item in fluxline_new["linestrings"]]
        fluxline_loc_new = [loc.lower() for loc in fluxline_loc_new]

        # Read flux.out
        flux_out_old = station.read_flux_out(
            target, fluxline_loc_old, reftime=None
        )

        # Identify missing columns and fill with NaN
        out_columns = fluxline_loc_new
        missing_cols = [col for col in out_columns if col not in flux_out_old.columns]
        for col in missing_cols:
            flux_out_old[col] = np.nan

        # Reorder to match column order
        flux_out_new = flux_out_old[out_columns]

        write_station_output(
            flux_out_new,
            os.path.join(out_dir, os.path.basename(file)),
            index_map="{:0.6f}",
            float_format="%.4e",
            na_rep="-999",
        )


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
