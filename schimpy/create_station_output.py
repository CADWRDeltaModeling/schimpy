#!/usr/bin/env python
import click
import os
import pandas as pd
import schimpy.station as station
from schimpy.util.yaml_load import yaml_from_file
import glob
import numpy as np
from datetime import datetime
import warnings


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


def generate_elapsed_time(run_start, hot_date, dt_in_sec, index_unit):

    if index_unit not in ["days", "seconds"]:
        raise ValueError(
            f"Invalid unit: '{index_unit}'. Must be one of ['days', 'seconds']"
        )

    run_start = datetime.strptime(run_start, "%Y-%m-%d")
    hot_date = datetime.strptime(hot_date, "%Y-%m-%d")

    total_seconds = (hot_date - run_start).total_seconds()
    elapsed_timestamps = np.arange(dt_in_sec, total_seconds + dt_in_sec, dt_in_sec)

    if index_unit == "days":
        elapsed_timestamps = elapsed_timestamps / 86400  # convert seconds to days
    elif index_unit == "seconds":
        pass

    return elapsed_timestamps


def create_station_output(
    hotstart,
    run_start,
    hot_date,
    input_type,
    old_input,
    new_input,
    out_dir,
    overwrite_existing,
    target,
):
    if hotstart:

        # Warn the user that some input arguments are not required with hotstart
        warnings.warn(
            "When --hotstart flag is present, `old_input` and `target` are ignored"
        )

        # if target is staout*, run_start and hot_date must be provided
        if run_start is None or hot_date is None:
            raise ValueError(
                "When --hotstart flag is present, `run_start` and `hot_date` must be provided"
            )

        # for hotstart, generate file list based on input_type
        if input_type == "station":
            file_list = [f"staout_{i}" for i in range(1, 10)]
        elif input_type in ["fluxflag", "fluxline"]:
            file_list = ["flux.out"]

    else:
        # if old_input is missing raise error
        if old_input is None:
            raise ValueError("`old_input` argument must be provided")

        # if target is missing raise error
        if target is None:
            if input_type == "station":
                raise ValueError(
                    "`target` argument pointing to staout must be provided. For example, staout*"
                )
            elif input_type in ["fluxflag", "fluxline"]:
                raise ValueError(
                    "`target` argument pointing to flux.out must be provided."
                )
        else:
            file_list = glob.glob(target)

    # check for overwriting
    for file in file_list:
        out_path = os.path.join(out_dir, os.path.basename(file))
        if os.path.exists(out_path):
            if overwrite_existing == True:
                warnings.warn(f"{out_path} already exists and will be overwritten.")
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

            if hotstart:

                # create empty DataFrame
                index = generate_elapsed_time(
                    run_start, hot_date, dt_in_sec=900, index_unit="seconds"
                )
                staout_new = pd.DataFrame(index=index, columns=out_columns)

            else:
                staout_old = read_staout_lite(file, old_input)

                # Identify columns in station_new that are missing from staout_old, fill with NaN
                missing_cols = [
                    col for col in out_columns if col not in staout_old.columns
                ]
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
        fluxflag_loc_new = [loc.lower() for loc in fluxflag_loc_new]
        out_columns = fluxflag_loc_new

        if hotstart:

            # create empty DataFrame
            index = generate_elapsed_time(
                run_start, hot_date, dt_in_sec=90, index_unit="days"
            )
            flux_out_new = pd.DataFrame(index=index, columns=out_columns)

        else:
            flux_out_old = station.read_flux_out(target, old_input, reftime=None)
            # Identify missing columns and fill with NaN
            missing_cols = [
                col for col in out_columns if col not in flux_out_old.columns
            ]
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

        fluxline_new = yaml_from_file(new_input)
        fluxline_loc_new = [item["station_id"] for item in fluxline_new["linestrings"]]
        fluxline_loc_new = [loc.lower() for loc in fluxline_loc_new]
        out_columns = fluxline_loc_new

        if hotstart:

            # create empty DataFrame
            index = generate_elapsed_time(
                run_start, hot_date, dt_in_sec=90, index_unit="days"
            )
            flux_out_new = pd.DataFrame(index=index, columns=out_columns)

        else:
            # extract station names
            fluxline_old = yaml_from_file(old_input)
            fluxline_loc_old = [
                item["station_id"] for item in fluxline_old["linestrings"]
            ]
            fluxline_loc_old = [loc.lower() for loc in fluxline_loc_old]

            # Read flux.out
            flux_out_old = station.read_flux_out(target, fluxline_loc_old, reftime=None)
            # Identify missing columns and fill with NaN
            missing_cols = [
                col for col in out_columns if col not in flux_out_old.columns
            ]
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
    "--hotstart",
    is_flag=True,
    help="Flag for creating desired SCHISM input files from scratch.",
)
@click.option(
    "--run_start",
    default=None,
    required=False,
    type=str,
    help="Run start date in the format YYYY-MM-DD. Required only when --hotstart flag is present.",
)
@click.option(
    "--hot_date",
    default=None,
    required=False,
    type=str,
    help="Hotstart date in the format YYYY-MM-DD. Required only when --hotstart flag is present.",
)
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
    required=False,
    type=click.Path(exists=True),
    help="Previous run's input file associated with `input_type`. Not required when --hotstart flag is present.",
)
@click.option(
    "--new_input",
    default=None,
    required=True,
    type=click.Path(exists=True),
    help="Current run's input file associated with `input_type`.",
)
@click.option(
    "--out_dir",
    default=".",
    required=False,
    type=click.Path(exists=True),
    help="Output location. Default is current directory.",
)
@click.option(
    "--overwrite_existing",
    default=None,
    required=False,
    type=bool,
    help="If true, existing output files will be overwritten. If false, warning given without generating file.",
)
@click.argument("target", required=False, metavar="target")
@click.help_option("-h", "--help")
def create_station_output_cli(
    hotstart,
    run_start,
    hot_date,
    input_type,
    old_input,
    new_input,
    out_dir,
    overwrite_existing,
    target,
):
    """
    Prepares station-related SCHISM inputs. \n
    target: SCHISM input associated with `input_type`. For example, use following input_type and target pair:\n
    station : staout*
    fluxflag : flux.out
    fluxline : flux.out


    Usage Examples:

    1. Create staout_* files from scratch (i.e., hotstart)
    create_station_output --hotstart --run_start 2021-01-01 --hot_date 2021-02-01 --input_type station --new_input station_new.in --out_dir . --overwrite_existing TRUE

    2. Create staout_* files based on a previous run
    create_station_output --input_type station --old_input station_old.in --new_input station_new.in --out_dir . --overwrite_existing TRUE staout*

    3. Create flux.out from scratch (i.e., hotstart) using fluxflag.prop
    create_station_output --hotstart --run_start 2021-01-01 --hot_date 2021-02-01 --input_type fluxflag --new_input fluxflag_new.prop --out_dir . --overwrite_existing TRUE

    4. Create flux.out based on a previous run using fluxflag.prop
    create_station_output --input_type fluxflag --old_input fluxflag_old.prop --new_input fluxflag_new.prop flux.out

    5. Create flux.out from scratch (i.e., hotstart) using flow_station_xsects.yaml
    create_station_output --hotstart --run_start 2021-01-01 --hot_date 2021-02-01 --input_type fluxline --new_input flow_station_xsects_new.yaml --out_dir . --overwrite_existing TRUE

    """

    create_station_output(
        hotstart,
        run_start,
        hot_date,
        input_type,
        old_input,
        new_input,
        out_dir,
        overwrite_existing,
        target,
    )

if __name__ == "__main__":
    create_station_output_cli()