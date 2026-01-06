#!/usr/bin/env python
# -*- coding: utf-8 -*

"""
Archive SCHISM time series outputs from station, flux and extracted output from a "study".

Goals:
1. Allow outputs from multiple scenarios to be collected in one directory by adding labels to the name.
2. Translate from fortranesque format to csv with # comments and datetime indexes.
3. Assure that the datetimes are not affected by floating point precision of the time columns
4. Add scenario info as data so that origins can be identified and data can be concatenated.

"""
from schimpy.station import *
from schimpy.param import *
import pandas as pd
import os
import yaml
import re
import glob
import click


def archive_time_series(
    rundir,
    ardir,
    runstart,
    scenario_label=None,
    scenario_data={},
    staouts={},
    stationfile="station.in",
    time_sharded=False,
):
    """Archive time series from rundir/outputs to ardir

    rundir/param.nml must point to a file with the correct run start
    """
    if not os.path.exists(ardir):
        raise ValueError(
            f"Archive directory {ardir} does not exist or permission problem"
        )

    if scenario_label is None:
        raise ValueError("Scenario label must be provided")

    do_flux = "flow" in stationfile or "flux" in stationfile
    if do_flux:
        archive_flux(
            rundir, ardir, scenario_label, stationfile, scenario_data, runstart=runstart
        )
    else:
        # archive_staout(rundir,ardir,scenario_label,scenario_data,runstart)
        # statouts={"fort.18": "sjrfrac"}
        # staouts={"sjrfrac.dat": "sjrfrac"}
        tunit = "d" if "bp" in stationfile else "s"
        archive_staout(
            rundir,
            ardir,
            scenario_label,
            scenario_data,
            staouts,
            stationfile,
            time_unit=tunit,
            runstart=runstart,
            time_sharded=time_sharded,
        )


def shard_number(x):
    a0 = os.path.splitext(x)[0]
    number_re = re.compile(r".+_(\d+)")
    aday = int(number_re.match(a0).group(1))
    return aday


def get_ordered_files(loc, pat, time_sharded):
    # if not time sharded, trivially return the one file
    # as a single member list
    if not time_sharded:
        all_files = [os.path.join(loc, pat)]
    else:
        searchpath = os.path.join(loc, pat)
        all_files = glob.glob(searchpath)
        all_files.sort(key=shard_number)
    return all_files


def infer_runstart(rundir):
    """Infer runstart based on param.nml in directory rundir"""
    paramfile = os.path.join(rundir, "param.nml")
    p = read_params(paramfile)
    return p.get_run_start()


def archive_staout(
    rundir,
    ardir,
    scenario_label,
    scenario_data,
    staouts=None,
    stationfile="station.in",
    float_format="%.3f",
    time_unit="s",
    multi=False,
    elim_default=True,
    do_flux=False,
    time_sharded=False,
    runstart=None,
):
    if runstart is None:
        runstart = infer_runstart(rundir)
    if staouts is None:
        staouts = {"staout_1": "elev", "staout_5": "temp", "staout_6": "salt"}
        if time_sharded:
            raise ValueError("Time sharding unexected for staout_* files")
    for s in staouts:
        print(f"Processing {s} in run directory {rundir} time_sharded={time_sharded}")
        loc = os.path.join(rundir, "outputs")
        ofiles = get_ordered_files(loc, s, time_sharded)
        dfs = []  # For concatenation in time in case there are more than one
        varlabel = staouts[s]
        if len(ofiles) == 0:
            print(f"No files found for pattern {s}")
            continue
        for fpath in ofiles:
            df = read_staout(
                fpath,
                station_infile=os.path.join(rundir, stationfile),
                reftime=runstart,
                time_unit=time_unit,
                multi=multi,
                elim_default=elim_default,
            )
            # df.pivot()
            for item in scenario_data:
                df[item] = scenario_data[item]
            df["variable"] = varlabel
            df["rundir"] = os.path.abspath(rundir).split("/")[-1]
            df.index.name = "datetime"
            dfs.append(df)

        dfout = pd.concat(dfs, axis=0)
        scenario_fname = f"{varlabel}_{scenario_label}.csv"
        outfpath = os.path.join(ardir, scenario_fname)
        dfout.to_csv(
            outfpath,
            sep=",",
            date_format="%Y-%m-%dT%H:%M:%S",
            float_format=float_format,
        )


def archive_flux(rundir, ardir, scenario_label, stationfile, scenario_data, runstart):
    print(f"Processing flux.out in directory {rundir}")
    if runstart is None:
        runstart = infer_runstart(rundir)
    fpath = os.path.join(rundir, "outputs", "flux.out")
    df = read_flux_out(
        os.path.join(rundir, "outputs", "flux.out"),
        names=os.path.join(rundir, stationfile),
        reftime=runstart,
    )
    print(df.head())
    df.index.name = "datetime"
    for item in scenario_data:
        df[item] = scenario_data[item]

    df["rundir"] = os.path.abspath(rundir).split("/")[-1]
    varlabel = "flow"
    scenario_fname = f"{varlabel}_{scenario_label}.csv"
    outfpath = os.path.join(ardir, scenario_fname)
    df.to_csv(outfpath, sep=",", date_format="%Y-%m-%dT%H:%M:%S")


def process_extracted_scalar(
    rundir, ardir, extract_data_file, variable, model_start, station_file, output_file
):
    """Process  extracted data into a time series"""
    process_extracted_scalar(rundir, ardir, "fort.18")
    # df = read_staout(archive_staout

    print(
        f"Processing extracted data file: {extract_data_file} model_start={model_start} far={variable} output_file={output_file}"
    )
    ts_out = pd.read_csv(sextract_data_file, sep="\s+", header=None, index_col=0)
    delta_t = (ts_out.index[1] - ts_out.index[0]) * 24 * 60
    freqstr = f"{int(delta_t)}min"
    # print(f"Detected frequency = {freqstr}")
    dr = pd.date_range(
        start=model_start + pd.Timedelta(days=ts_out.index[0]),
        periods=ts_out.shape[0],
        freq=freqstr,
    )
    ts_out.index = dr
    ts_out = ts_out.resample("1d").mean()

    if station_file.endswith("bp"):
        station_df = pd.read_csv(
            station_file,
            sep="\s+",
            index_col=0,
            skiprows=[1],
            header=0,
            usecols=range(1, 5),
            comment="!",
        )
    else:
        raise ValueError("Build point file expected")

    ncols = ts_out.shape[1]
    nstat = len(station_df)
    if ncols != nstat:
        raise ValueError(
            "Number of columns in salt output {ncols} must match number of locations in bp file {nroute}"
        )

    ts_out.columns = station_df.distance
    x2_prelim = ts_out.apply(find_x2, axis=1)
    x2_prelim.to_csv(output_file, float_format="%.1f")


def main_hardwire():
    rundir = "mss_base"
    ardir = "/home/eli/archive"
    archive_time_series(
        rundir, ardir, scenario_label="mss_base", scenario_data={"year": 2021}
    )


def ex():
    print(shard_sorter("hello_1", "hello_1"))
    print(shard_sorter("hello_10", "hello_2"))
    print(shard_sorter("hello_1.out", "hello_1.out"))
    print(shard_sorter("hello_10.out", "hello_2.out"))
    print(shard_sorter("hello_2.out", "hello_120.out"))


def archive_ts(
    rundir, ardir, label, scenario_data, stationfile, extracted, run_start, time_sharded
):
    """Archive time series from one simulation in a large study with many alternatives."""

    # Parse YAML strings
    scenario_data = yaml.safe_load(scenario_data)
    if stationfile and stationfile != "station.in":
        # Only parse as YAML if it looks like YAML, otherwise keep as string
        try:
            stationfile_parsed = yaml.safe_load(stationfile)
            if isinstance(stationfile_parsed, (dict, list)):
                stationfile = stationfile_parsed
        except:
            pass  # Keep as string if YAML parsing fails

    if extracted:
        staouts = yaml.safe_load(extracted)
    else:
        staouts = None

    # Parse datetime
    runstart = pd.to_datetime(run_start) if run_start else None

    print("time_sharded", time_sharded)

    for key, val in scenario_data.items():
        if val is None:
            raise ValueError(
                "Value in scenario_data is None. On the command line this may be an omitted space after colon"
            )
    if not staouts:
        print("None")
        staouts = {"staout_1": "elev", "staout_5": "temp", "staout_6": "salt"}
    for key, val in staouts.items():
        if val is None:
            raise ValueError(
                "Value in scenario_data is None. On the command line this may be an omitted space after colon"
            )

    archive_time_series(
        rundir,
        ardir,
        scenario_label=label,
        scenario_data=scenario_data,
        staouts=staouts,
        stationfile=stationfile,
        time_sharded=time_sharded,
        runstart=runstart,
    )

    print(rundir)
    print(ardir)
    print(scenario_data)
    print(label)
    print(staouts)


@click.command(
    help="""
    Archive time series from one simulation in a large with many alternatives 
    and store in a common location with better names and formats.
    
    Examples:
    
    $ archive_ts --ardir archive_ts --rundir mss_base --label base 
      --scenario_data '{year: 2001, sjr: 1400, exports: 1500}' 
      --extracted '{flux.out: flux}' --stationfile fluxflag.prop --run_start 2021-01-01

    $ archive_ts --ardir archive_ts --rundir mss_base --label base 
      --scenario_data '{year: 2001, sjr: 1400, exports:1500}'

    $ archive_ts --ardir archive_ts --rundir mss_base --label base 
      --scenario_data '{year: 2001, sjr: 1400, exports: 1500}' 
      --stationfile outputs/south_delta.bp 
      --extracted '{fracsjr_*.out: fracsjr, fracdelta_*.out: fracdelta}' 
      --time_sharded
    """
)
@click.option(
    "--rundir",
    default=None,
    help="location of the run launch dir where param.nml resides.",
)
@click.option(
    "--ardir",
    default=None,
    help="path to the archive receiving the data",
)
@click.option(
    "--label",
    default=None,
    help="scenario label for output. Final name in archive will be variable_label.csv",
)
@click.option(
    "--scenario_data",
    default="{}",
    help="dictionary of data that identifies scenario, which will be added as a column to the output. This is very helpful for stitching.",
)
@click.option(
    "--stationfile",
    default="station.in",
    help="name of annotated station.in, build pointfile or fluxflag.prop/flow_xsects.yaml",
)
@click.option(
    "--extracted",
    default=None,
    help='dictionary of extracted data in rundir/outputs in the form of "file_name: variable_name". '
    "You can only one station.in OR flux.out OR the files associated with a single build point per invocation. "
    "On the command line note that this is quoted and requires a space after the colons. "
    'Example: "{staout_1: elev, staout_5: temp, staout_6: salt}"',
)
@click.option(
    "--run_start",
    default=None,
    help="start date of run.",
)
@click.option(
    "--time_sharded",
    is_flag=True,
    default=False,
    help="if true, assume the extracted file is sharded in time.",
)
def archive_ts_cli(
    rundir, ardir, label, scenario_data, stationfile, extracted, run_start, time_sharded
):
    """Command Line Interface for archiving time series from one simulation in a large study with many alternatives."""

    archive_ts(
        rundir,
        ardir,
        label,
        scenario_data,
        stationfile,
        extracted,
        run_start,
        time_sharded,
    )


if __name__ == "__main__":
    archive_ts_cli()
