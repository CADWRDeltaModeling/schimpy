#!/usr/bin/env python
"""Script to make model date conversion convenient, converting elapsed model seconds to or from dates"""

import datetime
import pandas as pd
import re
import sys
import os.path
import click


@click.group()
@click.help_option("-h", "--help")
def model_time_cli():
    """Convert elapsed model seconds to or from dates"""
    pass


@model_time_cli.command(
    help="Interpret model times in elapsed seconds and translate between calendar time and elapsed. The script requires a subcommand like: $ model_time.py to_elapsed. You can also get subject-specific help on a subcommand by typing $ model_time.py subcommand --help"
)
@click.argument("dated_input", nargs=-1, required=True)
@click.option(
    "--start",
    required=False,
    type=str,
    help="Starting date and time basis for output if the input is a file.",
)
@click.option("--annotate", is_flag=True, default=False, help="Annotate output.")
@click.option(
    "--step",
    default=None,
    type=float,
    help="Model time step. If given, answer will be the integer time step.",
)
@click.option(
    "--out",
    default=None,
    type=str,
    help="Name of output file. If input is a *.th file the file will be converted and output to this file, otherwise printed to screen",
)
@click.option(
    "--skip_nan", is_flag=True, default=False, help="Skip a record with nan if True"
)
@click.help_option("-h", "--help")
def to_elapsed(dated_input, start, annotate, step, out, skip_nan):
    """Convert input datetime string or *.th file with datetime column to equivalent output in elapsed seconds."""
    s = start
    inputfile = list(dated_input)
    dt = step
    try:
        sdtime = datetime.datetime(*list(map(int, re.split(r"[^\d]", s))))
    except Exception:
        raise click.ClickException(f"Could not convert start time to datetime: {s}")
    click.echo(f"Model start time given: {sdtime}")
    if len(inputfile) > 1 and inputfile[0].endswith(".th"):
        raise click.ClickException("Only one file argument allowed at a time")
    th = len(inputfile) == 1 and os.path.exists(inputfile[0])
    outpath = out
    if outpath and not th:
        raise click.ClickException(
            "Outfile option only allowed if the input is a *.th file"
        )
    if th:
        infile = inputfile[0]
        file_to_elapsed(infile, s, outpath, annotate, skip_nan)
    else:
        describe_timestamps(inputfile, s, dt)


@model_time_cli.command()
@click.argument("elapsed_input", nargs=-1, required=True)
@click.option(
    "--start",
    required=True,
    type=str,
    help="Start time in ISO-like format 2009-03-12T00:00:00. Time part is optional.",
)
@click.option(
    "--step",
    default=None,
    type=float,
    help="Model time step in seconds. If given, answer will be the integer time step.",
)
@click.option(
    "--elapsed_unit",
    default="s",
    help="Time unit of input file. Must be either 's' for seconds or 'd' for days. Only used for files",
)
@click.option(
    "--time_format",
    default="%Y-%m-%dT%H:%M",
    help="Time format for output, e.g. the default is %%Y-%%m-%%dT%%H:%%M:%%S for 2009-03-14T22:40:00. Only used when converting fields.",
)
@click.option("--annotate", is_flag=True, default=False, help="Annotate output.")
@click.option(
    "--out",
    default=None,
    type=str,
    help="Name of output file. If input is a *.th file the file will be converted and output to this file, otherwise printed to screen",
)
@click.help_option("-h", "--help")
def to_date(elapsed_input, start, step, elapsed_unit, time_format, annotate, out):
    """Convert input elapsed seconds or *.th file with elapsed seconds as the time column to equivalent output with a datetime or annotated with datetimes."""
    s = start
    input = list(elapsed_input)
    dt = step
    try:
        sdtime = datetime.datetime(*list(map(int, re.split("[^\d]", s))))
    except Exception:
        raise click.ClickException(f"Could not convert start time to datetime: {s}")
    if len(input) > 1 and input[0].endswith(".th"):
        raise click.ClickException("Only one file argument allowed at a time")
    th = len(input) == 1 and os.path.exists(input[0])
    outpath = out
    if outpath and not th:
        raise click.ClickException(
            "Outfile option only allowed if the input is a *.th file"
        )
    if annotate and not th:
        raise click.ClickException(
            "Annotate option only allowed if the input is a *.th file"
        )
    if th:
        infile = input[0]
        file_to_timestamp(
            infile,
            s,
            outpath,
            annotate=annotate,
            elapsed_unit=elapsed_unit,
            time_format=time_format,
        )
    else:
        describe_elapsed(input, s, dt)


@model_time_cli.command()
@click.argument("elapsed_input", nargs=-1, required=True)
@click.option(
    "--start",
    required=True,
    help="Start time in ISO-like format 2009-03-12T00:00:00. Time part is optional.",
)
@click.option("--clip_start", required=True, help="Starting date for output.")
@click.option(
    "--out",
    default=None,
    type=str,
    help="Name of output file. If input is a *.th file the file will be converted and output to this file, otherwise printed to screen",
)
@click.help_option("-h", "--help")
def clip(elapsed_input, start, clip_start, out):
    """Clip (subset) an input file in elapsed time to a new, later, start date"""
    input = list(elapsed_input)
    try:
        start_dt = datetime.datetime(*list(map(int, re.split("[^\d]", start))))
        scliptime = datetime.datetime(*list(map(int, re.split("[^\d]", clip_start))))
    except Exception:
        raise click.ClickException("Could not convert start or clip_start to datetime.")
    th = len(input) == 1 and os.path.exists(input[0])
    if not th:
        raise click.ClickException("Clipping command requires file (.th) input")
    infile = input[0]
    outpath = out
    if not outpath:
        outfile = sys.stdout
    else:
        outfile = open(outpath, "w")
    with open(infile, "r") as thfile:
        prev_use = False
        prev_outline = None
        for line in thfile:
            if (
                line
                and len(line) > 1
                and not (
                    line.startswith("#")
                    or line.startswith("date")
                    or line.startswith("time")
                )
            ):
                splitline = line.split()
                timestr = splitline[0]
                msec_orig = float(timestr)
                mdtm = start_dt + datetime.timedelta(seconds=msec_orig)
                mdelta = mdtm - scliptime
                msec = mdelta.total_seconds()
                exact = msec == 0.0
                outline = line.replace(timestr, "%-10.1f " % msec)
                use = msec >= 0.0
                if use and not prev_use:
                    if prev_outline and not exact:
                        outfile.write(prev_outline)
                    outfile.write(outline)
                elif use and prev_use:
                    outfile.write(outline)
                prev_outline = outline
                prev_use = use
    if outfile != sys.stdout:
        outfile.close()


def prune_dated(name):
    return name.replace("_dated", "")


def multi_file_to_elapsed(input_files, output, start, name_transform="prune_dated"):
    import os
    import glob

    if name_transform is None:
        name_transform = lambda x: x
    elif name_transform == "prune_dated":
        name_transform = prune_dated

    if isinstance(input_files, list):
        inputs = input_files
        if isinstance(output, list):
            if not len(input_files) == len(output):
                raise ValueError(
                    "If inputs and outputs are both lists, they should be the same size"
                )
            outputs = output
        elif not os.path.isdir(output):
            raise ValueError(
                "output was a scalar, but not a valid directory name: {}".format(output)
            )
        else:
            outputs = [
                os.path.join(output, name_transform(y))
                for y in [os.path.split(x)[1] for x in input_files]
            ]
    else:
        is_glob = True
        is_dir = os.path.isdir(output)
        if isinstance(output, list) or not is_dir:
            raise ValueError(
                "If using blob search, output must be a directory not a list: {}".format(
                    output
                )
            )
        inputs = glob.glob(input_files)
        if len(inputs) == 0:
            raise ValueError("No files matched pattern: {}".format(input_files))
        outputs = [
            os.path.join(output, name_transform(y))
            for y in [os.path.split(x)[1] for x in inputs]
        ]
    print(inputs)
    print(outputs)
    for ifn, ofn in zip(inputs, outputs):
        print(ifn)
        file_to_elapsed(ifn, start, ofn)


def file_to_elapsed(infile, start, outpath=None, annotate=False, skip_nan=False):

    if not isinstance(start, datetime.datetime):
        start = datetime.datetime(*list(map(int, re.split(r"[^\d]", start))))

    if not outpath:
        outfile = sys.stdout
    else:
        outfile = open(outpath, "w")
    with open(infile, "r") as thfile:
        prev_use = False
        prev_outline = None
        no_record = True
        for iline, line in enumerate(thfile):
            if (
                line
                and len(line) > 1
                and not (
                    line.startswith("#")
                    or line.startswith("date")
                    or line.startswith("time")
                )
            ):
                splitline = line.split()
                if skip_nan and splitline[-1] == "nan":
                    continue
                if len(splitline) < 2:
                    continue
                timestr = splitline[0]
                if len(timestr) == 10:
                    # Only got the date,not the time
                    timestr += " %s" % splitline[1]
                use = True
                try:
                    mdtm = datetime.datetime(
                        *list(map(int, re.split("[^\d]", timestr)))
                    )
                except:
                    print(line)
                    raise ValueError(
                        "Could not parse time {} in line {}".format(timestr, iline)
                    )
                mdelta = mdtm - start
                mdtime = start + mdelta
                msec = mdelta.total_seconds()
                exact = msec == 0.0
                if annotate:
                    outline = "%s (%s,  %s)\n" % (line, mdtime, mdelta)
                    use = True
                else:
                    outline = line.replace(timestr, "%-10.1f " % msec)
                    use = msec >= 0.0

                if use and not prev_use:
                    if prev_outline and not exact:
                        outfile.write(prev_outline)
                    outfile.write(outline)
                    no_record = False
                elif use and prev_use:
                    outfile.write(outline)
                    no_record = False
                prev_outline = outline
                prev_use = use
        print("Last time string processed: {}".format(timestr))
        if no_record:
            if prev_outline is not None:
                outfile.write(prev_outline)

    if outfile != sys.stdout:
        outfile.close()


def file_to_timestamp(
    infile,
    start,
    outpath=None,
    annotate=False,
    elapsed_unit="s",
    time_format="%Y-%m-%dT%H:%M ",
):
    if elapsed_unit == "s":
        elapsed_fac = 1.0
    elif elapsed_unit == "d":
        elapsed_fac = 24 * 3600
    else:
        raise ValueError("elapsed_unit must be 's' or 'd'")

    if type(start) == str:
        start = datetime.datetime(*list(map(int, re.split("[^\d]", start))))
    if not outpath:
        outfile = sys.stdout
    else:
        outfile = open(outpath, "w")
    with open(infile, "r") as thfile:
        for line in thfile:
            if line and len(line) > 4:
                splitline = line.split()
                timestr = splitline[0]
                try:
                    msec = round(float(timestr) * elapsed_fac)
                except:
                    raise ValueError(
                        "elapsed time %s could not be coerced to float. \nDid you mean to convert a date to elapsed (use to_elapsed)?"
                        % timestr
                    )
                mdelta = datetime.timedelta(seconds=msec)
                mdtime = start + mdelta
                if annotate:
                    outline = "%s (%s,  %s)\n" % (line, mdtime, mdelta)
                else:
                    mdtimestr = mdtime.strftime(time_format)
                    outline = line.replace(timestr, mdtimestr, 1)
                outfile.write(outline)
    if outfile != sys.stdout:
        outfile.close()


def describe_elapsed(times, start, dt=None):
    if type(start) == str:
        start = datetime.datetime(*list(map(int, re.split("[^\d]", start))))
    for elapsed in times:
        elapsed = elapsed.lower()
        if elapsed.endswith("d") or elapsed.endswith("days") or elapsed.endswith("day"):
            mday = elapsed.split("d")[0]
            msec = float(mday) * (24.0 * 3600)  # mtime is in days
        elif elapsed.endswith("s") or elapsed.endswith("sec"):
            msec = elapsed.split("s")[0]
            msec = float(msec)
        else:
            # assuming that the input elapsed time is in seconds
            msec = float(elapsed)
        mdelta = datetime.timedelta(seconds=msec)
        print("Model seconds:  %s" % msec)
        print("Elapsed time:   %s" % mdelta)
        mdtime = start + mdelta
        print("Model datetime: %s" % mdtime)
        if dt:
            remain = abs(msec % dt)
            if abs(msec % dt) > 1.0e-6:
                print("Model step:    %s" % (msec / dt))
                print(
                    "\n Input time is %s seconds past the last even model time step\n"
                    % remain
                )
            else:
                print("Model step:    %s\n" % int(msec / dt))
        else:
            print("\n")


def elapsed_to_timestamp(input, time_basis, elapsed_unit="s"):
    """Take input dataframe or th file and convert to timestamp, returns dataframe"""

    if elapsed_unit == "s":
        elapsed_fac = 1.0
    elif elapsed_unit == "d":
        elapsed_fac = 24 * 3600
    else:
        raise ValueError("elapsed_unit must be 's' or 'd'")

    if type(time_basis) == str:
        time_basis = datetime.datetime(*list(map(int, re.split("[^\d]", time_basis))))

    if not isinstance(input, (pd.DataFrame, pd.Series)):
        in_df = pd.read_table(input, sep="\s+", comment="#", index_col=0, header=None)

    out_df = in_df.copy()
    out_df.index = pd.to_datetime(time_basis) + pd.to_timedelta(
        round(in_df.index.to_series() * elapsed_fac), unit="s"
    )

    return out_df


def is_elapsed(input):
    """Check whether the input file is 'elapsed' or 'timestamped'.
    Uses the first line of the file to determine this.
    If there are any non-numeric characters then that's a timestamped file and is_elapsed would return False
    """

    with open(input, "r") as f:
        first_line = f.readline().strip()
    is_elapsed = all(re.match(r"^[-+]?\d*\.?\d+$", s) for s in first_line.split())

    return is_elapsed


def read_th(input, time_basis=None, to_timestamp=True, elapsed_unit="s"):
    """Read input file and return dataframe.
    Automatically converts to timestamp if time_basis is supplied unless to_timestamp is False
    """

    # check if first line contains headers (and thus is elapsed)s
    if is_elapsed(input):
        # read elapsed file
        if time_basis is None:
            raise ValueError(
                f"The input file {input} is an elapsed .th file and there's no time_basis specified."
            )
        if to_timestamp:
            out_df = elapsed_to_timestamp(input, time_basis, elapsed_unit=elapsed_unit)
        else:
            out_df = pd.read_table(
                infile, sep="\s+", comment="#", index_col=0, header=None
            )

    else:
        # read timestamped file
        out_df = pd.read_table(input, sep="\s+", index_col="datetime", comment="#")
        out_df.index = pd.to_datetime(out_df.index, format="%Y-%m-%dT%H:%M")

    # monotonic increase check, finds any repeat datetimes and/or a mixup of order
    mon_inc = all(x < y for x, y in zip(out_df.index, out_df.index[1:]))
    if not mon_inc:
        # prints the row(s) where monotonicity is broken
        bad_rows = in_df.loc[
            in_df.index.to_series().diff() < pd.to_timedelta("0 seconds")
        ]
        raise ValueError(f"Non-monotonic datetime index found:\n{bad_rows}")

    return out_df


def get_headers(infile, no_index=True):
    if infile.split(".")[-1] == "th":
        with open(infile, "r") as headin:
            for line in headin:
                if line.startswith("#"):
                    continue
                headers = line.split()
                if len(headers) == 1:
                    headers = line.split(",")
                break

        if no_index:
            headers = headers[1:]
    elif infile.split(".")[-1] == "in":
        ss_in = read_source_sink_in(infile)[0]
        headers = ss_in.name.values

    return headers


def describe_timestamps(timestamps, start, dt=None):
    if type(start) == str:
        start = datetime.datetime(*list(map(int, re.split("[^\d]", start))))
    if not type(timestamps) == list:
        timestampse = [timestamps]
    for stamp in timestamps:
        # assume mtime argument is a date time and try to parse to datetime
        mdtime = datetime.datetime(*list(map(int, re.split("[^\d]", stamp))))
        mdelta = mdtime - start
        msec = mdelta.total_seconds()
        print("Datetime:      %s" % mdtime)
        print("Elapsed time:  %s" % mdelta)
        print("Model seconds: %s" % msec)
        if dt:
            remain = abs(msec % dt)
            if abs(msec % dt) > 1.0e-6:
                print("Model step:    %s" % (msec / dt))
                print(
                    "\n Input time is %s seconds past the last model time step\n"
                    % remain
                )
            else:
                print("Model step:    %s\n" % int(msec / dt))
        else:
            print("\n")


if __name__ == "__main__":
    model_time_cli()
