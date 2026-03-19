
#!/usr/bin/env python
"""Script to make model date conversion convenient, converting elapsed model seconds to or from dates"""

import datetime
import pandas as pd
import re
import sys
import os.path
import click

from schimpy.th_io import read_th, write_th, is_elapsed, get_column_formats


@click.group()
@click.help_option("-h", "--help")
def model_time_cli():
    """Convert elapsed model seconds to or from dates"""
    pass


def _parse_datetime_token(timestr):
    return datetime.datetime(*list(map(int, re.split(r"[^\d]", timestr))))


def _is_timestamp_token(token):
    try:
        _parse_datetime_token(token)
        return True
    except Exception:
        return False


def _is_timestamp_header_token(token):
    token = token.strip().lower()
    return token in {"datetime", "date", "time"}


def _iter_timestamp_data_lines(lines):
    """
    Yield only actual timestamped data lines from a dated .th file.

    This tolerates pandas-written MultiIndex column headers like:
        variable ...
        location ...
        datetime
        2000-01-01T00:00 ...

    Assumption: if a timestamped file has a formal index-name row, the
    index name is one of datetime/date/time.
    """
    saw_header = False

    for iline, raw_line in enumerate(lines):
        line = raw_line.rstrip("\n")
        stripped = line.strip()

        if not stripped:
            continue
        if stripped.startswith("#"):
            continue

        first = stripped.split()[0]

        if _is_timestamp_header_token(first):
            saw_header = True
            continue

        if _is_timestamp_token(first):
            yield iline, raw_line
            continue

        # Before the datetime row, tolerate pandas MultiIndex header rows
        # like "variable ..." and "location ...".
        if not saw_header:
            continue

        # After a datetime header, anything else is malformed.
        raise ValueError(
            f"Encountered non-data line after timestamp header at line {iline}: {line}"
        )


@model_time_cli.command(
    name="to_elapsed",
    help="Interpret model times in elapsed seconds and translate between calendar time and elapsed. The script requires a subcommand like: $ model_time.py to_elapsed. You can also get subject-specific help on a subcommand by typing $ model_time.py subcommand --help",
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


@model_time_cli.command(name="to_date")
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
        sdtime = datetime.datetime(*list(map(int, re.split(r"[^\d]", s))))
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


@model_time_cli.command(name="clip")
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
        start_dt = datetime.datetime(*list(map(int, re.split(r"[^\d]", start))))
        scliptime = datetime.datetime(*list(map(int, re.split(r"[^\d]", clip_start))))
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
        for _, line in _iter_timestamp_data_lines(thfile):
            splitline = line.split()
            timestr = splitline[0]
            mdtm = _parse_datetime_token(timestr)
            mdelta = mdtm - scliptime
            msec = mdelta.total_seconds()
            exact = msec == 0.0
            outline = line.replace(timestr, "%-10.1f " % msec, 1)
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
        last_timestr = None

        for iline, line in _iter_timestamp_data_lines(thfile):
            splitline = line.split()
            if skip_nan and splitline[-1] == "nan":
                continue
            if len(splitline) < 2:
                continue

            timestr = splitline[0]
            last_timestr = timestr

            try:
                mdtm = _parse_datetime_token(timestr)
            except Exception:
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
                outline = line.replace(timestr, "%-10.1f " % msec, 1)
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

        if last_timestr is not None:
            print("Last time string processed: {}".format(last_timestr))

        if no_record and prev_outline is not None:
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
        start = datetime.datetime(*list(map(int, re.split(r"[^\d]", start))))
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
        start = datetime.datetime(*list(map(int, re.split(r"[^\d]", start))))
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
        time_basis = datetime.datetime(*list(map(int, re.split(r"[^\d]", time_basis))))

    if not isinstance(input, (pd.DataFrame, pd.Series)):
        in_df = read_elapsed(input)
    else:
        in_df = input

    out_df = in_df.copy()
    out_df.index = pd.to_datetime(time_basis) + pd.to_timedelta(
        round(in_df.index.to_series() * elapsed_fac), unit="s"
    )

    return out_df


def read_elapsed(input):
    """Take input dataframe or th filename, returns dataframe
    
    If the file contains inline comments (after #), they are stored in a '__comment__' column.
    """

    if not isinstance(input, (pd.DataFrame, pd.Series)):
        # Extract inline comments manually
        comments_dict = {}
        data_lines = []
        header_tokens = None
        
        with open(input, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                # Capture header tokens from comment line if present
                if line.startswith('#') and header_tokens is None:
                    candidate = line.lstrip('#').strip()
                    if candidate:
                        tokens = candidate.split()
                        if any(not re.match(r"^[-+]?\d*\.?\d+$", t) for t in tokens):
                            header_tokens = tokens
                    continue
                # Process data lines
                if line.strip():
                    # Extract comment if present
                    if '#' in line:
                        data_part, comment_part = line.split('#', 1)
                        data_lines.append(data_part.rstrip())
                        # Store comment with leading # and space
                        comments_dict[len(data_lines) - 1] = '# ' + comment_part.lstrip()
                    else:
                        data_lines.append(line)
        
        # Write temporary file without comments for pandas to read
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tmp', delete=False) as tmp:
            for line in data_lines:
                tmp.write(line + '\n')
            tmp_path = tmp.name
        
        try:
            out_df = pd.read_table(
                tmp_path, sep=r"\s+", index_col=0, header=None, dtype=float
            )

            # Apply header tokens from comment line if it matches column count
            if header_tokens:
                data_col_count = len(out_df.columns)
                if len(header_tokens) == data_col_count + 1:
                    out_df.columns = header_tokens[1:]
                elif len(header_tokens) == data_col_count:
                    out_df.columns = header_tokens

            # Add comments column if any comments were found
            if comments_dict:
                out_df['__comment__'] = None
                for idx, comment in comments_dict.items():
                    if idx < len(out_df):
                        out_df.iloc[idx, out_df.columns.get_loc('__comment__')] = comment
        finally:
            import os as os_module
            os_module.remove(tmp_path)

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
        start = datetime.datetime(*list(map(int, re.split(r"[^\d]", start))))
    if not type(timestamps) == list:
        timestampse = [timestamps]
    for stamp in timestamps:
        # assume mtime argument is a date time and try to parse to datetime
        mdtime = datetime.datetime(*list(map(int, re.split(r"[^\d]", stamp))))
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
