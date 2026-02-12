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
                    #timestr += " %s" % splitline[1]
                    pass
                use = True
                try:
                    mdtm = datetime.datetime(
                        *list(map(int, re.split(r"[^\d]", timestr)))
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


def is_elapsed(input):
    """Check whether the input file is 'elapsed' or 'timestamped'.
    Uses the first non-comment, non-blank line of the file to determine this.
    Ignores any text after a '#' comment marker.
    """
    with open(input, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # Skip blank or comment lines
            # Ignore comments after '#'
            line = line.split("#", 1)[0].strip()
            is_elapsed = all(
                re.match(r"^[-+]?\d*\.?\d+$", s) for s in line.split() if s
            )
            return is_elapsed
    return False  # If no data lines found


def read_th(input, time_basis=None, to_timestamp=True, elapsed_unit="s", head=None):
    """Read input file and return dataframe.
    Automatically converts to timestamp if time_basis is supplied unless to_timestamp is False
    
    If the file contains inline comments (after #), they are stored in a '__comment__' column.
    """


    # check if first line contains headers (and thus is elapsed)s
    file_is_elapsed = is_elapsed(input)

    if file_is_elapsed:
        # read elapsed file
        if time_basis is None:
            raise ValueError(
                f"The input file {input} is an elapsed .th file and there's no time_basis specified."
            )
        if to_timestamp:
            out_df = elapsed_to_timestamp(input, time_basis, elapsed_unit=elapsed_unit)
        else:
            out_df = read_elapsed(input)
        
        if head is not None:
            headers = get_headers(head)
            # Exclude __comment__ column from validation
            num_data_cols = len([c for c in out_df.columns if c != '__comment__'])
            if len(headers) != num_data_cols:
                raise ValueError(
                    f"{os.path.basename(input)} has {num_data_cols} columns and the header calls for {len(headers)}!"
                )
            # Apply headers to data columns, preserving __comment__
            has_comment_col = '__comment__' in out_df.columns
            if has_comment_col:
                comment_col = out_df.pop('__comment__')
            
            out_df.columns = headers
            
            if has_comment_col:
                # Reset the index of comment_col to match the current dataframe
                comment_col.index = out_df.index
                out_df['__comment__'] = comment_col
        
        out_df.index.name = "datetime"

    else:
        # read timestamped file, extract comments (both inline and separate-line)
        comments_dict = {}
        with open(input, 'r') as f:
            data_lines = []
            header_line = None
            last_data_idx = None
            
            for line in f:
                line = line.rstrip('\n')
                # Store header line
                if header_line is None and not line.startswith('#') and not line.strip().startswith('#'):
                    if line and not any(c.isdigit() for c in line.split()[0]) if line.split() else False:
                        header_line = line
                        data_lines.append(line)
                        continue
                
                # Check if line is a comment line (starts with spaces and #)
                if line.strip() and line.lstrip().startswith('#'):
                    # This is a separate-line comment following a data line
                    if last_data_idx is not None:
                        # Store this comment for the last data line
                        comments_dict[last_data_idx] = line.strip()
                    continue
                
                # Skip pure comment lines at start
                if line.startswith('#'):
                    continue
                
                # Process data lines
                if line.strip():
                    # Extract inline comment if present
                    if '#' in line:
                        data_part, comment_part = line.split('#', 1)
                        data_lines.append(data_part.rstrip())
                        # Store inline comment
                        last_data_idx = len(data_lines) - 1
                        comments_dict[last_data_idx] = '# ' + comment_part.lstrip()
                    else:
                        data_lines.append(line)
                        last_data_idx = len(data_lines) - 1
        
        # Write temporary file without comments for pandas to read
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tmp', delete=False) as tmp:
            for line in data_lines:
                tmp.write(line + '\n')
            tmp_path = tmp.name
        
        try:
            out_df = pd.read_table(tmp_path, sep=r"\s+", index_col="datetime")
            out_df.index = pd.to_datetime(out_df.index, format="%Y-%m-%dT%H:%M")
            
            # Add comments column if any comments were found
            if comments_dict:
                out_df['__comment__'] = None
                for idx, comment in comments_dict.items():
                    # Adjust index for header line
                    if idx > 0 and idx - 1 < len(out_df):
                        out_df.iloc[idx - 1, out_df.columns.get_loc('__comment__')] = comment
        finally:
            import os as os_module
            os_module.remove(tmp_path)

    # monotonic increase check, finds any repeat datetimes and/or a mixup of order
    if not out_df.index.is_monotonic_increasing:
        # prints the row(s) where monotonicity is broken
        if file_is_elapsed:
            bad_rows = out_df.loc[
                out_df.index.to_series().diff() <= 0.
            ]
            raise ValueError(f"Non-monotonic elapsed index found:\n{bad_rows}")
        else:
            bad_rows = out_df.loc[
                out_df.index.to_series().diff() <= pd.to_timedelta("0 seconds")
            ]
            raise ValueError(f"Non-monotonic datetime index found:\n{bad_rows}")

    return out_df


def get_column_formats(filename):
    """Extract decimal format for each column from a .th file.
    
    Returns a dict mapping column names to format strings (e.g., "%.3f")
    """
    import re
    
    formats = {}
    with open(filename, 'r') as f:
        # Skip header lines
        header = None
        data_lines = []
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            if header is None and not line.startswith('datetime') and not line.startswith('date') and not line.startswith('time'):
                continue
            if header is None and (line.startswith('datetime') or line.startswith('date') or line.startswith('time')):
                header = line.split()
                continue
            if header is not None and line.strip() and not line.lstrip().startswith('#'):
                # Extract comment if present
                data_line = line.split('#')[0].rstrip() if '#' in line else line
                data_lines.append(data_line)
        
        if header is None or len(data_lines) == 0:
            return formats
        
        # Skip the time column (first column)
        for col_idx, col_name in enumerate(header[1:], start=1):
            decimals_list = []
            for data_line in data_lines:
                parts = data_line.split()
                if col_idx < len(parts):
                    value_str = parts[col_idx]
                    # Count decimals
                    if '.' in value_str:
                        decimals = len(value_str.split('.')[1])
                        decimals_list.append(decimals)
                    else:
                        decimals_list.append(0)
            
            # Use the maximum decimal places found
            if decimals_list:
                max_decimals = max(decimals_list)
                formats[col_name] = f"%.{max_decimals}f"
    
    return formats


def write_th(
    df,
    outpath,
    elapsed=None,
    ref_time=None,
    time_format="%Y-%m-%dT%H:%M",
    float_format="%.2f",
    elapsed_format="%.1f",
    include_header=None,
    index_name=None,
    source_file=None,
):
    """Write a dataframe to SCHISM .th format.

    Parameters
    ----------
    df : pandas.DataFrame or pandas.Series
        Time-indexed data. Datetime index produces a dated .th file.
        Numeric index produces an elapsed .th file.
    outpath : str or file-like
        Output .th filename or file object.
    elapsed : bool, optional
        Force elapsed (True) or dated (False) output. If None, inferred from
        the index dtype.
    ref_time : str or datetime, optional
        Reference time for converting datetime index to elapsed seconds (when
        elapsed=True) or elapsed index to datetimes (when elapsed=False).
    time_format : str
        Datetime format for dated output.
    float_format : str
        Format string for data values (e.g., "%.2f").
    elapsed_format : str
        Format string for elapsed time values (e.g., "%.1f").
    include_header : bool, optional
        Write column headers. Defaults to True for dated output and False for
        elapsed output.
    index_name : str, optional
        Name for the index column (defaults to "datetime" or "time").
    source_file : str, optional
        Path to source file to preserve comments from.
    """

    if isinstance(df, pd.Series):
        df = df.to_frame()

    if not isinstance(df, pd.DataFrame):
        raise ValueError("df must be a pandas DataFrame or Series")

    df_local = df.copy()
    is_datetime_index = pd.api.types.is_datetime64_any_dtype(df_local.index)

    if elapsed is None:
        elapsed = not is_datetime_index

    if include_header is None:
        include_header = not elapsed

    if ref_time is not None and isinstance(ref_time, str):
        ref_time = datetime.datetime(*list(map(int, re.split(r"[^\d]", ref_time))))

    # Extract comments from source file if provided
    comments = []
    col_formats = {}  # Will store per-column format strings
    if source_file and isinstance(source_file, str):
        try:
            with open(source_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        comments.append(line)
                    elif line.strip() and not line.startswith("date") and not line.startswith("time"):
                        # Stop at first data line
                        break
            # Extract column formats from source file
            col_formats = get_column_formats(source_file)
        except:
            pass  # If we can't read comments, just continue without them

    # Determine if outpath is a file object or path
    should_close = False
    if isinstance(outpath, str):
        outfile = open(outpath, mode="w", newline="\n")
        should_close = True
    else:
        outfile = outpath

    try:
        # Write comments first
        for comment in comments:
            outfile.write(comment)

        if elapsed:
            # If we have named columns, write them as a commented header line
            data_cols = [c for c in df_local.columns if c != '__comment__']
            if not include_header and data_cols:
                all_numeric = True
                for c in data_cols:
                    if isinstance(c, str):
                        if not re.match(r"^\d+$", c):
                            all_numeric = False
                            break
                    elif isinstance(c, (int, float)):
                        continue
                    else:
                        all_numeric = False
                        break
                if not all_numeric:
                    header_name = index_name or "time"
                    outfile.write("# " + header_name + " " + " ".join(map(str, data_cols)) + "\n")

            if is_datetime_index:
                if ref_time is None:
                    raise ValueError("ref_time is required to write elapsed output")
                elapsed_index = (df_local.index - ref_time).total_seconds()
            else:
                elapsed_index = df_local.index.astype(float)

            # Check if there's a __comment__ column to preserve inline comments
            has_comments = '__comment__' in df_local.columns
            if has_comments:
                comments_col = df_local.pop('__comment__')
            
            if include_header:
                header_df = df_local.copy()
                header_df.index = elapsed_index
                header_df.index.name = index_name or "time"
                
                if has_comments:
                    comments_col = header_df.pop('__comment__')
                    
                    # Write header
                    header_str = ' '.join(header_df.columns)
                    outfile.write(header_df.index.name + ' ' + header_str + '\n')
                    
                    # Write data rows with inline comments
                    for i, (idx, row) in enumerate(header_df.iterrows()):
                        tstr = elapsed_format % idx
                        values = []
                        for col_idx, (col_name, val) in enumerate(row.items()):
                            if pd.isna(val):
                                values.append('nan')
                            else:
                                try:
                                    # Use column-specific format if available
                                    format_str = None
                                    if col_name in col_formats:
                                        format_str = col_formats[col_name]
                                    elif isinstance(col_name, int) and col_formats:
                                        # Try to get format by position if columns are numeric
                                        col_list = list(col_formats.keys())
                                        if col_idx < len(col_list):
                                            format_str = col_formats[col_list[col_idx]]
                                    
                                    if format_str:
                                        values.append(format_str % float(val))
                                    else:
                                        values.append(float_format % float(val))
                                except Exception:
                                    values.append(str(val))
                        line = tstr + ' ' + ' '.join(values)
                        if i < len(comments_col) and pd.notna(comments_col.iloc[i]):
                            line = line + ' ' + comments_col.iloc[i]
                        outfile.write(line + '\n')
                else:
                    header_df.to_csv(
                        outfile,
                        sep=" ",
                        header=True,
                        float_format=float_format,
                    )
            else:
                data = df_local.to_numpy()
                for i, tval in enumerate(elapsed_index):
                    tstr = elapsed_format % tval
                    row = []
                    for val in data[i]:
                        if pd.isna(val):
                            row.append("nan")
                        else:
                            try:
                                row.append(float_format % float(val))
                            except Exception:
                                row.append(str(val))
                    line = f"{tstr} " + " ".join(row)
                    if has_comments and i < len(comments_col) and pd.notna(comments_col.iloc[i]):
                        line = line + ' ' + comments_col.iloc[i]
                    outfile.write(line + "\n")
        else:
            if is_datetime_index:
                out_df = df_local.copy()
            else:
                if ref_time is None:
                    raise ValueError("ref_time is required to write dated output")
                out_df = df_local.copy()
                out_df.index = pd.to_datetime(ref_time) + pd.to_timedelta(
                    df_local.index.astype(float), unit="s"
                )

            out_df.index.name = index_name or "datetime"
            
            # Check if there's a __comment__ column to preserve comments
            has_comments = '__comment__' in out_df.columns
            if has_comments:
                comments_col = out_df.pop('__comment__')
                
                # Write header
                if include_header:
                    header_str = ' '.join(out_df.columns)
                    outfile.write(out_df.index.name + ' ' + header_str + '\n')
                
                # Write data rows with inline comments
                for idx, row in out_df.iterrows():
                    # Format datetime
                    date_str = idx.strftime(time_format)
                    
                    # Format data values using column-specific formats if available
                    values = []
                    for col_idx, (col_name, val) in enumerate(row.items()):
                        if pd.isna(val):
                            values.append('nan')
                        else:
                            try:
                                # Use column-specific format if available, otherwise use default
                                # If col_name is numeric but col_formats has string keys, try to match by position
                                format_str = None
                                if col_name in col_formats:
                                    format_str = col_formats[col_name]
                                elif isinstance(col_name, int) and col_formats:
                                    # Try to get format by position if columns are numeric
                                    col_list = list(col_formats.keys())
                                    if col_idx < len(col_list):
                                        format_str = col_formats[col_list[col_idx]]
                                
                                if format_str:
                                    values.append(format_str % float(val))
                                else:
                                    values.append(float_format % float(val))
                            except Exception:
                                values.append(str(val))
                    
                    # Build line
                    line = date_str + ' ' + ' '.join(values)
                    
                    # Add comment if present
                    row_pos = out_df.index.get_loc(idx)
                    if row_pos < len(comments_col) and pd.notna(comments_col.iloc[row_pos]):
                        line = line + ' ' + comments_col.iloc[row_pos]
                    
                    outfile.write(line + '\n')
            else:
                out_df.to_csv(
                    outfile,
                    sep=" ",
                    header=include_header,
                    float_format=float_format,
                    date_format=time_format,
                )
    finally:
        if should_close:
            outfile.close()


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
