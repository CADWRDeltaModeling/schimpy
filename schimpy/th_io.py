"""I/O utilities for SCHISM .th time-history files.

Provides `read_th` and `write_th` as a shared module for model_time, th_calcs,
and merge_th.
"""

__all__ = [
    "is_elapsed",
    "get_column_formats",
    "read_th",
    "write_th",
]

import datetime
import os
import re
import tempfile

import pandas as pd


def is_elapsed(input_path):
    """Return True if .th file is elapsed seconds format.

    Parameters
    ----------
    input_path : str
        Path to a .th file.

    Returns
    -------
    bool
        True if the first non-empty, non-comment line is elapsed numeric time.
    """
    with open(input_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            line = line.split("#", 1)[0].strip()
            is_elapsed_line = all(
                re.match(r"^[-+]?\d*\.?\d+$", s) for s in line.split() if s
            )
            return is_elapsed_line
    return False


def get_column_formats(filename):
    """Extract per-column format strings (e.g., "%.3f") from .th input.

    Parameters
    ----------
    filename : str
        Path to a .th file.

    Returns
    -------
    dict
        Mapping column names to printf-style format specifiers.
    """
    formats = {}
    with open(filename, "r") as f:
        header = None
        data_lines = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if header is None and not line.startswith("datetime") and not line.startswith("date") and not line.startswith("time"):
                continue
            if header is None and (line.startswith("datetime") or line.startswith("date") or line.startswith("time")):
                header = line.split()
                continue
            if header is not None and line.strip() and not line.lstrip().startswith("#"):
                data_line = line.split("#")[0].rstrip() if "#" in line else line
                data_lines.append(data_line)

        if header is None or len(data_lines) == 0:
            return formats

        for col_idx, col_name in enumerate(header[1:], start=1):
            decimals_list = []
            for data_line in data_lines:
                parts = data_line.split()
                if col_idx < len(parts):
                    value_str = parts[col_idx]
                    if "." in value_str:
                        decimals = len(value_str.split(".")[1])
                        decimals_list.append(decimals)
                    else:
                        decimals_list.append(0)
            if decimals_list:
                max_decimals = max(decimals_list)
                formats[col_name] = f"%.{max_decimals}f"
    return formats


def read_th(input_path, time_basis=None, to_timestamp=True, elapsed_unit="s", head=None):
    """Read a SCHISM .th file and return a pandas DataFrame.

    Supports both timestamped (datetime index) and elapsed (seconds index) files.
    Preserves inline or trailing comments in a ``__comment__`` column.

    Parameters
    ----------
    input_path : str
        Path to .th file.
    time_basis : str or datetime, optional
        Reference time for converting elapsed series to timestamps.
    to_timestamp : bool, default True
        Convert elapsed data to timestamp index if True.
    elapsed_unit : str, default "s"
        Time units for elapsed files.
    head : str or None
        Optional header file to supply data column names for elapsed data.

    Returns
    -------
    pandas.DataFrame
        Time-series dataset with datetime index.
    """
    file_is_elapsed = is_elapsed(input_path)

    if file_is_elapsed:
        if time_basis is None:
            raise ValueError(
                f"The input file {input_path} is an elapsed .th file and there's no time_basis specified."
            )

        if to_timestamp:
            from schimpy.model_time import elapsed_to_timestamp

            out_df = elapsed_to_timestamp(input_path, time_basis, elapsed_unit=elapsed_unit)
        else:
            from schimpy.model_time import read_elapsed

            out_df = read_elapsed(input_path)

        if head is not None:
            from schimpy.model_time import get_headers

            headers = get_headers(head)
            num_data_cols = len([c for c in out_df.columns if c != "__comment__"])
            if len(headers) != num_data_cols:
                raise ValueError(
                    f"{os.path.basename(input_path)} has {num_data_cols} columns and the header calls for {len(headers)}!"
                )
            has_comment_col = "__comment__" in out_df.columns
            if has_comment_col:
                comment_col = out_df.pop("__comment__")
            out_df.columns = headers
            if has_comment_col:
                comment_col.index = out_df.index
                out_df["__comment__"] = comment_col

        out_df.index.name = "datetime"

    else:
        comments_dict = {}
        with open(input_path, "r") as f:
            data_lines = []
            header_line = None
            last_data_idx = None
            for line in f:
                line = line.rstrip("\n")
                if header_line is None and not line.startswith("#") and not line.strip().startswith("#"):
                    if line and not any(c.isdigit() for c in line.split()[0]) if line.split() else False:
                        header_line = line
                        data_lines.append(line)
                        continue
                if line.strip() and line.lstrip().startswith("#"):
                    if last_data_idx is not None:
                        comments_dict[last_data_idx] = line.strip()
                    continue
                if line.startswith("#"):
                    continue
                if line.strip():
                    if "#" in line:
                        data_part, comment_part = line.split("#", 1)
                        data_lines.append(data_part.rstrip())
                        last_data_idx = len(data_lines) - 1
                        comments_dict[last_data_idx] = "# " + comment_part.lstrip()
                    else:
                        data_lines.append(line)
                        last_data_idx = len(data_lines) - 1

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tmp", delete=False) as tmp:
            for line in data_lines:
                tmp.write(line + "\n")
            tmp_path = tmp.name

        try:
            out_df = pd.read_table(tmp_path, sep=r"\s+", index_col="datetime")
            out_df.index = pd.to_datetime(out_df.index, format="%Y-%m-%dT%H:%M")
            if comments_dict:
                out_df["__comment__"] = None
                for idx, comment in comments_dict.items():
                    if idx > 0 and idx - 1 < len(out_df):
                        out_df.iloc[idx - 1, out_df.columns.get_loc("__comment__")] = comment
        finally:
            os.remove(tmp_path)

    if not out_df.index.is_monotonic_increasing:
        if file_is_elapsed:
            bad_rows = out_df.loc[out_df.index.to_series().diff() <= 0.0]
            raise ValueError(f"Non-monotonic elapsed index found:\n{bad_rows}")
        else:
            bad_rows = out_df.loc[out_df.index.to_series().diff() <= pd.to_timedelta("0 seconds")]
            raise ValueError(f"Non-monotonic datetime index found:\n{bad_rows}")

    return out_df


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
    """Write time history data to a SCHISM .th file.

    Parameters
    ----------
    df : pandas.DataFrame or pandas.Series
        Time-indexed input data.
    outpath : str or file-like
        Output path or file-like object.
    elapsed : bool or None, optional
        Force elapsed format (True) or dated format (False). If None, inferred from index.
    ref_time : str or datetime, optional
        Reference datetime for elapsed<->timestamp conversions.
    time_format : str
        Datetime format string for dated output.
    float_format : str
        Number formatting for data values.
    elapsed_format : str
        Number formatting for elapsed timestamps.
    include_header : bool or None, optional
        Include header row. Defaults to True for dated output.
    index_name : str, optional
        Index name for output time column.
    source_file : str, optional
        Read source comments and formatting from this file.

    Returns
    -------
    None
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

    comments = []
    col_formats = {}
    if source_file and isinstance(source_file, str):
        try:
            with open(source_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        comments.append(line)
                    elif line.strip() and not line.startswith("date") and not line.startswith("time"):
                        break
            col_formats = get_column_formats(source_file)
        except Exception:
            pass

    should_close = False
    if isinstance(outpath, str):
        outfile = open(outpath, mode="w", newline="\n")
        should_close = True
    else:
        outfile = outpath

    try:
        for comment in comments:
            outfile.write(comment)

        if elapsed:
            data_cols = [c for c in df_local.columns if c != "__comment__"]
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

            has_comments = "__comment__" in df_local.columns
            if has_comments:
                comments_col = df_local.pop("__comment__")

            if include_header:
                header_df = df_local.copy()
                header_df.index = elapsed_index
                header_df.index.name = index_name or "time"

                if has_comments:
                    comments_col = header_df.pop("__comment__")
                    header_str = " ".join(header_df.columns)
                    outfile.write(header_df.index.name + " " + header_str + "\n")
                    for i, (idx, row) in enumerate(header_df.iterrows()):
                        tstr = elapsed_format % idx
                        values = []
                        for col_idx, (col_name, val) in enumerate(row.items()):
                            if pd.isna(val):
                                values.append("nan")
                            else:
                                try:
                                    format_str = None
                                    if col_name in col_formats:
                                        format_str = col_formats[col_name]
                                    elif isinstance(col_name, int) and col_formats:
                                        col_list = list(col_formats.keys())
                                        if col_idx < len(col_list):
                                            format_str = col_formats[col_list[col_idx]]
                                    if format_str:
                                        values.append(format_str % float(val))
                                    else:
                                        values.append(float_format % float(val))
                                except Exception:
                                    values.append(str(val))
                        line = tstr + " " + " ".join(values)
                        if i < len(comments_col) and pd.notna(comments_col.iloc[i]):
                            line = line + " " + comments_col.iloc[i]
                        outfile.write(line + "\n")
                else:
                    header_df.to_csv(outfile, sep=" ", header=True, float_format=float_format)
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
                        line = line + " " + comments_col.iloc[i]
                    outfile.write(line + "\n")
        else:
            if is_datetime_index:
                out_df = df_local.copy()
            else:
                if ref_time is None:
                    raise ValueError("ref_time is required to write dated output")
                out_df = df_local.copy()
                out_df.index = pd.to_datetime(ref_time) + pd.to_timedelta(df_local.index.astype(float), unit="s")

            out_df.index.name = index_name or "datetime"

            has_comments = "__comment__" in out_df.columns
            if has_comments:
                comments_col = out_df.pop("__comment__")
                if include_header:
                    header_str = " ".join(out_df.columns)
                    outfile.write(out_df.index.name + " " + header_str + "\n")
                for idx, row in out_df.iterrows():
                    date_str = idx.strftime(time_format)
                    values = []
                    for col_idx, (col_name, val) in enumerate(row.items()):
                        if pd.isna(val):
                            values.append("nan")
                        else:
                            try:
                                format_str = None
                                if col_name in col_formats:
                                    format_str = col_formats[col_name]
                                elif isinstance(col_name, int) and col_formats:
                                    col_list = list(col_formats.keys())
                                    if col_idx < len(col_list):
                                        format_str = col_formats[col_list[col_idx]]
                                if format_str:
                                    values.append(format_str % float(val))
                                else:
                                    values.append(float_format % float(val))
                            except Exception:
                                values.append(str(val))
                    line = date_str + " " + " ".join(values)
                    row_pos = out_df.index.get_loc(idx)
                    if row_pos < len(comments_col) and pd.notna(comments_col.iloc[row_pos]):
                        line = line + " " + comments_col.iloc[row_pos]
                    outfile.write(line + "\n")
            else:
                out_df.to_csv(outfile, sep=" ", header=include_header, float_format=float_format, date_format=time_format)
    finally:
        if should_close:
            outfile.close()
