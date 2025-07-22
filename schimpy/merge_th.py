import argparse

import pandas as pd
import numpy as np

from vtools.data.timeseries import datetime_elapsed
from schimpy import schism_yaml as syml


role_to_sign = {"source_tracer": 1, "source_flux": 1, "sink_flux": -1}
ambient = -9999


def read_th_spec(fname):
    with open(fname) as stream:
        try:
            return syml.load(stream)
        except yaml.YAMLError as exc:
            print("Fail")
            print(exc)


def read_locations_from_input(fname, role):
    with open(fname, "r") as stream:
        loc_yaml = syml.load(stream)
        if role == "source_flux" or role == "source_tracer":
            if "sources_sinks" in loc_yaml:
                section = loc_yaml["sources_sinks"]["sources"]
            elif "sources" in loc_yaml:
                section = loc_yaml["sources"]
            else:
                raise ValueError(
                    "sources section not found in file or in sources_sinks parent section"
                )
        elif role == "sink_flux":
            if "sources_sinks" in loc_yaml:
                section = loc_yaml["sources_sinks"]["sinks"]
            elif "sinks" in loc_yaml:
                section = loc_yaml["sinks"]
            else:
                raise ValueError(
                    "sinks section not found in file or in sources_sinks parent section"
                )

        return [item[0] for item in section.items()]
        # Assume sourc


def write_th(df, fname, elapsed, ref_time=None):
    print("writing ", fname)
    print(f"dataframe dimensions: {df.shape}")
    if elapsed:
        dfe = datetime_elapsed(df, reftime=ref_time, dtype=float)
        dfe.index.name = "time"
        with open(fname, mode="w", newline="\n") as outfile:
            dfe.to_csv(outfile, sep=" ", float_format="%.3f", header=False)
            # outfile.write("#"+str(dfe.columns.values.join(","))+"\n")
    else:
        sep = " " if fname.endswith("th") else ","
        df.to_csv(fname, sep=sep, float_format="%.3f", date_format="%Y-%m-%dT%H:%M")


def read_data(fname, variable):
    print(f"reading file: {fname}")
    if fname.endswith(".th"):
        sep = "\s+"
        comment = "#"
    else:
        sep = ","
        comment = "#"

    if variable is None:
        dfnew = pd.read_csv(
            fname,
            header=[0],
            sep=sep,
            index_col=0,
            parse_dates=[0],
            dtype=float,
            comment=comment,
        )
        return dfnew

    if variable.startswith("multi"):
        # Read directly as a multi index, appropriate for files with multiple variables
        dfnew = pd.read_csv(
            fname,
            header=[0, 1],
            sep=sep,
            index_col=0,
            parse_dates=[0],
            dtype=float,
            comment="#",
        )
    else:
        dfnew = pd.read_csv(
            fname, header=[0], sep=sep, index_col=0, parse_dates=[0], dtype=float
        )
        # Adds a level called variable to columns and sets it to current variable,
        # Resulting in a multiindex
        dfnew = pd.concat([dfnew], names=["variable"], keys=[variable], axis=1)
    return dfnew


def merge_th(th_spec):
    if isinstance(th_spec, str):
        spec = read_th_spec(th_spec)
    else:
        spec = th_spec

    spec = spec["merge_time_history"]
    start_time = pd.to_datetime(spec["start_time"])
    end_time = pd.to_datetime(spec["end_time"])
    th_files = spec["th_files"]
    vdf = {}

    for th_out in th_files:
        print("processing: ", th_out)
        th_config = th_files[th_out]
        role = th_config["role"]
        dated_output = (
            th_config["dated_output"] if "dated_output" in th_config else False
        )
        sign = role_to_sign[role]
        locfile = th_config["locs"]
        locs = read_locations_from_input(locfile, role)
        dfs = []
        for item in th_config["data"]:
            try:
                fname, variable = item["fname"], item["variable"]

            except:
                fname, variable = item, None

            dfnew = read_data(fname, variable)
            dfs.append(dfnew)

        if "tracer" in role:
            # Create column index representing desired output
            variables = th_config["variables"]  # establishes order and labels
            nvar = len(variables)
            nloc = len(locs)
            varlevel = np.repeat(variables, nloc)
            srclevel = np.tile(locs, nvar)
            multicol = pd.MultiIndex.from_arrays(
                [varlevel, srclevel], names=("variable", "source")
            )
            # mdf = pd.DataFrame(index=vdf.index,data=np.nan,columns=multicol,dtype=float)

            # horizontally stack, eliminate earlier copies
            vdf = pd.concat(dfs, axis=1).loc[start_time:end_time, :]
            vdf = vdf.loc[:, ~vdf.columns.duplicated(keep="last")].copy()

            mdf = vdf.reindex(multicol, axis=1)  # Choose the desired output

            # take care of anything assigned to a constant
            constantval = th_config["constant"] if "constant" in th_config else None
            if constantval is not None:
                vars = mdf.columns.get_level_values("variable")
                locs = mdf.columns.get_level_values(1)
                for var in constantval:
                    for item in constantval[var]:
                        pat = item["col_pattern"]
                        val = item["value"]
                        overwrite = item["overwrite"]
                        if val == "ambient":
                            if "source" in role:
                                val = ambient
                            else:
                                raise ValueError(
                                    "'ambient' keyword only works for sources (not, for instance, for boundaries)"
                                )  # todo
                        colmatch = (vars == var) & locs.str.match(pat)
                        cupdate = mdf.loc[:, colmatch].copy()  # for size and labels
                        cupdate.loc[:, :] = val  # for value
                        mdf.update(cupdate, overwrite=overwrite)

        else:
            mdf = pd.concat(dfs, axis=1).loc[start_time:end_time, :]
            mdf = mdf.loc[:, ~mdf.columns.duplicated(keep="last")].copy()
            mdf = mdf.reindex(locs, axis=1)
        if sign == 1:
            if ((mdf < 0) & (mdf > -9999.0)).any(axis=None):
                raise ValueError(f"In {th_out} variable signs do not meet convention")
        if sign == -1:
            if (mdf > 0).any(axis=None):
                raise ValueError(f"In {th_out} variable signs do not meet convention")

        colmissing = mdf.isnull().sum(axis=0)
        colmissing = colmissing[colmissing > 0]

        if colmissing.any():
            print("Missing data in columns:")
            print(colmissing)
            print(mdf.loc[mdf.isnull().any(axis=1), :])
            mdf.loc[mdf.isnull().any(axis=1), :].to_csv("missing.csv")
            raise ValueError(
                "Missing data. Check for gaps or non-matching labels between locations and time series."
            )

        _write_output(mdf, th_out, dated_output, start_time)


def _write_output(mdf, th_out, dated_output, ref_time):
    if dated_output is False:
        if not (th_out.endswith(".th")):
            raise ValueError("Elapsed file should have .th extension")
        write_th(mdf, th_out, elapsed=True, ref_time=ref_time)
    elif dated_output is True:
        write_th(mdf, th_out, elapsed=False)
    else:
        dated_name = dated_output
        write_th(mdf, th_out, elapsed=True, ref_time=ref_time)
        write_th(mdf, dated_name, elapsed=False)


def create_arg_parser():
    parser = argparse.ArgumentParser(
        description="Merge multicolumn th files to form union of columns and (for tracers) variables."
    )
    parser.add_argument("--input", help="Config file")
    return parser


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    configfile = args.input
    merge_th(configfile)


if __name__ == "__main__":
    merge_th("source_sink_data_test.yaml")
