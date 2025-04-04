#!/usr/bin/env python
import os.path
import os
import sys
import pandas as pd
import geopandas as gpd
import shapefile
from shapely.geometry import Point
from vtools.functions.unit_conversions import *
from dms_datastore.dstore_config import *


if sys.version_info[0] < 3:
    from pandas.compat import u
    from builtins import open, file, str
else:
    u = lambda x: x

import argparse
from vtools.data.timeseries import *

station_variables = [
    "elev",
    "air pressure",
    "wind_x",
    "wind_y",
    "temp",
    "salt",
    "u",
    "v",
    "w",
    "ssc",
]


def staout_name(var):
    try:
        ndx = station_variables.index(var)
        return "staout_{}".format(ndx + 1)
    except:
        raise ValueError(
            "Input variable is not standard station variable: {}".format(var)
        )


def read_staout(
    fname,
    station_infile,
    reftime,
    ret_station_in=False,
    multi=False,
    elim_default=False,
    time_unit="s",
):
    """Read a SCHISM staout_* file into a pandas DataFrame

    Parameters
    ----------
    fpath : fname
        Path to input staout file or a variable name in ["elev", "air pressure", "wind_x", "wind_y",  "temp", "salt", "u", "v", "w"] whose
        1-index will be mapped to a name like staout_1 for elev

    station_infile : str or DataFrame
        Path to station.in file or DataFrame from read_station_in

    reftime : Timestampe
        Start of simulation, time basis for staout file elapse time

    ret_station_in : bool
        Return station_in DataFrame for use, which may speed reading of a second file

    multi : bool
        Should the returned data have a multi index for the column with location and sublocation. If False the two are collapsed

    elim_default : bool
        If the MultiIndex is collapsed, stations with subloc "default" will be collapsed. Eg. ("CLC","default") becomes "CLC_default"

    time_unit : string
        Convertible to pandas frequency string, this is the timestamp of the file.

     Returns
     -------
     Result : DataFrame
         DataFrame with hierarchical index (id,subloc) and columns representing the staout data (collapsed as described above

    Examples
    --------

    >>> staout1,station_in = read_staout("staout_1","station.in",reftime=pd.Timestamp(2009,2,10),
                                 ret_station_in = True,multi=False,elim_default=True)
    >>> staout6 = read_staout("staout_6",station_in,reftime=pd.Timestamp(2009,2,10),multi=False,elim_default=True)

    """

    if isinstance(station_infile, str):
        station_in = read_station_in(station_infile)
    else:
        station_in = station_infile

    station_index = station_in.index.copy()

    if station_index.duplicated().any():
        print(station_index[station_index.duplicated()])
        raise ValueError(
            "Duplicate id/subloc pair in station.in file {}".format(station_infile)
        )
    staout = pd.read_csv(fname, index_col=0, sep="\s+", header=None)
    # todo: hardwire
    staout.mask(staout <= -999.0, inplace=True)
    staout.columns = station_index
    elapsed_datetime(staout, reftime=reftime, inplace=True, time_unit=time_unit)
    staout.index = staout.index.round("s")
    if not multi:
        if elim_default:
            staout.columns = [
                f"{loc}_{subloc}" if subloc != "default" else f"{loc}"
                for loc, subloc in staout.columns
            ]
        else:
            staout.columns = [f"{loc}_{subloc}" for loc, subloc in staout.columns]
    f = pd.infer_freq(staout.index)
    if f is None:
        # raise ValueError("Could not determine the time frequency of staoutfile")
        f2 = pd.infer_freq(staout.iloc[0:10, 0].index)
        newindex = pd.date_range(staout.index[0], freq=f2, periods=len(staout.index))
        staout.index = newindex
    else:
        staout = staout.resample(f).asfreq(f)
    return (staout, station_infile) if ret_station_in else staout


def convert_stations(input, output, write_sta_var="all"):
    """Read a station shapefile/.in file and write to a .in/shapefile

    Parameters
    ----------
    input : fname
       Path to input station.in style file or station.shp style file with station id, x, y, z, name and subloc fields

    output : fname
       Path to input station.in style file or station.shp style file with station id, x, y, z, name and subloc fields

    write_sta_var : 'all' or list(str)
       List of variables to put in output request from the choices 'elev', 'air pressure', 'wind_x', 'wind_y', 'temp', 'salt', 'u', 'v', 'w'
       or 'all' to include them all

    Returns
    -------
    Result : DataFrame
        DataFrame with hierarchical index (id,subloc) and columns x,y,z,name

    """
    stations = read_pointstrings(input)
    write_pointstrings(output, stations, request=write_sta_var)
    print(f"\tInput: {input} converted and written to Output: {output}")


def read_pointstrings(fpath):
    if fpath.endswith(".in"):
        return read_station_in(fpath)
    elif fpath.endswith(".shp"):
        return read_station_shp(fpath)
    else:
        raise ValueError("Not supported file type")


def write_pointstrings(fpath, station_in, request="all"):
    if fpath.endswith(".in"):
        return write_station_in(fpath, station_in, request=request)
    elif fpath.endswith(".shp"):
        return write_station_shp(fpath, station_in)
    else:
        raise ValueError("Not supported file type")


def read_station_shp(fpath, pop_xy=True):
    """Read a shapefile and convert into a pandas DataFrame

    Parameters
    ----------
    fpath : fname
       Path to input point shapefile - has station id, x, y, z, name and subloc labels (id is the station id, index will be autogenerated)

    pop_xy : bool
       Repopulate the x & y fields with point coordinates?

    Returns
    -------
    Result : DataFrame
        DataFrame that has station id, x, y, z, name and subloc labels (id is the station id, index will be autogenerated)

    """

    if os.path.exists(fpath):
        sf = shapefile.Reader(fpath)
        fields = [x[0] for x in sf.fields][1:]
        records = [y[:] for y in sf.records()]
        shps = [s.points[0] for s in sf.shapes()]

        df = pd.DataFrame(columns=fields, data=records)

        if pop_xy:
            for index, row in df.iterrows():
                df.loc[index, "x"] = round(shps[index][0], 2)
                df.loc[index, "y"] = round(shps[index][1], 2)

        return df

    else:
        raise ValueError("File not found")


def write_station_shp(fpath, station_in):
    """Write a point Shapefile file given a pandas DataFrame of metadata

    Parameters
    ----------
    fpath : fname
       Path to output station.in file

    station_in : DataFrame
       DataFrame that has station id, x, y, z, name and subloc labels (id is the station id, index will be autogenerated)
    """

    station_in["geometry"] = station_in.apply(
        lambda x: Point((float(x.x), float(x.y))), axis=1
    )

    gdf = gpd.GeoDataFrame(station_in, geometry="geometry")

    gdf.to_file(fpath, driver="ESRI Shapefile")


def read_station_in(fpath):
    """Read a SCHISM station.in file into a pandas DataFrame

    .. note::
        This only reads the tabular part, and assumes the BayDelta SCHISM format with columns:
        index x y z ! id subloc "Name"

        Note that there is no header and the delimiter is a space. Also note that the text beginning with !
        is extra BayDeltaSCHISM extra metadata, not required for vanilla SCHISM

     Parameters
     ----------
     fpath : fname
        Path to input station.in style file

     Returns
     -------
     Result : DataFrame
         DataFrame with hierarchical index (id,subloc) and columns x,y,z,name

    """

    with open(fpath, "r") as f:
        request = f.readline()
        n_entry = f.readline()
        stations = pd.read_csv(
            f,
            sep="\s+",
            header=None,
            names=["index", "x", "y", "z", "excl", "id", "subloc", "name"],
            usecols=["x", "y", "z", "id", "subloc", "name"],
            index_col=["id", "subloc"],
            na_values="-",
            keep_default_na=True,
        )
    return stations


def write_station_in(fpath, station_in, request=None):
    """Write a SCHISM station.in file given a pandas DataFrame of metadata

    Parameters
    ----------
    fpath : fname
       Path to output station.in file

    station_in : DataFrame
       DataFrame that has station id, x, y, z, name and subloc labels (id is the station id, index will be autogenerated)

    request :  'all' or list(str)
       List of variables to put in output request from the choices 'elev', 'air pressure', 'wind_x', 'wind_y', 'temp', 'salt', 'u', 'v', 'w'
       or 'all' to include them all
    """
    request_int = [0] * len(station_variables)
    if request == "all":
        request = ["all"]
    request_str = station_variables if request[0] == "all" else request
    request_int = [(1 if var in request_str else 0) for var in station_variables]
    request_explain = " !" + ",".join(station_variables)
    dfmerged = station_in.reset_index()
    dfmerged.index += 1
    dfmerged["excl"] = "!"
    nitem = len(dfmerged)
    # First two lines are a space delimited 1 or 0 for each request then the
    # total number of station requests
    buffer = " ".join([str(x) for x in request_int]) + "{}\n{}\n".format(
        request_explain, nitem
    )
    # Then the specific requests, here written to a string buffer
    buffer2 = dfmerged.to_csv(
        None,
        columns=["x", "y", "z", "excl", "id", "subloc", "name"],
        index_label="id",
        sep=" ",
        float_format="%.2f",
        header=False,
    )
    with open(fpath, "w", newline="") as f:
        f.write(buffer)
        f.write(u(buffer2))
        # f.write(u(buffer))
        # f.write(u(buffer2))


def read_station_subloc(fpath):
    """Read a BayDeltaSCHISM station_sublocs.csv  file into a pandas DataFrame

      The BayDelta SCHISM format has a header and uses "," as the delimiter and has these columns:
      id,subloc,z

      The id is the station id, which is the key that joins this file to the station database. 'subloc' is a label that describes
      the sublocation or subloc and z is the actual elevation of the instrument

      Example might be:
      id,subloc,z
      12345,upper,-0.5

      Other columns are allowed, but this will commonly merged with the station database file so we avoid column names like 'name' that might collide

    Parameters
    ----------
    fpath : fname
       Path to input station.in style file

    Returns
    -------
    Result : DataFrame
        DataFrame with hierarchical index (id,subloc) and data column z

    """

    df = pd.read_csv(fpath, sep=",", header=0, index_col=["id", "subloc"], comment="#")
    df["z"] = df.z
    return df[["z"]]


def read_station_dbase(fpath):
    """Read a BayDeltaSCHISM station data base csv  file into a pandas DataFrame

      The BayDelta SCHISM format is open, but expects these columns:
      index x y z ! id subloc "Name"


    Parameters
    ----------
    fpath : fname
       Path to input dbase style file

    Returns
    -------
    Result : DataFrame
        DataFrame with hierarchical index (id,subloc) and columns x,y,z,name

    """

    db = pd.read_csv(
        fpath, sep=",", comment="#", header=0, index_col="id", dtype={"agency_id": str}
    )
    db["agency_id"] = db["agency_id"].str.replace("'", "", regex=True)

    dup = db.index.duplicated()
    db.index = db.index.str.replace("'", "")
    if dup.sum(axis=0) > 0:
        print("Duplicates")
        print(db[dup])
        raise ValueError("Station database has duplicate id keys. See above")
    return db


def merge_station_subloc(station_dbase, station_subloc, default_z):
    """Merge BayDeltaSCHISM station database with subloc file, producing the union of all stations and sublocs including a default entry for stations with no subloc entry

    Parameters
    ----------
    station_dbase : DataFrame
       This should be the input that has only the station id as an index and includes other metadata like x,y,

    station_subloc : DataFrame
       This should have (id,subloc) as an index

    Returns
    -------
    Result : DataFrame
        DataFrame that links the information.

    """

    merged = station_dbase.reset_index().merge(
        station_subloc.reset_index(), left_on="id", right_on="id", how="left"
    )
    merged.fillna({"subloc": "default", "z": default_z}, inplace=True)
    merged.set_index(["id", "subloc"], inplace=True)

    return merged


def read_obs_links(fpath):
    """Read an obs_links csv file which has comma as delimiter and (id,subloc,variable) as index"""
    df = pd.read_csv(
        fpath,
        sep=",",
        header=0,
        index_col=["station_id", "subloc", "variable"],
        comment="#",
    )
    df.index = df.index.set_levels(df.index.levels[0].astype(str).str.lower(), level=0)
    dups = df.index.duplicated(keep="last")
    ndup = dups.sum(axis=0)
    if ndup > 0:
        print("{} duplicate index rows:".format(ndup))
        print(df.loc[dups, :])
        raise ValueError(
            "Duplicates not allowed in observation links file, see above list."
        )
    return df


def read_station_out(fpath_base, stationinfo, var=None, start=None):
    if var is None:
        fname = fpath_base
    else:
        try:
            fileno = station_variables.index(var)
        except ValueError:
            raise ValueError(
                "Variable name {} not on list: {}.format(var,station_variables"
            )
        fname = "{}_{:d}".format(fpath_base, fileno)
    data = pd.read_csv(
        fname,
        var,
        sep="\s+",
        index_col=0,
        header=None,
        names=stationinfo.index,
        dtype="d",
    )
    if start is not None:
        data = elapsed_datetime(data, reftime=start)
    f = pd.infer_freq(data.index)
    data = data.asfreq(f)
    return data


def flux_stations_from_yaml(inp):
    """Retrieve station id of fluxlines from yaml file or content"""
    import yaml

    if os.path.exists(inp):
        with open(inp) as f:
            content = yaml.full_load(f)
    else:
        content = inp
        if not "linestrings" in content:
            raise ValueError(
                "Could not fine 'linestrings' key. Was the input a string or filename? Valid file? For FLOW, make sure there is 'flow_station_input' parameter"
            )
    names = []
    linestrings = content["linestrings"]
    for ls in linestrings:
        if ls.get("Name") is not None:
            raise ValueError(
                "YAML file is case sensitive. Do not use Name as attribute label"
            )
        name = ls.get("name")
        if name in names:
            print("Duplicate name: {}".format(name))
        if not (name == name.lower()):
            # print("Station ids should be lower case, coercing")
            name = name.lower()
        names.append(name)
    return names


def station_names_from_file(fpath):
    ext = os.path.splitext(fpath)[1]
    if ext in (".yml", ".yaml"):
        station_names = flux_stations_from_yaml(fpath)
    elif ext == ".prop":
        station_names = []
        with open(fpath, "r") as f:
            for line in f:
                names = line.strip().split()
                if len(names) == 1:
                    station_names.append(names[0])
    else:
        raise ValueError(
            f"File type not recognized for harvesting station names: {fpath}"
        )

    return station_names


def read_flux_out(fpath, names, reftime):
    """Read fluxes from a SCHISM flux.out file

    Parameters
    ----------

    fpath : str
    Path to the file

    names : str
    name of file that contains  names of  flux areas, typically something like flow_xsects.yaml

    reftime : str
    start of simulation, against which relative times will be calculated
    """

    if isinstance(names, str):
        names = station_names_from_file(names)

    # Check uniqueness of names
    seen = set()
    uniq = []
    for x in names:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    if len(uniq) != len(names):
        raise ValueError("Duplicate station names.")
    names = [x.lower() for x in names]
    nstation = len(names)
    # probe = pd.read_csv(fpath,sep="\s+",index_col=0,header=None,dtype='d',nrows=2)
    # ncolfile = probe.shape[1]

    cols = [0, *(range(1, nstation + 1))]
    data = pd.read_csv(
        fpath,
        sep="\s+",
        index_col=0,
        header=None,
        usecols=cols,
        names=["time"] + names,
        dtype="d",
    )
    if reftime is not None:
        data = elapsed_datetime(data, reftime=reftime, time_unit="d")
        data.index = data.index.round(freq="s")
        f = pd.infer_freq(data.index)
        data = data.asfreq(f)

    # todo: freq when start is none?
    return data


def example():
    print(read_station_in("example_station.in"))
    stations_utm = read_station_dbase("stations_utm.csv")
    print(stations_utm)
    ssubloc = read_station_subloc("station_subloc.csv")
    stations_in = merge_station_subloc(stations_utm, ssubloc, default_z=-0.5)
    # stations_in = pd.merge(stations_utm,ssubloc,how='inner',left_index=True,right_index=True)
    # print(stations_in)
    station_request = ["salt", "elev"]
    write_station_in("station.in", stations_in, request=station_request)
    # stations_in = read_station_in("station.in")
    obs_links = read_obs_links("obs_links.csv")
    merged = stations_in.merge(obs_links, left_index=True, right_index=True, how="left")

    if True:
        print("**")
        print(obs_links)
        print("**")
        print(stations_in)
        print("**")
        print(stations_utm)
        print("**")
        print(merged)


uconversions = {"ft": m_to_ft, "ec": psu_ec_25c, "cfs": cms_to_cfs}


def station_subset(
    fpath,
    run_start,
    locs,
    extract_freq,
    convert=None,
    stationfile=None,
    isflux="infer",
    miss="raise",
):
    """Extract a subset of stations from an staout file or flux.out file

    Parameters
    ----------
    fpath : str
       Path to the output file to be read

    run_start : pd.Timestamp
       Start time (reference time) for the simulation elapsed time

    locs : pd.DataFrame or str
       A DataFrame with rows specifying the data to subset or a string that is the path to such a file.
       There are a few options. Minimally this file should have either one column called "station_id" or
       one called "id" and another called "subloc". If you use station_id, it should be an
       underscore-connected combination of id and subloc which should be unique, and this will be the treatment of the output.
       of the station. If you use "subloc" you can use "default" leave the subloc column blank. You can also use another optional
       column called "alias" and this will become the label used.

     extract_freq : str or pd.tseries.TimeOffset
       Frequency to extract ... this allows some economies if you want, say, 15min data. Use pandas freq string such as '15T' for 15min

     convert: str or function
        A small number of conversions are supported ("ft", "ec" for uS/cm, "cfs" for flow)

     stationfile : str
        Name of station file such as station.in. In the case of flow this will be a yaml file or fluxflag.prop file produced by our preprocessing system.
        If you leave this None, 'station.in' in the same directory as the output file will be assumed for staout files. For flow, a string must be supplied
        but will be tested first in the directory of execution and then side-by-side in that order.

     isflux : 'infer' | True | False
        Is the request for flux.out?

     miss : 'raise' | 'drop' | 'nan'
         What to do when a requested station_id does not exist. The default, raise, helps station lists from growing faulty. 'drop' will ignore the column and 'nan'
         will return nans for the column.

    Returns
    -------
    Result : DataFrame
        DataFrame that returns the converted and subsetted data. Column names will be the ids unless 'alias' is provided in locs, in which case those names will be
        swapped in.

    """

    locs.station_id = locs.station_id.str.lower()
    locs = locs.set_index("station_id")

    if isflux == "infer":
        if "staout" in fpath:
            isflux = False
        elif "flux" in fpath:
            isflux = True
        else:
            raise ValueError("Station output type (flux,staout) could not be inferred")

    if isflux:
        if not os.path.exists(stationfile):
            stationfile = os.path.join(os.path.split(fpath)[0], stationfile)
        staout = read_flux_out(fpath, stationfile, reftime=run_start)
    else:
        if stationfile is None:
            stationfile = os.path.join(os.path.split(fpath)[0], "station.in")
        print(stationfile)
        staout = read_staout(
            fpath,
            stationfile,
            reftime=run_start,
            ret_station_in=False,
            multi=False,
            elim_default=True,
        )
    staout.columns = [x.lower() for x in staout.columns]
    if extract_freq is not None:
        staout = staout.resample(extract_freq).interpolate()
    loc_ndx = locs.index
    missing = [label for label in locs.index if not label in staout.columns]
    if len(missing) > 0:
        if miss == "raise":
            print("Labels not in dataset:")
            print(missing)
            raise ValueError("Requested labels not in dataset")
        elif miss == "nan":
            subset = loc_ndx.values
        elif miss == "drop":
            subset = [label for label in locs.index if label in staout.columns]
        else:
            raise ValueError("miss must be one of raise | drop | nan")
    else:
        subset = locs.index.values

    def find_nondups(columns):
        nondups = list()
        cset = set(columns)
        for item in columns:
            try:
                cset.remove(item)
                nondups.append(True)
            except:
                print(
                    "Column name {} is duplcated in outputs, accepting first instance".format(
                        item
                    )
                )
                nondups.append(False)
                # raise ValueError("Column name {} is duplcated".format(item))
        return nondups

    nondups = find_nondups(staout.columns)
    staout = staout.loc[:, nondups]

    sub_df = staout.reindex(subset, axis="columns")
    if "alias" in locs.columns:
        sub_df.columns = locs.alias.values
    if convert is not None:
        cv = uconversions[convert] if isinstance(convert, str) else convert
        sub_df = cv(sub_df)
    return sub_df


def station_subset_multidir(
    dirs,
    staoutfile,
    run_start,
    locs,
    extract_freq,
    convert,
    stationfile=None,
    names=None,
    isflux="infer",
    miss="raise",
):
    """Extract a subset of stations from an staout file or flux.out file across a list of directories

    Parameters
    ----------

    dirs : list(str)
        List of directories. The output dataframe will have a column multindex (dir,station_id) where dir is the directory of the output.


    fpath : str
       Path to the output file to be read

    run_start : pd.Timestamp
       Start time (reference time) for the simulation elapsed time

    locs : pd.DataFrame or str
       A DataFrame with rows specifying the data to subset or a string that is the path to such a file.
       There are a few options for staout station files. Minimally this file should have either one column called "station_id" or
       one called "id" and another called "subloc". If you use station_id, it should be an
       underscore-connected combination of id and subloc and this will be the index treatment of the output.
       Flux files only have the station_id option. If you use "subloc" you can use "default" leave the subloc column blank.
       You can also include another optional column called "alias" and this will become the label used.

     extract_freq : str or pd.tseries.TimeOffset
       Frequency to extract ... this allows some economies if you want, say, 15min data. Use pandas freq string such as '15T' for 15min

     convert: str or function
        A small number of conversions are supported ("ft", "ec" for uS/cm, "cfs" for flow)

     stationfile : str
        Name or list of station file such as station.in. You can provide a list of the same length as dirs or a single value which will be assumed
        to be appropriate for all the directories. In the case of station.in, you can use None and 'station.in' in each directory
        will be assumed for staout files. For flow, a string (yaml file) must be supplied but will be tested first in the directory of
        execution and then side-by-side in that order.

     isflux : 'infer' | True | False
        Is the request for flux.out?

     miss : 'raise' | 'drop' | 'nan'
         What to do when a requested station_id does not exist. The default, raise, helps station lists from growing faulty. 'drop' will ignore the column and 'nan'
         will return nans for the column.

    Returns
    -------
    Result : DataFrame
        DataFrame that returns the converted and subsetted data. Column names will be the ids unless 'alias' is provided in locs, in which case those names will be
        swapped in.

    """

    dfs = []
    if stationfile is None or isinstance(stationfile, str):
        stationfile = [stationfile] * len(dirs)
    else:
        stationfile = list(stationfile)
        if len(stationfile) != len(dirs):
            raise ValueError(
                "stationfile must be same length as dirs if it is iterable"
            )
    for d, s in zip(dirs, stationfile):
        fpath = os.path.join(d, staoutfile)
        dfs.append(
            station_subset(
                fpath,
                run_start,
                locs,
                extract_freq,
                convert,
                stationfile=s,
                isflux=isflux,
                miss=miss,
            )
        )
    if names is None:
        names = dirs
    out = pd.concat(dfs, keys=names, names=["sim", "loc"], axis=1)
    return out


def convert_db_station_in(
    outfile="station.in",
    stationdb=None,
    sublocdb=None,
    station_request="all",
    default=-0.5,
):
    stations_utm = read_station_dbase(stationdb)
    ssubloc = read_station_subloc(sublocdb)
    stations_in = merge_station_subloc(stations_utm, ssubloc, default_z=-0.5)
    write_station_in(outfile, stations_in, request=station_request)


def create_arg_parser():
    """Create an argument parser"""
    parser = argparse.ArgumentParser(
        description="Create station.in file from station database (stations_utm.csv) and station subloc listing station_subloc.csv"
    )
    parser.add_argument(
        "--station_db",
        default=None,
        help="station database, otherwise station_dbase as configured in dms_datastore dstore_config file",
    )
    parser.add_argument(
        "--subloc_db",
        default=None,
        help="subloc listings for stations (otherwise default subloc from dms_datastore dstore_config file)",
    )
    parser.add_argument(
        "--request",
        default="all",
        nargs="+",
        help="requested variables or 'all' for all of them. Possibilities are: {}".format(
            ",".join(station_variables)
        ),
    )
    parser.add_argument(
        "--default_zcor",
        default="-0.5",
        help="z coordinate used when there is no listing for station id (z coordinate, not subloc from surface)",
    )
    parser.add_argument("--out", default="station.in", help="station.in formatted file")
    return parser


def main():
    """A main function to convert polygon files"""
    parser = create_arg_parser()
    args = parser.parse_args()
    stationdb = args.station_db
    sublocdb = args.subloc_db
    default = args.default_zcor
    request = args.request
    outfile = args.out
    if stationdb is None:
        stationdb = config_file("station_dbase")
    if sublocdb is None:
        sublocdb = config_file("sublocations")

    convert_db_station_in(outfile, stationdb, sublocdb, request, default)


if __name__ == "__main__":
    # example()
    main()
