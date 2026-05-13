"""Command line tool to splice two NetCDF time series of nudging data"""

import click
import numpy as np
import pandas as pd
import xarray as xr
from vtools.functions.merge import ts_splice


def load_nc_series(nc_file, start_date):
    """
    Load NetCDF tracer data and convert to pandas DataFrame.

    Parameters
    ----------
    nc_file : str
        Path to NetCDF file.
    start_date : str
        Start date corresponding to time=0 in the NetCDF file.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame indexed by datetime with flattened spatial dimensions.
    ds : xarray.Dataset
        Original dataset used as template.
    """

    ds = xr.open_dataset(nc_file)

    time_seconds = ds["time"].values
    start = pd.to_datetime(start_date)

    times = start + pd.to_timedelta(time_seconds, unit="s")

    data = ds["tracer_concentration"].values
    shape = data.shape

    flat = data.reshape(shape[0], -1)

    df = pd.DataFrame(flat, index=times)

    return df, ds


def rebuild_dataset(spliced_df, template_ds, start_date):
    """
    Convert spliced DataFrame back into xarray Dataset.

    Parameters
    ----------
    spliced_df : pandas.DataFrame
        Spliced time series data.
    template_ds : xarray.Dataset
        Template dataset for dimensions and metadata.
    start_date : str
        Start date corresponding to time=0.

    Returns
    -------
    xarray.Dataset
        Dataset ready to write to NetCDF.
    """

    start = pd.to_datetime(start_date)

    time_seconds = (spliced_df.index - start).total_seconds().astype(np.float32)

    node = template_ds.dims["node"]
    levels = template_ds.dims["nLevels"]
    one = template_ds.dims["one"]

    data = spliced_df.values.reshape(len(time_seconds), node, levels, one)

    ds_out = xr.Dataset()

    ds_out["time"] = xr.DataArray(
        time_seconds, dims=["time"], attrs=template_ds["time"].attrs
    )

    ds_out["map_to_global_node"] = template_ds["map_to_global_node"]

    ds_out["tracer_concentration"] = xr.DataArray(
        data,
        dims=["time", "node", "nLevels", "one"],
        attrs=template_ds["tracer_concentration"].attrs,
    )

    return ds_out


def splice_two_files(
    file1, file2, start_date1, start_date2, transition_date, transition
):
    """
    Splice two NetCDF time series.

    Parameters
    ----------
    file1 : str
        First NetCDF file.
    file2 : str
        Second NetCDF file.
    start_date1 : str
        Start date of first dataset.
    start_date2 : str
        Start date of second dataset.
    transition_date : str
        Timestamp where splice occurs.
    transition : {'prefer_first', 'prefer_last'}
        Transition rule passed to ts_splice.

    Returns
    -------
    xarray.Dataset
        Spliced dataset.
    """

    df1, ds1 = load_nc_series(file1, start_date1)
    df2, ds2 = load_nc_series(file2, start_date2)

    spliced = ts_splice([df1, df2], transition=transition)

    transition_ts = pd.to_datetime(transition_date)

    if transition == "prefer_first":
        spliced = spliced.loc[:transition_ts].combine_first(df2.loc[transition_ts:])

    if transition == "prefer_last":
        spliced = df1.loc[:transition_ts].combine_first(spliced.loc[transition_ts:])

    ds_out = rebuild_dataset(spliced, ds1, start_date1)

    return ds_out


@click.command()
@click.option("--file1", required=True, help="Path to first NetCDF file.")
@click.option("--file2", required=True, help="Path to second NetCDF file.")
@click.option(
    "--start-date1", required=True, help="Start date corresponding to time=0 in file1."
)
@click.option(
    "--start-date2", required=True, help="Start date corresponding to time=0 in file2."
)
@click.option("--transition-date", required=True, help="Date where splicing occurs.")
@click.option(
    "--transition",
    default="prefer_last",
    type=click.Choice(["prefer_first", "prefer_last"]),
    help="Transition rule for ts_splice.",
)
@click.option(
    "--output-dir", required=True, help="Directory to write spliced NetCDF file."
)
@click.option("--output-name", default=None, help="Optional output filename.")
def splice_netcdf_cli(
    file1,
    file2,
    start_date1,
    start_date2,
    transition_date,
    transition,
    output_dir,
    output_name,
):
    """
    Splice two NetCDF datasets containing tracer time series.

    The output starts at the beginning of file1 and ends at the
    end of file2. The transition_date marks the splice boundary.
    """

    ds_out = splice_two_files(
        file1, file2, start_date1, start_date2, transition_date, transition
    )

    import os

    os.makedirs(output_dir, exist_ok=True)

    if output_name:
        outfile = os.path.join(output_dir, output_name)
    else:
        outfile = os.path.join(output_dir, "spliced_output.nc")

    ds_out.to_netcdf(outfile)

    print("Wrote:", outfile)


if __name__ == "__main__":
    splice_netcdf_cli()
