from schimpy.model_time import read_th
import pandas as pd
import os


def calc_net_source_sink(
    vsource_file="vsource.th",
    vsink_file="vsink.th",
    time_basis=None,
    elapsed_unit="s",
    vsource_head=None,
    vsink_head=None,
    search_term=None,
    start_date=None,
    end_date=None,
):
    """Calculate the net source/sink terms supplied by the vsource.th and vsink.th files.
    Files can be specified, elapsed, or timestamped.
    By default it will look for vsource.th and vsink.th - the default source/sink inputs for SCHISM.
    """

    vsrc = read_th(
        vsource_file,
        time_basis=time_basis,
        elapsed_unit=elapsed_unit,
        head=vsource_head,
    )
    vsink = read_th(
        vsink_file, time_basis=time_basis, elapsed_unit=elapsed_unit, head=vsink_head
    )

    # subset by the search_term (eg: a region of source/sink terms)
    if search_term is not None:
        vsrc = vsrc.loc[:, [search_term in x for x in vsrc.columns]]
        vsink = vsink.loc[:, [search_term in x for x in vsink.columns]]

    # clip to time constraint
    if start_date is None:
        start_date = vsrc.index[0]
    if end_date is None:
        end_date = vsrc.index[-1]

    if (vsrc.index[0] > pd.to_datetime(end_date)) or (
        vsrc.index[-1] < pd.to_datetime(start_date)
    ):
        raise ValueError(
            f"File: {vsource_file} does not cover the dates requested. \n\tRequested dates are: {start_date} to {end_date}, \n\tand the file covers {vsrc.index[0]} to {vsrc.index[-1]}"
        )
    vsrc = vsrc[start_date:end_date]
    vsink = vsink[start_date:end_date]

    # Get the individual net source and net sink
    vsrc = vsrc.sum(axis=1)
    vsink = vsink.sum(axis=1)

    # Combine for the net timeseries
    net_df = vsrc + vsink

    return net_df, vsrc, vsink


def read_flux(
    flux="flux.th",
    flux_head=None,
    time_basis=None,
    elapsed_unit="s",
    start_date=None,
    end_date=None,
):
    """Read in the flux.th, use flux_head to assign headers to dataframe.
    If no flux_head, then the resultant dataframe will have no headers.
    Optionally subset by start_date and end_date.
    """

    flux_df = read_th(
        flux, time_basis=time_basis, elapsed_unit=elapsed_unit, head=flux_head
    )

    # Clip to time constraint
    if start_date is None:
        start_date = flux_df.index[0]
    if end_date is None:
        end_date = flux_df.index[-1]
    if (flux_df.index[0] > pd.to_datetime(end_date)) or (
        flux_df.index[-1] < pd.to_datetime(start_date)
    ):
        raise ValueError(
            f"File: {flux} does not cover the dates requested. \n\tRequested dates are: {start_date} to {end_date}, \n\tand the file covers {flux_df.index[0]} to {flux_df.index[-1]}"
        )
    flux_df = flux_df[start_date:end_date]

    return flux_df


def combine_flux(
    flux,
    comb_dict,
    time_basis=None,
    elapsed_unit="s",
    flux_head=None,
    start_date=None,
    end_date=None,
):
    """Combine flux data,
    comb_dict is a dictionary that combines the headers in the list into a single value with the header of the dictionary key.
    An example of this would be comb_dict={'Northern':['this','that','other']}
    Optionally subset by start_date and end_date.
    """

    if not isinstance(flux, (pd.DataFrame, pd.Series)):
        flux_df = read_flux(
            flux=flux,
            flux_head=flux_head,
            time_basis=time_basis,
            elapsed_unit=elapsed_unit,
            start_date=start_date,
            end_date=end_date,
        )
    else:
        flux_df = flux
        # Clip to time constraint if needed
        if start_date is None:
            start_date = flux_df.index[0]
        if end_date is None:
            end_date = flux_df.index[-1]
        flux_df = flux_df[start_date:end_date]

    out_df = pd.DataFrame(index=flux_df.index)

    # Read through comb_dict items and output that dataframe
    for col_name, flux_cols in comb_dict.items():
        print(f"\tCombining {col_name} from {','.join(flux_cols)}")

        combtemp = flux_df.loc[:, flux_cols].copy()
        out_df[col_name] = combtemp.sum(axis=1)

    return out_df


def normalize_df(df):
    normed = df.copy()
    for col in normed.columns:
        col_min = normed[col].min()
        col_max = normed[col].max()
        if col_max == col_min:
            if col_min != 0:
                normed[col] = 1
            else:
                normed[col] = 0
        else:
            normed[col] = (normed[col] - col_min) / (col_max - col_min)
    return normed


def struct_open_props(struct, th_data, out_freq="15min", datetime_idx=None):
    """Determine the structure header properties that are used to determine 'relative openness' of the structure
    For instance, a simple gate might just be (width * height)/(max(width) * max(height))
    """

    if struct["type"] == "weir":
        # width, elevation, op_down, op_up
        up_cols_to_multiply = ["install", "ndup", "op_up", "elev", "width"]
        down_cols_to_multiply = ["install", "ndup", "op_down", "elev", "width"]
    elif struct["type"] == "radial":
        # width, height, elevation, op_down, op_up
        up_cols_to_multiply = ["install", "ndup", "op_up", "elev", "width", "height"]
        down_cols_to_multiply = [
            "install",
            "ndup",
            "op_down",
            "elev",
            "width",
            "height",
        ]
    elif struct["type"] == "culvert":
        # rad, elevation, op_down, op_up
        up_cols_to_multiply = ["install", "ndup", "op_up", "elev", "rad"]
        down_cols_to_multiply = ["install", "ndup", "op_down", "elev", "rad"]
    elif struct["type"] == "radial_relheight":
        # width, height, elevation, op_down, op_up
        up_cols_to_multiply = ["install", "ndup", "op_up", "elev", "width", "height"]
        down_cols_to_multiply = [
            "install",
            "ndup",
            "op_down",
            "elev",
            "width",
            "height",
        ]
    elif struct["type"] == "weir_culvert":
        # elevation, width, op_down, op_up - weirs
        # elevation, width, op_up, op_down - culvert
        up_cols_to_multiply = [
            "install",
            "ndup_weir",
            "ndup_pipe",
            "elev_weir",
            "elev_pipe",
            "width_weir",
            "radius_pipe",
            "op_up_weir",
            "up_op_pipe",
        ]
        down_cols_to_multiply = [
            "install",
            "ndup_weir",
            "ndup_pipe",
            "elev_weir",
            "elev_pipe",
            "width_weir",
            "radius_pipe",
            "op_down_weir",
            "down_op_pipe",
        ]
    else:
        raise ValueError(f"structure type {struct['type']} is not yet supported!")

    # Take the dataframe and normalize to max/min of each column
    norm_df = normalize_df(th_data)
    norm_df["open_up"] = norm_df[up_cols_to_multiply].prod(axis=1)
    norm_df["open_down"] = -norm_df[down_cols_to_multiply].prod(axis=1)
    norm_df = norm_df.resample(out_freq).ffill()
    if datetime_idx is not None:
        norm_df = norm_df.reindex(datetime_idx)
        norm_df = norm_df.ffill()

    up_df = norm_df[["open_up"]].copy()
    up_df.columns = [f"{struct['name']}_up"]
    down_df = norm_df[["open_down"]].copy()
    down_df.columns = [f"{struct['name']}_down"]

    return up_df, down_df
