from schimpy.model_time import read_th, get_headers
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

    vsrc = read_th(vsource_file, time_basis=time_basis, elapsed_unit=elapsed_unit)
    vsink = read_th(vsink_file, time_basis=time_basis, elapsed_unit=elapsed_unit)

    # subset by the search_term (eg: a region of source/sink terms)
    if search_term is not None:
        src_head = get_headers(vsource_head)
        vsrc = vsrc.loc[:, [search_term in x for x in src_head]]
        sink_head = get_headers(vsink_head)
        vsink = vsink.loc[:, [search_term in x for x in sink_head]]

    # clip to time constraint
    if start_date is None:
        start_date = vsrc.index[0]
    if end_date is None:
        end_date = vsrc.index[-1]
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

    flux_df = read_th(flux, time_basis=time_basis, elapsed_unit=elapsed_unit)

    if flux_head is not None:
        # Assign headers
        headers = get_headers(flux_head)
        if len(headers) != len(flux_df.columns):
            raise ValueError(
                f"{os.path.basename(flux)} has {len(flux_df.columns)} columns and the header calls for {len(headers)}!"
            )
        flux_df.columns = headers

    # Clip to time constraint
    if start_date is None:
        start_date = flux_df.index[0]
    if end_date is None:
        end_date = flux_df.index[-1]
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
        print(f"Combining {col_name} from {','.join(flux_cols)}")

        combtemp = flux_df.loc[:, flux_cols].copy()
        out_df[col_name] = combtemp.sum(axis=1)

    return out_df
