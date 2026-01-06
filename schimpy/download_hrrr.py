# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 09:07:50 2023
Download NOAA High-Resolution Rapid Refresh (HRRR) Model using AWS bucket service
"""
import click
from schimpy.hrr3 import *
import datetime


def download_hrr(start_date, rnday, pscr, bbox):

    hr3 = HRRR(start_date=start_date, rnday=rnday, pscr=pscr, bbox=bbox)


def download_hrr_clip(
    start_date,
    destination,
    rnday,
    latitude_min,
    latitude_max,
    longitude_min,
    longitude_max,
):
    """Download NOAA High-Resolution Rapid Refresh (HRRR) Model using AWS bucket service.

    START_DATE: Starting date of HRRR data (format: MM/DD/YYYY, e.g., 09/19/2018)

    DESTINATION: Path to store downloaded HRRR data

    RNDAY: Number of days of data to be downloaded

    LATITUDE_MIN: Minimal latitude of bounding box

    LATITUDE_MAX: Maximal latitude of bounding box

    LONGITUDE_MIN: Minimal longitude of bounding box

    LONGITUDE_MAX: Maximal longitude of bounding box

    Example:

        download_hrrr 01/01/2023 g:\\temp 15 37.3 39.0 -123.15 -121.1
    """
    bbox = [
        longitude_min,
        latitude_min,
        longitude_max,
        latitude_max,
    ]
    pscr = destination
    start_date_parsed = datetime.datetime.strptime(start_date, "%m/%d/%Y")
    download_hrr(start_date_parsed, rnday, pscr, bbox)


@click.command()
@click.argument(
    "start_date",
    type=str,
)
@click.argument(
    "destination",
    type=str,
)
@click.argument(
    "rnday",
    type=int,
)
@click.argument(
    "latitude_min",
    type=float,
)
@click.argument(
    "latitude_max",
    type=float,
)
@click.argument(
    "longitude_min",
    type=float,
)
@click.argument(
    "longitude_max",
    type=float,
)
def download_hrr_clip_cli(
    start_date,
    destination,
    rnday,
    latitude_min,
    latitude_max,
    longitude_min,
    longitude_max,
):
    """Command Line Interface for download_hrr_clip function."""
    download_hrr_clip(
        start_date,
        destination,
        rnday,
        latitude_min,
        latitude_max,
        longitude_min,
        longitude_max,
    )


if __name__ == "__main__":
    download_hrr_clip_cli()


# bbox=[-123.15,37.3,-121.1,39.0]
# pscr="G:\\temp"
# rnday=5
# import datetime as dtm
# t0=dtm.datetime(2022,3,1)
# hr3=HRRR(start_date=t0,rnday=rnday,pscr=pscr,bbox=bbox)
