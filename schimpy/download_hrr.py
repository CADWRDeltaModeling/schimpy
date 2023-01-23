# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 09:07:50 2023
Download NOAA High-Resolution Rapid Refresh (HRRR) Model using AWS bucket service
"""
import argparse
from schimpy.hrr3 import *
import datetime

def create_arg_parser():
    """ Create an argument parser
        return: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    # Read in the input file
    parser = argparse.ArgumentParser(
        description="""Download Download NOAA High-Resolution Rapid Refresh (HRRR) Model using AWS bucket service
                       
                     Example:
                           download_hrr.py 01/01/2023  g:\temp 15 -123.15 37.3 -121.1 39.0
                     """)
    parser.add_argument(dest='start_date', default=None,
                        help='starting date of HRRR data, must be format of %m/%d/%y, like 09/19/18')
    parser.add_argument(dest='destination', default=None,
                        help='path to store downloaded HRRR data')
    parser.add_argument(dest='rnday', default=None,type=int,
                        help="number of days of data to be downloaded")
    parser.add_argument(dest='latitude_min', default=None,
                        help='Minimal latitude of bounding box for raw data to be downloaded')
    parser.add_argument(dest='latitude_max', default=None,
                        help='Maximal latitude of bounding box for raw data to be downloaded')
    parser.add_argument(dest='longititude_min', default=None,
                        help='Minimal longititude of bounding box for raw data to be downloaded')
    parser.add_argument(dest='longititude_max', default=None,
                        help='Maximal longititude of bounding box for raw data to be downloaded')
    
    
    return parser

def download_hrr(start_date,rnday,pscr,bbox):
    
    hr3=HRRR(start_date=start_date,rnday=rnday,pscr=pscr,bbox=bbox)


    
def main():
    """ Main function
    """
    parser = create_arg_parser()
    args = parser.parse_args()
    bbox=[args.longititude_min,args.latitude_min,args.longititude_max,args.latitude_max]
    pscr=args.destination
    rnday=args.rnday
    start_date= datetime.datetime.strptime(args.start_date, '%m/%d/%Y')
    download_hrr(start_date,rnday,pscr,bbox)

if __name__ == "__main__":
    main()

    
#bbox=[-123.15,37.3,-121.1,39.0]
#pscr="G:\\temp"
#rnday=5
#import datetime as dtm
#t0=dtm.datetime(2022,3,1)
#hr3=HRRR(start_date=t0,rnday=rnday,pscr=pscr,bbox=bbox)