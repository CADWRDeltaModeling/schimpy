#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import re
import datetime as dtm


def hotstart_inventory(start,dt=90):
    # convert start time string input to datetime
    if isinstance(start,str):
        start = dtm.date.fromisoformat(start)
        #start=list(map(int, re.split('[^\d]', start)))

    hots = glob.glob("hotstart_000000_*.nc")
    print(hots)
    print(start)
    if len(hots) == 0:
        hots = glob.glob("hotstart_0000_*.nc")
    hots.sort()
    iters = [int(x.split("_")[2].replace(".nc","")) for x in hots]

    iters.sort()
    times = [start + dtm.timedelta(seconds=x*dt) for x in iters]

    for it,t in zip(iters,times):
        print(f"{it}: {t}")
	


    
def create_arg_parser():
    parser = argparse.ArgumentParser("Lookup station metadata by partial string match on id or name")
    parser.add_argument('--dt',default=90,type=int,help="Time step in seconds of model")
    parser.add_argument('--start',default="",help = 'Start time in iso-like format, e.g. 2013-12-03')

    return parser    



def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    hotstart_inventory(args.start,args.dt)


if __name__ == "__main__":
    main()

