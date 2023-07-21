#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import re
import datetime as dtm
import os
import pandas as pd
from schimpy.param import *
from vtools import days



def hotstart_inventory(run_start=None,
                       dt=None,
                       nday=None,
                       workdir='.',
                       paramfile=None,
                       hot_freq=None
                       ):
    """ Create an inventory of existing hotstarts or expected hotstarts

    Existing vs expected depends on whether run in an outputs or study dir, which is detected.
    
    Parameters
    ----------

    start : Convertible to Datetime
       Start date of run. If None inferred from paramfile. Error if 
    
    dt : float
       dt in seconds, for instance perhaps 90s for clinic or 120s for tropic. This
       is needed to intepret the so-called iteration number in the hotstart labels 
       which are really time step numbers.

    workdir : str
        Directory to inventory. If it is an outputs directory, the inventory 
    comprise existing hotstarts. 
    
    paramfile : str
       Name of param.nml file, expected in workdir or workdir/.. If None, then
       both start and dt must be supplied. If all three are None, the name "param.nml" will
       be attempted

    
    Returns
    -------
        Dataframe, if programatic, listing hotstarts.  Should look like "Date, Iteration" pairs. CLI should print it out.
    
    Notes
    -----
        If the listing is done in the run dir, this will be the expected hotstarts. If 
        it is in the outputs dir, it will be an inventory of existing hotstarts.

    """
    run_start = pd.to_datetime(run_start)
    if type(hot_freq) == str:
        hot_freq = pd.tseries.frequencies.to_offset(hot_freq)
    print(run_start,dt,paramfile,workdir,hot_freq)
    hots = glob.glob(os.path.join(workdir,"hotstart_000000_*.nc"))
    is_existing = (len(hots) > 0)
    print("is_existing",is_existing)
    param_needed = (dt is None) or (run_start is None) or (hot_freq is None and not is_existing ) 

    if param_needed:
        if paramfile is None: 
            paramfile = 'param.nml'
        if is_existing or not os.path.exists(paramfile):
            paramfile = os.path.join(workdir,'..',paramfile)
        
        params = read_params(paramfile)
        dt_param = params['dt']
        run_start_param = params.run_start
        run_len_param = params['rnday']
        hot_freq_param = params.hotstart_freq
        print("run_start",run_start,run_start_param)
        if run_start is None or pd.isnull(run_start) : 
            print("yo",run_start_param)
            run_start = run_start_param
        if dt is None: dt = dt_param
        if nday is None or nday == 0: nday = run_len_param
        if hot_freq is None: hot_freq=hot_freq_param
    if not is_existing:
        print("hello",run_start,dt,paramfile,workdir,hot_freq)
        end = run_start + days(nday)
        t = run_start
        iters = []
        times = []
        dt_freq = pd.tseries.offsets.Second(dt)
        iters_per_hot = hot_freq/dt_freq
        itr = 0
        print(hot_freq,dt_freq,t,end,nday)
        while (t<end):
            times.append(t) 
            iters.append(itr)
            t = t + hot_freq
            itr = itr + int(iters_per_hot)
            print(f"t={t},itr={itr}")
        ndx = pd.DatetimeIndex(times)
        df = pd.DataFrame(index=ndx,data=iters)
        df.columns=["iteration"]
        df.index.name="datetime"
        print(df)


def hotstart_inventory2(start,dt=90):
    # convert start time string input to datetime
    if isinstance(start,str):
        start = dtm.date.fromisoformat(start)
        #start=list(map(int, re.split('[^\d]', start)))

    hots = glob.glob("hotstart_000000_*.nc")
    if len(hots) == 0:
        hots = glob.glob("hotstart_0000_*.nc")
    hots.sort()
    iters = [int(x.split("_")[2].replace(".nc","")) for x in hots]

    iters.sort()
    times = [start + dtm.timedelta(seconds=x*dt) for x in iters]

    for it,t in zip(iters,times):
        print("{}: {}".format(it,t))
   
def create_arg_parser():
    parser = argparse.ArgumentParser("Lookup station metadata by partial string match on id or name")
    parser.add_argument('--dt',default=90,type=int,help="Time step in seconds of model")
    parser.add_argument('--run_start',default="",help = 'Start time in iso-like format, e.g. 2013-12-03')
    parser.add_argument('--nday',default=0,type=int,help = 'Start time in iso-like format, e.g. 2013-12-03')
    parser.add_argument('--workdir',default='.',type=str,help="Time step in seconds of model")
    parser.add_argument('--paramfile',default="",type=str,help = 'Start time in iso-like format, e.g. 2013-12-03')
    parser.add_argument('--hot_freq',default=None,help = 'Start time in iso-like format, e.g. 2013-12-03')
    return parser    



def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    hotstart_inventory(args.run_start,args.dt,args.nday,args.workdir,args.paramfile,args.hot_freq)


if __name__ == "__main__":
    main()
