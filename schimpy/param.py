#!/usr/bin/env python
# -*- coding: utf-8 -*-

import schimpy.nml as nml
import pandas as pd



class Params(object):

    def __init__(self,fname,default=None):
        self._namelist = nml.parse(fname)
        self.cached_other={}
        self.default = self.process_default(default)

    def process_default(self,default):
        """Process default parameters

        Parameters
        ----------
        default : str or Param instance or 'repo'

            If default is 'repo', default is the GitHub repo version. If it is a filename, that name is evaluated. 
            If it is a templated version, that version is read. (config part of this not done yet)

        """
        
        return None

    def adjust_dt():
        """Adjust dt without perturbing variables that depend on it"""
        raise NotImplementedError("Adjusting dt safely not implemented")


    def set_run_start(self,run_start):
        """ Set start time
        
        Parameters
        ----------
            start : datetime
        Coercible to datetime

        """
        stime = run_start if type(run_start) == pd.Timestamp else pd.to_datetime(run_start)
        # parse them all before we start setting any
        syear= stime.year
        smonth= stime.month
        sday = stime.day
        shour= stime.hour

        self['start_year']  = syear
        self['start_month'] = smonth
        self['start_day']  = sday
        self['start_hour'] = shour

    def get_run_start(self):
        """Get start time as datetime"""
        syear= self['start_year']
        smonth= self['start_month']
        sday = self['start_day']
        shour= self['start_hour']
        return pd.Timestamp(year=syear,month=smonth,day=sday,hour=shour)

    run_start = property(get_run_start,set_run_start)


    def set_interval(self,name,freq):
        """Set binary output frequency using Pandas offset or string that evaluates as offset"""
        dt = int(self['dt'])
        if type(freq) in (str,pd.Timedelta): 
            freq=pd.tseries.frequencies.to_offset(freq)
            dt=pd.tseries.frequencies.to_offset(f"{dt}S")
            nspool = freq/dt           
        elif type(freq)== pd.tseries.offsets.DateOffset:
            dt=pd.tseries.frequencies.to_offset(f"{dt}S")
            nspool = freq/dt
            if abs(nspool - round(nspool)) > 0.01: 
                raise ValueError("Output interval not divisible by dt")
            else:
                nspool = round(nspool)
                print(nspool)
        else: 
            raise ValueError("Entry must be string or offset or something that is convertible to offset")
        self[name] = nspool

    def get_interval(self,name):
        dt = self['dt']
        sec=self[name]*dt
        freq = pd.Timedelta(sec,unit='S')

        return pd.tseries.frequencies.to_offset(freq)


    def set_hotstart_freq(self,freq):
        """Set hotstart frequency using Pandas offset or string that evaluates as offset
        
        Parameters
        ----------
        freq 

        Pandas offset or string that evaluates as offset. If None, frequency will be set using default (or 1 Hour) and station output disabled
        """
        self.setspool('nhot_write',freq)

    def get_hotstart_freq(self):
        return self.get_interval('nhot_write')

    hotstart_freq = property(get_hotstart_freq, set_hotstart_freq)



    def set_nc_out_freq(self,freq):
        """Set binary output frequency using Pandas offset or string that evaluates as offset"""
        self.set_interval('ihfskip',freq)


    def get_nc_out_freq(self):
        return self.get_interval("ihfskip")

    nc_out_freq = property(get_nc_out_freq, set_nc_out_freq)


    def set_nc_stack(self,freq):
        """Set binary output frequency using Pandas offset or string that evaluates as offset"""
        self.set_interval('ihfskip',freq)

    def get_nc_stack(self):
        return self.get_interval('ihfskip')

    nc_stack = property(get_nc_stack, set_nc_stack)



    def set_station_out_freq(self,freq):
        """Set station output frequency 
        
        Parameters
        ----------
        freq 

        Pandas offset or string that evaluates as offset. If None, frequency will be set using default (or 1 Hour) and station output disabled

        """
        self.set_interval('nspool_sta',freq)

    def get_station_out_freq(self):
        return self.get_interval('nspool_sta')

    station_out_freq = property(get_station_out_freq, set_station_out_freq)
    

    def diff(self,other,defaults=True):
        """Compare to another instance of Params 
        
        Parameters
        ----------
        other: str | Params
        
            Can be a string that evaluates to param filename, a Param object

        defaults : bool 
            Search includes defaults

        Returns
        -------
        diff : pd.DataFrame
            Data with multi index for section and parameter, including only parameters that 
            are different, including possibly one being absent from one Param set


        
        """
        raise NotImplementedError("Not sure what is most useful yet")

    def copy(self):
        """Return a copy of this Params set"""
        raise NotImplementedError


    def update(self,other,defaults=False):
        """Update from another instance of Params in-place
        
        Parameters
        ----------
        other: str | Params
        
            Can be a string that evaluates to param filename, a Param object

        defaults : bool 
            Search includes defaults


        
        """
        raise NotImplementedError("Not sure what is most useful yet")


    def searchfor(self,key,section=False):
        """ Search for key in all the sections 
        
        Parameters
        ----------
            key : str
            Key to search for

            section = bool
            If boolean, returns a tuple of section, value

            default = str
            Name of default to use for backup

        Returns
            Value cached under key

        Raises
            IndexError if key not present

        """
        for k in self._namelist.keys():
            sect = self._namelist[k]
            if key in sect.keys():
                if section:
                    return (k,sect[key])
                else:
                    return sect[key]

        # if we got here, the Param is not in the current file
        # search in (possibly cached) Params object named by default
        if self.default is not None:
            return self.default.searchfor(key,section=section)

        # If we get here, the key is not present in either this ParamSet or default
        raise IndexError(f"Key {key} not found in namespace")

    def __getitem__(self, key):
        item = self.searchfor(key)
        return item['value']

    def __setitem__(self, key, val):
        section, item = self.searchfor(key,section=True)
        self._namelist[section][key]['value']=val



def param_from_template(name):
    """Returns param based on named template files"""

def read_params(fname,default=None):
    with open(fname,"r") as fin:
        content=fin.read()
    p = Params(content,default)
    return p

def  test_param():
    test_param_file="C:/Delta/BayDeltaSCHISM/templates/bay_delta/param.nml.clinic"
    parms = read_params(test_param_file)
    print(parms)
    print(parms["rnday"])
    parms['rnday']=3000
    print(parms["rnday"])
    print(parms.run_start)
    parms.run_start='2010-02-04'
    print(parms.run_start)
    print("Stack")
    print(parms.nc_stack)
    parms.nc_stack = '8H'
    print(parms.nc_stack)
    print("IHFSKIP")
    print(parms["ihfskip"])
    print(parms.nc_out_freq)
    print(parms.station_out_freq)
    print(parms.hotstart_freq)


if __name__== '__main__':
    test_param()



