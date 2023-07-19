#!/usr/bin/env python
# -*- coding: utf-8 -*-

import namelist



class Params)(object):

    def __init__(self,fname,default=None):
        self.__namelist__ = naemelist.parse(fname)
        self.cached_other={}
        self.default = process_default



    def set_start_time(self,start):
        """ Set start time
        
        Parameters
        ----------
            start : datetime
        Coercible to datetime

        """
        pass

    def get_start_time(self):
        """Get start time as datetime"""
        pass

    start_time = property(get_start_time,get_start_time)

    def set_nday(self,datetime):
        pass

    def get_nday(self)
        pass

    nday = property(get_nday,get_nday)


    def hotstart_freq(self,freq):
        """Set hotstart frequency using Pandas offset or string that evaluates as offset
        
        Parameters
        ----------
        freq 

        Pandas offset or string that evaluates as offset. If None, frequency will be set using default (or 1 Hour) and station output disabled
        """
        pass

    def hotstart_freq(self)
        pass

    hotstart_interval = property(get_hotstart_interval, set_hotstart_interval)


    def set_nc_freq(self,freq):
        """Set binary output frequency using Pandas offset or string that evaluates as offset"""
        pass

    def get_nc_freq(self)
        pass

    nc_out_freq = property(get_nc_freq, set_nc_freq)


    def set_station_freq(self,freq):
        """Set station output frequency 
        
        Parameters
        ----------
        freq 

        Pandas offset or string that evaluates as offset. If None, frequency will be set using default (or 1 Hour) and station output disabled

        """
        pass

    def get_station_freq(self)
        pass

    nc_out_freq = property(get_nc_freq, set_nc_freq)
    

    def diff(self,other):
        #produces a list of differences
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
        for k in self.__namelist__.keys:
            if key in k.keys
                if complete:
                    return (k,key)
                else:
                    return k[key]

        # if we got here, the Param is not in the current file
        # search in (possibly cached) Params object named by default
        if self.default is not None:
            return self.default.searchfor(key,section=section)

        # If we get here, the key is not present in either this ParamSet or default
        raise IndexError(f"Key {key} not found in namespace")

    def __getitem__(self, key):
        item = searchfor(key)
        return item

    def __setitem__(self, key, val):
        section, item = searchfor(key,section=True)
        self.__namelist__[section][item]=val



def param_from_template(name):
    """Returns param based on named template files"""








