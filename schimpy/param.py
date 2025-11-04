#!/usr/bin/env python
# -*- coding: utf-8 -*-

import schimpy.nml as nml
import pandas as pd
from pathlib import Path
import os

class Params(object):

    def __init__(self, fname, default=None):
        # ...
        # OLD (problem): self._namelist = nml.parse(fname)
        # NEW (robust): read text if a path is given
        from . import nml  # wherever your nml lives

        # Accept either text or a path-like
        # if isinstance(fname, (str, os.PathLike, Path)):
        try:
            p = Path(fname)
            if p.exists():
                text = p.read_text(encoding="utf-8")
            else:
                # If it's a path-like string that doesn't exist, treat as text
                text = str(fname)
        except:
            text = fname  # already text

        self._namelist = nml.parse(text)  # feed TEXT to the parser

        if not isinstance(self._namelist, dict) or not self._namelist:
            raise ValueError(
                f"Parsed namelist is empty for input {fname!r}. "
                "Ensure you passed file contents or a valid file path."
            )

        self.default = self.process_default(default)

    def _require_section_keys(self, section: str, keys: tuple[str, ...]) -> dict:
        """
        Ensure a named section exists and that all required keys are present.
        Returns the section dict; raises KeyError with specifics otherwise.
        """
        sect = self._namelist.get(section)
        if not isinstance(sect, dict):
            raise KeyError(f"Section '{section}' not found in namelist.")
        missing = [k for k in keys if k not in sect]
        if missing:
            raise KeyError(f"Missing required keys in section '{section}': {missing}")
        return sect



    def process_default(self, default):
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

    # Get Properties --------------------------------------------------------------
    # All functions which retrieve properties from param.nml here:
    # -----------------------------------------------------------------------------

    def get_run_start(self):
        import pandas as pd
        opt = self._require_section_keys("OPT", ("start_year", "start_month", "start_day", "start_hour"))
        y = int(opt["start_year"]["value"])
        m = int(opt["start_month"]["value"])
        d = int(opt["start_day"]["value"])
        h = int(opt["start_hour"]["value"])
        return pd.Timestamp(y, m, d, h)


    def get_interval(self, name):
        dt = self["dt"]
        sec = self[name] * dt
        freq = pd.Timedelta(sec, unit="s")

        return pd.tseries.frequencies.to_offset(freq)

    def get_hotstart_freq(self):
        return self.get_interval("nhot_write")

    def get_nc_out_freq(self):
        return self.get_interval("ihfskip")

    def get_nc_stack(self):
        return self.get_interval("ihfskip")

    def get_station_out_freq(self):
        return self.get_interval("nspool_sta")

    def get_baro(self):
        ibc = self._namelist["CORE"]["ibc"]["value"]

        return "clinic" if ibc == 0 else "tropic"

    # Set Properties --------------------------------------------------------------
    # All functions which set properties to param.nml here:
    # -----------------------------------------------------------------------------

    def set_run_start(self, run_start):
        import pandas as pd
        ts = pd.to_datetime(run_start)
        opt = self._require_section_keys("OPT", ("start_year", "start_month", "start_day", "start_hour"))
        opt["start_year"]["value"]  = int(ts.year)
        opt["start_month"]["value"] = int(ts.month)
        opt["start_day"]["value"]   = int(ts.day)
        opt["start_hour"]["value"]  = int(ts.hour)

    run_start = property(get_run_start, set_run_start)

    def set_interval(self, name, freq):
        """Set binary output frequency using Pandas offset or string that evaluates as offset"""
        dt = int(self["dt"])
        if type(freq) in (str, pd.Timedelta):
            freq = pd.tseries.frequencies.to_offset(freq)
            dt = pd.tseries.frequencies.to_offset(f"{dt}s")
            nspool = freq / dt
        elif isinstance(freq, pd.offsets.DateOffset):
            dt = pd.tseries.frequencies.to_offset(f"{dt}s")
            nspool = freq / dt
            if abs(nspool - round(nspool)) > 0.01:
                raise ValueError("Output interval not divisible by dt")
            else:
                nspool = round(nspool)
        else:
            print(type(freq))
            raise ValueError(
                "Entry must be string or offset or something that is convertible to offset"
            )
        self[name] = int(nspool)

    def set_hotstart_freq(self, freq):
        """Set hotstart frequency using Pandas offset or string that evaluates as offset

        Parameters
        ----------
        freq

        If None, frequency will be set using default (or 1 Hour) and station output disabled
        """
        if freq is None:
            self["nhot"] = 0
            return
        else:
            self["nhot"] = 1
            self.set_interval("nhot_write", freq)

    hotstart_freq = property(get_hotstart_freq, set_hotstart_freq)

    def set_nc_out_freq(self, freq):
        """Set binary output frequency using Pandas offset or string that evaluates as offset"""
        self.set_interval("ihfskip", freq)

    nc_out_freq = property(get_nc_out_freq, set_nc_out_freq)

    def set_nc_stack(self, freq):
        """Set binary output frequency using Pandas offset or string that evaluates as offset"""
        self.set_interval("ihfskip", freq)

    nc_stack = property(get_nc_stack, set_nc_stack)

    def set_station_out_freq(self, freq):
        """Set station output frequency

        Parameters
        ----------
        freq : offset or string

            Sets output interval for staout files and ensures that output is enabled.
            If None, frequency will be set using default (or 1 Hour) and station output disabled

        """
        if freq is None:
            self["iout_sta"] = 0
        else:
            self.set_interval("nspool_sta", freq)
            self["iout_sta"] = 1

    station_out_freq = property(get_station_out_freq, set_station_out_freq)

    def sections(self, defaults=False):
        sections = self._namelist.keys()
        if defaults and self.default is not None:
            added = []
            for item in self.defaults.keys():
                if item not in sections:
                    added.append(item)
            sections.extend(added)
        return sections

    def to_dataframe(self, defaults=None):
        sections = self.sections(defaults)
        indices = []
        values = []

        for key in sections:
            items = self._namelist[key]
            for subkey in items.keys():
                if not isinstance(items[subkey], dict):
                    continue
                indices.append((key, subkey))
                values.append(items[subkey]["value"])
        index = pd.MultiIndex.from_tuples(indices, names=["section", "name"])
        df = pd.DataFrame(index=index, data=values)
        df.columns = ["values"]
        return df

    def diff(self, other, defaults=False):
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
        if defaults:
            raise NotImplementedError(
                "Diffing with defaults not supported at this time"
            )
        this = self.to_dataframe(None)
        this.columns = ["this"]
        otherdf = other.to_dataframe(None).convert_dtypes()
        otherdf.columns = ["other"]
        joined = this.merge(otherdf, how="outer", left_index=True, right_index=True)
        joined = joined.loc[joined.this != joined.other]
        return joined

    def copy(self):
        """Return a copy of this Params set"""
        raise NotImplementedError()

    def update(self, other, defaults=False):
        """Update from another instance of Params in-place

        Parameters
        ----------
        other: str | Params

            Can be a string that evaluates to param filename, a Param object

        defaults : bool
            Search includes defaults

        """
        raise NotImplementedError("Not sure what is most useful yet")

    def searchfor(self, key, section=False):
        """
        Search for key.
        - section == False (default): scan all dict sections; return value
        - section == True: scan all dict sections; return (section_name, value)
        - section is a str (e.g., 'OPT'/'CORE'): only search that section;
        return value when section == False, else (section_name, value)
        Raises IndexError if not found.
        """
        # section-specific search (deterministic; no global scan)
        if isinstance(section, str):
            sect = self._namelist.get(section)
            if isinstance(sect, dict) and key in sect:
                return (section, sect[key]) if section not in (False, None) else sect[key]
            raise IndexError(f"Key {key} not found in section '{section}'")

        # existing behavior: scan all sections
        for sect_name, sect in self._namelist.items():
            if isinstance(sect, dict) and key in sect:
                return (sect_name, sect[key]) if section else sect[key]

        if self.default is not None:
            return self.default.searchfor(key, section=section)

        raise IndexError(f"Key {key} not found in namespace")

    # 3) Hint map so __getitem__/__setitem__ don’t rely on global scans for known keys
    _SECTION_HINTS = {
        # OPT time quartet
        "start_year": "OPT", "start_month": "OPT", "start_day": "OPT", "start_hour": "OPT",
        # CORE writer validation fields
        "nhot_write": "CORE", "ihfskip": "CORE",
    }

    def __getitem__(self, key):
        sect_name = self._SECTION_HINTS.get(key)
        if sect_name:
            # deterministic, section-specific
            item = self.searchfor(key, section=sect_name)
            return item["value"]
        # fallback (legacy)
        item = self.searchfor(key)
        return item["value"]

    def __setitem__(self, key, val):
        sect_name = self._SECTION_HINTS.get(key)
        if sect_name:
            sect_name, _ = self.searchfor(key, section=sect_name)  # raises if absent
            self._namelist[sect_name][key]["value"] = val
            return
        sect_name, _ = self.searchfor(key, section=True)  # legacy path
        self._namelist[sect_name][key]["value"] = val

    # 5) Validate also reads from the known section (no global scan)
    def validate(self):
        core = self._require_section_keys("CORE", ("ihfskip",))
        schout = self._require_section_keys("SCHOUT", ("nhot_write",))
        nhotwrite = schout["nhot_write"]["value"]
        ihfskip   = core["ihfskip"]["value"]
        ratio = nhotwrite / ihfskip
        if abs(ratio - round(ratio)) > 1e-3:
            raise ValueError("nhot_write not divisible by ihfskip")

    def __getitem__(self, key):
        item = self.searchfor(key)
        return item["value"]

    def __setitem__(self, key, val):
        section, item = self.searchfor(key, section=True)
        self._namelist[section][key]["value"] = val

    def write(self, fname):
        self.validate()
        txt = nml.write(self._namelist)
        with open(fname, "w") as outfile:
            outfile.write(txt)


def param_from_template(name):
    """Returns param based on named template files"""


def read_params(fname, default=None):
    with open(fname, "r") as fin:
        content = fin.read()
    p = Params(content, default)
    return p


def test_param():
    test_param_file = "C:/Delta/BayDeltaSCHISM/templates/bay_delta/param.nml.clinic"
    parms = read_params(test_param_file)
    print(parms)
    print(parms["rnday"])
    parms["rnday"] = 3000
    print(parms["rnday"])
    print(parms.run_start)
    parms.run_start = "2010-02-04"
    print(parms.run_start)
    print("Stack")
    print(parms.nc_stack)
    parms.nc_stack = "8h"
    print(parms.nc_stack)
    print("IHFSKIP")
    print(parms["ihfskip"])
    print(parms.nc_out_freq)
    print(parms.station_out_freq)
    parms.station_out = None
    print(parms.station_out_freq)
    print(parms["iout_sta"])
    parms.station_out_freq = "15min"
    print(parms.station_out_freq)
    print(parms["iout_sta"])

    print(parms.hotstart_freq)
    print(parms["nhot"])
    parms.hotstart_freq = "10D"
    print(parms["nhot"])
    print(parms.hotstart_freq)
    parms.hotstart_freq = None
    print(parms["nhot"])
    parms.hotstart_freq = "5D"
    print(parms["nhot"])

    print(parms.sections())

    other_param_file = "C:/Delta/BayDeltaSCHISM/templates/bay_delta/param.nml.tropic"
    otherparms = read_params(other_param_file)
    df = parms.diff(otherparms)

    parms.nc_stack = pd.tseries.frequencies.to_offset("1D")
    parms.validate()
    parms.write("./junk.nml")


if __name__ == "__main__":
    test_param()
