#!/usr/bin/env python
# -*- coding: utf-8 -*-

import schimpy.nml as nml
import pandas as pd
from pathlib import Path
import os
import click

from dataclasses import dataclass
from typing import Any, Dict, Optional,Union

@dataclass(frozen=True)
class ParamAlias:
    """
    Metadata for a user-facing alias.

    Attributes
    ----------
    target : str
        Name of the underlying SCHISM parameter or Params property.
    kind : {"raw", "property"}
        "raw"      -> use Params.__getitem__/__setitem__ with SCHISM name.
        "property" -> use getattr/setattr on Params (e.g. run_nday, nc_out_freq).
    description : str
        Human-readable explanation for help/--help text.
    value_hint : str
        Short hint about the value format (e.g. "integer days", "pandas offset").
    """
    target: str
    kind: str
    description: str
    value_hint: str = ""


#: Registry of aliases that can be used from code or CLI (set_param)
PARAM_ALIASES: Dict[str, ParamAlias] = {
    # --- Run length -----------------------------------------------------------
    "run_nday": ParamAlias(
        target="run_nday",
        kind="property",
        description="Total simulation length in days (maps to CORE::rnday).",
        value_hint="integer days, e.g. 252",
    ),

    # --- Hotstart mode --------------------------------------------------------
    "hotstart_mode": ParamAlias(
        target="ihot",
        kind="raw",
        description=(
            "Hotstart option: 0 = cold start; 1 = restart from hotstart.nc but "
            "reset model time to 0; 2 = restart and continue from the time "
            "stored in hotstart.nc."
        ),
        value_hint="0, 1, or 2",
    ),

    # --- Global NetCDF output freq ---------------------------------------

    # --- Time span per NetCDF file -------------------------------------------
    "nc_out_file_span": ParamAlias(
        target="nc_stack",
        kind="property",
        description=(
            "Wall-clock time covered by each history file. Must be "
            "a multiple of nc_out_freq. Used to derive ihfskip, "
            "the number of (dt) time steps in the file, which is the parameter "
            " actually used by SCHISM. nc_out_steps_per_file is an alias " 
            " that more directly expresses the ihfskip parameter."
        ),
        value_hint="'1d', '10d', etc.",
    ),

    "nc_out_steps_per_file": ParamAlias(
        target="ihfskip",
        kind="raw",
        description=(
            "Number of model time steps per history file stack (ihfskip). "
            "Derived so that ihfskip is a multiple of nspool."
        ),
        value_hint="integer time steps, e.g. 960",
    ),

    # --- Hotstart writing frequency ------------------------------------------
    "hotstart_freq": ParamAlias(
        target="hotstart_freq",
        kind="property",
        description=(
            "Wall-clock interval between writing hotstart.nc files. Implemented "
            "via nhot and nhot_write; nhot_write must be a multiple of the nc_out_file_span "
            "or, more directly in schism-speak, nhot must be a multiple of ihfskip "            
        ),
        value_hint="'5d', '12h', or 'none'",
    ),

    # --- Station outputs: wall-clock interval --------------------------------
    "station_freq": ParamAlias(
        target="station_freq",
        kind="property",
        description=(
            "Wall-clock interval between samples written to station output files "
            "(staout_*). Use 'none' to disable station output."
        ),
        value_hint="'15min', '1h', etc. or 'none'",
    ),

    # --- Station outputs: step interval (low-level) --------------------------
    "station_freq_in_model_steps": ParamAlias(
        target="nspool_sta",
        kind="raw",
        description=(
            "Number of model time steps between station outputs (nspool_sta). "
            "Only used when station output is enabled (iout_sta != 0)."
        ),
        value_hint="integer time steps",
    ),
}


def _coerce_scalar_from_string(raw: str) -> Any:
    s = str(raw).strip()
    lower = s.lower()

    if lower in ("true", "t", "yes", "y"):
        return 1
    if lower in ("false", "f", "no", "n"):
        return 0

    try:
        return int(s)
    except ValueError:
        pass

    try:
        return float(s)
    except ValueError:
        pass

    return s

def _parse_alias_value(alias: ParamAlias, raw: str) -> Any:
    """
    Interpret CLI/programmatic values for aliases.

    For interval-like properties we normalize to lowercase freq strings
    (e.g. '1h', '30min', '3d') to avoid uppercase pandas warnings.
    """
    s = str(raw).strip()

    # Allow disabling for some interval-like aliases
    if alias.target in ("hotstart_freq", "station_out_freq"):
        if s.lower() in ("none", "off", "0", "disable", "disabled"):
            return None
        return _normalize_freq_string(s)

    # Core global output controls (intervals and spans)
    if alias.target in ("nc_out_freq", "nc_stack"):
        return _normalize_freq_string(s)

    if alias.target in ("run_nday",):
        return int(s)

    # Low-level step counts
    if alias.target in ("nspool", "ihfskip", "nspool_sta", "nhot_write"):
        return _coerce_scalar_from_string(s)

    # Fallback
    return _coerce_scalar_from_string(s)


def _normalize_freq_string(freq: Union[str , pd.tseries.offsets.BaseOffset , None]) -> Optional[str]:
    """
    Normalize a pandas-style frequency:

    - Accepts None
    - Accepts offsets or strings
    - Returns a lower-case string (e.g. '1h', '30min', '3d') to avoid pandas
      deprecation warnings about uppercase codes like 'H'.
    """
    if freq is None:
        return None
    if isinstance(freq, pd.tseries.offsets.BaseOffset):
        # The .rule_code is often something like '1H' or '15T'
        s = freq.rule_code
    else:
        s = str(freq)
    return s.strip().lower()


def _timedelta_from_offset_str(freq: str) -> pd.Timedelta:
    """
    Convert a (normalized) pandas offset string like '30min', '1h', '3d'
    to a Timedelta.
    """
    return pd.to_timedelta(freq)


##reserved SCHISM parameter names for tracer modules
## THOSE NAMES ARE ONLY FOR SCHISM TRACER MODULES THAT NEEDS
## TO BE SPECIFIED IN PARAM.NML
RESERVED_TRACER_PARAM_NAMES = {
    "ntracer_gen":"GEN",
    "ntracer_age":"AGE",
    "sed_class":"SED",
    "eco_class":"ECO",
    "ntrs_icm":"ICM"
}
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

    def get_tracers(self):
        """Get tracers types and number of tracers defined in param.nml

        Returns
        -------
        tracers : a dict
            dictionary whose keys are tracer names defined in param.nml, 
            and values are number of tracers for that tracer type.

        """
  
        tracers = {}
        for key in self._namelist["CORE"].keys():
            if key in RESERVED_TRACER_PARAM_NAMES.keys():
                tracer_name = RESERVED_TRACER_PARAM_NAMES[key]
                ntracer = self._namelist["CORE"][key]["value"]
                tracers[tracer_name] = ntracer
        return tracers
    # ------------------------------------------------------------------
    # Programmatic helper: set by SCHISM name or alias
    # ------------------------------------------------------------------



    def set_by_name_or_alias(self, name: str, value: Any) -> None:
        """
        Set a parameter value using either:

        - The original SCHISM name (e.g. 'ihot', 'rnday', 'ihfskip'), or
        - A friendly alias from PARAM_ALIASES (e.g. 'run_nday', 'nc_out_freq').

        Examples
        --------
        >>> p.set_by_name_or_alias("ihot", 0)
        >>> p.set_by_name_or_alias("run_nday", 365)
        >>> p.set_by_name_or_alias("nc_out_freq", "1h")
        """
        alias = PARAM_ALIASES.get(name)

        if alias is None:
            # No alias: treat 'name' as a raw SCHISM key
            if isinstance(value, str):
                value = _coerce_scalar_from_string(value)
            self[name] = value
            return

        # Alias: either go through underlying SCHISM key or property
        if alias.kind == "raw":
            target_key = alias.target
            if isinstance(value, str):
                value = _parse_alias_value(alias, value)
            self[target_key] = value
            return

        if alias.kind == "property":
            # property on Params, e.g. run_nday, nc_out_freq, station_out_freq
            parsed = _parse_alias_value(alias, value)
            setattr(self, alias.target, parsed)
            return

        raise ValueError(f"Unknown alias kind {alias.kind!r} for {name!r}")


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


    def get_freq(self, name):
        dt = self["dt"]
        sec = self[name] * dt
        freq = pd.Timedelta(sec, unit="s")

        return pd.tseries.frequencies.to_offset(freq)

    # 1) Hotstart
    def get_hotstart_freq(self):
        return self.get_freq("nhot_write")

    def set_hotstart_freq(self, freq):
        if freq is None:
            self["nhot"] = 0
        else:
            self["nhot"] = 1
            self.set_freq("nhot_write", freq)

    hotstart_freq = property(get_hotstart_freq, set_hotstart_freq)

    # 2) NetCDF history output *freq* (between writes)
    def get_nc_out_freq(self):
        # should use nspool * dt
        return self.get_freq("nspool")

    def set_nc_out_freq(self, freq):
        self.set_freq("nspool", freq)

    nc_out_freq = property(get_nc_out_freq, set_nc_out_freq)

    # 3) NetCDF history file *span* (time per file)
    def get_nc_out_file_span(self):
        # ihfskip * dt
        return self.get_freq("ihfskip")

    def set_nc_out_file_span(self, freq):
        self.set_freq("ihfskip", freq)

    nc_out_file_span = property(get_nc_out_file_span, set_nc_out_file_span)

    # 4) Station output freq
    def get_station_freq(self):
        return self.get_freq("nspool_sta")

    def set_station_freq(self, freq):
        if freq is None:
            self["iout_sta"] = 0
        else:
            self.set_freq("nspool_sta", freq)
            self["iout_sta"] = 1

    station_freq = property(get_station_freq, set_station_freq)


    def get_mode(self):
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


    # getter/setter for total run length (days), stored in CORE::rnday
    def get_run_nday(self):
        """
        Return total run length (in days) as an int, based on CORE::rnday.
        """
        core = self._require_section_keys("CORE", ("rnday",))
        return int(core["rnday"]["value"])

    def set_run_nday(self, nday):
        """
        Set total run length (in days), writing CORE::rnday as an int.
        """
        core = self._require_section_keys("CORE", ("rnday",))
        core["rnday"]["value"] = int(nday)

    run_nday = property(get_run_nday, set_run_nday)


    def set_freq(self, name, freq):
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
                raise ValueError("Output freq not divisible by dt")
            else:
                nspool = round(nspool)
        else:
            print(type(freq))
            raise ValueError(
                "Entry must be string or offset or something that is convertible to offset"
            )
        self[name] = int(nspool)

    def get_station_out_freq(self):
        return self.get_freq("nspool_sta")

    def set_station_out_freq(self, freq):
        """Set station output frequency

        Parameters
        ----------
        freq : offset or string

            Sets output freq for staout files and ensures that output is enabled.
            If None, frequency will be set using default (or 1 Hour) and station output disabled

        """
        if freq is None:
            self["iout_sta"] = 0
        else:
            self.set_freq("nspool_sta", freq)
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




    # 5) Validate also reads from the known section (no global scan)
    def validate(self):
        """
        Basic consistency checks for output controls:

        - nhot_write must be divisible by ihfskip when nhot == 1
        - if iout_sta != 0, nhot_write must also be divisible by nspool_sta
        """
        # CORE: ihfskip (stack steps)
        core = self._require_section_keys("CORE", ("ihfskip",))
        ihfskip = core["ihfskip"]["value"]

        # SCHOUT: nhot, nhot_write, station options
        schout = self._require_section_keys("SCHOUT", ("nhot_write",))
        nhot_write = schout["nhot_write"]["value"]
        nhot = schout.get("nhot", {}).get("value", 0)

        if nhot == 1:
            ratio = nhot_write / ihfskip
            if abs(ratio - round(ratio)) > 1e-3:
                raise ValueError("nhot_write not divisible by ihfskip when nhot=1")

        # Station relationship: mod(nhot_write, nspool_sta) must == 0 if iout_sta != 0
        iout_sta = schout.get("iout_sta", {}).get("value", 0)
        if iout_sta != 0 and "nspool_sta" in schout:
            nspool_sta = schout["nspool_sta"]["value"]
            ratio2 = nhot_write / nspool_sta
            if abs(ratio2 - round(ratio2)) > 1e-3:
                raise ValueError(
                    "nhot_write not divisible by nspool_sta when iout_sta != 0"
                )


    def __getitem__(self, key):
        item = self.searchfor(key)     # returns {"value": ...}
        return item["value"]

    def __setitem__(self, key, val):
        section, item = self.searchfor(key, section=True)
        self._namelist[section][key]["value"] = val

    def write(self, fname):
        self.validate()
        txt = nml.write(self._namelist)
        with open(fname, "w") as outfile:
            outfile.write(txt)
            
def format_param_alias_help() -> str:
    lines: list[str] = []
    lines.append("Alias names you can use for NAME:\n")

    for key in sorted(PARAM_ALIASES):
        a = PARAM_ALIASES[key]
        kind_label = "Params property" if a.kind == "property" else "raw SCHISM name"
        lines.append(f"  {key}")
        lines.append(f"    maps to : {a.target} ({kind_label})")
        if a.value_hint:
            lines.append(f"    value   : {a.value_hint}")
        lines.append(f"    meaning : {a.description}")
        lines.append("")

    return "\n".join(lines)



PARAM_ALIAS_HELP_TEXT = "\n" + format_param_alias_help()


PARAM_ALIAS_HELP_TEXT = format_param_alias_help()



def param_from_template(name):
    """Returns param based on named template files"""


def read_params(fname, default=None):
    with open(fname, "r") as fin:
        content = fin.read()
    p = Params(content, default)
    return p

def set_param(param_file: Path, pairs: tuple[str, ...], output: Path | None, dry_run: bool):
    """
    Set one or more parameters in a SCHISM param.nml file.

    PARAM_FILE is the path to the existing param.nml
    (e.g. param.nml, param.nml.clinic, retro_2013/param.nml).

    Remaining arguments must be NAME VALUE pairs:

        NAME1 VALUE1 [NAME2 VALUE2 [...]]

    NAME can be either:

      * the original SCHISM name (ihot, rnday, ihfskip, nspool_sta, ...)
      * an alias defined in PARAM_ALIASES
        (run_nday, nc_out_freq, station_out_freq, hotstart_freq, ...)

    VALUE is a scalar or freq:

      * integers / floats: 0, 252, 90.0
      * time freqs:    15min, 1H, 3D (pandas-style offsets)
      * 'none' for aliases like hotstart_freq / station_out_freq to disable them
      """
    changes = list(zip(pairs[::2], pairs[1::2]))

    # --- read once ---
    try:
        params = read_params(str(param_file))
    except Exception as e:
        raise click.ClickException(f"Failed to read {param_file}: {e}") from e

    click.echo(f"PARAM_FILE: {param_file}")

    # apply all changes in-memory
    for name, value in changes:
        # best-effort "before"
        before = None
        try:
            if name in PARAM_ALIASES:
                alias = PARAM_ALIASES[name]
                if alias.kind == "raw":
                    before = params[alias.target]
                else:
                    before = getattr(params, alias.target)
            else:
                before = params[name]
        except Exception:
            pass

        # apply via alias-aware helper
        try:
            params.set_by_name_or_alias(name, value)
        except Exception as e:
            raise click.ClickException(
                f"Failed to set {name} to {value!r}: {e}"
            ) from e

        # best-effort "after"
        after = None
        try:
            if name in PARAM_ALIASES:
                alias = PARAM_ALIASES[name]
                if alias.kind == "raw":
                    after = params[alias.target]
                else:
                    after = getattr(params, alias.target)
            else:
                after = params[name]
        except Exception:
            pass

        if before is not None or after is not None:
            click.echo(f"  {name}: {before!r} -> {after!r}")
        else:
            click.echo(f"  {name}: set to {value!r}")

    if dry_run:
        click.echo("Dry-run: no file written.")
        return

    dest = output or param_file

    # --- write once ---
    try:
        params.write(str(dest))
    except Exception as e:
        raise click.ClickException(f"Failed to write {dest}: {e}") from e

    click.echo(f"Wrote updated parameters to {dest}")

@click.command(
    name="set_param",
    context_settings={"help_option_names": ["-h", "--help"]},
    epilog=PARAM_ALIAS_HELP_TEXT,
)
@click.argument(
    "param_file",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
)
@click.argument("pairs", nargs=-1)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help=(
        "Write to this file instead of editing PARAM_FILE in place. "
        "If omitted, PARAM_FILE is overwritten."
    ),
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Parse and report changes but do not write anything.",
)
def set_param_cli(param_file: Path, pairs: tuple[str, ...], output: Path | None, dry_run: bool):
    """
    Command line tool to set SCHISM model parameters in PARAM_FILE.
    """
    if not pairs:
        raise click.ClickException(
            "No NAME VALUE pairs supplied.\n"
            "Usage: set_param PARAM_FILE NAME VALUE [NAME VALUE ...]"
        )

    if len(pairs) % 2 != 0:
        raise click.ClickException(
            "Expect an even number of arguments after PARAM_FILE: NAME VALUE pairs.\n"
            "Example: set_param param.nml ihot 0 run_nday 365"
        )

    set_param(param_file, pairs, output, dry_run)



if __name__ == "__main__":
    set_param_cli()
