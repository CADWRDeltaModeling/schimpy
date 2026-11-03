
# Reading model station output
Model station output for schism is read using the APIs below. Typically the data are on modeldir/outputs or lifted into Windows. The APIs require a reference time which is used to convert indices from elapsed to datetime indexes. The data are regular and the index need not be tested for integrety or possession of an index. For multiple outputs that overlap, vtools.ts splice() can be used (below). For specific conductance:
ts = psu_ec_25c(ts_in)

For filtering:
from vtools import cosine_lanczos
ts = ts.interpolate(limit=4)
tsfilt = cosine_lanczos(ts,'40h')

For caching read the caching notebook https://github.com/CADWRDeltaModeling/dms_datastore/blob/master/notebooks/cache.ipynb

Recommend the schism conda environment, which has all dependencies.

## Available names

station_variables = [
    "elev",
    "air pressure",
    "wind_x",
    "wind_y",
    "temp",
    "salt",
    "u",
    "v",
    "w",
    "ssc",
]

## Station data for variables other than flow 

from schimpy.station import read_staout
def read_staout(
    fname,
    station_infile,
    reftime,
    ret_station_in=False,
    multi=False,
    elim_default=False,
    time_unit="s",
):
    """Read a SCHISM staout_* file into a pandas DataFrame

    Parameters
    ----------
    fpath : fname
        Path to input staout file or a variable name in ["elev", "air pressure", "wind_x", "wind_y",  "temp", "salt", "u", "v", "w"] whose
        1-index will be mapped to a name like staout_1 for elev, staout_2 for air pressure, etc.

    station_infile : str or DataFrame
        Path to station.in file or DataFrame from read_station_in

    reftime : Timestampe
        Start of simulation, time basis for staout file elapse time

    ret_station_in : bool
        Return station_in DataFrame for use, which may speed reading of a second file

    multi : bool
        Should the returned data have a multi index for the column with location and sublocation. If False the two are collapsed

    elim_default : bool
        If the MultiIndex is collapsed, stations with subloc "default" will be collapsed. Eg. ("clc","default") becomes "clc"

    time_unit : string
        Convertible to pandas frequency string, this is the timestamp of the file.

     Returns
     -------
     Result : DataFrame
         DataFrame with hierarchical index (id,subloc) and columns representing the staout data (collapsed as described above

    Examples
    --------

    >>> staout1,station_in = read_staout("staout_1","station.in",reftime=pd.Timestamp(2009,2,10),
                                 ret_station_in = True,multi=False,elim_default=True)
    >>> staout6 = read_staout("staout_6",station_in,reftime=pd.Timestamp(2009,2,10),multi=False,elim_default=True)

    """

## Flow
from schimpy.station import read_staout
	def read_flux_out(fpath, names, reftime):
    """Read fluxes from a SCHISM flux.out file

    Parameters
    ----------

    fpath : str
    Path to the file

    names : str
    name of file that contains  names of  flux areas, typically something like flow_xsects.yaml or fluxflags.prop

    reftime : str
    start of simulation, against which relative times will be calculated
    """
	
from vtools import ts_splice
@align_inputs_strict(seq_arg=0, names_kw="names")
def ts_splice(
    series, names=None, transition="prefer_last", floor_dates=False, preserve_freq=True
):
    """
    Splice multiple time series together, prioritizing series in patches of time.

    Unlike `ts_merge`, which blends overlapping data points, `ts_splice` stitches
    together time series without overlap. The function determines when to switch
    between series based on a transition strategy.

    Parameters
    ----------
    series : tuple or list of pandas.DataFrame or pandas.Series
        A tuple or list of time series.
        Each series must have a `DatetimeIndex` and consistent column structure.

    names : None, str, or iterable of str, optional
        - If `None` (default), all input series must share common column names,
          and the output will merge common columns.
        - If a `str`, all input series must have a **single column**, and the
          output will be a DataFrame with this name as the column name.
        - If an iterable of `str`, all input DataFrames must have the same
          number of columns matching the length of `names`, and these will be
          used for the output.

    transition : {'prefer_first', 'prefer_last'} or list of pandas.Timestamp
        Defines how to determine breakpoints between time series:
        - `'prefer_first'`: Uses the earlier series on the list during until its valid timestamp.
        - `'prefer_last'`: Uses the later series starting from its first valid timestamp.
        - A list of specific timestamps can also be provided as transition points.

    floor_dates : bool, optional, default=False
        If `True`, inferred transition timestamps (`prefer_first` or `prefer_last`)
        are floored to the beginning of the day. This can introduce NaNs if the
        input series are regular with a `freq` attribute.

    preserve_freq : bool, default True
        If True, attempt to preserve a common regular frequency across inputs.
        When a common regular frequency can be established, the output is
        reindexed to that regular extent. When it cannot be established,
        ts_merge falls back to index-based merge behavior, which for strict_priority=True
        is a patchwise combination of the indexes and for False is a shuffle.

    Returns
    -------
    pandas.DataFrame or pandas.Series
        - If the input contains multi-column DataFrames, the output is a DataFrame
          with the same column structure.
        - If a collection of single-column `Series` is provided, the output will be
          a Series.
        - The output retains a `freq` attribute if all inputs share the same frequency.

    Notes
    -----
    - The output time index is the union of input time indices.
    - If `transition` is `'prefer_first'`, gaps may appear in the final time series.
    - If `transition` is `'prefer_last'`, overlapping data is resolved in favor
      of later series.

    See Also
    --------
    ts_merge : Merges series by filling gaps in order of priority.
    """	
	




# Vertical Grid

## vgrid.in Format vs schimpy Internal Representation (ivcor=1, LSC2)

### On-disk format (vgrid.in, version 5.10)

```
1                                    ← ivcor (1 = LSC2 localized sigma)
nvrt                                 ← max number of vertical levels across all nodes
kbp(1) kbp(2) ... kbp(n_nodes)      ← 1-based bottom level index per node
1  sigma(1,1) sigma(1,2) ...         ← level 1 values for all nodes (one row per level)
2  sigma(2,1) sigma(2,2) ...
...
nvrt  sigma(nvrt,1) sigma(nvrt,2)...
```

- **Sigma range**: −1 (bottom) to 0 (surface), monotonically increasing with level index.
- **Level ordering**: level 1 is the deepest possible; level `nvrt` is the surface.
- **Per-node bottom**: `kbp(i)` (1-based) is the first valid level for node i. Levels below `kbp` are padded with **−9.0** (sentinel for "not applicable").
- **Right-justified**: valid sigma values occupy levels `kbp` through `nvrt`; everything to the left is −9.0.
- **Output precision**: `%14.6f` — only 6 decimal places. Consecutive sigma values that differ by less than 5×10⁻⁷ will round to the same value and cause SCHISM to abort ("Weird side (3)").

### schimpy internal representation (`SchismLocalVerticalMesh`)

- **Left-justified**: `sigma[i, 0:Ni]` holds the `Ni` valid values; `sigma[i, Ni:]` is `NaN`.
- **kbps is 0-based**: `kbps[i] = nvrt − Ni`. File writes `kbps + 1`.
- **Same sigma convention**: values are in [−1, 0], bottom-first to surface-last within the valid slice.
- **Sentinel mapping**: on read, −9.0 → NaN. On write, NaN positions are filled with −9.0.

### Pipeline internal convention (lsc2_v2.py)

- **Sigma in [0, 1]**: 0 = surface, 1 = bottom. Left-justified: `sigma[i, 0] = 0`, `sigma[i, Ni-1] = 1`.
- **Conversion to SCHISM format**: negate (`-sigma`) then `flip_sigma` (reverses non-NaN slice). Result: [−1, ..., 0, NaN, ...] — bottom-first, matching `SchismLocalVerticalMesh` expectations.
- **Gotcha**: `flip_sigma` reverses the non-NaN portion. If the input is a palindrome (e.g., `[0, -1, 0]`), the flip is a no-op and the malformed array passes through unchanged. This happens when dry nodes (depth ≤ 0) have uninitialized interface values that clip to symmetric patterns.

### Key invariants

1. After the full pipeline, each node's sigma must be **strictly monotonically increasing** when written to disk.
2. Monotonicity must survive the `%14.6f` rounding — consecutive sigma must differ by at least ~1×10⁻⁶.
3. The z-coordinate check in SCHISM is at **side midpoints** (averaging two nodes), so both nodes of every edge must individually be monotone.


