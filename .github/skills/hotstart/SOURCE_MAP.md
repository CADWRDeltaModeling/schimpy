# SCHISM hotstart — source map

Only sources actually established during the work that produced this skill.
Paths are relative to their repository root unless noted. Line numbers are
omitted because they drift; symbol names are given instead.

---

## schimpy — hotstart implementation

### `schimpy/schism_hotstart.py`

Primary implementation. Reached from the `create_hotstart` console entry point.

| Symbol / key | Role |
|---|---|
| `wet_dry_check` | Completes node `idry`: matched nodes initialized from an elevation `hotstart_nc` retain source flags, while unmatched nodes use target `dp + eta <= h0`. It then derives sides and elements by taking the maximum node flag. |
| `self.h0` | Minimum depth. Defaults to `0.01`, overridden from `param.nml`, then from the hotstart yaml if present. |
| `_validated_max_blw_bed` | Enforces placement of `max_blw_bed`: required on elevation `hotstart_nc`, rejected on other variables, must be non-negative. |
| `_source_wet_mask` | Wet mask on the **source** mesh, `(dp_src + eta_src) > h0`. |
| `redirect_to_wet_donors` | Reselects donors for novel nodes from wet source nodes only. |
| `interp_from_mesh` | Where novel-node handling lives, scoped by `inpoly`; also applies the `max_blw_bed` floor `-dp - max_blw_bed` via `np.maximum`. |
| `simple_trend` | Constant or `x`/`y`/`z` expression; rewrites `max`/`min` to NumPy equivalents |
| `obs_points` | 1D station observations |
| `extrude_casts` | Vertical casts, e.g. USGS cruise / Polaris transects |
| `text_init` | Values from a gr3-format file |
| `patch_init` | Region dispatch to other initializers |
| `hotstart_nc` | Transfer from a prior hotstart |
| `schout_nc` | **Stub — body is `pass`.** Not usable. |
| `get_value` | Returns `list(dict_obj.values())[0]`, so the key name under an initializer is not validated |
| `create_dataarray` | Derives `tr_el` and `tr_nd0` from node tracers for prism centering. |
| variable-name maps | `eta2` ↔ `elevation`; `idry` ↔ `wetdry_node`; `idry_e` ↔ `wetdry_elem`; `elevation` uses `node2D` centering. |

**Answers:** how wet/dry is decided; what `max_blw_bed` does and where it is
legal; how novel nodes get donors; which initializers exist and which is a stub;
which variables are derived rather than transferred; how yaml variable names map
to netCDF names.

### Hotstart configuration yaml (per study, not in the repo)

Observed keys, with the meaning established here:

```
out_dir, hotstart.output_fn
hotstart.date          restart moment
hotstart.run_start     ORIGINAL run origin (not the restart date)
hotstart.time_step
hotstart.modules       e.g. [TEM, SAL]; sets tracer order
hotstart.hgrid_input_file, vgrid_input_file, vgrid_version
<variable>.initializer.patch_init:
    smoothing, regions_filename, allow_overlap, allow_incomplete, regions[]
<variable>.initializer.hotstart_nc:
    data_source, source_hgrid, source_vgrid, source_vgrid_version,
    max_blw_bed, novel_node_tol, distance_threshold, method
<variable>.initializer.text_init:
    data_source            e.g. a generated elev.ic
```

Variables seen configured independently: `elevation`, `temperature`, `salinity`,
`velocity_u`, `velocity_v`, `velocity_w`, `tke`.

**Answers:** what a working continuation configuration looks like; which keys
belong to which initializer.

---

## schimpy — supporting modules

| Path | Role | Answers |
|---|---|---|
| `schimpy/schism_mesh.py` (`read_mesh`) | Reads `.gr3`, including `elev.ic`-style files where the third column is the value | Node coordinates, depths, elements, edges |
| `schimpy/schism_setup.py` (`create_schism_setup`, `structure_node_paths`, `apply_polygons`) | Mesh-aware setup; resolves structure node paths; applies polygon attribute expressions | Which nodes a structure occupies; what depth enforcement produces |
| `schimpy/geo_tools.py` (`partition_check`) | Assigns mesh nodes to named regions; enforces coverage sufficiency | Whether a node resolves to the region you expect |
| `schimpy/schism_yaml.py` | Include-aware yaml with list-merge semantics: lists concatenate, dicts merge, scalars keep-first with a warning | Why two included files both setting a scalar do not sum |
| `schimpy/prepare_schism.py` | The preprocessor driver; holds the recognized structure configuration keys | Valid configuration key names |
| `schimpy/schism_structure.py` | Writes `hydraulics.in`; `use_time_series` becomes the "time series enabled" flag | Structure file layout and property names |

---

## schimpy — tests

| Path | Covers |
|---|---|
| `tests/test_hotstart_wet_dry.py` | The target-grid `H = dp + eta2 <= h0` rule, source-flag override behavior, and side/element derivation |
| `tests/test_hotstart_donors.py` | Wet-donor selection for novel nodes |
| `tests/test_hotstart_z_semantics.py` | `z` semantics; the `node2D` vs 3D distinction |
| `tests/test_structure_node_paths.py` | Node-pair validation for structures |
| `tests/test_breach_pool.py` | The `min_depth`-derived pool level identity |

**Answers:** the intended contract for each behavior, and which behaviors are
pinned against regression.

---

## schimpy — documentation

| Path | Role |
|---|---|
| `inundate_island_summary.md` | Workflow manual for grid modification and the continuation hotstart that follows it. Contains the `patch_init` precedence rules, the "no global initializer" pattern, the clock worked example, and a record of retracted findings. |
| `docsrc/notebooks/inundate_island.ipynb` | End-to-end synthetic worked case |
| `AGENTS.md` | Package scope and conventions; keeps Bay-Delta policy out of schimpy |

---

## Applied examples — BayDeltaSCHISM

`BayDeltaSCHISM/examples/hotstart/`. These are useful applied configuration
examples, but they are not self-contained test cases: the shared target grids,
vertical grids and source hotstarts are not distributed. Bay-Delta-specific
data choices belong here; schimpy stays domain-neutral.

| Path | Case |
|---|---|
| `examples/basic/` | Cold start, TEM and SAL, region-dispatched: ocean constant, bay/estuary from cruise casts, delta and marsh from station observations. `elev.ic: !include elev.yaml`. `run_start: default`. |
| `examples/flooded_island/` | Legacy newly flooded area. Demonstrates `patch_init` precedence, but its elevation `hotstart_nc` predates the required `max_blw_bed` key and its clock readme is obsolete. Use the generated inundation workflow for new work. |
| `examples/hotstart_from_previous_hotstart/` | Continuation from a prior restart |
| `examples/sed/`, `examples/bio/`, `examples/tracer_age/` | Additional module sets |
| `data_in/` | `usgs_cruise_station.csv`, `polaris_transect_<date>.csv`, `all_stations_temperature_<date>.csv`, `all_stations_salinity_<date>.csv`, plus sediment and ICM/cosine merged tables |
| `shapefile/` | Region shapefiles; the `region` attribute supplies the names matched in `patch_init.regions` |

**Answers:** what a cold start actually looks like; which initializer suits which
part of the estuary; the exact `variable` column names for cruise and station
data; how a flooded area is started dry without `max_blw_bed`.

**Caveats.**

- Four cases carry legacy Python drivers: `basic`, `flooded_island`, `sed`, and
  `tracer_age`. The drivers used to inject modules and write the dataset
  directly. That pattern is **superseded**: put required modules and inputs in
  yaml and run `create_hotstart <config.yaml>` or
  `sch create_hotstart <config.yaml>`.
- Module-specific yaml files currently lag the CLI contract in places. In
  particular, module lists formerly supplied by drivers must be explicit, and
  every elevation `hotstart_nc` must set `max_blw_bed`.
- `hotstart_from_previous_hotstart/` is a template for the downstream
  `hot_from_hot` wrapper, not a standalone `create_hotstart` example.
- Prefer implementation and tests for interface facts; treat these examples as
  demonstrations of initializer selection and regional policy until refreshed.

---

## SCHISM source — behavior that constrains the hotstart

Paths relative to the SCHISM source tree.

| Path | Symbol | Answers |
|---|---|---|
| `src/Hydro/solver_subs.F90` | `solve_jcg` | Why NaN goes global: `mpi_allreduce` dot products for `rdotr`/`rdotz`/`alpha`, and `x = x + alpha*p` over all nodes. Also why the `alpha == 0` guard misses NaN and why non-convergence is logged instead. |
| `src/Hydro/misc_subs.F90` | `ibilinear`, `quad_shape`, `zcoor` | Meaning of each `IBILINEAR` abort variant; the `ZCOOR: dry location` test `dp + eta2 <= h0` |
| `src/Hydro/bktrk_subs.F90` | backtracking `quad_shape` calls | Which `itag` corresponds to which call site |
| `src/Hydro/schism_step.F90` | block-face velocity update; "Save vel. at previous step (for hydraulics etc)" | That `block_nudge` relaxes toward the **previous step**, not a barrier-free solution; that block elements are excluded from the momentum solve, viscosity and the Shapiro filter |
| `src/Hydro/schism_init.F90` | hydraulics initialization | The `block_nudge` range guard `0 <= block_nudge < 1` |
| `src/Core/hydraulic_structures.F90` | `calc_struc_flow`, `load_structures`, `read_struct_ts`, `init_struct_time_series`, `irreg_time_history_advance` | That structures read `eta2` at reference nodes; which fields are schedulable per type; that `.th` files are opened as `<struct_name>.th` and read list-directed |

---

## Run configuration

| File | Keys established | Answers |
|---|---|---|
| `param.nml` (`param.nml.clinic` / `.tropic`) | `h0`, `dt`, `ihot`, `ihydraulics`, `nramp_elev`, `inunfl` | The wetting/drying threshold and time step the hotstart must be consistent with |
| `hydraulics.in` | structure count, `block_nudge` (nudging factor), per-structure node pairs and reference nodes | Which nodes must be wet; where structures sit |

**Run output for diagnosis**

| File | Answers |
|---|---|
| `fort.33` | JCG convergence log — the reliable NaN indicator |
| `mirror.out` | How many steps the run survived |
| `nonfatal_*`, `fort.12` | Warnings written before a fatal abort |

---

## Preprocessing outputs consumed by a hotstart

| Artifact | Role |
|---|---|
| `hgrid.gr3` (post-enforcement) | Target mesh and bathymetry |
| `vgrid.in.3d` / `vgrid.in` | Vertical grid; must correspond to the depths in use |
| `elev.ic` | Source for `text_init` region elevations |
| regions yaml | `regions_filename` for `patch_init`; contains a domain polygon plus one per new area |
| `diagnostics/` | Presence of `dem_misses` output indicates the bathymetry step ran |

---

## Ad hoc diagnostic scripts

Not part of the package; written per investigation and worth recreating rather
than preserving. The useful patterns were:

- probe one element or node by index, remembering **SCHISM is 1-based and the
  Python preprocessor is 0-based**;
- diff two `hgrid.gr3` files for node count, coordinate movement and depth
  change, then classify changed nodes by region;
- scan a hotstart for non-finite values, negative salinity, marginal water
  columns and extreme velocities;
- compare a rebuilt hotstart against its source on the same statistics;
- trace donor nodes for suspicious values;
- partition the mesh by region and assert that intended nodes resolve correctly.

A study-local `validate_hotstart.py` established the assertion style: clock
fields, finiteness, region elevations equal to `elev.ic`, and structure nodes
wet. See `EDGE_CASES.md` for why that set is necessary but not sufficient.

---

## Environment

Work was performed in the `schism` conda environment. The SCHISM environment is
the expected place to run schimpy for formal or informal testing.
