# Gated island inundation — summary

Tools for introducing a restored island to a SCHISM grid behind a hydraulic structure and
opening it gradually. Landed in `17dcf99` plus the core changes listed below.

## Workflow

The restoration is applied as a **modification of a finished grid**, not as a change to the
base preprocessing. Nothing in the original preprocessing configuration is edited, so the
unmodified grid stays reproducible and the restoration can be added, revised or dropped
without rebuilding it.

### 1. Work in a directory of its own

Create a working directory beside the finished preprocessing output and refer back to it
with relative paths. A worked layout:

```
prepro_out_2021/                       # first round: full preprocessing, already built
  hgrid.gr3                            #   corrected bathymetry and everything else
  vgrid.in.3d
  restoration_gate_mod/                # second round: this restoration, self contained
    inundate_breaches.yaml             #   the specification, hand written
    main_inundation.yaml               #   small prepare_schism config, hand written
    depth_enforce_inundate.yaml        #   generated
    elev_inundate.yaml                 #   generated
    hydraulic_structures_inundate.yaml #   generated
    inundate_regions.yaml              #   generated
    th_files/dated/*.th                #   generated, human readable
    th_files/elapsed/*.th              #   generated, what SCHISM reads
    prepro_out_inundate/               #   second round output
      hgrid.gr3  hgrid.ll  elev.ic  hydraulics.in
    hotstart_nudging/
      hotstart_2021grid.yaml           #   continuation config, hand written
      <date>/source_for_hotstart/      #   the prior run's restart and its grid
      hotstart.<date>.<iths>.nc        #   generated
```

Only three files are written by hand: the specification, the small `prepare_schism`
configuration, and the hotstart configuration.

### 2. Write the specification

One file describes every restoration area. Required per island: `name`, `polygon`, and a
non-empty `breaches` list. Required per breach: `name`, `min_depth`, `max_depth`, and
geometry as either `left`/`right` or a `pathway` (a `pathway` longer than two points also
needs `gate_span`). Everything else is optional.

Breach points are ordered **looking into the island**; that ordering is what tells the
generator which side is the interior. Depths are positive down, matching gr3.

### 3. Choose the gate type

The structure is declared per breach. In the Dutch Slough case each breach is a **throttled
weir**: crest at 1.5 m so it stays wet at most tides, `width` 20 m, and `op_downstream`
`1.0e-5` throttling the rating curve to a trickle that tops the interior pools up against
evaporation. `op_upstream` is `0.0`, so the pools never drain back to the channel. See
[Structure type and schedule](#structure-type-and-schedule) for why a weir is usually
preferable to a `transfer` and how the opening is scheduled.

### 4. Generate

```
sch inundate_island --config inundate_breaches.yaml --hgrid ../hgrid.gr3 --out-dir .
```

Writes the four yaml artifacts, plus a dated and an elapsed `.th` for every breach that
names a `breach_date`. Structures are validated against the mesh as they are built.

### 5. Second, simplified preprocessing round

```
prepare_schism main_inundation.yaml
```

Applies the dredge to the finished `hgrid.gr3`, writes `elev.ic`, and merges the generated
structures into the base `hydraulics.in`. See
[Simplified preprocessor integration](#simplified-preprocessor-integration).

### 6. Continuation hotstart from a prior run

```
create_hotstart hotstart_2021grid.yaml
```

Transfers a restart from a prior run onto the restoration grid: the domain keeps its prior
state through `hotstart_nc`, while each restoration region takes the generated `elev.ic`
levels through `text_init`. See [Hotstart integration](#hotstart-integration).

### 7. Stage the run

Copy `prepro_out_inundate/` outputs and `th_files/elapsed/*.th` into the run directory. The
`.th` files must sit beside `hydraulics.in`; SCHISM opens them by structure name.

### What each round is for

| round | input | output | rerun when |
|---|---|---|---|
| first, full preprocessing | DEMs, base configuration | `hgrid.gr3`, `vgrid.in.3d`, all gr3 | bathymetry or mesh changes |
| second, this tool | finished `hgrid.gr3` + specification | dredged grid, `elev.ic`, `hydraulics.in`, `.th` | breach geometry, gate or schedule changes |
| hotstart | second-round grid + a prior restart | `hotstart.nc` | either of the above changes |

The second round reads the first round's grid, so **any change to the base grid invalidates
everything downstream of it**. Regenerate in order.

## Entry points

```
sch inundate_island --config breaches.yaml --hgrid hgrid.gr3 --out-dir .
```

Geometry options: `--edge-length-factor`, `--freeboard`, `--ambient`, `--prefix`,
`--validate/--no-validate`.
Schedule options: `--breach-date`, `--run-start`, `--th-dir`, `--epoch-start`,
`--ramp-hours`, `--ramp-steps`, `--ramp-target`.

```python
from schimpy.inundate_island import (
    build_inundation_artifacts,   # worker: islands + mesh -> dict of yaml structures
    write_inundation_inputs,      # dict -> four yaml files
    build_breach_schedule,        # one schedule entry -> dated DataFrame
    write_breach_timeseries,      # schedules -> dated and elapsed .th
    validate_structures,          # build structures, check nodes against the polygons
    generate_inundation_inputs,   # file-facing: spec + hgrid -> written files
)
```

`build_inundation_artifacts` and `build_breach_schedule` do no file I/O. The CLI and
`generate_inundation_inputs` handle reading and writing.

## Input

```yaml
run_start: 2020-09-30          # model time origin; required if any breach_date is set

islands:
  - name: upper_andrus
    polygon: infer            # or a preprocessor polygon mapping with 'vertices'
    breaches:
      - name: main_1
        breach_date: 2022-04-15        # optional; enables the .th schedule
        left:  [625596.0, 4228143.9]   # ordered looking INTO the island
        right: [625856.1, 4228287.9]
        min_depth: -0.20               # undisturbed ground at the rim
        max_depth: 4.09                # dredged invert at the centre
        major_axis_len: 240.0          # optional; else 2 x factor x local edge length
        dredge_pad: 30.0               # optional; else one local edge length
        structure:
          type: weir
          configuration:
            n_duplicates: 1
            elevation: 1.5
            width: 20.0
            coefficient: 0.6
            op_downstream: 1.0e-5
            op_upstream: 0.0
```

Depths are positive down throughout, matching gr3.

### Required and optional keys

| level | required | optional |
|---|---|---|
| file | `islands` | `run_start`, `breach_date` |
| island | `name`, `polygon`, `breaches` | `breach_date` |
| breach | `name`, `min_depth`, `max_depth`, and `left`+`right` or `pathway` | `gate_span`, `major_axis_len`, `dredge_pad`, `structure`, `breach_date` |

`gate_span` is required when `pathway` has more than two points. `run_start` is required
once any breach resolves a `breach_date`, since the elapsed `.th` is measured from it.

`breach_date` belongs to the **breach**. A multi-year run can combine restoration sites
that come online years apart, so the date is carried per breach even when every breach at
one site shares it. It may also be set once per island or once at the top of the file;
the narrowest setting wins:

```
breach.breach_date  ->  island.breach_date  ->  file breach_date  ->  --breach-date
```

A breach that resolves no date gets no `.th` and no `use_time_series` flag, so scheduled
and unscheduled breaches can be mixed in one specification.

## Output

Four yaml files, all sharing one polygon per breach:

| file | contents |
|---|---|
| `hydraulic_structures_inundate.yaml` | `structures:` only — **no `nudging`**, see below |
| `depth_enforce_inundate.yaml` | `type: none`, `attribute: ellipse(...)` |
| `elev_inundate.yaml` | domain, then island, then breach entries **in that order** |
| `inundate_regions.yaml` | regions for hotstart `patch_init` |

Plus, for each breach carrying a `breach_date`:

| file | contents |
|---|---|
| `th_files/dated/<struct_name>.th` | schedule with an ISO datetime column and a header, for reading and diffing |
| `th_files/elapsed/<struct_name>.th` | the same rows in elapsed seconds, **no header** — this is what SCHISM reads |

`<struct_name>` is `<island>_<breach>`, matching the name in `hydraulics.in`. SCHISM opens
`<struct_name>.th` from the run directory, and only when `use_time_series: 1` appears in
the structure's configuration — the generator sets that flag automatically for scheduled
breaches.

## Structure type and schedule

The structure type defaults to `transfer` with `flow: 0.0`, which is a prescribed discharge
and nothing else. A **weir is usually the better instrument**, for three reasons visible in
`hydraulic_structures.F90`:

1. **It has a real closed state.** `depth_flow = min(height, max_elev - elev)`, and a
   non-positive `depth_flow` returns zero flow immediately. A crest above the current water
   level closes the structure by hydraulics rather than by prescription, and it reopens on
   its own when the tide tops the crest.
2. **It can be relaxed.** `read_struct_ts` reads `ttt, install, nduplicate, op_down, op_up,
   elev, width` for `WEIR` and `PIPE`, so crest elevation and width are schedulable. A
   `transfer` schedules only `ttt, install, prescribed_flow` — a discharge with no
   geometric meaning.
3. **It self-limits.** Once the interior rises past the crest the submergence term engages
   and attenuation falls to zero as the two sides equalise, so the interior cannot be
   overfilled.

`op_downstream` and `op_upstream` are pure linear multipliers applied last, independent of
the geometry, which makes them the clean throttle. Setting `op_upstream: 0.0` hard-blocks
the reverse direction, so an interior pool fills but never drains.

The generated schedule has three phases:

| when | rows | what |
|---|---|---|
| `epoch_start`, default 2005-01-01 | 1 | the configured setting, standing in for the infinite past |
| the `ramp_hours` before the breach, default 24 | `ramp_steps`, default 12 | `op_downstream` ramps to `ramp_target`, default 0.6 |
| `breach_date` | 1 | `install = 0` |

Deinstalling is what actually opens the breach: it zeroes `isblock_el` and `isblock_sd`, so
the elements rejoin the momentum solve and regain horizontal viscosity and the Shapiro
filter, which are all skipped on active block sides. The ramp exists to make that switch
small; `block_nudge` relaxes what remains.

Two cautions:

- **`coefficient` is not schedulable.** It is static in `hydraulics.in`. What ramps is the
  operating coefficient `op_downstream`.
- **Reference nodes must stay wet.** `elev_up` and `elev_down` are `eta2` at the structure's
  two reference nodes. On dry ground the generated initial condition sets `eta = -z -
  freeboard`, just below the bed, which is not a physical water level; a dry reference node
  would feed that artifact in as real head. The generator places the interior reference node
  inside the dredged pool, which is wet by construction, and `validate_structures` confirms
  the reference pair lands inside the dredge polygon.

## Simplified preprocessor integration

The generated files can be applied to an existing `hgrid.gr3` without rerunning the
full Bay-Delta preprocessing stack. Put a small `prepare_schism` configuration beside
the four generated files:

```yaml
prepro_output_dir: prepro_out_inundate

imports:
  - schimpy.ellipse.ellipse

mesh:
  mesh_inputfile: ../hgrid.gr3
  depth_enforcement:
    polygons:
      include:
        - depth_enforce_inundate.yaml
  ll_outputfile: hgrid.ll
  gr3_outputfile: hgrid.gr3

gr3:
  elev.ic: !include elev_inundate.yaml

hydraulics:
  include:
    - ../../bay_delta_template/hydraulic_structures.yaml
    - hydraulic_structures_inundate.yaml
  outputfile: hydraulics.in

# Optional: regenerate the vertical grid against the dredged bathymetry.
# vgrid: !include ../../bay_delta_template/vgrid_v2.yaml
```

Run it from that directory:

```
prepare_schism main_inundation.yaml
```

The integration has four important details:

1. Import `schimpy.ellipse.ellipse`, not only `schimpy.ellipse`. The generated depth
   expressions call `ellipse(...)`; importing only the module binds a non-callable name.
2. List the base hydraulic structures first and the generated structures second. The YAML
   include merge concatenates their `structures` lists. The generated file intentionally
   omits `nudging`, so the base file's scalar value is retained exactly once.
3. Use the generated `elev_inundate.yaml` as a complete `elev.ic` polygon definition. Its
   domain/island/breach ordering must not be rearranged.
4. Vertical-grid generation is optional. Include the existing vgrid recipe when the
   modified depths should receive a newly fitted vgrid; omit it when the modification run
   should retain the existing model vgrid.

The simplified run writes a modified `hgrid.gr3`, `hgrid.ll`, `elev.ic`, and combined
`hydraulics.in` under `prepro_output_dir`. `inundate_regions.yaml` is not consumed by this
step; use it later as the spatial regions input for hotstart `patch_init`. The `.th` files
are not consumed either \u2014 copy them to the run directory beside `hydraulics.in`.

## Hotstart integration

Use the modified `hgrid.gr3`, the generated `elev.ic`, and `inundate_regions.yaml` to
transfer a prior restart onto the restoration grid. `patch_init` accepts the generated
YAML polygon file directly; do not render it to a separate region GR3. It contains one
domain polygon plus one distinctly named polygon for each restoration area. The domain
guarantees complete coverage, satisfying `partition_check`'s sufficiency check. Each
restoration polygon overlaps the domain, so use `allow_overlap: true`.

There is no global or base initializer in this use of `patch_init`. Every region has a
peer initializer: the domain uses `hotstart_nc`, while each restoration region uses
`text_init` from the same generated `elev.ic`. This preserves the prior water level in
the domain and applies the generated dry-terrain level, `-dp-0.01`, plus the pooled breach
levels inside each restoration. Other variables can transfer directly from the prior
hotstart with the usual `hotstart_nc` interpolation settings.

Overlap precedence follows the order of `patch_init.regions`, not the feature order in
`inundate_regions.yaml`: the last matching configured region wins. List the domain first
and every restoration region after it. The generated region file carries the same warning.
The restoration polygons should remain mutually disjoint; `allow_overlap` is needed only
for their intentional overlap with the domain polygon.

```yaml
hotstart:
  date: 2021-10-05
  run_start: 2020-09-30
  time_step: 90
  hgrid_input_file: ../prepro_out_inundate/hgrid.gr3
  vgrid_input_file: ../../vgrid.in.3d
  vgrid_version: "5.10"

  elevation:
    initializer:
      patch_init:
        smoothing: false
        regions_filename: ../inundate_regions.yaml
        allow_overlap: true
        allow_incomplete: false
        regions:
          # Precedence is list order: domain first, restoration patches later.
          - region: domain
            initializer:
              hotstart_nc:
                data_source: ./source_for_hotstart/hotstart.nc
                source_hgrid: ./source_for_hotstart/hgrid.gr3
                source_vgrid: ./source_for_hotstart/vgrid.in.3d
                source_vgrid_version: "5.10"
                max_blw_bed: 0.01
                novel_node_tol: 0.001
          - region: dutch_emerson
            initializer:
              text_init:
                data_source: ../prepro_out_inundate/elev.ic
          - region: dutch_emerson_pool
            initializer:
              text_init:
                data_source: ../prepro_out_inundate/elev.ic
          - region: dutch_gilbert
            initializer:
              text_init:
                data_source: ../prepro_out_inundate/elev.ic
```

Repeat the `text_init` entry once for each distinctly named restoration region emitted by
the generator. `max_blw_bed` belongs inside the elevation `hotstart_nc` block; it applies
only to novel nodes supplied by that initializer. It is not valid at the top-level
`hotstart` block or for other variables.

Wet/dry state follows the elevation initializer. For elevation supplied by `hotstart_nc`,
nodes that coincide with source-grid nodes retain the source `idry` flag; other target
nodes are evaluated from the initialized elevation and target depth using
`H = dp + eta2 <= h0`. Elevation supplied by another initializer, such as the restoration
`text_init` patches above, uses that same target-grid test. Side and element flags are
then derived from the completed node flags. No separate wet/dry configuration is needed.

For an `ihot=2` continuation, retain the original run's nominal start in `run_start` and
put the restart date in `date`. `schism_hotstart` writes elapsed `time`, `iths`,
`nsteps_from_cold`, and `time_origin_of_simulation` from those values. For the dates and
90-second step above, the result is `time=31968000`, `iths=355200`,
`nsteps_from_cold=355200`, and origin `2020-09-30T00:00`; it does not reset the clock at
the restoration date.

## Design decisions worth knowing

**The pool level is derived, not specified.** `schimpy.ellipse.ellipse` clamps its radial
term at 1, so every node inside a breach polygon ends at `min_depth` or deeper. That makes
`-min_depth - freeboard` a constant water level meeting the island's `-z - freeboard`
exactly at the rim. It is applied as `type: min`, a **lower bound**, placed after the
island entry — it pools the dredged opening while leaving ambient water and dry ground
untouched. This removed the need for a `z_orig` (pre-dredge depth) mechanism entirely.

**Ordering in `elev.ic` is load bearing.** Domain, then island, then breach. The generated
file carries a comment saying so.

**One polygon per breach, reused verbatim.** Dredge, initial condition and hotstart region
all use the same vertices, so consistency is structural rather than conventional.

**`polygon: infer` walks the element graph** from each breach, treating elements a breach
line crosses as barriers. A set of breaches together pens in an island no single one
closes. The walk starts from the **island-side reference node the structure itself will
use** (`down_path` side), so the fill and the structure cannot disagree about which side is
the island. Every breach must reach the same connected component; disagreement is reported
with the coordinates of the offending node.

**Orientation rule.** `find_two_neighboring_node_paths` computes `norm = (-tan[1], tan[0])`
and sends nodes with `dot(norm, node - p1) > 0` to `down_path`. With `left`/`right` ordered
looking into the island, `down_path` is the interior. Pinned by tests; do not change
without them.

**The dredge runs past the opening** by `dredge_pad` at each end, and `major_axis_len`
defaults to `2 x edge_length_factor x local_edge_length` — the factor is the reach *each
way*. Both exist because the structure's end nodes and reference pair would otherwise land
on the taper, where the dredge has already returned to `min_depth`.

## Core changes made along the way

| commit | change |
|---|---|
| `12ea3b3` | `schism_yaml` list-include merge: concat lists, merge dicts, scalars keep-first with a warning. Previously two files each setting `nudging: 0.05` silently produced `0.1`. |
| `0593824`, `6546a27` | `create_structures` node-pair validation: equal length, crossing edge present, quad-or-two-triangles between consecutive pairs. Previously `zip` silently truncated unequal paths. |
| `37260f6`, `e03c552` | `wet_dry_check` now uses `H = dp + eta2 <= h0`. Previously it tested `build_z` output against `>= 0`, which never fired, so `idry`/`idry_s`/`idry_e` were never set. |
| `3a2853d` | Novel-node handling moved into `interp_from_mesh`, scoped by `inpoly`. Wet-donor selection now applies to **all** variables. `max_blw_bed` moved into the `hotstart_nc` block, required on elevation, rejected elsewhere, applied to novel nodes only. `fix_novel_elevation` deleted. |

Two findings raised during this work were **wrong and are retracted**:

- `simple_trend`'s `z` does **not** disagree with the gr3 `z` for elevation. `elevation` is
  `node2D` centering, whose `vgrid` is `node_z`, i.e. `dp`. They agree. The disagreement is
  real only for 3D variables, where `vgrid` is `build_z` output.
- A separate mask for novel-node extrapolation is unnecessary. The `patch_init` regions
  already carry that information.

## Limitations

- **`min_depth` must describe ground, not the levee crest.** The `min` bound only leaves
  the channel alone while the pool level stays below ambient.
- **Interior holes are lost.** `_fill_outline` keeps the longest ring, so a restoration
  area enclosing a borrow pit would have the hole discarded. Do not use the polygon for
  area.
- **The 50 percent fill guard is a heuristic**, and fires only after the walk has covered
  much of the mesh.
- **Inferred rings are not simplified** — one vertex per boundary node, so the files are
  awkward to read by hand. Shapely `simplify` was planned and omitted.
- **The dredge is specified, not computed.** Nothing checks `min_depth`/`max_depth` against
  surrounding bathymetry or the flow the structure should pass.
- **`_find_reference_node` is private** and called from outside `SchismSetup`.

## Not yet done

- bdschism wrapper and its `pyproject.toml` / `__main__.py` entries.
- A SCHISM run confirming the gates never trip the dry-structure-node rule.
- Nothing verifies that the dredge leaves the vertical grid valid. When the cuts move
  depths materially, regenerate the vgrid in the second round and point the hotstart
  `vgrid_input_file` at it.

## Map

- `schimpy/inundate_island.py` — generator, validation, CLI
- `tests/test_inundate_island.py` — geometry, fill, layering, validation, CLI
- `tests/test_breach_pool.py` — the `min_depth` pool identity on a four-node fixture
- `tests/testdata/inundate_island/breaches_example.yaml` — annotated sample spec
- `docsrc/notebooks/inundate_island.ipynb` — worked synthetic case, end to end
- `tests/test_hotstart_z_semantics.py`, `tests/test_hotstart_wet_dry.py`,
  `tests/test_hotstart_donors.py`, `tests/test_structure_node_paths.py` — core changes
