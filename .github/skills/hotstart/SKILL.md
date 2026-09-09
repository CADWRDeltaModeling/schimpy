---
name: schism-hotstart
description: >-
  Create, modify, extend, or debug SCHISM hotstart files with schimpy. Use for cold-start
  initialization, continuation from a prior run, transferring a restart onto a changed grid,
  initializing newly added terrain (restoration areas, breaches, levee removal), and for
  diagnosing NaN blowups or wet/dry problems that appear immediately after a hotstart.
---

# SCHISM hotstart creation and modification

## When to use this skill

Use it when the request involves the initial state of a SCHISM run: building
`hotstart.nc`, changing what a hotstart contains, moving a restart onto a
different grid, or explaining why a run fails at or near startup.

Triggering requests look like:

- "build a hotstart for this run" / "make a continuation from the prior run"
- "the grid changed, transfer the old restart onto it"
- "initialize the new restoration area / breach / flooded island"
- "the run dies immediately with NaN / IBILINEAR / ZCOOR / JCG did not converge"
- "why is this node wet (or dry) at the start?"
- "should this area start wet or dry?"
- "validate this hotstart before I burn a run"

Do **not** use it for gate operation schedules, breach geometry, or structure
type selection. Those belong to island-inundation work; this skill only covers
the initial state and the constraints structures place on it.

## Workflow variants to distinguish

Identify which one applies before doing anything. They have different required
inputs and different failure modes.

| Variant | Situation | Distinguishing feature |
|---|---|---|
| **Cold start** | No prior state | Values come from constants, coordinate formulas, cruise casts and station observations, usually dispatched per region |
| **Continuation, same grid** | Restart a run mid-stream | Source and target meshes identical; mostly a clock exercise |
| **Changed-grid extension** | Prior run exists, mesh or bathymetry changed | Needs spatial transfer plus a policy for nodes with no counterpart |
| **New terrain initialization** | Grid gained wettable area (restoration, breach, levee removal) | Needs an explicit wet/dry decision for the new area |

Cold start is at least as common as the others — do not treat it as a stub.
The last two are usually combined, and are where the subtlest failures live.

**Convention.** Drive `create_hotstart` from the yaml. A `create_hotstart.py`
driver script appears in older examples and is no longer the intended pattern;
do not add one, and treat its presence as a sign the example predates current
practice.

## Resolution rules

Resolve these before editing configuration. Classify each fact first; only the
third category justifies interrupting the user.

### Discoverable — go find it, do not ask

| Fact | Where |
|---|---|
| `h0` (min depth for wetting/drying) | `param.nml`; schimpy also accepts it in the hotstart yaml; default 0.01 |
| `dt`, `ihot`, `ihydraulics` | `param.nml` |
| Model time origin | `time_origin_of_simulation` attribute on the prior hotstart, or the prior run's config |
| Elapsed clock of the source restart | `time`, `iths`, `nsteps_from_cold` variables in the source hotstart |
| Node/element/side counts, depths | `hgrid.gr3` on both source and target |
| `nVert`, tracer count | `vgrid` and the `modules` list |
| Whether the mesh actually changed | Diff node count, node coordinates and depths between source and target `hgrid.gr3` |
| Whether structures exist and where | `hydraulics.in`, plus the structure yaml that generated it |
| Region names available to `patch_init` | The regions yaml referenced by `regions_filename` |
| Whether a prior preprocessing step ran | Presence and contents of its `diagnostics/` output |
| Which `patch_init` pattern a config uses | Whether `regions_filename` is a `.shp` or a polygon yaml |
| Available observation and cast inputs | The study's `data_in/` — station tables, cruise transects, and the `variable` column names inside them |

### Safely inferable — state the inference, then proceed

- `iths` and `nsteps_from_cold` follow from `(date - run_start) / time_step`.
- Descriptive output naming `hotstart.<yyyymmdd>.<seconds>.nc` is the local
  convention; follow it when extending an existing setup.
- In a changed-grid extension, the unchanged domain takes `hotstart_nc` and any
  newly added region takes an explicit initializer.
- Tracers, velocity components and TKE normally share the same source and the
  same interpolation settings as each other.
- A hotstart must be regenerated whenever `hgrid`, `vgrid`, or any file feeding
  a region initializer changes.

### Requires user input — genuine modeling decisions

Do not guess these. They change the physics, not just the plumbing.

- **`run_start` versus `date`.** For `ihot=2`, `run_start` is the *original*
  run's nominal origin, not the restart moment. Getting this wrong silently
  produces a wrong clock. Confirm it if it is not recoverable from the source.
- **Should new terrain start wet or dry?** A freshly excavated area can be
  initialized dry, as a shallow pool, or at ambient channel level. These are
  different experiments.
- **What elevation should new terrain take?** Ambient water level, a designed
  pond level, or "just below the bed" are all defensible and mutually exclusive.
- **Abrupt or gradual inundation?** Whether the area is open at t0 or opens
  later via a gate schedule determines whether the hotstart should show it
  connected.
- **Should the vertical grid be regenerated?** If depths moved materially, LSC2
  should generally be refit — but this invalidates the hotstart and is a
  workflow decision.

## Assumptions that must never be made silently

- **Do not assume the abort location identifies the cause.** SCHISM's elevation
  solver globalizes NaN across the whole domain in one iteration, so the
  reported element is usually innocent and often far away.
- **Do not assume differing node counts between source and target is an error.**
  It is expected when the grid gained terrain.
- **Do not assume anomalous values were introduced by interpolation.** Check
  whether the source already contained them before blaming the transfer.
- **Do not assume "all values finite" means the hotstart is usable.** Finiteness
  is necessary, not sufficient.
- **Do not assume `eta` at a dry node is a physical water level.** By convention
  it is placed just below the bed and is an artifact.
- **Do not assume wet/dry flags are always recomputed.** An elevation
  `hotstart_nc` preserves source `idry` at coincident nodes; target-grid
  `dp + eta` is used where no source flag is supplied.
- **Do not assume a region initializer covered every node you intended.** Verify
  coverage; polygons routinely clip a row short.
- **Do not assume the grid you were handed is correct.** Verify that bathymetry
  actually exists where the work depends on it.
- **Do not treat a tiny numerical change that "makes it run" as a fix.** It is
  usually masking a marginal state.
- **Do not assume `hotstart_nc` is the initializer being asked for.** A cold
  start uses formulas, casts and station observations instead.
- **Do not copy a `variable:` name between datasets.** Cast and observation files
  use their own column naming.

## Supporting files

- `DOMAIN_GUIDE.md` — concepts, configuration semantics, variable-by-variable
  treatment, clock arithmetic, wet/dry rules, and how NaN propagates.
- `EDGE_CASES.md` — failure modes and misleading symptoms, with the diagnosis
  each one actually warrants.
- `EXAMPLES.md` — worked judgment calls, including superficially similar cases
  that need different responses.
- `SOURCE_MAP.md` — the code, configuration and tests that answer specific
  questions.

## Default working method

1. Establish the workflow variant and resolve the facts above.
2. Read the existing hotstart yaml before editing it; these configurations are
   order-sensitive and key placement matters.
3. Prefer small probe scripts over reasoning from configuration alone. Measuring
   the mesh and the resulting file is fast and settles most disputes.
4. Validate the produced file against intent, not just against finiteness.
5. When something contradicts your model of the code, read the source. Several
   confident-sounding conclusions in this domain have turned out to be wrong.
