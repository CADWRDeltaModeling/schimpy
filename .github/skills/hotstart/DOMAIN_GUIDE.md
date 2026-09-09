# SCHISM hotstart — domain guide

Concepts, configuration semantics and behavior needed to build or repair a
hotstart. Organized by topic. Confidence is labelled throughout:

- **Established** — read from source or measured directly.
- **Convention** — how this project does it; not enforced by code.
- **Implementation** — specific to schimpy's current implementation.
- **Provisional** — plausible, not confirmed.
- **VERIFY AGAINST SOURCE** — check before relying on it.
- **OPEN QUESTION** — genuinely unresolved.

---

## 1. What a hotstart contains

**Established.** A schimpy-written hotstart for a `[TEM, SAL]` run carries:

```
time  iths  ifile  nsteps_from_cold          clock
eta2                                          free surface, per node
idry  idry_s  idry_e                          wet/dry flags: node, side, element
su2  sv2                                      horizontal velocity at sides
we                                            vertical velocity
q2  xl  dfv  dfh  dfq1  dfq2                  turbulence closure state
cumsum_eta  z                                 accumulated elevation, z coordinate
tr_nd  tr_nd0                                 tracers at nodes  (node, nVert, ntracers)
tr_el                                         tracers at elements (elem, nVert, ntracers)
```

Expectations on exact contents are extracted in `schism_hotstart.py`

Dimensions observed on a real Bay-Delta case: `node`, `elem`, `side`, `nVert`,
`ntracers`, `one`. Tracer order follows the `modules` list, so with
`modules: [TEM, SAL]`, tracer 0 is temperature and tracer 1 is salinity.

**Note.** `z` and `cumsum_eta` are written but SCHISM recomputes the vertical
coordinate from `eta2` and the bathymetry, so a "minimum layer thickness of 0"
seen in the raw `z` array is usually sub-bottom padding below `kbp`, not a real
degenerate column. **VERIFY AGAINST SOURCE** if this ever becomes load-bearing.

---

## 2. Clock semantics

**Established.** For an `ihot=2` continuation, the two dates mean different
things and are easy to transpose:

- `run_start` — the **original** run's nominal time origin. It does not move.
- `date` — the moment the restart represents.

`schism_hotstart` derives `time`, `iths`, `nsteps_from_cold` and the
`time_origin_of_simulation` attribute from these plus `time_step`.

Worked case: `run_start: 2020-09-30`, `date: 2021-10-05`, `time_step: 90` gives

```
time                  = 31968000      (elapsed seconds)
iths                  = 355200        (= 31968000 / 90)
nsteps_from_cold      = 355200
time_origin_of_simulation = 2020-09-30T00:00
```

The clock is **not** reset at the restoration or restart date for most use cases
although this would be possible by setting ihot=2 in `param.nml` while running schism.
Continuations stay on the original run's timeline so that time series inputs keyed to elapsed
seconds remain valid.

**Convention.** Output filename encodes both: `hotstart.<yyyymmdd>.<iths>.nc`. The run directory 
requires the file to be named or linked as `hotstart.nc`. The descriptive name is a staging convention
bit 

---

## 3. Wet and dry

**Established (implementation).** Node wet/dry flags depend on the elevation
initializer. When elevation comes from a prior `hotstart_nc`, target nodes that
coincide with source nodes retain the source `idry` value. Nodes without a
matched source flag, including novel nodes, are evaluated on the target grid:

```python
idry   = np.where(dp + eta <= h0, 1, 0)     # unmatched target node
idry_s = idry[edges].max(axis=1)            # side dry if either node is dry
idry_e = max(idry over element nodes)       # element dry if any node is dry
```

Consequences worth internalizing:

- A prior hotstart's node flag is preserved where source and target nodes match.
- Elsewhere, wetness is derived from target `dp + eta` against `h0`.
- The comparison is `<= h0`, so a column exactly at `h0` is **dry**.
- Dryness is contagious upward: one dry node dries its sides and its elements.
- `h0` comes from `param.nml` if available, may be overridden in the hotstart
  yaml, and defaults to `0.01`.

SCHISM's own dry test matches: `zcoor` aborts with `ZCOOR: dry location` when
`dp(inode) + eta2(inode) <= h0`.

### The dry-node elevation convention

**Convention, and a frequent source of confusion.** On dry ground the generated
initial condition places the surface *just below the bed*, typically
`eta = -z - freeboard` with `freeboard = 0.01`. That value is **not a physical
water level**. It exists only to make the node reliably dry.

This matters beyond the hotstart: anything that reads `eta2` at a node and
treats it as a water surface — notably hydraulic structures reading elevation at
their reference nodes — will consume the artifact as if it were real head. 

---

## 4. Initializers

**Established (implementation).** Every variable takes exactly one `initializer`.
Transferring from a prior run (`hotstart_nc`) is only one of several, and for a
cold start it is not used at all.

| Initializer | Source of values | Typical use |
|---|---|---|
| `simple_trend` | A constant, or an expression in `x`, `y`, `z` | Ocean constants, zero velocity, dry-terrain formulas |
| `obs_points` | Continuous station observations, `data` + `variable` | Delta and marsh tracers, where station coverage is dense |
| `extrude_casts` | Vertical CTD casts, `station` + `data` + `variable` | Bay and estuary tracers from USGS cruise / Polaris transects |
| `text_init` | A gr3-format file, `data_source` | Reusing a generated field such as `elev.ic` |
| `hotstart_nc` | A prior hotstart, plus its grids | Continuation and changed-grid transfer |
| `patch_init` | Dispatches to other initializers by region | Anything spatially heterogeneous |
| `schout_nc` | — | **Not implemented; the method body is `pass`.** Do not configure it. |

### `simple_trend` expressions

Accepts a scalar or a string expression evaluated over node coordinates, where
`z` is the node depth. `max` and `min` are rewritten to their NumPy equivalents,
so `max(0.97, -z-0.01)` is valid and means "ambient 0.97, but never above the
bed minus a centimetre".

This is the plain way to initialize new terrain dry: `values: -z-0.1` places the
surface a decimetre below the bed everywhere in the region. It is an alternative
to the `max_blw_bed` mechanism, not a companion to it — see §5.

**Gotcha (implementation).** The key under the initializer is read by
`get_value`, which returns `list(dict_obj.values())[0]` — the **first value,
whatever the key is called**. Both `value:` and `values:` appear in working
examples and both function. A misspelled key therefore fails silently rather than
erroring, and a second key would be ignored.

### Cast and observation initializers

`extrude_casts` takes a station table (`station:`) and a cast/transect table
(`data:`), plus the `variable:` column name as it appears in that file — the
names are dataset-specific, e.g. `"Temperature (Degrees Celsius)"` versus
`Salinity`. `obs_points` takes 1D station observations with `data:` and a
lower-case `variable:` such as `'temperature'`.

**Convention.** The usual Bay-Delta split is: ocean a constant, bay and estuary
from cruise casts, delta and marsh from continuous station data.

### Two `patch_init` patterns

Both are legitimate and they differ in how coverage is guaranteed. Recognize
which one a configuration is using before editing it.

| | Cold-start pattern | Changed-grid pattern |
|---|---|---|
| `regions_filename` | A **shapefile**, region names in a `region` attribute | A polygon **yaml** |
| Coverage | Only special regions need polygons | An explicit `domain` polygon covering everything |
| `allow_overlap` | Usually `False` | `True`, because each new region overlaps `domain` |
| Precedence | Regions are disjoint, so it does not arise | List order in `patch_init.regions`; last match wins |

**VERIFY AGAINST SOURCE.** Example documentation states that under the shapefile
pattern, areas not covered by any polygon are automatically assigned the region
name `other`. The corresponding handling in `geo_tools.partition_check` appears
commented out. Confirm before relying on an unlisted `other` region.

`regions_filename` accepting either a shapefile or a polygon yaml is the single
most useful thing to know here; a configuration written against one is not
mechanically portable to the other.

---

## 5. Transferring onto a changed grid

### `patch_init` and region precedence

**Established (implementation).** `patch_init` assigns a different initializer to
each spatial region.

- Regions come from a polygon yaml via `regions_filename`; it is consumed
  directly and should not be rendered to a region gr3 first.
- The file contains one **domain** polygon guaranteeing complete coverage, plus
  one distinctly named polygon per new area. Complete coverage is what satisfies
  `partition_check`'s sufficiency test.
- Each new-area polygon overlaps the domain, so `allow_overlap: true` is
  required.
- **Precedence follows the order of `patch_init.regions`, not the order of
  features in the polygon file.** The last matching configured region wins.
  List `domain` first, then the new areas.
- New-area polygons should stay mutually disjoint; overlap permission is only
  for their intended overlap with the domain.

There is no global or base initializer in this pattern. Every region has a peer.

### Novel nodes

**Established (implementation).** A target node with no counterpart in the source
mesh, within `novel_node_tol`, is *novel*. Novel nodes are handled inside
`interp_from_mesh`, scoped by the region polygon.

Two mechanisms apply, and they are distinct:

1. **Wet-donor redirection.** Novel nodes have their donor reselected from
   **wet source nodes only**. A newly excavated node beside a dry levee should
   not inherit that levee's state. This applies to *all* variables, not only
   elevation.
2. **Below-bed flooring** via `max_blw_bed`, described next.

**Historical note.** A separate mask for novel-node extrapolation was once
thought necessary. It is not — the `patch_init` regions already carry that
information. A function `fix_novel_elevation` existed and was deleted.

### `max_blw_bed`

**Established (implementation).** It is a **free-surface constraint**, and its
sign is the opposite of what the name suggests to most readers.

```python
floor = -dp[novel] - max_blw_bed
vout[novel] = np.maximum(vout[novel], floor)
```

So it is a **lower bound**: the interpolated surface may sit at most
`max_blw_bed` below the bed. With `max_blw_bed: 0.01` a floored node has
`H = dp + eta = -0.01`, i.e. dry by one centimetre.

Placement rules, enforced with explicit errors:

- **Required** on the elevation `hotstart_nc` initializer.
- **Rejected** on any other variable — it is meaningless for tracers.
- Not valid at the top-level `hotstart` block.
- Must be non-negative.
- Applies to **novel nodes only**, not to the whole field.

Do not confuse the floored value `H = -max_blw_bed` (dry, deliberate) with
`H = +h0` (marginally wet, usually accidental).

### Two ways to initialize new terrain

**Established.** These are alternatives, not companions, and mixing them is a
common source of confusion.

| Approach | How | When it fits |
|---|---|---|
| **Formula** | Give the new region a `simple_trend` such as `values: -z-0.1` | The whole region should be dry by construction; simplest and self-evident |
| **Transfer plus floor** | Give the region `hotstart_nc` or `text_init`, and floor novel nodes with `max_blw_bed` | The region should inherit real state where a counterpart exists, and only fall back to dry where none does |

The formula approach ignores the source entirely, so no donor question arises.
The transfer approach needs `max_blw_bed` precisely because some nodes have no
wet counterpart to inherit from.

### Interpolation settings

**Convention.** Tracers, velocity components and TKE typically use the same
source hotstart with `method: nearest` and a `distance_threshold`. Elevation
usually differs because it is the variable carrying the region policy.

---

## 6. Variables: what should and should not be transferred generically

**Established / Provisional mix.**

| Variable | Treatment | Note |
|---|---|---|
| `eta2` | Region-aware; the decision variable | Carries the wet/dry and new-terrain policy |
| `tr_nd`, `tr_nd0`, `tr_el` | Generic nearest-neighbour is normally fine | Watch for out-of-range values that break the equation of state |
| `su2`, `sv2`, `we` | Generic transfer acceptable | Large inherited velocities at a barely-wet new node are suspect |
| `q2`, `xl`, `dfv`, `dfh`, `dfq1`, `dfq2` | Generic transfer | `xl` must stay strictly positive; the k-kl closure divides by it |
| `idry` | Conditional | Preserve the source flag at nodes matched through elevation `hotstart_nc`; compute unmatched nodes from target `dp + eta` and `h0` |
| `idry_s`, `idry_e` | Derived | Recomputed from the completed target-node `idry` flags |
| `z`, `cumsum_eta` | Written, recomputed by SCHISM | Not a useful diagnostic target |

**Important.** Salinity must not go negative. The UNESCO equation of state takes
a square root of salinity, so a negative value produces NaN density and then NaN
velocity. This is a cheap and worthwhile check even though it was *not* the cause
in the case examined here. **VERIFY AGAINST SOURCE** for the specific EOS variant
compiled in your build.

Slightly negative *temperature* is not equivalent and is often benign — check
whether the source already had it before treating it as an interpolation defect.

### `z` semantics — a retracted confusion

**Established, previously misunderstood.** It was once believed that
`simple_trend`'s `z` disagreed with the gr3 `z` for elevation. It does not.
`elevation` uses `node2D` centering, whose `vgrid` is `node_z`, i.e. `dp`. They
agree. The disagreement is real only for 3D variables, where `vgrid` is
`build_z` output.

---

## 7. Vertical grid coupling

**Established in principle; case-specific in practice.** The hotstart is written
against a specific `vgrid`. With LSC2 (`ivcor=1`) each node carries its own level
count matched to its depth at generation time.

Therefore:

- If bathymetry changes materially, the vgrid should generally be refit against
  the new depths, and the hotstart regenerated against the new vgrid.
- Reading LSC2 sigma against a *different* depth does not by itself produce NaN,
  because the sigma values are still monotonic — it produces a badly resolved
  column, not an invalid one.
- A subtle trap: depth enforcement (dredging, breach cutting) usually happens in
  a **second** preprocessing pass, *after* the base grid and its vgrid were
  built. The vgrid then predates the cuts.

**OPEN QUESTION.** No automated check exists for whether a post-enforcement
depth change is large enough to require refitting the vgrid. Currently a
judgement call informed by comparing depths before and after enforcement.

---

## 8. How a bad initial state manifests at run time

This is the single most useful piece of knowledge in this skill.

### NaN is globalized by the elevation solver

**Established, read from source.** SCHISM solves the free-surface wave equation
with a Jacobi-preconditioned conjugate gradient (`solve_jcg`). The iteration
computes `rdotr`, `rdotz` and `alpha` as **global** `mpi_allreduce` sums, then
applies

```fortran
x = x + alpha*p        ! updates every node
r = r - alpha*sp
```

Consequences:

1. A single poisoned node anywhere makes the global dot products NaN.
2. `alpha` becomes NaN, and `x = x + alpha*p` makes **the entire elevation field**
   NaN in one iteration.
3. The guard `if(alpha==0.d0) call parallel_abort('JCG: division by zero')`
   does **not** catch NaN, because `NaN == 0` is false. It passes through
   silently.
4. The convergence test `rdotr <= rtol2*rdotr0` is also false for NaN, so the
   solver runs to `mxitn` and reports **"JCG did not converge"** before returning
   a NaN solution.

**Therefore: the element or node named in a downstream abort is essentially
random.** It is wherever the code first touched the poisoned field. In the case
examined, the aborting element was a clean rectangle roughly 9 km from the edited
region, wet, with entirely sane state.

**Diagnostic value:** `fort.33` reporting non-convergence, and the first step at
which it appears, is far more informative than the abort message.

### Reading the common aborts

| Abort | Means | Does *not* mean |
|---|---|---|
| `IBILINEAR: No roots` with NaN arguments | A NaN coordinate reached a quad shape-function solve | The element is degenerate |
| `IBILINEAR: Abnormal instances` / `Out of bound` with finite numbers | Genuinely a geometry or trajectory problem | — |
| `ZCOOR: dry location` | `dp + eta <= h0` where the code expected wet | — |
| `JCG did not converge` | Either genuine stiffness **or** NaN contamination | — |

The `itag` in an `IBILINEAR` message identifies the calling site, which
distinguishes backtracking from initialization.

**Rule of thumb:** NaN in the *arguments* means propagated garbage; wild but
finite values mean a real local problem.

---

## 9. Interaction with gates and structures

The hotstart does not configure structures, but structures constrain it.

**Established.** Blocked structure faces impose requirements on the initial
state:

- Structure reference nodes read `eta2` directly. They should be **wet**, or
  they feed the dry-node artifact into the flow calculation.
- Structure face nodes are commonly asserted to be wet in validation, since a
  dry structure face is usually a setup error.
- Elements inside an active block are excluded from the momentum solve, from
  horizontal viscosity and from the Shapiro filter. State there is not smoothed.

**Established, and a correction worth recording.** The `block_nudge` parameter is
**not** a blend between a prescribed velocity and a "barrier-free" momentum
solution. The source comment is explicit — the stored array is *"vel. at previous
step (for hydraulics etc)"*. The update is a first-order relaxation:

$$u^{n+1} = \beta\,v_{\text{th}} + (1-\beta)\,u^{n}$$

This is an IIR low-pass with fixed point exactly at $v_{\text{th}}$, converging
monotonically as $(1-\beta)^n$. With $\beta = 0.05$ that is an e-folding of about
20 steps. It is unconditionally stable for $0 \le \beta < 1$, which is precisely
the range enforced by the `block_nudge` guard at initialization.

So a small `block_nudge` is a **smoothing timescale, not weak enforcement**. An
earlier reading of this as "5% enforced, 95% free-running" was wrong, and any
instability argument built on it should be discarded.

**Practical implication for hotstarts:** a closed structure genuinely holds, so a
head difference across it will persist rather than blow through. That makes a
large initial head across a closed structure a *modeling* question rather than an
automatic numerical hazard.

---

## 10. Abrupt versus gradual inundation

**Convention.** Two different experiments, and the hotstart differs:

- **Abrupt** — the new area is connected at `t0`. Its initial elevation should be
  continuous with the water it connects to, or the run begins with a dam break.
- **Gradual** — the new area is isolated at `t0` behind a structure and opens
  later on a schedule. Its initial elevation is then free to differ from ambient,
  because the structure holds it.

In the gradual case the hotstart should show the area in its *pre-opening* state.
Confirm which is intended; they are not interchangeable and the choice changes
what "correct" means for every check downstream.

---

## 11. Validation: necessary versus sufficient

**Convention, learned the hard way.** A validation script that checks

- clock fields and time origin,
- all data variables finite,
- region elevations matching their source file,
- structure nodes wet,

can pass completely while the underlying grid is wrong. Every one of those
assertions held for a hotstart built on a grid whose restoration bathymetry was
absent.

Add checks that test **intent**, not just internal consistency:

- Does the depth-enforcement step report a plausible number of modified nodes?
- Does every node that should belong to a region actually resolve to it?
- Are there perched water columns — large `H` at nodes whose bed just moved?
- Did the count of marginally-wet nodes change sharply versus the source?
- Does the terrain being relied upon actually vary, or is it a constant fill?

---

## 12. Provisional and unresolved

**OPEN QUESTION — why an exactly-zero prescribed transfer flow behaved
differently from a tiny non-zero one.** A run that aborted with a globalized NaN
ran successfully after changing a prescribed structure flow from `0` to `1e-4`.
No mechanism was ever found. The implicit matrix is identical in both cases —
both flux terms are pure right-hand-side contributions — so only the RHS shifted,
by a physically negligible amount. Candidate explanations proposed and then
discarded: structurally singular matrix rows (tested, none), and weak structure
enforcement (based on a misreading of `block_nudge`). The case became moot when
the grid was found to lack bathymetry. **Treat a tiny numerical nudge that
"fixes" a blowup as masking, not repair.**

**OPEN QUESTION — significance of marginally-wet node counts.** A rebuilt
hotstart had roughly 4.5× as many wet nodes with `H < 0.02` as its source
(3710 vs 815). Suggestive of an ill-conditioned initial state but never shown to
be causal. Worth measuring; not yet worth acting on alone.

**Provisional.** Sub-bottom padding in the `z` array explains apparent
zero-thickness layers. Not confirmed against the writing code.
