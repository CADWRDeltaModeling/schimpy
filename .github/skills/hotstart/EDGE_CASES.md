# SCHISM hotstart — edge cases and failure modes

Corner cases that materially changed a diagnosis or a workflow. Ordered roughly
by how often they mislead.

---

## The abort location is not the problem location

**Situation**
A run dies shortly after a hotstart with something like
`IBILINEAR: No roots; 2 NaN 262636 NaN NaN`, naming a specific element.

**Tempting interpretation**
That element is degenerate. Go inspect its geometry, look for a pathological
quad, suspect the mesh near it.

**Why that can be wrong**
The arguments are `NaN`, not wild finite numbers. `delta`, `x` and `y` are all
NaN, which means a NaN *arrived*; the element did not create it. SCHISM's
elevation solver globalizes NaN across the entire domain in a single iteration,
because `alpha` comes from `mpi_allreduce` dot products and `x = x + alpha*p`
updates every node. Whichever element the code touches first afterwards is the
one that aborts. In the case examined, that element was a clean rectangle
(interior angles 91/84/93/92°, area 807 m²) about 9.4 km from the edited region,
wet, with entirely sane hotstart values.

**Preferred diagnosis or action**
Distinguish NaN arguments from finite-but-wild ones. For NaN, ignore the location
and look at `fort.33` for `JCG did not converge` and the first step it appears.
Then look for what could poison a single node: bad initial state, a degenerate
water column, or an out-of-range tracer. Inspect the element only to rule it out.

**Relevant implementation/docs**
`src/Hydro/solver_subs.F90` (`solve_jcg`); `src/Hydro/misc_subs.F90`
(`ibilinear`, `quad_shape`); `src/Hydro/bktrk_subs.F90` (`itag` values identify
the calling site).

---

## `if (alpha == 0.d0)` does not catch NaN

**Situation**
The JCG solver has an explicit divide-by-zero guard, so a NaN `alpha` should be
caught.

**Tempting interpretation**
The guard would have aborted with a clear message if the solve went bad.

**Why that can be wrong**
`NaN == 0` evaluates false in Fortran, so NaN passes straight through. The
convergence test `rdotr <= rtol2*rdotr0` is likewise false for NaN, so the solver
silently runs to `mxitn`, logs "did not converge", and returns a NaN solution.

**Preferred diagnosis or action**
Never treat the absence of a solver abort as evidence the solve was healthy.
Read `fort.33`.

**Relevant implementation/docs**
`src/Hydro/solver_subs.F90`, `solve_jcg`.

---

## Different node counts between source and target

**Situation**
The source hotstart has 400,606 nodes; the target grid has 420,435.

**Tempting interpretation**
Mismatch — something is wrong with the grid or the transfer.

**Why that can be wrong**
This is the *expected* signature of a changed-grid extension. The new grid gained
terrain the old one did not represent. Comparing raw counts tells you nothing by
itself.

**Preferred diagnosis or action**
Confirm the difference is confined to the region you intended to change. Compare
node coordinates and depths where the meshes *do* correspond, and check that the
added nodes fall inside the new-area footprint.

---

## Anomalous values that were already in the source

**Situation**
The rebuilt hotstart contains 81 negative temperature values, minimum −1.16 °C.

**Tempting interpretation**
Nearest-neighbour interpolation produced physically impossible values, so the
transfer is defective.

**Why that can be wrong**
Nearest-neighbour cannot invent a value — it copies one. Tracing the donors
showed the values were inherited verbatim from the source hotstart at donor
distance 0.00, at a handful of dry land nodes. The transfer was faithful.

**Preferred diagnosis or action**
Before blaming interpolation, trace suspicious values back to their donors and
check whether the source already contained them. If a method is a pure copy,
reason about it as a copy.

---

## Salinity sign versus temperature sign

**Situation**
Out-of-range tracer values in a hotstart.

**Tempting interpretation**
All out-of-range tracers are equally suspicious.

**Why that can be wrong**
They are not equivalent. Negative salinity is a genuine NaN generator because the
equation of state takes a square root of salinity. Slightly negative temperature
is often benign, especially at dry or sub-bottom levels.

**Preferred diagnosis or action**
Always check salinity for negatives explicitly; it is cheap. Treat negative
temperature as a question, not a defect, until you know where it came from.

---

## `max_blw_bed` has the opposite sign to most readings of its name

**Situation**
`max_blw_bed: 0.01` is set on the elevation initializer and the result is
`H = dp + eta = −0.01`.

**Tempting interpretation**
It caps how far below the bed the surface may go, or it forces a minimum water
depth.

**Why that can be wrong**
It is a **lower bound on the free surface**: `floor = -dp - max_blw_bed`, applied
as `np.maximum(vout, floor)`. It permits the surface to sit *at most*
`max_blw_bed` below the bed, producing a deliberately dry node with one
centimetre of negative depth. It applies to **novel nodes only**, and is
**required** on elevation `hotstart_nc` while being **rejected** elsewhere and at
the top level.

**Preferred diagnosis or action**
Read `H = -max_blw_bed` as "deliberately dry". Distinguish it sharply from
`H = +h0`, which is marginally wet and usually accidental.

**Relevant implementation/docs**
`schimpy/schism_hotstart.py`: `_validated_max_blw_bed`, and the flooring block
inside `interp_from_mesh`.

---

## Region polygons clipping one row short

**Situation**
A restoration area is initialized via `patch_init`, but some nodes inside the
feature received the domain's `hotstart_nc` value instead of the region's
`text_init` value.

**Tempting interpretation**
The initializer or the precedence order is broken.

**Why that can be wrong**
The precedence was correct; the polygon simply did not reach those nodes. In the
case examined, one row of each breach fell outside the region polygon and landed
in `domain`. The two rows were 15–25 m apart and ended up 1.40 m apart in
elevation — a step sitting exactly on a structure.

**Preferred diagnosis or action**
Verify coverage explicitly rather than trusting the polygon. Partition the mesh,
then assert that every node you care about resolves to the region you expect.
A quick tell: values that "should" be uniform within a feature coming out in two
distinct clusters.

**Relevant implementation/docs**
`schimpy.geo_tools.partition_check`; the regions yaml referenced by
`regions_filename`.

---

## Perched water columns from donor elevation inheritance

**Situation**
Nodes whose bed was lowered by excavation end up with 3–4.5 m water columns where
the source had 11 cm.

**Tempting interpretation**
Interpolation produced absurd depths.

**Why that can be wrong**
Elevation and depth were transferred independently and consistently. The donor
node had `dp = -0.883, eta = 1.594`, i.e. `H = 0.11 m` of water on a marsh. The
target node kept that surface but its bed had been cut to `dp = 2.823`, so
`H = dp + eta = 4.42 m`. Each half was right; the combination was not.

**Preferred diagnosis or action**
After any depth-changing step, check for nodes whose depth moved and whose
resulting `H` is implausible. The pattern is: inherited surface, new bed.

---

## Validation passing on a grid that has no bathymetry

**Situation**
The hotstart validation script passes — clock correct, all variables finite,
region elevations match `elev.ic`, all structure nodes wet — and the run still
dies.

**Tempting interpretation**
The hotstart is fine; the problem must be numerical, in the solver or the
structures.

**Why that can be wrong**
Those assertions test internal consistency, not intent. The grid in question had
no restoration bathymetry at all: the depth-enforcement step modified only 103
nodes out of 12,201 in the restoration polygons, and dry ground across 6,850
nodes spanned just 13 cm between the 5th and 95th percentiles — a flat plate, not
terrain. Every consistency check still passed.

**Preferred diagnosis or action**
Include intent checks. Look for: implausibly few modified nodes; near-zero
variance in terrain that should be graded; features (ponds, channels) represented
by a handful of nodes; and an empty `diagnostics/` directory where a previous run
produced DEM-miss output.

---

## An empty diagnostics directory is evidence

**Situation**
A preprocessing output folder exists and looks populated, but its `diagnostics/`
subdirectory is empty.

**Tempting interpretation**
Diagnostics are optional output; absence is uninteresting.

**Why that can be wrong**
The comparable full run wrote `dem_misses` files there. An empty directory means
the bathymetry step did not run at all in that pass — the strongest early signal
that the grid is not what was assumed.

**Preferred diagnosis or action**
Compare the output inventory of a partial pass against a known-good full pass.
Missing artifacts (no vgrid, no gr3 parameter files, no diagnostics) identify
which stages actually executed.

---

## A tiny numerical change that "fixes" the run

**Situation**
Changing a prescribed structure flow from `0` to `1e-4` turns a NaN abort into a
successful run.

**Tempting interpretation**
The zero value was the bug; the small non-zero value is the fix.

**Why that can be wrong**
The implicit matrix is identical for both — the relevant flux terms are pure
right-hand-side contributions — so only the RHS moved, by a physically negligible
amount (order `1e-6 m/s` across the face). A perturbation that small cannot
stabilize a genuinely unstable state; it can only shift a marginal one. Two
mechanisms were proposed and both failed: structurally singular matrix rows
(tested, zero found) and weak structure enforcement (based on a misreading of
`block_nudge`). The real defect was elsewhere.

**Preferred diagnosis or action**
Treat it as a signal that the state is marginal, and keep looking. Record it as
unexplained rather than resolved.

---

## `block_nudge` is a relaxation, not a blend

**Situation**
A structure face uses `block_nudge = 0.05`, and the update mixes 5% of a
prescribed velocity with 95% of a stored array.

**Tempting interpretation**
The structure is only 5% enforced, and 95% of the face flux is free-running
momentum — so a large head across a "closed" structure will blow through.

**Why that can be wrong**
The stored array is the velocity at the **previous step**, saved explicitly for
the hydraulics. The update is a first-order relaxation
`u^{n+1} = β·v_th + (1-β)·u^n`, an IIR low-pass whose fixed point is exactly
`v_th`, converging monotonically as `(1-β)^n` — about a 20-step e-folding at
β = 0.05. It is unconditionally stable for `0 ≤ β < 1`, which is exactly the
range the initialization guard enforces. Block faces are also excluded from the
momentum solve, so the stored value is pure filter state that the pressure
gradient never re-forces.

**Preferred diagnosis or action**
Read `block_nudge` as a smoothing timescale. A closed structure holds, so a large
initial head across one is a modeling question, not an automatic numerical
hazard.

**Relevant implementation/docs**
`src/Hydro/schism_step.F90`, the "Save vel. at previous step (for hydraulics
etc)" allocation and the block-face update; `src/Hydro/schism_init.F90`, the
`block_nudge` range guard.

---

## Marginally wet nodes at exactly `h0`

**Situation**
The minimum `H` over wet nodes is `0.010002`, essentially `h0`, and thousands of
nodes sit below `H = 0.02`.

**Tempting interpretation**
Either a bug, or definitely the cause of the blowup.

**Why that can be wrong**
It is consistent with the wet/dry rule (`idry = 1` when `dp + eta <= h0`), so
nodes just above `h0` are legitimately wet. The count rose about 4.5× versus the
source (3710 vs 815), which is suggestive, but it was never shown to be causal.

**Preferred diagnosis or action**
Measure it and compare against the source as a health indicator. Do not treat it
as a diagnosis on its own. **OPEN QUESTION.**

---

## Vertical grid predating depth enforcement

**Situation**
The hotstart references a vgrid built with the base grid, while breach or dredge
cuts were applied in a later preprocessing pass.

**Tempting interpretation**
Either it is fine because sigma is relative, or it is fatal because depths moved.

**Why that can be wrong**
Neither extreme. LSC2 sigma values remain monotonic when read against a different
depth, so this does not by itself produce NaN — it produces a poorly resolved
column. But each node's level count was fitted to its original depth, so a node
that moved from land to several metres of water is badly served.

**Preferred diagnosis or action**
Compare depths at the affected nodes before and after enforcement. If they moved
materially, regenerate the vgrid in the enforcement pass and repoint
`vgrid_input_file` — remembering that this invalidates any hotstart already built.

---

## Dry-node elevation consumed as real head

**Situation**
A structure's reference node sits on dry ground, and the structure computes flow
from elevations at its reference nodes.

**Tempting interpretation**
The elevation there is a water level like any other.

**Why that can be wrong**
On dry ground the generated initial condition sets `eta = -z - freeboard`, just
below the bed. It is an artifact chosen to guarantee dryness. If that exceeds the
real water level on the other side, the structure can be driven backwards — a dry
area appearing to drain.

**Preferred diagnosis or action**
Assert that structure reference nodes are wet. Where a generator places the
interior reference node inside a dredged pool, confirm that the pool is still wet
after any change to the pool depth or level.
