# SCHISM hotstart — worked examples

Examples chosen to teach judgment, especially between cases that look alike.

---

## 1. "The run dies immediately with an IBILINEAR NaN. Should I look for a pathological mesh?"

**Context.** A hotstart continuation on a grid that recently gained a restoration
area. Abort names a specific element.

**Distinction to recognize.** Whether the abort arguments are **NaN** or
**finite-but-extreme**. That single fact decides the entire investigation.

**Preferred reasoning.** NaN arguments mean the coordinate arrived poisoned. The
elevation solver globalizes NaN in one iteration, so the named element is
arbitrary. Confirm cheaply by inspecting it — if it is a well-shaped, wet element
far from the edits, the location is noise. Redirect to: `fort.33` for
non-convergence and the first failing step, then to what could poison a node.

**Bad response.** Spending effort on element quality metrics, then reporting that
the mesh is fine and the cause is unknown. The mesh check is a two-minute
rule-out, not the investigation.

---

## 2. "The source hotstart has 400,606 nodes, the target has 420,435."

**Context.** Changed-grid extension after adding a restoration area.

**Distinction to recognize.** Node-count mismatch is diagnostic in a *same-grid*
continuation and expected in a *changed-grid extension*.

**Preferred reasoning.** Establish which variant you are in first. Then verify
the difference is confined to the intended footprint: compare coordinates and
depths where the meshes correspond, and confirm the extra nodes fall inside the
new area.

**Bad response.** Flagging it as an error, or conversely waving it through
without checking that the difference is localized.

---

## 3. "There are negative temperatures. Nearest-neighbour shouldn't be able to do that."

**Context.** 81 negative temperature values in a rebuilt hotstart.

**Distinction to recognize.** Whether the method can *create* values or only
*copy* them.

**Preferred reasoning.** Nearest-neighbour copies. Trace the donors: if donor
distance is ~0 and the donor holds the same value, the source already had it and
the transfer is faithful. Separately, note that negative *temperature* is not the
NaN hazard that negative *salinity* is.

**Bad response.** Treating any out-of-range value as evidence of a broken
transfer, or conflating the two tracers' risk profiles.

---

## 4. "Should the new restoration area start wet or dry?"

**Context.** Grid gained a diked area with an excavated pond and a breach.

**Distinction to recognize.** This is a **modeling decision**, not something to
infer from the files.

**Preferred reasoning.** Ask, and frame the options concretely: dry with the
surface just below the bed; a designed pond at a specified level; or continuous
with ambient channel water. Tie the answer to whether the area is connected at
`t0` — if a structure isolates it, a non-ambient level is legitimate; if it is
open, a non-ambient level is a dam break.

**Bad response.** Picking whichever value makes the validation script pass, or
copying the previous case's number without checking that the connectivity story
is the same.

---

## 5. "Nodes in the new area have 4.5 m of water. That seems like a lot."

**Context.** The design intent was a shallow pond, "a couple of metres".

**Distinction to recognize.** *Which* nodes are deep. Pond nodes and cut/breach
nodes have different expected depths and different provenance.

**Preferred reasoning.** Separate them before answering. In the case examined the
pond was exactly as intended (median 1.5–1.75 m, max 2.4 m); every node deeper
than 3 m was a breach-cut node. Those split cleanly by region: nodes resolving to
`domain` carried the channel surface, nodes resolving to the new region carried
the pond surface, giving a 1.4 m step across a 15–25 m gap. The aggregate "max
depth" number concealed all of that.

**Bad response.** Reporting the global maximum and concluding the excavation
depth is wrong.

---

## 6. "Validation passes but the run still fails."

**Context.** Clock correct, all variables finite, region elevations match their
source, all structure nodes wet.

**Distinction to recognize.** Internal consistency versus correspondence to
intent.

**Preferred reasoning.** Those checks cannot detect a grid that never received
its bathymetry. Escalate to intent checks: how many nodes did depth enforcement
actually modify, relative to how many are in the footprint? Does terrain that
should be graded actually vary? Are features represented by a plausible number of
nodes? Did the earlier pass write its diagnostics?

**Bad response.** Concluding the initial state is fine and moving the
investigation into the solver. In the case examined, 103 modified nodes out of
12,201, and 13 cm of relief across 6,850, were the answer.

---

## 7. "I changed a prescribed flow from 0 to 0.0001 and it works now. Why?"

**Context.** A NaN blowup that disappears under a physically negligible change.

**Distinction to recognize.** "Runs" is not "fixed".

**Preferred reasoning.** Check whether the change alters the implicit matrix or
only the right-hand side. If only the RHS, and by a negligible amount, it cannot
stabilize an unstable state — it perturbs a marginal one. Say so plainly, keep
the real investigation open, and avoid inventing a mechanism to justify it.

**Bad response.** Constructing a plausible-sounding instability story to explain
it. This happened here, twice, and both explanations were wrong.

---

## 8. "Does switching the gate from a transfer to a weir change the hotstart?"

**Context.** Structure type changed in the preprocessing configuration.

**Distinction to recognize.** What lives in the hotstart versus what lives in the
runtime structure configuration.

**Preferred reasoning.** No configuration change is needed — the hotstart carries
elevation, velocity, tracers and wet/dry flags; structure type is read at runtime
from `hydraulics.in`. But the hotstart must still be **regenerated** if `hgrid`,
`vgrid`, or any file feeding a region initializer changed. Separate "needs
editing" from "needs rebuilding"; they are different answers.

**Bad response.** Either editing the hotstart yaml unnecessarily, or answering
"no change" and omitting the rebuild.

---

## 9. "The vgrid is from the base grid but the dredge came later."

**Context.** Two-pass preprocessing: base grid and vgrid first, depth enforcement
second.

**Distinction to recognize.** Between "produces NaN" and "poorly resolved".

**Preferred reasoning.** LSC2 sigma stays monotonic against a different depth, so
this is not a NaN source. It is a resolution problem: a node's level count was
fitted to its old depth. Measure how far the affected depths actually moved, then
decide. Flag the ordering consequence — regenerating the vgrid invalidates any
hotstart already built against the old one.

**Bad response.** Asserting confidently in either direction without measuring the
depth change.

---

## 10. "Set the new area's surface below the bed so it starts dry."

**Context.** Using `max_blw_bed` on the elevation initializer.

**Distinction to recognize.** The direction of the constraint, and its scope.

**Preferred reasoning.** It is a lower bound: `floor = -dp - max_blw_bed`, applied
with `np.maximum`. `max_blw_bed: 0.01` yields `H = -0.01` — dry by a centimetre.
It applies to **novel nodes only**, is **required** on elevation `hotstart_nc`,
and is **rejected** on other variables and at the top level. Verify the resulting
`H`, since `H = -0.01` (deliberately dry) and `H = +0.01` (marginally wet) look
almost identical in a summary but mean opposite things.

**Bad response.** Placing it at the top level or on tracers, or assuming it
guarantees a minimum water depth.

---

## 11. "Which regions should `patch_init` list, and in what order?"

**Context.** Domain plus several restoration areas.

**Distinction to recognize.** Precedence follows the **order of
`patch_init.regions`**, not the feature order in the polygon file.

**Preferred reasoning.** List `domain` first so later restoration regions win the
overlap. Set `allow_overlap: true`, because each restoration polygon overlaps the
domain by construction. Keep restoration polygons mutually disjoint. Give every
region a peer initializer — there is no global fallback in this pattern. Point
`regions_filename` at the polygon yaml directly rather than rendering it to a
region gr3.

**Bad response.** Reordering the polygon file expecting precedence to follow, or
adding a base initializer that silently competes.

---

## 12. "Build me a hotstart for a new run starting 2021-04-20."

**Context.** No prior run to continue from.

**Distinction to recognize.** This is a **cold start**, and `hotstart_nc` plays no
part in it. The work is choosing an initializer per region per variable.

**Preferred reasoning.** Reach for `patch_init` over a region shapefile and
dispatch by data availability rather than by geography for its own sake: a
constant via `simple_trend` where the ocean is effectively uniform, vertical
cruise casts via `extrude_casts` where transect data exist, and continuous
station data via `obs_points` where the station network is dense. Velocity
components are normally `simple_trend: 0.0`. Read the `variable:` column names
out of the actual data files rather than reusing them from another case.

**Bad response.** Assuming every hotstart is a transfer, and asking the user for
a source hotstart that does not exist.

---

## 13. "Start the new flooded area dry."

**Context.** The grid gained a flooded island or restoration area.

**Distinction to recognize.** Two mechanisms achieve this and they are
alternatives, not companions.

**Preferred reasoning.** If the region should simply be dry, give it a formula —
`simple_trend` with `values: -z-0.1` places the surface below the bed everywhere
in the region and asks no donor question at all. Use `hotstart_nc`/`text_init`
plus `max_blw_bed` only when the region should *inherit* real state where a
counterpart exists and fall back to dry only where none does. Choosing the
transfer path for a region that has no counterpart anywhere is needless
machinery.

**Bad response.** Configuring `max_blw_bed` alongside a `simple_trend` region, or
assuming the transfer path is required because the area is new.

---

## 14. "Copy the region configuration from the other study."

**Context.** Two hotstart configurations, both using `patch_init`.

**Distinction to recognize.** Whether `regions_filename` points at a **shapefile**
or a **polygon yaml**. The two patterns guarantee coverage differently.

**Preferred reasoning.** The shapefile pattern names regions through a `region`
attribute and typically runs with `allow_overlap: False`. The polygon-yaml
pattern used for changed-grid work carries an explicit `domain` polygon, runs
with `allow_overlap: True`, and depends on list-order precedence. Copying region
blocks between them without carrying the coverage strategy produces silently
uninitialized or wrongly-won nodes.

**Bad response.** Porting `regions:` entries verbatim and leaving `allow_overlap`
at whatever the donor file had.

---

## 15. "The example has a create_hotstart.py — should I write one?"

**Context.** An older worked example ships a Python driver script beside the yaml.

**Distinction to recognize.** Between the current interface and a superseded one.

**Preferred reasoning.** Drive `create_hotstart` from the yaml. The driver script
is legacy; its presence marks an example as predating current practice, which
also means its other contents deserve a second look rather than trust.

**Bad response.** Reproducing the script pattern in new work because an example
shows it.

## 16. "Is `run_start` the restart date?"

**Context.** Building an `ihot=2` continuation.

**Distinction to recognize.** `run_start` is the **original** run's origin;
`date` is the restart moment.

**Preferred reasoning.** Recover `run_start` from the source hotstart's
`time_origin_of_simulation` attribute where possible rather than asking. Then
sanity-check the derived clock: `iths = nsteps_from_cold = (date − run_start) /
time_step`. The clock does not reset at the restart date, which is what keeps
elapsed-time forcing files valid.

**Bad response.** Setting `run_start` to the restart date. It produces a
plausible-looking file with a silently wrong clock.
