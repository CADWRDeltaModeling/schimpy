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
