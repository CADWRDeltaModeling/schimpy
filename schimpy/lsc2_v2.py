
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LSC2 v2: Fully local + globally smooth vertical mesh construction for SCHISM.

This module implements the Step A–E pipeline:
  A) Compute continuous desired levels L* from a vertical size function.
  B) Quantize to integer levels with hysteresis and spatial gating.
  C) Build per-node target interface depths from the size function.
  D) Fit bottom & mid interfaces to targets with screened-Poisson smoothing,
     explicit bottom-thickness anchoring, and ordering/min-thickness guards.
  E) Final tiny stretch to match exact depth.

It preserves the idea of a general "vertical size function" abstraction.

Dependencies expected in the environment:
  - schimpy.schism_mesh: read_mesh, write_mesh
  - schimpy.schism_vertical_mesh: SchismLocalVerticalMesh, write_vmesh
  - laplace_smooth_data: smooth kernels (used here in wrappers)
"""

from __future__ import annotations
import logging
import pandas as pd
import numpy as np
import scipy.sparse as sp
import time
import math
from dataclasses import dataclass,field
from typing import Optional, Tuple, Callable, Dict

from schimpy.laplace_smooth_data import smooth_kernel_numba
from schimpy.schism_mesh import write_mesh

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Vertical size function abstraction
# -----------------------------------------------------------------------------

class VerticalSizeFunction:
    """Abstract base: maps level index -> cumulative depth for a given local depth.
    Implementations must be monotone in level index and approach D as level -> N.
    """

    def depth(self, levels: np.ndarray, depth: float, xycoords = None) -> np.ndarray:
        """Return cumulative depth (down from surface) at the provided level indices.
        Parameters
        ----------
        levels : np.ndarray (int or float)
            Level indices (0 at surface; last ~= total levels-1 at/near bed).
        depth : float
            Local total depth D (positive down).
        x : float, optional
            Extra parameter if needed by a custom function.
        """
        raise NotImplementedError("depth(levels, depth) must be implemented")


class WrapLegacyMeshFun(VerticalSizeFunction):
    """Wraps the legacy meshfun from lsc2.py which had signature depth(levels, depth, x)."""

    def __init__(self, legacy_meshfun):
        self._mf = legacy_meshfun

    def depth(self, levels: np.ndarray, depth: float, xycoords= None) -> np.ndarray:
        x = None
        return self._mf.depth(levels, depth, x)




class ScaledByFieldSizeFunction(VerticalSizeFunction):
    """Multiplies the level-to-depth mapping by a per-node scale factor.
    scale[i] > 1.0 => *smaller* dz (more levels) because cumulative depth grows slower.
    """

    def __init__(self, base: VerticalSizeFunction, scale_field: np.ndarray):
        self.base = base
        self.scale = scale_field

    def depth(self, levels: np.ndarray, depth: float, xycoords = None) -> np.ndarray:
        H = self.base.depth(levels, depth, coords=coords, node=node, context=context)
        s = 1.0
        if node is not None and 0 <= node < len(self.scale):
            s = float(self.scale[node])
        return H / max(s, 1e-9)


# --- New: σ-based size functions with dz cap and variants --------------------


class _ShapeAwareCountMixin:
    """Choose N (intervals) consistently with the shape ψ(s)."""
    count_spec = None  # {mode: max|top|bottom, value: meters}

    def set_count_spec(self, spec):
        self.count_spec = spec
        return self

    def _psi(self, s):
        raise NotImplementedError

    def required_intervals(self, D: float, count_spec: dict | None = None, Ncap: int = 200) -> int:
        import numpy as np, math
        if not np.isfinite(D) or D <= 0:
            return 1
        spec = count_spec or self.count_spec or {}
        mode = str(spec.get("mode", "max")).lower()
        val  = max(1e-6, float(spec.get("value", getattr(self, "dz_max", 1.5))))

        # >>> NEW: depth-aware ψ if the sizefun provides it
        psi_with_D = getattr(self, "_psi_with_depth", None)
        def psi_of(s):
            return psi_with_D(s, D) if psi_with_D is not None else self._psi(s)

        lo, hi = 1, max(1, Ncap)

        def violates(nint):
            s   = np.linspace(0.0, 1.0, nint + 1)
            psi = psi_of(s)
            dzs = D * np.diff(psi)
            if   mode == "top":    m = dzs[0]
            elif mode == "bottom": m = dzs[-1]
            elif mode == "max":    m = dzs.max()
            else:
                raise ValueError(f"Unknown count_spec.mode='{mode}' (expected: max|top|bottom)")
            # small tolerance so floating point can’t pass a barely-bad split
            return m > (val + 1e-9)

        while violates(hi) and hi < 10000:
            hi *= 2
        ans = hi
        while lo <= hi:
            mid = (lo + hi)//2
            if violates(mid): lo = mid + 1
            else:             ans = mid; hi = mid - 1
        return max(1, int(ans))
        
class SigmaCapSizeFunction(VerticalSizeFunction, _ShapeAwareCountMixin):
    """Pure-σ spacing with a per-node dz cap.
       N-1 = ceil(D/dz_max); H(ell) = D * (ell/(N-1))."""
    def __init__(self, dz_max: float = 1.5):
        self.dz_max = float(dz_max)

    
    # shape-aware hooks
    def _psi(self, s):
        import numpy as _np
        s = _np.asarray(s, float); return _np.clip(s, 0.0, 1.0)

    def set_count_spec(self, spec):
        self.count_spec = spec
        return self
    def depth(self, levels: np.ndarray, depth: float, xycoords=None) -> np.ndarray:
            D = max(float(depth), 1e-9)
            Nminus1 = max(1, int(math.ceil(D / self.dz_max)))
            # linear σ distribution
            return D * (np.asarray(levels, dtype=float) / Nminus1)


class SigmaPowerSizeFunction(VerticalSizeFunction, _ShapeAwareCountMixin):
    """σ spacing with dz cap and exponent p in (0,1) to bias toward bed."""
    def __init__(self, dz_max: float = 1.5, p: float = 0.7):
        self.dz_max = float(dz_max)
        self.p = float(p)

    
# shape-aware hooks
def _psi(self, s):
    import numpy as _np
    s = _np.clip(_np.asarray(s,float),0.0,1.0); return _np.power(s, self.p)

def set_count_spec(self, spec):
    self.count_spec = spec
    return self
def depth(self, levels: np.ndarray, depth: float, xycoords=None) -> np.ndarray:
        D = max(float(depth), 1e-9)
        Nminus1 = max(1, int(math.ceil(D / self.dz_max)))
        s = (np.asarray(levels, dtype=float) / Nminus1)
        return D * np.power(np.clip(s, 0.0, 1.0), self.p)


class SigmaTwoZoneSizeFunction(VerticalSizeFunction):
    """Reserve K bottom interfaces at dz_bottom (if possible), fill the rest linearly."""
    def __init__(self, dz_max: float = 1.5, dz_bottom: float = 0.8, K_bottom: int | None = None, frac_bottom: float = 0.25):
        self.dz_max = float(dz_max)
        self.dz_bottom = float(dz_bottom)
        self.K_bottom = K_bottom
        self.frac_bottom = float(frac_bottom)

    def depth(self, levels: np.ndarray, depth: float, xycoords=None) -> np.ndarray:
        D = max(float(depth), 1e-9)
        Nminus1 = max(1, int(math.ceil(D / self.dz_max)))
        K = int(self.K_bottom) if self.K_bottom is not None else max(2, int(self.frac_bottom * Nminus1))
        K = min(K, Nminus1)  # guard
        # bottom block occupies thickness B = min(K*dz_bottom, D)
        B = min(self.dz_bottom * K, D)
        # remaining interfaces above: M = Nminus1 - K (could be 0)
        M = max(0, Nminus1 - K)
        lev = np.asarray(levels, dtype=float)
        z = np.empty_like(lev, dtype=float)
        for i, ell in enumerate(lev):
            if ell <= M:
                # linear σ over the top segment of thickness (D-B) with M steps
                z[i] = 0.0 if M == 0 else (D - B) * (ell / M)
            else:
                j = ell - M  # 1..K bottom steps
                z[i] = (D - B) + min(j * self.dz_bottom, B)
        return z


class SigmaSBlendSizeFunction(VerticalSizeFunction, _ShapeAwareCountMixin):
    """
    σ→S blend (depth-from-surface return):
      - Choose N so nominal dz <= dz_max.
      - Build uniform σ levels s∈[0,1] (0=surface).
      - Compute S-stretch C(s; θb, θf) (normalized to [0,1]).
      - Blend: ψ(s) = (1−w)*s + w*S_norm(s), w = smoothstep((D−Hc)/W) clamped to [0,1], scaled by alpha.
      - Return depths = D * ψ(s), strictly monotone in [0,D].
    """
    def __init__(self, dz_max=1.5, theta_b=0.4, theta_f=3.0, alpha=1.0, Hc=6.0, W=4.0):
        self.dz_max  = float(dz_max)
        self.theta_b = float(theta_b)
        self.theta_f = float(theta_f)
        self.alpha   = float(alpha)
        self.Hc      = float(Hc)
        self.W       = float(max(W, 1e-3))

    @staticmethod
    def _C_raw(sigma, theta_b, theta_f):
        # sigma in [-1..0] (0=surface). Standard ROMS/SZ form.
        import numpy as _np
        tb, tf = float(theta_b), float(theta_f)
        return (1.0 - tb) * _np.sinh(tf * sigma) / (_np.sinh(tf) + 1e-12) + \
               tb * ( _np.tanh(tf * (sigma + 0.5)) - _np.tanh(tf * 0.5) ) * 0.5 / (_np.tanh(tf * 0.5) + 1e-12)

    def _C_norm(self, s01):
        import numpy as _np
        sigma = -_np.asarray(s01, dtype=float)  # [0..1] -> [0..-1]
        c = self._C_raw(sigma, self.theta_b, self.theta_f)
        c0 = self._C_raw(0.0, self.theta_b, self.theta_f)
        c1 = self._C_raw(-1.0, self.theta_b, self.theta_f)
        return (c - c0) / ((c1 - c0) + 1e-12)

    def _w(self, D):
        # smooth 0→1 around Hc with half-width W, then scale by alpha
        x = max(-1.0, min(1.0, (float(D) - self.Hc) / self.W))
        t = 0.5 * (x + 1.0)            # map [-1,1]→[0,1]
        w = t * t * (3 - 2 * t)        # smoothstep
        return max(0.0, min(1.0, self.alpha * w))

    def _psi_with_depth(self, s, D):
        import numpy as _np
        s = _np.clip(_np.asarray(s, float), 0.0, 1.0)
        s_ss = _np.clip(self._C_norm(s), 0.0, 1.0)
        w = self._w(D)
        return (1.0 - w) * s + w * s_ss

    def required_intervals(self, D: float) -> int:
        """
        Top-anchored count: choose smallest nint s.t. dz_top <= tau,
        where tau is specified (once) via count_spec = {'anchor':'top','value': tau}.
        """
        import numpy as _np, math
        if not _np.isfinite(D) or D <= 0.0:
            return 1
        spec = getattr(self, "count_spec", None) or {}
        if str(spec.get("anchor", "top")).lower() != "top":
            raise ValueError("S-blend requires count_spec={'anchor':'top','value': <meters>}")

        tau = float(spec.get("value", 0.75))  # desired top thickness
        tau = max(tau, 1e-6)

        # continuous slope of ψ at the surface (depth-aware)
        eps = 1e-6
        dpsi0 = float((self._psi_with_depth(eps, D) - self._psi_with_depth(0.0, D)) / eps)
        # top layer thickness ≈ D * dψ(0;D) / nint  ⇒ nint ≥ D*dψ(0;D)/tau
        kappa = (D * max(dpsi0, 1e-9)) / tau
        return max(1, int(math.ceil(kappa)))
        
    # shape-aware hooks
    def _psi(self, s):
        import numpy as _np
        s = _np.clip(_np.asarray(s,float),0.0,1.0); s_ss = _np.clip(self._C_norm(s),0.0,1.0); w = self._w(1.0); return (1.0-w)*s + w*s_ss

    def set_count_spec(self, spec):
        self.count_spec = spec
        return self

    def depth(self, levels, depth, xycoords=None):
        """
        Return cumulative depth H(levels; D) for the S-blend mapping.
        IMPORTANT: interval count N is taken from the provided `levels`
        (i.e., from the upstream Nlevels decision), not from dz_max.
        """
        import numpy as _np

        lev = _np.asarray(levels, dtype=float).reshape(-1)
        D = float(depth)

        # Trivial guards
        if lev.size == 0:
            return _np.array([], dtype=float)
        if not _np.isfinite(D) or D <= 0.0:
            return _np.zeros_like(lev, dtype=float)

        # Derive number of intervals from levels (robust to tiny float drift)
        lmin = _np.nanmin(lev)
        lmax = _np.nanmax(lev)
        nint = int(_np.rint(lmax - lmin))
        if nint < 1:
            nint = 1

        # Normalize to s in [0,1]
        s = _np.clip((lev - lmin) / nint, 0.0, 1.0)

        # Depth-aware S-blend ψ(s; D)
        # (identical to your current formula; you may also call self._psi_with_depth(s, D) if you added it)
        psi = self._psi_with_depth(s, D)

        # Enforce exact endpoints at the true min/max level positions
        i_min = int(_np.nanargmin(lev))
        i_max = int(_np.nanargmax(lev))
        psi[i_min] = 0.0
        psi[i_max] = 1.0

        # Monotonicity (levels are expected ascending; this keeps tiny numerical dips out)
        psi = _np.maximum.accumulate(psi)

        return D * psi


class BilinearDensitySizeFunction(object):

    def __init__(self):
        # self.params = np.array([0.6,0.066,0.002,1e-5])
        #self.params = np.array([0.5, 0.012, 0.1, 1e-5])
        #self.params = np.array([0.5, 0.005, 0.09, -8.5e-5])
        self.params = np.array([0.5, 0.005, 0.08, -3.0e-5])


    def density(self, z, h, x):
        (a, b, c, d) = self.params
        return a + b * h + c * z + d * z * h

    def depth(self, levels, depth, xycoords=None):
        # Raw bilinear cumulative depth (meters) at the provided integer levels
        t = np.asarray(levels, dtype=float)
        D = float(depth)
        (a, b, c, d) = [float(v) for v in self.params]
        m = a + b * D
        n = c + d * D
        # Guard near-linear regime
        if abs(n) < 1e-12:
            H_raw = (m * t)
        else:
            H_raw = (np.exp(n * t) - 1.0) * (m / n)
        # Ensure the last provided level lands exactly at the true depth D
        if H_raw.size:
            denom = float(H_raw[-1]) if np.isfinite(H_raw[-1]) else 0.0
            scale = (D / denom) if denom > 0.0 else (D / max(m, 1e-12))
            H = H_raw * scale
        else:
            H = H_raw
        return np.asarray(H).squeeze()

    def required_intervals(self, D: float) -> int:
        """
        Minimal nint such that H(nint; D) >= D for bilinear mapping
        H(ell; D) = (m/n) * (exp(n*ell) - 1),
        m = a + b*D, n = c + d*D.
        """
        import math, numpy as np
        if not np.isfinite(D) or D <= 0.0:
            return 1
        a, b, c, d = map(float, self._impl.params if hasattr(self, "_impl") else self.params)
        m = a + b * D
        n = c + d * D
        # Guard tiny |n| (treat as linear limit)
        if abs(n) < 1e-12:
            # H(ell) ~ m * ell; need m*ell >= D
            ell = D / max(m, 1e-12)
        else:
            ell = (1.0 / n) * math.log(1.0 + (D * n) / max(m, 1e-12))
        return max(1, int(math.ceil(ell)))
        
    def continuous_index_for_depth(self, D: float) -> float:
        import math
        # get coefficients whichever slot you store them in
        params = getattr(self, "params", None)
        if params is None and hasattr(self, "_impl"):
            params = getattr(self._impl, "params", None)
        a, b, c, d = [float(x) for x in params]
        m = a + b * D
        n = c + d * D
        if abs(n) < 1e-12:
            return D / max(m, 1e-12)   # linear limit
        return (1.0 / n) * math.log(1.0 + (D * n) / max(m, 1e-12))        
############################




DEBUG_LSC2 = False  # toggle at runtime as needed

def _first_new_nan(name, prev, cur, mask=None):
    """Raise on first entry that became NaN/Inf between prev and cur."""
    if not DEBUG_LSC2: 
        return
    m = ~np.isfinite(prev) & ~np.isfinite(cur)  # already NaN before, ignore
    became_bad = (~np.isfinite(cur)) & (np.isfinite(prev))
    if mask is not None:
        became_bad &= mask
    if became_bad.any():
        i, k = np.argwhere(became_bad)[0]
        raise RuntimeError(f"[{name}] first new NaN/Inf at node {i}, level {k}")

def _dbg(msg, *args):
    if DEBUG_LSC2:
        logger.debug(msg.format(*args))

def _assert_finite(name, arr, where_mask=None):
    if not DEBUG_LSC2: 
        return
    a = arr if where_mask is None else arr[where_mask]
    if a.size == 0:
        return
    bad = ~np.isfinite(a)
    if np.any(bad):
        # report a few examples to avoid spam
        idx = np.flatnonzero(bad)[:10]
        raise RuntimeError(f"[NaN/Inf] {name}: count={bad.sum()}, examples at {idx}")

def _assert_nonzero(name, arr, where_mask=None):
    if not DEBUG_LSC2:
        return
    a = arr if where_mask is None else arr[where_mask]
    if a.size == 0:
        return
    zeros = (a == 0)
    if np.any(zeros):
        idx = np.flatnonzero(zeros)[:10]
        raise RuntimeError(f"[Zero] {name}: zeros={zeros.sum()}, examples at {idx}")

def _summary_mask(name, m):
    if DEBUG_LSC2:
        _dbg("{}: true={}, false={}", name, int(m.sum()), int((~m).sum()))

# ---------------------------------------------------------------------
# Post-process σ: remove duplicates and fix bottom slivers
# ---------------------------------------------------------------------
def fix_sigma_pileups(
    sigma: np.ndarray,
    Nlevels: np.ndarray,
    depth: np.ndarray,
    tmin,
    mesh=None,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Inputs:
      - sigma: left-justified σ in [0,1], padded with NaN on the right
      - Nlevels: per-node level counts (includes bed & surface)
      - depth: total depth D (m). Only used to evaluate bottom thickness.
      - tmin: minimum bottom-cell thickness — scalar float or per-node array
      - mesh: optional, for x/y coordinates logging (mesh.nodes[:,0:2])
    Always checks/fixes; returns (sigma_fixed, Nlevels_fixed, log_df).
    Sweep rule: try to make the *bottom* cell ≥ tmin; if blocked by the
    next interface, drop that interface and try again; stop at first success.
    """
    # Accept scalar or per-node array
    tmin_arr = np.broadcast_to(np.asarray(tmin, dtype=float),
                               (sigma.shape[0],)).copy()
    TOL = 1.0e-10  # internal; not exposed
    n, maxN = sigma.shape
    S = sigma.copy()
    NL = Nlevels.copy()
    logs = []

    have_xy = hasattr(mesh, "nodes") and hasattr(mesh.nodes, "__getitem__")

    for i in range(n):
        Ni0 = int(NL[i])
        assert Ni0 >= 2
        s = S[i, :Ni0]
        assert np.all(np.isfinite(s))

        D = float(depth[i])
        tmin_i = float(tmin_arr[i])
        x = y = np.nan
        if have_xy:
            try:
                x, y = float(mesh.nodes[i, 0]), float(mesh.nodes[i, 1])
            except Exception:
                pass

        changed = False

        # 1) Remove duplicates anywhere (keep first of a run)
        if Ni0 >= 2:
            d = np.abs(np.diff(s))
            dup_mask = np.r_[True, d > TOL]  # True keeps; False removes
            dupes_removed = int(Ni0 - int(dup_mask.sum()))
            if dupes_removed:
                s = s[dup_mask]
                changed = True
        else:
            dupes_removed = 0

        Ni = int(s.size)
        if Ni <= 2:
            Ni = 2
            # degenerate; pad back to at least [0,1] if possible
            NL[i] = 2
            S[i, :] = np.nan
            S[i, :Ni] = [0., 1.]
            continue

        # Ensure exact surface at 1.0 if within tolerance
        if abs(s[-1] - 1.0) <= TOL:
            s[-1] = 1.0

        # 2a) Surface sliver fix: ensure the first cell (near σ=0, physical
        #     surface) is at least tmin thick.
        bottom_adjusted = False
        bottom_drops = 0
        bed_adjusted = False
        bed_drops = 0
        dz1_after = np.nan
        dz2_after = np.nan

        if D > 0.0 and Ni >= 3 and tmin_i > 0.0:
            # Keep trying until first success or run out of interfaces
            while True:
                # current surface-cell thickness
                dz = D * (s[1] - s[0])
                if dz + 1e-15 >= tmin_i:
                    bottom_adjusted = changed
                    break
                # try to raise s[1] to target
                target = s[0] + (tmin_i / max(D, 1e-12))
                if target < s[2] - TOL:
                    s[1] = target
                    changed = True
                    bottom_adjusted = True
                    break  # first success
                # blocked by s[2] → drop s[1] and retry
                s = np.delete(s, 1)
                Ni -= 1
                bottom_drops += 1
                changed = True
                if Ni < 3:
                    break

        # 2b) Bed sliver fix: ensure the last cell (near σ=1, physical
        #     bottom / bed) is at least tmin thick.
        if D > 0.0 and Ni >= 3 and tmin_i > 0.0:
            while True:
                dz = D * (s[-1] - s[-2])
                if dz + 1e-15 >= tmin_i:
                    bed_adjusted = changed
                    break
                # try to lower s[-2] to make room
                target = s[-1] - (tmin_i / max(D, 1e-12))
                if target > s[-3] + TOL:
                    s[-2] = target
                    changed = True
                    bed_adjusted = True
                    break
                # blocked by s[-3] → drop s[-2] and retry
                s = np.delete(s, -2)
                Ni -= 1
                bed_drops += 1
                changed = True
                if Ni < 3:
                    break

        # recompute dz1/dz2 after changes
        if D > 0.0 and Ni >= 3:
            dz1_after = D * (s[1] - s[0])
            dz2_after = D * (s[-1] - s[-2])
        elif D > 0.0 and Ni == 2:
            dz1_after = D  # single layer (bed→surface)

        if changed:
            # write back, pad with NaN to the right
            S[i, :] = np.nan
            S[i, :Ni] = s
            NL[i] = Ni
            logs.append({
                "node": i,
                "x": x, "y": y,
                "Ni_before": Ni0,
                "Ni_after": Ni,
                "dupes_removed": dupes_removed,
                "surface_drops": bottom_drops,
                "bed_drops": bed_drops,
                "dz1_after_m": dz1_after,
                "dz2_after_m": dz2_after,
                "depth_m": D,
            })

    # ------------------------------------------------------------------
    # Final pass: snap sigma to the output precision used by the vgrid
    # writer (6 decimal places, i.e. %14.6f).  Without this, two sigma
    # values that differ by less than 5e-7 survive the TOL=1e-10 check
    # above but become identical after rounding in vgrid.in, which causes
    # SCHISM to abort with "Weird side (3)".
    # ------------------------------------------------------------------
    OUTPUT_DECIMALS = 6
    for i in range(n):
        Ni = int(NL[i])
        if Ni < 2:
            continue
        s = np.round(S[i, :Ni], OUTPUT_DECIMALS)
        # Ensure endpoints are exact
        s[0] = 0.0
        s[-1] = 1.0
        # Remove consecutive duplicates introduced by rounding
        keep = np.r_[True, np.diff(s) > 0.0]
        # Always keep the last entry (surface = 1.0)
        keep[-1] = True
        if not keep.all():
            s = s[keep]
            S[i, :] = np.nan
            S[i, :len(s)] = s
            NL[i] = len(s)
            logs.append({
                "node": i,
                "x": float(mesh.nodes[i, 0]) if have_xy else np.nan,
                "y": float(mesh.nodes[i, 1]) if have_xy else np.nan,
                "Ni_before": Ni,
                "Ni_after": len(s),
                "dupes_removed": int(Ni - len(s)),
                "bottom_drops": 0,
                "dz1_after_m": np.nan,
                "dz2_after_m": np.nan,
                "depth_m": float(depth[i]),
            })

    df = pd.DataFrame(logs, columns=[
        "node","x","y","Ni_before","Ni_after",
        "dupes_removed","bottom_drops","dz1_after_m","dz2_after_m","depth_m"
    ])
    return S, NL, df


#############################

def _build_row_stochastic_W(mesh, n_nodes):
    """
    Build a CSR matrix W with row-stochastic averaging over *neighbors only*.
    Semantics:
      - (W @ x)[i] = mean( x[j] for j in N(i) ), i.e., excludes self i.
      - Each nonempty row is normalized to sum to 1. Empty rows are left as 0
        (the caller must decide how to treat nodes with no neighbors).
    Intended use:
      - Mid-column smoothing (and the one-pass interface smoother) where we
        combine neighbor information uniformly and, when a column lacks a
        given level k, we explicitly define a *ghost value* for that column
        (e.g., bed depth) **before** multiplying by W.
      - In other words: “neighbors-only + explicit ghost fill by the caller”.
    """
    indptr = [0]
    indices = []
    data = []
    for i in range(n_nodes):
        nbrs = mesh.get_neighbor_nodes(i)  # assume returns a list of neighbor node indices
        if not nbrs:
            # keep row empty; we'll treat row-sum as 1 during normalization
            indptr.append(len(indices))
            continue
        indices.extend(nbrs)
        data.extend([1.0] * len(nbrs))
        indptr.append(len(indices))

    W = sp.csr_matrix((np.asarray(data, float),
                       np.asarray(indices, int),
                       np.asarray(indptr, int)),
                      shape=(n_nodes, n_nodes))
    # Normalize rows to sum 1 (guard empty rows)
    row_sums = np.array(W.sum(axis=1)).ravel()
    row_sums[row_sums == 0.0] = 1.0
    Dinv = sp.diags(1.0 / row_sums)
    return Dinv @ W

def smooth_interfaces_onepass_with_bed_fill_sparse(mesh,
                                                   h: np.ndarray,
                                                   Nlevels: np.ndarray,
                                                   depth: np.ndarray,
                                                   alpha_top: float = 0.8,
                                                   alpha_bot: float = 0.1,
                                                   power: float = 1.0) -> None:
    """
    One Jacobi pass of interface-level smoothing across neighbors using a
    row-stochastic CSR averaging matrix (SciPy).

    Ghost-fill rule:
      If neighbor j lacks level k (k >= Nlevels[j]), use depth[j] (bed) as its value.
    Only nodes that *have* level k (k < Nlevels[i]) are updated. In-place update of h.
    """
    n, maxN = h.shape
    if n == 0 or maxN == 0:
        return

    W = _get_row_stochastic_W(mesh, n)
    hold = h.copy()  # Jacobi: read old, write new

    # Precompute safe denominators for sigma fraction
    Ni = np.maximum(Nlevels.astype(float), 2.0)

    for k in range(maxN):
        have_k = (Nlevels > k)
        if not np.any(have_k):
            continue

        # Nodes where this k is the bottom interface for that column
        bottom_mask = have_k & ( (Nlevels - 1) == k )
        # We update only non-bottom interfaces
        upd = have_k & (~bottom_mask)

        # Ghost fill: use bed where column lacks this level
        Hvirt = np.where(have_k, hold[:, k], depth)

        # Neighbor mean via sparse matvec
        nbr_mean = W @ Hvirt  # shape (n,)

        # Depth-weighted alpha (small/top .. larger/mid .. zero at bottom)
        sigma_frac = np.zeros(n, dtype=float)
        Ni = np.maximum(Nlevels.astype(float), 2.0)
        sigma_frac[have_k] = k / (Ni[have_k] - 1.0)
        alpha = alpha_top + (alpha_bot - alpha_top) * (sigma_frac ** power)

        # Do not move the bed interface
        alpha[bottom_mask] = 0.0

        # Jacobi blend; update only nodes that actually have this level and are not bottom
        h_col_new = (1.0 - alpha) * hold[:, k] + alpha * nbr_mean
        # Clamp: h must not exceed local depth (physically below bed)
        h_col_new = np.minimum(h_col_new, depth)
        h[upd, k] = h_col_new[upd]

def _get_row_stochastic_W(mesh, n_nodes):
    """
    Return a cached row-stochastic neighbor-averaging matrix (see
    _build_row_stochastic_W). W is used where we want:
      • neighbors-only means,
      • caller-controlled ghost fill (e.g., use bed for missing level k),
      • and no implicit inclusion of self in the average.
    """
    key = "_row_stochastic_W"
    W = getattr(mesh, key, None)
    if W is None or W.shape != (n_nodes, n_nodes):
        W = _build_row_stochastic_W(mesh, n_nodes)
        setattr(mesh, key, W)
    return W



# -----------------------------------------------------------------------------
# A) Continuous desired levels L* and light smoothing
# -----------------------------------------------------------------------------

def compute_Lstar(depth: np.ndarray,
                  sizefun: VerticalSizeFunction,
                  max_levels_cap: int = 80,
                  mesh=None) -> np.ndarray:
    """Compute real-valued desired number of levels L*(x) from the size function.

    Strategy: for each node i with depth D_i, evaluate cumulative depth H(ell; D_i)
    on a dense grid of ell (0..Lcap-1), then invert by linear interpolation
    to find real ell such that H(ell; D_i) == D_i.

    Returns
    -------
    Lstar : np.ndarray (float)
        Real-valued desired count of *levels* at each node (>=2).
    """
    n = depth.shape[0]
    # Clamp trivial shallows
    D = np.maximum(depth, 1e-6) # todo: this assumes that negative depth is invalid

    # Dense level grid to invert on (0..cap-1)
    ell = np.arange(max_levels_cap, dtype=float)
    Lstar = np.empty(n, dtype=float)

    # Vectorized evaluation by batching nodes with similar depths would be ideal.
    # For simplicity & clarity, loop; this is called once and is cheap relative to IO.
    for i in range(n):
        Di = float(D[i])

        # Prefer the size-function’s own count (e.g. S-blend with top anchor)
        if hasattr(sizefun, "required_intervals"):
            nint = int(sizefun.required_intervals(Di))
            # real-valued L* is “levels”, i.e., intervals+1
            Lstar[i] = max(2.0, float(nint + 1))
            continue

        # Fallback: old inversion on a dense grid
        H = np.asarray(sizefun.depth(ell, Di, 0.0)).squeeze()
        H = np.maximum.accumulate(H)
        if H[-1] < Di:
            dH = H[-1] - H[-2] if len(H) > 1 else 1.0
            Lstar[i] = (max_levels_cap - 1) + max(0.0, (Di - H[-1]) / max(dH, 1e-6))
        else:
            j = np.searchsorted(H, Di)
            if j == 0:
                Lstar[i] = 2.0
            else:
                Lstar[i] = (j - 1) + (Di - H[j-1]) / max(H[j] - H[j-1], 1e-9)

        if Lstar[i] < 2.0:
            Lstar[i] = 2.0


    return Lstar


def smooth_scalar_on_mesh(mesh, field: np.ndarray,
                          passes: int = 2,
                          kappa: float = 0.5,
                          dt: float = 1.0) -> np.ndarray:
    """Light isotropic Laplacian smoothing on a scalar field defined on mesh nodes.

    This is the legacy LSC2-v2 smoother.  It is intentionally retained so
    existing control files reproduce previous behavior unless the caller
    explicitly selects an edge-aware method.
    """
    neighbors = [np.array(mesh.get_neighbor_nodes(i), dtype=np.int32)
                 for i in range(mesh.nodes.shape[0])]
    return smooth_kernel_numba(field.astype(np.float64), neighbors, kappa, dt, passes)


def _smooth_scalar_on_mesh_bilateral(mesh,
                                     field: np.ndarray,
                                     depth: np.ndarray,
                                     passes: int = 2,
                                     kappa: float = 0.5,
                                     dt: float = 1.0,
                                     L_scale: float = 1.25,
                                     depth_scale: float = 4.0,
                                     length_power: float = 0.0) -> np.ndarray:
    """Edge-aware smoothing for scalar node fields such as L*.

    The update is a graph diffusion, but each edge conductance is reduced when
    neighboring nodes have a large contrast in the unsmoothed scalar field and/or
    bathymetric depth:

        w_ij = ell_ij**(-length_power)
               exp(-(dL/L_scale)**2)
               exp(-(dH/depth_scale)**2)

    Using the *initial* field for dL makes this a bilateral filter: narrow
    channels, pits, and submerged thalwegs can smooth along coherent features
    without being diffused as strongly across their margins.  Set either scale
    to <=0 to disable that factor.  length_power defaults to zero so the first
    attempt changes only anisotropy, not the historical neighbor averaging
    strength as a function of local mesh size.
    """
    x = np.asarray(field, dtype=np.float64).copy()
    ref = x.copy()
    H = np.asarray(depth, dtype=np.float64)
    n = x.size

    if passes <= 0 or kappa == 0.0:
        return x

    coords = np.asarray(mesh.nodes[:, 0:2], dtype=np.float64)
    use_L = L_scale is not None and float(L_scale) > 0.0
    use_H = depth_scale is not None and float(depth_scale) > 0.0
    lp = float(length_power or 0.0)

    neighbors = [np.asarray(mesh.get_neighbor_nodes(i), dtype=np.int32)
                 for i in range(n)]
    weights = []
    for i, nbrs in enumerate(neighbors):
        if nbrs.size == 0:
            weights.append(np.zeros(0, dtype=np.float64))
            continue
        w = np.ones(nbrs.size, dtype=np.float64)
        if lp != 0.0:
            dxy = coords[nbrs] - coords[i]
            ell = np.sqrt(np.sum(dxy * dxy, axis=1))
            w *= np.power(np.maximum(ell, 1.0e-12), -lp)
        if use_L:
            dL = ref[nbrs] - ref[i]
            w *= np.exp(-((dL / float(L_scale)) ** 2))
        if use_H:
            dH = H[nbrs] - H[i]
            w *= np.exp(-((dH / float(depth_scale)) ** 2))
        weights.append(w)

    # Jacobi iterations with row-normalized weighted neighbor mean.
    # Stability is comparable to the legacy smoother for kappa*dt <= 1.
    alpha = float(kappa) * float(dt)
    for _ in range(int(passes)):
        old = x.copy()
        for i, nbrs in enumerate(neighbors):
            if nbrs.size == 0:
                continue
            w = weights[i]
            sw = float(np.sum(w))
            if sw <= 0.0 or not np.isfinite(sw):
                continue
            nbr_mean = float(np.dot(w, old[nbrs]) / sw)
            x[i] = old[i] + alpha * (nbr_mean - old[i])
    return x


def smooth_Lstar_on_mesh(mesh,
                         Lstar: np.ndarray,
                         depth: np.ndarray,
                         method: str = "isotropic",
                         passes: int = 2,
                         kappa: float = 0.5,
                         dt: float = 1.0,
                         L_scale: float = 1.25,
                         depth_scale: float = 4.0,
                         length_power: float = 0.0) -> np.ndarray:
    """Smooth L* using either the legacy isotropic or edge-aware method."""
    method = str(method or "isotropic").strip().lower()
    if method in ("isotropic", "legacy", "laplace", "laplacian"):
        return smooth_scalar_on_mesh(mesh, Lstar, passes=passes, kappa=kappa, dt=dt)
    if method in ("bilateral", "edge_aware", "edge-aware", "anisotropic"):
        return _smooth_scalar_on_mesh_bilateral(
            mesh, Lstar, depth, passes=passes, kappa=kappa, dt=dt,
            L_scale=L_scale, depth_scale=depth_scale,
            length_power=length_power,
        )
    raise ValueError(
        f"Unknown Lsmooth_method='{method}'. Expected isotropic or bilateral."
    )


def _constraint_diffusion_weights(mesh,
                                  reference: np.ndarray,
                                  depth: np.ndarray,
                                  method: str = "isotropic",
                                  L_scale: float = 1.25,
                                  depth_scale: float = 4.0,
                                  length_power: float = 0.0):
    """Return neighbor lists and row-normalized weights for L-constraint diffusion.

    The weight semantics intentionally match ``smooth_Lstar_on_mesh``:
    ``isotropic`` gives an unweighted neighbor mean, while ``bilateral`` reduces
    conductance across large L/depth jumps and can optionally include edge length.
    """
    method = str(method or "isotropic").strip().lower()
    n = reference.size
    coords = np.asarray(mesh.nodes[:, 0:2], dtype=np.float64)
    ref = np.asarray(reference, dtype=np.float64)
    H = np.asarray(depth, dtype=np.float64)
    use_bilateral = method in ("bilateral", "edge_aware", "edge-aware", "anisotropic")
    if method not in ("isotropic", "legacy", "laplace", "laplacian",
                      "bilateral", "edge_aware", "edge-aware", "anisotropic"):
        raise ValueError(
            f"Unknown constraint_diffusion.method='{method}'. Expected isotropic or bilateral."
        )

    use_L = use_bilateral and L_scale is not None and float(L_scale) > 0.0
    use_H = use_bilateral and depth_scale is not None and float(depth_scale) > 0.0
    lp = float(length_power or 0.0) if use_bilateral else 0.0

    neighbors = []
    weights = []
    for i in range(n):
        nbrs = np.asarray(mesh.get_neighbor_nodes(i), dtype=np.int32)
        neighbors.append(nbrs)
        if nbrs.size == 0:
            weights.append(np.zeros(0, dtype=np.float64))
            continue
        w = np.ones(nbrs.size, dtype=np.float64)
        if lp != 0.0:
            dxy = coords[nbrs] - coords[i]
            ell = np.sqrt(np.sum(dxy * dxy, axis=1))
            w *= np.power(np.maximum(ell, 1.0e-12), -lp)
        if use_L:
            dL = ref[nbrs] - ref[i]
            w *= np.exp(-((dL / float(L_scale)) ** 2))
        if use_H:
            dH = H[nbrs] - H[i]
            w *= np.exp(-((dH / float(depth_scale)) ** 2))
        sw = float(np.sum(w))
        if sw > 0.0 and np.isfinite(sw):
            w /= sw
        weights.append(w)
    return neighbors, weights


def diffuse_Ltilde_with_constraints(mesh,
                                    Ltilde: np.ndarray,
                                    depth: np.ndarray,
                                    n_min: Optional[np.ndarray] = None,
                                    n_max: Optional[np.ndarray] = None,
                                    passes: int = 40,
                                    weight: float = 0.6,
                                    method: str = "isotropic",
                                    L_scale: float = 1.25,
                                    depth_scale: float = 4.0,
                                    length_power: float = 0.0,
                                    tol: float = 0.0) -> np.ndarray:
    """Diffuse continuous level counts while respecting polygon constraints.

    Fixed constraints (``n_min == n_max``) are treated as Dirichlet nodes and
    are re-imposed every iteration.  Range/min/max constraints are box
    constraints: after each screened diffusion update, constrained nodes are
    projected back into their allowed interval.  Unconstrained nodes are free
    except for the global two-level minimum.

    This stage is intended to run after the ordinary L* smoother.  Unlike the
    legacy inside-polygon clip, fixed constrained regions continue to act as
    sources during smoothing, so their influence can propagate into adjacent
    range polygons and unconstrained transition zones.
    """
    x0 = np.asarray(Ltilde, dtype=np.float64)
    x = x0.copy()
    n = x.size
    if passes <= 0 or weight <= 0.0 or (n_min is None and n_max is None):
        return x

    lo = np.full(n, 2.0, dtype=np.float64)
    hi = np.full(n, np.inf, dtype=np.float64)
    if n_min is not None:
        lo = np.maximum(lo, np.asarray(n_min, dtype=np.float64))
    if n_max is not None:
        hi = np.minimum(hi, np.asarray(n_max, dtype=np.float64))

    bad = lo > hi
    if np.any(bad):
        i = int(np.flatnonzero(bad)[0])
        raise ValueError(
            f"constraint diffusion got inconsistent bounds at node {i}: lo={lo[i]} hi={hi[i]}"
        )

    constrained = np.isfinite(hi) | (lo > 2.0)
    fixed = constrained & np.isfinite(hi) & np.isclose(lo, hi)

    # Project once before iteration so fixed nodes are already sources.
    x = np.minimum(np.maximum(x, lo), hi)
    x[fixed] = lo[fixed]

    alpha = min(max(float(weight), 0.0), 1.0)
    neighbors, weights = _constraint_diffusion_weights(
        mesh, x0, depth, method=method,
        L_scale=L_scale, depth_scale=depth_scale, length_power=length_power,
    )

    last_delta = 0.0
    for _ in range(int(passes)):
        old = x.copy()
        for i, nbrs in enumerate(neighbors):
            if fixed[i] or nbrs.size == 0:
                continue
            w = weights[i]
            if w.size == 0:
                continue
            nbr_mean = float(np.dot(w, old[nbrs]))
            # Screened diffusion: retain a pull to the pre-constraint smoothed field.
            x[i] = old[i] + alpha * (nbr_mean - old[i])

        x = np.minimum(np.maximum(x, lo), hi)
        x[fixed] = lo[fixed]
        last_delta = float(np.max(np.abs(x - old)))
        if tol and last_delta <= float(tol):
            break

    logger.info(
        "   constraint diffusion: passes=%d weight=%.3g method=%s constrained=%d fixed=%d max_delta=%.3g",
        int(passes), alpha, method, int(np.count_nonzero(constrained)),
        int(np.count_nonzero(fixed)), last_delta,
    )
    return x


def cleanup_integer_components(mesh,
                               Nlevels: np.ndarray,
                               max_component_nodes: int = 6,
                               max_passes: int = 5,
                               n_min: Optional[np.ndarray] = None,
                               n_max: Optional[np.ndarray] = None) -> np.ndarray:
    """Remove small strict-extremum components in integer level counts.

    This is a post-quantization topological simplification of the integer
    ``Nlevels`` field.  It repeatedly finds exact-value connected components
    whose node count is small, but removes them only when they are strict local
    extrema relative to their exterior boundary:

      * hill: component value > every boundary value
      * dip:  component value < every boundary value

    Small components that form a legitimate transition band are preserved.  For
    example, a small ``22`` component between exterior ``21`` and ``23`` values
    is not removed.  This avoids eroding stepped ramps into cliffs.

    The operation is scale-free: "small" is measured in mesh-node count, relying
    on the mesh adaptation to encode physical scale.  It is also constraint-aware:
    if polygon n_min/n_max arrays are present, candidate replacement values that
    would violate any node in the component are ignored.
    """
    N = np.asarray(Nlevels, dtype=np.int32).copy()
    max_nodes = int(max_component_nodes or 0)
    max_passes = int(max_passes or 0)
    if max_nodes <= 0 or max_passes <= 0:
        return N

    n = N.size
    neighbors = [np.asarray(mesh.get_neighbor_nodes(i), dtype=np.int32)
                 for i in range(n)]

    if n_min is None:
        nmin = np.full(n, -2**30, dtype=np.int32)
    else:
        nmin = np.asarray(n_min, dtype=np.int32)
    if n_max is None:
        nmax = np.full(n, 2**30, dtype=np.int32)
    else:
        nmax = np.asarray(n_max, dtype=np.int32)

    total_components = 0
    total_nodes = 0
    total_hills = 0
    total_dips = 0

    for ipass in range(max_passes):
        changed_components = 0
        changed_nodes = 0
        changed_hills = 0
        changed_dips = 0

        # Work high-to-low so unsupported peaks/ridges peel off naturally.
        # Dips are still removed in the same pass when their value is reached.
        values = np.unique(N).astype(np.int32)[::-1]

        for value in values:
            is_value = (N == value)
            if not np.any(is_value):
                continue
            visited = np.zeros(n, dtype=bool)
            seeds = np.flatnonzero(is_value)

            for seed in seeds:
                if visited[seed] or N[seed] != value:
                    continue

                # Flood-fill one exact-value component.
                stack = [int(seed)]
                visited[seed] = True
                comp = []
                boundary_vals = []

                while stack:
                    i = stack.pop()
                    comp.append(i)
                    for j in neighbors[i]:
                        jj = int(j)
                        if N[jj] == value:
                            if not visited[jj]:
                                visited[jj] = True
                                stack.append(jj)
                        else:
                            boundary_vals.append(int(N[jj]))

                if len(comp) > max_nodes or not boundary_vals:
                    continue

                comp_idx = np.asarray(comp, dtype=np.int32)
                bvals = np.asarray(boundary_vals, dtype=np.int32)

                # Only remove true local extrema.  If exterior values exist on
                # both sides of the component value, it is a transition band,
                # not a hill/dip artifact, and must be preserved.
                bmin = int(np.min(bvals))
                bmax = int(np.max(bvals))
                is_hill = int(value) > bmax
                is_dip = int(value) < bmin
                if not (is_hill or is_dip):
                    continue

                # Respect polygon constraints, if present, for every node in C.
                lo = int(np.max(nmin[comp_idx]))
                hi = int(np.min(nmax[comp_idx]))
                allowed = bvals[(bvals >= lo) & (bvals <= hi)]
                if allowed.size == 0:
                    continue

                # Choose boundary mode; tie-break by closeness to current value,
                # then by smaller value for deterministic behavior.
                u, counts = np.unique(allowed, return_counts=True)
                max_count = counts.max()
                tied = u[counts == max_count]
                repl = int(sorted(tied, key=lambda x: (abs(int(x) - int(value)), int(x)))[0])
                if repl == int(value):
                    continue

                N[comp_idx] = repl
                changed_components += 1
                changed_nodes += len(comp_idx)
                if is_hill:
                    changed_hills += 1
                else:
                    changed_dips += 1

        logger.info(
            "   integer cleanup: pass %d removed %d extrema components "
            "(%d nodes; hills=%d dips=%d)",
            ipass + 1, changed_components, changed_nodes,
            changed_hills, changed_dips,
        )
        total_components += changed_components
        total_nodes += changed_nodes
        total_hills += changed_hills
        total_dips += changed_dips
        if changed_components == 0:
            break

    if total_components:
        logger.info(
            "   integer cleanup: total removed %d extrema components "
            "(%d nodes; hills=%d dips=%d)",
            total_components, total_nodes, total_hills, total_dips,
        )
    return N


# -----------------------------------------------------------------------------
# B) Hysteretic quantization of levels (no regional bounds)
# -----------------------------------------------------------------------------

@dataclass
class HysteresisParams:
    delta_up: float = 0.7
    delta_down: float = 0.7
    Jmax: int = 1
    eps_depth: float = 0.75
    max_sweeps: int = 6
    do_median_cleanup: bool = True


def quantize_levels_with_hysteresis(mesh,
                                    Ltilde: np.ndarray,
                                    depth: np.ndarray,
                                    params: HysteresisParams = HysteresisParams()
                                    ) -> np.ndarray:
    """Convert real-valued desired levels to integers with hysteresis and spatial gate.

    Returns integer *levels* N (not layers).
    """
    N = np.rint(Ltilde).astype(np.int32)

    # Optional graph 2-coloring (naive: split by parity of index)
    colors = [np.where(np.arange(len(N)) % 2 == 0)[0],
              np.where(np.arange(len(N)) % 2 == 1)[0]]

    for _ in range(params.max_sweeps):
        changed = 0
        for color_nodes in colors:
            for i in color_nodes:
                k = int(N[i])
                d = float(Ltilde[i] - k)
                kprime = k
                if d > params.delta_up:
                    kprime = k + 1
                elif d < -params.delta_down:
                    kprime = k - 1

                if kprime != k:
                    ok = True
                    nbs = mesh.get_neighbor_nodes(i)
                    Di = depth[i]
                    for j in nbs:
                        if abs(kprime - int(N[j])) > params.Jmax and abs(Di - depth[j]) <= params.eps_depth:
                            ok = False
                            break
                    if ok:
                        N[i] = kprime
                        changed += 1
        if changed == 0:
            break

    if params.do_median_cleanup:
        # 1-ring median filter, then clip back toward Ltilde (floor/ceil)
        Nmed = N.copy()
        for i in range(len(N)):
            vals = [N[i]]
            vals += [N[j] for j in mesh.get_neighbor_nodes(i)]
            vals = sorted(vals)
            Nmed[i] = vals[len(vals)//2]
        N = np.minimum(np.maximum(Nmed, np.floor(Ltilde).astype(int)),
                       np.ceil(Ltilde).astype(int)).astype(np.int32)

    # Enforce minimum of 2 levels
    N[N < 2] = 2
    return N


# -----------------------------------------------------------------------------
# C) Targets for interfaces from size function
# -----------------------------------------------------------------------------

def build_target_interfaces(sizefun: VerticalSizeFunction,
                            depth: np.ndarray,
                            Nlevels: np.ndarray,
                            eta: float = 0.0,
                            mesh=None,
                            uniform_sigma_mask: Optional[np.ndarray] = None
                            ) -> Tuple[np.ndarray, np.ndarray]:
    """Compute target interface depths per node (down from surface).

    Parameters
    ----------
    uniform_sigma_mask : optional bool array (n_nodes,)
        If True for a node, use equal-spacing (uniform sigma) targets instead
        of the size-function profile.  Used for nodes where Nlevels exceeds
        the global-tmin feasibility cap.

    Returns
    -------
    hhat : (n_nodes, maxN) array with per-node valid range [0..N_i-1]
    valid_mask : boolean array same shape, True where defined
    """
    n = depth.shape[0]
    maxN = int(np.max(Nlevels))
    hhat = np.full((n, maxN), np.nan, dtype=float)
    valid = np.zeros((n, maxN), dtype=bool)

    levs_full = np.arange(maxN, dtype=float)

    for i in range(n):
        Ni = int(Nlevels[i])
        if Ni < 2:
            Ni = 2
        coords = None if mesh is None else (mesh.nodes[i,0], mesh.nodes[i,1])
        # If depth is <= 0 relative to η, treat as dry: zero interfaces, still 2 levels.
        Di = float(depth[i])
        if Di <= 0.0:
            hhat[i, :2] = (0.0, 0.0)
            valid[i, :2] = True
            continue

        # Uniform sigma for exception nodes (thin water columns with forced Nlevels)
        if uniform_sigma_mask is not None and uniform_sigma_mask[i]:
            for k in range(Ni):
                hhat[i, k] = k * Di / (Ni - 1)
            valid[i, :Ni] = True
            continue

        Hi = np.asarray(sizefun.depth(levs_full[:Ni], Di, xycoords=coords)).squeeze()

        Hi = np.maximum.accumulate(Hi)
        hhat[i, :Ni] = Hi
        valid[i, :Ni] = True

    return hhat, valid


# -----------------------------------------------------------------------------
# D) Fit interfaces with smoothing and bottom-thickness anchoring
# -----------------------------------------------------------------------------

@dataclass
class FitParams:
    # Smoothing & data-fit weights
    kappa: float = 0.35         # Laplacian weight per iteration
    dt: float = 0.5             # explicit step
    iters_bottom: int = 25      # iterations for lowest Kon
    iters_mid: int = 18         # iterations for mid column
    lam_D: float = 1.0          # pull to target hhat
    lam_B: float = 6.0          # bottom thickness anchoring weight
    K_bottom: int = 3           # number of lowest interfaces to treat as "bottom block"
    tmin: float = 0.35          # minimum layer thickness
    keep_top_thin: int = 2      # number of topmost interfaces with reduced smoothing
    top_relax: float = 0.2      # multiplier for kappa near surface
    no_collapse_ratio: float = 0.5  # don't shrink any layer below this fraction of pre-smooth
    dump_presmooth_npz: str ="pre_smooth.npz"  # Binary dump right before smoothing so smoothing can be debugged
    # --- DEBUG: interface-smoother-only mode & gains ---
    debug_smoother_only: bool = False
    sm_alpha_top: float = 0.35
    sm_alpha_bot: float = 0.02
    sm_power: float = 1.2
    sm_print_delta: bool = True
    sm_iters: int = 4


def smooth_bottom_thickness(depth: np.ndarray,
                            hhat: np.ndarray,
                            valid: np.ndarray,
                            tmin: float,
                            mesh) -> np.ndarray:
    """Build and lightly smooth the target bottom thickness field t_b(x)."""
    n, maxN = hhat.shape
    tb = np.zeros(n, dtype=float)
    for i in range(n):
        # find last two valid interfaces (N-1 and N-2)
        idx = int(np.sum(valid[i, :])) - 1
        if idx <= 0:
            tb[i] = tmin
            continue
        bottom = hhat[i, idx]
        lower = hhat[i, idx - 1]
        tb[i] = max(tmin, bottom - lower)
    tb = smooth_scalar_on_mesh(mesh, tb, passes=2, kappa=0.6, dt=1.0)
    tb = np.maximum(tb, tmin)
    return tb



def _neighbors_average(mesh, arr: np.ndarray) -> np.ndarray:
    out = arr.copy()
    for i in range(len(arr)):
        nbs = mesh.get_neighbor_nodes(i)
        if not nbs:
            continue
        out[i] = np.mean(arr[nbs])
    return out


# NOTE on A vs W:
#   A is *unnormalized* and includes self. It is used to reproduce the
#   “old bottom-block” semantics exactly:
#     vave[i] = mean( h[j,k] over {j in N(i) ∪ {i}} that actually *have* level k )
#   i.e., neighbors without level k are excluded from both sum and count.
#   We implement this by:
#     sum_k = A @ (mask_k * h[:,k])     # sum over contributors that have level k
#     cnt_k = A @ mask_k                # contributor count (no ghost fill)
#     vave  = sum_k / cnt_k             # per-node neighbor mean
#   This preserves the prior algorithm bit-for-bit and avoids any downward
#   bias that would arise from ghost-filling missing neighbors with the bed.

def _get_A_with_self(mesh, n_nodes):
    key = "_A_with_self"
    A = getattr(mesh, key, None)
    if A is not None and A.shape == (n_nodes, n_nodes):
        return A
    indptr = [0]; indices = []; data = []
    for i in range(n_nodes):
        nbrs = mesh.get_neighbor_nodes(i)
        # include self, and all neighbors
        all_idx = [i] + (nbrs or [])
        indices.extend(all_idx)
        data.extend([1.0] * len(all_idx))
        indptr.append(len(indices))
    A = sp.csr_matrix((np.asarray(data, float),
                       np.asarray(indices, int),
                       np.asarray(indptr, int)),
                      shape=(n_nodes, n_nodes))
    setattr(mesh, key, A)
    return A


def fit_interfaces(mesh,
                   depth: np.ndarray,
                   Nlevels: np.ndarray,
                   hhat: np.ndarray,
                   valid: np.ndarray,
                   tbottom: np.ndarray,
                   params: FitParams = FitParams(),
                   tmin_arr: Optional[np.ndarray] = None,
                   uniform_sigma_mask: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Fit interfaces h_k(x) to targets with smoothing, bottom anchoring, and guards.

    Parameters
    ----------
    tmin_arr : optional float array (n_nodes,)
        Per-node minimum layer thickness.  If None, uses params.tmin everywhere.
    uniform_sigma_mask : optional bool array (n_nodes,)
        If True for a node, reassert uniform sigma spacing in finalization
        (protects shallow exception nodes from smoother cliff-pressure).

    Returns
    -------
    h : (n_nodes, maxN) float array of final interfaces
    mask : boolean valid mask per node / per level
    """
    n, maxN = hhat.shape
    # Initialize h as targets (where defined) or simple sigma for safety
    h = np.where(valid, hhat, np.nan)

    # Per-node minimum thickness (scalar broadcast or caller-supplied array)
    if tmin_arr is None:
        tmin_arr = np.full(n, params.tmin, dtype=float)

    # Sequentially enforce constraints & smoothing: bottom block first
    # We'll work in-place on each interface level k
    kappa = params.kappa
    dt = params.dt
    lamD = params.lam_D
    lamB = params.lam_B
    K = params.K_bottom
    tmin = params.tmin          # global scalar kept for debug prints
    no_collapse = params.no_collapse_ratio

    # Build row-stochastic neighbor-averaging matrix once. Note that this will
    # be cached in the mesh object
    W = _get_row_stochastic_W(mesh, n)          # shape (n, n), CSR
    A = _get_A_with_self(mesh, n)   # sums over (neighbors + self)
    Ni_arr = Nlevels.astype(np.int32)

    # --- TEMP DEBUG ---
    DEBUG_NODE = 34067     # 0-based node index to trace
    DEBUG_MONO = False      # flip to False to silence prints
    # -------------------


    # Utility: enforce ordering & min thickness at a node
    def enforce_node_monotone(i, upto_k=None):
        Ni = int(Nlevels[i])
        if Ni < 2:
            return
        kk_max = Ni - 1 if upto_k is None else min(upto_k, Ni - 1)

        # --- TEMP DEBUG helper ---
        def _dbg_levels(tag):
            if DEBUG_MONO and i == DEBUG_NODE:
                arr = h[i, :Ni].copy()
                dif = np.diff(arr) if Ni >= 2 else np.array([])
                logger.debug("[monotone:%s] i=%d Ni=%d D=%.6f tmin=%.6f", tag, i, Ni, depth[i], tmin)
                logger.debug("  h  : %s", np.array2string(arr, precision=6, floatmode='fixed', max_line_width=140))
                if dif.size:
                    logger.debug("  Δh : %s  minΔ=%.6f",
                                 np.array2string(dif, precision=6, floatmode='fixed', max_line_width=140),
                                 np.nanmin(dif))
        # --- /TEMP DEBUG ---

        _dbg_levels("entry")

        # 1) Bottom-up (increasing k: surface→bed)
        tmin_i = tmin_arr[i]
        for k in range(1, kk_max + 1):
            prev = h[i, k-1]
            if not np.isfinite(prev):
                continue
            min_allowed = prev + tmin_i
            if (not np.isfinite(h[i, k])) or (h[i, k] < min_allowed):
                h[i, k] = min_allowed

        _dbg_levels("after_up")

        # 2) Pin bed exactly
        lastk = Ni - 1
        if (not np.isfinite(h[i, lastk])) or (h[i, lastk] > depth[i]):
            h[i, lastk] = depth[i]

        _dbg_levels("after_bedpin")

        # 3) Top-down (decreasing k: bed→surface)
        for k in range(lastk - 1, -1, -1):
            nextv = h[i, k+1]
            if not np.isfinite(nextv):
                continue
            max_allowed = nextv - tmin_i
            if (not np.isfinite(h[i, k])) or (h[i, k] > max_allowed):
                h[i, k] = max_allowed
            if h[i, k] < 0.0:
                h[i, k] = 0.0

        _dbg_levels("after_down")
	

    # Precompute per-node initial layer thickness for "no collapse" guard
    pre_thick = np.zeros((n, maxN), dtype=float)
    for i in range(n):
        Ni = int(Nlevels[i])
        for k in range(1, Ni):
            if valid[i, k] and valid[i, k-1]:
                pre_thick[i, k] = max(1e-6, h[i, k] - h[i, k-1])

    # -------------------------------------------------------------------------
    # Bottom block (A with self, *no* ghost fill):
    #   For each node i in the bottom band, its target level index is k_i = N[i]-block.
    #   The neighbor mean vave[i] must average *only* over contributors that
    #   actually have that level k (neighbors w/o level k do not contribute),
    #   and should *include self*. This matches the legacy “old way.”
    #   Implementation uses the unnormalized adjacency with self (A):
    #     sum_k = A @ (mask_k * h[:,k])
    #     cnt_k = A @ mask_k
    #     vave  = sum_k / cnt_k
    # -------------------------------------------------------------------------
    # Bottom block: k = N-1 down to N-K
    for i_iter in range(params.iters_bottom):
        # Iterate bottom-up within the band
        for block in range(1, params.K_bottom + 1):
            # Each node's target interface index for this block
            k_idx = Ni_arr - block
            have_k = (k_idx >= 1)                     # k must be >=1 to have a lower neighbor
            if not np.any(have_k):
                continue

            # Work by unique k to get correct neighbor means
            uniq_k = np.unique(k_idx[have_k])
            for k in uniq_k:
                # nodes using this k in this block
                S = (k_idx == k) & have_k
                if not np.any(S):
                    continue
                    
                # mask of nodes that actually have level k AND it's valid
                m_k = (Ni_arr > k) & valid[:, k]             
             
                # Debug block
                #_dbg("BOT it={} block={} k={} | S:{} m_k:{}", i_iter, block, k, S.sum(), m_k.sum())
                #_summary_mask("S mask", S)
                #_summary_mask("m_k (contributors)", m_k)
                # Debug

                # old semantics: average only over neighbors (and self) that have level k
                # no ghost-fill; exclude nodes without k from both sum and count
                hk = h[:, k]                         # (n,)
                

                
                # --- DEBUG (optional) ---
                #bad_contrib = (~np.isfinite(hk)) & m_k
                #if bad_contrib.any():
                #    _dbg("BOT k={} has {} NaN/Inf among would-be contributors", k, int(bad_contrib.sum()))
                #bad_noncontrib = (~np.isfinite(hk)) & (~m_k)
                #if bad_noncontrib.any():
                #    _dbg("BOT k={} has {} NaN/Inf outside m_k (must NOT leak into sums)", k, int(bad_noncontrib.sum()))
                # ------------------------      
                
                
                # Contributors are those that (1) have level k and valid, AND (2) are finite
                finite = np.isfinite(hk)
                contrib = m_k & finite

                # Sum over finite contributors only; non-contributors contribute 0
                hk_finite = np.where(contrib, hk, 0.0)
                sum_k = A @ hk_finite
                cnt_k = A @ contrib.astype(float)

                # Debug/asserts
                _assert_finite("sum_k", sum_k)
                _assert_finite("cnt_k", cnt_k)
                _assert_nonzero("cnt_k at S (bottom)", cnt_k, where_mask=S)

                # Neighbor mean like np.nanmean: divide where we have contributors
                vave = np.empty_like(sum_k)
                vave.fill(np.nan)
                np.divide(sum_k, cnt_k, out=vave, where=(cnt_k > 0.0))
                # end fix

                # Data-fit pull to hhat at k (only where valid)
                pull = np.zeros(n, dtype=float)
                valid_k = m_k                         # exactly the same validity
                pull[valid_k] = -lamD * (h[valid_k, k] - hhat[valid_k, k])

                # Bottom thickness anchoring (same formula as before)
                # (h_k - h_{k-1}) -> tbottom
                with_lower = S & (k >= 1) & valid[:, k-1]
                thick = np.zeros(n, dtype=float)
                thick[with_lower] = h[with_lower, k] - h[with_lower, k-1]
                pull[with_lower] += -lamB * (thick[with_lower] - tbottom[with_lower])

                # Explicit Euler update for this k, only on the S nodes
                delta = dt * (kappa * (vave - h[:, k]) + pull)
                h[S, k] = h[S, k] + delta[S]

        # Guards after each outer bottom iteration (unchanged)
        for i in range(n):
            Ni = int(Nlevels[i])
            if Ni < 2:
                continue
            # No-collapse per layer in bottom band
            lo = max(1, Ni - params.K_bottom)
            for k in range(lo, Ni):
                if valid[i, k] and valid[i, k-1]:
                    tpre = pre_thick[i, k]
                    tnow = h[i, k] - h[i, k-1]
                    tmin_guard = max(tmin_arr[i], no_collapse * tpre)
                    if tnow < tmin_guard:
                        h[i, k] = h[i, k-1] + tmin_guard
            enforce_node_monotone(i)

    # Mid-column smoothing (W neighbors-only + caller ghost fill):
    #   Here we use the row-stochastic W (neighbors only, no self) and the
    #   caller provides ghost values: Hvirt = h[:,k] where level k exists,
    #   otherwise the bed (depth) for that column. This matches the existing
    #   mid-column logic and the one-pass smoother behavior.
    # -------------------------------------------------------------------------
    # Mid-column smoothing & target fit
    # Cache/build once per call (above mid loop; place near other precomputes)
    # Mid-column smoothing & target fit
    # Cache/build once per call (above mid loop; place near other precomputes)
    W = _get_row_stochastic_W(mesh, n)
    Ni_arr = Nlevels.astype(np.int32)

    if DEBUG_LSC2:
        # mask of entries that exist: k < Nlevels[i]
        idx = np.arange(maxN)
        exists = idx[None, :] < Ni_arr[:, None]         # shape (n, maxN), bool
        _assert_finite("h before MID loop (existing levels only)", h, where_mask=exists)

    for it in range(params.iters_mid):
        for k in range(1, maxN - 1):
            # Nodes that have this level
            have_k = (Ni_arr > k)
            if not np.any(have_k):
                continue

            # Valid targets at this level
            valid_k = valid[:, k] & have_k


            # Debug 1) Inputs sanity
            _dbg("MID it={} k={} | have_k:{} valid_k:{}", it, k, int(have_k.sum()), int(valid_k.sum()))
            _assert_finite("h[:,k] pre (mid)", h[:, k], where_mask=have_k)
            _assert_finite("hhat[:,k] (mid)", hhat[:, k], where_mask=valid_k)

            # Ghost fill: use bed where column lacks level k
            # 2) Ghost fill for averaging and neighbor mean
            Hvirt = np.where(have_k, h[:, k], depth)
            _assert_finite("Hvirt (mid)", Hvirt)  # MUST be finite
            h_prev_k = h[:, k].copy()             # snapshot for nan-diff

            # Neighbor mean via sparse matvec
            nbr_mean = W @ Hvirt
            _assert_finite("nbr_mean (mid)", nbr_mean, where_mask=have_k)


            # Per-node zone scaling (vectorized)
            # near bed
            near_bed = have_k & (k >= (Ni_arr - params.K_bottom))
            # just above near-bed band
            near_bed_above = have_k & (k >= (Ni_arr - params.K_bottom - 3)) & (k < (Ni_arr - params.K_bottom))
            # near top
            near_top = have_k & (k <= params.keep_top_thin)

            zone_scale = np.ones(n, dtype=float)
            zone_scale[near_bed] = 0.6
            zone_scale[near_bed_above] = 0.8
            zone_scale[near_top] = params.top_relax

            # Target pull (only where valid)
            pull = np.zeros(n, dtype=float)
            pull[valid_k] = -params.lam_D * (h[valid_k, k] - hhat[valid_k, k])
            _assert_finite("pull (mid)", pull, where_mask=have_k)

            # Explicit Euler update (only nodes that truly have this level)
            alpha = dt * zone_scale * kappa
            delta = alpha * (nbr_mean - h[:, k]) + dt * pull
            # debug block
            _assert_finite("delta (mid)", delta, where_mask=have_k)
            
            # dEBUG Extra: catch first brand-new NaN introduced by this write
            # (define _first_new_nan helper once, as we discussed)
            # _first_new_nan(name, prev, cur, mask=None)
            upd = have_k
            h[upd, k] = h[upd, k] + delta[upd]
            _first_new_nan(f"MID it={it} k={k}", h_prev_k, h[:, k], mask=have_k)


        # Enforce ordering & guards after each outer mid-iter
        for i in range(n):
            Ni = int(Nlevels[i])
            if Ni < 2:
                continue
            enforce_node_monotone(i)

        # Enforce ordering & guards (only min thickness globally here)
        for i in range(n):
            Ni = int(Nlevels[i])
            if Ni < 2: continue
            enforce_node_monotone(i)

    if getattr(params, "dump_presmooth_npz", None):
        np.savez(
            params.dump_presmooth_npz,
            h=h.copy(),
            Nlevels=Nlevels.astype(np.int32),
            depth=depth.astype(np.float64),
        )
        logger.debug("[dump] wrote %s", params.dump_presmooth_npz)

    # --- after bottom & mid are done, right BEFORE finishing smoother ---
    if DEBUG_LSC2:
        idx = np.arange(maxN)
        Ni_arr = Nlevels.astype(np.int32)
        exists = idx[None, :] < Ni_arr[:, None]     # only real interfaces
        _assert_finite("h pre FINISH (exist)", h, where_mask=exists)


    # One or more Jacobi passes with ghost-bed fill
    for _ in range(int(getattr(params, "sm_iters", 1))):
        smooth_interfaces_onepass_with_bed_fill_sparse(
            mesh, h, Nlevels, depth,
            alpha_top=float(getattr(params, "sm_alpha_top", 0.1)),
            alpha_bot=float(getattr(params, "sm_alpha_bot", 0.40)),
            power=float(getattr(params, "sm_power", 1.0))
        )

    # --- immediately AFTER finishing smoother ---
    if DEBUG_LSC2:
        _assert_finite("h post FINISH (exist)", h, where_mask=exists)


    # Finalization: exact endpoints + guards (+ strict nudge)
    for i in range(n):
        Ni = int(Nlevels[i]); D = depth[i]
        if Ni < 2:
            continue

        # Dry / non-positive-depth nodes: assign uniform degenerate sigma
        if D <= 0.0:
            D_eff = max(D, 1e-6)
            for k in range(Ni):
                h[i, k] = D_eff * k / (Ni - 1)
            continue

        # Exception nodes (shallow, forced Nlevels by hysteresis): reassert
        # uniform sigma. The smoother pulls these toward deep neighbors;
        # uniform is the correct answer for them.
        if uniform_sigma_mask is not None and uniform_sigma_mask[i]:
            for k in range(Ni):
                h[i, k] = k * D / (Ni - 1)
            continue

        # Surface exactly 0
        h[i, 0] = 0.0

        # Forward guard: enforce minimum thickness everywhere
        for k in range(1, Ni):
            lb = h[i, k - 1] + tmin_arr[i]
            if h[i, k] < lb:
                h[i, k] = lb

        # Pin bed exactly to D
        h[i, Ni - 1] = D

        # Overflow fix: if smoother pushed levels above D, or the bed
        # cell would be thinner than tmin, redistribute proportionally
        # between a valid base and the bed.  Choose the split point so
        # that the redistributed spacing >= tmin.
        if Ni >= 3 and h[i, Ni - 2] > D - tmin_arr[i]:
            tmin_i = tmin_arr[i]
            # Scan downward: find the deepest level k whose value allows
            # (Ni-1-k) layers of at least tmin between h[k] and D.
            k_overflow = 1  # fallback: redistribute everything from level 1
            for k in range(Ni - 2, 0, -1):
                if h[i, k] <= D - (Ni - 1 - k) * tmin_i:
                    k_overflow = k + 1
                    break
            # Linearly redistribute from h[k_overflow-1] to D
            h_base = h[i, k_overflow - 1]
            n_redist = Ni - k_overflow  # number of levels to redistribute
            for j in range(n_redist):
                h[i, k_overflow + j] = h_base + (j + 1) * (D - h_base) / n_redist

        # Strict separation nudge (fp safety)
        for k in range(1, Ni):
            if not (h[i, k] > h[i, k - 1]):
                h[i, k] = np.nextafter(h[i, k - 1], np.float64('inf'))

    # --- final tripwire before return ---
    if DEBUG_LSC2:
        nan_exist = (~np.isfinite(h)) & exists
        if nan_exist.any():
            i, k = np.argwhere(nan_exist)[0]
            raise RuntimeError(f"[EXIT fit_interfaces] NaN/Inf on existing interface at (node {i}, level {k})")

    return h, valid


# -----------------------------------------------------------------------------
# E) Build sigma and ensure SCHISM formatting expectations
# -----------------------------------------------------------------------------

def build_sigma_from_interfaces(depth: np.ndarray,
                                h: np.ndarray,
                                valid: np.ndarray,
                                Nlevels: np.ndarray) -> np.ndarray:
    """Compute per-node sigma values (0..1) with right padding for SCHISM.

    sigma[i, :Ni] = h[i, :Ni] / depth[i]
    Remaining entries are filled with NaN; writers will convert to -9.0 as needed.
    """
    n, maxN = h.shape
    sigma = np.full((n, maxN), np.nan, dtype=float)
    D = np.maximum(depth, 1e-6)
    for i in range(n):
        Ni = int(Nlevels[i])
        if Ni < 2: 
            raise ValueError(f"Bad skip at node {i}")
            #continue
        sigma[i, :Ni] = h[i, :Ni] / D[i]
        # Guard against tiny negative or >1 due to numerics
        sigma[i, :Ni] = np.clip(sigma[i, :Ni], 0.0, 1.0)
    return sigma

def apply_interface_smoother_only(mesh,
                                  h_in: np.ndarray,
                                  Nlevels: np.ndarray,
                                  depth: np.ndarray,
                                  alpha_top: float = 0.25,
                                  alpha_bot: float = 0.05,
                                  power: float = 1.0,
                                  print_delta: bool = True):
    """
    Run ONLY the interface-level smoother on a frozen snapshot (h_in, Nlevels, depth),
    and return the smoothed h plus a per-level max|Δh| diagnostic.

    No fitting, no guards. Use this to iterate gains quickly.
    """
    h = h_in.copy()
    h_before = h.copy()

    # Ensure endpoints for stable ghost-fill behavior
    n, maxN = h.shape
    for i in range(n):
        Ni = int(Nlevels[i])
        if Ni < 2:
            continue
        if depth[i] <= 0.0:
            h[i, :Ni] = 0.0
            continue
        h[i, 0] = 0.0
        h[i, Ni-1] = depth[i]

    # One Jacobi pass (CSR, ghost-bed fill)
    smooth_interfaces_onepass_with_bed_fill_sparse(
        mesh, h, Nlevels, depth,
        alpha_top=alpha_top,
        alpha_bot=alpha_bot,
        power=power
    )

    if print_delta:
        with np.errstate(invalid="ignore"):
            dk = np.nanmax(np.abs(h - h_before), axis=0)
        logger.info("smoother-only max|Δh| by level (first 16): %s",
                     ", ".join(f"{v:.3f}" for v in dk[:min(16, len(dk))]))

    return h

# -----------------------------------------------------------------------------
# Optional boundary priors (soft control near selected boundaries)
# -----------------------------------------------------------------------------

@dataclass
class BoundaryPriors:
    upstream_nodes: Optional[np.ndarray] = None
    upstream_target_levels: int = 2   # => 1 layer
    upstream_rings: int = 5           # number of edge-rings over which to fade
    ocean_nodes: Optional[np.ndarray] = None
    ocean_target_levels: Optional[int] = None  # if None, no ocean prior on levels
    ocean_rings: int = 8

def _rings_from_seeds(mesh, seeds: np.ndarray, max_rings: int) -> np.ndarray:
    """Return integer ring distance from the seed set per node (0 at seeds, inf if unreachable)."""
    import collections
    n = mesh.nodes.shape[0]
    dist = np.full(n, np.inf, dtype=float)
    q = collections.deque()
    for s in np.asarray(seeds, dtype=int):
        if 0 <= s < n:
            dist[s] = 0.0
            q.append(s)
    while q:
        i = q.popleft()
        if dist[i] >= max_rings:
            continue
        for j in mesh.get_neighbor_nodes(i):
            if dist[j] == np.inf:
                dist[j] = dist[i] + 1.0
                if dist[j] < max_rings:
                    q.append(j)
    return dist

def apply_boundary_priors_to_Ltilde(mesh,
                                    Ltilde: np.ndarray,
                                    priors: Optional[BoundaryPriors]) -> np.ndarray:
    """Blend Ltilde toward boundary target levels over a few rings (soft, shock-free)."""
    if priors is None:
        return Ltilde
    out = Ltilde.copy()
    n = len(out)

    def blend_toward(nodes, rings, target_levels):
        if nodes is None or len(nodes) == 0 or rings <= 0 or target_levels is None:
            return
        dist = _rings_from_seeds(mesh, np.asarray(nodes, dtype=int), rings)
        w = np.clip((rings - dist) / max(rings, 1), 0.0, 1.0)
        # Smooth taper (cosine) to avoid kinks
        w = 0.5 * (1.0 - np.cos(np.pi * w))
        out[:] = (1.0 - w) * out + w * float(target_levels)

    # Upstream collapse to 1 layer (2 levels)
    blend_toward(priors.upstream_nodes, priors.upstream_rings, priors.upstream_target_levels)
    # Ocean control (optional)
    blend_toward(priors.ocean_nodes, priors.ocean_rings, priors.ocean_target_levels)

    return out

@dataclass
class PipelineParams:
    # Smoothing & quantization
    priors: Optional[BoundaryPriors] = None
    # Smoothing & quantization
    L_smooth_passes: int = 2
    L_smooth_kappa: float = 0.5
    L_smooth_method: str = "isotropic"
    L_smooth_L_scale: float = 1.25
    L_smooth_depth_scale: float = 4.0
    L_smooth_length_power: float = 0.0
    constraint_diffusion_enable: bool = False
    constraint_diffusion_passes: int = 40
    constraint_diffusion_weight: float = 0.6
    constraint_diffusion_method: str = "same"
    constraint_diffusion_L_scale: Optional[float] = None
    constraint_diffusion_depth_scale: Optional[float] = None
    constraint_diffusion_length_power: Optional[float] = None
    constraint_diffusion_tol: float = 0.0
    integer_cleanup_enable: bool = False
    integer_cleanup_max_nodes: int = 6
    integer_cleanup_max_passes: int = 5
    hysteresis: HysteresisParams = field(default_factory=HysteresisParams)
    fit: FitParams = field(default_factory=FitParams)
    # Polygon region constraints (per-node arrays, or None)
    n_min: Optional[np.ndarray] = None   # int32, minimum Nlevels per node
    n_max: Optional[np.ndarray] = None   # int32, maximum Nlevels per node




def _write_debug_node_attr(mesh, debug: Optional[Dict[str, str]], key: str, values) -> None:
    """Write a debug GR3 node-attribute field if requested.

    Kept deliberately small so diagnostics can be emitted immediately after
    the corresponding quantity is computed, rather than at the end of the
    full vgrid pipeline.
    """
    if not debug or key not in debug:
        return
    write_mesh(mesh, debug[key], node_attr=np.asarray(values, dtype=float))
    logger.info("   diagnostic %s → %s", key, debug[key])

def run_pipeline(mesh,
                 depth: np.ndarray,
                 sizefun: VerticalSizeFunction,
                 eta: float,
                 pp: PipelineParams = PipelineParams(),
                 debug: Optional[Dict[str, str]] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run the full A–E pipeline and return sigma, Nlevels, h, and tmin_arr.

    Returns
    -------
    sigma : (n, maxN) float
    Nlevels : (n,) int
    h : (n, maxN) float
    tmin_arr : (n,) float  — per-node minimum layer thickness

    If debug is a dict of filenames, optional GR3s will be written for inspection.
    Diagnostics are emitted as soon as their source quantity is available:
    Lstar/Ltilde/Nmin/Nmax/Nlevels/Nlayers/uniform_sigma/tbottom are written
    before interface fitting; bottom_thickness is written after fitted interfaces
    exist.
    """

    t0 = time.perf_counter()
    logger.info("A) Computing L* ...")

    # A) L*
    dpos = np.maximum(depth, 1e-6)
    Lstar = compute_Lstar(dpos, sizefun, mesh=mesh)
    #Lstar = compute_Lstar(depth, sizefun, mesh=mesh)
    Ltilde = smooth_Lstar_on_mesh(
        mesh, Lstar, dpos,
        method=pp.L_smooth_method,
        passes=pp.L_smooth_passes,
        kappa=pp.L_smooth_kappa,
        dt=1.0,
        L_scale=pp.L_smooth_L_scale,
        depth_scale=pp.L_smooth_depth_scale,
        length_power=pp.L_smooth_length_power,
    )
    logger.info(
        "   L* smoothing: method=%s passes=%d kappa=%.3g L_scale=%.3g depth_scale=%.3g length_power=%.3g",
        pp.L_smooth_method, pp.L_smooth_passes, pp.L_smooth_kappa,
        pp.L_smooth_L_scale, pp.L_smooth_depth_scale, pp.L_smooth_length_power,
    )


    # Apply optional boundary priors softly (no shocks)
    Ltilde = apply_boundary_priors_to_Ltilde(mesh, Ltilde, pp.priors)

    # Apply polygon constraints to the continuous level field.  The legacy path
    # only clips inside constrained polygons.  The optional constraint-diffusion
    # path treats fixed constraints as Dirichlet nodes and range/min/max
    # constraints as projected box bounds during smoothing, so constrained
    # regions continue to influence adjacent transition zones.
    if pp.n_min is not None or pp.n_max is not None:
        if pp.constraint_diffusion_enable:
            method = str(pp.constraint_diffusion_method or "same")
            if method.strip().lower() == "same":
                method = pp.L_smooth_method
            L_scale = (pp.L_smooth_L_scale if pp.constraint_diffusion_L_scale is None
                       else pp.constraint_diffusion_L_scale)
            depth_scale = (pp.L_smooth_depth_scale if pp.constraint_diffusion_depth_scale is None
                           else pp.constraint_diffusion_depth_scale)
            length_power = (pp.L_smooth_length_power if pp.constraint_diffusion_length_power is None
                            else pp.constraint_diffusion_length_power)
            Ltilde = diffuse_Ltilde_with_constraints(
                mesh, Ltilde, dpos,
                n_min=pp.n_min, n_max=pp.n_max,
                passes=pp.constraint_diffusion_passes,
                weight=pp.constraint_diffusion_weight,
                method=method,
                L_scale=L_scale,
                depth_scale=depth_scale,
                length_power=length_power,
                tol=pp.constraint_diffusion_tol,
            )
        else:
            out = Ltilde.copy()
            n_raised = 0
            n_lowered = 0

            if pp.n_min is not None:
                m = out < pp.n_min
                n_raised = int(np.count_nonzero(m))
                out[m] = pp.n_min[m].astype(float)

            if pp.n_max is not None:
                m = out > pp.n_max
                n_lowered = int(np.count_nonzero(m))
                out[m] = pp.n_max[m].astype(float)

            Ltilde = out
            n_constrained = int((pp.n_min is not None and (pp.n_min > 2).sum()) or 0) + \
                            int((pp.n_max is not None and (pp.n_max < 99999).sum()) or 0)
            logger.info(
                "   polygon constraints: applied inside polygons only at %d constrained nodes "
                "(raised Ltilde at %d, lowered at %d)",
                n_constrained, n_raised, n_lowered,
            )

    _write_debug_node_attr(mesh, debug, "Lstar", Lstar)
    _write_debug_node_attr(mesh, debug, "Ltilde", Ltilde)
    _write_debug_node_attr(mesh, debug, "Lconstrained", Ltilde)
    if pp.n_min is not None:
        _write_debug_node_attr(mesh, debug, "Nmin", pp.n_min)
    if pp.n_max is not None:
        _write_debug_node_attr(mesh, debug, "Nmax", pp.n_max)

    logger.info("   done in %.2f s", time.perf_counter() - t0)

    # B) hysteretic integerization (levels)
    t1 = time.perf_counter()
    logger.info("B) Hysteretic integerization (levels) ...")
    Nlevels = quantize_levels_with_hysteresis(mesh, Ltilde, dpos, pp.hysteresis)

    # Hard-clip Nlevels to polygon constraints (after quantization)
    if pp.n_min is not None:
        clipped_up = int((Nlevels < pp.n_min).sum())
        Nlevels = np.maximum(Nlevels, pp.n_min)
        if clipped_up:
            logger.info("   hard-clip: raised Nlevels at %d nodes to meet n_min", clipped_up)
    if pp.n_max is not None:
        clipped_dn = int((Nlevels > pp.n_max).sum())
        Nlevels = np.minimum(Nlevels, pp.n_max)
        if clipped_dn:
            logger.info("   hard-clip: lowered Nlevels at %d nodes to meet n_max", clipped_dn)

    _write_debug_node_attr(mesh, debug, "Nlevels_raw", Nlevels)

    if pp.integer_cleanup_enable:
        N_before_cleanup = Nlevels.copy()
        Nlevels = cleanup_integer_components(
            mesh, Nlevels,
            max_component_nodes=pp.integer_cleanup_max_nodes,
            max_passes=pp.integer_cleanup_max_passes,
            n_min=pp.n_min,
            n_max=pp.n_max,
        )
        n_changed = int(np.count_nonzero(Nlevels != N_before_cleanup))
        if n_changed:
            delta = Nlevels.astype(np.int32) - N_before_cleanup.astype(np.int32)
            logger.info(
                "   integer cleanup: changed %d nodes (min Δ=%d max Δ=%d)",
                n_changed, int(delta.min()), int(delta.max()),
            )
        else:
            logger.info("   integer cleanup: no Nlevels changed")

    # --- Per-node tmin: relax where hysteresis asks for more levels than
    #     the global tmin allows, instead of capping Nlevels down. ---
    feas_cap = (np.floor(dpos / max(pp.fit.tmin, 1e-9)).astype(int) + 1)
    feas_cap = np.maximum(feas_cap, 2)
    # Nodes where hysteresis wants more levels than global tmin can provide
    uniform_sigma_mask = Nlevels > feas_cap
    # Build per-node tmin: global value everywhere, relaxed for exception nodes
    tmin_arr = np.full(len(dpos), pp.fit.tmin, dtype=float)
    exc_nodes = np.where(uniform_sigma_mask)[0]
    for i in exc_nodes:
        tmin_arr[i] = dpos[i] / max(int(Nlevels[i]) - 1, 1)
    n_exc = int(uniform_sigma_mask.sum())
    if n_exc > 0:
        exc_depths = dpos[uniform_sigma_mask]
        logger.info(
            "   %d nodes exceed global tmin feasibility "
            "(depth range %.2f-%.2f m); relaxing tmin locally, using uniform sigma targets.",
            n_exc, exc_depths.min(), exc_depths.max(),
        )
    Nlayers = (Nlevels - 1).astype(int)
    maxN = int(np.max(Nlevels))
    _write_debug_node_attr(mesh, debug, "Nlevels", Nlevels)
    _write_debug_node_attr(mesh, debug, "Nlayers", Nlayers)
    _write_debug_node_attr(mesh, debug, "uniform_sigma", uniform_sigma_mask)
    logger.info("   done in %.2f s", time.perf_counter() - t1)
    
    t2 = time.perf_counter()
    # C) targets
    logger.info("C) Building target interfaces ...")

    maxN = int(np.max(Nlevels))
    hhat, valid = build_target_interfaces(sizefun, depth, Nlevels, eta=eta, mesh=mesh,
                                          uniform_sigma_mask=uniform_sigma_mask)

    # bottom thickness target & smoothing
    tb = smooth_bottom_thickness(depth, hhat, valid, pp.fit.tmin, mesh)    # dpos to depth
    _write_debug_node_attr(mesh, debug, "tbottom", tb)
    logger.info("   done in %.2f s", time.perf_counter() - t2)


    # D) fit interfaces
    t3 = time.perf_counter()
    logger.info("D) Fitting interfaces ...")
    h, valid2 = fit_interfaces(mesh, depth, Nlevels, hhat, valid, tb, pp.fit,
                               tmin_arr=tmin_arr,
                               uniform_sigma_mask=uniform_sigma_mask)
    logger.info("   done in %.2f s", time.perf_counter() - t3)

    # E) sigma
    t4 = time.perf_counter()
    logger.info("E) Building sigma ...")
    sigma = build_sigma_from_interfaces(depth, h, valid2, Nlevels)              # todo dpos to depth

    logger.info("   done in %.2f s", time.perf_counter() - t4)

    logger.info("Total pipeline time: %.2f s", time.perf_counter() - t0)

    # Optional debug output that depends on fitted interfaces
    if debug and "bottom_thickness" in debug:
        bt = np.zeros(len(Nlevels), dtype=float)
        for i in range(len(bt)):
            Ni = int(Nlevels[i])
            if Ni >= 2:
                bt[i] = h[i, Ni-1] - h[i, Ni-2]
        _write_debug_node_attr(mesh, debug, "bottom_thickness", bt)

    return sigma, Nlevels, h, tmin_arr
