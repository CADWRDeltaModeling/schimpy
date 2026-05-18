#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Validate a SCHISM LSC2 vgrid.in file.

Checks
------
1. Padding consistency: -9.0 sentinel appears exactly at levels < kbp.
2. Strict monotonicity: sigma values increase from kbp to nvrt.
3. Output-precision monotonicity: sigma survives %14.6f rounding.
4. Thin-layer detection: flag layers thinner than a fraction of uniform.
5. Side-midpoint monotonicity: the z-check SCHISM actually performs.

Designed for use as both a CLI (``check_vgrid``) and in test suites::

    from schimpy.check_vgrid import check_vgrid
    issues = check_vgrid("vgrid.in", hgrid="hgrid.gr3", eta=1.5)
    assert len(issues) == 0, issues
"""

import logging
from dataclasses import dataclass
from typing import List

import click
import numpy as np

from schimpy.logging_config import resolve_loglevel, configure_logging

logger = logging.getLogger(__name__)


@dataclass
class VgridIssue:
    """Single validation failure."""

    check: str
    node: int  # 0-based
    level: int = -1  # 0-based internal, or -1 if node-level
    detail: str = ""
    x: float = float("nan")
    y: float = float("nan")

    def __str__(self):
        loc = f"node {self.node}"
        if self.level >= 0:
            loc += f" level {self.level}"
        if np.isfinite(self.x) and np.isfinite(self.y):
            loc += f" ({self.x:.2f}, {self.y:.2f})"
        return f"[{self.check}] {loc}: {self.detail}"


def check_vgrid(
    vgrid_path,
    *,
    vgrid_version="5.10",
    hgrid=None,
    eta=0.0,
    thin_fraction=0.1,
    max_report=50,
) -> List[VgridIssue]:
    """Validate a vgrid.in file and return a list of issues.

    Parameters
    ----------
    vgrid_path : str or Path
        Path to vgrid.in.
    vgrid_version : str
        '5.8' or '5.10'.
    hgrid : str or Path, optional
        Path to hgrid.gr3.  Required for thin-layer and side-midpoint checks.
    eta : float
        Reference free-surface elevation for depth computation.
    thin_fraction : float
        Flag layers thinner than ``thin_fraction * depth / nlayers``.
    max_report : int
        Cap on issues per check category.

    Returns
    -------
    list of VgridIssue
    """
    from schimpy.schism_vertical_mesh import read_vmesh

    vmesh = read_vmesh(str(vgrid_path), vgrid_version=vgrid_version)
    sigma = vmesh.sigma  # (n_nodes, nvrt), left-justified, NaN-padded
    kbps = vmesh.kbps  # (n_nodes,) 0-based
    nvrt = vmesh.n_vert_levels()
    n_nodes = sigma.shape[0]

    issues: List[VgridIssue] = []

    def _cap(check_name):
        return sum(1 for iss in issues if iss.check == check_name) >= max_report

    # ------------------------------------------------------------------
    # Load hgrid (coordinates + depth) if provided
    # ------------------------------------------------------------------
    mesh = None
    depth = None
    node_xy = None
    if hgrid is not None:
        from schimpy.schism_mesh import read_mesh

        mesh = read_mesh(str(hgrid))
        node_xy = mesh.nodes[:, :2]  # (n_nodes, 2) x,y
        h0 = mesh.nodes[:, 2]
        depth = eta + h0

    def _xy(i):
        """Return (x, y) for node *i*, or (NaN, NaN) if no hgrid."""
        if node_xy is not None and i < len(node_xy):
            return float(node_xy[i, 0]), float(node_xy[i, 1])
        return float("nan"), float("nan")

    # ------------------------------------------------------------------
    # 1. Padding / kbp consistency
    # ------------------------------------------------------------------
    logger.info(
        "Check 1: padding and kbp consistency (%d nodes, nvrt=%d)", n_nodes, nvrt
    )
    for i in range(n_nodes):
        if _cap("padding"):
            break
        kbp = kbps[i]
        n_levels = nvrt - kbp
        valid = sigma[i, :n_levels]
        padding = sigma[i, n_levels:]

        n_nan_in_valid = int(np.isnan(valid).sum())
        if n_nan_in_valid > 0:
            issues.append(
                VgridIssue(
                    "padding",
                    i,
                    detail=f"kbp={kbp} implies {n_levels} valid levels, "
                    f"but {n_nan_in_valid} are NaN",
                    x=_xy(i)[0],
                    y=_xy(i)[1],
                )
            )

        if padding.size > 0:
            n_finite_in_pad = int(np.isfinite(padding).sum())
            if n_finite_in_pad > 0:
                issues.append(
                    VgridIssue(
                        "padding",
                        i,
                        detail=f"kbp={kbp} implies {kbp} padding slots, "
                        f"but {n_finite_in_pad} are finite",
                        x=_xy(i)[0],
                        y=_xy(i)[1],
                    )
                )

    n_pad = sum(1 for iss in issues if iss.check == "padding")
    logger.info("  padding issues: %d", n_pad)

    # ------------------------------------------------------------------
    # 2. Strict monotonicity (full precision)
    # ------------------------------------------------------------------
    logger.info("Check 2: strict monotonicity (full precision)")
    for i in range(n_nodes):
        if _cap("monotone"):
            break
        n_levels = nvrt - kbps[i]
        if n_levels < 2:
            continue
        s = sigma[i, :n_levels]
        if not np.all(np.isfinite(s)):
            continue
        d = np.diff(s)
        bad = np.where(d <= 0.0)[0]
        for b in bad:
            issues.append(
                VgridIssue(
                    "monotone",
                    i,
                    level=int(b),
                    detail=f"sigma[{b}]={s[b]:.10f}, sigma[{b+1}]={s[b+1]:.10f}, "
                    f"diff={d[b]:.2e}",
                    x=_xy(i)[0],
                    y=_xy(i)[1],
                )
            )
            if _cap("monotone"):
                break

    n_mono = sum(1 for iss in issues if iss.check == "monotone")
    logger.info("  monotonicity issues: %d", n_mono)

    # ------------------------------------------------------------------
    # 3. Output-precision monotonicity (6 decimal places)
    # ------------------------------------------------------------------
    logger.info("Check 3: output-precision monotonicity (6 decimal places)")
    for i in range(n_nodes):
        if _cap("precision"):
            break
        n_levels = nvrt - kbps[i]
        if n_levels < 2:
            continue
        s = sigma[i, :n_levels]
        if not np.all(np.isfinite(s)):
            continue
        sr = np.round(s, 6)
        d = np.diff(sr)
        bad = np.where(d <= 0.0)[0]
        for b in bad:
            issues.append(
                VgridIssue(
                    "precision",
                    i,
                    level=int(b),
                    detail=f"rounded sigma[{b}]={sr[b]:.6f} == "
                    f"sigma[{b+1}]={sr[b+1]:.6f}",
                    x=_xy(i)[0],
                    y=_xy(i)[1],
                )
            )
            if _cap("precision"):
                break

    n_prec = sum(1 for iss in issues if iss.check == "precision")
    logger.info("  precision issues: %d", n_prec)

    # ------------------------------------------------------------------
    # 4. Thin-layer detection
    # ------------------------------------------------------------------
    if mesh is not None:
        logger.info(
            "Check 4: thin-layer detection (threshold=%.0f%% of uniform)",
            thin_fraction * 100,
        )
        for i in range(n_nodes):
            if _cap("thin_layer"):
                break
            n_levels = nvrt - kbps[i]
            nlayers = n_levels - 1
            if nlayers < 1 or i >= len(depth):
                continue
            D = depth[i]
            if D <= 0.0:
                continue
            s = sigma[i, :n_levels]
            if not np.all(np.isfinite(s)):
                continue
            uniform_dz = D / nlayers
            threshold = thin_fraction * uniform_dz
            # sigma in [-1,0]; physical dz = -D * dsigma
            ds = np.diff(s) * (-D)
            bad = np.where(np.abs(ds) < threshold)[0]
            for b in bad:
                issues.append(
                    VgridIssue(
                        "thin_layer",
                        i,
                        level=int(b),
                        detail=f"dz={abs(ds[b]):.6f} m < {threshold:.6f} m "
                        f"({thin_fraction:.0%} of uniform "
                        f"{uniform_dz:.4f} m), "
                        f"depth={D:.2f} m, nlayers={nlayers}",
                        x=_xy(i)[0],
                        y=_xy(i)[1],
                    )
                )
                if _cap("thin_layer"):
                    break

        n_thin = sum(1 for iss in issues if iss.check == "thin_layer")
        logger.info("  thin-layer issues: %d", n_thin)
    else:
        logger.info("Check 4: thin-layer detection skipped (no hgrid)")

    # ------------------------------------------------------------------
    # 5. Side-midpoint monotonicity (the SCHISM "Weird side" check)
    # ------------------------------------------------------------------
    if mesh is not None:
        logger.info("Check 5: side-midpoint monotonicity")
        # mesh.edges is (n_edges, 5): col 0,1 are node indices
        edges = mesh.edges

        for edge_i in range(edges.shape[0]):
            if _cap("side_midpoint"):
                break
            n1, n2 = int(edges[edge_i, 0]), int(edges[edge_i, 1])
            if n1 >= n_nodes or n2 >= n_nodes:
                continue
            kbp1, kbp2 = kbps[n1], kbps[n2]
            nl1, nl2 = nvrt - kbp1, nvrt - kbp2
            # Side bottom index (SCHISM convention: min of two kbps)
            kbs = min(kbp1, kbp2)
            n_side_levels = nvrt - kbs
            if n_side_levels < 2:
                continue

            # Reproduce SCHISM's zs(k) = avg of sigma at max(k, kbp)
            zs_prev = None
            for k in range(n_side_levels):
                schism_k = kbs + k
                k1_int = min(max(schism_k - kbp1, 0), nl1 - 1)
                k2_int = min(max(schism_k - kbp2, 0), nl2 - 1)
                s1 = sigma[n1, k1_int]
                s2 = sigma[n2, k2_int]

                if not (np.isfinite(s1) and np.isfinite(s2)):
                    continue

                zs = (s1 + s2) / 2.0
                if zs_prev is not None and zs <= zs_prev:
                    issues.append(
                        VgridIssue(
                            "side_midpoint",
                            n1,
                            level=schism_k,
                            detail=f"side ({n1},{n2}) level {schism_k}: "
                            f"zs={zs:.10f} <= prev={zs_prev:.10f}, "
                            f"s1={s1:.8f} s2={s2:.8f}",
                            x=_xy(n1)[0],
                            y=_xy(n1)[1],
                        )
                    )
                    if _cap("side_midpoint"):
                        break
                zs_prev = zs

        n_side = sum(1 for iss in issues if iss.check == "side_midpoint")
        logger.info("  side-midpoint issues: %d", n_side)
    else:
        logger.info("Check 5: side-midpoint monotonicity skipped (no hgrid)")

    logger.info("Total issues: %d", len(issues))
    return issues


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------


@click.command()
@click.argument("vgrid", type=click.Path(exists=True))
@click.option(
    "--vgrid_version",
    type=str,
    default="5.10",
    help="vgrid format version ('5.8' or '5.10')",
)
@click.option(
    "--hgrid",
    type=click.Path(exists=True),
    default=None,
    help="Path to hgrid.gr3 (enables thin-layer and side checks)",
)
@click.option("--eta", type=float, default=0.0, help="Reference free-surface elevation")
@click.option(
    "--thin-fraction",
    type=float,
    default=0.1,
    help="Flag layers thinner than this fraction of uniform dz",
)
@click.option(
    "--max-report", type=int, default=50, help="Max issues reported per category"
)
@click.option("--debug", is_flag=True, help="Debug-level logging")
@click.option("--quiet", is_flag=True, help="Suppress console output")
@click.help_option("-h", "--help")
def check_vgrid_cli(vgrid, vgrid_version, hgrid, eta, thin_fraction, max_report, debug, quiet):
    """Validate a SCHISM vgrid.in file.

    Reports padding inconsistencies, non-monotone sigma, precision
    collapses, and optionally thin layers and side-midpoint monotonicity
    (when --hgrid is given).  Node x,y coordinates are included in the
    output when --hgrid is provided.

    Node IDs in the output are 0-based (Python convention).

    Exit code 0 if no issues found, 1 otherwise.
    """
    level, console = resolve_loglevel(debug=debug, quiet=quiet)
    configure_logging(
        package_name="schimpy",
        level=level,
        console=console,
    )

    issues = check_vgrid(
        vgrid,
        vgrid_version=vgrid_version,
        hgrid=hgrid,
        eta=eta,
        thin_fraction=thin_fraction,
        max_report=max_report,
    )

    if issues:
        for iss in issues:
            click.echo(str(iss))
        raise SystemExit(1)
    else:
        click.echo("vgrid OK")


if __name__ == "__main__":
    check_vgrid_cli()
