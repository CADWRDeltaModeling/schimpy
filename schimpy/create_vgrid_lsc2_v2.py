#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create LSC2 v2 vertical grid (BETA).

This module implements the v2 vertical grid generator for SCHISM,
using the LSC2 v2 pipeline (locally adaptive, globally smooth).

It can be called from prepare_schism (dict-based YAML config) or
standalone via the CLI entry point.

The v2 generator is currently in beta. For the stable legacy generator,
use vgrid_generator_version: v1.
"""

import logging
import os
from pathlib import Path

import click
import numpy as np

from schimpy.lsc2 import flip_sigma
from schimpy.logging_config import resolve_loglevel, configure_logging
from schimpy.lsc2_v2 import (
    BilinearDensitySizeFunction,
    BoundaryPriors,
    FitParams,
    HysteresisParams,
    PipelineParams,
    ScaledByFieldSizeFunction,
    SigmaCapSizeFunction,
    SigmaPowerSizeFunction,
    SigmaSBlendSizeFunction,
    SigmaTwoZoneSizeFunction,
    fix_sigma_pileups,
    run_pipeline,
)
from schimpy.schism_mesh import read_mesh, write_mesh
from schimpy.schism_polygon import SchismPolygonDictConverter
from schimpy.schism_setup import ensure_outdir
from schimpy.schism_vertical_mesh import SchismLocalVerticalMesh, write_vmesh

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Size function registry
# ---------------------------------------------------------------------------
_SIZEFUN_REGISTRY = {
    "bilinear": lambda p: _make_bilinear(p),
    "sigma_cap": lambda p: SigmaCapSizeFunction(
        **{k: p[k] for k in ("dz_max",) if k in p}
    ),
    "sigma_power": lambda p: SigmaPowerSizeFunction(
        **{k: p[k] for k in ("dz_max", "p") if k in p}
    ),
    "sigma_twozone": lambda p: SigmaTwoZoneSizeFunction(
        **{k: p[k] for k in ("dz_max", "dz_bottom", "K_bottom", "frac_bottom") if k in p}
    ),
    "sigma_sblend": lambda p: SigmaSBlendSizeFunction(
        **{k: p[k] for k in ("dz_max", "theta_b", "theta_f", "alpha", "Hc", "W") if k in p}
    ),
}


def _make_bilinear(params):
    sf = BilinearDensitySizeFunction()
    if params:
        keys = ("a", "b", "c", "d")
        vals = [float(params[k]) if k in params else float(sf.params[i])
                for i, k in enumerate(keys)]
        sf.params = np.array(vals)
    return sf


# ---------------------------------------------------------------------------
# Region constraint builder (works with polygon dicts from YAML include:)
# ---------------------------------------------------------------------------

def _build_node_constraints_from_polygons(mesh, polygons):
    """Build per-node n_min / n_max arrays from a list of SchismPolygon objects.

    Parameters
    ----------
    mesh : SchismMesh
    polygons : list of SchismPolygon
        Each polygon has .type ('min' or 'max') and .attribute (layer count).

    Returns
    -------
    n_min, n_max : ndarray or None
    """
    from shapely.geometry import Point
    from shapely.prepared import prep as shapely_prep

    n = mesh.nodes.shape[0]
    n_min = np.full(n, 2, dtype=np.int32)
    n_max = np.full(n, 99999, dtype=np.int32)
    has_min = False
    has_max = False

    for poly in polygons:
        ptype = str(poly.type).strip().lower()
        # External convention: attribute = number of layers
        # Internal convention: Nlevels = layers + 1
        value = int(float(poly.attribute)) + 1
        if ptype not in ("min", "max"):
            logger.warning(
                "region_constraints: skipping polygon '%s' with type='%s'",
                poly.name, ptype,
            )
            continue
        prepared = shapely_prep(poly)
        for i in range(n):
            pt = Point(float(mesh.nodes[i, 0]), float(mesh.nodes[i, 1]))
            if prepared.contains(pt):
                if ptype == "min":
                    n_min[i] = max(n_min[i], value)
                    has_min = True
                else:
                    n_max[i] = min(n_max[i], value)
                    has_max = True

    conflict = n_min > n_max
    if conflict.any():
        idx = np.where(conflict)[0]
        raise ValueError(
            f"Polygon min > max at {conflict.sum()} nodes "
            f"(first: node {idx[0]}, min={n_min[idx[0]]}, max={n_max[idx[0]]})"
        )

    n_min_out = n_min if has_min else None
    n_max_out = n_max if has_max else None
    n_constrained = int((n_min > 2).sum()) + int((n_max < 99999).sum())
    logger.info(
        "region_constraints: %d polygons, %d node constraints applied",
        len(polygons), n_constrained,
    )
    return n_min_out, n_max_out


def _load_region_constraints(mesh, rc_section):
    """Parse the region_constraints section (dict with 'polygons' key)
    and return (n_min, n_max) arrays.

    The dict comes from YAML after include: expansion, so it contains
    the polygon list inline.
    """
    if rc_section is None:
        return None, None
    if not isinstance(rc_section, dict):
        raise TypeError(
            "region_constraints must be a dict with a 'polygons' key "
            "(use 'include: <file>' to reference an external polygon YAML)."
        )
    polygons = SchismPolygonDictConverter().read(rc_section)
    if not polygons:
        logger.warning("region_constraints section present but contains no polygons")
        return None, None
    return _build_node_constraints_from_polygons(mesh, polygons)


# ---------------------------------------------------------------------------
# Config dict → dataclass helpers
# ---------------------------------------------------------------------------

def _build_pipeline_params(section):
    """Build PipelineParams from the vgrid YAML section dict.

    Only keys present in the YAML override defaults; the dataclass
    is the single source of truth for fallback values.
    """
    alg = section.get("algorithm", {}) or {}
    hyst_d = alg.get("hysteresis", {}) or {}
    fit_d = alg.get("fit", {}) or {}

    # Pass through only keys that are valid dataclass fields.
    # YAML already handles type coercion; the dataclass supplies defaults
    # for anything not specified.
    _hyst_fields = set(HysteresisParams.__dataclass_fields__)
    _fit_fields = set(FitParams.__dataclass_fields__)

    hyst = HysteresisParams(**{k: v for k, v in hyst_d.items() if k in _hyst_fields})
    fit = FitParams(**{k: v for k, v in fit_d.items() if k in _fit_fields})

    pp_kwargs = {}
    if "Lsmooth_passes" in alg:
        pp_kwargs["L_smooth_passes"] = int(alg["Lsmooth_passes"])
    if "Lsmooth_kappa" in alg:
        pp_kwargs["L_smooth_kappa"] = float(alg["Lsmooth_kappa"])
    if "constraint_taper_rings" in section:
        pp_kwargs["constraint_taper_rings"] = int(section["constraint_taper_rings"])

    pp = PipelineParams(hysteresis=hyst, fit=fit, **pp_kwargs)
    return pp


def _build_sizefun(section):
    """Instantiate a VerticalSizeFunction from the depth_function section."""
    df = section.get("depth_function", {}) or {}
    name = str(df.get("name", "bilinear")).lower()
    params = dict(df.get("params", {}) or {})

    try:
        base = _SIZEFUN_REGISTRY[name](params)
    except KeyError:
        raise ValueError(
            f"Unknown depth_function.name='{name}'. "
            f"Valid choices: {sorted(_SIZEFUN_REGISTRY)}"
        )

    # Anchor wiring (e.g. S-blend top anchor)
    anchor = params.get("anchor")
    if anchor:
        if "top" in anchor:
            spec = {"anchor": "top", "value": float(anchor["top"])}
            if hasattr(base, "set_count_spec"):
                base.set_count_spec(spec)
        else:
            raise ValueError("Unsupported anchor; expected {'top': <meters>}.")

    return base


# ---------------------------------------------------------------------------
# Main entry point (dict-based, called from prepare_schism)
# ---------------------------------------------------------------------------

def vgrid_gen_v2(
    hgrid,
    vgrid_out,
    vgrid_version,
    eta,
    *,
    depth_function=None,
    algorithm=None,
    region_constraints=None,
    constraint_taper_rings=3,
    dz_scale_gr3=None,
    debug_prefix=None,
    pileup_log=None,
    **kwargs,
):
    """Generate a SCHISM vgrid using the LSC2 v2 pipeline (BETA).

    This is the v2 counterpart of ``vgrid_gen`` in ``create_vgrid_lsc2.py``.
    It accepts the same top-level parameters (hgrid, vgrid_out, vgrid_version,
    eta) plus v2-specific nested dicts for depth_function, algorithm, and
    region_constraints.

    Parameters
    ----------
    hgrid : str or SchismMesh
        Path to hgrid.gr3 or an already-loaded mesh object.
    vgrid_out : str
        Output vgrid filename.
    vgrid_version : str
        SCHISM version string ('5.8' or '5.10').
    eta : float
        Reference free-surface elevation.
    depth_function : dict, optional
        ``{name: ..., params: {...}}``
    algorithm : dict, optional
        Nested dict with ``hysteresis`` and ``fit`` sub-dicts.
    region_constraints : dict, optional
        Dict with ``polygons`` key (from YAML ``include:`` expansion).
    constraint_taper_rings : int
        Soft-blend rings at constraint polygon edges.
    dz_scale_gr3 : str, optional
        Path to a scalar GR3 for per-node dz scaling.
    debug_prefix : str, optional
        Prefix path for debug GR3 snapshots.
    pileup_log : str, optional
        Path for CSV log of fixed pileup nodes.
    """
    if kwargs:
        unknown = sorted(kwargs.keys())
        logger.warning(
            "vgrid v2 (beta): ignoring unknown keys: %s", ", ".join(unknown)
        )

    if vgrid_version not in ("5.8", "5.10"):
        raise ValueError(f"vgrid_version must be '5.8' or '5.10', got '{vgrid_version}'")

    logger.info(
        "vgrid v2 (beta): generating vertical grid. "
        "Note: v2 is experimental; for the stable legacy generator use vgrid_generator_version: v1"
    )

    # --- Mesh ---
    if hasattr(hgrid, "n_nodes"):
        mesh = hgrid
    else:
        logger.info("Reading mesh: %s", hgrid)
        mesh = read_mesh(hgrid)
    h0 = mesh.nodes[:, 2]
    depth = eta + h0

    # --- Build section dict for helper functions ---
    section = {
        "depth_function": depth_function,
        "algorithm": algorithm,
        "constraint_taper_rings": constraint_taper_rings,
    }

    # --- Size function ---
    sizefun = _build_sizefun(section)

    # Optional per-node dz scaling
    if dz_scale_gr3:
        scale_mesh = read_mesh(dz_scale_gr3)
        scale = scale_mesh.nodes[:, 2].astype(float)
        sizefun = ScaledByFieldSizeFunction(sizefun, scale_field=scale)

    # --- Pipeline params ---
    pp = _build_pipeline_params(section)

    # --- Region constraints ---
    n_min, n_max = _load_region_constraints(mesh, region_constraints)
    pp.n_min = n_min
    pp.n_max = n_max

    # --- Debug outputs ---
    debug = None
    if debug_prefix:
        debug = {
            "Lstar": f"{debug_prefix}Lstar.gr3",
            "Ltilde": f"{debug_prefix}Lstar_smooth.gr3",
            "Nlevels": f"{debug_prefix}nlevels.gr3",
            "Nlayers": f"{debug_prefix}nlayers.gr3",
            "tbottom": f"{debug_prefix}tbottom_target.gr3",
            "bottom_thickness": f"{debug_prefix}bottom_thickness_final.gr3",
            "uniform_sigma": f"{debug_prefix}uniform_sigma.gr3",
            "Nmin": f"{debug_prefix}Nmin.gr3",
            "Nmax": f"{debug_prefix}Nmax.gr3",
        }
        for p in debug.values():
            d = os.path.dirname(p)
            if d:
                os.makedirs(d, exist_ok=True)

    # --- Run pipeline ---
    sigma, Nlevels, h, tmin_arr = run_pipeline(
        mesh, depth, sizefun, eta, pp, debug=debug
    )

    # Post-process: fix pileups
    sigma, Nlevels, pileup_df = fix_sigma_pileups(
        sigma=sigma,
        Nlevels=Nlevels,
        depth=depth,
        tmin=tmin_arr,
        mesh=mesh,
    )

    # Guarantee no all-NaN sigma rows
    for i in range(sigma.shape[0]):
        Ni = int(Nlevels[i])
        if Ni < 2 or not np.isfinite(sigma[i, :Ni]).any():
            sigma[i, :] = np.nan
            sigma[i, 0] = 0.0
            sigma[i, 1] = 1.0
            Nlevels[i] = 2

    if len(pileup_df):
        logger.info("vgrid v2 (beta): fixed %d pileup nodes", len(pileup_df))
        if pileup_log:
            pileup_df.to_csv(pileup_log, index=False)
            logger.info("vgrid v2 (beta): pileup log → %s", pileup_log)
    else:
        logger.info("vgrid v2 (beta): no pileups detected")

    # --- Write output ---
    vmesh = SchismLocalVerticalMesh(flip_sigma(-sigma))
    out_dir = os.path.dirname(vgrid_out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    write_vmesh(vmesh, vgrid_out, vgrid_version)
    logger.info("vgrid v2 (beta): wrote %s", vgrid_out)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

@click.command()
@click.option("--hgrid", required=True, help="Path to hgrid.gr3")
@click.option("--eta", type=float, default=1.5,
              help="Reference free-surface elevation")
@click.option("--vgrid_version", type=str, default="5.10",
              help="SCHISM version for vgrid output ('5.8' or '5.10')")
@click.option("--out_vgrid", type=str, default="vgrid.in",
              help="Output vgrid filename")
@click.option("--config", type=click.Path(exists=True), default=None,
              help="YAML config with depth_function, algorithm, "
              "region_constraints sections")
@click.option("--debug_prefix", type=str, default=None,
              help="Prefix for debug GR3 snapshots")
@click.option("--pileup-log", type=str, default=None,
              help="CSV log of fixed pileup nodes")
@click.option("--logdir", type=click.Path(path_type=Path), default=None,
              help="Directory for log files")
@click.option("--debug", is_flag=True, help="Enable debug-level logging")
@click.option("--quiet", is_flag=True, help="Suppress console logging")
@click.help_option("-h", "--help")
def create_vgrid_lsc2_v2_cli(hgrid, eta, vgrid_version, out_vgrid, config,
                              debug_prefix, pileup_log, logdir, debug, quiet):
    """Create SCHISM vgrid using LSC2 v2 pipeline (BETA).

    For the stable legacy generator, use create_vgrid_lsc2.
    """
    from schimpy.schism_yaml import load as schism_yaml_load

    level, console = resolve_loglevel(debug=debug, quiet=quiet)
    configure_logging(
        package_name="schimpy",
        level=level,
        console=console,
        logdir=logdir,
        logfile_prefix="create_vgrid_lsc2_v2",
    )

    cfg = {}
    if config:
        with open(config, "r") as f:
            cfg = schism_yaml_load(f) or {}

    vgrid_gen_v2(
        hgrid=hgrid,
        vgrid_out=out_vgrid,
        vgrid_version=vgrid_version,
        eta=eta,
        depth_function=cfg.get("depth_function"),
        algorithm=cfg.get("algorithm"),
        region_constraints=cfg.get("region_constraints"),
        constraint_taper_rings=int(cfg.get("constraint_taper_rings", 3)),
        debug_prefix=debug_prefix,
        pileup_log=pileup_log,
    )


if __name__ == "__main__":
    create_vgrid_lsc2_v2_cli()
