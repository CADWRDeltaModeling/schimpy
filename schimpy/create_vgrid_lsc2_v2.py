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
from schimpy.schism_setup import SchismSetup, ensure_outdir
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
# Region constraint builder (uses SchismSetup.apply_polygons for node selection)
# ---------------------------------------------------------------------------

def _layers_to_levels(value):
    """Convert a user-facing layer count to internal Nlevels."""
    return int(float(value)) + 1


def _parse_range_attribute(attr, *, polygon_name):
    """Parse a two-value layer-count range from YAML/GIS attributes."""
    if isinstance(attr, str):
        s = attr.strip()
        if not s:
            raise ValueError(
                f"region_constraints polygon '{polygon_name}' has type='range' "
                "but attribute is an empty string"
            )
        try:
            import ast
            parsed = ast.literal_eval(s)
        except (ValueError, SyntaxError):
            parsed = [part.strip() for part in s.strip("[]()").split(",")]
        attr = parsed

    if isinstance(attr, np.ndarray):
        vals = attr.tolist()
    elif isinstance(attr, (list, tuple)):
        vals = list(attr)
    else:
        raise ValueError(
            f"region_constraints polygon '{polygon_name}' has type='range' "
            f"but attribute={attr!r} is not a two-value sequence or string"
        )

    if len(vals) != 2:
        raise ValueError(
            f"region_constraints polygon '{polygon_name}' has type='range' "
            f"but attribute={attr!r} parsed to {len(vals)} values; "
            "expected [min_layers, max_layers]"
        )
    return vals


def _constraint_level_bounds_from_polygon_dict(polygon):
    """Return (lo, hi) internal level-count bounds for one vgrid polygon dict.

    User-facing region constraint attributes are layer counts.  The vgrid
    algorithm works in level/interface counts, so conversion happens here and
    only here.
    """
    name = polygon.get("name") or polygon.get("Name")
    ptype = str(polygon.get("type", "none")).strip().lower()
    attr = polygon.get("attribute")

    if attr is None:
        raise ValueError(f"region_constraints polygon {name!r} is missing attribute")

    if ptype == "min":
        lo = _layers_to_levels(attr)
        hi = None
    elif ptype == "max":
        lo = None
        hi = _layers_to_levels(attr)
    elif ptype in ("none", "fixed", "set", ""):
        lo = hi = _layers_to_levels(attr)
    elif ptype == "range":
        vals = _parse_range_attribute(attr, polygon_name=name)
        lo = _layers_to_levels(vals[0])
        hi = _layers_to_levels(vals[1])
    else:
        raise ValueError(
            f"region_constraints polygon {name!r} has unsupported type={ptype!r}. "
            "Expected one of: min, max, none, range."
        )

    if lo is not None and hi is not None and lo > hi:
        raise ValueError(
            f"region_constraints polygon {name!r} has min > max after conversion "
            f"to levels (min={lo}, max={hi}). Constraint attributes are layer "
            "counts in [min_layers, max_layers] order."
        )
    return lo, hi


def _clone_polygon_with_type_and_attribute(polygon, ptype, attribute):
    """Return a plain polygon dict with vgrid bounds converted to levels."""
    out = dict(polygon)
    out["type"] = ptype
    out["attribute"] = attribute
    return out

def _build_node_constraints_from_polygon_dicts(mesh, polygon_dicts):
    """Build vgrid n_min/n_max arrays using SchismSetup.apply_polygons.

    Reuses common polygon infrastructure for node selection, but keeps
    source tracking so min/max conflicts identify the responsible polygons.
    """
    n = mesh.nodes.shape[0]
    min_polygons = []
    max_polygons = []
    min_sources = []
    max_sources = []

    for ipoly, polygon in enumerate(polygon_dicts):
        lo, hi = _constraint_level_bounds_from_polygon_dict(polygon)
        name = polygon.get("name", f"polygon_{ipoly}")
        src = f"{ipoly}: {name}"

        if lo is not None:
            min_polygons.append(
                _clone_polygon_with_type_and_attribute(polygon, "min", lo)
            )
            min_sources.append(src)

        if hi is not None:
            max_polygons.append(
                _clone_polygon_with_type_and_attribute(polygon, "max", hi)
            )
            max_sources.append(src)

    setup = SchismSetup(logger)
    setup.mesh = mesh

    n_min = None
    n_max = None

    min_source = np.full(n, "", dtype=object)
    max_source = np.full(n, "", dtype=object)

    if min_polygons:
        n_min_work = np.full(n, 2, dtype=np.int32)

        for poly, src in zip(min_polygons, min_sources):
            vals = setup.apply_polygons([poly], default=2).astype(np.int32)
            changed = vals > n_min_work
            n_min_work[changed] = vals[changed]
            min_source[changed] = src

        n_min = n_min_work

    if max_polygons:
        n_max_work = np.full(n, 99999, dtype=np.int32)

        for poly, src in zip(max_polygons, max_sources):
            vals = setup.apply_polygons([poly], default=99999).astype(np.int32)
            changed = vals < n_max_work
            n_max_work[changed] = vals[changed]
            max_source[changed] = src

        n_max = n_max_work

    n_min_check = n_min if n_min is not None else np.full(n, 2, dtype=np.int32)
    n_max_check = n_max if n_max is not None else np.full(n, 99999, dtype=np.int32)

    conflict = n_min_check > n_max_check
    if conflict.any():
        idx = np.where(conflict)[0]
        i = idx[0]

        raise ValueError(
            f"Polygon min > max at {conflict.sum()} nodes.\n"
            f"First conflict:\n"
            f"  node        : {i}\n"
            f"  min level   : {n_min_check[i]} from {min_source[i]}\n"
            f"  max level   : {n_max_check[i]} from {max_source[i]}\n"
            f"  min layers  : {n_min_check[i] - 1}\n"
            f"  max layers  : {n_max_check[i] - 1}\n"
            "Constraint attributes in the user interface are layer counts; internal diagnostics "
            "Nmin/Nmax here are levels. Check for overlapping polygons with conflicting min/max layer counts. "
            "Use type='range' with attribute=[min_layers, max_layers] to specify a range of allowed layers."
            "Avoid complex polygons, MULTIPOLYGONs with holes or self-intersections, as they may produce unexpected node selections."
        )

    n_constrained = int((n_min_check > 2).sum()) + int((n_max_check < 99999).sum())
    logger.info(
        "region_constraints: %d polygons, %d node constraints applied",
        len(polygon_dicts), n_constrained,
    )

    return n_min, n_max

def _load_region_constraints(mesh, rc_section):
    """Parse the region_constraints section and return (n_min, n_max) arrays."""
    if rc_section is None:
        return None, None
    if not isinstance(rc_section, dict):
        raise TypeError(
            "region_constraints must be a dict with a 'polygons' key "
            "(use 'include: <file>' to reference an external polygon YAML)."
        )
    polygons = rc_section.get("polygons") or []
    if not polygons:
        logger.warning("region_constraints section present but contains no polygons")
        return None, None
    return _build_node_constraints_from_polygon_dicts(mesh, polygons)


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
    if "Lsmooth_method" in alg:
        pp_kwargs["L_smooth_method"] = str(alg["Lsmooth_method"])
    if "Lsmooth_L_scale" in alg:
        pp_kwargs["L_smooth_L_scale"] = float(alg["Lsmooth_L_scale"])
    if "Lsmooth_depth_scale" in alg:
        pp_kwargs["L_smooth_depth_scale"] = float(alg["Lsmooth_depth_scale"])
    if "Lsmooth_length_power" in alg:
        pp_kwargs["L_smooth_length_power"] = float(alg["Lsmooth_length_power"])

    # Optional projected/Dirichlet-like diffusion for polygon constraints.
    # This is intentionally nested and default-off so legacy configs reproduce.
    cd = alg.get("constraint_diffusion", {}) or {}
    if not isinstance(cd, dict):
        raise TypeError("algorithm.constraint_diffusion must be a mapping if provided")
    if "constraint_diffusion_enable" in alg:
        pp_kwargs["constraint_diffusion_enable"] = bool(alg["constraint_diffusion_enable"])
    if "enabled" in cd:
        pp_kwargs["constraint_diffusion_enable"] = bool(cd["enabled"])
    if "constraint_diffusion_passes" in alg:
        pp_kwargs["constraint_diffusion_passes"] = int(alg["constraint_diffusion_passes"])
    if "passes" in cd:
        pp_kwargs["constraint_diffusion_passes"] = int(cd["passes"])
    if "constraint_diffusion_weight" in alg:
        pp_kwargs["constraint_diffusion_weight"] = float(alg["constraint_diffusion_weight"])
    if "weight" in cd:
        pp_kwargs["constraint_diffusion_weight"] = float(cd["weight"])
    if "method" in cd:
        pp_kwargs["constraint_diffusion_method"] = str(cd["method"])
    if "constraint_diffusion_method" in alg:
        pp_kwargs["constraint_diffusion_method"] = str(alg["constraint_diffusion_method"])
    if "L_scale" in cd:
        pp_kwargs["constraint_diffusion_L_scale"] = float(cd["L_scale"])
    if "depth_scale" in cd:
        pp_kwargs["constraint_diffusion_depth_scale"] = float(cd["depth_scale"])
    if "length_power" in cd:
        pp_kwargs["constraint_diffusion_length_power"] = float(cd["length_power"])
    if "tol" in cd:
        pp_kwargs["constraint_diffusion_tol"] = float(cd["tol"])

    # Integer topological cleanup can be given as flat keys or as a nested
    # integer_cleanup: {enabled, max_component_nodes, max_passes} section.
    ic = alg.get("integer_cleanup", {}) or {}
    if not isinstance(ic, dict):
        raise TypeError("algorithm.integer_cleanup must be a mapping if provided")
    if "integer_cleanup_enable" in alg:
        pp_kwargs["integer_cleanup_enable"] = bool(alg["integer_cleanup_enable"])
    if "enabled" in ic:
        pp_kwargs["integer_cleanup_enable"] = bool(ic["enabled"])
    if "integer_cleanup_max_nodes" in alg:
        pp_kwargs["integer_cleanup_max_nodes"] = int(alg["integer_cleanup_max_nodes"])
    if "max_component_nodes" in ic:
        pp_kwargs["integer_cleanup_max_nodes"] = int(ic["max_component_nodes"])
    if "integer_cleanup_max_passes" in alg:
        pp_kwargs["integer_cleanup_max_passes"] = int(alg["integer_cleanup_max_passes"])
    if "max_passes" in ic:
        pp_kwargs["integer_cleanup_max_passes"] = int(ic["max_passes"])


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
    dz_scale_gr3=None,
    diagnostics=False,
    diagnostics_dir=None,
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
    dz_scale_gr3 : str, optional
        Path to a scalar GR3 for per-node dz scaling.
    diagnostics : bool or dict, optional
        If True, write the standard LSC2 v2 diagnostic GR3 files and pileup CSV.
        A dict may be used as ``{enabled: true, dir: <path>}``.
    diagnostics_dir : str, optional
        Directory for standard diagnostics when diagnostics is True.
    debug_prefix : str, optional
        Legacy prefix path for debug GR3 snapshots. Still supported.
    pileup_log : str, optional
        Legacy path for CSV log of fixed pileup nodes. Still supported.
    """
    if kwargs:
        unknown = sorted(kwargs.keys())
        logger.warning(
            "vgrid v2 (beta): ignoring unknown keys: %s", ", ".join(unknown)
        )

    if vgrid_version not in ("5.8", "5.10"):
        raise ValueError(f"vgrid_version must be '5.8' or '5.10', got '{vgrid_version}'")

    # --- Diagnostics output routing ---
    # New config form is diagnostics: true|false plus diagnostics_dir.  Keep the
    # legacy debug_prefix/pileup_log knobs intact; explicit legacy paths win.
    if isinstance(diagnostics, dict):
        diagnostics_dir = diagnostics.get("dir", diagnostics_dir)
        diagnostics = bool(diagnostics.get("enabled", diagnostics.get("write", True)))
    else:
        diagnostics = bool(diagnostics)

    if diagnostics:
        diag_dir = diagnostics_dir or "vgrid_lsc2_v2_diagnostics"
        os.makedirs(diag_dir, exist_ok=True)
        if debug_prefix is None:
            debug_prefix = os.path.join(diag_dir, "")
        if pileup_log is None:
            pileup_log = os.path.join(diag_dir, "pileups.csv")

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
            "Lconstrained": f"{debug_prefix}Lstar_constrained.gr3",
            "Nlevels_raw": f"{debug_prefix}nlevels_raw.gr3",
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
@click.option("--diagnostics/--no-diagnostics", default=None,
              help="Write standard diagnostic GR3 files and pileup CSV")
@click.option("--diagnostics-dir", type=str, default=None,
              help="Directory for standard diagnostics")
@click.option("--logdir", type=click.Path(path_type=Path), default=None,
              help="Directory for log files")
@click.option("--debug", is_flag=True, help="Enable debug-level logging")
@click.option("--quiet", is_flag=True, help="Suppress console logging")
@click.help_option("-h", "--help")
def create_vgrid_lsc2_v2_cli(hgrid, eta, vgrid_version, out_vgrid, config,
                              debug_prefix, pileup_log, diagnostics, diagnostics_dir,
                              logdir, debug, quiet):
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

    cfg_diagnostics = cfg.get("diagnostics", False)
    if diagnostics is not None:
        cfg_diagnostics = diagnostics

    vgrid_gen_v2(
        hgrid=hgrid,
        vgrid_out=out_vgrid,
        vgrid_version=vgrid_version,
        eta=eta,
        depth_function=cfg.get("depth_function"),
        algorithm=cfg.get("algorithm"),
        region_constraints=cfg.get("region_constraints"),
        diagnostics=cfg_diagnostics,
        diagnostics_dir=diagnostics_dir or cfg.get("diagnostics_dir"),
        debug_prefix=debug_prefix,
        pileup_log=pileup_log,
    )


if __name__ == "__main__":
    create_vgrid_lsc2_v2_cli()
