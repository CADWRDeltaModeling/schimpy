#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Driver module to prepares input files for a SCHISM run."""
from schimpy import __version__ as schimpy_version
from schimpy.schism_setup import create_schism_setup, check_and_suggest, ensure_outdir
from schimpy.grid_opt import GridOptimizer
from schimpy.stacked_dem_fill import stacked_dem_fill
from schimpy.small_areas import small_areas
from schimpy.split_quad import split_quad
from schimpy import schism_yaml
from schimpy.create_vgrid_lsc2 import vgrid_gen
from schimpy.mesh_volume_tvd import (
    refine_volume_tvd,
    ShorelineOptions,
    FloorOptions,
    TVDOptions,
)
from packaging.version import Version, InvalidVersion
import datetime
import importlib
import numpy as np
import subprocess
import os
import argparse
import logging
import warnings
import shutil


__all__ = ["prepare_schism"]


def create_arg_parser():
    """Create ArgumentParser"""
    parser = argparse.ArgumentParser(description="Prepare SCHISM input files.")
    parser.add_argument(
        dest="main_inputfile", default=None, help="main input file name"
    )
    return parser


def check_min_schimpy_version(inputs) -> None:
    """
    Check the optional 'min_schimpy_version' key in a configuration dict
    against the installed schimpy version.

    Parameters
    ----------
    inputs : dict
        Parsed YAML configuration (top-level dict). If it contains the
        key 'min_schimpy_version', this function validates that the
        installed schimpy version is greater or equal.

    Raises
    ------
    ValueError
        If the installed version is older than the required minimum.
    """
    required = inputs.get("min_schimpy_version")
    if not required:
        return  # nothing to check

    try:
        current = Version(schimpy_version)
        required_v = Version(str(required))
    except InvalidVersion as e:
        raise ValueError(f"Invalid version string in check: {e}") from None

    if current < required_v:
        raise ValueError(
            f"schimpy {required_v} or newer required, but found {current}."
        )


def create_hgrid(s, inputs, logger):
    """Preprocess the hgrid file"""
    section_name = "mesh"
    section = inputs.get(section_name)
    if section is not None:
        split_param = section.get("split_quad")
        if split_param is not None:
            # This should be an in-place
            split_quad(
                s.mesh, inputs["prepro_output_dir"], logger=logger, **split_param
            )
        small_area_param = section.get("small_areas")
        if small_area_param is not None:
            small_area_param["prepro_output_dir"] = inputs["prepro_output_dir"]
            # This just emits warnings unless the fail threshold is met
            small_areas(s.mesh, logger=logger, **small_area_param)

        open_boundaries = section.get("open_boundaries")
        if open_boundaries is not None:
            logger.info("Processing open boundaries...")
            check_and_suggest(open_boundaries, ("linestrings",))
            s.create_open_boundaries(open_boundaries)

        # Fill the missing land and island boundary information
        logger.info("Filling missing land and island boundaries...")
        s.mesh.fill_land_and_island_boundaries()

        # Volumetric optimization of mesh
        option_name = "depth_optimization"

        # if option_name in section.keys():
        opt_params = section.get(option_name)
        default_depth_for_missing_dem = 2.0
        if opt_params is not None:
            method = opt_params.get("method")
            if method is None:
                raise ValueError(
                    "Depth optimization parameters must include method = [volume | volume_tvd]. "
                    "The former is the traditional method which is now deprecated. "
                    "volume_tvd is the new method which requires different parameters which you can "
                    "find listed in a more current template"
                )

            if method == "volume":
                # Legacy GridOptimizer path (deprecated)
                dem_list = section.get("dem_list")
                if dem_list is None:
                    raise ValueError("dem_list must be provided for the mesh optimization")
                expected_items = ("damp", "damp_shoreline", "face_coeff", "volume_coeff")
                check_and_suggest(opt_params, expected_items)
                logger.info("Start optimizing the mesh (legacy 'volume').")
                optimizer = GridOptimizer(
                    mesh=s.mesh,
                    demfiles=dem_list,
                    na_fill=default_depth_for_missing_dem,
                    logger=logger,
                    out_dir=inputs["prepro_output_dir"],
                )
                optimized_elevation = optimizer.optimize(opt_params)
                s.mesh.nodes[:, 2] = np.negative(optimized_elevation)
            elif method=="volume_tvd":
                dem_list = section.get("dem_list")
                if dem_list is None:
                    raise ValueError("For method=volume_tvd, provide dem_list")
                s.dem_list = dem_list
                logger.info(f"Using stacked DEM fill for depth initialization. Length of DEM list: {len(dem_list)}")
                s.mesh.nodes[:, 2] = stacked_dem_fill(
                            dem_list,
                            s.mesh.nodes[:, :2],
                            inputs["prepro_output_dir"],
                            require_all=False,
                            na_fill=default_depth_for_missing_dem,
                            negate=True
                        )
                sl = opt_params.get("shoreline")
                shoreline = None
                if sl is not None:
                    href_spec = sl.get("href", 0.0)
                    # NEW: support dict form to build an href.gr3 via polygons
                    if isinstance(href_spec, dict):
                        href_name = href_spec.get("gr3_name", "href.gr3")   # default name
                        href_path = ensure_outdir(inputs["prepro_output_dir"], href_name)
                        polygons = href_spec.get("polygons", [])
                        default  = href_spec.get("default", None)
                        smooth   = href_spec.get("smooth", None)
                        # reuse the exact machinery used by the 'gr3' section:
                        s.create_node_partitioning(href_path, polygons, default, smooth)  # writes GR3
                        href_arg = href_path
                    else:
                        href_arg = href_spec  # float or string — existing behavior                    
                    
                    logger.debug(f"Using {href_arg} as reference for shoreline discovery")
                    shoreline = ShorelineOptions(
                        href=href_arg,
                        deep_delta=float(sl.get("deep_delta", 1.0)),
                        shore_delta=float(sl.get("shore_delta", 3.0)),
                        seeds=sl.get("seeds"),
                        use_default_seeds=bool(sl.get("use_default_seeds", True)),
                        smooth_relax_factor=float(sl.get("smooth_relax_factor", 0.4)),
                        smooth_strip_max=int(sl.get("smooth_strip_max", 24)),
                        smooth_eps_deg=float(sl.get("smooth_eps_deg", 55.0)),
                        filter_deep=bool(sl.get("filter_deep", True)),
                        epsg=int(sl.get("epsg", 26910)),
                        shore_csv=sl.get("shore_csv"),
                        shore_shp=sl.get("shore_shp"),
                    )

                fl = opt_params.get("floor")
                floor = None
                if fl is not None:
                    floor = FloorOptions(
                        filter=str(fl.get("filter", "max")).lower(),
                        window=int(fl.get("window", 5)),
                        enforce_floor=bool(fl.get("enforce_floor", False)),
                        reg_weight=float(fl.get("reg_weight", 0.0)),
                    )

                tv = opt_params.get("tvd")
                tvd = None
                if tv is not None:
                    tvd = TVDOptions(
                        steps=int(tv.get("steps", 32)),
                        dt=float(tv.get("dt", 1.0)),
                        mu=float(tv.get("mu", 1.0)),
                        lambda_l2=float(tv.get("lambda_l2", 500.0)),
                        tv_weight=float(tv.get("tv_weight", 1.0)),
                        l2_weight=float(tv.get("l2_weight", 5e-5)),
                        cfl_target=float(tv.get("cfl_target", 1.0)),
                        clip_eps=tv.get("clip_eps", None),
                        rel_eta=float(tv.get("rel_eta", 2e-3)),
                        max_backtracks=int(tv.get("max_backtracks", 4)),
                    )

                logger.info("Start optimizing the mesh ('volume_tvd').")
                _ = refine_volume_tvd(
                    s.mesh,
                    dem_spec=dem_list,
                    out_dir=inputs["prepro_output_dir"],
                    shoreline=shoreline,
                    floor=floor,
                    tvd=tvd,
                    cache_dir=os.path.join(inputs["prepro_output_dir"], ".dem_cache"),
                    logger=logger
                )
            else:
                raise ValueError(f"Unknown depth_optimization.method: {method!r}")
        else:
            dem_list = section.get("dem_list")
            s.dem_list = dem_list
            if dem_list is not None:
                s.mesh.nodes[:, 2] = np.negative(
                    stacked_dem_fill(
                        dem_list,
                        s.mesh.nodes[:, :2],
                        inputs["prepro_output_dir"],
                        require_all=False,
                        na_fill=default_depth_for_missing_dem,
                    )
                )
        # This deprecation can be removed at some point
        try:
            ee = section.get("elev_enforcement")
        except:
            ee = None
        if ee is not None:
            raise ValueError("elev_enforcement has been renamed depth_enforcement")

        depth_enforce_params = section.get("depth_enforcement")
        if depth_enforce_params is not None:
            if "polygons" in depth_enforce_params:
                s.mesh.nodes[:, 2] = s.apply_polygons(
                    default=None, polygons=depth_enforce_params["polygons"]
                )
            if "linestrings" in depth_enforce_params:
                s.mesh.nodes[:, 2] = s.apply_linestring_ops(
                    default=None, linestrings=depth_enforce_params["linestrings"]
                )

        # Write hgrid.gr3
        output_dir = inputs["prepro_output_dir"]
        option_name = "gr3_outputfile"
        if option_name in section:
            logger.info("Writing hgrid file...")
            hgrid_out_fpath = ensure_outdir(
                inputs["prepro_output_dir"], section[option_name]
            )
            s.write_hgrid(hgrid_out_fpath, boundary=True)

        # Write hgrid.ll
        option_name = "ll_outputfile"
        if option_name in section:
            logger.info("Creating hgrid.ll file...")
            hgrid_ll_fpath = ensure_outdir(
                inputs["prepro_output_dir"], section[option_name]
            )
            s.write_hgrid_ll(hgrid_ll_fpath, boundary=True)


def create_vgrid(s, inputs, logger):
    section_name = "vgrid"
    section = inputs.get(section_name)
    if section is not None:
        output_dir = inputs["prepro_output_dir"]
        if "hgrid" in section:
            hgrid = section["hgrid"]
        else:
            hgrid = s.mesh
            #
            # msection = inputs.get('mesh')
            # if 'gr3_outputfile' in msection:
            #    logger.info('Using gr3_outputfile from mesh section for generating vgrid')
            #    hgrid = msection['gr3_outputfile']
            # elif 'mesh_inputfile' in msection:
            #    logger.warning('Using mesh_inputfile from mesh section for generating vgrid. This makes sense if you are doing an abbreviated or follow-up preprocessing job')
            #    hgrid  = msection['mesh_inputfile']
        vgrid_out = section["vgrid_out"]
        section["vgrid_out"] = ensure_outdir(inputs["prepro_output_dir"], vgrid_out)
        vgrid_gen(hgrid, **section)


def create_source_sink(s, inputs, logger):
    """Create source_sink.in"""
    dict_ss = inputs.get("sources_sinks")
    if dict_ss is None:
        return
    logger.info("Processing sources/sinks inputs...")
    expected_items = ("sources", "sinks", "outputfile")
    check_and_suggest(dict_ss, expected_items)
    sources = dict_ss.get("sources")
    sources_sinks = {}
    if sources is not None:
        sources_sinks["sources"] = sources
    sinks = dict_ss.get("sinks")
    if sinks is not None:
        sources_sinks["sinks"] = sinks
    fname = dict_ss.get("outputfile")

    if fname is not None:
        output_dir = inputs["prepro_output_dir"]
        fname = ensure_outdir(output_dir, fname)
        logger.info("Creating %s..." % fname)
        s.create_source_sink_in(sources_sinks, fname)


def create_gr3_with_polygons(s, inputs, logger):
    """Create GR3 files with polygons"""
    dict_gr3 = inputs.get("gr3")
    if dict_gr3 is None:
        return
    logger.info("Processing gr3 outputs...")
    expected_items = ("polygons", "default", "smooth")
    output_dir = inputs["prepro_output_dir"]
    for fname, item in dict_gr3.items():
        check_and_suggest(item, expected_items)
        polygons = item.get("polygons", [])
        if polygons is not None:
            polygon_items = (
                "name",
                "vertices",
                "type",
                "attribute",
                "imports",
                "smooth",
            )
            for polygon in polygons:
                check_and_suggest(polygon, polygon_items)
    for fname, item in dict_gr3.items():
        if fname is None:
            logger.warning("No filename is given in one of the gr3 specs")
            continue
        fname = ensure_outdir(inputs["prepro_output_dir"], fname)
        polygons = item.get("polygons", [])
        default = item.get("default")
        logger.info("Creating %s..." % fname)
        smooth = item.get("smooth")
        s.create_node_partitioning(fname, polygons, default, smooth)


def create_prop_with_polygons(s, inputs, logger):
    """Create prop files with polygons"""
    dict_prop = inputs.get("prop")
    if dict_prop is None:
        return
    logger.info("Processing prop outputs...")
    expected_items = ("default", "polygons")
    output_dir = inputs["prepro_output_dir"]
    for fname, item in dict_prop.items():
        check_and_suggest(item, expected_items)
        polygons = item.get("polygons", [])
        if polygons is not None:
            polygon_items = ("name", "vertices", "type", "attribute")
            for polygon in polygons:
                check_and_suggest(polygon, polygon_items)
    for fname, item in dict_prop.items():
        if fname is None:
            logger.warning("No filename is given in one of prop")
            continue
        fname = ensure_outdir(inputs["prepro_output_dir"], fname)
        polygons = item.get("polygons", [])
        default = item.get("default")
        logger.info("Creating %s..." % fname)
        s.create_prop_partitioning(fname, polygons, default)


def create_structures(s, inputs, logger):
    """Create a structure file"""
    dict_struct = inputs.get("hydraulics")
    if dict_struct is None:
        return
    logger.info("Processing structures...")
    expected_items = ("nudging", "structures", "outputfile")
    output_dir = inputs["prepro_output_dir"]
    check_and_suggest(dict_struct, expected_items)
    structures = dict_struct.get("structures")
    if structures is None:
        logger.error("No structures in hydraulics section")
        raise ValueError("No structures in hydraulics section")
    nudging = dict_struct.get("nudging")
    structures = get_structures_from_yaml(structures)
    s.create_structures(structures, nudging)
    fname = dict_struct.get("outputfile")
    fname = ensure_outdir(inputs["prepro_output_dir"], fname)
    if fname is not None:
        fname = os.path.expanduser(fname)
        logger.info("Creating %s..." % fname)
        s.write_structures(fname)


def get_structures_from_yaml(inputs):
    """Get structures from hydraulic_structures.yaml file"""
    structures = inputs.copy()
    structure_items = ("name", "type", "end_points", "configuration", "reference")
    configuration_items = (
        "n_duplicates",
        "elevation",
        "width",
        "height",
        "radius",
        "coefficient",
        "op_downstream",
        "op_upstream",
        "use_time_series",
        "coefficient_height",
        "culvert_n_duplicates",
        "culvert_elevation",
        "culvert_radius",
        "culvert_coefficient",
        "culvert_op_downstream",
        "culvert_op_upstream",
    )
    for structure in structures:
        check_and_suggest(structure, structure_items)
        conf = structure.get("configuration")
        if conf is not None:
            check_and_suggest(conf, configuration_items)

    return structures


def create_fluxflag(s, inputs, logger):
    """Create fluxflag.gr3"""
    dict_flow = inputs.get("flow_outputs")
    if dict_flow is None:
        return
    logger.info("Processing fluxflag outputs...")
    expected_items = ("linestrings", "outputfile")
    output_dir = inputs["prepro_output_dir"]
    check_and_suggest(dict_flow, expected_items)
    flowlines = dict_flow.get("linestrings")
    if flowlines is None:
        raise ValueError("No flowlines in flow_outputs")
    fname = dict_flow.get("outputfile")
    if fname is None:
        logger.info("outputfile not given for flow_outputs. Using fluxflag.prop")
        fname = "fluxflag.prop"
    fname = ensure_outdir(inputs["prepro_output_dir"], fname)
    logger.info("Creating %s..." % fname)
    s.create_flux_regions(flowlines, fname)
    with open(fname, "a") as f:
        for line in flowlines:
            buf = "{}\n".format(line["name"])
            f.write(buf)


def update_spatial_inputs(s, inputs, logger):
    """Create SCHISM grid inputs.

    Parameters
    ----------
    s: SchismSetup
        schism setup object
    inputs: dict
        inputs from an input file
    """
    create_hgrid(s, inputs, logger)
    create_vgrid(s, inputs, logger)
    create_gr3_with_polygons(s, inputs, logger)
    create_source_sink(s, inputs, logger)
    create_prop_with_polygons(s, inputs, logger)
    create_structures(s, inputs, logger)
    create_fluxflag(s, inputs, logger)


def update_temporal_inputs(s, inputs):
    """Create temporal inputs. Under development"""
    # create in interpolated tide file
    output_dir = inputs["prepro_output_dir"]
    sf_tide_out_fpath = os.path.join(output_dir, sf_tide_out_fname)
    s.interpolate_tide(time_start, time_end, dt, sf_tide_in_fpath, sf_tide_out_fpath)
    # Run the FORTRAN code to create elev2D.th
    hgrid_out_fpath = os.path.join(output_dir, hgrid_out_fname)
    webtide_grid_fpath = os.path.join(input_dir, webtide_grid_fname)
    webtide_fpath = os.path.join(input_dir, webtide_fname)
    elev2d_fpath = os.path.join(output_dir, elev2d_fname)
    p = subprocess.Popen(
        [
            "./gen_elev2D_4_NAVD88",
            sf_tide_out_fpath,
            hgrid_out_fpath,
            webtide_grid_fpath,
            webtide_fpath,
            elev2d_fpath,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return_code = p.wait()
    if return_code != 0:
        for l in p.stdout:
            print(l)
        for l in p.stderr:
            print(l)


def check_nested_match(node):
    for key, item in node.items():
        if isinstance(item, dict):
            if key in item.keys():
                node[key] = item[key]
            check_nested_match(item)


def item_exist(inputs, name):
    return True if name in inputs else False


def setup_logger(outdir):
    """Set up a logger"""
    logging_level = logging.INFO
    logging_fname = os.path.join(outdir, "prepare_schism.log")
    logging.basicConfig(level=logging_level, filename=logging_fname, filemode="w")
    console = logging.StreamHandler()
    console.setLevel(logging_level)
    formatter = logging.Formatter("%(message)s")
    console.setFormatter(formatter)
    logging.getLogger("").addHandler(console)


def main():
    """main function"""
    parser = create_arg_parser()
    args = parser.parse_args()
    prepare_schism(args)


def process_output_dir(inputs):
    """Identify or create dir for results and diagnostics and returns the name"""
    if "prepro_output_dir" in inputs:
        outdir = inputs["prepro_output_dir"]
        force = True
    else:
        raise ValueError("Inplace preprocessing not allowed. prepro_output_dir must be specified in launch yaml file. \nConvention is 'prepro_out'")
    if os.path.exists(outdir):
        created = False
    else:
        if force:
            os.mkdir(outdir)
        else:
            raise ValueError("Output directory (output_dir) does not exist")
    return outdir

echo_header = None
def echo_file_header():
    """Return a string with processing timestamp and software versions."""
    global echo_header
    if echo_header is not None: 
        return echo_header
    header_lines = []
    # Timestamp
    now = datetime.datetime.now().isoformat(timespec="seconds")
    header_lines.append(f"# Processed {now}")

    # Version info
    header_lines.append("# Software versions:")

    for pkg in ["bdschism", "schimpy", "vtools"]:
        try:
            mod = importlib.import_module(pkg)
            version = getattr(mod, "__version__", "Not Available")
        except Exception:
            version = "Not Available"
        header_lines.append(f"# {pkg}: {version}")

    echo_header = "\n".join(header_lines) + "\n"
    return echo_header


def process_prepare_yaml(in_fname, use_logging=True):
    """Process the main input YAML file and return the inputs dict without processing any SCHISM inputs."""

    if not os.path.exists(in_fname):
        raise ValueError("Main input file not found")
    with open(in_fname, "r") as f:
        inputs = schism_yaml.load(f)
        outdir = process_output_dir(inputs)

    if use_logging is True:
        setup_logger(outdir)
        logger = logging.getLogger("SCHISM")
    else:
        logger = logging.getLogger("")
    logger.info("Start pre-processing SCHISM inputs...")

    keys_top_level = [
        "config",
        "env",
        "min_schimpy_version",
        "prepro_output_dir",
        "mesh",
        "gr3",
        "vgrid",
        "prop",
        "hydraulics",
        "sources_sinks",
        "flow_outputs",
        "copy_resources",
    ] + schism_yaml.include_keywords
    logger.info("Processing the top level...")
    check_and_suggest(list(inputs.keys()), keys_top_level, logger)

    logger.info("Checking for matching nested keys...")
    check_nested_match(inputs)

    out_fname = os.path.splitext(in_fname)[0] + "_echo" + os.path.splitext(in_fname)[1]

    out_fname = os.path.join(inputs["prepro_output_dir"], out_fname)
    with open(out_fname, "w") as f:
        f.write(echo_file_header())
        f.write(schism_yaml.safe_dump(inputs))

    return inputs, logger


def prepare_schism(args, use_logging=True):
    inputs, logger = process_prepare_yaml(args.main_inputfile, use_logging)

    # Mesh section
    if item_exist(inputs, "mesh"):
        logger.info("Processing mesh section...")
        mesh_items = inputs["mesh"]
        keys_mesh_section = [
            "mesh_inputfile",
            "dem_list",
            "open_boundaries",
            "split_quad",
            "small_areas",
            "depth_optimization",
            "depth_enforcement",
            "gr3_outputfile",
            "ll_outputfile",
        ] + schism_yaml.include_keywords
        check_and_suggest(list(mesh_items.keys()), keys_mesh_section)
        check_min_schimpy_version(inputs)
        if item_exist(inputs["mesh"], "mesh_inputfile"):
            # Read the grid file to be processed
            mesh_input_fpath = os.path.expanduser(mesh_items["mesh_inputfile"])
            s = create_schism_setup(mesh_input_fpath, logger)
            update_spatial_inputs(s, inputs, logger)
        else:
            raise ValueError("No mesh input file in the mesh section.")
    else:
        raise ValueError("No mesh section in the main input.")

    if item_exist(inputs, "copy_resources"):
        logger.info("Copying resources to output dir")
        copy_spec = inputs["copy_resources"]
        for key, item in copy_spec.items():
            outpath = os.path.join(outdir, item)
            if os.path.normpath(key) == os.path.normpath(outpath):
                continue
            logger.info(f"copy_resources: copying {key} to {outpath}")
            shutil.copyfile(key, outpath)

    logger.info(echo_file_header())
    logger.info("Done.")


if __name__ == "__main__":
    main()
