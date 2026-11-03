# -*- coding: utf-8 -*-
"""Generate SCHISM ``bctides.in`` files from a YAML boundary description.

All elevation, velocity, temperature, salinity and tracer boundary types read
by SCHISM (``schism_init.F90``) are supported. Each open boundary lists its
``variables``; each variable names a ``source`` keyword that maps one-to-one
onto a SCHISM type code. A numeric ``source`` (or a list of numbers for a
tracer module) means a constant (code 2). Omitting a variable means code 0.

Source keywords
---------------
elevation (``iettype``)
    ``th_text`` (1, ``elev.th``), number (2), ``tidal`` (3),
    ``th_nc`` (4, ``elev2D.th.nc``), ``tidal_th_nc`` (5, tides added to
    ``elev2D.th.nc``).
velocity (``ifltype``)
    ``th_text`` (1, ``flux.th``), number (2, discharge, negative is inflow),
    ``tidal`` (3), ``th_nc`` (4, ``uv3D.th.nc``), ``tidal_th_nc`` (5),
    ``th_nc_relaxed`` (-4), ``tidal_th_nc_relaxed`` (-5), ``flather_1`` (-1),
    ``discharge_relation`` (-2).
temperature, salinity (``itetype``, ``isatype``) and tracers (``itrtype``)
    ``th_text`` (1), number or list of numbers (2), ``initial_profile`` (3),
    ``th_nc`` (4, ``TEM_3D.th.nc``, ``SAL_3D.th.nc`` or ``[MOD]_3D.th.nc``).
    Each requires ``nudge`` in [0, 1].

See the "Boundary conditions (bctides)" chapter of the schimpy documentation
for the complete YAML reference and examples.
"""
from schimpy.schism_yaml import load
from schimpy.schism_mesh import SchismMeshGr3Reader, BoundaryType
from schimpy.logging_config import configure_logging, resolve_loglevel
import datetime
import logging
import numbers
from pathlib import Path
import numpy as np
import click


__all__ = ["bctide_cli"]

logger = logging.getLogger(__name__)

#: SCHISM type code used for a numeric (constant) ``source``.
CONSTANT_CODE = 2

#: Elevation ``source`` keywords and their SCHISM ``iettype`` codes.
ELEVATION_SOURCES = {
    "th_text": 1,
    "tidal": 3,
    "th_nc": 4,
    "tidal_th_nc": 5,
}

#: Velocity ``source`` keywords and their SCHISM ``ifltype`` codes.
VELOCITY_SOURCES = {
    "th_text": 1,
    "tidal": 3,
    "th_nc": 4,
    "tidal_th_nc": 5,
    "th_nc_relaxed": -4,
    "tidal_th_nc_relaxed": -5,
    "flather_1": -1,
    "discharge_relation": -2,
}

#: Temperature, salinity and tracer ``source`` keywords and their SCHISM codes.
SCALAR_SOURCES = {
    "th_text": 1,
    "initial_profile": 3,
    "th_nc": 4,
}


@click.command()
@click.argument("main_yaml", type=click.Path(exists=True))
@click.option("--logdir", type=click.Path(path_type=Path), default=None,
              help="Directory for a log file.")
@click.option("--debug", is_flag=True, help="Enable debug logging.")
@click.option("--quiet", is_flag=True, help="Suppress console logging.")
@click.help_option("-h", "--help")
def bctide_cli(main_yaml, logdir, debug, quiet):
    """Create SCHISM bctides.in files from MAIN_YAML.

    MAIN_YAML names the mesh and lists one entry per bctides file to write:

    \b
        mesh:
          mesh_inputfile: ./hgrid.gr3
          open_boundaries: !include open_boundary.yaml   # optional
        bctides:
          include:
            - bctides.in.2d.yaml
            - bctides.in.3d.yaml

    Each bctides entry is written to the file given by its ``name``. The
    boundary assignment chosen for every open boundary is logged.

    \b
    Example:
      bctide main_bctides.yaml
    """
    level, console = resolve_loglevel(debug=debug, quiet=quiet)
    configure_logging(
        package_name="schimpy",
        level=level,
        console=console,
        logdir=logdir,
        logfile_prefix="bctide",
    )
    with open(main_yaml, "r") as fn:
        main_yaml_dict = load(fn)

    mesh_in = main_yaml_dict["mesh"]["mesh_inputfile"]
    sr = SchismMeshGr3Reader()

    hgrid = sr.read(mesh_in)
    open_boundary_segments = None
    if "open_boundaries" in main_yaml_dict["mesh"].keys():
        open_boundary_segments = main_yaml_dict["mesh"]["open_boundaries"]
    bctides_dic = main_yaml_dict.get("bctides")
    for bctides_yaml in bctides_dic:
        out = bctides_yaml["name"]
        by = load_boundary(hgrid, bctides_yaml, open_boundary_segments)
        by.write_bctides(out)

def update_mesh_open_boundaries(mesh, segments):
    """Overtide mesh open boundaries with linestrings from segments dict
       borrowed from schism_setup.py
       segments: dict containing "linestrings" key
    """
    mesh.clear_boundaries()
    boundary_only = True
    linestrings = segments.get("linestrings")
    if linestrings is None:
        raise ValueError("Linestrings is required for open boundaries specification")
    for item in linestrings:
        name = item.get("name")
        if name is None:
            name = item.get("Name")
        p = np.array(item["coordinates"])
        start = mesh.find_closest_nodes(p[0], 1, boundary_only)
        end = mesh.find_closest_nodes(p[1], 1, boundary_only)
        path = mesh._build_boundary_node_string(start, end, boundary_only)
        comment = '! Open boundary "%s"' % name
        mesh.add_boundary(path, BoundaryType.OPEN, comment)
    mesh.fill_land_and_island_boundaries()




def load_boundary(hgrid,fn,open_boundary_segments):
    return boundary(hgrid,fn,open_boundary_segments)


def _source_code(variable, spec, table, label):
    """Resolve the ``source`` of a boundary variable to a SCHISM type code.

    Parameters
    ----------
    variable : str
        Variable name used in error messages, e.g. ``"velocity"``.
    spec : dict
        Variable specification containing a ``source`` entry.
    table : dict
        Mapping of accepted ``source`` keywords to SCHISM type codes.
    label : str
        Boundary name used in error messages.

    Returns
    -------
    code : int
        SCHISM type code. A numeric source gives :data:`CONSTANT_CODE`.
    keyword : str
        The keyword, or ``"constant"`` for a numeric source.

    Raises
    ------
    ValueError
        If ``source`` is missing or is not a number or a keyword in `table`.
    """
    if "source" not in spec:
        raise ValueError(f"Boundary '{label}': {variable} is missing 'source'.")
    source = spec["source"]
    if isinstance(source, numbers.Number) and not isinstance(source, bool):
        return CONSTANT_CODE, "constant"
    if isinstance(source, str) and source in table:
        return table[source], source
    raise ValueError(
        f"Boundary '{label}': {variable} source {source!r} is not supported. "
        f"Use a number for a constant or one of: {', '.join(table)}."
    )


def _unit_interval(spec, key, what, label):
    """Return required factor ``spec[key]`` after checking it lies in [0, 1].

    Parameters
    ----------
    spec : dict
        Specification holding the factor.
    key : str
        Key of the factor, e.g. ``"nudge"`` or ``"inflow_relax"``.
    what : str
        Description used in error messages.
    label : str
        Boundary name used in error messages.

    Returns
    -------
    float or int
        The factor as given in the YAML.

    Raises
    ------
    ValueError
        If the factor is missing, not numeric or outside [0, 1].
    """
    if key not in spec:
        raise ValueError(f"Boundary '{label}': {what} requires '{key}'.")
    value = spec[key]
    if (not isinstance(value, numbers.Number) or isinstance(value, bool)
            or not 0.0 <= value <= 1.0):
        raise ValueError(
            f"Boundary '{label}': {what} '{key}' must be a number in [0, 1], got {value!r}."
        )
    return value


class boundary(object):
    """
    A class to generate boundary condition input file for SCHISM
    """

    def __init__(self, hgrid,bc_yaml=None,boundary_segments=None):
        """ initialize the boundary condition from yaml file"""

        #main_id = "bctides"
        # "date" is optional, default to Jan 1, 2020
        # ensure bc_yaml is a dict and normalize/insert a datetime object so the
        # later `self.date = bc_yaml["date"]` works reliably
        if bc_yaml is None:
            bc_yaml = {}
        if "date" not in bc_yaml or bc_yaml["date"] is None:
            bc_yaml["date"] = datetime.datetime(2020, 1, 1)
        else:
            d = bc_yaml["date"]
            if isinstance(d, str):
            # try ISO format first, then common YYYY-MM-DD
                try:
                    bc_yaml["date"] = datetime.datetime.fromisoformat(d)
                except Exception:
                    try:
                       bc_yaml["date"] = datetime.datetime.strptime(d, "%Y-%m-%d")
                    except Exception:
                       raise ValueError(f"Unrecognized date format for bctides date: {d}")
            elif isinstance(d, datetime.date) and not isinstance(d, datetime.datetime):
                bc_yaml["date"] = datetime.datetime(d.year, d.month, d.day)
            else:
                raise ValueError(f"Unrecognized date format for bctides date: {d}")
        self.date = bc_yaml["date"]
        if "mode" not in bc_yaml or bc_yaml["mode"] is None: ## raise error if mode is not given, since it is required for determining the boundary condition type
            raise ValueError(" 'mode' is required but not specified in bctides YAML. Please specify 'mode' as either 'barotropic' or 'baroclinic'.")
        self.mode = bc_yaml.get("mode")
        if self.mode not in ["barotropic", "baroclinic"]:
            raise ValueError(f"Invalid mode specified in bctides YAML: {self.mode}. Expected 'barotropic' or 'baroclinic'.")
        self.computation_temperature_salinity = True
        if self.mode == "barotropic":
            if not("temperature_salinity_computation" in bc_yaml.keys()):
                raise ValueError("temperature_salinity_computation is required for barotropic mode but not specified in bctides YAML. Please specify 'temperature_salinity_computation' as true or false.")
            self.computation_temperature_salinity = bc_yaml.get("temperature_salinity_computation")

        if "earth_tides" in bc_yaml.keys():
            self.earth_tides = bc_yaml["earth_tides"]
        else:
            self.earth_tides = None
        if "boundary_forcing_tides" in bc_yaml.keys():
            self.boundary_tides = bc_yaml["boundary_forcing_tides"]
        else:
            self.boundary_tides = None
        
        self.open_boundaries = bc_yaml["open_boundaries"]
        self.hgrid = hgrid  

        ## if boundary grid is provided, use it to override the hgrid boundary 
        ## grid
        self.boundary_segments = boundary_segments
        if self.boundary_segments is not None:
            ## also check if the number of open boundaries and names match
            ## with bctides open boundaries
            num_boundary_segments = len(self.boundary_segments["linestrings"])
            if num_boundary_segments != len(self.open_boundaries):
                raise ValueError(
                    "number of boundary segments from included open_boundary yaml does not match number of open boundaries in bctides YAML: %s vs %s"
                    % (num_boundary_segments, len(self.open_boundaries))
                )
            for i in range(num_boundary_segments):
                bname_yaml = self.open_boundaries[i]["name"]
                bname_segment = self.boundary_segments["linestrings"][i]["name"]
                if bname_yaml != bname_segment:
                    raise ValueError(
                        "boundary name mismatch between included open_boundary YAML and bctides YAML: %s vs %s"
                        % (bname_yaml, bname_segment)
                    )       
            update_mesh_open_boundaries(self.hgrid, self.boundary_segments)

        ## tracer order
        self.tracer_lst = {
            "GEN": 1,
            "AGE": 2,
            "SED": 3,
            "ECO": 4,
            "ICM": 5,
            "COSINE": 6,
            "FECO": 7,
            "TIMOR": 8,
            "FABM": 9,
            "DVD": 10,
        }

##   also read in tracer modules defined in the bctides yaml here
        self.tracer_modules = bc_yaml.get("modules", [])
##   find out number of tracer modules and set up tracer order list
##   accroding to the self.tracer_lst defined above, if module not found
##   module section in yaml file looks like:
##    modules:
##      - name: GEN
##        num_tracers: 1
##      - name: SED
##        num_tracers: 5

        self.tracer_mod_pos = {}
        self.tracer_mod_num = {}
        tracer_mods = []
        for mod in self.tracer_modules:
            mod_id = mod.get("name")
            num_tracers = mod.get("num_tracers", 0) 
            if mod_id in self.tracer_lst.keys():
                tracer_mods.append(mod_id)
                self.tracer_mod_num[mod_id] = num_tracers
            else:
                raise ValueError("tracer module %s not recognized" % mod) 
        tracer_mods.sort(key=lambda d: self.tracer_lst[d])
        for kk in range(len(tracer_mods)):
            self.tracer_mod_pos[tracer_mods[kk]] = kk  
            

    def _norm_earth(self,consts):
        out = []
        for item in consts:
            if not isinstance(item, dict):
                # unsupported format -> skip
                continue

            # case: explicit dict with "name" key
            if "name" in item:
                name = item.get("name")
                amp = item.get("amplitude", item.get("amp", 0.0))
                node = item.get("node_factor", item.get("node", 1.0))
                freq = item.get("angular_frequency", item.get("frequency", 0.0))
                eqa = item.get(
                    "equilibrium_argument",
                    item.get("eqa", item.get("phase", 0.0)),
                )
                out.append(
                    {
                        "name": name,
                        "amplitude": amp,
                        "node_factor": node,
                        "angular_frequency": freq,
                        "equilibrium_argument": eqa,
                    }
                )
                continue

            # case: mapping like {K1: {...}} or {K1: [..]}
            key = next(iter(item))
            vals = item[key]
            name = key

            if isinstance(vals, dict):
                amp = vals.get("amplitude", vals.get("amp", 0.0))
                node = vals.get("node_factor", vals.get("node", 1.0))
                freq = vals.get("angular_frequency", vals.get("frequency", 0.0))
                eqa = vals.get(
                    "equilibrium_argument",
                    vals.get("eqa", vals.get("phase", 0.0)),
                )
                out.append(
                    {
                        "name": name,
                        "amplitude": amp,
                        "node_factor": node,
                        "angular_frequency": freq,
                        "equilibrium_argument": eqa,
                    }
                )
                continue

            if isinstance(vals, (list, tuple)):
                # support common list formats used historically
                if len(vals) == 5:
                    # example: [?, amplitude, freq, node_factor, eqa]
                    amp = vals[1]
                    freq = vals[2]
                    node = vals[3]
                    eqa = vals[4]
                elif len(vals) == 4:
                    # example: [amp, freq, node, eqa]
                    amp = vals[0]
                    freq = vals[1]
                    node = vals[2]
                    eqa = vals[3]
                else:
                    amp = vals[0] if len(vals) > 0 else 0.0
                    freq = vals[2] if len(vals) > 2 else 0.0
                    node = vals[3] if len(vals) > 3 else 1.0
                    eqa = vals[-1] if len(vals) > 0 else 0.0
                out.append(
                    {
                        "name": name,
                        "amplitude": amp,
                        "node_factor": node,
                        "angular_frequency": freq,
                        "equilibrium_argument": eqa,
                    }
                )
                continue

            # unsupported nested format -> skip
            continue

        return out


    def _norm_boundary(self,consts):
        """
        Normalize boundary (forcing) tidal constituent specifications into list of dicts with keys:
          name, angular_frequency, node_factor, earth_equilibrium_argument

        Supported input shapes (from YAML examples):
          - explicit dict with "name": {...}
          - mapping like {O1: {...}} or {O1: [freq, node, eqa]}
          - list/tuple forms [freq, node, eqa] when provided as the value for a key
          - simple string -> defaults
        """
        out = []
        for item in consts or []:
            try:
                # explicit dict with 'name'
                if isinstance(item, dict) and "name" in item:
                    name = item.get("name")
                    freq = item.get("angular_frequency", item.get("frequency", 0.0))
                    node = item.get("node_factor", item.get("node", 1.0))
                    eqa = item.get(
                        "equilibrium_argument", item.get("eqa", item.get("phase", 0.0))
                    )
                    out.append(
                        {
                            "name": name,
                            "angular_frequency": freq,
                            "node_factor": node,
                            "equilibrium_argument": eqa,
                        }
                    )
                    continue

                # mapping {NAME: {...}} or {NAME: [...]}
                if isinstance(item, dict):
                    key = next(iter(item))
                    vals = item[key]
                    name = key

                    if isinstance(vals, dict):
                        freq = vals.get("angular_frequency", vals.get("frequency", 0.0))
                        node = vals.get("node_factor", vals.get("node", 1.0))
                        eqa = vals.get(
                            "equilibrium_argument", vals.get("eqa", vals.get("phase", 0.0))
                        )
                        out.append(
                            {
                                "name": name,
                                "angular_frequency": freq,
                                "node_factor": node,
                                "equilibrium_argument": eqa,
                            }
                        )
                        continue

                    if isinstance(vals, (list, tuple)):
                        if len(vals) >= 3:
                            freq, node, eqa = vals[0], vals[1], vals[2]
                        else:
                            # partial lists: fill with defaults
                            freq = vals[0] if len(vals) > 0 else 0.0
                            node = vals[1] if len(vals) > 1 else 1.0
                            eqa = vals[-1] if len(vals) > 0 else 0.0
                        out.append(
                            {
                                "name": name,
                                "angular_frequency": freq,
                                "node_factor": node,
                                "equilibrium_argument": eqa,
                            }
                        )
                        continue

                    # vals is scalar (frequency) or unsupported -> try to use as frequency
                    if isinstance(vals, (int, float)):
                        out.append(
                            {
                                "name": name,
                                "angular_frequency": vals,
                                "node_factor": 1.0,
                                "equilibrium_argument": 0.0,
                            }
                        )
                        continue

                    # otherwise skip
                    continue

                # plain string -> treat as name with defaults
                if isinstance(item, str):
                    out.append(
                        {
                            "name": item,
                            "angular_frequency": 0.0,
                            "node_factor": 1.0,
                            "equilibrium_argument": 0.0,
                        }
                    )
                    continue

            except Exception:
                # skip malformed entries
                continue

        return out
    
    def _handle_tidal_constituents_val(self,vals,consts_name,const_type_str,boundary_id,node_id):
        """
        Determine and return the numeric amplitude (value) for a tidal constituent given a
        variety of accepted input formats.

        The function accepts a scalar numeric value, an indexable sequence (list/tuple/np.ndarray)
        from which a per-node value will be selected, or a string expression which will be
        evaluated to produce a numeric result.

        Args:
            vals: add argument meaning comment here
                The input value(s) for the tidal constituent. Accepted types:
                - numbers.Number: treated as a scalar amplitude and returned as-is.
                - list, tuple, np.ndarray: treated as a per-node collection; the element at
                  index `node_id` is returned. A mismatched length will raise ValueError.
                - str: treated as a Python expression and evaluated with eval(); the
                  evaluated result is returned if successful.
            consts_name: 
                Human-readable name or identifier of the tidal constituent (used in error
                messages).
            const_type_str: a
                Short description of the constituent paramter (such as elev_phase
                elev_amp,u_amp,u_phase) type (used in error messages).
            boundary_id: 
                Identifier for the boundary to which the constituent applies (used in
                error messages).
            node_id: 
                Integer index of the node for which a per-node value should be selected
                when `vals` is an indexable collection.

        Returns:
            float: The resolved amplitude value for the constituent at the requested node.

        Raises:
            ValueError: If
              - `vals` is an indexable collection but indexing with `node_id` fails
                (e.g., length mismatch).
              - `vals` is a string but its evaluation via eval() fails.
              - `vals` is of an unsupported type.

        Notes:
            - The function uses eval() to evaluate string expressions. For safety, ensure
              inputs are trusted or consider replacing eval() with a safer parser if
              untrusted input must be supported.
            - The returned value is a float-compatible numeric; callers may want to cast
              or validate the result further if required by downstream code.
        """
        ##
        # amplitude value determination
        val = 0.0
        if isinstance(vals, numbers.Number):
            val = vals
        elif isinstance(vals, (list, tuple, np.ndarray)):
            try:
                val = vals[node_id]
            except Exception:
                raise ValueError(
                    f"{const_type_str} list length does not match number of nodes for \
                    tidal constituent {consts_name} at boundary {boundary_id} \n"
                )
        elif isinstance(vals, str):
            try:
                val = eval(vals)
            except Exception:
                raise ValueError(f"{val} is not a supported tidal expressionfor tidal constituent \
                                 {const_type_str} for tidal constituent {consts_name} at boundary \
                                    {boundary_id} \n "
                )
        else:
            raise ValueError(f"Unsupported tidal type {vals} for tidal constituent {const_type_str} for\
                              tidal constituent {consts_name} at boundary {boundary_id} \n")
        return val

    def check_inflow_boundary_spec(self,var_dic, bound_name,tracer_name):

        if self.mode == "barotropic":
            if not self.computation_temperature_salinity:
                if var_dic is not None:
                    raise ValueError(
                        f"Boundary '{bound_name}' has specified '{tracer_name}' boundary conditions, but 'temperature_salinity_computation' \
                        is set to false in barotropic mode. Please remove '{tracer_name}' specification or set 'temperature_salinity_computation' \
                        to true."
                    )
            else:
                 if var_dic is None and (tracer_name in ["temperature", "salinity"]):
                    raise ValueError(
                        f"Boundary '{bound_name}' is missing required '{tracer_name}' specification for inflow/bidirectional boundary in \
                          barotropic mode with 'temperature_salinity_computation' enabled. Please specify '{tracer_name}' boundary conditions."
                    )
                 else:
                    if var_dic is not None and not (tracer_name in ["temperature", "salinity"]):
                        raise ValueError(
                            f"Boundary '{bound_name}' has specified '{tracer_name}' boundary conditions in barotropic mode. \
                              Please remove '{tracer_name}' specification" )
        
        else: # baroclinic mode,         
            if var_dic is None:
                raise ValueError(
                    f"boundary '{bound_name}' is missing variable specification. Please specify {tracer_name} for inflow/bidirectional boundary."
                )

    _VARIABLES = ("elevation", "velocity", "temperature", "salinity", "tracers")

    def _resolve_open_boundary(self, index, spec, num_nodes):
        """Validate one open boundary and resolve its SCHISM type codes.

        The checks mirror those SCHISM makes when reading ``bctides.in`` so that
        problems surface during preprocessing rather than at model start.

        Parameters
        ----------
        index : int
            Zero-based position of the boundary in hgrid order.
        spec : dict
            The open boundary entry from the bctides YAML.
        num_nodes : int
            Number of nodes on the boundary.

        Returns
        -------
        dict
            Resolved assignment with keys ``label``, ``flow_direction``,
            ``num_nodes``, one entry per variable in ``elevation``,
            ``velocity``, ``temperature`` and ``salinity`` (each a dict with
            ``code``, ``keyword``, ``spec`` and, for scalars, ``nudge``) and
            ``tracers`` (a list ordered by SCHISM module order whose items are
            None or dicts with ``module``, ``code``, ``keyword``, ``values``
            and ``nudge``).

        Raises
        ------
        ValueError
            If the specification is incomplete, inconsistent or uses an
            unsupported keyword.
        """
        bound_name = spec.get("name", "")
        label = bound_name if bound_name else f"boundary_{index}"

        if "flow_direction" not in spec:
            raise ValueError(
                f"Boundary '{label}' is missing required 'flow_direction'. "
                "Use 'inflow', 'outflow' or 'bidirection'."
            )
        flow_direction = spec["flow_direction"]

        variables = {}
        for entry in spec.get("variables") or []:
            if not isinstance(entry, dict) or len(entry) != 1:
                raise ValueError(
                    f"Boundary '{label}': each item under 'variables' must be a single "
                    f"variable such as 'velocity: {{source: th_text}}'; got {entry!r}. "
                    "Check the YAML indentation."
                )
            ((vname, vspec),) = entry.items()
            if vname not in self._VARIABLES:
                raise ValueError(
                    f"Boundary '{label}': unknown variable '{vname}'. "
                    f"Expected one of: {', '.join(self._VARIABLES)}."
                )
            if vspec is None:
                raise ValueError(
                    f"Boundary '{label}': '{vname}' has no specification. "
                    "Check the YAML indentation."
                )
            if vname in variables:
                raise ValueError(f"Boundary '{label}': '{vname}' is specified more than once.")
            variables[vname] = vspec

        if flow_direction == "outflow":
            for vname in ("temperature", "salinity", "tracers"):
                if vname in variables:
                    raise ValueError(
                        f"Boundary '{label}' is configured as an outflow; a {vname} "
                        "boundary specification is not allowed."
                    )
        elif flow_direction in ("inflow", "bidirection"):
            check_list = ["temperature", "salinity"]
            if len(self.tracer_modules) > 0:
                check_list.append("tracers")
            elif "tracers" in variables:
                raise ValueError(
                    f"Boundary '{label}' has specified 'tracers' boundary conditions, but no "
                    "tracer modules are defined in the bctides YAML. Please add tracer modules "
                    "to the top-level 'modules' section or remove 'tracers' from the boundary."
                )
            for var in check_list:
                self.check_inflow_boundary_spec(variables.get(var), label, var)
        else:
            raise ValueError(
                f"Boundary '{label}' has invalid 'flow_direction' value: {flow_direction}. "
                "Expected 'inflow', 'outflow', or 'bidirection'."
            )

        unset = {"code": 0, "keyword": None, "spec": None}
        bc = {
            "label": label,
            "flow_direction": flow_direction,
            "num_nodes": num_nodes,
        }

        elev = variables.get("elevation")
        bc["elevation"] = dict(unset)
        if elev is not None:
            code, keyword = _source_code("elevation", elev, ELEVATION_SOURCES, label)
            if code in (3, 5) and "constituents" not in elev:
                raise ValueError(
                    f"Boundary '{label}': elevation source '{keyword}' requires 'constituents'."
                )
            bc["elevation"] = {"code": code, "keyword": keyword, "spec": elev}
        elev_id = bc["elevation"]["code"]

        vel = variables.get("velocity")
        bc["velocity"] = dict(unset)
        if vel is not None:
            code, keyword = _source_code("velocity", vel, VELOCITY_SOURCES, label)
            if code in (3, 5, -5) and "constituents" not in vel:
                raise ValueError(
                    f"Boundary '{label}': velocity source '{keyword}' requires 'constituents'."
                )
            if code in (-1, -2) and elev_id != 0:
                raise ValueError(
                    f"Boundary '{label}': velocity source '{keyword}' requires that no "
                    "elevation boundary be specified."
                )
            if code in (-4, -5):
                if elev_id == 0:
                    raise ValueError(
                        f"Boundary '{label}': velocity source '{keyword}' requires an "
                        "elevation boundary."
                    )
                for key in ("inflow_relax", "outflow_relax"):
                    _unit_interval(vel, key, f"velocity source '{keyword}'", label)
            if code == -2:
                coefs = vel.get("coefficients")
                if (not isinstance(coefs, list) or len(coefs) != 4
                        or not all(isinstance(c, numbers.Number) for c in coefs)):
                    raise ValueError(
                        f"Boundary '{label}': velocity source 'discharge_relation' requires "
                        f"'coefficients' as a list of 4 numbers, got {coefs!r}."
                    )
            if code == -1:
                for key in ("eta_m0", "vnormal_nodes"):
                    vals = vel.get(key)
                    if not isinstance(vals, list) or len(vals) != num_nodes:
                        raise ValueError(
                            f"Boundary '{label}': velocity source 'flather_1' requires "
                            f"'{key}' as a list with one entry per boundary node ({num_nodes})."
                        )
            bc["velocity"] = {"code": code, "keyword": keyword, "spec": vel}

        for vname in ("temperature", "salinity"):
            vspec = variables.get(vname)
            bc[vname] = dict(unset)
            if vspec is not None:
                code, keyword = _source_code(vname, vspec, SCALAR_SOURCES, label)
                nudge = _unit_interval(vspec, "nudge", vname, label)
                bc[vname] = {"code": code, "keyword": keyword, "spec": vspec, "nudge": nudge}

        tracer_specs = [None] * len(self.tracer_mod_pos)
        tracers = variables.get("tracers")
        if tracers is not None:
            if not isinstance(tracers, list):
                raise ValueError(
                    f"Boundary '{label}': 'tracers' must be a list of module entries."
                )
            for tracer in tracers:
                tracer_mod = tracer.get("module")
                if tracer_mod not in self.tracer_mod_pos:
                    raise ValueError(
                        f"Tracer module '{tracer_mod}' is not defined in the bctides 'modules' section. "
                        f"Please add it to the top-level modules list. "
                        f"Valid modules: {', '.join(sorted(self.tracer_mod_pos.keys()))}"
                    )
                pos = self.tracer_mod_pos[tracer_mod]
                if tracer_specs[pos] is not None:
                    raise ValueError(
                        f"Boundary '{label}': tracer module '{tracer_mod}' is specified more than once."
                    )
                if "relax" in tracer:
                    raise ValueError(
                        f"Boundary '{label}': tracer module '{tracer_mod}' uses 'relax'; "
                        "the tracer factor is named 'nudge'."
                    )
                what = f"tracer module '{tracer_mod}'"
                source = tracer.get("source")
                values = None
                if isinstance(source, (list, numbers.Number)) and not isinstance(source, bool):
                    values = source if isinstance(source, list) else [source]
                    if not all(isinstance(x, numbers.Number) for x in values):
                        raise ValueError(
                            f"Boundary '{label}': {what} constant source must contain only "
                            f"numbers, got {source!r}."
                        )
                    if len(values) != self.tracer_mod_num[tracer_mod]:
                        raise ValueError(
                            f"Boundary '{label}': {what} boundary source length {len(values)} "
                            "does not match number of tracers defined in module section: "
                            f"{self.tracer_mod_num[tracer_mod]}"
                        )
                    code, keyword = CONSTANT_CODE, "constant"
                else:
                    code, keyword = _source_code(what, tracer, SCALAR_SOURCES, label)
                nudge = _unit_interval(tracer, "nudge", what, label)
                tracer_specs[pos] = {
                    "module": tracer_mod,
                    "code": code,
                    "keyword": keyword,
                    "values": values,
                    "nudge": nudge,
                }
        bc["tracers"] = tracer_specs
        return bc

    @staticmethod
    def _describe_assignment(bc):
        """Return a one-line summary of a resolved boundary for logging.

        Parameters
        ----------
        bc : dict
            Resolved assignment from :meth:`_resolve_open_boundary`.

        Returns
        -------
        str
            Summary listing each assigned variable as ``keyword(code)``.
        """
        parts = []
        for vname in ("elevation", "velocity", "temperature", "salinity"):
            var = bc[vname]
            if var["code"] == 0:
                continue
            text = f"{vname}={var['keyword']}({var['code']})"
            if "nudge" in var:
                text += f" nudge={var['nudge']}"
            if var["code"] in (-4, -5):
                text += (f" inflow_relax={var['spec']['inflow_relax']}"
                         f" outflow_relax={var['spec']['outflow_relax']}")
            parts.append(text)
        for tracer in bc["tracers"]:
            if tracer is not None:
                parts.append(
                    f"{tracer['module']}={tracer['keyword']}({tracer['code']}) "
                    f"nudge={tracer['nudge']}"
                )
        summary = ", ".join(parts) if parts else "no variables assigned"
        return f"Boundary '{bc['label']}' ({bc['flow_direction']}): {summary}"

    def _write_tidal_constituents(self, outf, spec, fields, num_nodes, what, label):
        """Write per-node tidal amplitude/phase blocks for one boundary.

        Parameters
        ----------
        outf : file-like
            Open ``bctides.in`` file.
        spec : dict
            Variable specification holding ``constituents``.
        fields : sequence of str
            Constituent fields written per node, e.g. ``("amplitude", "phase")``.
        num_nodes : int
            Number of nodes on the boundary.
        what : str
            Prefix used in error messages, e.g. ``"elev"``.
        label : str
            Boundary name used in error messages.
        """
        for tidal_constituent in spec["constituents"]:
            name = tidal_constituent.get("name")
            if name is None:
                if isinstance(tidal_constituent, dict) and len(tidal_constituent) > 0:
                    name = next(iter(tidal_constituent))
                else:
                    raise ValueError(
                        f" {what} tidal constituent name is empty for open boundary {label} \n"
                    )
            outf.write(name + "\n")
            entries = tidal_constituent[name]
            arrays = [entries[k][field] for k, field in enumerate(fields)]
            for kk in range(num_nodes):
                vals = [
                    self._handle_tidal_constituents_val(
                        arr, name, f"{what}_{field}", label, kk
                    )
                    for arr, field in zip(arrays, fields)
                ]
                sep = "   " if len(vals) == 2 else " "
                outf.write(sep.join(str(v) for v in vals) + "\n")

    def _write_open_boundary(self, outf, bc):
        """Write the type-code line and parameters for one open boundary.

        Parameters
        ----------
        outf : file-like
            Open ``bctides.in`` file.
        bc : dict
            Resolved assignment from :meth:`_resolve_open_boundary`.
        """
        label = bc["label"]
        num_nodes = bc["num_nodes"]
        elev = bc["elevation"]
        vel = bc["velocity"]
        codes = [num_nodes, elev["code"], vel["code"],
                 bc["temperature"]["code"], bc["salinity"]["code"]]
        codes += [t["code"] if t is not None else 0 for t in bc["tracers"]]
        outf.write(" ".join(str(c) for c in codes) + " !  " + label + "\n")

        if elev["code"] == CONSTANT_CODE:
            outf.write(str(elev["spec"]["source"]) + "\n")
        elif elev["code"] in (3, 5):
            self._write_tidal_constituents(
                outf, elev["spec"], ("amplitude", "phase"), num_nodes, "elev", label
            )

        vel_id = vel["code"]
        vspec = vel["spec"]
        if vel_id == CONSTANT_CODE:
            outf.write(str(vspec["source"]) + "\n")
        elif vel_id in (-4, -5):
            outf.write(f"{vspec['inflow_relax']} {vspec['outflow_relax']}\n")
        elif vel_id == -1:
            outf.write("eta_mean\n")
            for eta_m0 in vspec["eta_m0"]:
                outf.write(str(eta_m0) + "\n")
            outf.write("vn_mean\n")
            for vn_node in vspec["vnormal_nodes"]:
                vn_node = vn_node if isinstance(vn_node, list) else [vn_node]
                outf.write(" ".join(str(vn) for vn in vn_node) + "\n")
        elif vel_id == -2:
            outf.write(" ".join(str(c) for c in vspec["coefficients"]) + "\n")
        if vel_id in (3, 5, -5):
            self._write_tidal_constituents(
                outf, vspec, ("u_amplitude", "u_phase", "v_amplitude", "v_phase"),
                num_nodes, "vel", label,
            )

        for vname in ("temperature", "salinity"):
            var = bc[vname]
            if var["code"] == CONSTANT_CODE:
                outf.write(str(var["spec"]["source"]) + "\n")
            if var["code"] > 0:
                outf.write(str(var["nudge"]) + "\n")

        for tracer in bc["tracers"]:
            if tracer is None:
                continue
            if tracer["code"] == CONSTANT_CODE:
                outf.write(" ".join(str(v) for v in tracer["values"]) + "\n")
            outf.write(str(tracer["nudge"]) + "\n")

    def write_bctides(self, bctides_file):
        logger.info("Writing bctides file %s", bctides_file)
 
        with open(bctides_file, "w") as outf:
            tt = self.date.strftime("%m/%d/%Y %H:%M")
            outf.write(tt)
            outf.write("\n")
            # normalize tidal-constituent input formats so downstream code can
            # normalize tidal-constituent input formats so downstream code can
            # expect a list of dicts with explicit keys.
 
            if self.earth_tides and "constituents" in self.earth_tides:
                self.earth_tides["constituents"] = self._norm_earth(
                    self.earth_tides["constituents"]
                )
            if self.boundary_tides and "constituents" in self.boundary_tides:
                self.boundary_tides["constituents"] = self._norm_boundary(
                    self.boundary_tides["constituents"]
                )
            ##out earth tidal potential if any
            if self.earth_tides:
                num_tidal = len(self.earth_tides["constituents"])
                cutoff_depth = self.earth_tides["cutoff_depth"]
                outf.write(str(num_tidal) + " " + str(cutoff_depth) + 
                           "! # of constituents used in earth tidal potential, \
                           cut-off depth for applying tidal potential\n")
                for i in range(num_tidal):
                    outf.write(
                        self.earth_tides["constituents"][i]["name"]
                    )
                    outf.write("\n")
                    amp = self.earth_tides["constituents"][i]["amplitude"]
                    node_factor = self.earth_tides["constituents"][i][
                        "node_factor"
                    ]
                    freq = self.earth_tides["constituents"][i][
                        "angular_frequency"
                    ]
                    eqa = self.earth_tides["constituents"][i][
                        "equilibrium_argument"
                    ]
                    outf.write(
                        str(amp)
                        + " "
                        + str(node_factor)
                        + " "
                        + str(freq)
                        + " "
                        + str(eqa)
                    )
                    outf.write("\n")
            else:
                num_tidal = 0
                cutoff_depth = 40
                outf.write(str(num_tidal) + " " + str(cutoff_depth) + 
                           " ! # of constituents used in earth tidal potential," +
                           "cut-off depth for applying tidal potential\n")

            ## out boundary forcing tidal
            if self.boundary_tides:
                num_tidal = len(self.boundary_tides["constituents"])
                outf.write(str(num_tidal) + " ! total # of tidal boundary forcing frequencies  \n")
                for i in range(num_tidal):
                    outf.write(
                        self.boundary_tides["constituents"][i]["name"]
                    )
                    outf.write("\n")
                    freq = self.boundary_tides["constituents"][i][
                        "angular_frequency"
                    ]
                    node_factor = self.boundary_tides["constituents"][i][
                        "node_factor"
                    ]
                    eqa = self.boundary_tides["constituents"][i][
                        "equilibrium_argument"
                    ]
                    outf.write(str(freq) + " " + str(node_factor) + " " + str(eqa))
                    outf.write("\n")
            else:
                num_tidal = 0
                outf.write(str(num_tidal) + " ! total # of tidal boundary forcing frequencies  \n")

            
            hgrid = self.hgrid
        
            
            hgrd_open_boundaries = list(
            boundary
            for boundary in hgrid.boundaries
            if boundary.btype == BoundaryType.OPEN
            )

            if len(hgrd_open_boundaries) != len(self.open_boundaries):
                raise ValueError(
                    "boundary YAML has different number of openboundary than model hgrid, YAML:%s hgrid:%s"
                    % (len(self.open_boundaries), len(hgrd_open_boundaries))
                )
            else:
                num_open_boundaries = len(self.open_boundaries)
                ## check to make sure all open_boundaries name match
                ## with hgrid boundary names
                ## summary all the bounary names from self.open_boundaries and
                ## hgrid_open_boundaries, then compare them
                bname_yaml_lst = []
                bname_hgrid_lst = []
                for i in range(num_open_boundaries):
                    bname_yaml_lst.append(self.open_boundaries[i].get("name",""))
                    if  hgrd_open_boundaries[i].comment is None:
                        bname_hgrid_lst.append("")  
                    else:
                        bname_hgrid =  hgrd_open_boundaries[i].comment.strip().split('"')[1]
                        bname_hgrid_lst.append(bname_hgrid)
                ## if both boundary names are same, skip the check
                if bname_yaml_lst != bname_hgrid_lst:
                    ## if not same, check each boundary name defined in ymal can
                    ## be found in hgrid boundary names, and reorder self.open_boundaries
                    ## according to hgrid boundary order
                    reordered_open_boundaries_tmp = [None]*num_open_boundaries
                    order_indices = []
                    for i in range(num_open_boundaries):
                        bname_yaml = self.open_boundaries[i].get("name","")
                        if bname_yaml not in bname_hgrid_lst:
                            raise ValueError(
                                "boundary name %s from YAML not found in hgrid boundary names %s"
                                % (bname_yaml,bname_hgrid_lst)
                            )
                        else:
                            index_in_hgrid = bname_hgrid_lst.index(bname_yaml)   
                            reordered_open_boundaries_tmp[index_in_hgrid] = self.open_boundaries[i] 
                            order_indices.append(index_in_hgrid)
                    ## reorered self.open_boundaries should be same as hgrid boundary names now
                    ## if not, raise error
                    for i in range(num_open_boundaries):
                        bname_yaml = reordered_open_boundaries_tmp[i].get("name","")
                        bname_hgrid = bname_hgrid_lst[i]
                        if bname_yaml != bname_hgrid:
                            raise ValueError(
                                "after reordering, boundary name mismatch between YAML and hgrid: %s vs %s"
                                % (bname_yaml, bname_hgrid)
                            )
                    self.open_boundaries = reordered_open_boundaries_tmp

                    
            
                # for i in range(num_open_boundaries):
                #     bname_yaml = self.open_boundaries[i]["name"]
                #     if  hgrd_open_boundaries[i].comment is None:
                #         print("Warning: hgrid boundary %s has no name, skipping boundary name check" % i)
                #         continue ## skip if hgrid boundary name is empty
                #     bname_hgrid =  hgrd_open_boundaries[i].comment.strip().split('"')[1]
                #     ## if bname in hgrid is empty, skip the check   
                #     if bname_hgrid == "":
                #         continue
                #     if bname_yaml != bname_hgrid:
                #         raise ValueError(
                #             "boundary name mismatch between YAML and hgrid: %s vs %s"
                #             % (bname_yaml, bname_hgrid)
                #         )

                # num_tracer_mod = 0
                # tracer_mods = []
                # tracer_mod_pos = {}
                # for i in range(num_open_boundaries):
                #     if "tracers" in self.open_boundaries[i].keys():
                #         boundary_tracer_mod_num = len(
                #             self.open_boundaries[i]["tracers"]
                #         )
                #         if boundary_tracer_mod_num > num_tracer_mod:
                #             num_tracer_mod = boundary_tracer_mod_num
                #         for j in range(boundary_tracer_mod_num):
                #             tracer_mod = self.open_boundaries[i]["tracers"][j]["tracer"]
                #             if not (tracer_mod in self.tracer_lst.keys()):
                #                 raise ValueError(
                #                     tracer_mod + " is not a supported tracer module\n"
                #                 )
                #             else:
                #                 if not (tracer_mod in tracer_mods):
                #                     tracer_mods.append(tracer_mod)
                # tracer_mods.sort(key=lambda d: self.tracer_lst[d])
                # for kk in range(len(tracer_mods)):
                #     tracer_mod_pos[tracer_mods[kk]] = kk

                ## write out number of boundary here
                outf.write(str(num_open_boundaries) + " ! # of open boundary \n")

                for i in range(num_open_boundaries):
                    bc = self._resolve_open_boundary(
                        i, self.open_boundaries[i], len(hgrid.boundaries[i].nodes)
                    )
                    logger.info(self._describe_assignment(bc))
                    self._write_open_boundary(outf, bc)


if __name__ == "__main__":
    by = bctide_cli()



