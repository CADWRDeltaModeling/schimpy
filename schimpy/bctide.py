# -*- coding: utf-8 -*-
"""

This class will read in a yaml file which defines SCHISM boundary and generate bctides.in.
It supports  all type of elevation, velocity, salinity, temperature and tracer boundaries.
At the end of this script, there is a synthetic example which demonstrates the format of bctide YAML file.

"""
from schimpy.schism_yaml import load
from schimpy.schism_mesh import SchismMeshGr3Reader, BoundaryType
import datetime
from zoneinfo import ZoneInfo
import numbers
import numpy as np
import click


__all__ = ["bctide_cli"]


@click.command(
    help=(
        "The tool create SCHISM bctide inputs by reading a number of open boundary condition yaml and mesh \
         from a  main input yaml \n\n"
        "Example:\n"
        "  bctide bctide_main.yaml\n "
    )
)
@click.argument('main_yaml', type=click.Path(exists=True))

def bctide_cli(main_yaml):
     """ read in main yaml file which contians the grid topology and
         boundary condtion yaml files, this input file looks like below
         mesh:
         mesh_inputfile: ./hgrid.gr3  
         open_boundaries: !include open_boundary.yaml
  
         bctides:
         bctides.in.2d: !include bctides.in.2d.yaml
         bctides.in.3d: !include bctides.in.3d.yaml
         bctides.in.sed.3d: !include bctides.in.sed.3d.yaml

         for each boundary condition yaml file, it will generate the corresponding boundary
         class and write out the SCHISM boundary condition btdides input file.
         """
     with open(main_yaml, "r") as fn:
         main_yaml_dict = load(fn)

     mesh_in = main_yaml_dict["mesh"]["mesh_inputfile"]
     sr = SchismMeshGr3Reader()

     hgrid = sr.read(mesh_in)
     open_boundary_segments = None
     if "open_boundaries" in main_yaml_dict["mesh"].keys():
         open_boundary_segments = main_yaml_dict["mesh"]["open_boundaries"]

     bctides_dic = main_yaml_dict.get("bctides")
#     for bctides_key in bctides_dic.keys():
     for out in bctides_dic:
         print("---------processing bctides file %s----------" % out)
         bctides_yaml = bctides_dic[out]
         by = load_boundary(hgrid,bctides_yaml,open_boundary_segments)
  #       by.write_bctides(bctides_key)
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

        self.elev_source = {
            "elev history": 1,
            "constant": 2,
            "tidal": 3,
            "2D elev history": 4,
            "tidal and 2D elev history": 5,
            "unspecified": 0,
        }
        self.vel_source = {
            "flux history": 1,
            "constant": 2,
            "tidal": 3,
            "3D uv history": 4,
            "tidal and 3D uv history": 5,
            "flather 1": -1,
            "unspecified": 0,
        }

        self.temp_source = {
            "temperature history": 1,
            "constant": 2,
            "initial profile": 3,
            "3D temperature history": 4,
            "unspecified": 0,
        }

        self.salt_source = {
            "salinity history": 1,
            "constant": 2,
            "initial profile": 3,
            "3D salinity history": 4,
            "unspecified": 0,
        }

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
        self.tracer_source = {
            "time history": 1,
            "constant": 2,
            "initial profile": 3,
            "unspecified": 0,
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

    def write_bctides(self, bctides_file):
 
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
                outf.write(str(num_tidal) + " " + str(cutoff_depth) + "\n")
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
                outf.write(str(num_tidal) + " " + str(cutoff_depth) + "\n")

            ## out boundary forcing tidal
            if self.boundary_tides:
                num_tidal = len(self.boundary_tides["constituents"])
                outf.write(str(num_tidal) + "\n")
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
                outf.write(str(num_tidal) + "\n")

            
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
                                "boundary name %s from YAML not found in hgrid boundary names"
                                % bname_yaml
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
                outf.write(str(num_open_boundaries) + "\n")

                for i in range(num_open_boundaries):
                    bound_name = self.open_boundaries[i]["name"]
                    print("processing open boundary %s " % bound_name)
                    num_nodes = len(hgrid.boundaries[i].nodes)
                    node_id_lst = hgrid.boundaries[i].nodes
                    elev_id = 0
                    elev_boundary = None
                    ## if flow_direction is given and is outflow, then
                    ## this boundary shouldn't have any tracer boundary, such
                    ## as temperature, salinity or other scalars
                    

                    ## flow_direction is must for open boundary, if not given, raise error
                    if not("flow_direction" in self.open_boundaries[i].keys()):
                        raise ValueError(
                            f"Boundary '{bound_name}' is missing required 'flow_direction' specification. Please specify 'flow_direction' as either 'inflow' or 'outflow'."
                        )

                    # variables in YAML are provided as a list of dicts; helper to find a variable entry
                    vars_list = self.open_boundaries[i].get("variables", [])
                    def _get_var(name):
                        for vv in vars_list:
                            if isinstance(vv, dict) and name in vv:
                                return vv[name]
                        return None
                    
                    flow_direction = self.open_boundaries[i]["flow_direction"]
                    if flow_direction == "outflow":
                        if "temperature" in self.open_boundaries[i].keys():
                            raise ValueError(
                                f"Boundary '{bound_name}' is configured as an outflow; a temperature boundary "
                                "specification is not allowed."
                            )
                        if "salinity" in self.open_boundaries[i].keys():
                            raise ValueError(
                                f"Boundary '{bound_name}' is configured as an outflow; a salinity boundary "
                                "specification is not allowed."
                            )
                        if "tracers" in self.open_boundaries[i].keys():
                            raise ValueError(
                                f"Boundary '{bound_name}' is configured as an outflow; tracer "
                                "specification is not allowed."
                            )
                    elif flow_direction == "inflow" or flow_direction == "bidirection":

                        check_list = ["temperature", "salinity","tracers"]
                        if len(self.tracer_modules) == 0: ## if no tracer module is defined, then tracers should not be specified in the boundary
                            check_list = ["temperature", "salinity"]
                            ## if tracers is specified in the boundary but no tracer module is defined, raise error
                            if "tracers" in self.open_boundaries[i].keys():
                                raise ValueError(
                                    f"Boundary '{bound_name}' has specified 'tracers' boundary conditions, but no tracer modules are defined in the bctides YAML. \
                                        Please add tracer modules to the top-level 'modules' section or remove 'tracers' specification from the boundary."
                                )
                        for var in check_list:
                            self.check_inflow_boundary_spec(_get_var(var), bound_name,var)
                    else:
                        raise ValueError(
                            f"Boundary '{bound_name}' has invalid 'flow_direction' value: {flow_direction}. "
                            "Expected 'inflow', 'outflow', or 'bidirection'."
                        )   




                    if _get_var("elevation") is not None:
                        elev_boundary = _get_var("elevation")
                        elev_source = elev_boundary["source"]
                        elev_key = elev_source
                        if isinstance(elev_source, numbers.Number):
                            elev_key = "constant"
                        if ("tidal" in elev_key) and ("elev2D.th.nc" in elev_key):
                            elev_key = "tidal elev2D.th.nc"

                        try:
                            elev_id = self.elev_source[elev_key]
                        except:
                            raise ValueError(
                                elev_key + " elevation boundary is not supported"
                            )

                    vel_id = 0
                    vel_boundary = None
                    if _get_var("velocity") is not None:
                        vel_boundary = _get_var("velocity")
                        vel_source = vel_boundary["source"]
                        vel_key = vel_source

                        if isinstance(vel_source, numbers.Number):
                            vel_key = "constant"
                        if ("tidal" in vel_key) and ("uv3D.th.nc" in vel_key):
                            vel_key = "tidal uv3D.th.nc"

                        if ("flather " in vel_key) and ("1 " in vel_key):
                            vel_key = "flather 1"

                        try:
                            vel_id = self.vel_source[vel_key]
                        except:
                            raise ValueError(
                                vel_key + " velocity boundary is not supported"
                            )

                    ## output temperature and salinity boundary
                    temp_id = 0
                    temp_boundary = None
                    if _get_var("temperature") is not None:
                        temp_boundary = _get_var("temperature")
                        temp_source = temp_boundary["source"]
                        temp_key = temp_source
                        if isinstance(temp_source, numbers.Number):
                            temp_key = "constant"

                        try:
                            temp_id = self.temp_source[temp_key]
                        except:
                            raise ValueError(
                                temp_key + " temperature boundary is not supported"
                            )
                    salt_id = 0
                    salt_boundary = None
                    if _get_var("salinity") is not None:
                        salt_boundary = _get_var("salinity")
                        salt_source = salt_boundary["source"]
                        salt_key = salt_source
                        if isinstance(salt_source, numbers.Number):
                            salt_key = "constant"

                        try:
                            salt_id = self.salt_source[salt_key]
                        except:
                            raise ValueError(
                                salt_key + " temperature boundary is not supported"
                            )

                    ## output tracer boundary
                    tracer_boundary_sources = [0] * len(self.tracer_mod_pos)
                    ## this list save sorted tracer boundary index according to SCHISM code order
                    tracer_boundary_lst_sorted = []
                    if _get_var("tracers") is not None:
                        tracers_var = _get_var("tracers")
                        boundary_tracer_mod_num = len(tracers_var)
                        for j in range(boundary_tracer_mod_num):
                            tracer_boundary = tracers_var[j]["source"]
                            tracer_mod = tracers_var[j]["module"]
                            ## if tracer mod not in self.tracer_mod_pos, raise error
                            if not (tracer_mod in self.tracer_mod_pos.keys()):  
                                raise ValueError(
                                    f"Tracer module '{tracer_mod}' is not defined in the bctides 'modules' section. "
                                    f"Please add it to the top-level modules list. "
                                    f"Valid modules: {', '.join(sorted(self.tracer_mod_pos.keys()))}"
                                )

                            tracer_boundary_key = tracer_boundary
                            ## set boundary key to "constant" if source is a number or list of numbers
                            if isinstance(tracer_boundary, numbers.Number):
                                tracer_boundary_key = "constant"
                            elif isinstance(tracer_boundary, list):
                                boundary_tracer_mod_num = len(tracer_boundary)
                                ## if curent tracer mod has multiple tracers, check if the
                                ## length of tracer_boundary match the number of tracers defined
                                ## in self.tracer_mod_num
                                if boundary_tracer_mod_num != self.tracer_mod_num[tracer_mod]:
                                    raise ValueError( self.open_boundaries[i]["name"] + ": "+
                                        tracer_mod + " boundary source length "
                                        + str(boundary_tracer_mod_num)
                                        + " does not match number of tracers defined in module section: "
                                        + str(self.tracer_mod_num[tracer_mod])
                                    )

                                if all(
                                    isinstance(x, numbers.Number)
                                    for x in tracer_boundary
                                ):
                                    tracer_boundary_key = "constant"
                                else:
                                    raise ValueError(
                                        tracer_mod
                                        + " "
                                        + str(tracer_boundary)
                                        + " boundary is not supported"
                                    )
                            tracer_boundary_source = self.tracer_source[tracer_boundary_key]
                            pos = self.tracer_mod_pos[tracer_mod]
                            tracer_boundary_sources[pos] = tracer_boundary_source
                            tracer_boundary_lst_sorted.append(j)

                        tracer_boundary_lst_sorted.sort(
                            key=lambda jj: self.tracer_mod_pos[
                                tracers_var[jj]["module"]
                            ]
                        )
                    outf.write(str(num_nodes) + " ")
                    outf.write(str(elev_id) + " ")
                    outf.write(str(vel_id) + " ")
                    outf.write(str(temp_id) + " ")
                    outf.write(str(salt_id) + " ")
                    for ii in range(len(self.tracer_mod_pos)):
                        outf.write(str(tracer_boundary_sources[ii]) + " ")
                    outf.write("\n")

                    ## output all the boundary parameters
                    ## elevation boundary parameters
                    if elev_id == 2:
                        ## const
                        outf.write(str(elev_source))
                    elif (elev_id == 3) or (elev_id == 5):
                        ## tidal forcing
                        # write each tidal constituent: name line then amplitude/phase per node
                        for tidal_constituent in elev_boundary["constituents"]:
                            # name must exist (YAML examples use 'name' after normalization)
                            name = tidal_constituent.get("name")
                            if name is None:
                                # fallback to first key if mapping style was used
                                if isinstance(tidal_constituent, dict) and len(tidal_constituent) > 0:
                                    name = next(iter(tidal_constituent))
                                else:
                                    ## raise error  
                                    raise ValueError(f" elev tidal constituent name is empty for open boundary {i} \n")
                            outf.write(name + "\n")
                            amp = tidal_constituent[name][0]["amplitude"]
                            phase = tidal_constituent[name][1]["phase"]

                            for kk in range(num_nodes):
                                node_id = node_id_lst[kk]
                                #x = hgrid.nodes[node_id, 0]
                                #y = hgrid.nodes[node_id, 1]
                                amp_val = self._handle_tidal_constituents_val(amp,name,"elev_amplitude",i,kk)
                                phase_val = self._handle_tidal_constituents_val(phase,name,"elev_phase",i,kk)

                                outf.write(str(amp_val) + "   " + str(phase_val) + "\n")

                    ## velocity boundary parameters
                    if vel_id == 2:
                        ## constp 
                        outf.write(str(vel_source))
                    elif (vel_id == 3) or (vel_id == 5):
                        ## tidal forcing
                        for tidal_constituent in vel_boundary["constituents"]:
                            # name must exist (YAML examples use 'name' after normalization)
                            name = tidal_constituent.get("name")
                            if name is None:
                                # fallback to first key if mapping style was used
                                if isinstance(tidal_constituent, dict) and len(tidal_constituent) > 0:
                                    name = next(iter(tidal_constituent))
                                else:
                                    ## raise error  
                                    raise ValueError(f" vel tidal constituent name is empty for open boundary {i} \n")
                            outf.write(name + "\n")
                            u_amp = tidal_constituent[name][0]["u_amplitude"]
                            u_phase = tidal_constituent[name][1]["u_phase"]
                            v_amp = tidal_constituent[name][2]["v_amplitude"]
                            v_phase = tidal_constituent[name][3]["v_phase"]

                            for kk in range(num_nodes):

                                node_id = node_id_lst[kk]


                                u_amp_val = self._handle_tidal_constituents_val(u_amp,name,"u_amplitude",i,kk)
                                u_phase_val = self._handle_tidal_constituents_val(u_phase,name,"u_phase",i,kk)
                                v_amp_val = self._handle_tidal_constituents_val(v_amp,name,"v_amplitude",i,kk)
                                v_phase_val = self._handle_tidal_constituents_val(v_phase,name,"v_phase",i,kk)
                                
                                
                                outf.write(
                                    str(u_amp_val)
                                    + " "
                                    + str(u_phase_val)
                                    + " "
                                    + str(v_amp_val)
                                    + " "
                                    + str(v_phase_val)
                                    + "\n"
                                )
                    elif vel_id == -4:
                        relax_inflow = vel_boundary["inflow_relax"]
                        relax_outflow = vel_boundary["outflow_relax"]
                        outf.write(str(relax_inflow) + " " + str(relax_outflow) + "\n")
                    elif vel_id == -1:
                        eta_mean = vel_boundary["eta_mean"]
                        outf.write(str(eta_mean) + "\n")
                        for kk in range(num_nodes):
                            eta_m0 = vel_boundary["eta_m0"][kk]
                            outf.write(str(eta_m0) + "\n")
                        vn_mean = vel_boundary["vnormal_mean"]
                        outf.write(str(vn_mean) + "\n")
                        for kk in range(num_nodes):
                            vn_node = vel_boundary["vnormal_nodes"][kk]
                            for vn in vn_node:
                                outf.write(str(vn) + " ")
                            outf.write("\n")

                    ## temperature bc parameters
                    if temp_id == 2:
                        temp_bc = temp_boundary["source"]
                        outf.write(str(temp_bc) + "\n")
                    if temp_id > 0:
                        nudge = temp_boundary["nudge"]
                        outf.write(str(nudge) + "\n")

                    ## salt bc parameters
                    if salt_id == 2:
                        salt_bc = salt_boundary["source"]
                        outf.write(str(salt_bc) + "\n")
                    if salt_id > 0:
                        nudge = salt_boundary["nudge"]
                        outf.write(str(nudge) + "\n")

                    ## tracer bc paraemters
                    for ii in range(len(tracer_boundary_lst_sorted)):
                        tracer_index = tracer_boundary_lst_sorted[ii]
                        tracer = _get_var("tracers")[tracer_index]
                        tracer_bc_source = tracer_boundary_sources[ii]
                        if tracer_bc_source == 2:
                            tracer_bc_const = tracer["source"]
                            if isinstance(tracer_bc_const, list):
                                for val in tracer_bc_const:
                                    outf.write(str(val) + " ")
                            else:
                                outf.write(str(tracer_bc_const))
                            outf.write("\n")
                        if tracer_bc_source > 0:
                            relax = tracer["relax"]
                            outf.write(str(relax) + "\n ")


if __name__ == "__main__":
    by = bctide_cli()



