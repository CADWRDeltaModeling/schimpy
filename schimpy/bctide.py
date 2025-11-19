# -*- coding: utf-8 -*-
"""

This class will read in a yaml file which defines SCHISM boundary and generate bctides.in.
It supports  all type of elevation, velocity, salinity, temperature and tracer boundaries.
At the end of this script, there is a synthetic example which demonstrates the format of bctide YAML file.

"""
from schimpy.schism_yaml import load
from schimpy.schism_mesh import SchismMeshGr3Reader
import numbers


__all__ = ["load_boundary"]


def load_boundary(fn):
    return boundary(fn)


class boundary(object):
    """
    A class to generate boundary condition input file for SCHISM
    """

    def __init__(self, yaml_fn=None):

        with open(yaml_fn, "r") as fn:
            bc_yaml = load(fn)

            main_id = "bctides"
            self.date = bc_yaml[main_id]["date"]
            if "earth_tidals" in bc_yaml[main_id].keys():
                self.earth_tidals = bc_yaml[main_id]["earth_tidals"]
            else:
                self.earth_tidals = None
            if "bounary_forcing_tidals" in bc_yaml[main_id].keys():
                self.boundary_tidals = bc_yaml[main_id]["bounary_forcing_tidals"]
            else:
                self.boundary_tidals = None
            
            self.open_boundaries = bc_yaml[main_id]["open_boundaries"]
            self.hgrid = bc_yaml[main_id]["hgrid_input_file"]

            self.elev_type = {
                "elev.th": 1,
                "constant": 2,
                "tidal": 3,
                "elev2D.th.nc": 4,
                "tidal elev2D.th.nc": 5,
                "none": 0,
            }
            self.vel_type = {
                "flux.th": 1,
                "constant": 2,
                "tidal": 3,
                "uv3D.th.nc": 4,
                "tidal uv3D.th.nc": 5,
                "flather 1": -1,
                "none": 0,
            }

            self.temp_type = {
                "TEM_1.th": 1,
                "constant": 2,
                "initial profile": 3,
                "TEM_3D.th.nc": 4,
                "none": 0,
            }

            self.salt_type = {
                "SAL_1.th": 1,
                "constant": 2,
                "initial profile": 3,
                "SAL_3D.th.nc": 4,
                "none": 0,
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
            self.tracer_type = {
                "time history": 1,
                "constant": 2,
                "initial profile": 3,
                "none": 0,
            }

    def write_bctides(self, bctides_file):

        with open(bctides_file, "w") as outf:
            tt = self.date.strftime("%Y-%m-%d %H:%M")
            outf.write(tt)
            outf.write("\n")

            ##out earth tidal potential if any
            if self.earth_tidals:
                num_tidal = len(self.earth_tidals["tidal_constituents"])
                cutoff_depth = self.earth_tidals["tidal_cutoff_depth"]
                outf.write(str(num_tidal) + " " + str(cutoff_depth) + "\n")
                for i in range(num_tidal):
                    outf.write(
                        self.earth_tidals["tidal_constituents"][i]["name"]
                    )
                    outf.write("\n")
                    amp = self.earth_tidals["tidal_constituents"][i]["amplitude"]
                    node_factor = self.earth_tidals["tidal_constituents"][i][
                        "node_factor"
                    ]
                    freq = self.earth_tidals["tidal_constituents"][i][
                        "angular_frequency"
                    ]
                    eqa = self.earth_tidals["tidal_constituents"][i][
                        "earth_equilibrium_argument"
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
            if self.boundary_tidals:
                num_tidal = len(self.boundary_tidals["tidal_constituents"])
                outf.write(str(num_tidal) + "\n")
                for i in range(num_tidal):
                    outf.write(
                        self.boundary_tidals["tidal_constituents"][i]["name"]
                    )
                    outf.write("\n")
                    freq = self.boundary_tidals["tidal_constituents"][i][
                        "angular_frequency"
                    ]
                    node_factor = self.boundary_tidals["tidal_constituents"][i][
                        "node_factor"
                    ]
                    eqa = self.boundary_tidals["tidal_constituents"][i][
                        "earth_equilibrium_argument"
                    ]
                    outf.write(str(freq) + " " + str(node_factor) + " " + str(eqa))
                    outf.write("\n")
            else:
                num_tidal = 0
                outf.write(str(num_tidal) + "\n")

            ## output open boundaries
            ## read in mesh grid file
            sr = SchismMeshGr3Reader()
            hgrid = sr.read(self.hgrid)

            if len(hgrid.boundaries) < len(self.open_boundaries):
                raise ValueError(
                    "boundary YAML has more number of openbounary than model hgrid, YAML:%s hgrid:%s"
                    % (len(self.open_boundaries), len(hgrid.boundaries))
                )
            else:
                num_open_boundaries = len(self.open_boundaries)
                num_tracer_mod = 0
                tracer_mods = []
                tracer_mod_pos = {}
                for i in range(num_open_boundaries):
                    if "tracers" in self.open_boundaries[i].keys():
                        boundary_tracer_mod_num = len(
                            self.open_boundaries[i]["tracers"]
                        )
                        if boundary_tracer_mod_num > num_tracer_mod:
                            num_tracer_mod = boundary_tracer_mod_num
                        for j in range(boundary_tracer_mod_num):
                            tracer_mod = self.open_boundaries[i]["tracers"][j]["tracer"]
                            if not (tracer_mod in self.tracer_lst.keys()):
                                raise ValueError(
                                    tracer_mod + " is not a supported tracer module\n"
                                )
                            else:
                                if not (tracer_mod in tracer_mods):
                                    tracer_mods.append(tracer_mod)
                tracer_mods.sort(key=lambda d: self.tracer_lst[d])
                for kk in range(len(tracer_mods)):
                    tracer_mod_pos[tracer_mods[kk]] = kk

                for i in range(num_open_boundaries):
                    num_nodes = len(hgrid.boundaries[i].nodes)
                    node_id_lst = hgrid.boundaries[i].nodes
                    elev_id = 0
                    elev_boundary = None
                    if "elevation_boundary" in self.open_boundaries[i].keys():
                        elev_boundary = self.open_boundaries[i]["elevation_boundary"]
                        elev_source = elev_boundary["type"]
                        elev_key = elev_source
                        if isinstance(elev_source, numbers.Number):
                            elev_key = "constant"
                        if ("tidal" in elev_key) and ("elev2D.th.nc" in elev_key):
                            elev_key = "tidal elev2D.th.nc"

                        try:
                            elev_id = self.elev_type[elev_key]
                        except:
                            raise ValueError(
                                elev_key + " elevation boundary is not supported"
                            )

                    vel_id = 0
                    vel_boundary = None
                    if "velocity_boundary" in self.open_boundaries[i].keys():
                        vel_boundary = self.open_boundaries[i]["velocity_boundary"]
                        vel_source = vel_boundary["type"]
                        vel_key = vel_source

                        if isinstance(vel_source, numbers.Number):
                            vel_key = "constant"
                        if ("tidal" in vel_key) and ("uv3D.th.nc" in vel_key):
                            vel_key = "tidal uv3D.th.nc"

                        if ("flather " in vel_key) and ("1 " in vel_key):
                            vel_key = "flather 1"

                        try:
                            vel_id = self.vel_type[vel_key]
                        except:
                            raise ValueError(
                                vel_key + " velocity boundary is not supported"
                            )

                    ## output temperature and salinity boundary
                    temp_id = 0
                    temp_boundary = None
                    if "temperature_boundary" in self.open_boundaries[i].keys():
                        temp_boundary = self.open_boundaries[i]["temperature_boundary"]
                        temp_source = temp_boundary["type"]
                        temp_key = temp_source
                        if isinstance(temp_source, numbers.Number):
                            temp_key = "constant"

                        try:
                            temp_id = self.vel_type[vel_key]
                        except:
                            raise ValueError(
                                temp_key + " temperature boundary is not supported"
                            )
                    salt_id = 0
                    salt_boundary = None
                    if "salinity_boundary" in self.open_boundaries[i].keys():
                        salt_boundary = self.open_boundaries[i]["salinity_boundary"]
                        salt_source = salt_boundary["type"]
                        salt_key = salt_source
                        if isinstance(salt_source, numbers.Number):
                            salt_key = "constant"

                        try:
                            salt_id = self.salt_type[salt_key]
                        except:
                            raise ValueError(
                                salt_key + " temperature boundary is not supported"
                            )

                    ## output tracer boundary
                    tracer_boundary_types = [0] * num_tracer_mod
                    ## this list save sorted tracer boundary index according to SCHISM code order
                    tracer_boundary_lst_sorted = []
                    if "tracers" in self.open_boundaries[i].keys():
                        boundary_tracer_mod_num = len(
                            self.open_boundaries[i]["tracers"]
                        )
                        for j in range(boundary_tracer_mod_num):
                            tracer_boundary = self.open_boundaries[i]["tracers"][j]["type"]
                            tracer_mod = self.open_boundaries[i]["tracers"][j]["tracer"]
                            tracer_boundary_key = tracer_boundary
                            if isinstance(tracer_boundary, numbers.Number):
                                tracer_boundary_key = "constant"
                            elif isinstance(tracer_boundary, list):
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
                            tracer_boundary_type = self.tracer_type[tracer_boundary_key]
                            pos = tracer_mod_pos[tracer_mod]
                            tracer_boundary_types[pos] = tracer_boundary_type
                            tracer_boundary_lst_sorted.append(j)

                        tracer_boundary_lst_sorted.sort(
                            key=lambda jj: tracer_mod_pos[
                                self.open_boundaries[i]["tracers"][jj]["tracer"]
                            ]
                        )
                    outf.write(str(num_nodes) + " ")
                    outf.write(str(elev_id) + " ")
                    outf.write(str(vel_id) + " ")
                    outf.write(str(temp_id) + " ")
                    outf.write(str(salt_id) + " ")
                    for ii in range(num_tracer_mod):
                        outf.write(str(tracer_boundary_types[ii]) + " ")
                    outf.write("\n")

                    ## output all the boundary parameters
                    ## elevation boundary parameters
                    if elev_id == 2:
                        ## const
                        outf.write(str(elev_source))
                    elif (elev_id == 3) or (elev_id == 5):
                        ## tidal forcing
                        num_tidal_constituents = elev_boundary["tidal_constituents"]
                        for tidal_constituent in elev_boundary["tidal_constituents"]:
                            outf.write(tidal_constituent["name"] + "\n")
                            for kk in range(num_nodes):
                                amp = tidal_constituent["amplitude"]
                                phase = tidal_constituent["phase"]
                                node_id = node_id_lst[kk]
                                x = hgrid.nodes[node_id, 0]
                                y = hgrid.nodes[node_id, 1]
                                if isinstance(amp, numbers.Number):
                                    amp_val = amp
                                elif isinstance(amp, str):  # try to evalate expression
                                    try:
                                        amp_val = eval(amp)
                                    except:
                                        raise ValueError(
                                            amp
                                            + " is not a supported tidal amp expression\n"
                                        )
                                if isinstance(phase, numbers.Number):
                                    phase_val = phase
                                elif isinstance(
                                    phase, str
                                ):  # try to evalate expression
                                    try:
                                        phase_val = eval(phase)
                                    except:
                                        raise ValueError(
                                            phase
                                            + " is not a supported tidal phase expression\n"
                                        )
                                outf.write(str(amp_val) + "   " + str(phase_val) + "\n")

                    ## velocity boundary parameters
                    if vel_id == 2:
                        ## const
                        outf.write(str(vel_source))
                    elif (vel_id == 3) or (vel_id == 5):
                        ## tidal forcing
                        for tidal_constituent in vel_boundary["tidal_constituents"]:
                            outf.write(tidal_constituent["name"] + "\n")
                            for kk in range(num_nodes):
                                u_amp = tidal_constituent["u_amplitude"]
                                u_phase = tidal_constituent["u_phase"]
                                v_amp = tidal_constituent["v_amplitude"]
                                v_phase = tidal_constituent["v_phase"]
                                node_id = node_id_lst[kk]
                                x = hgrid.nodes[node_id, 0]
                                y = hgrid.nodes[node_id, 1]
                                
                                if isinstance(u_amp, numbers.Number):
                                    u_amp_val = u_amp
                                elif isinstance(
                                    u_amp, str
                                ):  # try to evalate expression
                                    try:
                                        u_amp_val = eval(u_amp)
                                    except:
                                        raise ValueError(
                                            u_amp
                                            + " is not a supported tidal amp expression\n"
                                        )
                                if isinstance(u_phase, numbers.Number):
                                    u_phase_val = u_phase
                                elif isinstance(
                                    u_phase, str
                                ):  # try to evalate expression
                                    try:
                                        u_phase_val = eval(u_phase)
                                    except:
                                        raise ValueError(
                                            u_phase
                                            + " is not a supported tidal phase expression\n"
                                        )

                                if isinstance(v_amp, numbers.Number):
                                    v_amp_val = v_amp
                                elif isinstance(
                                    v_amp, str
                                ):  # try to evalate expression
                                    try:
                                        v_amp_val = eval(v_amp)
                                    except:
                                        raise ValueError(
                                            u_amp
                                            + " is not a supported tidal amp expression\n"
                                        )
                                if isinstance(v_phase, numbers.Number):
                                    v_phase_val = v_phase
                                elif isinstance(
                                    v_phase, str
                                ):  # try to evalate expression
                                    try:
                                        v_phase_val = eval(v_phase)
                                    except:
                                        raise ValueError(
                                            v_phase
                                            + " is not a supported tidal phase expression\n"
                                        )

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
                        temp_bc = temp_boundary["type"]
                        outf.write(str(temp_bc) + "\n")
                    if temp_id > 0:
                        nudge = temp_boundary["nudge"]
                        outf.write(str(nudge) + "\n")

                    ## salt bc parameters
                    if salt_id == 2:
                        salt_bc = salt_boundary["type"]
                        outf.write(str(salt_bc) + "\n")
                    if salt_id > 0:
                        nudge = salt_boundary["nudge"]
                        outf.write(str(nudge) + "\n")

                    ## tracer bc paraemters
                    for ii in range(len(tracer_boundary_lst_sorted)):
                        tracer_index = tracer_boundary_lst_sorted[ii]
                        tracer = self.open_boundaries[i]["tracers"][tracer_index]
                        tracer_bc_type = tracer_boundary_types[ii]
                        if tracer_bc_type == 2:
                            tracer_bc_const = tracer["type"]
                            if isinstance(tracer_bc_const, list):
                                for val in tracer_bc_const:
                                    outf.write(str(val) + " ")
                            else:
                                outf.write(str(tracer_bc_const))
                            outf.write("\n")
                        if tracer_bc_type > 0:
                            relax = tracer["relax"]
                            outf.write(str(relax) + "\n ")


if __name__ == "__main__":
    bcyaml = "G:\\schism\\bctides\\bctides.yaml"
    bctides = "G:\\schism\\bctides\\bctides.in"
    by = load_boundary(bcyaml)
    by.write_bctides(bctides)


##################### YAML EXMAPLE #####################################################

# bctides:
#     date: 2011-01-28
#     hgrid_input_file: G:\schism\bctides\hgrid.gr3
#     earth_tidals:
#         tidal_cutoff_depth: 40
#         tidal_constituents:
#           - name: Sa
#             amplitude: 0.05
#             node_factor: 1.0
#             angular_frequency: 0.0027
#             earth_equilibrium_argument: 0.01
#           - name: O1
#             amplitude: 0.2
#             node_factor: 1.1
#             angular_frequency: 4.001
#             earth_equilibrium_argument: 0.0
#           - name: M2
#             amplitude: 0.5
#             node_factor: 0.99
#             angular_frequency: 1.9323
#             earth_equilibrium_argument: 0.0
#     bounary_forcing_tidals:
#         tidal_constituents:
#           - name: Z0
#             angular_frequency: 0.0
#             node_factor: 1
#             earth_equilibrium_argument: 0
#           - name: O1
#             angular_frequency: 0.675977E-04
#             node_factor: 1.11945
#             earth_equilibrium_argument: 7.47302
#           - name: K1
#             angular_frequency: 0.729212E-04
#             node_factor: 1.07391
#             earth_equilibrium_argument: 206.78674
#           - name: M2
#             angular_frequency: 0.140519E-04
#             node_factor: 0.97907
#             earth_equilibrium_argument: 217.04138
#     open_boundaries:
#       - open_boundary:
#         name: ocean
#         elevation_boundary:
#             boundary: tidal
#             tidal_constituents:
#               - name : Z0
#                 amplitude: 1.0
#                 phase: 0.0
#               - name: O1
#                 amplitude: 0.226
#                 phase: 206
#               - name: K1
#                 amplitude: 0.369
#                 phase: 220
#               - name: M2
#                 amplitude: 0.578
#                 phase: 190
#         velocity_boundary:
#             boundary: tidal
#             tidal_constituents:
#               - name: Z0
#                 u_amplitude: 0.20
#                 u_phase: 0.0
#                 v_amplitude: 0.10
#                 v_phase: 0.0
#               - name: O1
#                 u_amplitude: 0.0226
#                 u_phase: 206
#                 v_amplitude: 0.0226
#                 v_phase: 206
#               - name: K1
#                 u_amplitude: 0.0369
#                 u_phase: 220
#                 v_amplitude: 0.0226
#                 v_phase: 206
#               - name: M2
#                 u_amplitude: 0.0578
#                 u_phase: 190
#                 v_amplitude: 0.0226
#                 v_phase: 206

#         temperature_boundary:
#             boundary: initial profile
#             nudge: 0.5
#         salinity_boundary:
#             boundary: initial profile
#             nudge: 0.5
#         tracer_boundaries:
#             boundary: none
#       - open_boundary:
#         name: coyote
#         elevation_boundary:
#             boundary: elev2D.th.nc
#         velocity_boundary:
#             boundary: uv3D.th.nc
#         temperature_boundary:
#             boundary: initial profile
#             nudge: 0.5
#         salinity_boundary:
#             boundary: initial profile
#             nudge: 0.5
#       - open_boundary:
#         name: ccc_rock_slough
#         velocity_boundary:
#             boundary: flux.th
#         temperature_boundary:
#             boundary: initial profile
#             nudge: 1.0
#         salinity_boundary:
#             boundary: initial profile
#             nudge: 1.0
#         tracers:
#           - tracer: GEN
#             number_of_tracers: 2
#             boundary: time history
#             nudge: 0.5
#           - tracer: AGE
#             number_of_tracers: 2
#             boundary: [10, 20]
#             nudge: 1.0
#           - tracer: SED
#             number_of_tracers: 5
#             boundary: [0.1, 0.2, 0.1, 0.12, 0.23]
#             nudge: 0.3
#       - open_boundary:
#         name: ccc_old
#         velocity_boundary:
#             boundary: flux.th
#         temperature_boundary:
#             boundary: initial profile
#             nudge: 0.5
#         salinity_boundary:
#             boundary: initial profile
#             nudge: 0.5
#         tracers:
#           - tracer: GEN
#             number_of_tracers: 1
#             boundary: 10
#             nudge: 0.5
#       - open_boundary:
#         name: swp
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: cvp
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: sjr
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: calaveras
#         elevation_boundary:
#              boundary: 0.3
#         velocity_boundary:
#              boundary: 0.3
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: moke
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: american
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: sac
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: yolo_toedrain
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: yolo
#         velocity_boundary:
#              boundary: flux.th
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: north bay
#         elevation_boundary:
#             boundary: tidal elev2D.th.nc
#             tidal_constituents:
#               - name: Z0
#                 amplitude: 1.0
#                 phase: 0.0
#               - name: O1
#                 amplitude: 0.226
#                 phase: 206
#               - name: K1
#                 amplitude: 0.369
#                 phase: 220
#               - name: M2
#                 amplitude: 0.578
#                 phase: 190
#         velocity_boundary:
#              boundary: tidal uv3D.th.nc
#              tidal_constituents:
#               - name: Z0
#                 u_amplitude: 0.20
#                 u_phase: 0.0
#                 v_amplitude: 0.10
#                 v_phase: 0.0
#               - name: O1
#                 u_amplitude: 0.0226
#                 u_phase: 206
#                 v_amplitude: 0.0226
#                 v_phase: 206
#               - name: K1
#                 u_amplitude: 0.0369
#                 u_phase: 220
#                 v_amplitude: 0.0226
#                 v_phase: 206
#               - name: M2
#                 u_amplitude: 0.0578
#                 u_phase: 190
#                 v_amplitude: 0.0226
#                 v_phase: 206
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: napa
#         velocity_boundary:
#              boundary: flather 1
#              eta_mean: 0.23
#              eta_m0: [0.23,0.23,0.23]
#              vnormal_mean: 0.01
#              vnormal_nodes:
#                - [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.01, 0.01, 0.01]
#                - [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.01, 0.01, 0.01]
#                - [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.01, 0.01, 0.01]
#         temperature_boundary:
#             boundary: 10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#       - open_boundary:
#         name: ccc_vic
#         elevation_boundary:
#             boundary: elev.th
#         velocity_boundary:
#             boundary: flux.th
#         temperature_boundary:
#             boundary:  10.0
#             nudge: 1.0
#         salinity_boundary:
#             boundary: SAL_3D.th.nc
#             nudge: 1.0
#############################################################################################################################
