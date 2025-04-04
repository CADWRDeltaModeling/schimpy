import param

from collections import namedtuple

Params = namedtuple("Params", ["core", "opt", "schout"])


class Core(param.Parameterized):
    ipre = param.Integer(
        default=0, bounds=(0, 1), doc="Pre-processing flag (1: on; 0: off)"
    )
    ibc = param.Integer(default=0, bounds=(0, 1), doc="Barotropic/baroclinic flags")
    ibtp = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="Transport equation flag (0: not solved, 1: solved)",
    )
    rnday = param.Number(
        default=0.0, bounds=(0, None), doc="Total simulation time in days"
    )
    dt = param.Number(default=0.0, bounds=(0, None), doc="Time step in seconds")
    msc2 = param.Integer(
        default=0, bounds=(0, None), doc="Spectral resolution for WWM module"
    )
    mdc2 = param.Integer(
        default=0, bounds=(0, None), doc="Spectral resolution for WWM module"
    )
    eco_class = param.Integer(
        default=0, bounds=(0, None), doc="# of tracer classes for EcoSim module"
    )
    ntracer_gen = param.Integer(
        default=0, bounds=(0, None), doc="# of tracer classes for GEN module"
    )
    ntracer_age = param.Integer(
        default=0, bounds=(0, None), doc="# of tracer classes for AGE module"
    )
    sed_class = param.Integer(
        default=0, bounds=(0, None), doc="# of tracer classes for SED3D module"
    )
    nspool = param.Integer(
        default=0, bounds=(0, None), doc="Output interval for global netcdf outputs"
    )
    ihfskip = param.Integer(
        default=0,
        bounds=(0, None),
        doc="New output stack interval for global netcdf outputs",
    )


class Opt(param.Parameterized):
    gen_wsett = param.Boolean(default=True, doc="Generate wave setup")
    flag_fib = param.Integer(
        default=0, bounds=(0, 1), doc="Fecal Indicator Bacteria model flag"
    )
    ics = param.Integer(
        default=0, bounds=(0, 1), doc="Cohesive sediment transport flag"
    )
    rearth_pole = param.Number(
        default=6356.750, bounds=(0, None), doc="Earth's polar radius [m]"
    )
    rearth_eq = param.Number(
        default=6378.137, bounds=(0, None), doc="Earth's equatorial radius [m]"
    )
    indvel = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="Include non-hydrostatic pressure gradient effects in velocity calculations",
    )
    imm = param.Integer(
        default=0,
        bounds=(0, 2),
        doc="Options for sediment bed morph, 0: no bed deformation, 1: bed deformation, 2: 3D bed deformation",
    )
    ibdef = param.Integer(
        default=10, bounds=(1, 100), doc="# of steps during which deformation occurs"
    )
    ihot = param.Integer(default=0, bounds=(0, 1), doc="Model hot start flag")
    ihydraulics = param.Integer(default=0, bounds=(0, 1), doc="Hydraulics module flag")
    izonal5 = param.Integer(
        default=0, bounds=(0, 1), doc="5th-order horizontal advection scheme flag"
    )
    slam0 = param.Number(
        default=-124,
        bounds=(-360, 360),
        doc="Centers of projection used to convert lat/lon to Cartesian coordinates. These are used if a variable Coriolis parameter is employed (ncor=1)",
    )
    sfea0 = param.Number(
        default=45,
        bounds=(-360, 360),
        doc="Centers of projection used to convert lat/lon to Cartesian coordinates. These are used if a variable Coriolis parameter is employed (ncor=1)",
    )
    iupwind_mom = param.Integer(
        default=0, bounds=(0, 1), doc="Upwinding flag for momentum advection"
    )
    ihorcon = param.Integer(
        default=0, bounds=(0, 1), doc="Horizontal advection of the Coriolis force flag"
    )
    hvis_coef0 = param.Number(
        default=0.1, bounds=(0, None), doc="Horizontal eddy viscosity coefficient"
    )
    ishapiro = param.Integer(default=0, bounds=(0, 1), doc="Shapiro filtering flag")
    shapiro0 = param.Number(
        default=0.8, bounds=(0, 1), doc="Shapiro filtering parameter"
    )
    niter_shap = param.Integer(
        default=0,
        bounds=(0, None),
        doc="Maximum number of iterations for Shapiro filter",
    )
    ihdif = param.Integer(default=0, bounds=(0, 1), doc="Horizontal diffusion flag")
    thetai = param.Number(
        default=0.7, bounds=(0, 1), doc="Angle of horizontal diffusion"
    )
    drampbc = param.Number(
        default=1.0, bounds=(0, None), doc="Damping factor for boundary conditions"
    )
    dramp = param.Number(
        default=1.0, bounds=(0, None), doc="Damping factor for free surface"
    )
    nadv = param.Integer(default=3, bounds=(0, None), doc="Advection scheme flag")
    dtb_min = param.Number(
        default=1.0, bounds=(0, None), doc="Minimum allowable time step size [s]"
    )
    dtb_max = param.Number(
        default=100.0, bounds=(0, None), doc="Maximum allowable time step size [s]"
    )
    h0 = param.Number(
        default=1.0, bounds=(0, None), doc="Initial water depth in meters"
    )
    nchi = param.Integer(
        default=50, bounds=(0, None), doc="Number of horizontal grid cells"
    )
    dzb_min = param.Number(
        default=1e-5, bounds=(0, None), doc="Minimum depth below the free surface"
    )
    hmin_man = param.Number(
        default=0.05, bounds=(0, None), doc="Minimum water depth for nodes"
    )
    ncor = param.Integer(
        default=1, bounds=(0, None), doc="Number of stages in the coriolis integration"
    )
    rlatitude = param.Number(
        default=47.5,
        bounds=(-90, 90),
        doc="Latitude of the origin of the coriolis coordinate system",
    )
    coricoef = param.Number(
        default=1.0, bounds=(0, None), doc="Coriolis coefficient (1/f)"
    )
    nws = param.Integer(default=0, bounds=(0, None), doc="Switch to enable wind")
    wtiminc = param.Number(
        default=3600.0, bounds=(0, None), doc="Time increment for wind forcing"
    )
    iwind_form = param.Integer(
        default=1, bounds=(1, 2), doc="Form of the wind input data"
    )
    drampwind = param.Number(
        default=86400.0, bounds=(0, None), doc="Ramp time for wind forcing"
    )
    iwindoff = param.Integer(
        default=0, bounds=(0, 1), doc="Switch to turn off wind forcing"
    )
    ihconsv = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="Switch to turn on horizontal mass conservation correction",
    )
    isconsv = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="Switch to turn on vertical mass conservation correction",
    )
    itur = param.Integer(default=0, bounds=(0, 3), doc="Turbulence closure model")
    dfv0 = param.Number(
        default=1.0, bounds=(0, None), doc="Vertical eddy viscosity/diffusivity at z=0"
    )
    dfh0 = param.Number(
        default=0.0,
        bounds=(0, None),
        doc="Horizontal eddy viscosity/diffusivity at z=0",
    )
    h1_pp = param.Number(
        default=0.2,
        bounds=(0, None),
        doc="1st parameter for piece-wise linear bottom friction",
    )
    h2_pp = param.Number(
        default=1.0,
        bounds=(0, None),
        doc="2nd parameter for piece-wise linear bottom friction",
    )
    vdmax_pp1 = param.Number(
        default=0.02,
        bounds=(0, None),
        doc="Maximum vertical eddy viscosity for piece-wise linear bottom friction",
    )
    vdmax_pp2 = param.Number(
        default=0.01,
        bounds=(0, None),
        doc="Maximum vertical eddy viscosity for piece-wise linear bottom friction",
    )
    vdmin_pp1 = param.Number(
        default=1e-5,
        bounds=(0, None),
        doc="Minimum vertical eddy viscosity for piece-wise linear bottom friction",
    )
    vdmin_pp2 = param.Number(
        default=1e-5,
        bounds=(0, None),
        doc="Minimum vertical eddy viscosity for piece-wise linear bottom friction",
    )
    tdmin_pp1 = param.Number(
        default=0.01,
        bounds=(0, None),
        doc="Minimum time step size for predictor-corrector scheme for 1st order turbulence closure.",
    )
    tdmin_pp2 = param.Number(
        default=0.001,
        bounds=(0, None),
        doc="Minimum time step size for predictor-corrector scheme for 2nd order turbulence closure.",
    )
    mid = param.Selector(
        objects=["MY", "KL", "KE", "KW", "UB"], default="KL", doc="Turbulence model."
    )
    stab = param.Selector(
        objects=["GA", "KC"], default="GA", doc="Turbulence stability model."
    )
    xlsc0 = param.Number(
        default=1.0,
        bounds=(0, None),
        doc="Initial value for cross-front mixing parameter",
    )
    ibcc_mean = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="Boundary condition for mean flow at solid boundaries",
    )
    flag_ic = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="Flag for initialization of free surface and velocity field at start of simulation",
    )
    start_year = param.Integer(default=2000, bounds=(1900, None), doc="Start year")
    start_month = param.Integer(default=1, bounds=(1, 12), doc="Start month")
    start_day = param.Integer(default=1, bounds=(1, 31), doc="Start day")
    start_hour = param.Integer(default=0, bounds=(0, 23), doc="Start hour")
    utc_start = param.Integer(
        default=8,
        bounds=(-24, 24),
        doc="UTC flag for start time (0: Local time, 1: UTC)",
    )
    itr_met = param.Integer(
        default=3,
        bounds=(1, 3),
        doc="Type of meteorological forcing input (1:WRF, 2:GFS, 3:ECMWF)",
    )
    h_tvd = param.Number(
        default=2.0, bounds=(0, None), doc="TVD limiter parameter for advection scheme"
    )
    eps1_tvd_imp = param.Number(
        default=1e-3,
        bounds=(0, None),
        doc="TVD limiter parameter for implicit vertical advection",
    )
    eps2_tvd_imp = param.Number(
        default=1e-4,
        bounds=(0, None),
        doc="TVD limiter parameter for implicit vertical advection",
    )
    ip_weno = param.Integer(
        default=0,
        bounds=(0, 3),
        doc="Order of WENO scheme for horizontal advection order of accuracy: 0- upwind; 1- linear polynomial, 2nd order; 2- quadratic polynomial, 3rd order",
    )
    courant_weno = param.Number(
        default=0.8,
        bounds=(0, 1),
        doc="Courant number limit for WENO scheme for horizontal advection",
    )
    ntd_weno = param.Integer(
        default=1,
        bounds=(1, None),
        doc="Number of time derivatives for WENO scheme for horizontal advection",
    )
    nquad = param.Integer(
        default=5, bounds=(1, None), doc="Order of quadrature for vertical advection"
    )
    epsilon1 = param.Number(
        default=1e-5,
        bounds=(0, None),
        doc="Criteria for relaxation of vertical velocity in the CFL limiter",
    )
    epsilon2 = param.Number(
        default=1e-3,
        bounds=(0, None),
        doc="Criteria for relaxation of vertical velocity in the CFL limiter",
    )
    epsilon3 = param.Number(
        default=1e-4,
        bounds=(0, None),
        doc="Criteria for relaxation of vertical velocity in the CFL limiter",
    )
    ielad_weno = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="Switch for ENO/WENO advection scheme for elevations ielad, if ielad=1, use ELAD method to suppress dispersion (inactive at the moment)",
    )
    small_elad = param.Number(
        default=1e-6,
        bounds=(1e-8, 1e-3),
        doc="Small number to prevent negative water depths",
    )
    i_prtnftl_weno = param.Integer(
        default=1, bounds=(0, 2), doc="Switch for partitions for WENO advection scheme"
    )
    inu_tr = param.Integer(
        default=0, bounds=(0, 4), doc="Switch for turbulence closure scheme"
    )
    step_nu_tr = param.Number(
        default=0.001, bounds=(0, None), doc="Time step for the turbulence closure"
    )
    vnh1 = param.Number(
        default=1, bounds=(0, None), doc="Horizontal viscosity coefficient for momentum"
    )
    vnh2 = param.Number(
        default=0,
        bounds=(0, None),
        doc="Horizontal viscosity coefficient for turbulent quantities",
    )
    vnf1 = param.Number(
        default=1, bounds=(0, None), doc="Vertical viscosity coefficient for momentum"
    )
    vnf2 = param.Number(
        default=1,
        bounds=(0, None),
        doc="Vertical viscosity coefficient for turbulent quantities",
    )
    moitn0 = param.Number(
        default=0.1, bounds=(0, None), doc="Time step for sea ice modeling"
    )
    mxitn0 = param.Integer(
        default=5, bounds=(0, None), doc="Number of subcycles for sea ice modeling"
    )
    rtol0 = param.Number(
        default=1e-6,
        bounds=(1e-12, 1e-2),
        doc="Relative tolerance for iterative solvers",
    )
    iflux = param.Integer(default=2, bounds=(0, 2), doc="Switch for flux type")
    iflux_out_format = param.Integer(
        default=0, bounds=(0, 2), doc="Format for outputting fluxes"
    )
    inter_mom = param.Integer(
        default=1, bounds=(0, 2), doc="Interpolation for momentum advection"
    )
    h_bcc1 = param.Number(
        default=1.5, bounds=(0, None), doc="Constant for BCC turbulence closure scheme"
    )
    inu_elev = param.Integer(
        default=1, bounds=(0, 2), doc="Switch for viscosity in elevation advection"
    )
    inu_uv = param.Integer(
        default=1, bounds=(0, 2), doc="Switch for viscosity in velocity advection"
    )
    ihhat = param.Integer(
        default=0, bounds=(0, 1), doc="Switch for parameterization of mean circulation"
    )
    kr_co = param.Number(
        default=0.5, bounds=(0, None), doc="Coefficient for vertical mixing"
    )
    rmaxvel = param.Number(
        default=2.0, bounds=(1, 10), doc="Maximum value of the absolute velocity"
    )
    velmin_btrack = param.Number(
        default=1e-4, bounds=(1e-8, 1e-2), doc="Minimum velocity for bottom tracking"
    )
    btrack_nudge = param.Number(
        default=9.013e-3, bounds=(0, 0.5), doc="Switch for nudging of bottom tracking"
    )
    ibtrack_test = param.Integer(
        default=1, bounds=(0, 1), doc="Switch for testing bottom tracking"
    )
    irouse_test = param.Integer(
        default=0, bounds=(0, 1), doc="Switch for testing Rouse number"
    )
    inunfl = param.Boolean(True, doc="If .true., the inundation from ocean is included")
    shorewafo = param.Number(
        default=0.0, bounds=(0, None), doc="Shoreline wave friction coefficient"
    )
    ic_elev = param.Number(
        default=0.0, bounds=(None, None), doc="Initial condition of elevation (m)"
    )
    nramp_elev = param.Number(
        default=0.0, bounds=(None, None), doc="Duration of ramping for ic_elev (s)"
    )
    inv_atm_bnd = param.Boolean(
        False, doc="If .true., atmospheric pressure is forced at the boundary"
    )
    prmsl_ref = param.Number(
        default=1013.25,
        bounds=(None, None),
        doc="Reference pressure for atmospheric pressure (Pa)",
    )
    s1_mxnbt = param.Number(
        default=0.0,
        bounds=(0, None),
        doc="s1 * maximum boundary stress for wind stress calculation",
    )
    s2_mxnbt = param.Number(
        default=0.0,
        bounds=(0, None),
        doc="s2 * maximum boundary stress for wind stress calculation",
    )
    iharind = param.Integer(default=1, bounds=(0, 1), doc="Index of harmonic analysis")
    icou_elfe_wwm = param.Integer(
        default=0, bounds=(0, 1), doc="Coupling option between elfe and WWM models"
    )
    drampwafo = param.Number(
        default=0.0, bounds=(None, None), doc="Duration of ramping for shorewafo (s)"
    )
    nstep_wwm = param.Integer(
        default=1, bounds=(1, None), doc="Number of steps between WWM model calls"
    )
    hmin_radstress = param.Number(
        default=0.1,
        bounds=(0, None),
        doc="Minimum depth for radiation stress calculation (m)",
    )
    turbinj = param.Boolean(False, doc="If .true., the turbine injects turbulence")
    turbinjds = param.Number(
        default=0.0, bounds=(0, None), doc="Turbine injection diameter size (m)"
    )
    alphaw = param.Number(
        default=0.0, bounds=(None, None), doc="WWM swell attenuation coefficient"
    )
    fwvor_advxy_stokes = param.Boolean(
        False, doc="If .true., advective terms in the Stokes drift are calculated"
    )
    fwvor_advz_stokes = param.Boolean(
        False,
        doc="If .true., vertical advective terms in the Stokes drift are calculated",
    )
    fwvor_gradpress = param.Boolean(
        False,
        doc="If .true., the pressure gradient force in the FWVOR model is calculated",
    )
    fwvor_breaking = param.Boolean(
        False,
        doc="If .true., the breaking parameterization in the FWVOR model is turned on",
    )
    fwvor_streaming = param.Boolean(
        False,
        doc="If .true., the streaming parameterization in the FWVOR model is turned on",
    )
    fwvor_wveg = param.Boolean(
        False,
        doc="If .true., the WVEG wave generation term is calculated in the FWVOR model",
    )
    fwvor_wveg_NL = param.Boolean(
        False, doc="If .true., the nonlinear WVEG term is calculated in the FWVOR model"
    )
    wafo_obcramp = param.Number(
        1, bounds=(0, None), doc="WAFO parameter for offshore boundary ramp"
    )
    iwbl = param.Integer(
        0, bounds=(0, 2), doc="WBL model flag (0=Off, 1=JK model, 2=TW model)"
    )
    cur_wwm = param.Boolean(False, doc="Use wave-current interaction?")
    if_source = param.Integer(
        2, bounds=(0, 2), doc="Infragravity wave source (0=none, 1=roller, 2=eddy)"
    )
    dramp_ss = param.Number(
        60.0, bounds=(0, None), doc="Ramp time for sea surface forcing (hours)"
    )
    ieos_type = param.Integer(
        0,
        bounds=(0, 1),
        doc="Equation of state 0: original UNESCO; 1: linear function of T ONLY",
    )
    ieos_pres = param.Number(
        0.0, bounds=(0, None), doc="Reference pressure for nonlinear EOS (dbar)"
    )
    eos_a = param.Number(-0.1, bounds=(None, 0), doc="EOS coefficient A")
    eos_b = param.Number(1001.0, bounds=(None, None), doc="EOS coefficient B")
    slr_rate = param.Number(0.0, bounds=(0, None), doc="Sea level rise rate (mm/year)")
    rho0 = param.Number(1000.0, bounds=(0, None), doc="Reference density (kg/m^3)")
    shw = param.Number(0.025, bounds=(0, None), doc="Surface wave height (m)")
    isav = param.Integer(0, bounds=(0, 1), doc="Save hotstart files (0=no, 1=yes)")
    nstep_ice = param.Integer(
        0, bounds=(0, None), doc="Number of time steps per ice cycle"
    )
    iunder_deep = param.Integer(
        0, bounds=(0, 1), doc="Allow subgrid depth < hmin (0=no, 1=yes)"
    )
    h1_bcc = param.Number(1.0, bounds=(0, None), doc="Boundary condition parameter 1")
    h2_bcc = param.Number(0.5, bounds=(0, None), doc="Boundary condition parameter 2")
    hw_depth = param.Number(15.0, bounds=(0, None), doc="Wave depth parameter (m)")
    hw_ratio = param.Number(1.2, bounds=(0, None), doc="Wave height/depth ratio")
    level_age = param.List(
        item_type=float,
        default=[-999.0],
        bounds=(-999.0, None),
        doc="Minimum age of a level (days)",
    )
    vclose_surf_frac = param.Number(
        0.0,
        bounds=(0, 1),
        doc="Fraction of initial cell depth to set vclose near surface",
    )
    iadjust_mass_consv0 = param.Integer(
        1, bounds=(0, 1), doc="Adjust mass conservation (0=no, 1=yes)"
    )
    ipre2 = param.Integer(
        0, bounds=(0, 1), doc="Solve pressure-Poisson equation? (0=no, 1=yes)"
    )
    ielm_transport = param.Integer(
        default=0,
        bounds=(0, 1),
        doc="""Hybrid ELM-FV transport for performance; used only with itr_met>=3. 
    If ielm_transport=1, the hybrid scheme is invoked and max_subcyc represents the max # of subcycling per time step in transport allowed; 
    if the actual # of subcycling in a prism at a time step exceeds this threshold, more efficient ELM transport is used locally 
    (at the expense of mass conservation, so make sure this option is used sparingly).""",
    )
    max_subcyc = param.Integer(
        default=10, bounds=(1, None), doc="Maximum number of subcycles"
    )
    i_hmin_airsea_ex = param.Integer(
        default=1,
        bounds=(1, 2),
        doc="Enforce minimum water depth for air-sea exchange? (0: no, 1: yes) Option to locally turn off heat exchange.",
    )
    hmin_airsea_ex = param.Number(
        default=0.01, bounds=(0, None), doc="Minimum depth (m) for air-sea exchange"
    )
    itransport_only = param.Integer(
        default=0, bounds=(0, 1), doc="Whether to do only transport"
    )
    meth_sink = param.Integer(
        default=2, bounds=(1, 2), doc="Method for resolving horizontal mass divergence"
    )
    iloadtide = param.Integer(
        default=0, bounds=(0, 1), doc="Whether to include tidal forcing"
    )
    loadtide_coef = param.Number(
        default=1.0, bounds=(0, None), doc="Scaling factor for tidal forcing"
    )
    nu_sum_mult = param.Number(
        default=1.0, bounds=(0, None), doc="Multiplier for horizontal viscosity"
    )
    i_hmin_salt_ex = param.Integer(
        default=1,
        bounds=(0, 2),
        doc="Enforce minimum water depth for salt flux computation? (0: no, 1: yes)",
    )
    hmin_salt_ex = param.Number(
        default=0.01,
        bounds=(0, None),
        doc="Minimum depth (m) for salt flux computation",
    )
    h_massconsv = param.Number(
        default=0.001,
        bounds=(0, None),
        doc="Threshold for enforcing mass conservation in meters",
    )
    lev_tr_source = param.Integer(
        default=1, bounds=(1, None), doc="Level of tracers (1: constant, 2: linear)"
    )
    rinflation_icm = param.Number(
        default=2.0,
        bounds=(0, None),
        doc="Ratio of inflation coefficient in ICM to target coefficient",
    )
    iprecip_off_bnd = param.Integer(
        default=1,
        bounds=(0, 1),
        doc="Turn off precipitation in dry cells? (0: no, 1: yes)",
    )


class OutControls(param.Parameterized):
    """Parameters controlling output"""

    def get_iof_array_names(self):
        out_params = self.param.params().values()
        out_var_names = [param.name for param in out_params]
        return out_var_names[1:]  # first is always name

    def to_iof_array(self):
        out_params = self.param.params()
        valuemap = dict(self.param.get_param_values())
        iof_array = []
        for i, name in enumerate(out_params):
            if i == 0:  # skip
                pass
            else:
                iof_array.append(int(valuemap[name]))
        return iof_array

    def from_iof_array(self, iof_array):
        bool_values = [bool(val) for val in iof_array]
        # Enumerate the variables and create a param_dict
        param_dict = {
            param_name: bool_value
            for param_name, bool_value in zip(self.param.params(), [0] + bool_values)
        }
        self.set_param(**param_dict)


class HydroOutput(OutControls):
    """Class to specify which variables to output from the SCHISM model"""

    dry_flag_node = param.Boolean(
        True, doc="Whether to output dry flag nodes.", constant=True
    )
    elevation = param.Boolean(False, doc="Whether to output water surface elevation.")
    air_pressure = param.Boolean(False, doc="Whether to output air pressure.")
    air_temperature = param.Boolean(False, doc="Whether to output air temperature.")
    specific_humidity = param.Boolean(False, doc="Whether to output specific humidity.")
    solar_radiation = param.Boolean(
        False, doc="Whether to output solar radiation."
    )  # 5
    sensible_heat = param.Boolean(False, doc="Whether to output sensible heat.")
    latent_heat = param.Boolean(False, doc="Whether to output latent heat.")
    upward_longwave = param.Boolean(
        False, doc="Whether to output upward longwave radiation."
    )
    downward_longwave = param.Boolean(
        False, doc="Whether to output downward longwave radiation."
    )
    total_heat = param.Boolean(False, doc="Whether to output total heat flux.")  # 10
    evaporation_rate = param.Boolean(False, doc="Whether to output evaporation rate.")
    precipitation_rate = param.Boolean(
        False, doc="Whether to output precipitation rate."
    )
    bottom_stress_xy = param.Boolean(
        False, doc="Whether to output bottom stress in x,y direction."
    )
    wind_speed_xy = param.Boolean(
        False, doc="Whether to output wind speed in x,y direction."
    )
    wind_stress_xy = param.Boolean(
        False, doc="Whether to output wind stress in x,y direction."
    )  # 15
    depth_avg_vel_xy = param.Boolean(
        False, doc="Whether to output depth-averaged velocity in x,y direction."
    )
    vertical_velocity = param.Boolean(False, doc="Whether to output vertical velocity.")
    temperature = param.Boolean(False, doc="Whether to output water temperature.")
    salinity = param.Boolean(False, doc="Whether to output water salinity.")
    water_density = param.Boolean(False, doc="Whether to output water density.")  # 20
    diffusivity = param.Boolean(
        False, doc="Whether to output vertical eddy diffusivity."
    )
    viscosity = param.Boolean(False, doc="Whether to output vertical eddy viscosity.")
    turbulent_kinetic_ener = param.Boolean(
        False, doc="Whether to output turbulent kinetic energy."
    )
    mixing_length = param.Boolean(False, doc="Whether to output mixing length.")
    z_coordinates = param.Boolean(False, doc="Whether to output z coordinates.")  # 25
    horizontal_vel_xy = param.Boolean(
        False, doc="Whether to output horizontal velocity in x,y direction."
    )
    horizontal_side_vel_xy = param.Boolean(
        False, doc="Whether to output horizontal side velocity in x,y direction."
    )
    vertical_vel_at_element = param.Boolean(
        False, doc="Whether to output vertical velocity at element."
    )
    temperature_at_element = param.Boolean(
        False, doc="Whether to output water temperature at element."
    )
    salinity_at_element = param.Boolean(
        False, doc="Whether to output salinity at element."
    )  # 30
    barotropic_pres_grad_xy = param.Boolean(
        False, doc="Whether to output barotropic pressure gradient in x,y direction."
    )


class WindWaveOutput(OutControls):
    """A collection of parameters that define the wind wave outputs."""

    sig_wave_height = param.Boolean(
        default=True, doc="Flag to output significant wave height"
    )
    mean_wave_period = param.Boolean(
        default=True, doc="Flag to output mean wave period"
    )
    zero_downcross_period = param.Boolean(
        default=True, doc="Flag to output zero downcrossing period"
    )
    tm_10 = param.Boolean(
        default=True, doc="Flag to output spectral moments of order 10"
    )
    mean_wave_number = param.Boolean(
        default=True, doc="Flag to output mean wave number"
    )
    mean_wave_length = param.Boolean(
        default=True, doc="Flag to output mean wave length"
    )
    mean_wave_direction = param.Boolean(
        default=True, doc="Flag to output mean wave direction"
    )
    mean_dir_spreading = param.Boolean(
        default=True, doc="Flag to output mean directional spreading"
    )
    peak_period = param.Boolean(default=True, doc="Flag to output peak wave period")
    continuous_peak_period = param.Boolean(
        default=True, doc="Flag to output continuous peak wave period"
    )
    peak_phase_vel = param.Boolean(
        default=True, doc="Flag to output peak phase velocity"
    )
    peak_n_factor = param.Boolean(default=True, doc="Flag to output peak N-factor")
    peak_group_vel = param.Boolean(
        default=True, doc="Flag to output peak group velocity"
    )
    peak_wave_number = param.Boolean(
        default=True, doc="Flag to output peak wave number"
    )
    peak_wave_length = param.Boolean(
        default=True, doc="Flag to output peak wave length"
    )
    dominant_direction = param.Boolean(
        default=True, doc="Flag to output dominant wave direction"
    )
    peak_spreading = param.Boolean(
        default=True, doc="Flag to output peak directional spreading"
    )
    discrete_peak_direction = param.Boolean(
        default=True, doc="Flag to output discrete peak direction"
    )
    orbital_velocity = param.Boolean(
        default=True, doc="Flag to output orbital velocity"
    )
    rms_orbital_velocity = param.Boolean(
        default=True, doc="Flag to output root mean square of orbital velocity"
    )
    bottom_excursion_period = param.Boolean(
        default=True, doc="Flag to output bottom excursion period"
    )
    bottom_wave_period = param.Boolean(
        default=True, doc="Flag to output bottom wave period"
    )
    u_resell_number = param.Boolean(default=True, doc="Flag to output U-Resell number")
    frictional_velocity = param.Boolean(
        default=True, doc="Flag to output frictional velocity"
    )
    charnock_coeff = param.Boolean(
        default=True, doc="Flag to output Charnock coefficient"
    )
    roughness_length = param.Boolean(
        default=True, doc="Flag to output roughness length"
    )
    roller_diss_rate = param.Boolean(
        default=True, doc="Flag to output roller dissipation rate"
    )
    diss_rate_dep_breaking = param.Boolean(
        default=True,
        doc="Flag to output dissipation rate due to depth-induced breaking",
    )
    diss_rate_bott_friction = param.Boolean(
        default=True, doc="Flag to output dissipation rate due to bottom friction"
    )
    diss_rate_whitecapping = param.Boolean(
        default=True, doc="Flag to output dissipation rate due to whitecapping"
    )
    diss_rate_vegetation = param.Boolean(
        default=True, doc="Flag to output dissipation rate due to vegetation"
    )
    energy_input_atmos = param.Boolean(
        default=True, doc="Flag to output energy input from atmosphere"
    )
    wave_energy_dir_x = param.Boolean(
        default=True, doc="Flag to output wave energy flux direction in X direction"
    )
    wave_energy_dir_y = param.Boolean(
        default=True, doc="Flag to output wave energy flux direction in Y direction"
    )


class TracerGenOutput(param.Parameterized):
    ntracer_gen = param.Integer(default=0, doc="Number of Tracer gens")

    def get_tracer_names(self):
        return [f"GEN_{i+1}" for i in range(self.ntracer_gen)]


class Schout(param.Parameterized):
    nc_out = param.Boolean(default=True, doc="If True, netCDF output is enabled.")
    nhot = param.Integer(
        default=0, bounds=(0, None), doc="Number of hot-start files to output."
    )
    nhot_write = param.Integer(
        default=0, bounds=(0, None), doc="Number of hot-start files to write."
    )
    iout_sta = param.Integer(
        default=0, bounds=(0, None), doc="Output frequency for STA files."
    )
    nspool_sta = param.Integer(
        default=0,
        bounds=(0, None),
        doc="Number of time steps between writes to STA files.",
    )
    iof_ugrid = param.Boolean(
        default=False, doc="If True, output unstructured grid variables."
    )
    iof_hydro = HydroOutput()
    iof_wwm = WindWaveOutput()
    iof_gen = (
        TracerGenOutput()
    )  # has equal to Core.ntracer_gen files if USE_GEN is true
    iof_age = param.Boolean(
        default=False, doc="If True, output age variables."
    )  # if USE_AGE is true
    iof_sed = param.Boolean(default=False, doc="If True, output sediment variables.")
    iof_eco = param.Boolean(default=False, doc="If True, output ecosystem variables.")
    iof_icm_core = param.Boolean(
        default=False, doc="If True, output ICM core variables."
    )
    iof_icm_silica = param.Boolean(
        default=False, doc="If True, output ICM silica variables."
    )
    iof_icm_zb = param.Boolean(
        default=False, doc="If True, output ICM zoobenthos variables."
    )
    iof_icm_ph = param.Boolean(default=False, doc="If True, output ICM pH variables.")
    iof_icm_cbp = param.Boolean(
        default=False,
        doc="If True, output ICM carbon, biomass, and production variables.",
    )
    iof_icm_sav = param.Boolean(default=False, doc="If True, output ICM SAV variables.")
    iof_icm_veg = param.Boolean(
        default=False, doc="If True, output ICM vegetation variables."
    )
    iof_icm_sed = param.Boolean(
        default=False, doc="If True, output ICM sediment variables."
    )
    iof_icm_ba = param.Boolean(
        default=False, doc="If True, output ICM blue crab variables."
    )
    iof_icm_dbg = param.Boolean(
        default=False, doc="If True, output ICM debugging variables."
    )
    iof_cos = param.Boolean(default=False, doc="If True, output COS variables.")
    iof_fib = param.Boolean(default=False, doc="If True, output FIB variables.")
    iof_sed2d = param.Boolean(
        default=False, doc="If True, output 2D sediment variables."
    )
    iof_ice = param.Boolean(default=False, doc="If True, output ice variables.")
    iof_ana = param.Boolean(default=False, doc="If True, output analytical fields.")
    iof_marsh = param.Boolean(default=False, doc="If True, output marsh variables.")
    iof_dvd = param.Boolean(default=False, doc="If True, output DVD variables.")


def get_type_dict(cls):
    return {
        k: v["type"] if "type" in v else "str" for k, v in cls.param.schema().items()
    }


def coerce_to_type(schema, value):
    typename = "str"
    if "type" in schema:
        typename = schema["type"]
    if typename == "integer":
        return int(value)
    elif typename == "number":
        return float(
            value.replace("d", "e")
        )  # added replacement for fortran values of type 1000.d0
    elif typename == "boolean":
        return bool(value)
    elif typename == "str":
        return value.replace("'", "")  # fortran code has string single quotes
    elif typename == "array":
        values = value.split(",")
        tvalues = []
        for v in values:
            if "items" in schema:
                tv = coerce_to_type(schema["items"], v)
            else:
                tv = coerce_to_type({"type": "str"}, v)
            tvalues.append(tv)
        return tvalues
    else:
        raise BaseException(f"Unknown typename {typename}")


def get_value_dict(map, cls):
    """Create a dict from namelist from name: value elements and coerce to the type in the schema"""
    xdict = {key: map[key]["value"] for key in map if key != "comment"}
    tdict = cls.param.schema()  # get_type_dict(cls)
    xkdict = {}  # keep only those defined in the class
    for p in cls.param:
        if p in xdict:
            try:
                xkdict[p] = coerce_to_type(tdict[p], xdict[p])
            except Exception as err:
                print(f"Invalid value for {p} : {xdict[p]}")
    return xkdict


def create_params(namelists):
    params = []
    for section, cls in {"CORE": Core, "OPT": Opt, "SCHOUT": Schout}.items():
        inst = cls.__new__(cls)
        if section in namelists:
            inst.__init__(**get_value_dict(namelists[section], cls))
        params.append(inst)
    p = Params(*tuple(params))
    p.schout.iof_hydro = HydroOutput(
        **build_iof_dict("iof_hydro", namelists["SCHOUT"], HydroOutput)
    )
    return p


def build_iof_dict(iof_name, vdict, cls):
    """Builds dictionary of names from the class cls and the values
    that have the form  iof_name(n) = 0 or 1
    where n is an integer and the values are 0 or 1 which map onto False (0) or True (1)
    """
    vars = list(cls.param.schema().keys())[1:]  # first one is always name
    position_values = {
        key.split("(")[1].split(")")[0]: vdict[key]["value"]
        for key in vdict
        if key.startswith(iof_name)
    }
    return {vars[int(k)]: bool(int(v)) for k, v in position_values.items()}
