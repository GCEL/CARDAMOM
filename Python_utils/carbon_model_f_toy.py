import numpy as np

# explicit publics
CARBON_MODEL = None
soil_frac_clay = None
soil_frac_sand = None
nos_soil_layers = None
extracted_C = None
cica_time = None
gs_demand_supply_ratio = None
gs_total_canopy = None
gb_total_canopy = None
canopy_par_MJday_time = None
dim_1 = None
dim_2 = None
nos_trees = None
nos_inputs = None
leftDaughter = None
rightDaughter = None
nodestatus = None
xbestsplit = None
nodepred = None
bestvar = None

# Parameters
# useful technical parameters
xacc = 1e-4  # accuracy parameter for zbrent bisection proceedure
vsmall = np.finfo(float).tiny * 1e3

nos_root_layers = 2
nos_soil_layers = nos_root_layers + 1
pi = 3.1415927
pi_1 = 0.3183099  # pi**(-1)
two_pi = 6.283185  # pi*2
deg_to_rad = 0.01745329  # pi/180
sin_dayl_deg_to_rad = 0.3979486  # sin(23.45 * deg_to_rad)
boltz = 5.670400e-8  # Boltzmann constant (W.m-2.K-4)
emissivity = 0.96
emiss_boltz = 5.443584e-08  # emissivity * boltz
sw_par_fraction = 0.5  # fraction of short-wave radiation which is PAR
freeze = 273.15
gs_H2O_CO2 = 1.646259  # The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
gs_H2Ommol_CO2mol_day = 142.2368  # The ratio of H20:CO2 diffusion for gs, including seconds per day correction
gb_H2O_CO2 = 1.37  # The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
mmol_to_kg_water = 1.8e-5  # milli mole conversion to kg
umol_to_gC = 1.2e-5  # conversion of umolC -> gC
gC_to_umol = 83333.33  # conversion of gC -> umolC; umol_to_gC**(-1)
Rcon = 8.3144  # Universal gas constant (J.K-1.mol-1)
vonkarman = 0.41  # von Karman's constant
vonkarman_1 = 2.439024  # 1 / von Karman's constant
cpair = 1004.6  # Specific heat capacity of air; used in energy balance J.kg-1.K-1

# photosynthesis / respiration parameters
kc_half_sat_25C = 310  # CO2 half saturation, saturation value
kc_half_sat_gradient = 23.956  # CO2 half sat, half sat
co2comp_sat_25C = 36.5  # CO2 compensation point, saturation
co2comp_gradient = 9.46  # CO2 comp point, half sat
# Each of these are temperature sensitivity

# hydraulic parameters
gplant = 4  # plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
root_resist = 25  # Root resistivity (MPa s g mmol−1 H2O)
max_depth = 2  # max root depth (m)
root_k = 100  # root biomass needed to reach 50% depth (gbiomass/m2)
root_radius = 0.00029  # root radius (m) Bonen et al 2014 = 0.00029
root_radius_1 = 1 / root_radius
root_cross_sec_area = np.pi * root_radius * root_radius  # root cross-sectional area (m2)
root_density = 0.31e6  # root density (g biomass m-3 root)
# 0.5e6 Williams et al 1996
# 0.31e6 Bonan et al 2014
root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-1)
const_sfc_pressure = 101325  # (Pa)  Atmospheric surface pressure
head = 0.009807  # head of pressure (MPa/m)
head_1 = 101.968  # inverse head of pressure (m/MPa)

# structural parameters
canopy_height = 9  # canopy height assumed to be 9 m
tower_height = canopy_height + 2  # tower (observation) height assumed to be 2 m above canopy
min_wind = 0.2  # minimum wind speed at canopy top
min_layer = 0.03  # minimum thickness of the third rooting layer (m)
top_soil_depth = 0.30  # thickness of the topsoil layer (m)
min_root = 5  # minimum root biomass (gBiomass.m-2)
min_lai = 0.1  # minimum LAI assumed for aerodynamic conductance calculations (m2/m2)

# timing parameters
seconds_per_hour = 3600  # Number of seconds per hour
seconds_per_day = 86400  # Number of seconds per day
seconds_per_day_1 = 1.157407e-05  # Inverse of seconds per day

# ACM-GPP-ET parameters
pn_max_temp = 6.842942e+01  # Maximum daily max temperature for photosynthesis (oC)
pn_min_temp = -1e+06  # Minimum daily max temperature for photosynthesis (oC)
pn_opt_temp = 3.155960e+01  # Optimum daily max temperature for photosynthesis (oC)
pn_kurtosis = 1.889026e-01  # Kurtosis of photosynthesis temperature response
# bespoke
# pn_max_temp = 59  # Maximum daily max temperature for photosynthesis (oC)
# pn_min_temp = -4  # Minimum daily max temperature for photosynthesis (oC)
# pn_opt_temp = 30  # Optimum daily max temperature for photosynthesis (oC)
# pn_kurtosis = 0.07  # Kurtosis of photosynthesis temperature response
e0 = 3.661204  # Quantum yield gC/MJ/m2/day PAR
minlwp_default = -1.808224  # minimum leaf water potential (MPa)
soil_iso_to_net_coef_LAI = -2.717467  # Coefficient relating soil isothermal net radiation to net.
iWUE = 6.431150e-03  # Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
soil_swrad_absorption = 9.989852e-01  # Fraction of SW rad absorbed by soil
max_lai_lwrad_release = 9.516639e-01  # 1-Max fraction of LW emitted from the canopy to be released
lai_half_lwrad_release = 4.693329e+00  # LAI at which LW emitted from the canopy to be released at 50%
soil_iso_to_net_coef_SW = -3.500964e-02  # Coefficient relating soil isothermal net radiation to net.
soil_iso_to_net_const = 3.455772e+00  # Constant relating soil isothermal net radiation to net
max_par_transmitted = 1.628077e-01  # Max fraction of canopy incident PAR transmitted to soil
max_nir_transmitted = 2.793660e-01  # Max fraction of canopy incident NIR transmitted to soil
max_par_reflected = 1.629133e-01  # Max fraction of canopy incident PAR reflected to the sky
max_nir_reflected = 4.284365e-01  # Max fraction of canopy incident NIR reflected to the sky
canopy_iso_to_net_coef_SW = 1.480105e-02  # Coefficient relating SW to the adjustment between isothermal and net LW
canopy_iso_to_net_const = 3.753067e-03  # Constant relating canopy isothermal net radiation to net
canopy_iso_to_net_coef_LAI = 2.455582e+00  # Coefficient relating LAI to the adjustment between isothermal and net LW

minlwp = minlwp_default

# forest rotation specific info
extracted_C = np.empty(0)
# Metrics on photosynthetic activity
gs_demand_supply_ratio = np.empty(0)
gs_total_canopy = np.empty(0)
gb_total_canopy = np.empty(0)
cica_time = np.empty(0)
canopy_par_MJday_time = np.empty(0)

# arrays for the emulator
dim_1 = None
dim_2 = None
nos_trees = None
nos_inputs = None
leftDaughter = np.empty((0, 0))
rightDaughter = np.empty((0, 0))
nodestatus = np.empty(0)
xbestsplit = np.empty(0)
nodepred = np.empty(0)
bestvar = np.empty(0)

# hydraulic model variables
water_retention_pass = None
soil_layer = None
layer_thickness = np.zeros(nos_soil_layers + 1)
soil_frac_clay = np.zeros(nos_soil_layers)
soil_frac_sand = np.zeros(nos_soil_layers)
uptake_fraction = np.zeros(nos_root_layers)
demand = np.zeros(nos_root_layers)
water_flux_mmolH2Om2s = None

root_reach = None
root_biomass = None
total_water_flux = None
canopy_wind = None
ustar = None
ustar_Uh = None
air_density_kg = None
ET_demand_coef = None
roughl = None
displacement = None
max_supply = None
meant = None
canopy_swrad_MJday = None
canopy_par_MJday = None
soil_swrad_MJday = None
canopy_lwrad_Wm2 = None
soil_lwrad_Wm2 = None
sky_lwrad_Wm2 = None
stomatal_conductance = None
potential_conductance = None
minimum_conductance = None
aerodynamic_conductance = None
convert_ms1_mol_1 = None
lambda_ = None
psych = None
slope = None
water_vapour_diffusion = None
dynamic_viscosity = None
kinematic_viscosity = None

# Module level variables for ACM_GPP_ET parameters
delta_gs = None
Ceff = None
iWUE_step = None
metabolic_limited_photosynthesis = None
light_limited_photosynthesis = None
ci = None
gb_mol = None
rb_mol_1 = None
pn_airt_scaling = None
co2_half_sat = None
co2_comp_point = None

# Module level variables for step-specific met drivers
mint = None
maxt = None
leafT = None
swrad = None
co2 = None
doy = None
wind_spd = None
vpd_kpa = None
lai_1 = None
lai = None

# Module level variables for step-specific timing information
cos_solar_zenith_angle = None
seconds_per_step = None
days_per_step = None
dayl_hours_fraction = None
dayl_seconds = None
dayl_seconds_1 = None
dayl_hours = None

def zbrent(called_from, func, x1, x2, tol, toltol):
    """
    This is a bisection routine. When ZBRENT is called, we provide a reference to a particular function
    and also two values which bound the arguments for the function of interest. ZBRENT finds a root of
    the function (i.e. the point where the function equals zero), that lies between the two bounds.
    There are five exit conditions:
    1) The first proposal for the root of the function equals zero
    2) The proposal range has been reduced to less than tol
    3) The magnitude of the function is less than toltol
    4) The maximum number of iterations has been reached
    5) The root of the function does not lie between supplied bounds
    For a full description, see Press et al. (1986).
    """
    # local variables
    iter = 0
    ITMAX = 8
    EPS = 6e-8

    # calculations
    a = x1
    b = x2
    fa = func(a)
    fb = func(b)
    tol0 = tol * 0.5

    # Check that we haven't (by fluke) already started with the root...
    if abs(fa) < toltol:
        return a
    elif abs(fb) < toltol:
        return b

    c = b
    fc = fb

    for iter in range(1, ITMAX + 1):
        # If the new value (f(c)) doesn't bracket the root with f(b) then adjust it
        if fb * fc > 0.0:
            c = a
            fc = fa
            d = b - a
            e = d

        if abs(fc) < abs(fb):
            a, b, c = b, c, a
            fa, fb, fc = fb, fc, fa

        tol1 = EPS * abs(b) + tol0
        xm = 0.5 * (c - b)

        if abs(xm) <= tol1 or abs(fb) < toltol:
            return b

        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)

            if p > 0.0:
                q = -q

            p = abs(p)

            if 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q)):
                e = d
                d = p / q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d

        a = b
        fa = fb

        if abs(d) > tol1:
            b += d
        else:
            b += sign(tol1, xm)

        fb = func(b)

    return b


def sign(a, b):
    """Return the sign of b with the sign of a."""
    return abs(a) * (1 if b >= 0 else -1)


def meteorological_constants(input_temperature, input_temperature_K, input_vpd_kPa):
    """
    Determine some multiple-use constants used by a wide range of functions.
    All variables here are linked to air temperature and thus invariant between
    iterations and can be stored in memory.
    """
    # Density of air (kg/m3)
    air_density_kg = 353 / input_temperature_K
    # Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / (input_temperature_K * Rcon)
    # Latent heat of vaporization, function of air temperature (J.kg-1)
    lambda_val = 2501000 - 2364 * input_temperature

    # Psychrometric constant (kPa K-1)
    psych = 0.0646 * exp(0.00097 * input_temperature)
    # Straight line approximation of the true slope; used in determining
    # the relationship slope
    mult = input_temperature + 237.3
    # Rate of change of saturation vapor pressure with temperature (kPa.K-1)
    slope = (2502.935945 * exp(17.269 * input_temperature / mult)) / (mult * mult)

    # Estimate frequently used atmospheric demand component
    ET_demand_coef = air_density_kg * cpair * input_vpd_kPa

    # Determine diffusion coefficient (m2.s-1), temperature-dependent (pressure dependence neglected).
    # Jones p51; appendix 2
    # Temperature adjusted from standard 20oC (293.15 K), NOTE that 1/293.15 = 0.003411223
    # 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    water_vapour_diffusion = 0.0000242 * ((input_temperature_K / 293.15) ** 1.75)

    # Calculate the dynamic viscosity of air (kg.m-2.s-1)
    dynamic_viscosity = ((input_temperature_K ** 1.5) / (input_temperature_K + 120)) * 1.4963e-6
    # Kinematic viscosity (m2.s-1)
    kinematic_viscosity = dynamic_viscosity / air_density_kg



# # --------------------------------------------------------------------------------------------------------------------------------------------------------------------


# import numpy as np

# def opt_max_scaling(max_val, min_val, optimum, kurtosis, current):
#     """
#     Estimates a 0-1 scaling based on a skewed Gaussian distribution with a
#     given optimum, maximum, and kurtosis. Minimum is assumed to be at infinity
#     (or near enough).
#     """
    
#     # Code with implicit assumption of min bound at infinity
#     # if current >= max_val:
#     #     return 0.0
#     # else:
#     #     dummy = np.exp(np.log((max_val - current) / (max_val - optimum)) * kurtosis * (max_val - optimum))
#     #     return dummy * np.exp(kurtosis * (current - optimum))
    
#     # Code with explicit min bound
#     opt_max_scaling = np.exp(kurtosis * np.log((max_val - current) / (max_val - optimum)) * (max_val - optimum)) \
#                     * np.exp(kurtosis * np.log((current - min_val) / (optimum - min_val)) * (optimum - min_val))
    
#     # Sanity check, allows for overlapping parameter ranges
#     if np.isnan(opt_max_scaling):
#         opt_max_scaling = 0.0
    
#     return opt_max_scaling

# import numpy as np

# def calculate_shortwave_balance(lai, sw_par_fraction, swrad, max_par_transmitted, max_nir_transmitted,
#                                 max_nir_reflected, max_par_reflected, soil_swrad_absorption):
    
#     """
#     Subroutine estimates the canopy and soil absorbed shortwave radiation (MJ/m2/day).
#     Radiation absorption is paritioned into NIR and PAR for canopy, and NIR + PAR for soil.

#     SPA uses a complex multi-layer radiative transfer scheme including
#     reflectance, transmittance any absorption. However, for a given
#     canopy vertical profiles, the LAI absorption relationship is readily
#     predicted via Michaelis-Menten or non-rectangular hyperbola as done here.
#     """

#     # Local parameters
#     clump = 1.0  # Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75); He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
#     decay = -0.5  # Decay coefficient for incident radiation

#     # Determine canopy absorption, reflectance and transmittance as function of LAI

#     # First, we consider how much radiation is likely to be incident on the
#     # canopy, or put another way what fraction passes straight through the canopy?
#     # Local variables
#     transmitted_fraction = np.exp(decay * lai * clump)
#     # Second, of the radiation which is incident on the canopy what fractions
#     # are transmitted through, reflected from or absorbed by the canopy
#     canopy_transmitted_fraction = np.exp(decay * lai * 0.5 * clump)

#     # Canopy transmitted of PAR & NIR radiation towards the soil
#     trans_par_fraction = canopy_transmitted_fraction * max_par_transmitted
#     trans_nir_fraction = canopy_transmitted_fraction * max_nir_transmitted
#     # Canopy reflected of near infrared and photosynthetically active radiation
#     reflected_nir_fraction = canopy_transmitted_fraction * max_nir_reflected
#     reflected_par_fraction = canopy_transmitted_fraction * max_par_reflected
#     # Canopy absorption of near infrared and photosynthetically active radiation
#     absorbed_nir_fraction = 1.0 - reflected_nir_fraction - trans_nir_fraction
#     absorbed_par_fraction = 1.0 - reflected_par_fraction - trans_par_fraction

#     # Estimate canopy absorption of incoming shortwave radiation
#     # Estimate multiple use par and nir components
#     par = sw_par_fraction * swrad
#     nir = (1.0 - sw_par_fraction) * swrad

#     # Estimate the radiation which directly bypasses the canopy...
#     trans_par_MJday = par * transmitted_fraction
#     trans_nir_MJday = nir * transmitted_fraction
#     # ...and update the canopy intercepted radiation
#     par -= trans_par_MJday
#     nir -= trans_nir_MJday

#     # Estimate incoming shortwave radiation absorbed, transmitted and reflected by the canopy (MJ.m-2.day-1)
#     canopy_par_MJday = par * absorbed_par_fraction
#     canopy_nir_MJday = nir * absorbed_nir_fraction
#     trans_par_MJday += par * trans_par_fraction
#     trans_nir_MJday += nir * trans_nir_fraction
#     refl_par_MJday = par * reflected_par_fraction
#     refl_nir_MJday = nir * reflected_nir_fraction

#     # Estimate soil absorption of shortwave passing through the canopy

#     # Then the radiation incident and ultimately absorbed by the soil surface itself (MJ.m-2.day-1)
#     soil_par_MJday = trans_par_MJday * soil_swrad_absorption
#     soil_nir_MJday = trans_nir_MJday * soil_swrad_absorption
#     # combine totals for use is soil evaporation
#     soil_swrad_MJday = soil_nir_MJday + soil_par_MJday


#     # Estimate canopy absorption of soil reflected shortwave radiation
#     # This additional reflection / absorption cycle is needed to ensure > 0.99
#     # of incoming radiation is explicitly accounted for in the energy balance.
#     # calculate multiple use variables
#     par = trans_par_MJday - soil_par_MJday
#     nir = trans_nir_MJday - soil_nir_MJday
#     # how much of the reflected radiation directly bypasses the canopy...
#     refl_par_MJday += par * transmitted_fraction
#     refl_nir_MJday += nir * transmitted_fraction
#     # ...and update the canopy on this basis
#     par *= (1.0 - transmitted_fraction)
#     nir *= (1.0 - transmitted_fraction)

#     # Update the canopy radiation absorption based on the reflected radiation (MJ.m-2.day-1)
#     canopy_par_MJday += par * absorbed_par_fraction
#     canopy_nir_MJday += nir * absorbed_nir_fraction
#     # Update the total radiation reflected back into the sky, i.e. that which is
#     # now transmitted through the canopy
#     refl_par_MJday += par * trans_par_fraction
#     refl_nir_MJday += nir * trans_nir_fraction

#     # Combine to estimate total shortwave canopy absorbed radiation
#     canopy_swrad_MJday = canopy_par_MJday + canopy_nir_MJday

#     # Check energy balance (uncomment if needed)
#     # balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
#     # if (np.abs((balance - swrad) / swrad) > 0.01):
#     #     print("SW residual frac =", (balance - swrad) / swrad, "SW residual =", balance, "SW in =", swrad)

#     return par, absorbed_par_fraction#, canopy_swrad_MJday, soil_swrad_MJday

# def arrhenious(a, b, t):
#     # The equation is simply a * exp(b * (t - 25.0) / (t + 273.15))
#     # However, precision in this routine matters as it affects many others.
#     # To maximize precision, the calculations have been split.

#     # Local variables
#     numerator = t - 25.0
#     denominator = t + 273.15
#     answer = a * np.exp(b * numerator / denominator)

#     return answer


# # constants: 
# gC_to_umol = 83333.33e0    # conversion of gC -> umolC; umol_to_gC**(-1d0)

# pn_max_temp = 6.842942e1   # Maximum daily max temperature for photosynthesis (oC)
# pn_min_temp = -1e6         # Minimum daily max temperature for photosynthesis (oC)
# pn_opt_temp = 3.155960e1   # Optimum daily max temperature for photosynthesis (oC)
# pn_kurtosis = 1.889026e-1  # Kurtosis of photosynthesis temperature response

# e0 = 3.661204e0  # Quantum yield gC/MJ/m2/day PAR
# sw_par_fraction = 0.5e0 # fraction of short-wave radiation which is PAR
# max_par_transmitted = 1.628077e-1 # Max fraction of canopy incident PAR transmitted to soil
# max_nir_transmitted = 2.793660e-01 # Max fraction of canopy incident NIR transmitted to soil
# max_par_reflected = 1.629133e-01  # Max fraction of canopy incident PAR reflected to sky
# max_nir_reflected = 4.284365e-01  # Max fraction of canopy incident NIR reflected to sky
# soil_swrad_absorption = 9.989852e-01 # Fraction of SW rad absorbed by soil

# mint = met(2,n)  # minimum temperature (oC) # in
# maxt = met(3,n)  # maximum temperature (oC) # in
# swrad = met(4,n) # incoming short wave radiation (MJ/m2/day)
# leafT = (maxt*0.75) + (mint*0.25)   # initial day time canopy temperature (oC)


# # Estimate incoming shortwave radiation absorbed, transmitted and reflected by the canopy (MJ.m-2.day-1)
# canopy_par_MJday = par * absorbed_par_fraction


# def acm_gpp_stage_1(lai, ceff, leafT):
#     '''
#     Estimate the light and temperature limited photosynthesis components.
#     See acm_gpp_stage_2() for estimation of CO2 supply limitation and combination of light, temperature and CO2 co-limitation

#     Metabolic limited photosynthesis

#     maximum rate of temperature and nitrogen (canopy efficiency) limited
#     photosynthesis (gC.m-2.day-1 -> umolC/m2/day)
#     '''
#     # ceff: Canopy efficiency (gC/m2leaf/day)

#     metabolic_limited_photosynthesis = gC_to_umol*lai*ceff*opt_max_scaling(pn_max_temp, pn_min_temp, pn_opt_temp, pn_kurtosis,leafT)
#     # Light limited photosynthesis
#     # calculate light limted rate of photosynthesis (gC.m-2.day-1)
#     light_limited_photosynthesis = e0 * canopy_par_MJday


# #   subroutine acm_gpp_stage_1