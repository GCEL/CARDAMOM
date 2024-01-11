# from C code in DALEC_CDEA_GLEAM

import math

def ACM(pars, consts):
    gc = (abs(pars[8]) ** consts[9]) / (consts[5] * pars[9] + 0.5 * (pars[1] - pars[2]))
    pp = pars[0] * pars[3] / gc * consts[0] * math.exp(consts[7] * pars[1])
    qq = consts[2] - consts[3]
    ci = 0.5 * (pars[4] + qq - pp + math.sqrt((pars[4] + qq - pp) ** 2 - 4 * (pars[4] * qq - pp * consts[2])))
    e0 = consts[6] * (pars[0] ** 2) / ((pars[0] ** 2) + consts[8])
    dec = -23.4 * math.cos((360. * (pars[5] + 10.) / 365.) * pars[10] / 180.) * pars[10] / 180.
    mult = math.tan(pars[6] * pars[10] / 180) * math.tan(dec)

    if mult >= 1:
        dayl = 24.
    elif mult <= -1:
        dayl = 0.
    else:
        dayl = 24. * math.acos(-mult) / pars[10]

    cps = e0 * pars[7] * gc * (pars[4] - ci) / (e0 * pars[7] + gc * (pars[4] - ci))
    GPP = cps * (consts[1] * dayl + consts[4])
    return GPP

def offset(L, w):
    mxc = [0.000023599784710, 0.000332730053021, 0.000901865258885, -0.005437736864888,
           -0.020836027517787, 0.126972018064287, -0.188459767342504]
    lf = math.log(L - 1)
    os = mxc[0] * (lf ** 6) + mxc[1] * (lf ** 5) + mxc[2] * (lf ** 4) + mxc[3] * (lf ** 3) + \
         mxc[4] * (lf ** 2) + mxc[5] * lf + mxc[6]

    os = os * w
    return os

def CARBON_MODEL(MET, pars, deltat, nr, lat, LAI, NEE, FLUXES, POOLS):
    pi = 3.1415927
    timestep = 0.0  # Placeholder for the dynamic timestep

    gpppars = [0.0] * 11
    gpppars[3] = 1
    gpppars[6] = lat
    gpppars[8] = -2.0
    gpppars[9] = 1.0
    gpppars[10] = pi

    constants = [pars[10], 0.0156935, 4.22273, 208.868, 0.0453194, 0.37836, 7.19298, 0.011136, 2.1001, 0.789798]

    nopools = 6
    nomet = 10
    nofluxes = 28

    wf = pars[15] * math.sqrt(2) / 2
    wl = pars[13] * math.sqrt(2) / 2
    ff = (math.log(pars[4]) - math.log(pars[4] - 1)) / 2
    fl = (math.log(1.001) - math.log(0.001)) / 2
    osf = offset(pars[4], wf)
    osl = offset(1.001, wl)
    sf = 365.25 / pi

    CF = [0.1, 0.9, 0.1, 0.1, 0.5, 0.01]
    rfac = 0.5

    for n in range(nr):
        p = nopools * n
        nxp = nopools * (n + 1)
        m = nomet * n
        f = nofluxes * n

        LAI[n] = POOLS[p + 1] / pars[16]

        gpppars[0] = LAI[n]
        gpppars[1] = MET[m + 2]
        gpppars[2] = MET[m + 1]
        gpppars[4] = MET[m + 4]
        gpppars[5] = MET[m + 5]
        gpppars[7] = MET[m + 3]

        FLUXES[f + 0] = ACM(gpppars, constants) * MET[m + 9]
        FLUXES[f + 1] = math.exp(pars[9] * 0.5 * (MET[m + 2] + MET[m + 1]))
        FLUXES[f + 2] = pars[1] * FLUXES[f + 0]
        FLUXES[f + 3] = (FLUXES[f + 0] - FLUXES[f + 2]) * pars[2]
        FLUXES[f + 4] = (FLUXES[f + 0] - FLUXES[f + 2] - FLUXES[f + 3]) * pars[13 - 1]
        FLUXES[f + 5] = (FLUXES[f + 0] - FLUXES[f + 2] - FLUXES[f + 3] - FLUXES[f + 4]) * pars[4 - 1]
        FLUXES[f + 6] = FLUXES[f + 0] - FLUXES[f + 2] - FLUXES[f + 3] - FLUXES[f + 5] - FLUXES[f + 4]
        FLUXES[f + 8] = (2 / math.sqrt(pi)) * (ff / wf) * math.exp(-pow(math.sin((MET[m + 0] - pars[14] + osf) / sf) * sf / wf, 2))
        FLUXES[f + 15] = (2 / math.sqrt(pi)) * (fl / wl) * math.exp(-pow(math.sin((MET[m + 0] - pars[11] + osl) / sf) * sf / wl, 2))

        timestep = MET[m + 8]

        FLUXES[f + 7] = POOLS[p + 0] * (1 - pow(1 - FLUXES[f + 15], timestep)) / timestep
        FLUXES[f + 9] = POOLS[p + 1] * (1 - pow(1 - FLUXES[f + 8], timestep)) / timestep
        FLUXES[f + 10] = POOLS[p + 3] * (1 - pow(1 - pars[6 - 1], timestep)) / timestep
        FLUXES[f + 11] = POOLS[p + 2] * (1 - pow(1 - pars[7 - 1], timestep)) / timestep
        FLUXES[f + 12] = POOLS[p + 4] * (1 - pow(1 - FLUXES[f + 1] * pars[8 - 1], timestep)) / timestep
        FLUXES[f + 13] = POOLS[p + 5] * (1 - pow(1 - FLUXES[f + 1] * pars[9 - 1], timestep)) / timestep
        FLUXES[f + 14] = POOLS[p + 4] * (1 - pow(1 - pars[1 - 1] * FLUXES[f + 1], timestep)) / timestep

        POOLS[nxp + 0] = POOLS[p + 0] + (FLUXES[f + 4] - FLUXES[f + 7]) * timestep
        POOLS[nxp + 1] = POOLS[p + 1] + (FLUXES[f + 3] - FLUXES[f + 9] + FLUXES[f + 7]) * timestep
        POOLS[nxp + 3] = POOLS[p + 3] + (FLUXES[f + 6] - FLUXES[f + 10]) * timestep
        POOLS[nxp + 2] = POOLS[p + 2] + (FLUXES[f + 5] - FLUXES[f + 11]) * timestep
        POOLS[nxp + 4] = POOLS[p + 4] + (FLUXES[f + 9] + FLUXES[f + 11] - FLUXES[f + 12] - FLUXES[f + 14]) * timestep
        POOLS[nxp + 5] = POOLS[p + 5] + (FLUXES[f + 14] - FLUXES[f + 13] + FLUXES[f + 10]) * timestep

        NEE[n] = -FLUXES[f + 0] + FLUXES[f + 2] + FLUXES[f + 12] + FLUXES[f + 13]

        if MET[m + 6] > 0:
            POOLS[nxp + 0] *= (1 - MET[m + 6])
            POOLS[nxp + 1] *= (1 - MET[m + 6])
            POOLS[nxp + 3] *= (1 - MET[m + 6])

        FLUXES[f + 16] = 0

        if MET[m + 7] > 0:
            for nn in range(6):
                FLUXES[f + 17 + nn] = POOLS[nxp + nn] * MET[m + 7] * CF[nn] / timestep
            for nn in range(5):
                FLUXES[f + 23 + nn] = POOLS[nxp + nn] * MET[m + 7] * (1 - CF[nn]) * (1 - rfac) / timestep

            for nn in range(4):
                POOLS[nxp + nn] -= (FLUXES[f + 17 + nn] + FLUXES[f + 23 + nn]) * timestep

            POOLS[nxp + 4] += (FLUXES[f + 23] + FLUXES[f + 23 + 1] + FLUXES[f + 23 + 2] - FLUXES[f + 17 + 4] - FLUXES[f + 23 + 4]) * timestep
            POOLS[nxp + 5] += (FLUXES[f + 23 + 3] + FLUXES[f + 23 + 4] - FLUXES[f + 17 + 5]) * timestep

            for nn in range(6):
                FLUXES[f + 16] += FLUXES[f + 17 + nn]

        NEE[n] += FLUXES[f + 16]

# Example usage:
MET = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
pars = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2]
deltat = 0.1
nr = 10
lat = 45.0
LAI = [0.0] * nr
NEE = [0.0] * nr
FLUXES = [0.0] * (28 * nr)
POOLS = [0.0] * (6 * (nr + 1))

CARBON_MODEL(MET, pars, deltat, nr, lat, LAI, NEE, FLUXES, POOLS)

print(LAI)
print(NEE)
print(FLUXES)
print(POOLS)



# import numpy as np

# # gC_to_umol = 83333.33d0,     & ! conversion of gC -> umolC; umol_to_gC**(-1d0)


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


# def arrhenious(a, b, t):
#     # The equation is simply a * exp(b * (t - 25.0) / (t + 273.15))
#     # However, precision in this routine matters as it affects many others.
#     # To maximize precision, the calculations have been split.

#     # Local variables
#     numerator = t - 25.0
#     denominator = t + 273.15
#     answer = a * np.exp(b * numerator / denominator)

#     return answer



# def acm_gpp_stage_1():
#     # Estimate the light and temperature limited photosynthesis components.
#     # See acm_gpp_stage_2() for estimation of CO2 supply limitation and
#     # combination of light, temperature and CO2 co-limitation

#     # Metabolic limited photosynthesis

#     # maximum rate of temperature and nitrogen (canopy efficiency) limited
#     # photosynthesis (gC.m-2.day-1 -> umolC/m2/day)
#     metabolic_limited_photosynthesis = gC_to_umol*lai*ceff*opt_max_scaling(pn_max_temp,pn_min_temp,pn_opt_temp,pn_kurtosis,leafT)


# #   subroutine acm_gpp_stage_1