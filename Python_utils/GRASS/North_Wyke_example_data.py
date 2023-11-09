'''
Authorship: Songyan Zhu (szhu4@ed.ac.uk)
Converting Python standalone CARDAMOM npy input data into csv for the common CARADAMOM framework 
'''


import numpy as np
import pandas as pd
from copy import deepcopy
from pathlib import Path

root = Path.cwd() 
workingdir = root.joinpath('workspace/project_data/CARDAMOM/MDF_DALEC_GRASS').as_posix()
sitename = 'greatfield'

meti = np.array(np.load('%s/%s_M.npy' %(workingdir,sitename)),order="F")
obs_laii = np.load("%s/%s_O.npy" %(workingdir,sitename))
met       = deepcopy(meti)
obs_lai   = deepcopy(obs_laii)

nodays   = met.shape[1]

colnames = [
    'run_day', 'minT', 'maxT', 'srad', 'co2', 'DOY', 'lagged_precip', 'Lai_loss', 'burnt_area_fraction',
    '21d_minT', '21d_photoperiod', '21d_vpd', 'Forest_mgmt_after_clearing', 'Mean_T'
]

df_met = pd.DataFrame(met.transpose(), columns = colnames, index = pd.date_range('2017-01-01', periods = nodays, freq='7D'))
df_met['doy'] = df_met.index.dayofyear
df_met['mint_C'] = df_met['minT']
df_met['maxt_C'] = df_met['maxT']

df_obs_lai = pd.DataFrame(obs_lai, columns = ['LAI_m2m2'], index = pd.date_range('2017-01-01', periods = nodays, freq='7D'))
df_obs_lai['LAI_unc_m2m2'] = df_obs_lai['LAI_m2m2'] * 0.0921436 + 0.2124204
df_obs_lai['doy'] = df_obs_lai.index.dayofyear
df_obs_lai = df_obs_lai[['doy', 'LAI_m2m2', 'LAI_unc_m2m2']]

df_obs_lai.to_csv('North_Wyke_initial_obs.csv')
df_met.to_csv('North_Wyke_timeseries_met.csv')