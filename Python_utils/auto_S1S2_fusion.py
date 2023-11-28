'''
Sentinel-1/-2 data fusion to generate high-resolution LAI for the UK EOCIS project
Author: Songyan Zhu
Contact: szhu4@ed.ac.uk
---------Log of Changes-------------
Created: 2023-10-20
Updated: 2023-10-25
    |-> Reordered the VPD processing chain (sel time then interp) to save some time
    |-> Added support of Xgboost and catboost to further speed up for the future when GPU is available
Updated: 2023-11-04
    |-> Added trash collection to save memory
    |-> Added exception catching
    !-> Surpress non-fatal warnings 
Updated: 2023-11-05
    |-> Added argparse for operational runs
Updated: 2023-11-06
    |-> Integrate daily fused LAI netCDF files into one and resample into weekly
Updated: 2023-11-06
    |-> Rename to auto_S1_S2_fusion.py
    |-> Rename config file to auto_S1_S2_fusion.yaml
'''

import yaml
import pickle
import warnings
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
try:
    from xgboost import XGBRegressor
except Exception as e:
    print("xgboost is not installed but it's not essential...")
try:
    from catboost import CatBoostRegressor
except Exception as e:
    print("catboost is not installed but it's not essential...")

warnings.simplefilter('ignore')


def load_yaml_config(p):
    with open(p, "r") as stream:
        try:
            yamlfile = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            assert(exc)
    return yamlfile

def load_pickle(p):
    with open(p, "rb") as f:
        ds = pickle.load(f)
    return ds

def dump_pickle(ds, p, large = False):
    with open(p, "wb") as f:
        if large:
            pickle.dump(ds, f, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            pickle.dump(ds, f)

def load_s1s2(p, dt, names, usenames = None, epsg_code = '4326', shp = None, buffer_val = 0):
    if usenames == None: usenames = names
    rnc = rxr.open_rasterio(p, band_as_variable = True)
    rnc = rnc.rio.write_crs(f"EPSG: {epsg_code}", inplace = False)
    if epsg_code != '4326':
        rnc = rnc.rio.reproject("EPSG:4326")
    d = dict(zip(list(rnc.keys()), names))
    d.update({'x': 'longitude', 'y': 'latitude'})
    rnc = rnc.rename(d)[usenames]
    rnc = rnc.expand_dims(time=[dt])
    if not (shp is None): 
        rnc = rnc.rio.clip(shp.geometry.buffer(buffer_val).values, shp.crs)#.drop_vars('spatial_ref')
    return rnc

def train(verbose = 0):
    cfg = load_yaml_config('auto_S1S2_fusion.yaml')
    s1_folder = Path(cfg['s1_folder'])
    s2_folder = Path(cfg['s2_folder'])
    vpd_folder = Path(cfg['vpd_folder'])
    shp_file = Path(cfg['shp_file'])
    epsg1 = cfg['epsg1']
    epsg2 = cfg['epsg2']

    names1 = ['VV', 'VH']
    names2 = ['LAI', 'LAIunc']
    df_path_s1 = pd.DataFrame([[pd.to_datetime(p.stem.split('S1_D')[1], format = '%Y%m%d'), p] for p in s1_folder.glob('*.tif')], columns = ['DATETIME', 'PATH_S1']).set_index('DATETIME').sort_index()
    df_path_s2 = pd.DataFrame([[pd.to_datetime(p.stem.split('_')[-1].split('T')[0], format = '%Y%m%d'), p] for p in s2_folder.glob('*.tif')], columns = ['DATETIME', 'PATH_S2']).set_index('DATETIME').sort_index()
    df_path_vpd = pd.DataFrame([[pd.to_datetime(p.stem.split('_')[-1], format = '%Y%m'), p] for p in vpd_folder.glob('vpd_daily_mean*.nc')], columns = ['DATETIME', 'PATH_VPD']).set_index('DATETIME').sort_index()

    train_dts = df_path_s1.index.drop_duplicates().intersection(df_path_s2.index.drop_duplicates())

    if shp_file:
        shp = gpd.read_file(shp_file)
        if shp.crs.to_epsg() != 4326: shp = shp.to_crs(4326)


    nc_train = []
    for dt in train_dts:
        try:
            # if dt != pd.to_datetime('2017-07-14'): continue
            s1_ps = df_path_s1.loc[dt, 'PATH_S1']
            s2_ps = df_path_s2.loc[dt, 'PATH_S2']
            if verbose: print(s1_ps, s2_ps)
            
            if type(s1_ps) == pd.core.series.Series:
                s1 = xr.merge([load_s1s2(p, dt, names1, epsg_code = epsg1, shp = shp) for p in s1_ps])
            else:
                s1 = load_s1s2(s1_ps, dt, names1, epsg_code = epsg1, shp = shp)
            
            if type(s2_ps) == pd.core.series.Series:
                s2 = xr.merge([load_s1s2(p, dt, names2, epsg_code = epsg2, shp = shp) for p in s2_ps])
            else:
                s2 = load_s1s2(s2_ps, dt, names2, epsg_code = epsg2, shp = shp)
            s2 = s2.where(
                (abs(s2['LAI']) < 1000) & (s2['LAI'] > 0)
            )
            s1 = s1.interp(longitude = s2.longitude, latitude = s2.latitude)
            
            nct = xr.merge([s1, s2])
                
            vpd_p = df_path_vpd.loc[pd.to_datetime(f'{dt.year}-{str(dt.month).zfill(2)}-01', format = '%Y-%m-%d'), 'PATH_VPD']
            vpd = xr.load_dataset(vpd_p, drop_variables = ['tstep', 'lat_dim', 'long_dim']).set_coords(['timestp', 'Latitude', 'Longitude'])
            vpd = vpd.rename({
                'timestp': 'time', 'tstep': 'time', 
                'lat_dim': 'latitude', 'Latitude': 'latitude',
                'long_dim': 'longitude', 'Longitude': 'longitude'
            })
            vpd = vpd.where(vpd.time == dt, drop=True)
            vpd = vpd.interp(longitude = nct.longitude, latitude = nct.latitude)['vpd_mean']

            nct = xr.merge([nct, vpd])
            
            nc_train.append(nct)
            del(s1); del(s2); del(vpd); del(nct)
        except Exception as e:
            print(e)
        
    nc_train = xr.merge(nc_train)
    df_train = nc_train[['LAI', 'VV', 'VH', 'vpd_mean']].drop_vars('spatial_ref').to_dataframe().dropna().reset_index().set_index('time')
    df_train['year'] = df_train.index.year
    df_train['month'] = df_train.index.month
    df_train['day'] = df_train.index.day

    X_train, X_test, y_train, y_test = train_test_split(df_train[['VV', 'VH', 'vpd_mean', 'year', 'month', 'day']], df_train['LAI'], test_size=0.2,random_state=0) # , 'latitude', 'longitude'
    regr = RandomForestRegressor(n_estimators=100) 
    regr.fit(X_train, y_train)

    savefile = Path(cfg['model_file'])
    dump_pickle(regr, savefile)

    if verbose:
        model_score = (regr.score(X_test, y_test))
        print(model_score)


def run(verbose = 0):
    cfg = load_yaml_config('auto_S1S2_fusion.yaml')
    s1_folder = Path(cfg['s1_folder'])
    vpd_folder = Path(cfg['vpd_folder'])
    shp_file = Path(cfg['shp_file'])
    epsg1 = cfg['epsg1']
    savefile = Path(cfg['model_file'])
    regr = load_pickle(savefile)

    names1 = ['VV', 'VH']
    df_path_s1 = pd.DataFrame([[pd.to_datetime(p.stem.split('S1_D')[1], format = '%Y%m%d'), p] for p in s1_folder.glob('*.tif')], columns = ['DATETIME', 'PATH_S1']).set_index('DATETIME').sort_index()
    df_path_vpd = pd.DataFrame([[pd.to_datetime(p.stem.split('_')[-1], format = '%Y%m'), p] for p in vpd_folder.glob('vpd_daily_mean*.nc')], columns = ['DATETIME', 'PATH_VPD']).set_index('DATETIME').sort_index()

    if shp_file:
        shp = gpd.read_file(shp_file)
        if shp.crs.to_epsg() != 4326: shp = shp.to_crs(4326)
            
    for dt in df_path_s1.index:
        try:
            if verbose: print(dt)
            s1_ps = df_path_s1.loc[dt, 'PATH_S1']
            if type(s1_ps) == pd.core.series.Series:
                s1 = xr.merge([load_s1s2(p, dt, names1, epsg_code = epsg1, shp = shp) for p in s1_ps])
            else:
                s1 = load_s1s2(s1_ps, dt, names1, epsg_code = epsg1, shp = shp)
            nct = s1

            vpd_p = df_path_vpd.loc[pd.to_datetime(f'{dt.year}-{str(dt.month).zfill(2)}-01', format = '%Y-%m-%d'), 'PATH_VPD']
            vpd = xr.load_dataset(vpd_p, drop_variables = ['tstep', 'lat_dim', 'long_dim']).set_coords(['timestp', 'Latitude', 'Longitude'])
            vpd = vpd.rename({
                'timestp': 'time', 'tstep': 'time', 
                'lat_dim': 'latitude', 'Latitude': 'latitude',
                'long_dim': 'longitude', 'Longitude': 'longitude'
            })
            vpd = vpd.where(vpd.time == dt, drop=True)
            vpd = vpd.interp(longitude = nct.longitude, latitude = nct.latitude)['vpd_mean']

            nct = xr.merge([nct, vpd])
            df_app = nct[['VV', 'VH', 'vpd_mean']].drop_vars('spatial_ref').to_dataframe().dropna().reset_index().set_index('time')
            df_app['year'] = df_app.index.year
            df_app['month'] = df_app.index.month
            df_app['day'] = df_app.index.day

            lai = regr.predict(df_app[['VV', 'VH', 'vpd_mean', 'year', 'month', 'day']]) # , 'latitude', 'longitude'
            df_app['LAI'] = lai
            nco = df_app[['latitude', 'longitude', 'LAI']].reset_index().set_index(['time', 'latitude', 'longitude']).to_xarray()

            savefile = Path(cfg['savefolder']).joinpath('FLAI_' + dt.strftime('%Y%m%d') + '.nc')

            if not savefile.exists(): nco.to_netcdf(savefile)

            del(s1); del(nct); del(vpd); del(df_app); del(lai); del(nco)
        except Exception as e:
            print(e)

def integrate_resample_weekly(verbose = 0):
    cfg = load_yaml_config('auto_S1S2_fusion.yaml')
    savefolder = cfg['savefolder']
    savefolder = Path(savefolder)

    # nc = xr.merge([xr.open_dataset(p) for p in savefolder.glob('*.nc')])
    df_path = []
    for p in list(savefolder.glob('*.nc')):
        dt = pd.to_datetime(p.stem.split('_')[1], format = '%Y%m%d')
        df_path.append([dt, p])
    df_path = pd.DataFrame(df_path, columns = ['DATETIME', 'PATH']).set_index('DATETIME').sort_index()

    nc0 = xr.open_dataset(df_path['PATH'][0])
    nc = []
    for dt in df_path.index:
        p = df_path.loc[dt, 'PATH']
        nct = xr.open_dataset(p)
        nct = nct.interp(latitude = nc0.latitude, longitude = nc0.longitude)
        nc.append(nct)
        if verbose: print(dt)

    nc = xr.merge(nc)

    start_day = pd.to_datetime(nc.time.values[0])
    start_day = pd.to_datetime(f'{start_day.year}-{start_day.month}-01', format = '%Y-%m-%d')

    nc = nc.resample(time = '7D', origin = start_day).mean()
    nc.to_netcdf('FLAI_weekly.nc')


if __name__ == '__main__':
    # Example: auto_S1S2_fusion.py -m t -v 0
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode')
    parser.add_argument('-v', '--verbose')
    args = parser.parse_args()
    run_mode = args.mode
    verbose = args.verbose
    if run_mode.lower() in ['t', 'train', 'training']:
        train(verbose = verbose)
    elif run_mode.lower() in ['r', 'run', 'running']:
        run(verbose = verbose)
    elif run_mode.lower() in ['w']:
        integrate_resample_weekly(verbose = verbose)