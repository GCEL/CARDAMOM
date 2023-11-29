'''
Auto-generating Sentinel-2 LAI
Author: Songyan Zhu
Contact: szhu4@ed.ac.uk
---------Log of Changes-------------
Make sure you have SNAP installed, 
snappy installed and configured, 
and Python 3.6 virtual environment activated!

Created: 2023-11-16
Updated: 2023-11-16
    |-> Make it operational
Updated: 2023-11-17
    |-> Add error catching for unzip
Updated: 2023-11-29
    |-> Remove failed files for both unzipping and biophysical
    |-> Save nc to root_proj.joinpath(data_collection)
'''

import os
import sys
import yaml
import zipfile
import argparse
import shutil
from pathlib import Path

def load_yaml_config(p):
    with open(p, "r") as stream:
        try:
            yamlfile = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            assert(exc)
    return yamlfile

def init_snappy(snappy_dir):
    sys.path.append(snappy_dir)
    # from snappy import ProductIO
    # from snappy import HashMap
    # from snappy import jpy
    # from snappy import GPF
    # # Example: 
    # file_path = '/exports/csce/datastore/geos/groups/gcel/EOCIS/SNAP_LAI/SnapPy.py/snappy/testdata/MER_FRS_L1B_SUBSET.dim'
    # p = ProductIO.readProduct(file_path)
    # list(p.getBandNames())

def init_env(p_config):
    cfg = load_yaml_config(p_config)
    init_snappy(cfg['snappy_dir'])
    root_proj = Path(cfg['workspace'])
    data_collection = cfg['data_collection']
    return root_proj, data_collection

def unzip_sentinel(root_proj, data_collection):
    paths = root_proj.glob(f'{data_collection}*.zip')
    for src in paths:
        try:
            dest = src.parent.joinpath(data_collection)
            # if dest.joinpath(src.stem + '.SAFE').exists():
            #     print('yes')
            print(dest, src.stem)
            with zipfile.ZipFile(src.as_posix(),"r") as zip_ref:
                zip_ref.extractall(dest.as_posix())
            os.remove(src)
        except Exception as e:
            print(e)
            with open('auto_SNAP_LAI_logging.txt', 'a') as f:
                f.write(src.as_posix() + '\n')
            os.remove(src)
            
def generate_biophysics(root_proj, data_collection):
    from snappy import ProductIO
    from snappy import HashMap
    from snappy import jpy
    from snappy import GPF
    paths = list(root_proj.joinpath(data_collection).glob('*.SAFE'))
    n_files = len(paths)
    print(f'Beginning to retrieve {n_files} {data_collection} files...')
    for cnt, p in enumerate(paths):
        try:
            print(p.stem)
            p = list(p.glob('MTD*.xml'))[0]
            savefile = root_proj.joinpath(data_collection).joinpath(p.parent.stem + '.nc')
            if savefile.exists(): 
                print(f'{cnt + 1} done, {n_files - cnt} remaining...')
                print('-' * 100)
                continue

            # -----------------------------------------------------------------------
            prod = ProductIO.readProduct(p.as_posix())

            GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
            HashMap = jpy.get_type('java.util.HashMap')

            parameters = HashMap()
            parameters.put('targetResolution', 20)
            prod_resampled = GPF.createProduct('Resample', parameters, prod)
            # -----------------------------------------------------------------------
            # https://forum.step.esa.int/t/biophysical-parameter-snappy/11507/3
            parameters = HashMap()
            parameters.put('computeCab', True)
            parameters.put('computeCw', True)
            parameters.put('computeFapar', True)
            parameters.put('computeFcover', True)
            parameters.put('computeLAI', True)
            prod_bio = GPF.createProduct('BiophysicalOp',parameters,prod_resampled)
            # -----------------------------------------------------------------------
            ProductIO.writeProduct(prod_bio, savefile.as_posix(), 'NetCDF4-CF') # "GeoTIFF-BigTIFF"
            print(f'{cnt + 1} done, {n_files - cnt} remaining...')
            print('-' * 100)
            shutil.rmtree(p.parent)
        except Exception as e:
            print(e)
            shutil.rmtree(p.parent)

if __name__ == '__main__':
    # Example: python auto_SNAP_LAI.py -m unzip
    # Example: python auto_SNAP_LAI.py -m lai
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", nargs = "?", default = 'auto_SNAP_LAI.yaml', type = str)
    parser.add_argument('-m', '--mode')
    args = parser.parse_args()
    p_config = args.p
    run_mode = args.mode

    # conda activate snap # Python 3.6 virtual environment 
    root_proj, data_collection = init_env(p_config)
    if run_mode.lower() in ['unzip', 'u']:
        unzip_sentinel(root_proj, data_collection)
    elif run_mode.lower() in ['lai']:
        generate_biophysics(root_proj, data_collection)
    else:
        raise Exception('Wrong mode, unzip or generate LAI!')