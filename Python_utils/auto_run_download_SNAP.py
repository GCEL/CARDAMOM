'''
Automatically run the whole processing chain
Author: Songyan Zhu
Contact: szhu4@ed.ac.uk
---------Log of Changes-------------
Created: 2023-11-17
Updated: 2023-11-28
    |-> Add splitting roi function
Updated: 2023-11-28
    |-> Check if the downloaded log file exists
'''

import os
import subprocess
import yaml
import numpy as np
import pandas as pd
from split_roi import split_roi
from auto_Sentinel_download import load_yaml_config
from pathlib import Path

if __name__ == '__main__':
    # UK small

    minlon = -8.5
    maxlon = 2
    minlat = 49.8
    maxlat = 59.6

    nrow = np.ceil(maxlat - minlat) # 1 deg
    ncol = np.ceil(maxlon - minlon) # 1 deg

    nrow = int(nrow)
    ncol = int(ncol)
    rois = split_roi(minlon, maxlon, minlat, maxlat, nrow, ncol)
    rois = [f"POLYGON(({roi[0][0]} {roi[0][1]},{roi[1][0]} {roi[1][1]},{roi[2][0]} {roi[2][1]},{roi[3][0]} {roi[3][1]},{roi[0][0]} {roi[0][1]}))'" for roi in rois]


    run_dates = pd.date_range('2017-01-01', '2023-01-01', freq = '1D')
    for rn in range(len(run_dates) - 1):
        start_date = run_dates[rn].strftime('%Y-%m-%d')
        end_date = run_dates[rn + 1].strftime('%Y-%m-%d')
        print(start_date, end_date)
        print('=' * 100)
        download_config = load_yaml_config('auto_Sentinel_download.yaml')
        download_config.update({
            'start_date': start_date,
            'end_date': end_date
        })
        with open('auto_Sentinel_download.yaml', 'w') as file:
            _ = yaml.dump(download_config, file, sort_keys=False)

        download_config = load_yaml_config('auto_Sentinel_download.yaml')
        for roi in rois:
            download_config.update({
                'roi_str': roi
            })

            with open('auto_Sentinel_download.yaml', 'w') as file:
                _ = yaml.dump(download_config, file, sort_keys=False)
        
            # subprocess.run('conda activate base', shell = True)
            print('Downloading...')
            subprocess.run('python auto_Sentinel_download.py', shell = True)

            # recorded downloaded files...
            if Path('auto_Sentinel_downloaded_files.txt').exists():
                with open('auto_Sentinel_downloaded_files.txt', 'r') as f:
                    lines = f.readlines()
                    lines = [l.strip() for l in lines]
                if 'downloaded_files' in download_config.keys():
                    lines = download_config['downloaded_files'] + lines
                    lines = np.array(lines)
                    lines = [str(l) for l in np.unique(lines)]
                    download_config['downloaded_files'] = lines
                else:
                    download_config['downloaded_files'] = lines

                with open('auto_Sentinel_download.yaml', 'w') as file:
                    _ = yaml.dump(download_config, file, sort_keys=False)
                os.remove('auto_Sentinel_downloaded_files.txt')

            # subprocess.run('conda activate /exports/csce/datastore/geos/groups/gcel/EOCIS/SNAP_LAI/conda_env/snap', shell = True)
            # subprocess.run('python auto_SNAP_LAI.py -m unzip', shell = True)
            # subprocess.run('python auto_SNAP_LAI.py -m lai', shell = True)
            # subprocess.run('conda deactivate', shell = True)
            print('Unzipping...')
            try:
                subprocess.run("conda run -n snap python auto_SNAP_LAI.py -m unzip", shell = True)
            except Exception as e:
                print(e)
            print('Generating LAI...')
            try:
                subprocess.run("conda run -n snap python auto_SNAP_LAI.py -m lai", shell = True)
            except Exception as e:
                print(e)
