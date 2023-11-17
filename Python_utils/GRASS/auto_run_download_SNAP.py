import subprocess

if __name__ == '__main__':
    # subprocess.run('conda activate base', shell = True)
    print('Downloading...')
    subprocess.run('python auto_Sentinel_download.py', shell = True)
    # subprocess.run('conda activate /exports/csce/datastore/geos/groups/gcel/EOCIS/SNAP_LAI/conda_env/snap', shell = True)
    # subprocess.run('python auto_SNAP_LAI.py -m unzip', shell = True)
    # subprocess.run('python auto_SNAP_LAI.py -m lai', shell = True)
    # subprocess.run('conda deactivate', shell = True)
    print('Unzipping...')
    subprocess.run("conda run -n snap python auto_SNAP_LAI.py -m unzip", shell = True)
    print('Generating LAI...')
    subprocess.run("conda run -n snap python auto_SNAP_LAI.py -m lai", shell = True)