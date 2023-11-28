import numpy as np
import xarray as xr
from pathlib import Path


def check_tile_coverage(data_folder):
    paths = data_folder.glob('*.nc')
    for p in paths:
        nct = xr.load_dataset(p)
        lons = nct.lon.values
        lats = nct.lat.values
        # print(f'Lon: {lons[0], lons[-1], np.abs(lons[-1] - lons[0])}')
        # print(f'Lat: {lats[0], lats[-1], np.abs(lats[-1] - lats[0])}')
        print(p.stem)
        print('LON:', lons.min(), lons.max(), np.abs(lons.min() - lons.max()))
        print('LAT:', lats.min(), lats.max(), np.abs(lats.min() - lats.max()))
        nct.close()
        del(nct)

if __name__ == '__main__':
    data_folder = Path('/exports/csce/datastore/geos/groups/gcel/EOCIS/S2_LAI_workflow/S2_LAI_SNAP')
    check_tile_coverage(data_folder)