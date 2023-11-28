'''
Split regions into smaller rois for retrieving data from Copernicus Data Space Ecosystem 
Author: Songyan Zhu
Contact: szhu4@ed.ac.uk
---------Log of Changes-------------
Created: 2023-11-21
Updated: 2023-11-27
    |-> Make it standalone from the scigeo (also by szhu4) package
Updated: 2023-11-28
    |-> Generating roi in the format of Copernicus Data Space Ecosystem, polygons by 1deg x 1deg
To do:
    |-> Check the area of a tile to decide how many rois are needed.
'''

import numpy as np

# UK small

minlon = -8.5
maxlon = 2
minlat = 49.8
maxlat = 59.6

nrow = np.ceil(maxlat - minlat) # 1 deg
ncol = np.ceil(maxlon - minlon) # 1 deg

nrow = int(nrow)
ncol = int(ncol)

def split_roi(minlon, maxlon, minlat, maxlat, nrow, ncol):
    lons = np.linspace(minlon, maxlon, ncol)
    lats = np.linspace(minlat, maxlat, nrow)
    lons, lats = np.meshgrid(lons, lats)

    sub_coords = []
    for i in range(nrow - 1):
        for j in range(ncol - 1):
            rect_lons = [
                lons[i, j], lons[i, j + 1],
                lons[i + 1, j + 1], lons[i + 1, j]
            ]
            rect_lats = [
                lats[i, j], lats[i, j + 1],
                lats[i + 1, j + 1], lats[i + 1, j]
            ]
            pairs = [list(pair) for pair in zip(rect_lons, rect_lats)]
            sub_coords.append(pairs)
    return sub_coords

rois = split_roi(minlon, maxlon, minlat, maxlat, nrow, ncol)
# rois = [roi[0] + roi[2] for roi in rois]
rois = [f"POLYGON(({roi[0][0]} {roi[0][1]},{roi[1][0]} {roi[1][1]},{roi[2][0]} {roi[2][1]},{roi[3][0]} {roi[3][1]},{roi[0][0]} {roi[0][1]}))'" for roi in rois]
print(len(rois))
print(rois)