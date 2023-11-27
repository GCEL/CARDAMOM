import numpy as np

# UK small


minlon = -8.5
maxlon = 2
minlat = 49.8
maxlat = 59.6

nrow = 3
ncol = 3

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
rois = [roi[0] + roi[2] for roi in rois]
print(len(rois))
print(rois)