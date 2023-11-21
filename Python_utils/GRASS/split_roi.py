# UK small

minlon = -8.5
maxlon = 2
minlat = 49.8
maxlat = 59.6

nrow = 3
ncol = 3

rois = geobox.split_roi(minlon, maxlon, minlat, maxlat, nrow, ncol)
rois = [roi[0] + roi[2] for roi in rois]
print(len(rois))
print(rois)