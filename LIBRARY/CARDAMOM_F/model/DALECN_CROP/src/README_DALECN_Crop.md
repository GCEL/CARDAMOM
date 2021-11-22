# DALEC-CropN
DALEC-Crop model developments for simulating the impacts of N dilution on the ACM equation used to determine the leaf N and temperature limited rate of photosynthesis.

Code used for inferring crop leaf nitrogen per area (LNA) based on the calibration of a semi-empirical critical N dilution function derived from the ATEC project experimental field trials (see plot design). Code adapted from DALEC_GSI_DFOL_FR_CROP.f90 - key changes are the use of N dilution slope and intercept parameters applied in lines 376 to 389. From a calibration of parameters 1-34 based on a prior calibration step, the CARDAMOM MDF framework optimises the slope and intercept of the N dilusion (parameters 35 and 36) in order to find an N dilution fit consistent with an LAI observational time-series.

![TRIALS_Github](https://user-images.githubusercontent.com/43847496/97121237-9dd99800-1714-11eb-8e95-5672db870db2.png)

![N_fit_Github](https://user-images.githubusercontent.com/43847496/97122637-3d038d00-171f-11eb-96b7-2966736ecfe4.png)

A DALEC-CropN parameter file (see: crosstreatment_calibration.csv) that was generated based on a calibration of key parameters determining C availability and phenology (see table below) across a gradient of N fertiliser treatment rates for two winter wheat seasons under the constraints of LAI observations.

![calibrated_parameters](https://user-images.githubusercontent.com/43847496/97156962-5a604780-176f-11eb-8ada-18ce864c9ce0.png)
