# softgating-recon

Description:
Softgating recon matlab function allows to perform UTE MR image reconstruction using one of the five available soft-gating functions. There is a matlab script 'run_soft_recon.m' as well a python script 'Softgating_multiRecon.py'

Requirements:
Matlab 2020 or later
GE recon functions (proprietary; not included in this repository)

Usage:
[f] = run_soft_recon(Patient,PfileName,functype,MAGorPHS,const_width,sigma,cal)

Inputs:
Patient: Name of the patient (eg. 'P120')
P-file: Name of the P-file (eg. 'P00000.7')
functype: Name of the soft-gating function you need to reconstruct with (eg. 'Exponential')
MAGorPHS: Name of the the soft-gating weight assignments file (eg. 'SoftGating_8Bins_Phs.mat')
const_width: Threshold limit (eg.1)
sigma: Slope (eg. 0.5)
cal: Name of the calibration file ('P00000.7')

Author:
Deep B. Gandhi
