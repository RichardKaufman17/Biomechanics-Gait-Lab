# Biomechanics-Gait-Lab
In this lab, motion tracking in combination with a force plate was utilized to model a subject's gait. The data and code here were utilized to plot the subject's joint angles, joint angular velocities, and joint moments over the course of his stance phase. these metrics were compared for three different trials, where the subject walked, walked fast, and jogged. A full report for the lab, with the resuls, can be found in Biomechanics Gait Analysis Lab.pdf in this repository.
## Datasets
This repository contains three datasets in the form of .mat files:
- **normal03.mat:** Data acquired when the subject walked
- **fast03.mat:** Data acquired when the subject fast walked
- **jogging09.mat:** Data acquired when the subject jogged
<br />
Each file is the same format, and contains tables with the force plate and motion capture data

## Running the Code
This repository contains one MATLAB scrip to run the analysis: biomehcnaicsGaitLab. To run this script, make sure the three .mat files made available are in the same folder. 
