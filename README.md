# Biomechanics-Gait-Lab
In this lab, motion tracking in combination with a force plate was utilized to model a subject's gait. The data and code here were utilized to plot the subject's joint angles, joint angular velocities, and joint moments over the course of his stance phase. these metrics were compared for three different trials, where the subject walked, walked fast, and jogged. A full report for the lab, with the resuls, can be found in Biomechanics Gait Analysis Lab.pdf in this repository.
## Datasets
This repository contains three datasets in the form of .mat files:
- **normal03.mat:** Data acquired when the subject walked
- **fast03.mat:** Data acquired when the subject fast walked
- **jogging09.mat:** Data acquired when the subject jogged

### Tables 
Each dataset has the same format, and contains tables with the force plate and motion capture data. Below I describe the tables in each .mat file that are needed for the assignment:
<br />

<br />**digitaldata:** Contains the raw moment and force data acquired from the force plate. Note the subject walked in the X-direction with the Z-axis pointing up. This table has 6 columns:
1. Force in the X direction
2. Force in the Y direction
3. Force in the Z direction
4. Moment in the X direction
5. Moment in the Y direction
6. Moment in the Z direction

<br />**markerposition:** Contains the raw motion capture data, with each column representing the x, y or z position of a marker. The table has 39 columns:
- 1.-3. X,Y and Z positions of the right ankle
- 4.-6. X,Y and Z position of the tip of the subject's right foot
- 7.-9. X,Y and Z position of the right heel
- 10.-12. X,Y and Z position of the right knee
- 13.-15. X,Y and Z position of the right trochanter
- 16.-18. X,Y and Z positions of the left ankle
- 19.-21. X,Y and Z position of the tip of the subject's left foot
- 22.-24. X,Y and Z position of the left heel
- 25.-27. X,Y and Z position of the left knee
- 28.-30. X,Y and Z position of the left trochanter
- 31.-33. X,Y and Z position of the right shoulder
- 34.-36. X,Y and Z position of the left shoulder
- 37.-39. X,Y and Z position of the left Back

<br />**COP:** Contains the location of the center of pressure of the subject when he moved over the force plate. It is calculated using the force and moment data from the force plate. This table has two colums:
1. X-position of the COP
2. Y-Postion of the COP


## Running the Code
This repository contains one MATLAB scrip to run the analysis: biomehcnaicsGaitLab. To run this script, make sure the three .mat files made available are in the same folder. 
