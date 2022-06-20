Before running any scripts add the folder this is contained within to your MATLAB path...

Cover figure:
"PATH\\SpeedForSubmissionAdditionalData\CodesToDrawFigures\DrawPyFracMeshesCoverFigure.m"
Fig.1 
"PATH\\SpeedForSubmissionAdditionalData\CodesToDrawFigures\DrawPyFracMeshes.m"
Fig.2 
"PATH\\SpeedForSubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\VelocityGravityFracturingResults2c.m" (flag on line 12: TenVol=0 and TenVol=1)
Fig.3 and Fig.4
"PATH\\SpeedForSubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\PyFracLargeScaleData\HeightHeightAndVelocityVelocityPlot.m"
Fig.5
"PATH\\SpeedForSubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\AnalogVelocityData\PlottingData1900s_VVplot.m""

######################################################
PyFrac experiments
######################################################


See PyFrac documentation on how to read the output files. 
To run simulations. Install PyFrac (check the tests run) then run: "PATH\\SpeedForSubmissionAdditionalData\PyFracScripts\TimVerticalPropagationCaseNrmLong.py"
To export the meshes for use in figures such as Fig.1 use:
"PATH\\SpeedForSubmissionAdditionalData\PyFracScripts\DrawSim_PyFrac.py"
and to get the speed at a given time use:
"PATH\\SpeedForSubmissionAdditionalData\PyFracScripts\GetFractureSpeed_PyFrac.py"
This saves a csv to the cwd with the heights, tip velocity etc.
