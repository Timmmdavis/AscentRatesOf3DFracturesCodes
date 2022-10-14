Before running any scripts add the folder this is contained within to your MATLAB path...

Cover figure:
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\DrawPyFracMeshesCoverFigure.m"

Main text:
Fig.1 
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\DrawPyFracMeshes.m"
Fig.2a
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\VelocityGravityFracturingResults2c.m" 
Fig.2b
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\VelocityGravityFracturingResults2cVsTime.m" 
Fig.3 and Fig.4
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\PyFracLargeScaleData\HeightHeightAndVelocityVelocityPlot.m"
Fig.5a
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\AnalogVelocityData\VelocityGravityGelatinTimeVsLen.m""
Fig.5b
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\AnalogVelocityData\VelocityGravityGelatinResultsHHplotInjTimes.m""

Appendix:
Fig.6a
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\VelocityGravityFracturingResults2cComparison.m" 
Fig.6b
"PATH\\SpeedForResubmissionAdditionalData\CodesToDrawFigures\AscentVelocity\VelocityGravityFracturingResults2cVsTimeComparison.m" 

######################################################
PyFrac experiments
######################################################


See PyFrac documentation on how to read the output files. 
To run simulations. Install PyFrac (check the tests run) then run: "PATH\\SpeedForResubmissionAdditionalData\PyFracScripts\TimVerticalPropagationCaseNrmLong.py"
To export the meshes for use in figures such as Fig.1 use:
"PATH\\SpeedForResubmissionAdditionalData\PyFracScripts\DrawSim_PyFrac.py"
and to get the speed at a given time use:
"PATH\\SpeedForResubmissionAdditionalData\PyFracScripts\GetFractureSpeed_PyFrac.py"
This saves a csv to the cwd with the heights, tip velocity etc.
