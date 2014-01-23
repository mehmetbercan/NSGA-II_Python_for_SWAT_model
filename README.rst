####
NSGA-II for SWAT Watershed Model
####
nsga2 1.0.0

Released: 23-January-2014

************
Introduction
************
Soil and Water Assessment Tool (SWAT) is a conceptual distributed continuous time model that has capability
 of running on a daily and sub-daily time step.

Non-Dominated Sorting Genetic Algorithm II (NSGA-II)is a multi-objective optimization algorithm used as an 
automatic calibration tool in various disciplines.

There are Three pythoncodes in this folder: 
* Extract_rch.py---------> Extracts output.rch file based on ./NSGA2.IN/observed_rch.txt file.
* SWAT_ParameterEdit.py--> Edits SWAT parameters based on ./model.in file using Backup folder for default
 SWAT parameter values.
* NSGA-II_Bekele_Cloud.py> Runs the NSGA-II for SWAT model. When running, it uses ./nsga2_mid.cmd file and
 produces ./model.in file for parameter values.



To Setup Model:
*  Copy SWATtxtinout content in this folder and ./Backup folder.
*  Setup NSGA-II inputs in ./NSGA2.IN folder.
  * nsga2.def file--------> Defines NSGA-II methods setting.
  * nsga2_par.def file----> Defines parameters to be calibrated and their constrains.
  Currently there are 24 flow parameters can be calibrated.
  * observed_rch.txt file-> This file contains observation data. When there is missing data,
  it should not be included in this file. Column one is the squential numbers which should
  jump when there is missing data (e.g. 3,4,7,8,9,10 where the data for 1,2,5 and 6'th days
  are missing).



..Now, you can start code by running NSGA-II_Bekele_Cloud.py.
..For more information please contact with **Mehmet Ercan** at mehmetbercan@gmail.com or ercanm@engr.sc.edu

