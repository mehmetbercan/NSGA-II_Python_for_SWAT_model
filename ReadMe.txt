-- Non-Dominated Sorting Geneting Algorithm II (NSGA-II) for SWAT model --

There are Three pythoncodes in this folder: 
1. Extract_rch.py---------> Extracts output.rch file based on ./NSGA2.IN/observed_rch.txt file.
2. SWAT_ParameterEdit.py--> Edits SWAT parameters based on ./model.in file using Backup folder
	for default SWAT parameter values
3. NSGA-II_Bekele_Cloud.py> Runs the NSGA-II for SWAT model. When running, it uses ./nsga2_mid.cmd
	file and produces ./model.in file for parameter values.



To Setup Model:
1. Copy SWATtxtinout content in this folder and ./Backup folder.
2. Setup NSGA-II inputs in ./NSGA2.IN folder.
	2.1. nsga2.def file--------> Defines NSGA-II methods setting.
	2.2. nsga2_par.def file----> Defines parameters to be calibrated and their constrains.
		Currently there are 24 flow parameters can be calibrated.
	2.3. observed_rch.txt file-> This file contains observation data. When there is missing data,
		it should not be included in this file. Column one is the squential numbers which should
		jump when there is missing data (e.g. 3,4,7,8,9,10 where the data for 1,2,5 and 6'th days
		are missing).



Now, you can start code by running NSGA-II_Bekele_Cloud.py.
For more information please contact with Mehmet Ercan at mehmetbercan@gmail.com or ercanm@engr.sc.edu



