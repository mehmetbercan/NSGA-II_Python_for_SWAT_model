#################################
NSGA-II for SWAT Watershed Model
#################################
nsga2sw 3.0

Released: 31-October-2014


************
Introduction
************
**Soil and Water Assessment Tool (SWAT)** is a conceptual distributed continuous
time model that has capability of running on a daily and sub-daily time step.

**Non-Dominated Sorting Genetic Algorithm II (NSGA-II)** is a multi-objective
optimization algorithm used as an automatic calibration tool in various disciplines.

************
Setup
************  
 
**Install the Python module:**

*  Python setuptools are required for installation
*  Download the repository
*  Open a command prompt and "cd" to "./nsga2lib"
*  Then command "Python setup.py install"
*  Alternatively, you can install directly from github using the command in the command prompt ```pip install https://github.com/mehmetbercan/NSGA-II_Python_for_SWAT_model/zipball/master```

**Setup SWAT Model:** 

*  Create *"./Backup"* folder in  *"SWATtxtinout"* folder and copy content in *"SWATtxtinout"* to the *"SWATtxtinout/Backup"* folder.
*  Copy *"ExampleModel/swatTest/NSGA2.IN"* in to the *"SWATtxtinout"* folder.
*  Setup NSGA-II inputs in *"./NSGA2.IN"* folder. 

  * **nsga2.def file:** Defines NSGA-II methods setting.
  * **nsga2_par.def file:** Defines parameters to be calibrated and their constrains.Currently there are 24 flow parameters can be calibrated.
  * **observed_rch.txt file:** This file contains observation data. When there is missing data,it should not be included in this file. Column one is the squential numbers which should jump when there is missing data (e.g. 3,4,7,8,9,10 where the data for 1,2,5 and 6'th days are missing).

*  NSGA-II calibration can be started using *"ExampleModel/ExampleTest.py"*.
*  *"ExampleModel/ExampleTest.py"* can be edited for certain purposes.


*  Definition of model Outputs in *"./NSGA2.OUT"* folder.

  * Once the the *"ExampleModel/ExampleTest.py"* finishes execution, it creates output files bellow.
  * **Output.out file:** Contains all results starting from first generation to last generation.
  * **Plot.out file:** Contains fitness results from the population of last generation (pareto front).
  * **g_rank_record.out file:** Related with NSGA-II recods.


**Running NSGA-II in SWAT-CUP:** 

Although above procedure is enough to run NSGA-II, this section explains how to run it in SWAT-CUP. This calibration tool is designed in a way to migrate into SWAT-CUP with an only need to add the NSGA-II calibration method to SWAT-CUP interface. The procedure bellow is not final intent but it is a way around to run NSGA-II in SWAT-CUP. 

*  Follow the above procedure but do not run (*"ExampleModel/ExampleTest.py"* )

*  Open SWAT-CUP and create a new GLUE project by following instruction within the interface

  * The "SWATtxtinout" used during this process is only for creating SWAT-CUP GLUE project and will not be used for any other purposes. Thus, using a small size model here would be faster.

*  Under SWAT-CUP "Project Explorer", extend "Executable Files" and right Click on "GLUE_run.bat" to open "item location"

*  Right click "GLUE_run.bat" file and click "edit"

  * Delete everything and paste the entire directory for *"ExampleTest.py"*, *"CreateGlueFiles4VisualizationInSWAT-CUP.py"* and *"MoveGlueFilesInCorrespondingSWAT-CUPfolder.py"*. This will tell the SWAT-CUP to run the NSGA-II instead of GLUE  in the NSGA-II directories, turn the results in to the GLUE format and move them into SWAT-CUP GLUE directories. 
  * Make sure the Python inputs (the SWAT and SWAT-CUP directories) for NSGA-II and GLUE folders are correctly set in *"ExampleTest.py"*, *"CreateGlueFiles4VisualizationInSWAT-CUP.py"* and *"MoveGlueFilesInCorrespondingSWAT-CUPfolder.py"* files (the last two are located under "ExternalPythonScripts" folder).

*  Click on "Calibrate..." (green gear shape) on the upper menu and click on "GLUE_run.bat"

  * Again, this starts NSGA-2 within its own directories and brings results into SWAT-CUP GLUE directories

*  To generate more results in SWAT-CUP, click on "Calibrate..." on the upper menu and click on "GLUE_Post.bat" 

*  Plots and outputs can be viewed under "Calibration Outputs" and "Sensitivity analysis" in the SWAT-CUP "Project Explorer" window to the left.
	
	
************
Notes
************ 


*  Read papers bellow to understand the process behing the scripts
*  Visit `my website <http://mehmetbercan.com/research/researchCal.html>`_ for more information
*  If you encounter any problems or have suggestions for the future development, please contact **Mehmet B. Ercan** at mehmetbercan@gmail.com or ercanm@engr.sc.edu.

**Credit:** 

Please cite one of the bellow articles if you use this code:

Ercan, M. B. and J. L. Goodall(2016), Design and implementation of a general software library for using NSGA-II with SWAT for multi-objective model calibration., *Environmental Modelling & Software*, 84, 112-120. doi:`10.1016/j.envsoft.2016.06.017 <http://www.sciencedirect.com/science/article/pii/S1364815216302547>`_.

Ercan, M. B. and J. L. Goodall (2014), A Python tool for multi-gage calibration of SWAT models using the NSGA-II algorithm., In: Ames, D.P., Quinn, N.W.T., Rizzoli, A.E. (Eds.), 2014. *Proceedings of the 7th International Congress on Environmental Modelling and Software, June 15-19, San Diego, California, USA*. (4):2325-2331, 2014. doi:`10.13140/2.1.3865.4407 <http://www.iemss.org/sites/iemss2014/papers/iemss2014_submission_212.pdf>`_. 


