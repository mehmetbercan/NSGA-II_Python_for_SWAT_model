#################################
NSGA-II for SWAT Watershed Model
#################################
nsga2sw 2.0

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
*  Open a command prompt and "cd" to "./nsga2lib"
*  Then command "Python setup.py install"
 
**Setup SWAT Model:** 

*  Create *"./Backup"* folder in  *"SWATtxtinout"* folder and copy content in *"SWATtxtinout"* to the *"SWATtxtinout/Backup"* folder.
*  Coppy *"ExampleModel\swatTest\NSGA2.IN"* in to the *"SWATtxtinout"* folder.
*  Setup NSGA-II inputs in *"./NSGA2.IN"* folder. 

  * **nsga2.def file:** Defines NSGA-II methods setting.
  * **nsga2_par.def file:** Defines parameters to be calibrated and their constrains.Currently there are 24 flow parameters can be calibrated.
  * **observed_rch.txt file:** This file contains observation data. When there is missing data,it should not be included in this file. Column one is the squential numbers which should jump when there is missing data (e.g. 3,4,7,8,9,10 where the data for 1,2,5 and 6'th days are missing).

*  NSGA-II calibration can be started using *"ExampleModel\ExampleTest.py"*.
*  *"ExampleModel\ExampleTest.py"* can be edited for certain purposes.


*  Definition of model Outputs in *"./NSGA2.OUT"* folder.

  * Once the the *"ExampleModel\ExampleTest.py"* finishes execution, it creates output files bellow.
  * **Output.out file:** Contains all results starting from first generation to last generation.
  * **Plot.out file:** Contains fitness results from the population of last generation (pareto front).
  * **g_rank_record.out file:** Related with NSGA-II recods.

Visit `my website <http://mehmetbercan.com/research/researchCal.html>`_ for more information.

If you encounter any problems or have suggestions for the future development, 
please contact **Mehmet B. Ercan** at mehmetbercan@gmail.com or ercanm@engr.sc.edu.




