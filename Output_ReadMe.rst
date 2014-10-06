**Model Output:**

In this section, the output files will be explained to make users familiar with them. At the end of this section, the example output files will be used for demonstration purposes.

*  Output.out:

  * This file contains all results starting from first generation to last generation (Last generation is the final result that included in next output file, *"plot.out"*).
  * This ouput is in same format with Prof. Deb's C code output. In our case, varibale columns indicate SWAT calibration parameter values.



*  Plot.out:

  * This file contains fitness results from the population of last generation (pareto front).
  * Each line, displaying objective function values, is a member of the pareto front. If two objective defined, the order in definition file is applied. if more than one site defined for objective functions, the order in observation file is applied. 



*  g_rank_record.out: 

  * This file contains only results related with NSGA-II recods and does not have any data related with SWAT model outputs.



*  Example output demonstration: 

  * **output.out file:** The first eight columns represent the calibration parameter values sperated with single space. The calibration paramers are in same order defined in *"nsga2_par.def"* file. The next four columns represent fitness values. There are two objective function defined in *"nsga2.def"* file (ObjFunc=5 --> E, PB) and two observation sites defined in *"observed_rch.txt"*. Therefore, the order of four fitness values are E for the first site, PB for the first site, E for the second site and PB for the second site. Next two columns are related with NSGA-II methods. *"|**|"* sign seperates previous population results from current population results. After *"|**|"* sign, the order of columns are same as just defined for the current population. 
  * **plot.out file:** There are four columns each representing a fitness value for a pareto front member. The order is same as defined above for fitness columns of *"output.out"* file. Therefore, the order of columns are as E for the first site, PB for the first site, E for the second site and PB for the second site.