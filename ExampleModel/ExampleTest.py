from nsga2lib import nsga2main, SWATnsga2funcs, nsga2funcs
import copy

#---- Input ----(space on directory may cause problems during execution)
SWATtxtinoutDirectory = r"C:\Users\ercanm\Desktop\GitHupRepositories\NSGA-II_Python_for_SWAT_model-version_2.0\ExampleModel\swatTest"
#---------------

ns=nsga2main.nsga2(SWATtxtinoutDirectory) 
ns.DetermineInitialPopulation()

#Loop through generations
TotalNumGenerations = ns.gener
i=0
while i < TotalNumGenerations:
    ns.Selection(); #select solutions from old population to mate population
    ns.old_pop_ptr = copy.deepcopy(ns.new_pop_ptr); #Now, old population becames new population before new population changes
    ns.Crossover() #New population changes based on mate population
    ns.Mutation()  #New population mutates
    nsga2funcs.decode(ns.new_pop_ptr, ns.vlen, ns.lim_b); #Turn binary calibration parameters into normal numbers
    SWATnsga2funcs.CalculateObjectiveFunctions(ns.new_pop_ptr,ns.Outlet_Obsdata,ns.FuncOpt,ns.FuncOptAvr,ns.parname,i+1,ns.SWATdir);
    ns.CreateMatePopulation(i+1) # Old and New populations goes throuth Elitism, crowding distances, nondominated sorting
    nsga2funcs.decode(ns.mate_pop_ptr, ns.vlen, ns.lim_b); 
    #print
    nsga2funcs.report(ns.old_pop_ptr,ns.mate_pop_ptr,i+1,ns.gener,ns.SWATdir,ns.ncross,ns.nmut);
    #copy mate to new population
    ns.new_pop_ptr = copy.deepcopy(ns.mate_pop_ptr);
    i+=1


print "The NSGA-II execution finished. Look at the results in NSGA2.OUT folder.";

