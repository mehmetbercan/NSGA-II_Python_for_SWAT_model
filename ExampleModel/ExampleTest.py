from nsga2lib import nsga2, SWATutilities, nsga2utilities


#---- Input ----(space on directory may cause problems during execution)
SWATtxtinoutDirectory = r"C:\Users\ercanm\Desktop\GitHupRepositories\NSGA-II_Python_for_SWAT_model\ExampleModel\swatTest"
#---------------

NSGAII=nsga2.nsga2(SWATtxtinoutDirectory) 
NSGAII.CreateInitialPopulation()

#Loop through generations
TotalNumGenerations = NSGAII.ngener
i=0
while i < TotalNumGenerations:
    NSGAII.CreateChildPopulation() #Thorough selection, crossover and mutation child population created from old population
    
    nsga2utilities.decode(NSGAII.new_pop_ptr, NSGAII.vlen, NSGAII.lim_b); #Turn binary calibration parameters into normal numbers.
    #new_pop_ptr=child population. vlen=the no.of bits assigned to the each calibration parameters. lim_b=range of calibration parameters.
    
    SWATutilities.CalculateObjectiveFunctions(NSGAII.new_pop_ptr,NSGAII.Outlet_Obsdata,NSGAII.FuncOpt,NSGAII.FuncOptAvr,NSGAII.parname,i+1,NSGAII.SWATdir);
    
    NSGAII.CreateParentPopulation(i+1) # Old and New populations goes throuth Elitism, crowding distances, nondominated sorting
    #and create the old population for next generation. Report is printed during this function process.
    i+=1


print "The NSGA-II execution finished. Look at the results in NSGA2.OUT folder.";

