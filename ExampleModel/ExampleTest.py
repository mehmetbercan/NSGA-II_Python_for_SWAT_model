import os#, sys
#sys.path.append(os.path.join(os.getcwd(),"SWATnsga2Libs")) #Use this if you do not want to install library
from nsga2lib import nsga2, SWATutilities, nsga2utilities


#---- Input ----(space on directory may cause problems during execution)
SWATtxtinoutDirectory = os.path.join(os.getcwd(),"swatTest")
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
    
    if NSGAII.is_parallel:
        SWATutilities.CalculateObjectiveFunctionsinParallel(NSGAII.new_pop_ptr,NSGAII.Outlet_Obsdata,NSGAII.FuncOpt,NSGAII.FuncOptAvr,NSGAII.parname,i+1,NSGAII.SWATdir, NSGAII.cpu_number);
    else:
        SWATutilities.CalculateObjectiveFunctions(NSGAII.new_pop_ptr,NSGAII.Outlet_Obsdata,NSGAII.FuncOpt,NSGAII.FuncOptAvr,NSGAII.parname,i+1,NSGAII.SWATdir);
    
    NSGAII.CreateParentPopulation(i+1) # Old and New populations goes throuth Elitism, crowding distances, nondominated sorting
    #and create the old population for next generation. Report is printed during this function process.
    i+=1


print("The NSGA-II execution finished. Look at the results in NSGA2.OUT folder.");

