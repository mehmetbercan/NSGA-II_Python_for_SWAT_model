
#-------------------------------------------------------------------------------
# Name:        NSGA-II
# Purpose:     Calculates objective functions for NSGA-II
#
# Author:      Mehmet B. Ercan (mehmetbercan@gmail.com)
#
# Created:     12/29/2020
# Edited:      
# Copyright:   (c) Mehmet B. Ercan 2014
# Licence:     MIT
#-------------------------------------------------------------------------------


import os, copy, random
#sys.path.append(os.path.join(os.getcwd(),"SWATnsga2Libs")) #Use this if you do not want to install library
from nsga2lib import nsga2, nsga2utilities

HERE = os.path.dirname(os.path.realpath(__file__))

# Objective Function Calculator: Equvalent of Model Run and Calculate Objective Function
def ZTD2_problem_obj_func(parvals):
    '''
    This is a problem that is used to test GA
    '''
    n = len(parvals)
    x1 = parvals[0]
    g = 1+(9*sum(parvals[1:]))/(n-1)
    f1 = x1
    f2 = g*(1-(x1/g)**2)
    return [f1, f2]

def CalculateObjectiveFunctions(population):
    """
    This function must be modified by user. Fitness values in population dictionary has to be updated here.
    Objective function values are considered to be best once they get closer to zero and worst when they are closer to +1 (or +infinity)
    """
    #/*Initializing the max rank to zero*/
    population["maxrank"]=0
    popsize = len(population["ind"])
    nfunc = len(population["ind"][0]["fitness"])
    nchrom = len(population["ind"][0]["xbin"])

    ParameterValues=[]
    for i in range(popsize):
        ParameterValues.append(population["ind"][i]["xbin"]) #/* problem variables */
        

#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#This part should be edited by the user based on the specific model to be calibrated.
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    #Population Loop: Paralization can be  applied here
    for i in range(popsize): 
        #Edit model input parameters using ParameterValues[i]
        parvals=ParameterValues[i]
  
        
        #Calculate Objective functions
        objectivefuncs = ZTD2_problem_obj_func(parvals)
                      
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#The user should stop editing after this line.
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #Add objective functions to population
        population["ind"][i]['fitness'] = objectivefuncs
    return;
#-------------------------------------------------------------------------------
    
# --------------------------------------------------------------------
# --------------------------------------------------------------------
# Define NSGA2 Settings
setting_dict = {'PopSize': 100,
                'GenNumber': 250,
                'CrossPrb': 0.9,
                'CrossTyp': 2,
                'Bits': 6,
                'MutPrb': 0.5,
                'seed': 0.5,
                'ObjFuncNum': 2,
                'M': 100,
                'ReadMFrmOut': 0}
# Define Parameters and their limits
para_dict = {'PARAM1':[0,1]}
for i in range(2, 31):
    para_dict['PARAM{}'.format(i)] = [0,1]
# --------------------------------------------------------------------
# --------------------------------------------------------------------



NSGAII=nsga2.nsga2(setting_dict, para_dict, HERE) 
NSGAII.CreateInitialPopulation(CalculateObjectiveFunctions)

#Loop through generations
TotalNumGenerations = NSGAII.ngener
i=0
while i < TotalNumGenerations:
    '''
    INFO:
        NSGAII.new_pop_ptr=child population; 
        NSGAII.vlen=the no.of bits assigned to the each calibration parameters; 
        NSGAII.lim_b=range of calibration parameters (upper and lower bounds).
    '''
    
    
    print ('Running generation: {} ...'.format(i))
    
    #Thorough selection, crossover and mutation child population created from old population
    NSGAII.CreateChildPopulation() 
    
    #Turn binary calibration parameters into normal numbers
    nsga2utilities.decode(NSGAII.new_pop_ptr, NSGAII.vlen, NSGAII.lim_b); 
    
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # update non unique child population (new_pop_ptr) members
    matched_indices, matched_2_indices, matched_indices_historic, matched_2_indices_historic = NSGAII.Update_nonUnique_childpop_members()
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    
    # calculate fitness
    CalculateObjectiveFunctions(NSGAII.new_pop_ptr);

    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # put original pop member back in if updated pop member did not improve results (did not generate smaller summed objective function)
    for p1 in matched_indices:
        p2 = matched_2_indices[matched_indices.index(p1)]
        sum_fitness1 = sum(NSGAII.new_pop_ptr["ind"][p1]["fitness"])
        sum_fitness2 = sum(NSGAII.new_pop_ptr["ind"][p2]["fitness"])
        if sum_fitness1 >= sum_fitness2:
            NSGAII.new_pop_ptr["ind"][p1] = copy.deepcopy(NSGAII.new_pop_ptr["ind"][p2])
    for p1 in matched_indices_historic:
        p2 = matched_2_indices_historic[matched_indices_historic.index(p1)]
        sum_fitness1 = sum(NSGAII.new_pop_ptr["ind"][p1]["fitness"])
        sum_fitness2 = sum(NSGAII.historic_record["Fitnesses"][p2])
        if sum_fitness1 >= sum_fitness2:
            NSGAII.new_pop_ptr["ind"][p1]['xbin'] = copy.deepcopy(NSGAII.historic_record["Parameters"][p2])
            NSGAII.new_pop_ptr["ind"][p1]['fitness'] = copy.deepcopy(NSGAII.historic_record["Fitnesses"][p2])

    # record unique values into historic dictionary
    NSGAII.Record_Unique_Population_Sets(NSGAII.new_pop_ptr)
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------

    # ranking based on calculated objective functions
    nsga2utilities.rankcon(NSGAII.new_pop_ptr);
    
    NSGAII.CreateParentPopulation(i+1) # Old and New populations goes throuth Elitism, crowding distances, nondominated sorting
    #and create the old population for next generation. Report is printed during this function process.

    # plot
    if (i+1)%25==0:
        df = nsga2utilities.CreatePopulationDataframe(NSGAII.mate_pop_ptr, NSGAII.parname)
        df[['f1', 'f2']].plot.scatter('f1','f2', title='Generation {}'.format(i+1), figsize=(6.5,6))
    
    i+=1


print("The NSGA-II execution finished. Look at the results in NSGA2.OUT folder.");

