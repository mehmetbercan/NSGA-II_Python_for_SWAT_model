import os#, sys
#sys.path.append(os.path.join(os.getcwd(),"SWATnsga2Libs")) #Use this if you do not want to install library




#-------------------------------------------------------------------------------
# Name:        NSGA-II
# Purpose:     Calculates objective functions for NSGA-II
#
# Author:      Mehmet B. Ercan (mehmetbercan@gmail.com)
#
# Created:     10/29/2014
# Edited:      01/11/2017
# Copyright:   (c) Mehmet B. Ercan 2014
# Licence:     MIT
#-------------------------------------------------------------------------------


import os
from nsga2lib import nsga2, nsga2utilities

HERE = os.path.dirname(os.path.realpath(__file__))

# Objective Function Calculator: Equvalent of Model Run and Calculate Objective Function
def SCHproblem_obj_func(parvals):
    '''
    SCH is a one parameter and two fitness problem. Optimal solution is within
    0 and 2. Upper and lower limits are -10**3 and 10**3, respectively.
    '''
    x = parvals[0]
    f1 = x**2
    f2 = (x - 2)**2
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
        objectivefuncs = SCHproblem_obj_func(parvals)
                      
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
setting_dict = {'PopSize': 150,
                'GenNumber': 15,
                'CrossPrb': 0.9,
                'CrossTyp': 2,
                'Bits': 15,
                'MutPrb': 0.5,
                'seed': 0.5,
                'ObjFuncNum': 2,
                'M': 1000,
                'ReadMFrmOut': 0}
# Define Parameters and their limits
para_dict = {"PARAM1": [-10**3, 10**3]}
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

    
    # calculate fitness
    CalculateObjectiveFunctions(NSGAII.new_pop_ptr);
    
    # ranking based on calculated objective functions
    nsga2utilities.rankcon(NSGAII.new_pop_ptr);
    
    NSGAII.CreateParentPopulation(i+1) # Old and New populations goes throuth Elitism, crowding distances, nondominated sorting
    #and create the old population for next generation. Report is printed during this function process.

    # plot
    df = nsga2utilities.CreatePopulationDataframe(NSGAII.mate_pop_ptr, NSGAII.parname)
    df[['f1', 'f2']].plot.scatter('f1','f2', title='Generation {}'.format(i))
    
    i+=1


print("The NSGA-II execution finished. Look at the results in NSGA2.OUT folder.");

