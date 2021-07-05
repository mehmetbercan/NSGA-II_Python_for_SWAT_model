
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
                'GenNumber': 20,
                'CrossPrb': 0.9,
                'CrossTyp': 2,
                'Bits': 16,
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
    
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------

    # check same members of child population / update
    population = copy.deepcopy(NSGAII.new_pop_ptr)
    nchrom = len(population["ind"][0]["xbin"])

    popsize = len(population["ind"])
    ParameterValues=[]
    for p in range(popsize):
        ParameterValues.append(population["ind"][p]["xbin"]) 
    
    # determine non unique pop members that needs to be updated (if matches two (or more), the last of them stays as unique)
    isUnique = []
    matched_indices = []; matched_2_indices = []
    for p1 in range(popsize):
        for p2 in range(p1+1, popsize):
            issame = False
            if p1 != p2:
                  issame = all(ParameterValues[p1] == ParameterValues[p2])
            if issame:
                matched_indices.append(p1)
                matched_2_indices.append(p2)
                break
        isUnique.append(not issame)
    indices_4_popmember_that_needs_update = matched_indices # matched_indices is same as [i for i, x in enumerate(isUnique) if x == False]
    # print ('non-unique replaced: {}'.format(len(indices_4_popmember_that_needs_update)))
    
    # update child population if necessary
    for nonunique_i in indices_4_popmember_that_needs_update:
        # randomly increase or decrese values 5% of its limit range
        for p in range(nchrom):
            par = NSGAII.new_pop_ptr["ind"][nonunique_i]["xbin"][p]
            if random.randint(0, 1) == 0:
                par = par + ((NSGAII.lim_b[p][1] - NSGAII.lim_b[p][0]) * 0.05)
            else:
                par = par - ((NSGAII.lim_b[p][1] - NSGAII.lim_b[p][0]) * 0.05)
            #make sure the new value is in parameter value limits
            if par < NSGAII.lim_b[p][0]:
                par = NSGAII.lim_b[p][0]
            if par > NSGAII.lim_b[p][1]:
                par = NSGAII.lim_b[p][1]
            NSGAII.new_pop_ptr["ind"][nonunique_i]["xbin"][p] = par
    
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------

    
    # calculate fitness
    CalculateObjectiveFunctions(NSGAII.new_pop_ptr);


    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------

    # put original pop member back in if updated pop member did not improve results (did not generate smaller summed objective function)
    for p1 in indices_4_popmember_that_needs_update:
        p2 = matched_2_indices[indices_4_popmember_that_needs_update.index(p1)]
        sum_fitness1 = sum(NSGAII.new_pop_ptr["ind"][p1]["fitness"])
        sum_fitness2 = sum(NSGAII.new_pop_ptr["ind"][p2]["fitness"])
        if sum_fitness1 >= sum_fitness2:
            NSGAII.new_pop_ptr["ind"][p1] = NSGAII.new_pop_ptr["ind"][p2]

    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    
    # ranking based on calculated objective functions
    nsga2utilities.rankcon(NSGAII.new_pop_ptr);
    
    NSGAII.CreateParentPopulation(i+1) # Old and New populations goes throuth Elitism, crowding distances, nondominated sorting
    #and create the old population for next generation. Report is printed during this function process.

    # plot
    df = nsga2utilities.CreatePopulationDataframe(NSGAII.mate_pop_ptr, NSGAII.parname)
    df[['f1', 'f2']].plot.scatter('f1','f2', title='Generation {}'.format(i+1), figsize=(6.5,6))
    
    i+=1


print("The NSGA-II execution finished. Look at the results in NSGA2.OUT folder.");

