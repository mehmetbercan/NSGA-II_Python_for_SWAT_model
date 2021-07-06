
#-------------------------------------------------------------------------------
# Name:        NSGA-II
# Purpose:     Multi-Objective calibration of SWAT model
#
# Author:      Mehmet B. Ercan (mehmetbercan@gmail.com) at USC/Columbia, SC.
#
# Created:     10/29/2014
# Edited:      12/29/2020
# Copyright:   (c) Mehmet B. Ercan 2014
# Licence:     MIT
#-------------------------------------------------------------------------------

'''
P.S: Old population (old_pop_ptr) is a parent population, new population (new_pop_ptr)
is a child population and mate population (mate_pop_ptr) is an intermediate population between
transactions. General process starts with an old population which goes through selection to
create mate population which is then used to form a new population through crossover. This
new population goes through mutation to get final new population. Then, new and old populations
go through elitism, crowding distances, nondominated sorting to create a mate population.
This mate population then copied as an old population for the next generation. The same porcess
repeats for every generation.
'''

import copy, sys, numpy, os, random, shutil
from nsga2lib import nsga2utilities
class nsga2:
    def __init__(self,setting_dict, para_dict, TestDir):

        popsize = setting_dict['PopSize'] #Population size (an even no.)
        ngener = setting_dict['GenNumber'] #the no.of generations
        pcross = setting_dict['CrossPrb'] #the cross-over probability (between 0.5 and 1)
        optype = setting_dict['CrossTyp'] #Crossover type 1 for Simple one & 2 for Uniform X-over
        nbits = setting_dict['Bits']  #No.of bits assigned to each variable(parameters)
        pmutprb = setting_dict['MutPrb'] #Mutation prblty proportion range(between 0-1)
        seed = setting_dict['seed'] #random seed(between 0 and 1)
        M = setting_dict['M']      #Number of Latin Hypercube Sampling intervals
        nfunc = setting_dict['ObjFuncNum'] #Number of Objective Functions
        ReadMFrmOut = setting_dict['ReadMFrmOut'] #1= Read last population from output.out (use "1" when you want to re-start with same parameters defined in parameter file)
        UniqueParSetSize = 1000 #number of unique parameter set to be saved / used to replace duplicates (default = 1000 - will record last 1000 unique pop members)


        nchrom=len(para_dict);
        if nchrom <= 0: sys.exit("ERROR: 'nsga2_par.txt' files does not have prameters (or paramters doesn't start with 'a','r' or 'v')")
        chrom = 0#Chromosome length (Total Sum of the bit value)
        vlen=[];lim_b=[]; parname = []
        for key, val in para_dict.items():
            vlen.append(nbits) #the no.of bits assigned to the %d variable\n"%(i+1)
            chrom += nbits;
            parname.append(key)
            lim_b.append(val) #lower & the upper limits of the %d variable\n"%(i+1)
        pmut_b = pmutprb * 1.0/chrom #the mutation probability for binary strings (between 0 and 1.0/chrom)
        
        
        self.Modeldir=TestDir
        
        self.popsize=popsize
        self.ngener=ngener
        self.pcross=pcross
        self.optype=optype
        self.M=M
        self.ReadMFrmOut=ReadMFrmOut
        self.UniqueParSetSize=UniqueParSetSize

        self.chrom=chrom
        self.nchrom=nchrom
        self.nfunc=nfunc
        self.pmut_b=pmut_b
        self.lim_b=lim_b
        self.vlen=vlen
        self.parname=parname
        
        #---
        self.nmut=0
        self.ncross=0
        self.old_pop_ptr=nsga2utilities.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        self.new_pop_ptr=nsga2utilities.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        self.mate_pop_ptr=nsga2utilities.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        self.historic_record={'Parameters':[], 'Fitnesses':[]}
        #/*Initialize the random no generator*/
        self.warmup_random = random_(seed); #
    #-------------------------------------------------------------------------------


    #-------------------------------------------------------------------------------
    def CreateInitialPopulation(self, CalculateObjectiveFunctions):
        old_pop_ptr=self.old_pop_ptr
        #Check if NSGA2.OUT exist
        if not os.path.exists(self.Modeldir+"/NSGA2.OUT"): os.makedirs(self.Modeldir+"/NSGA2.OUT")
       
        #Define the the initital population
        if self.ReadMFrmOut == 1:#Read Last population from output.out
            #Read outputout
            f = open(self.Modeldir+"/NSGA2.OUT/output.out","r")
            lines = f.readlines()
            prmtrno = int(lines[5].split(")")[0].split("binary")[1]) #Number of parameters (binary)
            if self.nchrom != prmtrno: sys.exit("ERROR: parameter number is not equal to parameter numbers defined in output.out file")
            #Loop through lines
            i=6; LastGenPars=[]; gennum=0; previousgennum=0
            while i < len(lines):
                i += 2
                #Deal with generation title part and get generation number
                if lines[i] == "\n":
                    try:gennum = int(lines[i+5].split("\n")[0].split("->")[1])
                    except:pass
                    i += 10
                    if i >= len(lines):#break if the loop is completed
                        break
                if gennum > previousgennum:
                    previousgennum = gennum
                    LastGenPars=[]
                Parameters = [x for x in lines[i].split("\n")[0].split("|**|")[1].split(" ") if x!=""][:prmtrno] #Get parameters on line i after |**| (mate population)
                LastGenPars.append(Parameters)  
            f.close()
            if self.popsize != len(LastGenPars): sys.exit("ERROR: poulation size is not equal to population size in output.out file")
            #Write values in old_pop_ptr
            for i in range(0,self.popsize):
                for j in range(0,self.nchrom):
                    old_pop_ptr["ind"][i]["xbin"][j] = LastGenPars[i][j]
            #Copy old output file
            shutil.copy2(self.Modeldir+"/NSGA2.OUT/output.out", self.Modeldir+"/NSGA2.OUT/output_previous.out")            
            #/*Function Calculaiton*/
            CalculateObjectiveFunctions(old_pop_ptr)
            self.Record_Unique_Population_Sets(old_pop_ptr)
        else:
            #Defining Latin Hypercube Sampling population
            InitialLHSpop = nsga2utilities.CreateDefaultPopulation(self.M,self.chrom,self.nchrom,self.nfunc)
            #Latin Hypercube Samples
            LHSamples = []
            for j in range(0,len(self.lim_b)):
                LowBound = self.lim_b[j][0]
                UpBound = self.lim_b[j][1]
                parVals = []
                for i in range(0,self.M+1):
                    point = i * 1. / self.M
                    parValue = (point * (UpBound - LowBound)) + LowBound
                    parVals.append(parValue)
                LHSamples.append(parVals)
            #Random selection of values from each interval of LHS
            for i in range(0,self.M):
                for j in range(0,self.nchrom):
                    rndinteger = int(round(self.M*random.random(),0))
                    InitialLHSpop["ind"][i]["xbin"][j] = LHSamples[j][rndinteger]         
            #/*Function Calculaiton*/
            CalculateObjectiveFunctions(InitialLHSpop)
            self.Record_Unique_Population_Sets(InitialLHSpop)
            #Select the first population from InitialLHSpop
            n=0
            for rank in range(1,self.M+1):
                for i in range(0,self.M):
                    if rank >= InitialLHSpop["ind"][i]["rank"] and InitialLHSpop["ind"][i]["flag"] != 99:
                        old_pop_ptr["ind"][n] = copy.deepcopy(InitialLHSpop["ind"][i])
                        n+=1
                        InitialLHSpop["ind"][i]["flag"] = 99
                        if n == self.popsize:
                            break
                else:continue
                break
        nsga2utilities.reverse_decode(old_pop_ptr,self.vlen,self.lim_b)
        nsga2utilities.rankcon(old_pop_ptr)
        self.old_pop_ptr=old_pop_ptr
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    def Record_Unique_Population_Sets(self, population):
        # check non-unique members of child population / update
        ParameterValues=[]; Objectives=[]
        for p in range(self.popsize):
            ParameterValues.append(copy.deepcopy(population["ind"][p]["xbin"])) 
            Objectives.append(copy.deepcopy(population["ind"][p]["fitness"])) 
        
        # determine non unique pop members 
        isUnique = []
        for p1 in range(self.popsize):
            for p2 in range(p1+1, self.popsize):
                issame = False
                if p1 != p2:
                    issame = all(ParameterValues[p1] == ParameterValues[p2])
                if issame:
                    break
            isUnique.append(not issame)
        Unique_indices = [i for i, x in enumerate(isUnique) if x == True]

        for p in Unique_indices:
            self.historic_record['Parameters'].append(ParameterValues[p])
            self.historic_record['Fitnesses'].append(Objectives[p])
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    def Update_nonUnique_childpop_members(self):
        population = copy.deepcopy(self.new_pop_ptr)

        # check non-unique members of child population / update
        ParameterValues=[]; 
        for p in range(self.popsize):
            ParameterValues.append(population["ind"][p]["xbin"])
        
        # determine non unique pop members 
        matched_indices = []; matched_2_indices = []
        matched_indices_historic = []; matched_2_indices_historic = []
        historic_ParameterValues = self.historic_record['Parameters']
        for p1 in range(self.popsize):
            for p2 in range(len(historic_ParameterValues)):
                issame = False
                if p1 != p2:
                    issame = all(ParameterValues[p1] == historic_ParameterValues[p2])
                if issame:
                    matched_indices_historic.append(p1)
                    matched_2_indices_historic.append(p2)
                    break
            if not issame:
                for p2 in range(p1+1, self.popsize):
                    issame = False
                    if p1 != p2:
                        issame = all(ParameterValues[p1] == ParameterValues[p2])
                    if issame:
                        matched_indices.append(p1)
                        matched_2_indices.append(p2)
                        break

        indices_4_popmember_that_needs_update = matched_indices + matched_indices_historic
        
        # update child population if necessary
        for nonunique_i in indices_4_popmember_that_needs_update:
            # randomly change one parameter
            par_2_change_index = random.randint(0, self.nchrom-1)
            upper_threshold = self.lim_b[par_2_change_index][1]
            lower_threshold = self.lim_b[par_2_change_index][0]
            random_par_value = random.uniform(lower_threshold, upper_threshold)
            self.new_pop_ptr["ind"][nonunique_i]["xbin"][par_2_change_index] = random_par_value

        return matched_indices, matched_2_indices, matched_indices_historic, matched_2_indices_historic
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    def CreateChildPopulation(self):
        #Selection
        nsga2utilities.Selection(self.old_pop_ptr,self.mate_pop_ptr,self.warmup_random)
        #Crossover
        if(self.optype == 1): #/*Binary Cross-over*/
            self.ncross = nsga2utilities.crossover(self.new_pop_ptr,self.mate_pop_ptr,self.warmup_random,self.pcross,self.ncross);

        if(self.optype == 2):#/*Binary Uniform Cross-over*/
            self.ncross = nsga2utilities.unicross(self.new_pop_ptr ,self.mate_pop_ptr,self.warmup_random,self.pcross,self.ncross);
        #Mutation
        self.nmut = nsga2utilities.Mutation(self.new_pop_ptr,self.warmup_random,self.pmut_b,self.nmut)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    def CreateParentPopulation(self,generationNo):
        nsga2utilities.CreateMatePopFromNewandOldPops(self.old_pop_ptr ,self.new_pop_ptr ,self.mate_pop_ptr,generationNo, self.Modeldir)
        nsga2utilities.decode(self.mate_pop_ptr, self.vlen, self.lim_b);
        nsga2utilities.report(self.old_pop_ptr,self.mate_pop_ptr,generationNo,self.ngener,self.Modeldir,self.ncross,self.nmut); #Print Report: old_pop_ptr is old population and mate_pop_ptr is here the created old population for next generation
        self.new_pop_ptr = copy.deepcopy(self.mate_pop_ptr); 
        self.old_pop_ptr = copy.deepcopy(self.new_pop_ptr); #Parent Population
    #-------------------------------------------------------------------------------















#-------------------------------------------------------------------------------
def advance_random(oldrand):
    #/* Create next batch of 55 random numbers */
    for j1 in range(24):
        new_random = oldrand[j1] - oldrand[j1+31];
        if(new_random < 0.0): new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    for j1 in range(24,55):
        new_random = oldrand [j1] - oldrand [j1-24];
        if(new_random < 0.0): new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    return oldrand

class random_:
    #/* Get random off and running */
    def __init__(self,random_seed): #Warmup_random(random_seed):
        oldrand = numpy.zeros(55,float)
        oldrand[54]=random_seed
        new_random = 0.000000001;
        prev_random = random_seed;
        for j1 in range(1,55):
            ii = (21*j1)%54;
            oldrand[ii] = new_random;
            new_random = prev_random-new_random;
            if(new_random<0.0): new_random = new_random + 1.0;
            prev_random = oldrand[ii];
        self.oldrand = advance_random(advance_random(advance_random(oldrand)));
        self.jrand = 0

    def randomperc(self):
    #/* Fetch a single random number between 0.0 and 1.0 - Subtractive Method */
    #/* See Knuth, D. (1969), v. 2 for details */
    #/* name changed from random() to avoid library conflicts on some machines*/
        self.jrand +=1
        if(self.jrand >= 55):
            self.jrand = 1;
            self.oldrand = advance_random(self.oldrand);
        return self.oldrand[self.jrand]
#-------------------------------------------------------------------------------
