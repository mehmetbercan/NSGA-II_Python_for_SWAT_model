
#-------------------------------------------------------------------------------
# Name:        NSGA-II
# Purpose:     Multi-Objective calibration of SWAT model
#
# Author:      Mehmet B. Ercan (mehmetbercan@gmail.com) at USC/Columbia, SC.
#
# Created:     10/29/2014
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

import copy, sys, numpy, os, random, shutil, nsga2utilities, SWATutilities
from math import floor

class nsga2:
    def __init__(self,SWATtxtinoutFolderDirectory):
        """nsga2 calibration functions"""
        """SWATtxtinout folder should have 'nsga.in' subfolder with nsga2 input files"""
        libpath = os.path.dirname(nsga2utilities.__file__)
        #Copy necessary files to SWAT directory (Operating Platform Specific)
        shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","Extract_rch.py"), SWATtxtinoutFolderDirectory)
        shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","SWAT_ParameterEdit.py"), SWATtxtinoutFolderDirectory)
        if "win" in sys.platform.lower():
            print ("Operating System is {0}".format(sys.platform))
            shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","nsga2_mid.cmd"), SWATtxtinoutFolderDirectory)
            shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","swat.exe"), SWATtxtinoutFolderDirectory)
        elif "linux" in sys.platform.lower():
            print ("Operating System is {0}".format(sys.platform))
            shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","nsga2_mid.sh"), SWATtxtinoutFolderDirectory)
            shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","swat2012_627"), SWATtxtinoutFolderDirectory)
            shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","Makefile"), SWATtxtinoutFolderDirectory)
        else: #goes with windows for now
            shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","nsga2_mid.cmd"), SWATtxtinoutFolderDirectory)
            shutil.copy2(os.path.join(libpath,"ScriptsForSWATtxt","swat.exe"), SWATtxtinoutFolderDirectory)
            
            
        #Read ('nsga2.def') NSGA-II binary options input
        nsga2def = SWATtxtinoutFolderDirectory+"/NSGA2.IN/nsga2.def"
        nsga2pardef = SWATtxtinoutFolderDirectory+"/NSGA2.IN/nsga2_par.def"
        observed_rch = SWATtxtinoutFolderDirectory+"/NSGA2.IN/observed_rch.txt"
        f = open(nsga2def, "r")
        lines = f.readlines()
        popsize = int(lines[1].split()[1]) #Population size (an even no.)
        ngener = int(lines[2].split()[1]) #the no.of generations
        pcross = float(lines[3].split()[1]) #the cross-over probability (between 0.5 and 1)
        optype = int(lines[4].split()[1]) #Crossover type 1 for Simple one & 2 for Uniform X-over
        bits = int(lines[5].split()[1])  #No.of bits assigned to each variable(parameters)
        pmutprp = float(lines[6].split()[1]) #Mutation prblty proportion range(between 0-1)
        seed = float(lines[7].split()[1]) #random seed(between 0 and 1)
        M = int(lines[8].split()[1])       #Number of Latin Hypercube Sampling intervals
        FuncOpt = int(lines[9].split()[1]) #1=E; 2=R^2; 3=E,log E; 4=E,R^2,log E; 5=E,PB
        FuncOptAvr = int(lines[10].split()[1]) #0=Do not average; 1=Average Objective sites; 2=Average Objective Functions
        ReadMFrmOut = int(lines[11].split()[1]) #1= Read last population from output.out (use "1" when you want to re-start with same parameters defined in parameter file)
        f.close()
        #Read 'nsga2_par.def'
        f = open(nsga2pardef, "r")
        lines = f.readlines(); nchrom=0;
        for i in range(1,len(lines)):
            if lines[i][0]=="\n" or lines[i][0]=="-":break;
            nchrom = i #no. of binary-coded variables (--number of parameters--)
        if nchrom <= 0: sys.exit("ERROR: 'nsga2_par.def' files does not have prameters (or paramters doesn't start with 'a','r' or 'v')")
        chrom = 0#Chromosome length (Total Sum of the bit value)
        vlen=[];lim_b=[]; parname = []
        for i in xrange(nchrom):
            vlen.append(bits) #the no.of bits assigned to the %d variable\n"%(i+1)
            chrom += bits;
            parname.append(lines[i+1].split()[0])
            lim_b.append([float(lines[i+1].split()[1]),float(lines[i+1].split()[2])]) #lower & the upper limits of the %d variable\n"%(i+1)
        pmut_b = pmutprp * 1.0/chrom #the mutation probability for binary strings (between 0 and 1.0/chrom)
        f.close()
        #Read 'observed_rch.txt'
        f = open(observed_rch, "r")
        lines = f.readlines()
        #Read Observed Streamflow
        Outlet_Obsdata = {}; outlet = -99; nofdatapoint=-99;
        for i in range(0,len(lines)):
            try:
                if lines[i][0:10:]=='output_rch':
                    outlet = int(lines[i][11:16].split(" ")[0])
                    nofdatapoint = int(lines[i+1].split(" ")[0])
                    Obsdata = []
                    for j in range((i+3), (i+3+nofdatapoint)):
                        Obsdata.append(float((lines[j].split("\t")[2]).split("\n")[0]))
                    Outlet_Obsdata[outlet] = Obsdata
            except: sys.exit("ERROR: check the 'observed_rch.txt' file (the gage on line:"+str(i)+")")
        #Calculate number of objective functions
        nfunc = 0 #no. of objective functions
        nsites = int(lines[0].split("  ")[0])
        if FuncOpt==1 or FuncOpt==2:
            nfunc = nsites            
        if FuncOpt==3 or FuncOpt==5:
            nfunc = 2*nsites
        if FuncOpt==4:
            nfunc = 3*nsites
        if FuncOptAvr==1:
            nfunc=nfunc/nsites
        if FuncOptAvr==2:
            nfunc=nsites
        f.close()
        self.parname=parname
        self.popsize=popsize
        self.ngener=ngener
        self.pcross=pcross
        self.optype=optype
        self.seed=seed
        self.M=M
        self.FuncOpt=FuncOpt
        self.FuncOptAvr=FuncOptAvr
        self.ReadMFrmOut=ReadMFrmOut
        self.nchrom=nchrom
        self.vlen=vlen
        self.chrom=chrom
        self.lim_b=lim_b
        self.pmut_b=pmut_b
        self.Outlet_Obsdata=Outlet_Obsdata
        self.nfunc=nfunc
        self.SWATdir=SWATtxtinoutFolderDirectory
        #---
        self.nmut=0
        self.ncross=0
        self.old_pop_ptr=nsga2utilities.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        self.new_pop_ptr=nsga2utilities.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        self.mate_pop_ptr=nsga2utilities.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        #/*Initialize the random no generator*/
        self.warmup_random = random_(seed); #
    #-------------------------------------------------------------------------------


    #-------------------------------------------------------------------------------
    def CreateInitialPopulation(self):
        old_pop_ptr=self.old_pop_ptr
        #Check if NSGA2.OUT exist
        if not os.path.exists(self.SWATdir+"/NSGA2.OUT"): os.makedirs(self.SWATdir+"/NSGA2.OUT")
       
        #Define the the initital population
        if self.ReadMFrmOut == 1:#Read Last population from output.out
            #Read outputout
            f = open(self.SWATdir+"/NSGA2.OUT/output.out","r")
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
            for i in xrange(self.popsize):
                for j in xrange(self.nchrom):
                    old_pop_ptr["ind"][i]["xbin"][j] = LastGenPars[i][j]
            #Copy old output file
            shutil.copy2(self.SWATdir+"/NSGA2.OUT/output.out", self.SWATdir+"/NSGA2.OUT/output_previous.out")            
            #/*Function Calculaiton*/
            SWATutilities.CalculateObjectiveFunctions(old_pop_ptr,self.Outlet_Obsdata,self.FuncOpt,self.FuncOptAvr,self.parname,"Previous NSGA-II run Last Population",self.SWATdir)
            
        else:
            #Defining Latin Hypercube Sampling population
            InitialLHSpop = nsga2utilities.CreateDefaultPopulation(self.M,self.chrom,self.nchrom,self.nfunc)
            #Latin Hypercube Samples
            LHSamples = []
            for j in xrange(len(self.lim_b)):
                LowBound = self.lim_b[j][0]
                UpBound = self.lim_b[j][1]
                parVals = []
                for i in xrange(self.M+1):
                    point = i * 1. / self.M
                    parValue = (point * (UpBound - LowBound)) + LowBound
                    parVals.append(parValue)
                LHSamples.append(parVals)
            #Random selection of values from each interval of LHS
            for i in xrange(self.M):
                for j in xrange(self.nchrom):
                    rndinteger = int(round(self.M*random.random(),0))
                    InitialLHSpop["ind"][i]["xbin"][j] = LHSamples[j][rndinteger]         
            #/*Function Calculaiton*/
            SWATutilities.CalculateObjectiveFunctions(InitialLHSpop,self.Outlet_Obsdata,self.FuncOpt,self.FuncOptAvr,self.parname,"LHS",self.SWATdir)
            #Select the first population from InitialLHSpop
            n=0
            for rank in range(1,self.M+1):
                for i in xrange(self.M):
                    if rank >= InitialLHSpop["ind"][i]["rank"] and InitialLHSpop["ind"][i]["flag"] != 99:
                        old_pop_ptr["ind"][n] = copy.deepcopy(InitialLHSpop["ind"][i])
                        n+=1
                        InitialLHSpop["ind"][i]["flag"] = 99
                        if n == self.popsize:
                            break
                else:continue
                break
        nsga2utilities.reverse_decode(old_pop_ptr,self.vlen,self.lim_b)
        SWATutilities.rankcon(old_pop_ptr)
        self.old_pop_ptr=old_pop_ptr
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
        nsga2utilities.CreateMatePopFromNewandOldPops(self.old_pop_ptr ,self.new_pop_ptr ,self.mate_pop_ptr,generationNo,self.SWATdir)
        nsga2utilities.decode(self.mate_pop_ptr, self.vlen, self.lim_b);
        nsga2utilities.report(self.old_pop_ptr,self.mate_pop_ptr,generationNo,self.ngener,self.SWATdir,self.ncross,self.nmut); #Print Report: old_pop_ptr is old population and mate_pop_ptr is here the created old population for next generation
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


