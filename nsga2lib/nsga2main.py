
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

import copy, sys, numpy, os, random, shutil, nsga2funcs, SWATnsga2funcs
from math import floor

class nsga2:
    def __init__(self,SWATtxtinoutFolderDirectory):
        """nsga2 calibration functions"""
        """SWATtxtinout folder should have 'nsga.in' subfolder with nsga2 input files"""
        libpath = os.path.dirname(nsga2funcs.__file__)
        #Copy necessary files to SWAT directory
        shutil.copy2(libpath+"/ScriptsForSWATtxt/Extract_rch.py", SWATtxtinoutFolderDirectory)
        shutil.copy2(libpath+"/ScriptsForSWATtxt/nsga2_mid.cmd", SWATtxtinoutFolderDirectory)
        shutil.copy2(libpath+"/ScriptsForSWATtxt/swat.exe", SWATtxtinoutFolderDirectory)
        shutil.copy2(libpath+"/ScriptsForSWATtxt/SWAT_ParameterEdit.py", SWATtxtinoutFolderDirectory)
        #Read ('nsga2.def') NSGA-II binary options input
        nsga2def = SWATtxtinoutFolderDirectory+"/NSGA2.IN/nsga2.def"
        nsga2pardef = SWATtxtinoutFolderDirectory+"/NSGA2.IN/nsga2_par.def"
        observed_rch = SWATtxtinoutFolderDirectory+"/NSGA2.IN/observed_rch.txt"
        f = open(nsga2def, "r")
        lines = f.readlines()
        popsize = int(lines[1].split()[1]) #Population size (an even no.)
        gener = int(lines[2].split()[1]) #the no.of generations
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
        #print "popsize,gener,pcross,optype,bits,pmutprp,seed,M,FuncOpt,FuncOptAvr,ReadMFrmOut"
        #print popsize,gener,pcross,optype,bits,pmutprp,seed,M,FuncOpt,FuncOptAvr,ReadMFrmOut
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

        #---Thinking if observations are not all rch file. I will incorporate the others later.-----------------------------------------------------------------------------------------------------------------
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
        self.gener=gener
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
        self.old_pop_ptr=nsga2funcs.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        self.new_pop_ptr=nsga2funcs.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        self.mate_pop_ptr=nsga2funcs.CreateDefaultPopulation(self.popsize,self.chrom,self.nchrom,self.nfunc)
        #/*Initialize the random no generator*/
        self.warmup_random = random_(seed); #
        #-------------------------------------------------------------------------------

    def DetermineInitialPopulation(self):
        old_pop_ptr=self.old_pop_ptr
        #Check if NSGA2.OUT exist
        if not os.path.exists(self.SWATdir+"/NSGA2.OUT"): os.makedirs(self.SWATdir+"/NSGA2.OUT")
        
        #-------------------------------- INITIAL POPULATION DETERMINATION ----------------------------------------    
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
                    old_pop_ptr["ind"][i]["xbin"][j] = LastGenPars[0][j]
            #Copy old output file
            shutil.copy2(self.SWATdir+"/NSGA2.OUT/output.out", self.SWATdir+"/NSGA2.OUT/output_previous.out")            
            #/*Function Calculaiton*/
            SWATnsga2funcs.CalculateObjectiveFunctions(old_pop_ptr,self.Outlet_Obsdata,self.FuncOpt,self.FuncOptAvr,self.parname,"Previous NSGA-II run Last Population",self.SWATdir)
            
        else:
            #Defining Latin Hypercube Sampling population
            InitialLHSpop = nsga2funcs.CreateDefaultPopulation(self.M,self.chrom,self.nchrom,self.nfunc)
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
            SWATnsga2funcs.CalculateObjectiveFunctions(InitialLHSpop,self.Outlet_Obsdata,self.FuncOpt,self.FuncOptAvr,self.parname,"LHS",self.SWATdir)
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
        nsga2funcs.reverse_decode(old_pop_ptr,self.vlen,self.lim_b)
        SWATnsga2funcs.rankcon(old_pop_ptr)
        self.old_pop_ptr=old_pop_ptr
        #-------------------------------- INITIAL POPULATION DETERMINATION ----------------------------------------


    #-------------------------------------------------------------------------------
    #/*This is the file to get the different individuals selected*/
    def Selection(self):
        old_pop_ptr=self.old_pop_ptr; pop2_ptr=self.mate_pop_ptr; warmup_random=self.warmup_random
        popsize = len(old_pop_ptr["ind"])
        chrom = len(old_pop_ptr["ind"][0]['genes'])

        r = popsize;
        s = chrom;

        indZeros = nsga2funcs.indzeros(chrom)

        k=-1;
        for n in xrange(popsize):
            k+=1
            j=0;j1 = 0;

            rnd2 = warmup_random.randomperc();
            rnd2 = popsize*rnd2;
            rnd = int(floor(rnd2));
            if(rnd == 0):rnd = popsize - k;
            if(rnd == popsize):rnd = abs(popsize-2)/2;

            # /*Select first parent randomly*/
            if rnd <= 0:j = indZeros #The population has max individual members in the c code but not here so last item there is zero
            else:j = old_pop_ptr["ind"][rnd-1];

            rnd2 = warmup_random.randomperc();
            rnd2 = popsize * rnd2;
            rnd1 = int(floor(rnd2));
            if (rnd1 == 0):rnd1 = popsize - n;
            if(rnd1 == popsize):rnd1 = abs(popsize - 4)/2;

            #/*Select second parent randomly*/
            if rnd1 <= 0:j1 = indZeros #The population has max individual members in the c code but not here so last item there is zero
            else:j1 = old_pop_ptr["ind"][rnd1-1];

            s1_ptr = j["genes"][:];
            fit_ptr1 = j["rank"];
            f1_ptr = j["cub_len"];

            s2_ptr = j1["genes"][:];
            fit_ptr2 = j1["rank"];
            f2_ptr = j1["cub_len"];

        #/*---SELECTION PROCEDURE---*/
          #/*Comparing the fitnesses*/
            if(fit_ptr1 > fit_ptr2):pop2_ptr["ind"][k]["genes"][:] = s2_ptr[:];
            elif(fit_ptr1 < fit_ptr2):pop2_ptr["ind"][k]["genes"][:]=s1_ptr[:];
            elif(f1_ptr < f2_ptr):pop2_ptr["ind"][k]["genes"][:] = s2_ptr[:];
            else:pop2_ptr["ind"][k]["genes"][:] = s1_ptr[:];
        self.old_pop_ptr=old_pop_ptr; self.mate_pop_ptr=pop2_ptr; self.warmup_random=warmup_random
        return

    #/*CROSSOVER----------------------------*/
    def Crossover(self):
        if(self.optype == 1): #/*Binary Cross-over*/
            self.ncross = nsga2funcs.crossover(self.new_pop_ptr,self.mate_pop_ptr,self.warmup_random,self.pcross,self.ncross);

        if(self.optype == 2):#/*Binary Uniform Cross-over*/
            self.ncross = nsga2funcs.unicross(self.new_pop_ptr ,self.mate_pop_ptr,self.warmup_random,self.pcross,self.ncross);
    #-------------------------------------------------------------------------------
    #/* This is the module used to formulate the mutation routine*/
    def Mutation(self):
        rand1=self.warmup_random.randomperc();
        j=0
        while j < self.popsize:
            #/*Select bit */
            i=0
            while i < self.chrom:
                rand1 = self.warmup_random.randomperc();

                #/*Check whether to do mutation or not*/
                if(rand1 <= self.pmut_b):
                    if(self.new_pop_ptr['ind'][j]['genes'][i] == 0):
                        self.new_pop_ptr['ind'][j]['genes'][i] =1;
                    else:
                        self.new_pop_ptr['ind'][j]['genes'][i]=0;
                    self.nmut+=1;
                i+=1
            j+=1
        return;
    #-------------------------------------------------------------------------------
    
    #/*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
    #/*Elitism And Sharing Implemented*/ #### Nondominated sorting, crowding distances
    def CreateMatePopulation(self,generationNo):
        nsga2funcs.CreateMatePopFromNewandOldPops(self.old_pop_ptr ,self.new_pop_ptr ,self.mate_pop_ptr,generationNo)
















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


