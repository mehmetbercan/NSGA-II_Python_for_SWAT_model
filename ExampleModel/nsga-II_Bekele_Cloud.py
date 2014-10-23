# -*- coding: cp1252 -*-
#-------------------------------------------------------------------------------
# Name:        NSGA-II
# Purpose:     Multi-Objective calibration of SWAT model
#
# Author:      Mehmet B. Ercan (mehmetbercan@gmail.com) at USC/Columbia, SC.
#
# Created:     11/14/2013
# Copyright:   (c) Mehmet B. Ercan 2013
# Licence:
#-------------------------------------------------------------------------------

import copy, datetime, sys, numpy, os, math, random, shutil
from math import floor

def MainNSGA2():
    print "\n"*5,"                                   NSGA-II\n"
    print "    Author: Mehmet B. Ercan (mehmetbercan@gmail.com) at USC/Columbia, SC USA"
    
    #Read Input file
    parname,popsize,gener,pcross,optype,seed,M,FuncOpt,FuncOptAvr,ReadMFrmOut,nchrom,vlen,chrom,lim_b,pmut_b,Outlet_Obsdata,nfunc=ReadInput()

    #Open Plot file
    if not os.path.exists(r'./NSGA2.OUT' ): os.makedirs(r'./NSGA2.OUT')
    Plot = open("./NSGA2.OUT/plot.out","w"); #("plot.out","w");
    Plot.writelines("# Feasible and Non-dominated Objective Vector\n");

    #counters for mutation and crossover
    nmut = 0;
    ncross = 0;

    #typedef struct    /*Popuation Structure*/
    def population(popsize):
        popltn={}
        popltn["maxrank"]=0        #/*Maximum rank present in the population*/
        popltn["rankno"]=numpy.zeros(2*popsize,int)#/*Individual at different ranks*/
        popltn["ind"]=[]    #/*Different Individuals*/
        for i in xrange(popsize):
            indvdl = {}
            indvdl["genes"]=numpy.zeros(chrom,int)#/*bianry chromosome*/
            indvdl["rank"]=0      #/*Rank of the individual*/
            indvdl["flag"]=0      #/*Flag for ranking*/
            indvdl["xbin"]=numpy.zeros(nchrom,float)#/*list of decoded value of the chromosome */
            indvdl["fitness"]=numpy.zeros(nfunc,float)#/*Fitness values */
            indvdl["cub_len"]=0.0 #/*crowding distance of the individual*/
            popltn["ind"].append(indvdl)
        return popltn

    #/*Defining the population Structures*/
    old_pop_ptr = population(popsize)
    new_pop_ptr = population(popsize)
    mate_pop_ptr = population(popsize)
    
    #/*Initialize the random no generator*/
    warmup_random = random_(seed); #
    
#-------------------------------- INITIAL POPULATION DETERMINATION ----------------------------------------    
    #Define the the initital population
    if ReadMFrmOut == 1:#Read Last population from output.out
        #Read outputout
        f = open("./NSGA2.OUT/output.out","r")
        lines = f.readlines()
        prmtrno = int(lines[5].split(")")[0].split("binary")[1]) #Number of parameters (binary)
        if nchrom != prmtrno: sys.exit("ERROR: parameter number is not equal to parameter numbers defined in output.out file")
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
        if popsize != len(LastGenPars): sys.exit("ERROR: poulation size is not equal to population size in output.out file")
        #Write values in old_pop_ptr
        for i in xrange(popsize):
            for j in xrange(nchrom):
                old_pop_ptr["ind"][i]["xbin"][j] = LastGenPars[0][j]
        #Open output file
        shutil.copy2("./NSGA2.OUT/output.out", "./NSGA2.OUT/output_previous.out")
        Output = open("./NSGA2.OUT/output.out","w");
        #/*Function Calculaiton*/
        func(old_pop_ptr,Outlet_Obsdata,FuncOpt,FuncOptAvr,parname,"Previous NSGA-II run Last Population")
        
    else:
        #Open output file
        Output = open("./NSGA2.OUT/output.out","w");
        #Defining Latin Hypercube Sampling population
        InitialLHSpop = population(M)
        #Latin Hypercube Samples
        LHSamples = []
        for j in xrange(len(lim_b)):
            LowBound = lim_b[j][0]
            UpBound = lim_b[j][1]
            parVals = []
            for i in xrange(M+1):
                point = i * 1. / M
                parValue = (point * (UpBound - LowBound)) + LowBound
                parVals.append(parValue)
            LHSamples.append(parVals)
        #Random selection of values from each interval of LHS
        for i in xrange(M):
            for j in xrange(nchrom):
                rndinteger = int(round(M*random.random(),0))
                InitialLHSpop["ind"][i]["xbin"][j] = LHSamples[j][rndinteger]

        #/*Function Calculaiton*/
        func(InitialLHSpop,Outlet_Obsdata,FuncOpt,FuncOptAvr,parname,"LHS")

        #Select the first population from InitialLHSpop
        n=0
        for rank in range(1,M+1):
            for i in xrange(M):
                if rank >= InitialLHSpop["ind"][i]["rank"] and InitialLHSpop["ind"][i]["flag"] != 99:
                    old_pop_ptr["ind"][n] = copy.deepcopy(InitialLHSpop["ind"][i])
                    n+=1
                    InitialLHSpop["ind"][i]["flag"] = 99
                    if n == popsize:
                        break
            else:continue
            break
    reverse_decode(old_pop_ptr,vlen,lim_b)
    rankcon(old_pop_ptr)
#-------------------------------- INITIAL POPULATION DETERMINATION ----------------------------------------
    
    #/********************************************************************/
    #/*----------------------GENERATION STARTS HERE----------------------*/
    i=0
    while i < gener:
        #/*--------SELECT----------------*/
        nselect(old_pop_ptr,mate_pop_ptr,warmup_random);

        #/*CROSSOVER----------------------------*/
        if(optype == 1): #/*Binary Cross-over*/
            ncross = crossover(new_pop_ptr,mate_pop_ptr,warmup_random,pcross,ncross);

        if(optype == 2):#/*Binary Uniform Cross-over*/
            ncross = unicross(new_pop_ptr ,mate_pop_ptr,warmup_random,pcross,ncross);
        #/*------MUTATION-------------------*/
        nmut = mutate(new_pop_ptr,warmup_random,pmut_b,nmut); #/*Binary Mutation */

        #/*-------DECODING----------*/
        decode(new_pop_ptr,vlen,lim_b); #/*Decoding for binary strings*/

        #/*----------FUNCTION EVALUATION-----------*/
        #-----------------------------------------------------
        func(new_pop_ptr,Outlet_Obsdata,FuncOpt,FuncOptAvr,parname,i+1); #------ MAIN PART FOR SWAT RUN ----
        #-----------------------------------------------------

        #/*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
        #/*Elitism And Sharing Implemented*/
        keepalive(old_pop_ptr ,new_pop_ptr ,mate_pop_ptr,i+1); # Elitism, crowding distances
        decode(mate_pop_ptr,vlen,lim_b);

        #/*------------------REPORT PRINTING--------------------------------*/
        report(i,old_pop_ptr,mate_pop_ptr,Output,Plot,gener);

        #/*==================================================================*/
        new_pop_ptr = copy.deepcopy(mate_pop_ptr);
        old_pop_ptr = copy.deepcopy(new_pop_ptr);

        i+=1

    #/*                   Generation Loop Ends                                */
    #/************************************************************************/

    Output.writelines("NO. OF CROSSOVER = %d\n"%ncross);
    Output.writelines("NO. OF MUTATION = %d\n"%nmut);
    Output.writelines("------------------------------------------------------------\n");
    Output.writelines("---------------------------------Thanks---------------------\n");
    Output.writelines("-------------------------------------------------------------\n");
    print ("NOW YOU CAN LOOK IN THE FILE output.out\n");

    #/*Closing the files*/
    Output.close()
    Plot.close()






















######################################################################################################
# Author:Mehmet B. Ercan (University of South Carolina/CEE/Water Resources, mehmetbercan@gmail.com)
#Includes all the functions for binary nsga2 calibration
######################################################################################################
#-------------------------------------------------------------------------------
def ReadInput():
    #Read ('nsga2.def') NSGA-II binary options input
    f = open("./NSGA2.IN/nsga2.def", "r")
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
    f = open("./NSGA2.IN/nsga2_par.def", "r")
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
    f = open("./NSGA2.IN/observed_rch.txt", "r")
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
    return parname,popsize,gener,pcross,optype,seed,M,FuncOpt,FuncOptAvr,ReadMFrmOut,nchrom,vlen,chrom,lim_b,pmut_b,Outlet_Obsdata,nfunc
#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#This is the program used to evaluate the value of the function
def func(population,Outlet_Obsdata,FuncOpt,FuncOptAvr,parname, generation):
    #/*Initializing the max rank to zero*/
    population["maxrank"]=0

    popsize = len(population["ind"])
    nfunc = len(population["ind"][0]["fitness"])

    i=0;
    ParameterValues=[]
    while i < popsize:
        ParameterValues.append(population["ind"][i]["xbin"]) #/* problem variables */
        i+=1

    #&&&& without pralel SWAT run &&&&&
    outlets = Outlet_Obsdata.keys()
    outlets.sort()
    for i in xrange(len(ParameterValues)): #population loop
        print "\n"*5,"-"*45,"\nGeneration: ", generation, "  Simulation: ", i+1, "\n", "-"*45
        #Print parameter set in model.in file
        modelinf = open(".\model.in","w")
        writeline =''
        for j in xrange(len(ParameterValues[i])): #parameter loop
            writeline += parname[j]+'\t'+str(ParameterValues[i][j])+'\n'       
        modelinf.writelines(writeline)
        modelinf.close()
        
        #Run command file (SWATedit, SWAT and extract exe files)
        os.system(r'.\nsga2_mid.cmd')

        #Read 'model.out' file
        modelrchf = open('./model.out','r')
        lines = modelrchf.readlines()
        Outlet_Modeldata = {}
        k=0; Modeldata=[]
        for outlet in outlets:
            nofdatapoints = len(Outlet_Obsdata[outlet])
            for j in range(k,k+nofdatapoints):
                Modeldata.append(float(lines[j].split()[1]))
            Outlet_Modeldata[outlet] = Modeldata
            k = j+1
            Modeldata=[]
            
        #Calculate Objective functions for each site (gage)
        objectivefuncs = []
        for outlet in outlets:
            outflowSWAT = Outlet_Modeldata[outlet]
            outflowUSGS = Outlet_Obsdata[outlet]
            #Define x and y for model efficiency coefficients
            x = outflowSWAT #Simulated parameters
            y = outflowUSGS #Measured parameters

            if FuncOpt == 1:
                E = Nash_Sutcliffe(x,y) #Nash-Sutcliffe model efficiency coefficient
                E0best = 1 - E #0 is the best and +infinity is the worst
                objectivefuncs.append(E0best)
            if FuncOpt == 2:
                R = numpy.corrcoef(x, y)[0,1] 
                R2 = math.pow(R,2) #Corelation coefficient
                R20best = 1 - R2 #0 is the best and + 1 is the worst
                objectivefuncs.append(R20best)
            if FuncOpt == 3:
                E = Nash_Sutcliffe(x,y) #Nash-Sutcliffe model efficiency coefficient
                E0best = 1 - E #0 is the best and +infinity is the worst
                LE = Log_Nash_Sutcliffe(x,y) #Log Nash-Sutcliffe model efficiency coefficient
                LE0best = 1 - LE #0 is the best and +infinity is the worst
                objectivefuncs.append(E0best)
                objectivefuncs.append(LE0best)
            if FuncOpt == 5:
                E = Nash_Sutcliffe(x,y) #Nash-Sutcliffe model efficiency coefficient
                E0best = 1 - E #0 is the best and +infinity is the worst
                PB = PercentBias(x,y) #Log Nash-Sutcliffe model efficiency coefficient
                PB0best = abs(PB/100.0) #0 is the best and +infinity is the worst 
                objectivefuncs.append(E0best)
                objectivefuncs.append(PB0best)
            if FuncOpt == 4:
                E = Nash_Sutcliffe(x,y) #Nash-Sutcliffe model efficiency coefficient
                E0best = 1 - E #0 is the best and +infinity is the worst
                LE = Log_Nash_Sutcliffe(x,y) #Log Nash-Sutcliffe model efficiency coefficient
                LE0best = 1 - LE #0 is the best and +infinity is the worst
                R = numpy.corrcoef(x, y)[0,1] 
                R2 = math.pow(R,2) #Corelation coefficient
                R20best = 1 - R2 #0 is the best and + 1 is the worst
                objectivefuncs.append(E0best)
                objectivefuncs.append(LE0best)
                objectivefuncs.append(R20best)  
        #Average objective functions
        nobjfunc_=1
        if FuncOpt==3 or FuncOpt==5:
            nobjfunc_=2
        if FuncOpt==4:
            nobjfunc_=3
        nobjsite_=len(outlets)
        newobjectivefuncs=[]
        if FuncOptAvr==1: #Average Objective sites
            for i in range(0,nobjfunc_):
                objfuncav=0
                for j in range(0,nobjsite_):
                    objfuncav+=objectivefuncs[(i+j*nobjfunc_)]/nobjsite_
                newobjectivefuncs.append(objfuncav)  
        elif FuncOptAvr==2: #Average Objective Functions
            for i in range(0,nobjsite_):
                objsiteav=0
                for j in range(0,nobjfunc_):
                    objsiteav+=objectivefuncs[(i*nobjfunc_+j)]/nobjfunc_
                newobjectivefuncs.append(objsiteav)             
        else: #Do not average
            newobjectivefuncs = objectivefuncs
        #Add objective functions to population
        population["ind"][i]['fitness'] = newobjectivefuncs
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    #/*---------------------------* RANKING *------------------------------*/
    rankcon(population);
    return;

##Nash–Sutcliffe model efficiency coefficient
def Nash_Sutcliffe(SimulatedStreamFlow, ObservedStreamFlow):
    '''(SimulatedStreamFlow, ObservedStreamFlow)''' 
    x=SimulatedStreamFlow
    y=ObservedStreamFlow
    A=0.0 #dominator
    B=0.0 #deminator
    tot = 0.0
    for i in range(0, len(y)):
        tot = tot + y[i]
    average = tot / len(y)
    for i in range(0, len(y)):
        A = A + math.pow((y[i] - x[i]), 2)
        B = B + math.pow((y[i] - average), 2)
    E = 1 - (A/B) # Nash-Sutcliffe model eficiency coefficient
    return E

##Logaritmic Nash–Sutcliffe model efficiency coefficient
def Log_Nash_Sutcliffe(SimulatedStreamFlow, ObservedStreamFlow):
    '''(SimulatedStreamFlow, ObservedStreamFlow)''' 
    x=SimulatedStreamFlow
    y=ObservedStreamFlow
    A=0.0 #dominator
    B=0.0 #deminator
    tot = 0.0
    for i in range(0, len(y)):
        tot = tot + y[i]
    average = tot / len(y)
    for i in range(0, len(y)):
        X = x[i]
        Y = y[i]
        if X == 0:
            X = 1e-3
        if Y == 0:
            Y = 1e-3
        A = A + math.pow((math.log(Y) - math.log(X)), 2) #log = ln or log_e
        B = B + math.pow((math.log(Y) - math.log(average)), 2)
    E = 1 - (A/B) # Nash-Sutcliffe model eficiency coefficient
    return E

#Define definition for Percent Bias model efficiency coefficient---used up in the class
def PercentBias(SimulatedStreamFlow, ObservedStreamFlow):
    '''(SimulatedStreamFlow, ObservedStreamFlow)'''  
    x=SimulatedStreamFlow
    y=ObservedStreamFlow
    A=0.0 #dominator
    B=0.0 #deminator
    for i in range(0, len(y)):
        A = A + (y[i] - x[i])
        B = B + y[i]
    PB = 100.0*(A/B) # Percent Bias model eficiency coefficient
    return PB
#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#/*This program subroutine is used to print the report*/
def report(t,pop1_ptr,pop2_ptr,rep_ptr,lastit,gener):
    popsize = len(pop1_ptr["ind"])
    nchrom = len(pop1_ptr["ind"][0]["xbin"])
    nfunc = len(pop1_ptr["ind"][0]["fitness"])

    rep_ptr.writelines("\n\n---------------------------------------------------\n");
    rep_ptr.writelines("Generation No.     ->%d\n"%(t+1));
    rep_ptr.writelines("------------------------------------------------------\n");

    rep_ptr.writelines(" variables (binary %d)  fitness (%d)  rank cublen || variables  fitness rank cublen\n"%(nchrom,nfunc));

    i=0
    while i < popsize:
        rep_ptr.writelines("\n------------------------------------------------\n");

        ptr1_b = pop1_ptr['ind'][i]['xbin'];
        ptr2_b = pop2_ptr['ind'][i]['xbin'];

        fptr = pop1_ptr['ind'][i]['fitness'];
        fptr1 = pop2_ptr['ind'][i]['fitness'];

        rptr = pop1_ptr['ind'][i]['rank'];
        rptr1 = pop2_ptr['ind'][i]['rank'];

        j=0
        while j < nchrom:
            rep_ptr.writelines("%f "%ptr1_b[j]);
            j+=1

        if (t == gener-1):
            j=0
            while j < nfunc:
                if (rptr1 == 1):
                    lastit.writelines("%f\t"%(fptr1[j]));
                j+=1
            if (rptr1 == 1):
                lastit.writelines("\n");

        fptr =  pop1_ptr['ind'][i]['fitness'];
        fptr1 = pop2_ptr['ind'][i]['fitness'];

        j=0
        while j < nfunc:
            rep_ptr.writelines("  %.4f"%fptr[j]);
            j+=1

        rep_ptr.writelines(" %d "%rptr);

        rep_ptr.writelines("%f "%pop1_ptr['ind'][i]['cub_len']);
        rep_ptr.writelines("|**|");

        j=0
        while j < nchrom:
            rep_ptr.writelines("%f "%ptr2_b[j]);
            j+=1
        j=0
        while j < nfunc:
            rep_ptr.writelines("  %f"%fptr1[j]);
            j+=1
        rep_ptr.writelines(" %d "%rptr1);

        rep_ptr.writelines(" %f "%pop2_ptr['ind'][i]['cub_len']);

        i+=1

    rep_ptr.writelines("\n--------------------------------------------------\n\n");
    rep_ptr.writelines("-------------------------------------------------------\n");
    return;
#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#/*This is the file for formulating the crossover process*/
def crossover(new_pop_ptr,mate_pop_ptr,warmup_random,pcross,ncross):
    popsize = len(new_pop_ptr["ind"])
    chrom = len(new_pop_ptr["ind"][0]["genes"])

    rnd=warmup_random.randomperc();

    i=0;y=0;n=0;
    while i < popsize/2:
        chld1=new_pop_ptr['ind'][n]['genes'];
        n = n+1;

        chld2=new_pop_ptr['ind'][n]['genes'];
        n = n+1;

        par1 = mate_pop_ptr['ind'][y]['genes'];
        y = y+1;

        par2 = mate_pop_ptr['ind'][y]['genes'];
        y = y+1;

        rnd = warmup_random.randomperc();
        iptr=0
        if (rnd < pcross):
            ncross+=1;
            rnd = warmup_random.randomperc();
            c = math.floor(rnd*(chrom+10));
            mating_site = c;
            if(mating_site >= chrom):
                mating_site = mating_site/2.0;

            k=0
            while k < chrom:
                if(k > mating_site-1):
                    chld1[iptr] = par2[iptr];
                    chld2[iptr] = par1[iptr];
                else:
                    chld1[iptr] = par1[iptr];
                    chld2[iptr] = par2[iptr];
                k+=1
                iptr+=1
        else:
            k=0
            while k < chrom:
                chld1[iptr] = par1[iptr];
                chld2[iptr] = par2[iptr];
                k+=1
                iptr+=1
        i+=1
    return ncross;

#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#/*This is the program to decode the chromosome to get real values*/
def decode(pop_ptr,vlen,lim_b):
    popsize = len(pop_ptr["ind"])
    nchrom = len(pop_ptr["ind"][0]["xbin"])

    i=0
    while i < popsize:
        m=0
        coef = numpy.zeros(nchrom,float)
        gptr = 0
        while m < nchrom:
            #/*finding out the co-efficient 2 to the power of
            #(l-1) where l is the no of bits assigned to this variable
            #For More Info Study DEB's Book*/
            sum_ = 0;
            k=0
            while k < vlen[m]:
                b = pop_ptr['ind'][i]['genes'][gptr];
                d = vlen[m] - k - 1;
                c = pow(2,d);
                sum_ =sum_ + c * b;
                k+=1
                gptr+=1

            x = vlen[m];
            coef[m] = pow(2,x) - 1;
            pop_ptr['ind'][i]['xbin'][m] =lim_b[m][0] + (sum_/coef[m])*(lim_b[m][1]-lim_b[m][0]);
            m+=1
        i+=1
    return;

def reverse_decode(pop_ptr,vlen,lim_b):#Calculates binary (genes) values from the xbin values
    popsize = len(pop_ptr["ind"])
    nchrom = len(pop_ptr["ind"][0]["xbin"])

    i=0
    while i < popsize:
        m=0
        while m < nchrom:
            sum_ = int((pop_ptr['ind'][i]['xbin'][m] - lim_b[m][0])*(pow(2,vlen[m]) - 1)/(lim_b[m][1]-lim_b[m][0]))
            k=vlen[m]-1
            while k > -1:
                if sum_%2 == 1:
                    pop_ptr['ind'][i]['genes'][(m*vlen[0] + k)] = 1
                    sum_ = sum_ - 1
                    sum_ = sum_ / 2
                else:
                    pop_ptr['ind'][i]['genes'][(m*vlen[0] + k)] = 0
                    sum_ = sum_ / 2
                k-=1
            m+=1
        i+=1
    return;
#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#/*This is a routine to keep the fronts alive (caring the end problem)*/
def globpop(popsize,chrom,nchrom,nfunc):
    popltn = {}
    popltn["maxrank"] = 0   #/*Max rank of the global population*/
    popltn["rankar"] = numpy.zeros((2*popsize,2*popsize), int) #/*record of array of individual numbers at a particular rank */
    popltn["rankno"] = numpy.zeros(2*popsize, int);           #/*record of no. of individuals at a particular rank*/
    popltn["genes"] = numpy.zeros((2*popsize,chrom), int)
    popltn["rank"] = numpy.zeros(2*popsize, int)            #/*rank of different individuals*/
    popltn["flag"] = numpy.zeros(2*popsize, int)               #/*Setting the flag */
    popltn["fitness"] = numpy.zeros((2*popsize,nfunc), float) #/*Fitness function values for the different	                          # individuals*/
    popltn["cub_len"] = numpy.zeros(2*popsize, float);              #/*Dummyfitness*/
    popltn["xbin"] = numpy.zeros((2*popsize,nchrom) , float);    #/* binray-coded variables */
    return popltn

def keepalive(pop1_ptr,pop2_ptr,pop3_ptr,gen):
    popsize = len(pop1_ptr["ind"])
    chrom = len(pop1_ptr["ind"][0]["genes"])
    nchrom = len(pop1_ptr["ind"][0]["xbin"])
    nfunc = len(pop1_ptr["ind"][0]["fitness"])

    globalpop = globpop(popsize,chrom,nchrom,nfunc)
    global_pop_ptr = globpop(popsize,chrom,nchrom,nfunc)

    fpara1 = numpy.zeros((2*popsize,2),float);
    Lastrank = 0

    #/*Forming the global mating pool*/
    i=0
    while i<popsize:
        if(nchrom > 0):
            #/*Binary Coded GA genes are copied*/
            k=0
            while k<chrom:
                globalpop["genes"][i][k]=pop1_ptr["ind"][i]["genes"][k];
                globalpop["genes"][i+popsize][k] = pop2_ptr["ind"][i]["genes"][k];
                k+=1
            k=0
            while k<nchrom:
                globalpop["xbin"][i][k] = pop1_ptr["ind"][i]["xbin"][k];
                globalpop["xbin"][i+popsize][k] = pop2_ptr["ind"][i]["xbin"][k];
                k+=1

        #/*Fitness is copied to the global pool */
        l=0
        while l<nfunc:
            globalpop["fitness"][i][l] = pop1_ptr["ind"][i]['fitness'][l];
            globalpop["fitness"][i+popsize][l] = pop2_ptr["ind"][i]["fitness"][l];
            l+=1

        #/*Initialising the dummyfitness to zero */
        globalpop["cub_len"][i] = 0;
        globalpop["cub_len"][i+popsize] = 0;
        i+=1

    global_pop_ptr = globalpop;

    #/*Finding the global ranks */
    grank(gen,popsize,global_pop_ptr,nfunc,globalpop);

    m = globalpop['maxrank'];

    #/* Sharing the fitness to get the dummy fitness */
    i=0
    while i<m:
        gshare(i+1,popsize,nfunc,globalpop,fpara1);
        i+=1

    poolf = popsize;
    pool = 0;
    
    #/*Initializing the flags of population to zero */
    i=0
    while i<(2*popsize):
        globalpop['flag'][i] = 0;
        i+=1
    #// decide which all solutions belong to the pop3
    rec = 0;
    st = 0;
    i=0
    while i<m:
        #/*    Elitism Applied Here     */
        st = pool;
        pool += globalpop["rankno"][i];

        if(pool <= popsize):
            k=0
            while k<(2*popsize):
                if(globalpop['rank'][k] == i+1):
                    globalpop['flag'][k] = 1;
                k+=1
            pop3_ptr['rankno'][i] = globalpop['rankno'][i];
        else:
            sel = popsize - st;
            Lastrank = i+1;
            pop3_ptr['rankno'][i] = sel;
            gsort(i+1,sel,popsize,globalpop);
            break;
        i+=1

    k = 0;i=0;
    while(i < 2*popsize and k < popsize):
        if(nchrom > 0):
            if(globalpop['flag'][i] == 1):
                gene1_ptr = globalpop['genes'][i];
                xbin1_ptr = globalpop['xbin'][i];
                gene2_ptr = 0
                gene2_ptr = pop3_ptr['ind'][k]['genes'];
                xbin2_ptr = pop3_ptr['ind'][k]['xbin'];

                j=0
                while j<chrom:
                    gene2_ptr[j] = gene1_ptr[j];
                    j+=1
                j=0
                while j<nchrom:
                    xbin2_ptr[j] = xbin1_ptr[j];
                    j+=1

        if(globalpop['flag'][i] == 1):
            j=0
            while j<(nfunc):
                pop3_ptr['ind'][k]['fitness'][j] = globalpop['fitness'][i][j];
                j+=1
            pop3_ptr['ind'][k]['cub_len'] = globalpop['cub_len'][i];

            jj=0

            pop3_ptr['ind'][k]['rank'] = globalpop['rank'][i];
            k+=1;  #// increment the pop3 counter
        i+=1

    pop3_ptr['maxrank'] = Lastrank;

    return;

def grank(gen,popsize,global_pop_ptr,nfunc,globalpop):
    gr = open("./NSGA2.Out/g_rank_record.out","a");
    gr.writelines("Genration no. = %d\n"%gen);
    #/*----------------------------* RANKING *---------------------------------*/
    rnk = 0;
    nondom = 0;
    popsize1 = 2*popsize;
    gflg = numpy.zeros(2*popsize,int)

    i=0
    while i<(popsize1):
        gflg[i] = 2;
        i+=1

    k=0
    while k<(popsize1):
        q =  0;
        j=0
        while j<(popsize1):
            if (gflg[j] != 1): break;
            j+=1
        if(j == (popsize1)): break;
        rnk = rnk +1;
        j=0
        while j<(popsize1):
            if(gflg[j] == 0): gflg[j] = 2;
            j+=1
        i=0
        while i<(popsize1):
            if(gflg[i] != 1 and gflg[i] != 0):
                ptr1 = global_pop_ptr["fitness"][i];
                j=0
                while j<(popsize1):
                    if( i!= j):
                        if(gflg[j] != 1):
                            ptr2 = global_pop_ptr["fitness"][j];
                            val = indcmp1(ptr1,ptr2,nfunc);
                            if( val == 2):
                                gflg[i] = 0;#/* individual 1 is dominated */
                                break;
                            if(val == 1):
                                gflg[j] = 0;#/* individual 2 is dominated */
                            if(val == 3):
                                nondom+=1;#/* individual 1 & 2 are non dominated */
                                if(gflg[j] != 0):gflg[j] = 3;
                    j+=1
                if( j == (popsize1)):
                    global_pop_ptr["rank"][i] = rnk;
                    gflg[i] = 1;
                    global_pop_ptr["rankar"][rnk-1][q] =  i;
                    q+=1;
            i+=1
        global_pop_ptr["rankno"][rnk-1] = q;
        k+=1
    global_pop_ptr["maxrank"] = rnk;
    gr.writelines("   RANK     No Of Individuals\n");
    i=0
    while i<(rnk):
        gr.writelines("\t%d\t%d\n"%(i+1,globalpop["rankno"][i]));
        i+=1

    gr.close();
    return;

def indcmp1(ptr1,ptr2,nfunc):
    fit1=numpy.zeros(nfunc,float)
    fit2=numpy.zeros(nfunc,float)
    i=0
    value = 3 # Mehmet: I added as in some cases value is not defined
    while i<(nfunc):
        fit1[i] = ptr1[i];
        fit2[i] = ptr2[i];
        i+=1
    m = 0;n=0;
    while(m < nfunc and fit1[m] <= fit2[m]):
        if((fit2[m] - fit1[m]) < 1e-7): n+=1;
        m+=1;
    if(m == nfunc):
        if(n == nfunc): value = 3;
        else: value = 1;                    #/*value = 1 for dominating*/
    else:
        m = 0;n = 0;
        while(m < nfunc and fit1[m] >= fit2[m]):
            if((fit1[m] - fit2[m]) < 1e-7): n+=1;
            m+=1;
        if(m == nfunc):
            if(n != nfunc):
                value = 2;                       #/*value =  2 for dominated */
        else: value = 3;                   #/*value = 3 for incomparable*/
    return value;

#/* This is the file used to sort the dummyfitness arrays */
def gsort(rnk,sel,popsize,globalpop):
    array = numpy.zeros((2*popsize,2),float);

    q = globalpop['rankno'][rnk-1];

    i=0
    while i<q:
        array[i][0] = globalpop['rankar'][rnk-1][i];
        a = globalpop['rankar'][rnk-1][i];
        array[i][1] = globalpop['cub_len'][a];
        i+=1
    i=0
    while i<q:
        j=i+1
        while j<q:
            if(array[i][1] < array[j][1]):
                temp = array[i][1];
                temp1 = array[i][0];
                array[i][1] = array[j][1];
                array[i][0] = array[j][0];

                array[j][1] = temp;
                array[j][0] = temp1;
            j+=1
        i+=1

    i=0
    while i<(sel):
        a = array[i][0];
        globalpop['flag'][a] = 1;
        i+=1
    return;

def gshare(rnk,popsize,nfunc,globalpop,fpara1):
    length = numpy.zeros((2*popsize,2),float)

    m1 = globalpop['rankno'][rnk-1];

    j=0
    while j<(nfunc):
        i=0
        while i<(m1):
            fpara1[i][0] = 0;
            fpara1[i][1] = 0;
            i+=1
        i=0
        while i<(m1):
            a = globalpop['rankar'][rnk-1][i];
            fpara1[i][0] = float(a) ;
            fpara1[i][1] = globalpop['fitness'][a][j];
            i+=1

        sort(m1,fpara1); #/*Sort the arrays in ascending order of the fitness*/

        max = fpara1[m1-1][1];
        min = fpara1[0][1];  #// Added 18.08.2003
        Diff = max-min;      #// Added 18.08.2003 and 5 subsequent lines
        if (Diff < 0.0):
            print("Something wrong in keepaliven.h (gshare)\n");
            exit(1);
        i=0
        while i<(m1):
            if(i == 0 or i == (m1-1)):
                length[i][0] = fpara1[i][0];
                length[i][1] = 100*max;
            else:
                length[i][0] = fpara1[i][0];
                if Diff == 0.0: #Mehmet: Added in case Diff=0 error
                    Diff = 1e-10
                length[i][1] = abs(fpara1[i+1][1]- fpara1[i-1][1])/Diff; #// crowding distances are normalized 18.08.2003
            i+=1
        i=0
        while i<(m1):
            a = length[i][0];
            globalpop['cub_len'][a] += length[i][1]
            i+=1;
        j+=1

    return;

def sort(m1,fpara1):
    k1=0
    while k1<(m1-1):
        i1=k1+1
        while i1<m1:
            if(fpara1[k1][1] > fpara1[i1][1]):
                temp = fpara1[k1][1];
                temp1 = fpara1[k1][0];
                fpara1[k1][1] = fpara1[i1][1];
                fpara1[k1][0] = fpara1[i1][0];
                fpara1[i1][1] = temp;
                fpara1[i1][0] = temp1;
            i1+=1
        k1+=1
    return;
#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#/* This is the module used to formulate the mutation routine*/
def mutate(new_pop_ptr,warmup_random,pmut_b,nmut):
    popsize = len(new_pop_ptr["ind"])
    chrom = len(new_pop_ptr["ind"][0]["genes"])

    rand1=warmup_random.randomperc();

    j=0
    while j < popsize:
        #/*Select bit */
        i=0
        while i < chrom:
            rand1 = warmup_random.randomperc();

            #/*Check whether to do mutation or not*/
            if(rand1 <= pmut_b):
                if(new_pop_ptr['ind'][j]['genes'][i] == 0):
                    new_pop_ptr['ind'][j]['genes'][i] =1;
                else:
                    new_pop_ptr['ind'][j]['genes'][i]=0;
                nmut+=1;
            i+=1
        j+=1
    return nmut;
#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#This functions will be used when ncons != 0 (number of Constraints is not zero).
#/*This also demarkates the different Pareto Fronts*/
def rankcon(population): #uses "indcmp3(ptr1,ptr2)" function
    #/*---* RANKING *---*/
    #/*Initializing the ranks to zero*/
    rnk = 0 ; #/*rank*/

    nondom = 0 ;#/*no of non dominated members*/
    maxrank1 = 0;#/*Max rank of the population*/

    #/*Initializing all the flags to 2*/
    popsize = len(population["ind"])
    j=0;
    while j < popsize:
        population["ind"][j]["flag"] = 2;
        j+=1

    k=0;
    while k < popsize:
        q = 0;
        j=0;
        while j < popsize:
            if (population["ind"][j]["flag"] != 1):break;
            j+=1
            #/*Break if all the individuals are assigned a rank*/
        if(j == popsize): break;

        rnk = rnk + 1;

        j=0;
        while j < popsize:
            if(population["ind"][j]["flag"] == 0): population["ind"][j]["flag"] = 2;
            j+=1
            #/*Set the flag of dominated individuals to 2*/

        i=0;
        while i < popsize:
            if(population["ind"][i]["flag"] != 1 and population["ind"][i]["flag"] != 0):
                j=0;
                while j < popsize:
                    #/*Select the other individual which has not got a rank*/
                    if( i!= j):
                        if(population["ind"][j]["flag"] != 1):
                            #/*Compare the two individuals for
                            #fitness*/
                            val = indcmp3(population,i,j); #/*value obtained after comparing two individuals*/

                            #/*VAL = 2 for dominated individual
                            #which rank to be given*/

                            #/*VAL = 1 for dominating individual
                            #which rank to be given*/

                            #/*VAL = 3 for non comparable
                            #individuals*/
                            if( val == 2):
                                population["ind"][i]["flag"] = 0;
                                #/* individual 1 is dominated */
                                break;
                            if(val == 1):
                                population["ind"][j]["flag"] = 0;
                                #/* individual 2 is dominated */
                            if(val == 3):
                                nondom+=1;
                                #/* individual 1 & 2 are
                                #non dominated */
                                if(population["ind"][j]["flag"] != 0):
                                    population["ind"][j]["flag"] = 3;

                    j+=1
                #/*loop over j ends*/
                if( j == popsize):
                    #/*Assign the rank and set the flag*/
                    population["ind"][i]["rank"] = rnk;
                    population["ind"][i]["flag"] = 1;
                    q+=1
            i+=1
            #/*Loop over flag check ends*/
        k+=1
        #/*Loop over i ends */
        population["rankno"][rnk-1] = q ;
    maxrank1 = rnk;

    #/*     Find Max Rank of the population    */
    i=0;
    while i < popsize:
        rnk = population["ind"][i]["rank"];
        if(rnk > maxrank1):maxrank1 = rnk;
        i+=1

    population["maxrank"] = maxrank1;

    return

def indcmp3(population,ii,jj):
    #/*Routine Comparing the two individuals*/

    nfunc = len(population["ind"][0]["fitness"])
    fit1=numpy.zeros(nfunc,float)
    fit2=numpy.zeros(nfunc,float)
    i=0;
    while i < nfunc:
        fit1[i] = population["ind"][ii]["fitness"][i]#*ptr1++; population["ind_ptr"] = population["ind"][i];ptr1 = population["ind_ptr"]["fitness"][0];
        fit2[i] = population["ind"][jj]["fitness"][i]#*ptr2++;population["ind_ptr"] = population["ind"][j];ptr2 = population["ind_ptr"]["fitness"][0];
        i+=1

    m = 0;
    n = 0;
    while(m < nfunc and fit1[m] <= fit2[m]):
        if(fit1[m]== fit2[m]): n+=1;
        m+=1;

    value = 0
    if(m == nfunc):
        if(n == nfunc): value = 3;
        else: value = 1;             #/*value = 1 for dominationg*/
    else:
        m = 0;
        n = 0;
        while(m < nfunc and fit1[m] >= fit2[m]):
            if(fit1[m]== fit2[m]): n+=1;
            m+=1;
        if(m == nfunc):
            if(n != nfunc):value = 2;   #/*value =  2 for dominated */
            else: value =3;
        else: value = 3;   #/*value = 3 for incomparable*/

    return value;
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























#-------------------------------------------------------------------------------
#/*This is the file to get the different individuals selected*/
def nselect(old_pop_ptr,pop2_ptr,warmup_random):
    popsize = len(old_pop_ptr["ind"])
    chrom = len(old_pop_ptr["ind"][0]['genes'])

    r = popsize;
    s = chrom;

    indZeros = indzeros(chrom)

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
    return

def indzeros(chrom):
    indvdl = {}
    indvdl["genes"]=numpy.zeros(chrom,int)#/*bianry chromosome*/
    indvdl["rank"]=0      #/*Rank of the individual*/
    indvdl["cub_len"]=0.0 #/*crowding distance of the individual*/
    return indvdl
#-------------------------------------------------------------------------------























#-------------------------------------------------------------------------------
#/* This is the header file to do the uniform crossover */
def unicross(new_pop_ptr, mate_pop_ptr,warmup_random,pcross,ncross):
    popsize = len(new_pop_ptr["ind"])
    chrom = len(new_pop_ptr["ind"][0]["genes"])

    i=0;y=0;n=0;
    while i < popsize/2:
        j=0;
        while j < chrom:
            rnd = warmup_random.randomperc();

            #/*Checking whether to do cross-over or not*/
            if(rnd <= pcross):
                ncross+=1;
                new_pop_ptr['ind'][y]['genes'][j] = mate_pop_ptr['ind'][n+1]['genes'][j]
            else:
                new_pop_ptr['ind'][y]['genes'][j] = mate_pop_ptr['ind'][n]['genes'][j]

            new_pop_ptr['ind'][y+1]['genes'][j] = mate_pop_ptr['ind'][n+1]['genes'][j]

            j+=1
        y = y+2;
        n = n+2;
        i+=1

    return ncross;
#-------------------------------------------------------------------------------



if __name__ == '__main__':
    MainNSGA2()




