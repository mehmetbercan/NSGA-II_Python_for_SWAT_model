#Runs SWAT model with all solutions and calculate objective functions
import numpy, os, math, sys
from nsga2lib import nsga2utilities

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
##Nash-Sutcliffe model efficiency coefficient
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

##Logaritmic Nash-Sutcliffe model efficiency coefficient
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

def CalculateObjectiveFunctions(population,Outlet_Obsdata,FuncOpt,FuncOptAvr,parname, generation,SWATdir):
    nsga2utilities.round_parameters(population)
    os.chdir(SWATdir)
    #/*Initializing the max rank to zero*/
    population["maxrank"]=0

    popsize = len(population["ind"])
    nchrom = len(population["ind"][0]["xbin"])

    ParameterValues=[]
    for i in range(popsize):
        ParameterValues.append(population["ind"][i]["xbin"]) #/* problem variables */

    #&&&& without pralel SWAT run &&&&&
    outlets = list(Outlet_Obsdata.keys())
    outlets.sort()
    for i in range(popsize): #population loop
        print ("\n"*5,"-"*45,"\nGeneration: ", generation, "  Simulation: ", i+1, "\n", "-"*45)
        #Print parameter set in model.in file
        modelinf = open(os.path.join(os.getcwd(),"model.in"),"w")
        writeline =''
        for j in range(nchrom): #parameter loop
            writeline += parname[j]+'\t'+str(ParameterValues[i][j])+'\n'       
        modelinf.writelines(writeline)
        modelinf.close()
        
        #Run command file (SWATedit, SWAT and extract exe files)
        if "win" in sys.platform.lower():
            os.system(SWATdir+'/nsga2_mid.cmd')
        elif "linux" in sys.platform.lower():
            os.system(SWATdir+'/nsga2_mid.sh')
        else:
            os.system(SWATdir+'/nsga2_mid.cmd')

        #Read 'model.out' file
        modelrchf = open(os.path.join(SWATdir,"model.out"),'r')
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
            for k in range(0,nobjfunc_):
                objfuncav=0
                for j in range(0,nobjsite_):
                    objfuncav+=objectivefuncs[(k+j*nobjfunc_)]/nobjsite_
                newobjectivefuncs.append(objfuncav)  
        elif FuncOptAvr==2: #Average Objective Functions
            for k in range(0,nobjsite_):
                objsiteav=0
                for j in range(0,nobjfunc_):
                    objsiteav+=objectivefuncs[(k*nobjfunc_+j)]/nobjfunc_
                newobjectivefuncs.append(objsiteav)             
        else: #Do not average
            newobjectivefuncs = objectivefuncs
        #Add objective functions to population
        population["ind"][i]['fitness'] = newobjectivefuncs
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    #/*---------------------------* RANKING *------------------------------*/
    nsga2utilities.round_fitness(population)
    rankcon(population);
    return;
#-------------------------------------------------------------------------------
















def CalculateObjectiveFunctionsinParallel(population,Outlet_Obsdata,FuncOpt,FuncOptAvr,parname, generation,SWATdir, cpu_number):
    nsga2utilities.round_parameters(population)
    os.chdir(SWATdir)
    #/*Initializing the max rank to zero*/
    population["maxrank"]=0

    popsize = len(population["ind"])
    nchrom = len(population["ind"][0]["xbin"])

    ParameterValues=[]
    for i in range(popsize):
        ParameterValues.append(population["ind"][i]["xbin"]) #/* problem variables */

    #&&&& without pralel SWAT run &&&&&
    outlets = list(Outlet_Obsdata.keys())
    outlets.sort()
    for i in range(popsize): #population loop
        print ("\n"*5,"-"*45,"\nGeneration: ", generation, "  Simulation: ", i+1, "\n", "-"*45)
        #Print parameter set in model.in file
        modelinf = open(os.path.join(os.getcwd(),"model.in"),"w")
        writeline =''
        for j in range(nchrom): #parameter loop
            writeline += parname[j]+'\t'+str(ParameterValues[i][j])+'\n'       
        modelinf.writelines(writeline)
        modelinf.close()
        
        #Run command file (SWATedit, SWAT and extract exe files)
        if "win" in sys.platform.lower():
            os.system(SWATdir+'/nsga2_mid_prll.cmd')
        elif "linux" in sys.platform.lower():
            os.system(SWATdir+'/nsga2_mid.sh')
        else:
            os.system(SWATdir+'/nsga2_mid.cmd')

        #Read 'model.out' file
        modelrchf = open(os.path.join(SWATdir,"model.out"),'r')
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
            for k in range(0,nobjfunc_):
                objfuncav=0
                for j in range(0,nobjsite_):
                    objfuncav+=objectivefuncs[(k+j*nobjfunc_)]/nobjsite_
                newobjectivefuncs.append(objfuncav)  
        elif FuncOptAvr==2: #Average Objective Functions
            for k in range(0,nobjsite_):
                objsiteav=0
                for j in range(0,nobjfunc_):
                    objsiteav+=objectivefuncs[(k*nobjfunc_+j)]/nobjfunc_
                newobjectivefuncs.append(objsiteav)             
        else: #Do not average
            newobjectivefuncs = objectivefuncs
        #Add objective functions to population
        population["ind"][i]['fitness'] = newobjectivefuncs
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    #/*---------------------------* RANKING *------------------------------*/
    nsga2utilities.round_fitness(population)
    rankcon(population);
    return;
#-------------------------------------------------------------------------------
