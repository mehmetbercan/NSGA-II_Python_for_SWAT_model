#Creates Mutated population (non-dominating sorting and crowding distance functions)
import numpy, os
from math import floor


def round_parameters(pop_ptr, decimal = 6):
    popsize = len(pop_ptr["ind"])
    nchrom = len(pop_ptr["ind"][0]["xbin"])
    i=0
    while i < popsize:
        m = 0
        while m < nchrom:
            pop_ptr["ind"][i]["xbin"][m] = round(pop_ptr["ind"][i]["xbin"][m],decimal)
            m += 1
        i += 1
def round_fitness(pop_ptr, decimal = 3):
    popsize = len(pop_ptr["ind"])
    nfunc = len(pop_ptr["ind"][0]["fitness"])
    i=0
    while i < popsize:
        m = 0
        while m < nfunc:
            pop_ptr["ind"][i]["fitness"][m] = round(pop_ptr["ind"][i]["fitness"][m],decimal)
            m += 1
        i += 1
        
#-------------------------------------------------------------------------------
#/*This program subroutine is used to print the report*/
def report(pop1_ptr,pop2_ptr,igen,ngen,SWATdir,ncross,nmut):
    outputoutfile=os.path.join(SWATdir,"NSGA2.OUT","output.out")
    plotoutfile=os.path.join(SWATdir,"NSGA2.OUT","plot.out")
    outputfile=open(outputoutfile, "a")
    if (igen == ngen):
        plotfile=open(plotoutfile, "w")
        plotfile.writelines("# Feasible and Non-dominated Objective Vector\n");
            
    popsize = len(pop1_ptr["ind"])
    nchrom = len(pop1_ptr["ind"][0]["xbin"])
    nfunc = len(pop1_ptr["ind"][0]["fitness"])

    outputfile.writelines("\n\n---------------------------------------------------\n");
    outputfile.writelines("Generation No.     ->%d\n"%(igen));
    outputfile.writelines("------------------------------------------------------\n");

    outputfile.writelines("OldPop:  variables (binary %d)  fitness (%d)  rank cublen || MatePop: variables  fitness rank cublen\n"%(nchrom,nfunc));

    i=0
    while i < popsize:
        outputfile.writelines("\n------------------------------------------------\n");

        ptr1_b = pop1_ptr['ind'][i]['xbin'];
        ptr2_b = pop2_ptr['ind'][i]['xbin'];

        fptr = pop1_ptr['ind'][i]['fitness'];
        fptr1 = pop2_ptr['ind'][i]['fitness'];

        rptr = pop1_ptr['ind'][i]['rank'];
        rptr1 = pop2_ptr['ind'][i]['rank'];

        j=0
        while j < nchrom:
            outputfile.writelines("%f "%ptr1_b[j]);
            j+=1

        if (igen == ngen):
            j=0
            while j < nfunc:
                if (rptr1 == 1):
                    plotfile.writelines("%f\t"%(fptr1[j]));
                j+=1
            if (rptr1 == 1):
                plotfile.writelines("\n");

        fptr =  pop1_ptr['ind'][i]['fitness'];
        fptr1 = pop2_ptr['ind'][i]['fitness'];

        j=0
        while j < nfunc:
            outputfile.writelines("  %.4f"%fptr[j]);
            j+=1

        outputfile.writelines(" %d "%rptr);

        outputfile.writelines("%f "%pop1_ptr['ind'][i]['cub_len']);
        outputfile.writelines("|**|");

        j=0
        while j < nchrom:
            outputfile.writelines("%f "%ptr2_b[j]);
            j+=1
        j=0
        while j < nfunc:
            outputfile.writelines("  %f"%fptr1[j]);
            j+=1
        outputfile.writelines(" %d "%rptr1);

        outputfile.writelines(" %f "%pop2_ptr['ind'][i]['cub_len']);

        i+=1

    outputfile.writelines("\n--------------------------------------------------\n\n");
    outputfile.writelines("-------------------------------------------------------\n");
    
    if (igen == ngen):
        plotfile.close()
        outputfile.writelines("NO. OF CROSSOVER = %d\n"%ncross);
        outputfile.writelines("NO. OF MUTATION = %d\n"%nmut);
        outputfile.writelines("------------------------------------------------------------\n");
        outputfile.writelines("---------------------------------Thanks---------------------\n");
        outputfile.writelines("-------------------------------------------------------------\n");
    outputfile.close()
    return;
#-------------------------------------------------------------------------------




#decode--------------------------------------------------------------------------
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
#decode--------------------------------------------------------------------------









#typedef struct    /*Popuation Structure*/
def CreateDefaultPopulation(popsize,chrom,nchrom,nfunc):
    popltn={}
    popltn["maxrank"]=0        #/*Maximum rank present in the population*/
    popltn["rankno"]=numpy.zeros(2*popsize,int)#/*Individual at different ranks*/
    popltn["ind"]=[]    #/*Different Individuals*/
    for i in range(popsize):
        indvdl = {}
        indvdl["genes"]=numpy.zeros(chrom,int)#/*bianry chromosome*/
        indvdl["rank"]=0      #/*Rank of the individual*/
        indvdl["flag"]=0      #/*Flag for ranking*/
        indvdl["xbin"]=numpy.zeros(nchrom,float)#/*list of decoded value of the chromosome */
        indvdl["fitness"]=numpy.zeros(nfunc,float)#/*Fitness values */
        indvdl["fitness"][:] = 1e15 #To avoid being best fit in case of objective function error during initialization
        indvdl["cub_len"]=0.0 #/*crowding distance of the individual*/
        popltn["ind"].append(indvdl)
    return popltn




#Selection----------------------------------------------------------------------
#/*This is the function to get the different individuals selected*/
"""All mate_pop_ptr members are updated. Selection() selects two random member from old_pop_ptr and
puts best one based on ranking ("rank") and crowding distance ("cub_len") to mate_pop_ptr untill all mate_pop_ptr filled."""
def indzeros(chrom): #Creates individual with "genes", "rank" and "cub_len" with zero values. genes= a parameter set in a binary form, cub_len=crowding distance
    indvdl = {}
    indvdl["genes"]=numpy.zeros(chrom,int)#/*bianry chromosome*/
    indvdl["rank"]=0      #/*Rank of the individual*/
    indvdl["cub_len"]=0.0 #/*crowding distance of the individual*/
    return indvdl

def Selection(old_pop_ptr,pop2_ptr,warmup_random): #nselect() in deb's c code #select solutions from old population to the mate population. All members on mate population updated.
    popsize = len(old_pop_ptr["ind"])
    chrom = len(old_pop_ptr["ind"][0]['genes']) #chrom=number of bits of a parameter set

    r = popsize;
    s = chrom;

    indZeros = indzeros(chrom)

    k=-1;
    for n in range(popsize):
        k+=1
        j=0;j1 = 0;

        rnd2 = warmup_random.randomperc();
        rnd2 = popsize*rnd2;
        rnd = int(floor(rnd2));
        if(rnd == 0):rnd = popsize - k;
        if(rnd == popsize):rnd = abs(popsize-2)/2;

        # /*Select first parent randomly*/
        if rnd <= 0:j = indZeros #The population has max individual members in the c code but not here so last item there is zero
        else:j = old_pop_ptr["ind"][int(rnd-1)];

        rnd2 = warmup_random.randomperc();
        rnd2 = popsize * rnd2;
        rnd1 = int(floor(rnd2));
        if (rnd1 == 0):rnd1 = popsize - n;
        if(rnd1 == popsize):rnd1 = abs(popsize - 4)/2;

        #/*Select second parent randomly*/
        if rnd1 <= 0:j1 = indZeros #The population has max individual members in the c code but not here so last item there is zero
        else:j1 = old_pop_ptr["ind"][int(rnd1-1)];

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
#Selection----------------------------------------------------------------------



#Crossover------------------------------------------------------------------------
def crossover(new_pop_ptr,mate_pop_ptr,warmup_random,pcross,ncross): #New population changes based on mate population
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

#/* uniform crossover */
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
#Crossover-------------------------------------------------------------------------



#Mutation----------------------------------------------------------------------
#/* This is the module used to formulate the mutation routine*/
def Mutation(new_pop_ptr,warmup_random,pmut_b,nmut): #mutate() in deb's c code #New population mutates
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
#Mutation-----------------------------------------------------------------------




#Keepaliven---------------------------------------------------------------------------
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

def indcmp1(ptr1,ptr2,nfunc): #Used in NonDominatedSorting
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

def NonDominatedSorting(gen,popsize,global_pop_ptr,nfunc,globalpop,SWATdir):#grank
    grankoutfile=os.path.join(SWATdir,"NSGA2.OUT","g_rank_record.out")
    gr = open(grankoutfile,"a");
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



#/* This is the file used to sort the dummyfitness arrays */
def gsort(rnk,sel,popsize,globalpop): #Used in CreateMatePopulation
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
        globalpop['flag'][int(a)] = 1;
        i+=1
    return;

def sort(m1,fpara1): #Used in CrowdingDistance
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

def CrowdingDistance(rnk,popsize,nfunc,globalpop,fpara1):#gshare
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
            print("Something wrong in keepaliven.h (CrowdingDistance)\n");
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
            globalpop['cub_len'][int(a)] += length[i][1]
            i+=1;
        j+=1

    return;


#/*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
#/*Elitism And Sharing Implemented*/ #### Nondominated sorting, crowding distances
def CreateMatePopFromNewandOldPops(pop1_ptr,pop2_ptr,pop3_ptr,gen,SWATdir): #keepaliven
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

    #####-------- Non-dominated sorting --------#####  #/*Finding the global ranks */ 
    NonDominatedSorting(gen,popsize,global_pop_ptr,nfunc,globalpop,SWATdir);

    m = globalpop['maxrank'];

    #####-------- Crowding Distance --------##### #/* Sharing the fitness to get the dummy fitness */
    i=0
    while i<m:
        CrowdingDistance(i+1,popsize,nfunc,globalpop,fpara1);
        i+=1

    poolf = popsize;
    pool = 0;
    
    #/*Initializing the flags of population to zero */
    i=0
    while i<(2*popsize):
        globalpop['flag'][i] = 0;
        i+=1
    #####-------- Creating Mate Population --------##### #// decide which all solutions belong to the pop3
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
    #####-------- Rest: Copying selected solutins to mate population --------#####
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
#Keepaliven-------------------------------------------------------------------------
