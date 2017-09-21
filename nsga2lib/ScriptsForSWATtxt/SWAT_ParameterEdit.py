#Changes SWAT files based on the values defined in model.in
import os, glob, shutil


DefaultDir = os.path.join(os.getcwd(),"Backup")
modelin = os.path.join(os.getcwd(),"model.in")

#------------------------------------ Change parameter values-------------------------------------
#Changes SWAT input files based on imet and parameter value
#-------------------------------------------------------------------------------------------------
def SWATparameterChange(TheFileDir, TheValue, IMET, TheLineNo, TheStartSpaceNo, TheSpaceLenght, multiopt="notmulti"):
    #Read, save in a library, delete
    File = open(TheFileDir, "r") # open the file
    lines = File.readlines() # Read lines
    Lines = {}
    for i in range(1, len(lines)+1):
        Lines[i] = lines[i-1]
    File.close()
    os.remove(TheFileDir) #remove the file

    # get the original value
    if IMET == 2 or IMET == 3 or multiopt=="multi":
        originalfile = DefaultDir + "/" + TheFileDir.split("/")[len(TheFileDir.split("/"))-1]
        File = open(originalfile, "r") # open the file
        lines = File.readlines() # Read lines
        oLines = {}
        for i in range(1, len(lines)+1):
            oLines[i] = lines[i-1]
        File.close()
    else:
        oLines = Lines
    originalvalue = []
    if multiopt=="multi":
        theline = oLines[TheLineNo][TheStartSpaceNo:-1:]
        for i in range(0, len(theline), TheSpaceLenght):
            originalvalue.append(float(theline[i:(i+TheSpaceLenght):]))
    else:
        originalvalue.append(float(oLines[TheLineNo][TheStartSpaceNo:(TheSpaceLenght+TheStartSpaceNo):]))

    # calculate the new value based on imet value (1-Replace, 2-Add, 3-%Add)
    TheNewValue = ""
    for oval in originalvalue:
        if IMET == 1: #replace
            TheNewValue += str(round((TheValue),4)).rjust(TheSpaceLenght," ")
        if IMET == 2: #add
            TheNewValue += str(round((TheValue + oval),4)).rjust(TheSpaceLenght," ")
        if IMET == 3: #percent add
            TheNewValue += str(round((TheValue/100.0 + 1.0)* oval,4)).rjust(TheSpaceLenght," ")
    #Change the line
    firstpart = Lines[TheLineNo][:TheStartSpaceNo] + TheNewValue
    secondpart = Lines[TheLineNo][len(firstpart):]
    Lines[TheLineNo] = firstpart + secondpart
    #rewrite the file
    File = open(TheFileDir, "w")
    for i in range(1, len(Lines)+1):
        File.write(Lines[i])
    File.close()
#-------------------------------------------------------------------------------------------------

#Read model.in
par_ValImet={}
f=open(modelin,"r")
lines=f.readlines()
for line in lines:
    if line=="\n":break
    columns=line.split("\n")[0].split("\t")
    imetc=columns[0].split("__")[0]
    if imetc=="v":imet=1
    if imetc=="a":imet=2
    if imetc=="r":imet=3
    par=columns[0].split("__")[1]
    value=float(columns[1])
    par_ValImet[par]=[value,imet]
f.close()
    
#Remove output files
try:os.remove("./output.hru"); os.remove("./output.sub")
except:pass

#change parameter
for par in par_ValImet.keys():
    value = par_ValImet[par][0]
    imet = par_ValImet[par][1]
    if par == "ALPHA_BF.gw": ##AlphaBf# `ok E=0.9997
        for path in glob.glob( os.path.join("./", '*.gw') ):
            SWATparameterChange(path, value, imet, 5, 0, 16)
        print("ALPHA_BF.gw--> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "BIOMIX.mgt": ##Biomix# `ok E=0.9998
        for path in glob.glob( os.path.join("./", '*.mgt') ):
            SWATparameterChange(path, value, imet, 10, 0, 16)
        print("BIOMIX.mgt---> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "CANMX.hru": ##Canmx#`ok E=0.9996
        for path in glob.glob( os.path.join("./", '*.hru') ):
            SWATparameterChange(path, value, imet, 9, 0, 16)
        print("CANMX.hru----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "CH_K2.rte": ##Ch_K2# `ok E=0.9990
        for path in glob.glob( os.path.join("./", '*.rte') ):
            SWATparameterChange(path, value, imet, 7, 0, 14)
        print("CH_K2.rte---->  Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "CH_N2.rte": ##Ch_N2# `ok E=0.9768 at chn2=1.0
        for path in glob.glob( os.path.join("./", '*.rte') ):
            SWATparameterChange(path, value, imet, 6, 0, 14)
        print("CH_N2.rte---->  Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "CN2.mgt": ##Cn2# `ok E=0.9990
        for path in glob.glob( os.path.join("./", '*.mgt') ):
            SWATparameterChange(path, value, imet, 11, 0, 16)
        print("CN2.mgt------> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "EPCO.hru": ##Epco# `ok E=0.9997
        for path in glob.glob( os.path.join("./", '*.hru') ):
            SWATparameterChange(path, value, imet, 11, 0, 16)
        print("EPCO.hru-----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "ESCO.hru": ##Esco# `ok E=0.9998
        for path in glob.glob( os.path.join("./", '*.hru') ):
            SWATparameterChange(path, value, imet, 10, 0, 16)
        print("ESCO.hru-----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "GW_DELAY.gw": ##GwDelay# `ok E=0.9997
        for path in glob.glob( os.path.join("./", '*.gw') ):
            SWATparameterChange(path, value, imet, 4, 0, 16)
        print("GW_DELAY.gw--> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "GW_REVAP.gw": ##GwRevap# `ok E=0.9993
        for path in glob.glob( os.path.join("./", '*.gw') ):
            SWATparameterChange(path, value, imet, 7, 0, 16)
        print("GW_REVAP.gw--> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "GWQMN.gw": ##Gwqmn# `ok E=0.9993
        for path in glob.glob( os.path.join("./", '*.gw') ):
            SWATparameterChange(path, value, imet, 6, 0, 16)
        print("GWQMN.gw-----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    #Not in ChangeparAll2009flowpar.dat (blai and slope missing in it)
    if par == "RCHRG_DP.gw": ##Rchrg_Dp# `NA E=0.8768
        for path in glob.glob( os.path.join("./", '*.gw') ):
            SWATparameterChange(path, value, imet, 9, 0, 16)
        print("RCHRG_DP.gw--> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "REVAPMN.gw": ##Revapmn# `ok E=0.9996
        for path in glob.glob( os.path.join("./", '*.gw') ):
            SWATparameterChange(path, value, imet, 8, 0, 16)
        print("REVAPMN.gw---> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SFTMP.bsn": ##Sftmp# `ok E=0.9996 (Not that No snow fall occurs in Eno --results were same in different values of sftmp)
        path = "./basins.bsn"
        SWATparameterChange(path, value, imet, 4, 0, 16)
        print("SFTMP.bsn----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SLSUBBSN.hru": ##Slsubbsn# `ok E=0.9996
        for path in glob.glob( os.path.join("./", '*.hru') ):
            SWATparameterChange(path, value, imet, 3, 0, 16)
        print("SLSUBBSN.hru-> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SMFMN.bsn": ##Smfmn# `ok E=0.9996 (Not that No snow fall occurs in Eno --results were same in different values of sftmp for testoriginal)
        path = "./basins.bsn"
        SWATparameterChange(path, value, imet, 7, 0, 16)
        print("SMFMN.bsn----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")
        
    if par == "SMFMX.bsn": ##Smfmx#  `ok E=0.9996 (Not that No snow fall occurs in Eno --results were same in different values of sftmp for testoriginal)
        path = "./basins.bsn"
        SWATparameterChange(path, value, imet, 6, 0, 16)
        print("SMFMX.bsn----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")
        
    if par == "SMTMP.bsn": ##Smtmp# `ok E=0.9996 (Not that No snow fall occurs in Eno --results were same in different values of sftmp)
        path = "./basins.bsn"
        SWATparameterChange(path, value, imet, 5, 0, 16)
        print("SMTMP.bsn----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SOL_ALB().sol": #SolAlb# #It does not change any result with change in sol file for test.
        for path in glob.glob( os.path.join("./", '*.sol') ):
            SWATparameterChange(path, value, imet, 17, 27, 12, "multi")
        print("SOL_ALB().sol> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SOL_AWC().sol": ##Sol_Awc# `ok E=0.9996
        for path in glob.glob( os.path.join("./", '*.sol') ):
            SWATparameterChange(path, value, imet, 10, 27, 12, "multi")
        print("SOL_AWC().sol> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SOL_K().sol": ##SolK# `ok E=0.9996
        for path in glob.glob( os.path.join("./", '*.sol') ):
            SWATparameterChange(path, value, imet, 11, 27, 12, "multi")
        print("SOL_K().sol--> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SOL_Z().sol": #SolZ# `ok E=0.9900--
        for path in glob.glob( os.path.join("./", '*.sol') ):
            SWATparameterChange(path, value, imet, 8, 27, 12, "multi")
        print("SOL_Z().sol--> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "SURLAG.bsn": ##Surlag# `ok E=0.9995
        path = "./basins.bsn"
        SWATparameterChange(path, value, imet, 20, 0, 16)
        print("SURLAG.bsn---> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "TIMP.bsn": ##Timp# `ok E=0.9996 (Not that No snow fall occurs in Eno --results were same in different values of sftmp for testoriginal)
        path = "./basins.bsn"
        SWATparameterChange(path, value, imet, 8, 0, 16)
        print("TIMP.bsn-----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")

    if par == "TLAPS.sub": ##Tlaps#    #It does not change any result with change in sol file for test.
        for path in glob.glob( os.path.join("./", '*.sub') ):
            SWATparameterChange(path, value, imet, 22, 0, 16)
        print("TLAPS.sub----> Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)")


###--------------Nitrogen Parameters----------------
##    if par == -1: ##CMN # `ok E=
##        path = SWATtxtinout + "/" + "basins.bsn"
##        SWATparameterChange(path, value, imet, 27, 0, 16)
##        #print "--CMN--- Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)"
##    if par == -2: ##N_UPDIS  # `ok E=
##        path = SWATtxtinout + "/" + "basins.bsn"
##        SWATparameterChange(path, value, imet, 28, 0, 16)
##        #print "--N_UPDIS--- Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)"
##    if par == -3: ##RSDCO # `ok E=
##        path = SWATtxtinout + "/" + "basins.bsn"
##        SWATparameterChange(path, value, imet, 34, 0, 16)
##        #print "--RSDCO--- Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)"
##    if par == 21: #Soil NO3# `ok E=
##        for path in glob.glob( os.path.join(SWATdir, '*.mgt') ):
##            SWATparameterChange(path, value, imet, 4, 27, 12, "multi")
##        #print "--Soil NO3----- Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)"
##    if par == 19: #Soil organic N# `ok E=
##        for path in glob.glob( os.path.join(SWATdir, '*.mgt') ):
##            SWATparameterChange(path, value, imet, 5, 27, 12, "multi")
##        #print "--Soil organic N----- Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)"
##    if par == -4: ##ERORGN # `ok E=
##        for path in glob.glob( os.path.join(SWATdir, '*.mgt') ):
##            SWATparameterChange(path, value, imet, 13, 0, 16)
##        #print "--ERORGN ----- Value: ", "%.5e"%value, "IMET: ", imet, " (1-Replace, 2-Add, 3-%Add)"
###--------------Nitrogen Parameters----------------
