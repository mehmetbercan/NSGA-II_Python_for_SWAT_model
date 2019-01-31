#Extract output.rch similar to observed_rch.txt struction in nsga2.in folder
import sys,datetime,os

#----------------- Class ------------------
class SWAT:
    def __init__(self,SWATtxtinoutFolderDirectory):
        self.__dir = SWATtxtinoutFolderDirectory
        #Read file. cio to get begining date and end date
        SWAT = open(os.path.join(SWATtxtinoutFolderDirectory,"file.cio"), "r")
        lines = SWAT.readlines()
        skipyear = int(lines[59][12:16])
        FCbeginyear = int(lines[8][12:16]) + skipyear #begining year
        FCendyear = FCbeginyear + int(lines[7][12:16])-1 - int(lines[59][12:16])#ending year
        if skipyear == 0:
            FCbeginday = int(lines[9][12:16])  #begining julian day
        else:
            FCbeginday = 1  #begining julian day
        FCendday = int(lines[10][12:16])  #ending julian day
        SWAT.close()
        FCbegindate = datetime.datetime(FCbeginyear, 1, 1) + datetime.timedelta(FCbeginday - 1)
        FCenddate = datetime.datetime(FCendyear, 1, 1) + datetime.timedelta(FCendday - 1)
        self.__fcbegindate = FCbegindate
        self.__fcenddate = FCenddate

    def DefineTimePeriod(self, BeginDate, EndDate):
        '''Date Format: "m/d/yyyy" '''
        self.__begindate = (datetime.datetime.strptime(BeginDate, "%m/%d/%Y"))
        self.__enddate = (datetime.datetime.strptime(EndDate, "%m/%d/%Y"))

    def ReadDailySWATrchfile(self, OutletSubbasinNumber,columnNumber):
        ''' saves the rch file column and corresponding date'''
        '''columnNumber=0 for area and 1 for the next and so on''' 
        SWAT_Directory = self.__dir        
        #Get file.cio dates
        FCbegindate = self.__fcbegindate
        FCenddate = self.__fcenddate
        #Get defined dates
        BD = self.__begindate
        ED = self.__enddate
        #Check if the date range of the SWAT run covers the defined date range
        if FCbegindate>BD:
            sys.exit("Error: The begin date in observed_rch.txt is earlier than SWAT run time span.")
        if FCenddate<ED:
            sys.exit("Error: The end date in observed_rch.txt is later than SWAT run time span.")
        #Read output.rch to get outflow
        SWAT = open(os.path.join(SWAT_Directory,"output.rch"), "r")
        SWATlines = SWAT.readlines()
        Date_Value = {}
        date = FCbegindate
        for line in range(9, len(SWATlines)):
            if float(SWATlines[line][6:11:]) == OutletSubbasinNumber: #get the values for outlet subbasin
                outflow = float(SWATlines[line][(columnNumber*12+26):(columnNumber*12+26+12):])
                Date_Value[date] = outflow
                date = date + datetime.timedelta(days=1)
        SWAT.close()
        # Get the values between defined date
        thedate=BD
        ddate=[]; dailystreamflow=[]
        while thedate != ED+datetime.timedelta(days=1):
            ddate+=[thedate]
            dailystreamflow+=[Date_Value[thedate]]
            thedate=thedate + datetime.timedelta(days=1)
        #Get array values for daily streamflow
        self.get_DailyDate = ddate
        self.get_DailyStreamflow = dailystreamflow
#Convert dictionary (Date_value) to array (date and value)---used up in the class
def DictionaryofDate_valuetoArrays(Date_value):
    '''Returns (array): date, value '''
    date = list(Date_value.keys())
    date.sort()
    value = []
    for d in date:
        value.append(Date_value[d])
    return date, value
#----------------- Class ------------------



#Read 'observed_rch.txt'
f = open(os.path.join(os.getcwd(),"NSGA2.IN","observed_rch.txt"), "r")
lines = f.readlines()
#Read Observed Streamflow
Outlet_Obsorder = {}; outlet = -99; nofdatapoint=-99; begindate=-99; enddate=-99; ColumnNo=-99
for i in range(0,len(lines)):
    try:
        if lines[i][0:10:]=='output_rch':
            outlet = int(lines[i][11:16].split(" ")[0])
            nofdatapoint = int(lines[i+1].split(" ")[0])
            begindate = lines[i+2].split(" ")[0].split("|")[0]
            enddate = lines[i+2].split(" ")[0].split("|")[1]
            ColumnNo = int(lines[i+2].split(" ")[0].split("|")[2])
            Obsorder = []
            for j in range((i+3), (i+3+nofdatapoint)):
                Obsorder.append(float(lines[j].split("\t")[0]))
            Outlet_Obsorder[outlet] = Obsorder
    except: sys.exit("ERROR: check the 'observed_rch.txt' file (the gage on line:"+str(i)+")")

outlets=list(Outlet_Obsorder.keys())
outlets.sort()
mlines=""
for outlet in outlets:
    #Read output.rch
    sw = SWAT("./")
    sw.DefineTimePeriod(begindate,enddate)
    sw.ReadDailySWATrchfile(outlet,ColumnNo) #(OutletSubbasinNumber,columnNumber)
    simstrflw=sw.get_DailyStreamflow
    #Get model.out lines
    for i1 in Outlet_Obsorder[outlet]:
        i=int(i1)-1
        mlines+=str(i+1)+"    "+str(simstrflw[i])+"\n"
#Print the model.out file
f=open(os.path.join(os.getcwd(),"model.out"),"w")
f.writelines(mlines)
f.close()

print ("Excraction is done for outlets: " , outlets)


















