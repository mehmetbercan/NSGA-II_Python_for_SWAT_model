#Moves Created GLUE files from NSGA-II to SWAT-CUP folders

import shutil


#---- Input ----(space on directory may cause problems during execution)
SWATtxtinoutDirectory = r"copy complete directory to SWAT folder that was set for NSGA-II"
GLUEdirectory = r'copy and paste complete directory create by SWAT-CUP (PROJECTNAME.Glue.SwatCup)'
#---------------

#Define paths
nsga2sourcePath=SWATtxtinoutDirectory+r'/SWAT-CUP_GLUE_files'
glueinPaht=GLUEdirectory+r'/GLUE.IN'
glueoutPaht=GLUEdirectory+r'/GLUE.OUT'

#Copy files
shutil.copy2(nsga2sourcePath+r'/glue.inf', glueinPaht+r'/glue.inf')
shutil.copy2(nsga2sourcePath+r'/glue_obs.dat', glueinPaht+r'/glue_obs.dat')
shutil.copy2(nsga2sourcePath+r'/var_file_name.txt', glueinPaht+r'/var_file_name.txt')
shutil.copy2(nsga2sourcePath+r'/modelpara.beh', glueoutPaht+r'/modelpara.beh')
shutil.copy2(nsga2sourcePath+r'/modelpara.beh', glueoutPaht+r'/modelpara.out')
shutil.copy2(nsga2sourcePath+r'/modelres.beh', glueoutPaht+r'/modelres.beh')

print 'GLUE files successfully copied into SWAT-CUP GLUE directories.'
