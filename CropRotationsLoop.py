import sys, string, os, random, time, pdb, math, operator, arcgisscripting, csv, glob

analysis_dir = 'C:\\Users\\ritvik\\Documents\\PhD\\Projects\\CropRotationsUSA\\src\\'
base_dir     = 'C:\\Users\\ritvik\\Documents\\PhD\\Projects\\CropIntensity\\'
first_yr = 2006
last_yr  = 2010
num_rot = 4

print "Started processing at :", time.strftime("%H:%M:%S", time.localtime())
for subdir, dir_list, files in os.walk(base_dir):
    break

range_of_yrs = []
for i in range(last_yr-first_yr+1):
    range_of_yrs.append(first_yr+i)

state_file = open(analysis_dir+'states.txt', 'rb')
lines = state_file.readlines()
frst = False
for line in lines:
    state_name = '_'+line.split()[0]+'_'
    #print state_name
    state_ras_files = []
    for j in range(len(range_of_yrs)):
        for position, item in enumerate(dir_list):
            if (str(range_of_yrs[j]) in item):
                list_files = glob.glob(base_dir+dir_list[position]+os.sep+'*'+state_name+'*'+str(range_of_yrs[j])+'*.tif')
                if list_files:
                    if frst == False:
                        first_yr = range_of_yrs[j]
                        frst = True 
                    state_ras_files.append(list_files)

    if len(state_ras_files) > num_rot:
        state_ras_files = state_ras_files[:num_rot]
        last_yr = first_yr + num_rot
        

print "Finished processing at :", time.strftime("%H:%M:%S", time.localtime())