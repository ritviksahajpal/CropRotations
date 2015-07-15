################################################################
# May 3, 2010
# CropRotations.py
#
# Produces crop rotation patterns from CDL data
#################################################################
# Import system modules
import sys, string, os, random, time, pdb, math, operator, arcgisscripting, csv, glob, logging, random

gp = arcgisscripting.create()   # Create the Geoprocessor object
gp.CheckOutExtension("spatial") # Check out any necessary licenses
gp.AddToolbox("C:/Program Files (x86)/ArcGIS/ArcToolBox/Toolboxes/Spatial Analyst Tools.tbx")
gp.AddToolbox("C:/Program Files (x86)/ArcGIS/ArcToolBox/Toolboxes/Data Management Tools.tbx")
gp.overwriteoutput = True

base_dir     = 'C:\\Users\\ritvik\\Documents\\PhD\\Projects\\CropIntensity\\'
analysis_dir = 'C:\\Users\\ritvik\\Documents\\PhD\\Projects\\CropRotationsUSA\\src\\'

#######################################################################
# USER SPECIFIED PARAMETERS
# The following parameters are user specified and might need to be changed
# prj_dir: Top level directory where all the code, data and outputs are stored
# data_dir: Contains the input data
# out_dir: Contains the output and intermediate analysis stuff
# inp_raster_files: The list of Crop data layer (CDL) files for each year
# maxCropRotations: The total number of crop rotations we want to determine
#######################################################################
# prj_dir defines the project space
#prj_dir    = 'C:\\Documents and Settings\\Ritvik\\My Documents\\Projects\\CropRotations\\'
state = ''
max_crop_rot = 10
min_area_rot = 0.80# data_dir contains the input data
#data_dir   = prj_dir + 'Data\\land_use_land_cover_NASS_CDL_'+state.lower()+'\\land_use_land_cover\\NASS_'+state+'\\'
#data_dir = prj_dir+os.sep+'Data'+os.sep+state+os.sep+state+os.sep
# inp_raster_files is the list of CDL files for each year
# IL: 0.02, SD: 0.05, ND: 0.05, WI:0.05; IN: 0.05; MI: 0.05 NY = 0.05 PA: 0.02
min_perc_for_rot = 0.1
min_perc_for_comb = 0.1

RAS = 'ras_'
OUT_RAS = 'out_ras_'
COMB = 'comb_'
YRS = 'yrs'
FRST = 'frst_'
ROT_FRST_GRASS = 'all_'
CDL_ROT = 'cdl_rot_'
ROT = 'rot_'
BASE_YR_WOUT_URBAN = 'wout_urb_'
EXISTING_ROTATIONS = 'existing_rotations_CSite.csv'
#EXISTING_ROTATIONS = 'latest_rot.csv'
HUMAN_READABLE_ROTATIONS = 'human_readable_rotations.csv'
PRODUCT = 'product_'
RAND_MIN = 10001
RAND_MAX = 20000
ROTATION_STATS = 'rot_stats.csv'
STATES = 'states.txt'
CDL_VAL_NAME_MAP = 'CDL_Value_Name_Mapping.csv'

# Determine the crop rotation sequence
#############################################################################
# CROPS                                 HERBS                                   WOODY
# 1    Corn                            58    Clover/Wildflowers              63    Wood lands
# 3    Rice                            61    Fallow/ Idle crop land          141    NLCD - Deciduous Forest
# 4    Sorghum                         62    Grass/Pasture/Non-Ag            142    NLCD - Evergreen Forest
# 5    Soybean                         171    NLCD - Grassland Herbaceous     143    NLCD - Mixed Forest
# 6    Sunflowers                      181    NLCD - Pasture/Hay              152    NLCD - Shrubland
# 21    Barley                          37      Other Hays
# 22    Durum Wheat
# 23    Spring Wheat
# 24    Winter Wheat
# 25    Other Small Grains
# 26    Win. Wht./Soyb. Dbl. cropped
# 27    Rye
# 28    Oats
# 29    Millet
# 31    Canola
# 32    Flaxseed
# 33    Safflower
# 36    Alfalfa
# 41    Sugerbeets
# 42    Dry Beans
# 52    Lentils
# 53    Peas
# CATEGORIES
# All the crops are categorized separately. All the herbs (57, 59, 62, 152 and 171)
# are lumped together as  'grassland' and all the woody areas (141, 142, 143, 63)
# are lumped together as 'forest'
# RULES
# 1. Any value (say 1: Corn) occuring more than 70% of the years is taken to be continuous (e.g. continuous corn in this case)
# 2. Two values occuring in tandem (say 1: Corn, and 5: Soybean) constitute a rotation if they occur more than 70% of the years
# 3. Three values occuring in tandem (say 1, 5 and 24: W. Wheat) constitute a rotation if they occur more than 70% of the years
# 4. 0: Background; 121,122, 131: Developed/Open Space; 111: Open water; 190: Woody wetland; These categories are to be excluded entirely
# 5. Remaining areas should be covered by the dominant cover type from 1,2,3
##############################################################################
# Find raster names
# Find the right file type
# Find the year name in the string
#http://stackoverflow.com/questions/2089036/sorting-csv-in-python
def convert(types, values):
    return [t(v) for t, v in zip(types, values)]

#http://stackoverflow.com/questions/2089036/sorting-csv-in-python
def sort_csv(csv_filename, types, sort_key_columns):
    """sort (and rewrite) a csv file.
    types:  data types (conversion functions) for each column in the file
    sort_key_columns: column numbers of columns to sort by"""
    data = []
    f = open(csv_filename, 'rb')
    for row in csv.reader(f):
        data.append(convert(types, row))
    data.sort(key=operator.itemgetter(sort_key_columns))
    f = open(csv_filename, 'wb')
    csv.writer(f).writerow(['FROM', 'TO', 'VALUE'])
    csv.writer(f).writerows(data)
    
# Returns the number of string positions that are different in two crop rotations a and b
# E.g. a = [1, 2, 3, 5]; b = [2, 4, 3, 5] differ in 2 positions (1st and 2nd)
def lev(a, b):
        cnt = 0
        for i in range(len(a)):
                if(a[i] <> b[i]):
                        cnt = cnt + 1
            
        return cnt

# Returns true if two strings are permutations.
def checkPermutations(a, b):
        return (sorted(a) == sorted(b))
    
def writeCDLData(cdl_comb, fileName):
        fileHandle = open(fileName, 'w')
        for i in range(len(cdl_comb)):
                for j in range(len(cdl_comb[i])):
                    fileHandle.write(str(cdl_comb[i][j]))
                    if (j < len(cdl_comb[i])-1):
                        fileHandle.write(', ')
                    else:
                        fileHandle.write('\n')        
        fileHandle.close()        

def returnLandUses(file_name):
    cdl_file = open(file_name, 'r')
    lines = cdl_file.readlines()
    elem = []
    for line in lines:
        #print line.split()[0]
        elem.append(line.split()[0])
    
    cdl_file_str = ""
    val = "\"Value\" = "
    
    for i in range(len(elem)):
        cdl_file_str+=(val)
        cdl_file_str+=elem[i]
        if i < len(elem)-1:
            cdl_file_str+=" OR "
    
    return cdl_file_str

def writeReClassFile(file_hndl, from_str, value_str):
    # Assume that the to_str is the same as from_str
    file_hndl.write(str(from_str)+','+str(from_str)+','+str(value_str)+'\n')

def synchronizeRotation(vec, yrs_in_rot, first_yr, prev_first_yr):
    new_vec = []
    count = 0
    if(first_yr <> prev_first_yr):
        for j in range(prev_first_yr, first_yr):
            new_vec.append(vec[yrs_in_rot-(first_yr-prev_first_yr)+count])
            count += 1
        for i in range(num_rot - (first_yr-prev_first_yr)):
                new_vec.append(vec[i])
        #new_vec.append(vec[:(yrs_in_rot - (first_yr-prev_first_yr))])
    else:
        new_vec = vec
            
    return new_vec

def computeCropRotations(crp_rot_files, first_yr, last_yr):    
    rasters_to_delete = []
    # numYears is the number of years for which we have crop rotation data
    yrs_in_rot = len(crp_rot_files)

    # Extract the simplified raster data
    all_inp_rasters = '"'
    print 'Processing ' + str(yrs_in_rot) + ' rasters located in ' + data_dir
    
    for i in range(yrs_in_rot):
        if (i == 0):
            all_inp_rasters = '"'+out_dir+OUT_RAS+state+str(i)
        elif (i == yrs_in_rot -1):
            all_inp_rasters = all_inp_rasters+'; '+out_dir+OUT_RAS+state+str(i)+'"'
        else:
            all_inp_rasters = all_inp_rasters+'; '+out_dir+OUT_RAS+state+str(i)
    
        print '\tExtracting attributes for raster ' + os.path.split(crp_rot_files[i])[1] 
        # Simplify each of the rasters: i.e.extract the crop information
        cdl_2_str = returnLandUses(analysis_dir+'crops.txt')
        gp.ExtractByAttributes_sa(crp_rot_files[i], cdl_2_str, out_dir+RAS+state+str(i)) 
        rasters_to_delete.append(out_dir+OUT_RAS+state+str(i))
        rasters_to_delete.append(out_dir+RAS+state+str(i))
    
    base_yr_data  = crp_rot_files[base_yr-first_yr]
    # Remove the urban areas from raster for base_yr
    print '\tRemoving urban areas for '+os.path.split(base_yr_data)[1]
    urban_str = returnLandUses(analysis_dir+'urban.txt')
    gp.SetNull_sa(base_yr_data, base_yr_data, out_dir+BASE_YR_WOUT_URBAN+state, urban_str)
    rasters_to_delete.append(out_dir+BASE_YR_WOUT_URBAN+state)

    # For each raster replace the NoData values with cell values from raster for base_yr        
    for i in range(yrs_in_rot):
        if i <> ( (base_yr-first_yr) + ((yrs_in_rot-1)-(last_yr-first_yr)) ):
            print '\tReplacing the NoData cells for ' + os.path.split(crp_rot_files[i])[1]
            soma_exp = "CON(isnull("+out_dir+RAS+state+str(i) + "),"+base_yr_data+","+out_dir+RAS+state+str(i)+")"
            gp.SingleOutputMapAlgebra_sa(soma_exp, out_dir+OUT_RAS+state+str(i))
        else:
            gp.CopyRaster_management(out_dir+RAS+state+str(i), out_dir+OUT_RAS+state+str(i))
    
    # Combine the simplified rasters
    comb_rasters = out_dir+COMB+state+'_'+str(yrs_in_rot)+YRS
    print 'Merging all CDLs'
    gp.Combine_sa(all_inp_rasters, comb_rasters)
    # For each line of the attribute table in the comb_rasters file, determine the appropriate crop rotation
    rows       = gp.searchcursor(comb_rasters)
    row        = rows.next()
    cdl_comb = []
    while row <> None:
            cdl_info_for_pixel = []
            for i in range(yrs_in_rot):
                cdl_info_for_pixel.append(row.GetValue(OUT_RAS+state+str(i)))
    
            # Create a tuple containing cdlInfoForPixel and the number of pixels with that crop rotation type
            cdl_list = list(cdl_info_for_pixel)
            cdl_list.append(row.COUNT)
            # A '1' in the right most column(USED bit) indicates that the crop rotation will be used.In subsequent simplifications,
            # most of the 1's will be turned to 0's
            # E.g. cdlInfoForPixel can look like: 1, 5, 1, 5, 300, 1 (corn, soy, corn, soy, COUNT, 1)
            cdl_list.append(1)
            cdl_list.append(row.VALUE)
            cdl_comb.append(cdl_list)
            row = rows.next()
    
    # Sort the data based on the COUNT which is in column index yrs_in_rot (0 based counting)
    cdl_comb = sorted(cdl_comb, key=operator.itemgetter(yrs_in_rot), reverse = True)
    writeCDLData(cdl_comb, out_dir + 'cdl_comb_'+state+'.csv')
    
    # Find the forested Land in base_yr
    forest_data = out_dir+'forest'+str(base_yr)
    print 'Extracting forest land in '+str(base_yr)
    cdl_f_str = returnLandUses(analysis_dir+'forests.txt')
    gp.ExtractByAttributes_sa(crp_rot_files[yrs_in_rot-2], cdl_f_str, forest_data)
    # Find the grassland/pasture/herbaceus vegetation land in base_yr
    grasslands_data = out_dir+'grass'+str(base_yr)
    print 'Extracting grassland in '+str(base_yr)
    cdl_g_str = returnLandUses(analysis_dir+'grasslands.txt')
    gp.ExtractByAttributes_sa(crp_rot_files[yrs_in_rot-2], cdl_g_str, grasslands_data)
    
    # Retain a maximum of maxCropRotations number of rotations. For this:
    # 1. Remove crop rotations which occupy less than 5% of the area
    # 2. Combine crop rotations which are very similar to each other i.e. by looking at the Levenshtein distance
    # 3. Combine crop rotations which differ from each other temporally e.g. [corn,soybean,corn] and [soybean,corn,soybean] are combined.
    # Simplifications:
    
    # 2. Loop through the file, find all the lines which are offset by less than 25% and combine them
    print 'Simplifying crop rotations'
    # Also create the reclassification table
    recl_1 = open(out_dir+'recl_1.csv', 'w')
    #most_common_lu = cdl_comb[0][yrs_in_rot+2] # Basically the value field from the raster attribute table
    most_pixels = cdl_comb[0][yrs_in_rot]
    total_pixels = 0
    for i in range(len(cdl_comb)):
        total_pixels += cdl_comb[i][yrs_in_rot]
    
    cur_num_rot = 0
    cur_perc_land_in_rot = 0.0
    perc_simplified_land = 0.0
    
    for i in range(len(cdl_comb)-2):
        if (cdl_comb[i][yrs_in_rot+1]!=0):
            cur_num_rot += 1
            cur_perc_land_in_rot += float(cdl_comb[i][yrs_in_rot])/float(total_pixels)
            
            if (cur_num_rot<max_crop_rot and cur_perc_land_in_rot<min_area_rot):
                for j in range(i+1, len(cdl_comb)):
                    if (cdl_comb[j][yrs_in_rot+1] != 0):
                        # Compute the levenshtein distance between crop rotation types i and j
                        diff_btwn_rot = float(lev(cdl_comb[i][:yrs_in_rot], \
                                            cdl_comb[j][:yrs_in_rot]))/float(yrs_in_rot)
                        # Check if the strings are permutations of each other
                        is_permutation = checkPermutations(cdl_comb[i][:yrs_in_rot], \
                                            cdl_comb[j][:yrs_in_rot])
                        # If the levenshtein distance is less than 25% of the rotation length, then
                        # the rotations are similar enough to combine. Also combine if the rotations
                        # are just permutations of one another
                        if (((diff_btwn_rot <= 1.0/yrs_in_rot) and (float(cdl_comb[j][yrs_in_rot])/float(most_pixels) <= min_perc_for_comb)) \
                                or (is_permutation and \
                                    float(cdl_comb[j][yrs_in_rot])/float(most_pixels) <= min_perc_for_rot)):
                            cdl_comb[i][yrs_in_rot] = cdl_comb[i][yrs_in_rot] +cdl_comb[j][yrs_in_rot]
                            # Set the USED bit of rotation 2 to 0
                            cdl_comb[j][yrs_in_rot+1] = 0
                            # In the reclassification file, set the reclassification bit to the rotation 1
                            writeReClassFile(recl_1, cdl_comb[j][yrs_in_rot+2], cdl_comb[i][yrs_in_rot+2])
                            cdl_comb[j][yrs_in_rot+2] = cdl_comb[i][yrs_in_rot+2]
                            perc_simplified_land += (1.0-diff_btwn_rot)*cdl_comb[j][yrs_in_rot]/total_pixels
                            #cur_perc_land_in_rot += cdl_comb[j][yrs_in_rot]/total_pixels
            else:
                for j in range(i+1, len(cdl_comb)):
                    if (cdl_comb[j][yrs_in_rot+1]!=0):
                        cdl_comb[0][yrs_in_rot] = cdl_comb[0][yrs_in_rot] +cdl_comb[j][yrs_in_rot]
                        # Set the USED bit of rotation 2 to 0
                        cdl_comb[j][yrs_in_rot+1] = 0
                        # In the reclassification file, set the reclassification bit to the rotation 1              
                        writeReClassFile(recl_1, cdl_comb[j][yrs_in_rot+2], cdl_comb[0][yrs_in_rot+2])
                        cdl_comb[j][yrs_in_rot+2] = cdl_comb[0][yrs_in_rot+2]
                        diff_btwn_rot = float(lev(cdl_comb[0][:yrs_in_rot], cdl_comb[j][:yrs_in_rot]))/float(yrs_in_rot)
                        perc_simplified_land += (1.0-diff_btwn_rot)*cdl_comb[j][yrs_in_rot]/total_pixels
                        #cur_perc_land_in_rot += cdl_comb[j][yrs_in_rot]/total_pixels                                  
    
    float_str1 = "%.2f" % (cur_perc_land_in_rot*100)
    float_str2 = "%.2f" % (perc_simplified_land*100)
    float_str3 = "%.2f" % ((perc_simplified_land+cur_perc_land_in_rot)*100)
    logger.info('% of rotations which match CDL data exactly: '+str(float_str1)+'%')
    logger.info('% of rotations which have been simplified: '+str(float_str2)+'%')
    logger.info('% net accuracy: '+str(float_str3)+'%')
    stat_writer.write(str(state)+', '+str(cur_num_rot)+', '+str(first_yr)+', '+str(last_yr)+', '+str(float_str1)+', '+str(float_str3)+'\n')
    stat_writer.flush()


    cdl_comb = sorted(cdl_comb, key=operator.itemgetter(yrs_in_rot+1), reverse = True)                
    writeCDLData(cdl_comb,out_dir+CDL_ROT+state+'.csv')
    recl_1.close()
    
    sort_csv(out_dir+'recl_1.csv', (int, int, int), 0)
    
    # Reclassify
    rot_data = out_dir+ROT+state+'_'+str(yrs_in_rot)+YRS
    gp.ReclassByTable_sa(comb_rasters, out_dir+'recl_1.csv', "FROM", "TO", "VALUE", rot_data, "DATA")
        
    # read existing rotations
    prev_rot_file = csv.reader(open(analysis_dir+EXISTING_ROTATIONS,'a+'))
    add_more_rot = []
    prev_rot = {}
    for str_row in prev_rot_file:
        row = [int(x) for x in str_row]
        #yrs_in_prev_rot = len(row)-1
        for i in range(1,len(row)-1):
            prev_rot[row[0]] = row[1:len(row)] 
                                                                                                                
    # Extract the existing rotations
    cur_rot = {}
    for i in range(len(cdl_comb)):                                                                            
        if(cdl_comb[i][yrs_in_rot+1]>0):
            #for j in range(yrs_in_rot):                                                                                                                                                                 
            sync_vec = synchronizeRotation(cdl_comb[i][0:yrs_in_rot], yrs_in_rot, first_yr, prev_first_yr)
            cur_rot[cdl_comb[i][yrs_in_rot+2]] = sync_vec
    
    # Compare the newly created rotations with the existing rotations and reclassify
    # 1. Iteratre through cur_rot and compare each cur_rot with the prev_rot
    # 2. If a cur_rot matches a prev_rot then put an entry in the reclass file
    # 3. If not, then add the cur_rot to the prev_rot file
    recl_cur_to_prev = open(out_dir+'recl_cur_to_prev.csv', 'wb')
    k1=k2=v1=v2=[]

    remove_cur_rot_items = []
    for k1,v1 in cur_rot.iteritems():
        for k2,v2 in prev_rot.iteritems():
            if(v1 == v2):
                # write in reclass file
                recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(k2)+'\n')
                remove_cur_rot_items.append(k1)
                
    for i in range(len(remove_cur_rot_items)):
        cur_rot.pop(remove_cur_rot_items[i])
        
    for k1,v1 in cur_rot.iteritems():
        add_cur_rot = True
        for k2,v2 in prev_rot.iteritems():
            if(k1==k2):
                tmp = k1 + random.randrange(RAND_MIN,RAND_MAX)
                lis = tmp, v1
                recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(tmp)+'\n')
                add_more_rot.append(lis)
                add_cur_rot = False
        if add_cur_rot == True:
            lis = k1, v1
            recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(k1)+'\n')
            add_more_rot.append(lis)
                        

#    for k1,v1 in cur_rot.iteritems():
#        v1_in_v2 = False
#        for k2,v2 in prev_rot.iteritems():
#            if(v1 == v2):
#                # write in reclass file
#                recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(k2)+'\n')
#                v1_in_v2 = True
#            elif(v1 <> v2 and k1 == k2):
#                #write in prev_rot_file
#                #tmp = k1+random.randrange(RAND_MIN,RAND_MAX)
#                tmp = k1 + max_prev_rot_id
#                lis = tmp, v1
#                recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(tmp)+'\n')
#                print '------'
#                print k1, tmp, v1
#                add_more_rot.append(lis)
#        if (v1_in_v2 == False):
#            # write in prev_rot_file
#            lis = k1, v1
#            recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(k1)+'\n')
#            add_more_rot.append(lis)
#        
    
    recl_cur_to_prev.close()
    sort_csv(out_dir+'recl_cur_to_prev.csv', (int, int, int), 0)
    
    append_to_prev_rot_file = open(analysis_dir+EXISTING_ROTATIONS,'a+')
    for j in range(len(add_more_rot)):
        # e.g. add_more_rot = (23, [24, 61, 61, 24])
        tmp = add_more_rot[j]
        append_to_prev_rot_file.write(str(tmp[0])+', ')
        for i in range(len(tmp[1])):
            append_to_prev_rot_file.write(str(tmp[1][i]))
            if i < (len(tmp[1])-1):
                append_to_prev_rot_file.write(', ')
            else:
                append_to_prev_rot_file.write('\n')
    append_to_prev_rot_file.close()
    
    # Reclassify
    final_rot_data = out_dir+PRODUCT+state+'_'+str(yrs_in_rot)
    gp.ReclassByTable_sa(rot_data, out_dir+'recl_cur_to_prev.csv', "FROM"   , "TO", "VALUE", final_rot_data, "DATA")
    rot_and_frst = out_dir+FRST+state+'_'+str(yrs_in_rot)+YRS
    rot_frst_grass = out_dir+ROT_FRST_GRASS+state+'_'+str(yrs_in_rot)+YRS
    # Append the forest data to the remaining set
    soma_exp = "CON(isnull(" + forest_data + ")," + final_rot_data + "," + forest_data + ")"
    gp.SingleOutputMapAlgebra_sa(soma_exp, rot_and_frst)
    rasters_to_delete.append(forest_data)
    rasters_to_delete.append(grasslands_data)
    rasters_to_delete.append(rot_and_frst)
    # Append the grasslands data to the remaining set
    soma_exp = "CON(isnull(" + grasslands_data + ")," + rot_and_frst + "," + grasslands_data + ")"
    gp.SingleOutputMapAlgebra_sa(soma_exp, rot_frst_grass)
    
    for i in range(len(rasters_to_delete)):
        print 'Deleting raster '+os.path.split(rasters_to_delete[i])[1]
        gp.Delete_management(rasters_to_delete[i], "")


print "Started processing at :", time.strftime("%H:%M:%S", time.localtime())
for subdir, dir_list, files in os.walk(base_dir):
    break

state_file = open(analysis_dir+STATES, 'rb')
lines = state_file.readlines()
data_dir = base_dir
stat_writer = open(data_dir+'Output'+os.sep+ROTATION_STATS,'wb')
stat_writer.write('State, # Rot, Start Yr, End Yr, % CDL Match, % Overall Match\n')
stat_writer.flush()

for line in lines:
    frst = False
    prev_first_yr = 2007
    first_yr = 2007
    last_yr  = 2010
    base_yr  = 2008 # This is the raster which will be used to fill in the gaps
    num_rot = 4
    state = line.split()[0]
    state_ras_files = []
    # out_dir contains the output and intermediate analysis stuff
    out_dir  = data_dir + 'Output'+os.sep+state+os.sep
    if not os.path.exists (out_dir):
        os.makedirs(out_dir)
    cur_time = str(time.strftime("%H:%M:%S",time.localtime())).replace(':','_')
    LOG_FILENAME = data_dir + os.sep+'Output'+os.sep+'Log_CropRotations'+cur_time+'.log'
    logging.basicConfig(filename = LOG_FILENAME, level=logging.DEBUG)
    # Logging levels are DEBUG, INFO, WARNING, ERROR, and CRITICAL
    logger = logging
    range_of_yrs = []
    for i in range(last_yr-first_yr+1):
        range_of_yrs.append(first_yr+i)
    print 'Evaluating '+state
    logger.info('CROP ROTATIONS FOR '+state)
    
    for j in range(len(range_of_yrs)):
        for position, item in enumerate(dir_list):
            if (str(range_of_yrs[j]) in item):
                list_files = glob.glob(base_dir+dir_list[position]+os.sep+'*_'+state+'_*'+str(range_of_yrs[j])+'*.tif')
                if list_files:
                    if frst == False:
                        first_yr = range_of_yrs[j]
                        frst = True 
                    state_ras_files.append(''.join(list_files))

    if len(state_ras_files) > num_rot:
        state_ras_files = state_ras_files[:num_rot]
        last_yr = first_yr + num_rot - 1
    
    logger.info('List of raster files ')    
    for k in range(len(state_ras_files)):
        logger.info(state_ras_files[k])
    logger.info('First Year '+str(first_yr))
    logger.info('Last Year '+str(last_yr))
    
    computeCropRotations(state_ras_files, first_yr, last_yr)
            
cdl_map_file = csv.reader(open(base_dir+'Output'+os.sep+CDL_VAL_NAME_MAP, 'r'))
cdl_map = {}
for row in cdl_map_file:
    cdl_map[int(row[0])] = row[1]

# Write a use friendly version of the EXISTING_ROTATIONS file
prev_rot_file = csv.reader(open(analysis_dir+EXISTING_ROTATIONS,'r'))
human_readable_rot = open(analysis_dir+HUMAN_READABLE_ROTATIONS,'wb')
    
for str_row in prev_rot_file:
    line = []
    row = [int(x) for x in str_row]
    #yrs_in_prev_rot = len(row)-1
    line.append(row[0])

    for i in range(1,len(row)):
        line.append(cdl_map[row[i]])

    print str_row
    print line
    csv.writer(human_readable_rot).writerow(line)

print "Finished processing at :", time.strftime("%H:%M:%S", time.localtime())