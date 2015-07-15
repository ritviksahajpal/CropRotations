################################################################
# May 3, 2010
# CropRotations.py
#
# Produces crop rotation patterns from CDL data
#################################################################
# Import system modules
import sys, string, os, random, time, pdb, math, operator, arcgisscripting, csv, glob, logging, random, shutil, arcpy

gp = arcgisscripting.create(9.3)   # Create the Geoprocessor object
gp.CheckOutExtension("spatial") # Check out any necessary licenses
gp.AddToolbox("C:/Program Files (x86)/ArcGIS/Desktop10.0/ArcToolbox/Toolboxes/Spatial Analyst Tools.tbx")
gp.AddToolbox("C:/Program Files (x86)/ArcGIS/Desktop10.0/ArcToolbox/Toolboxes/Data Management Tools.tbx")
gp.overwriteoutput = True
gp.extent="MAXOF"

base_dir     = 'C:\\Users\\ritvik\\Documents\\PhD\\Projects\\CropIntensity\\'
date = 'oct7_2011'
tag = '_'+date+'_all_2008_2010'
max_crop_rot = 15
min_area_rot = 0.70# base_dir contains the input data
STATES = 'states_remaining.txt'
#EXISTING_ROTATIONS = 'existing_rotations_CSite.csv'
EXISTING_ROTATIONS = date+'.csv'
output_dir   = base_dir+'Output'+tag
analysis_dir = base_dir+os.sep+'src\\'
CDL_VAL_NAME_MAP = 'CDL_Value_Name_Mapping.csv'
DELETE = True

#######################################################################
# USER SPECIFIED PARAMETERS
# The following parameters are user specified and might need to be changed
# prj_dir: Top level directory where all the code, data and outputs are stored
# base_dir: Contains the input data
# out_dir: Contains the output and intermediate analysis stuff
# inp_raster_files: The list of Crop data layer (CDL) files for each year
# maxCropRotations: The total number of crop rotations we want to determine
#######################################################################
# prj_dir defines the project space
#prj_dir    = 'C:\\Documents and Settings\\Ritvik\\My Documents\\Projects\\CropRotations\\'
state = ''
#base_dir   = prj_dir + 'Data\\land_use_land_cover_NASS_CDL_'+state.lower()+'\\land_use_land_cover\\NASS_'+state+'\\'
#base_dir = prj_dir+os.sep+'Data'+os.sep+state+os.sep+state+os.sep
# inp_raster_files is the list of CDL files for each year
# IL: 0.02, SD: 0.05, ND: 0.05, WI:0.05; IN: 0.05; MI: 0.05 NY = 0.05 PA: 0.02
min_perc_for_rot = 0.1
min_perc_for_comb = 0.1

RAS = 'ras_'
OUT_RAS = 'out_ras_'
COMB = 'comb_'
YRS = 'yrs'
FRST = 'frst_'
ROT_FRST_GRASS = 'rfg_'
ROT_FRST_GRASS_URB = 'all_'
CDL_ROT = 'cdl_rot_'
ROT = 'rot_'
STATS = '_stats'
BASE_YR_WOUT_URBAN = 'wout_urb_'
#EXISTING_ROTATIONS = 'latest_rot.csv'
HUMAN_READABLE_ROTATIONS = 'human_readable_rotations'+date+'.csv'
PRODUCT = 'product_'
RAND_MIN = 300
ROTATION_STATS = 'rot_stats'+tag+'.csv'
check_rotation = False # If false, then fold each crop rotation with the crop rotation it is most similar to within the max_crop_rot list
remove_urban_and_wtlnd = True #Should be false if you want to create a wall to wall product
                          
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

# Returns true if list a is a rotation of list b    
def checkRotation(a, b):
    if (len(a) <> len(b)):
        return False
    
    tmp_a = str(a).strip('[]')
    tmp_b = str(b+b).strip('[]')
    
    return (tmp_a in tmp_b)

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
    print 'Processing ' + str(yrs_in_rot) + ' rasters located in ' + base_dir
    
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
        dir_path = os.path.dirname(crp_rot_files[i])
        file_name = os.path.basename(crp_rot_files[i])[:-11][-7:]
        # Copy the TIF file to a GRID file to avoid a ArcGIS bug
        print file_name
        gp.CopyRaster_management(crp_rot_files[i], dir_path+os.sep+file_name)
        crp_rot_files[i] = dir_path+os.sep+file_name
        gp.ExtractByAttributes_sa(crp_rot_files[i], cdl_2_str, out_dir+RAS+state+str(i))
        rasters_to_delete.append(out_dir+OUT_RAS+state+str(i))
        rasters_to_delete.append(out_dir+RAS+state+str(i))
    
    base_yr_data  = crp_rot_files[base_yr-first_yr]

    # Remove the urban areas from raster for base_yr
    if remove_urban_and_wtlnd==True:
        print '\tRemoving urban areas for '+os.path.split(base_yr_data)[1]
        urban_str = returnLandUses(analysis_dir+'urban_wetlands.txt')
        gp.SetNull_sa(base_yr_data, base_yr_data, out_dir+BASE_YR_WOUT_URBAN+state, urban_str)
        rasters_to_delete.append(out_dir+BASE_YR_WOUT_URBAN+state)
    else:
        pass
        #gp.SetNull_sa(base_yr_data, base_yr_data, out_dir+BASE_YR_WOUT_URBAN+state, "\"VALUE\" = 0")
        #gp.CopyRaster_management(base_yr_data, out_dir+BASE_YR_WOUT_URBAN+state, "", "", "", "NONE", "NONE", "")

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
    ras_cell_size_x = 56#gp.GetRasterProperties (comb_rasters, 'CELLSIZEX')
    ras_cell_size_y = 56#gp.GetRasterProperties (comb_rasters, 'CELLSIZEY')
    # Area is computed in ha (convert 56m*56m into ha)
    area_raster = ras_cell_size_x*ras_cell_size_y*(0.0001)

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
    cdl_f_str = returnLandUses(analysis_dir+'forests_all.txt')
    gp.ExtractByAttributes_sa(base_yr_data, cdl_f_str, forest_data)
    # Find the grassland/pasture/herbaceus vegetation land in base_yr
    grasslands_data = out_dir+'grass'+str(base_yr)
    print 'Extracting grassland in '+str(base_yr)
    cdl_g_str = returnLandUses(analysis_dir+'grasslands_all.txt')
    gp.ExtractByAttributes_sa(base_yr_data, cdl_g_str, grasslands_data)
    # Find the urban/wetland land in base_yr
    urban_data = out_dir+'urb'+str(base_yr)
    print 'Extracting urban and wetland in '+str(base_yr)
    cdl_u_str = returnLandUses(analysis_dir+'urban_wetlands.txt')
    gp.ExtractByAttributes_sa(base_yr_data, cdl_u_str, urban_data)
    # Retain a maximum of maxCropRotations number of rotations. For this:
    # 1. Remove crop rotations which occupy less than 5% of the area
    # 2. Combine crop rotations which are very similar to each other i.e. by looking at the Levenshtein distance
    # 3. Combine crop rotations which differ from each other temporally e.g. [corn,soybean,corn] and [soybean,corn,soybean] are combined.
    # Simplifications:
    
    # 2. Loop through the file, find all the lines which are offset by less than 25% and combine them
    print 'Simplifying crop rotations'
    # Also create the reclassification table
    recl_1 = open(out_dir+'recl_1.csv', 'w')
    state_stats = open(out_dir+state+STATS+'.csv', 'w') # Contains information on how many pixels each rotation occupies
    state_rots = open(out_dir+ROT+state+'.csv', 'w') # Contains information on how many pixels each rotation occupies
    
    #most_common_lu = cdl_comb[0][yrs_in_rot+2] # Basically the value field from the raster attribute table
    most_pixels = cdl_comb[0][yrs_in_rot]
    total_pixels = 0
    for i in range(len(cdl_comb)):
        total_pixels += cdl_comb[i][yrs_in_rot]
    
    cur_num_rot = 0
    cur_perc_land_in_rot = 0.0
    perc_simplified_land = 0.0
    
    if (check_rotation == False):
        # The first max_crop_rot rotations are selected by default
        cur_max_crop_rot = max_crop_rot
        for i in range(max_crop_rot):
            cur_num_rot += 1
            cur_perc_land_in_rot += float(cdl_comb[i][yrs_in_rot])/float(total_pixels)
            if(cur_perc_land_in_rot>min_area_rot):
                cur_max_crop_rot = i+1
                break
    
        for i in range(cur_max_crop_rot, len(cdl_comb)):
            cur_diff = 1.0
            rot_to_match_with =  0
            for j in range(cur_max_crop_rot):   
                # Compute the levenshtein distance between crop rotation types i and j
                diff_btwn_rot = float(lev(cdl_comb[i][:yrs_in_rot], cdl_comb[j][:yrs_in_rot]))/float(yrs_in_rot)
                # Check if the strings are permutations of each other
                is_permutation = checkPermutations(cdl_comb[i][:yrs_in_rot], cdl_comb[j][:yrs_in_rot])
                if diff_btwn_rot < cur_diff:
                    rot_to_match_with = j
    
                    cur_diff = diff_btwn_rot
       
            cdl_comb[rot_to_match_with][yrs_in_rot] += cdl_comb[i][yrs_in_rot]
            # Set the USED bit of rotation 2 to 0
            cdl_comb[i][yrs_in_rot+1] = 0
            # In the reclassification file, set the reclassification bit to the rotation 1
            writeReClassFile(recl_1, cdl_comb[i][yrs_in_rot+2], cdl_comb[rot_to_match_with][yrs_in_rot+2])
            cdl_comb[i][yrs_in_rot+2] = cdl_comb[rot_to_match_with][yrs_in_rot+2]
            perc_simplified_land += (1.0-cur_diff)*cdl_comb[i][yrs_in_rot]/total_pixels
    else: #check_rotation == true
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
                            if (checkRotation(cdl_comb[j][:yrs_in_rot], cdl_comb[i][:yrs_in_rot]) and \
                                        float(cdl_comb[j][yrs_in_rot])/float(most_pixels) <= min_perc_for_rot):
                                cdl_comb[i][yrs_in_rot] += cdl_comb[j][yrs_in_rot]
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
    
    # In the statename_stats.csv file print the crop rotation ID and the number of pixels it contains
    # In the statename_rots.csv file print the crop rotation ID and the rotation list
    state_stats.write('ID, Num_Pixels, Area_ha\n')
    state_rots.write('ID, Rotation\n')
    for i in range(len(cdl_comb)):
        if(cdl_comb[i][yrs_in_rot+1]>0):
            state_stats.write(str(cdl_comb[i][yrs_in_rot+2])+', '+ str(cdl_comb[i][yrs_in_rot])+\
                              ', '+ str(float(cdl_comb[i][yrs_in_rot])*area_raster)+'\n')
            sync_vec = synchronizeRotation(cdl_comb[i][:yrs_in_rot], yrs_in_rot, first_yr, prev_first_yr)
            state_rots.write(str(cdl_comb[i][yrs_in_rot+2])+', '+str(sync_vec).strip('[]')+'\n')
    state_stats.close()
    state_rots.close()
    
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
    # 1. Iterate through cur_rot and compare each cur_rot with the prev_rot
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
        
    global RAND_MIN
    for k1,v1 in cur_rot.iteritems():
        add_cur_rot = True
        for k2,v2 in prev_rot.iteritems():
            if(k1==k2):
                tmp = k1 + RAND_MIN#random.randrange(RAND_MIN,RAND_MAX)
                lis = tmp, v1
                recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(tmp)+'\n')
                add_more_rot.append(lis)
                add_cur_rot = False
                RAND_MIN += 300
        if add_cur_rot == True:
            lis = k1, v1
            recl_cur_to_prev.write(str(k1)+', '+str(k1)+', '+str(k1)+'\n')
            add_more_rot.append(lis)
    
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
    try:
        gp.ReclassByTable_sa(rot_data, out_dir+'recl_cur_to_prev.csv', "FROM"   , "TO", "VALUE", final_rot_data, "DATA")
    except:
         print gp.GetMessages()
    rot_and_frst = out_dir+FRST+state+'_'+str(yrs_in_rot)+YRS
    rot_frst_grass = out_dir+ROT_FRST_GRASS+state+'_'+str(yrs_in_rot)+YRS
    rot_frst_grass_urb = out_dir+ROT_FRST_GRASS_URB+state+'_'+str(yrs_in_rot)+YRS
    # Append the forest data to the remaining set
    soma_exp = "CON(isnull(" + forest_data + ")," + final_rot_data + "," + forest_data + ")"
    gp.SingleOutputMapAlgebra_sa(soma_exp, rot_and_frst)
    rasters_to_delete.append(forest_data)
    rasters_to_delete.append(grasslands_data)
    rasters_to_delete.append(rot_and_frst)
    rasters_to_delete.append(rot_frst_grass)
    # Append the grasslands data to the remaining set
    soma_exp = "CON(isnull(" + grasslands_data + ")," + rot_and_frst + "," + grasslands_data + ")"
    gp.SingleOutputMapAlgebra_sa(soma_exp, rot_frst_grass)
    # Append the grasslands data to the remaining set
    soma_exp = "CON(isnull(" + urban_data + ")," + rot_frst_grass + "," + urban_data + ")"
    gp.SingleOutputMapAlgebra_sa(soma_exp, rot_frst_grass_urb)
    
    for i in range(len(rasters_to_delete)):
        print 'Deleting raster '+os.path.split(rasters_to_delete[i])[1]
        if DELETE:
            gp.Delete_management(rasters_to_delete[i], "")
        
    return rot_frst_grass_urb


print "Started processing at :", time.strftime("%H:%M:%S", time.localtime())
for subdir, dir_list, files in os.walk(base_dir):
    break

state_file = open(analysis_dir+STATES, 'rb')
lines = state_file.readlines()
if not os.path.exists (output_dir):
    os.makedirs(output_dir)
stat_writer = open(output_dir+os.sep+ROTATION_STATS,'wb')
stat_writer.write('State, # Rot, Start Yr, End Yr, % CDL Match, % Overall Match\n')
stat_writer.flush()
to_merge_files = []
merged_states_ras = output_dir+os.sep+'merged'

for line in lines:
    frst = False
    prev_first_yr = 2008
    first_yr = 2008
    last_yr  = 2010
    base_yr  = 2008 # This is the raster which will be used to fill in the gaps
    num_rot = last_yr-first_yr+1
    state = line.split()[0]
    max_crop_rot = 15
    min_area_rot = 0.70# base_dir contains the input data
    ####################################
    if (state=='nd'):
        max_crop_rot = 30
    elif (state=='sd'):
        max_crop_rot = 20
    elif (state=='nv' or state=='az'):
        min_area_rot = 0.70
    elif (state=='id' or state=='ks' or state=='la' or state=='ms'):
        max_crop_rot = 20
    ####################################
    # out_dir contains the output and intermediate analysis stuff
    out_dir  = output_dir+os.sep+state+os.sep
    if not os.path.exists (out_dir):
        os.makedirs(out_dir)
    state_ras_files = []
    
    cur_time = str(time.strftime("%H:%M:%S",time.localtime())).replace(':','_')
    LOG_FILENAME = output_dir+os.sep+'Log_CropRotations'+cur_time+'.log'
    logging.basicConfig(filename = LOG_FILENAME, level=logging.DEBUG)
    # Logging levels are DEBUG, INFO, WARNING, ERROR, and CRITICAL
    logger = logging
    # Copy the EXISTING_ROTATIONS file to the current output directory
    shutil.copyfile(analysis_dir+EXISTING_ROTATIONS, output_dir+os.sep+EXISTING_ROTATIONS)

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
    
    ras = computeCropRotations(state_ras_files, first_yr, last_yr)
    to_merge_files.append(ras)     

print 'Merging all crop rotation rasters...'
soma_exp = "merge("
for i in range(len(to_merge_files)):
    soma_exp+=to_merge_files[i]
    if i < len(to_merge_files)-1:
        soma_exp+=", "
soma_exp+=")"
#print soma_exp
gp.SingleOutputMapAlgebra_sa(soma_exp, merged_states_ras)

cdl_map_file = csv.reader(open(output_dir+os.sep+CDL_VAL_NAME_MAP, 'r'))
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