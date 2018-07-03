import numpy as np 
import astropy 
from astropy.table import Table, Column, join 
import matplotlib.pyplot as plt
from astropy.io import ascii 
from astroML.crossmatch import crossmatch_angular 
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.plotting import hist
import csv 


######---------------- Variables-----------------#############



filter_band1 = 'i'
filter_band2 = 'r'
filter_band3 = 'g'
field_name = 'test'
date =  '000'






#--------------Import the obs cataloug in first filter band----------------------
MAG_APP_OBS_band1, MAGERR_APP_OBS_band1, MAG_AUTO_OBS_band1, MAGERR_AUTO_OBS_band1, XPEAK_OBS_band1, YPEAK_OBS_band1, X_IMG_OBS_band1, Y_IMG_OBS_band1, RA_OBS_band1, DEC_OBS_band1 = np.loadtxt(field_name +'_' + filter_band1 + '_band_' + date +'.cat', unpack = True) 

Coords_obs_band1 = np.empty((len(MAG_APP_OBS_band1), 4), dtype = np.float64)
Coords_obs_band1[:, 0] = RA_OBS_band1
Coords_obs_band1[:, 1] = DEC_OBS_band1
Coords_obs_band1[:, 2] = MAG_AUTO_OBS_band1
Coords_obs_band1[:, 3] = MAGERR_AUTO_OBS_band1


#--------------Import the obs cataloug in 2second filter band----------------------
MAG_APP_OBS_band2, MAGERR_APP_OBS_band2, MAG_AUTO_OBS_band2, MAGERR_AUTO_OBS_band2, XPEAK_OBS_band2, YPEAK_OBS_band2, X_IMG_OBS_band2, Y_IMG_OBS_band2, RA_OBS_band2, DEC_OBS_band2 = np.loadtxt(field_name +'_' + filter_band2 + '_band_' + date +'.cat', unpack = True) 

Coords_obs_band2 = np.empty((len(MAG_APP_OBS_band2), 4), dtype = np.float64)
Coords_obs_band2[:, 0] = RA_OBS_band2
Coords_obs_band2[:, 1] = DEC_OBS_band2
Coords_obs_band2[:, 2] = MAG_AUTO_OBS_band2
Coords_obs_band2[:, 3] = MAGERR_AUTO_OBS_band2


#--------------Import the obs cataloug in 3rd filter band----------------------
MAG_APP_OBS_band3, MAGERR_APP_OBS_band3, MAG_AUTO_OBS_band3, MAGERR_AUTO_OBS_band3, XPEAK_OBS_band3, YPEAK_OBS_band3, X_IMG_OBS_band3, Y_IMG_OBS_band3, RA_OBS_band3, DEC_OBS_band3 = np.loadtxt(field_name +'_' + filter_band3 + '_band_' + date +'.cat', unpack = True) 

Coords_obs_band3 = np.empty((len(MAG_APP_OBS_band3), 4), dtype = np.float64)
Coords_obs_band3[:, 0] = RA_OBS_band3
Coords_obs_band3[:, 1] = DEC_OBS_band3
Coords_obs_band3[:, 2] = MAG_AUTO_OBS_band3
Coords_obs_band3[:, 3] = MAGERR_AUTO_OBS_band3

print(Coords_obs_band3)


#---------------Cross match between sources between the two bands -------------- #

max_radius = 1./3600 # 1 arcsec

dist_between, ind_row = crossmatch_angular(Coords_obs_band1, Coords_obs_band2, max_radius)
match = ~np.isinf(dist_between)
#print(match)


match_table = Table() 
match_table['band2_match_TF'] = match 
match_table['band2_match_ID'] = ind_row
match_table['band1_obs_mag'] = MAG_AUTO_OBS_band1
match_table['band1_obs_mags_err'] = MAGERR_AUTO_OBS_band1
match_table['band1_obs_ra'] = RA_OBS_band1
match_table['band1_obs_dec'] = DEC_OBS_band1

band2_match_true = []
band2_row_matched = []
band1_band1_mags_matched = []
band1_obs_mags_error_matched = []
band1_band1_obs_ra_matched = []
band1_band1_obs_dec_matched = []

for row in match_table: 
	if row['band2_match_TF'] == True:
		band2_match_true.append(row['band2_match_TF'])
		band2_row_matched.append(row['band2_match_ID'])
		band1_band1_mags_matched.append(row['band1_obs_mag'])
		band1_obs_mags_error_matched.append(row['band1_obs_mags_err'])
		band1_band1_obs_ra_matched.append(row['band1_obs_ra'])
		band1_band1_obs_dec_matched.append(row['band1_obs_dec'])


band2_RA = []
band2_DEC = []
band2_MAG = []
band2_MAGerr = []


for i in band2_row_matched: 
	RA = Coords_obs_band2[i, 0]
	DEC = Coords_obs_band2[i, 1]
	MAG =Coords_obs_band2[i, 2]
	MAGERR = Coords_obs_band2[i, 3]
	band2_RA.append(RA)
	band2_DEC.append(DEC)
	band2_MAG.append(MAG)
	band2_MAGerr.append(MAGERR)
	

source_match_table = Table()
source_match_table['band2_match_true'] = band2_match_true
source_match_table['band2_match_rowID'] = band2_row_matched
source_match_table['band2_mag'] = band2_MAG
source_match_table['band2_mag_error'] = band2_MAGerr
source_match_table['band2_ra'] = band2_RA
source_match_table['band2_dec'] = band2_DEC
source_match_table['band1_mags_matched'] = band1_band1_mags_matched
source_match_table['band1_mag_error'] = band1_obs_mags_error_matched
source_match_table['band1_obs_ra_matched'] = band1_band1_obs_ra_matched
source_match_table['band1_obs_dec_matched'] = band1_band1_obs_dec_matched	



###---------------------------------- clean out unreadable magnitudes (anything greather then 30 ------------------------------###########

clean_band2_match_true = []
clean_band2_row_matched = []
clean_band2_MAG = []
clean_band2_MAGerr = []
clean_band2_RA = []
clean_band2_DEC = []
clean_band1_band1_mags_matched = []
clean_band1_obs_mags_error_matched = []
clean_band1_band1_obs_ra_matched = []
clean_band1_band1_obs_dec_matched = []


bad_band2_mag =[]
bad_band2_magerr =[]
bad_band2_ra =[]
bad_band2_dec =[]
bad_band1_mag =[]
bad_band1_magerr =[]
bad_band1_ra =[]
bad_band1_dec =[]


for row in source_match_table:
	if row['band1_mags_matched'] < 30:
		clean_band2_match_true.append(row['band2_match_true'])
		clean_band2_row_matched.append(row['band2_match_rowID'])
		clean_band2_MAG.append(row['band2_mag'])
		clean_band2_MAGerr.append(row['band2_mag_error'])
		clean_band2_RA.append(row['band2_ra'])
		clean_band2_DEC.append(row['band2_dec'])
		clean_band1_band1_mags_matched.append(row['band1_mags_matched'])
		clean_band1_obs_mags_error_matched.append(row['band1_mag_error'])
		clean_band1_band1_obs_ra_matched.append(row['band1_obs_ra_matched'])
		clean_band1_band1_obs_dec_matched.append(row['band1_obs_dec_matched'])
		
	else: 
		bad_band2_mag.append(row['band2_mag'])
		bad_band2_magerr.append(row['band2_mag_error'])
		bad_band2_ra.append(row['band2_ra'])
		bad_band2_dec.append(row['band2_dec'])
		bad_band1_mag.append(row['band1_mags_matched'])
		bad_band1_magerr.append(row['band1_mag_error'])
		bad_band1_ra.append(row['band1_obs_ra_matched'])
		bad_band1_dec.append(row['band1_obs_dec_matched'])



b1_cleaned_table = Table()
b1_cleaned_table['band2_match_true'] = clean_band2_match_true
b1_cleaned_table['band2_match_rowID'] = clean_band2_row_matched
b1_cleaned_table['band2_mag'] = clean_band2_MAG
b1_cleaned_table['band2_mag_error'] = clean_band2_MAGerr
b1_cleaned_table['band2_ra'] = clean_band2_RA
b1_cleaned_table['band2_dec'] = clean_band2_DEC
b1_cleaned_table['band1_mags_matched'] = clean_band1_band1_mags_matched
b1_cleaned_table['band1_mag_error'] = clean_band1_obs_mags_error_matched
b1_cleaned_table['band1_obs_ra_matched'] = clean_band1_band1_obs_ra_matched
b1_cleaned_table['band1_obs_dec_matched'] = clean_band1_band1_obs_dec_matched



clean2_band2_match_true = []
clean2_band2_row_matched = []
clean2_band2_MAG = []
clean2_band2_MAGerr = []
clean2_band2_RA = []
clean2_band2_DEC = []
clean2_band1_band1_mags_matched = []
clean2_band1_obs_mags_error_matched = []
clean2_band1_band1_obs_ra_matched = []
clean2_band1_band1_obs_dec_matched = []

for row in b1_cleaned_table:
	if row['band2_mag'] < 30:
		clean2_band2_match_true.append(row['band2_match_true'])
		clean2_band2_row_matched.append(row['band2_match_rowID'])
		clean2_band2_MAG.append(row['band2_mag'])
		clean2_band2_MAGerr.append(row['band2_mag_error'])
		clean2_band2_RA.append(row['band2_ra'])
		clean2_band2_DEC.append(row['band2_dec'])
		clean2_band1_band1_mags_matched.append(row['band1_mags_matched'])
		clean2_band1_obs_mags_error_matched.append(row['band1_mag_error'])
		clean2_band1_band1_obs_ra_matched.append(row['band1_obs_ra_matched'])
		clean2_band1_band1_obs_dec_matched.append(row['band1_obs_dec_matched'])
		
	else: 
		bad_band2_mag.append(row['band2_mag'])
		bad_band2_magerr.append(row['band2_mag_error'])
		bad_band2_ra.append(row['band2_ra'])
		bad_band2_dec.append(row['band2_dec'])
		bad_band1_mag.append(row['band1_mags_matched'])
		bad_band1_magerr.append(row['band1_mag_error'])
		bad_band1_ra.append(row['band1_obs_ra_matched'])
		bad_band1_dec.append(row['band1_obs_dec_matched'])	
	
final_cleaned_matches = Table()
final_cleaned_matches['band2_match_true'] = clean2_band2_match_true
final_cleaned_matches['band2_match_rowID'] = clean2_band2_row_matched
final_cleaned_matches['band2_mag'] = clean2_band2_MAG
final_cleaned_matches['band2_mag_error'] = clean2_band2_MAGerr
final_cleaned_matches['band2_ra'] = clean2_band2_RA
final_cleaned_matches['band2_dec'] = clean2_band2_DEC
final_cleaned_matches['band1_mags_matched'] = clean2_band1_band1_mags_matched
final_cleaned_matches['band1_mag_error'] = clean2_band1_obs_mags_error_matched
final_cleaned_matches['band1_obs_ra_matched'] = clean2_band1_band1_obs_ra_matched
final_cleaned_matches['band1_obs_dec_matched'] = clean2_band1_band1_obs_dec_matched



####------------------------------------ Create final table of band 1 and band 2 -------------------------------------------------- #########



match_i_r_table = Table()
match_i_r_table['band1_mag'] = clean2_band1_band1_mags_matched
match_i_r_table['band1_mag_err'] = clean2_band1_obs_mags_error_matched
match_i_r_table['band1_ra'] = clean2_band1_band1_obs_ra_matched
match_i_r_table['band1_dec'] = clean2_band1_band1_obs_dec_matched
match_i_r_table['band2_mag'] = clean2_band2_MAG
match_i_r_table['band2_mag_err'] = clean2_band2_MAGerr
match_i_r_table['band2_ra'] = clean2_band2_RA
match_i_r_table['band2_dec'] = clean2_band2_DEC

matches_2_band = np.empty((len(match_i_r_table['band2_ra']), 2), dtype = np.float64)
matches_2_band[:, 0] = match_i_r_table['band2_ra']
matches_2_band[:, 1] = match_i_r_table['band2_dec']
#matches_2_band[:, 2] = match_i_r_table['band2_mag']
#matches_2_band[:, 3] = match_i_r_table['band2_mag_err']

#print(matches_2_band)
#print(len(matches_2_band))


#plt.plot(match_i_r_table['band2_ra'], match_i_r_table['band1_ra'], 'bo')
#plt.show()

#####----------------------------------- cross match to band 3 ---------------------------------------------------


max_radius_2 = 1./3600 # 1 arcsec

dist_between_2, ind_row_2 = crossmatch_angular( matches_2_band, Coords_obs_band3,max_radius_2)
match2 = ~np.isinf(dist_between_2)
print(match2)

match2_table = Table() 
match2_table['band3_match_TF'] = match2 
match2_table['band3_match_ID'] = ind_row_2
match2_table['band1_obs_mag'] = match_i_r_table['band1_mag']
match2_table['band1_obs_mags_err'] = match_i_r_table['band1_mag_err']
match2_table['band1_obs_ra'] = match_i_r_table['band1_ra']
match2_table['band1_obs_dec'] = match_i_r_table['band1_dec']
match2_table['band2_obs_mag'] = match_i_r_table['band2_mag']
match2_table['band2_obs_mags_err'] = match_i_r_table['band2_mag_err']
match2_table['band2_obs_ra'] = match_i_r_table['band2_ra']
match2_table['band2_obs_dec'] = match_i_r_table['band2_dec']

#print(match2_table)
#print(len(match2_table))

all_true_band1_mag = []
all_true_band1_mag_err = []
all_true_band1_ra = []
all_true_band1_dec = []

all_true_band2_mag = []
all_true_band2_mag_err = []
all_true_band2_ra = []
all_true_band2_dec = []

all_true_band3_id = []
all_true_band3_mag = []
all_true_band3_mag_err = []
all_true_band3_ra = []
all_true_band3_dec = []


for row in match2_table: 
	if row['band3_match_TF'] == True:
		all_true_band3_id.append(row['band3_match_ID'])
		all_true_band1_mag.append(row['band1_obs_mag'])
		all_true_band1_mag_err.append(row['band1_obs_mags_err'])
		all_true_band1_ra.append(row['band1_obs_ra'])
		all_true_band1_dec.append(row['band1_obs_dec'])
		all_true_band2_mag.append(row['band2_obs_mag'])
		all_true_band2_mag_err.append(row['band2_obs_mags_err'])
		all_true_band2_ra.append(row['band2_obs_ra'])
		all_true_band2_dec.append(row['band2_obs_dec'])

band3_RA = []
band3_DEC = []
band3_MAG = []
band3_MAGerr = []

#print(all_true_band1_mag)

for i in all_true_band3_id: 
	RA_3 = Coords_obs_band3[i, 0]
	DEC_3 = Coords_obs_band3[i, 1]
	MAG_3 =Coords_obs_band3[i, 2]
	MAGERR_3 = Coords_obs_band3[i, 3]
	band3_RA.append(RA_3)
	band3_DEC.append(DEC_3)
	band3_MAG.append(MAG_3)
	band3_MAGerr.append(MAGERR_3)
	

###---------------------------------------- create final table with all band color information ----------------------------#########

final_table = Table()
final_table['band1_mag'] = all_true_band1_mag
final_table['band1_mag_err'] = all_true_band1_mag_err
final_table['band1_ra'] = all_true_band1_ra
final_table['band1_dec'] = all_true_band1_dec
final_table['band2_mag'] = all_true_band2_mag
final_table['band2_mag_err'] = all_true_band2_mag_err
final_table['band2_ra'] = all_true_band2_ra
final_table['band2_dec'] = all_true_band2_dec
final_table['band3_mag'] = band3_MAG
final_table['band3_mag_err'] = band3_MAGerr
final_table['band3_ra'] = band3_RA
final_table['band3_dec'] = band3_DEC

print(final_table)

outputfile = 'crossmatched_'+ field_name+ '_'+ filter_band1 + '_' + filter_band2 + '_' +filter_band3+ '_' + date + '.csv'

with open(outputfile, 'w') as output:
	writer = csv.writer(output, lineterminator = '\n')
	for  val in final_table:
		writer.writerow([val])

