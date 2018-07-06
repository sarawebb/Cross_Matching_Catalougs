import numpy as np 
import astropy 
from astropy.table import Table, Column, join 
import matplotlib.pyplot as plt
from astropy.io import ascii 
from astroML.crossmatch import crossmatch_angular 
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.plotting import hist
import csv 
import sys 
from astropy.io import fits 

#---------------------- Enter values below------------------------# 


#--- RA and Dec of the position of canidate you want to crossmamtch---- 
RA = 157.9972
DEC = -36.52447

#--- Please enter the filter bands in order used in previous scripts! 
filter_band1 = 'i'
filter_band2 = 'r'
filter_band3 = 'g'


#--- the seperation of the search area to cross match with 
sep_arcsec = 5 







can_entry =  np.empty((1, 2), dtype = np.float64)
can_entry[:,0] = RA
can_entry[:,1] = DEC

print(can_entry)








crossed_matched_cat = 'crossmatched_test_Bandsi_r_g_000.fits'
hdul = fits.open(crossed_matched_cat)
data = hdul[1].data
print(data)



data_table = Table()
data_table['band_' +filter_band1+'_mag'] = data['band_' +filter_band1+'_mag']
data_table['band_' +filter_band1+'_mag_err'] = data['band_' +filter_band1+'_mag_err']
data_table['band_' +filter_band1+'_ra'] = data['band_' +filter_band1+'_ra']
data_table['band_' +filter_band1+'_dec'] = data['band_' +filter_band1+'_dec']

data_table['band_' +filter_band2+'_mag'] = data['band_' +filter_band2+'_mag']
data_table['band_' +filter_band2+'_mag_err'] = data['band_' +filter_band2+'_mag_err']
data_table['band_' +filter_band2+'_ra'] = data['band_' +filter_band2+'_ra']
data_table['band_' +filter_band2+'_dec'] = data['band_' +filter_band2+'_dec']

data_table['band_' +filter_band3+'_mag'] = data['band_' +filter_band3+'_mag']
data_table['band_' +filter_band3+'_mag_err'] = data['band_' +filter_band3+'_mag_err']
data_table['band_' +filter_band3+'_ra'] = data['band_' +filter_band3+'_ra']
data_table['band_' +filter_band3+'_dec'] = data['band_' +filter_band3+'_dec']


#----------- make coords_cat for the the first band to cross match with the entred coords ----------------# 
coords_cat =  np.empty((len(data), 12), dtype = np.float64)
coords_cat[:,0] = data['band_' +filter_band1+'_ra']
coords_cat[:,1] = data['band_' +filter_band1+'_dec']
coords_cat[:,2] = data['band_' +filter_band1+'_mag']





max_radius = sep_arcsec/3600. # 5 arcsec






dist_between_can, ind_row_can = crossmatch_angular(coords_cat, can_entry, max_radius)
match_can = ~np.isinf(dist_between_can)
print(match_can)

coords_matched_band1_mag = []
coords_matched_band1_mag_err = []
coords_matched_band2_mag = []
coords_matched_band2_mag_err = []
coords_matched_band3_mag = []
coords_matched_band3_mag_err = []

coords_matched_table = Table()
coords_matched_table['match_true_false'] = match_can
coords_matched_table['band_' +filter_band1+'_mag'] = coords_cat[:, 2]
coords_matched_table['band_' +filter_band1+'_mag_err'] = data_table['band_' +filter_band1+'_mag_err']
coords_matched_table['band_' +filter_band1+'_ra'] = coords_cat[:,0]
coords_matched_table['band_' +filter_band1+'_dec'] = coords_cat[:,1]

coords_matched_table['band_' +filter_band2+'_mag'] =  data_table['band_' +filter_band2+'_mag']
coords_matched_table['band_' +filter_band2+'_mag_err'] =  data_table['band_' +filter_band2+'_mag_err']
coords_matched_table['band_' +filter_band2+'_ra'] = data_table['band_' +filter_band2+'_ra']
coords_matched_table['band_' +filter_band2+'_dec'] = data_table['band_' +filter_band2+'_dec']

coords_matched_table['band_' +filter_band3+'_mag'] = data_table['band_' +filter_band3+'_mag']
coords_matched_table['band_' +filter_band3+'_mag_err'] =  data_table['band_' +filter_band3+'_mag_err']
coords_matched_table['band_' +filter_band3+'_ra'] = data_table['band_' +filter_band3+'_ra']
coords_matched_table['band_' +filter_band3+'_dec'] = data_table['band_' +filter_band3+'_dec']



for row in coords_matched_table: 
	if row['match_true_false'] == True: 
		print('test')
		coords_matched_band1_mag.append(row['band_' +filter_band1+'_mag'])
		coords_matched_band1_mag_err.append(row['band_' +filter_band1+'_mag_err'])
		coords_matched_band2_mag.append(row['band_' +filter_band2+'_mag'])
		coords_matched_band2_mag_err.append(row['band_' +filter_band2+'_mag_err'])
		coords_matched_band3_mag.append(row['band_' +filter_band3+'_mag'])
		coords_matched_band3_mag_err.append(row['band_' +filter_band3+'_mag_err'])


print('Cand coords ' + filter_band1 + 'band magnitude is ', coords_matched_band1_mag, 'error', coords_matched_band1_mag_err) 
print('Cand coords ' + filter_band2 + 'band magnitude is ', coords_matched_band2_mag, 'error ',  coords_matched_band2_mag_err)
print('Cand coords ' + filter_band3 + 'band magnitude is ', coords_matched_band3_mag, 'error ' , coords_matched_band3_mag_err)


