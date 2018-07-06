import numpy as np 
import astropy 
from astropy.table import Table, Column, join 
import matplotlib.pyplot as plt
from astropy.io import ascii 
from astroML.crossmatch import crossmatch_angular 
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.plotting import hist
import csv 
from astropy.io import fits 

### Please enter filter bands exactly how they were entred in Crosd_Match_Photom_Cats.py 

filter_band1 = 'i'
filter_band2 = 'r'
filter_band3 = 'g'

crossed_matched_cat = 'crossmatched_test_Bandsi_r_g_000.fits'
hdul = fits.open(crossed_matched_cat)
data = hdul[1].data

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

print(data_table) 


g_r_mag = data_table['band_' +filter_band3+'_mag'] - data_table['band_' +filter_band2+'_mag']
r_i_mag = data_table['band_' +filter_band2+'_mag'] - data_table['band_' +filter_band1+'_mag']

plt.plot(g_r_mag, r_i_mag, 'r.') 
plt.xlabel('g-r')
plt.ylabel('r-i')
plt.show() 
