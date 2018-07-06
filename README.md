# For cross matching script (Cross_Match_Photom_Cats.py) 

Please note this script is structured to take in SExtractor catalog outputs with the following columns: 
 
  1 MAG_APER		   Fixed aperture magnitude vector	     
  2 MAGERR_APER 	   RMS error vector for fixed aperture mag.  
  3 MAG_AUTO		   Kron-like elliptical aperture magnitude   
  4 MAGERR_AUTO 	   RMS error for AUTO magnitude 	     
  5 XPEAK_IMAGE 	   x-coordinate of the brightest pixel       
  6 YPEAK_IMAGE 	   y-coordinate of the brightest pixel       
  7 X_IMAGE		   Object position along x		     
  8 Y_IMAGE		   Object position along y		     
  9 ALPHA_J2000 	   Right ascension of barycenter (J2000)     
 10 DELTA_J2000 	   Declination of barycenter (J2000)	     
 
## 1. Structure of Catalouges for input

Please insure all catalogs are in this format or edit script to accept alturnate formate 
**field_name +'_' + filter_band1 + '_band_+' + date +'.cat'**
eg. antlia_i_band_170318.cat 
Data can be ob observations or of processing. 

## 2. Set Parameters at start of script

filter_band1 = 'i'
filter_band2 = 'r'
filter_band3 = 'g'
field_name = 'antlia'
date =   1703017 

## 3. Run from commandline or in IDE 

## 4. The full cross matched table will be outputed as a fits file: 
The fits file will be in the following format: 
**crossmatched_+ field_name + _bands + filter_band1 + _ + filter_band2 + _ +filter_band3+ _ + date + .fits**

#  For plotting color color plots of crossmatched cats ( Plot_Color_Color.py) 

## 1. In script change the filter bands 1, 2, 3 to the same input used with Cross_Match_Photom_Cats.py

## 2. Run from Command line for standard r-i vs g-r plot or edit newly created table for speicfic plots 

# For cross matching candidate to catalog

## 1. Edit input params in script and enter the candidate coords and search radius in arc seconds. 

## 2. Run Script and and matches in the three bands with errors will be produced 
