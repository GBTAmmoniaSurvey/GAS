import GAS_first_look as GAS_FL
import numpy as np
# 
# OrionA
#
a_rms = [ 0, 136, 300, 410, 515, 680]
b_rms = [78, 198, 340, 450, 610, 795]
index_rms=GAS_FL.create_index( a_rms, b_rms)
index_peak=np.arange(335,420)

file_in='OrionA/OrionA_0102_NH3_11_all.fits'
file_out=file_in.replace('.fits','_base1.fits')
file_new=GAS_FL.GAS_baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
GAS_FL.GAS_peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

# 
# NGC 1333
#
a_rms = [  0, 165, 239, 325, 435, 540, 700]
b_rms = [100, 200, 270, 370, 475, 630, 795]
index_rms=GAS_FL.create_index( a_rms, b_rms)
index_peak=np.arange(370,450)

file_in='NGC1333/NGC1333_ABCDEFGH_NH3_11_all.fits'
file_out=file_in.replace('.fits','_base1.fits')
file_new=GAS_FL.GAS_baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
GAS_FL.GAS_peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

# 
# Taurus 01-02
#
a_rms = [  0, 145, 230, 310, 420, 520, 690]
b_rms = [100, 205, 270, 375, 475, 640, 795]
index_rms=GAS_FL.create_index( a_rms, b_rms)
index_peak=np.arange(370,415)

file_in='Taurus/Taurus_12_NH3_11_all.fits'
file_out=file_in.replace('.fits','_base1.fits')
file_new=GAS_FL.GAS_baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
GAS_FL.GAS_peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

# 
# Oph L1689 01-02
#
a_rms = [  0, 145, 230, 310, 420, 525, 690]
b_rms = [ 95, 210, 265, 360, 470, 640, 795]
index_rms=GAS_FL.create_index( a_rms, b_rms)
index_peak=np.arange(370,415)

file_in='Oph_L1698/L1698_12_NH3_11_all.fits'
file_out=file_in.replace('.fits','_base1.fits')
file_new=GAS_FL.GAS_baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
GAS_FL.GAS_peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)


# 
# Serpens_Main NH3 (1,1)
#
a_rms = [  0, 255, 575, 825, 1000, 1300, 1619]
b_rms = [140, 455, 640, 860, 1100, 1490, 1750]
index_rms=GAS_FL.create_index( a_rms, b_rms)
index_peak=np.arange(600,800)

file_in='Serpens_Main/Serpens_Main_NH3_11.fits'
file_out=file_in.replace('.fits','_base1.fits')
file_new=GAS_FL.GAS_baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
GAS_FL.GAS_peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)


# 
# L1544 NH3 (1,1)
#
a_rms = [  0,  55, 135, 185, 235, 310]
b_rms = [ 40, 115, 165, 215, 290, 350]
index_rms=GAS_FL.create_index( a_rms, b_rms)
index_peak=np.arange(168,180)

file_in='L1544/L1544_NH3_11.fits'
file_out=file_in.replace('.fits','_base1.fits')
file_new=GAS_FL.GAS_baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
GAS_FL.GAS_peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

# 
# L1544 NH3 (2,2)
#
a_rms = [  0,  20]
b_rms = [ 10,  25]
index_rms=GAS_FL.create_index( a_rms, b_rms)
index_peak=np.arange(10,18)

file_in='L1544/L1544_NH3_22.fits'
file_out=file_in.replace('.fits','_base1.fits')
file_new=GAS_FL.GAS_baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
GAS_FL.GAS_peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
