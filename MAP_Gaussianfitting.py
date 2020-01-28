import os
import glob

from pathlib import Path
import numpy as np

import matplotlib as mp
from matplotlib import pyplot as plt
from matplotlib import gridspec, rcParams, rc

import astropy.units as u
from astropy.nddata import CCDData
from astropy.io import fits

from skimage.feature import peak_local_max
from ccdproc import trim_image

from skimage.feature import peak_local_max

from ccdproc import Combiner, combine

from photutils import aperture_photometry as AP
from photutils.aperture import CircularAnnulus as CAn
from photutils.aperture import CircularAperture as CAp

from astropy.stats import sigma_clip, gaussian_fwhm_to_sigma
from astropy.modeling.models import Gaussian1D, Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter

from scipy.spatial.distance import euclidean as ds
#%%

def disable_mplkeymaps():
    rc('keymap', 
       fullscreen='',
       home='',
       back='',
       forward='',
       pan='',
       zoom='',
       save='',
       quit='',
       grid='',
       yscale='',
       xscale='',
       all_axes=''
       )
FONTSIZE = 12 # Change it on your computer if you wish.
rcParams.update({'font.size': FONTSIZE})
#%%
def trim1(img):
    
#def trim1(x1=256,x2=256,y1=256,y2=256,img):
#    if (x1 >= x2) or (y1 >= y2):
#        return print("interval error")
    
    
    


    if type(img) == list:
        for fpath in img:
            name = os.path.basename(fpath)
            ccd = trim_image( CCDData.read( fpath),
                            # fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
                            # fits_section = "[1921:2176, 1921:2176]"
                             fits_section = "[1792:2047, 1792:2047]"
                             )
            if "cgan" in str(fpath):
                tpath = Path(toppath/"t_prep_cgan")
                tpath.mkdir(exist_ok=True)
                ccd.write(tpath/name,overwrite=True)
            else:
                tpath = Path(toppath/"t_prep")
                tpath.mkdir(exist_ok=True)
                ccd.write(tpath/name,overwrite=True)
    else:
        name = os.path.basename(img)
        ccd = trim_image( CCDData.read( img),
#                         fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
                         fits_section = "[1792:2047, 1792:2047]"
#                         fits_section = "[1:4096, 1:4096]"
                         )
        if "cgan" in str(img):
            tpath = Path(toppath/"t_prep_cgan")
            tpath.mkdir(exist_ok=True)
            ccd.write(tpath/name,overwrite=True)
        else:
            tpath = Path(toppath/"t_prep")
            tpath.mkdir(exist_ok=True)
            ccd.write(tpath/name,overwrite=True)
            

#    return
#%%
        
def find_dir8( new , oldlist ):
    if ( [new[0]+1, new[1]] in oldlist ):
        return "pos6"
    elif ( [new[0]-1, new[1]] in oldlist ):
        return "pos4"
    elif ( [new[0], new[1]+1] in oldlist ):
        return "pos8"
    elif ( [new[0], new[1]-1] in oldlist ):
        return "pos2"
    elif ( [new[0]+1, new[1]+1] in oldlist ):
        return "pos9"
    elif ( [new[0]+1, new[1]-1] in oldlist ):
        return "pos3"
    elif ( [new[0]-1, new[1]+1] in oldlist ):
        return "pos7"
    elif ( [new[0]-1, new[1]-1] in oldlist ):
        return "pos1"

#%%
        
def sub_NT_1st(peak_pos, img, threshold = 40 ):
    del_ind = []

    for i in range(len(peak_pos)):
        if np.abs(img[ peak_pos[i][0], peak_pos[i][1]]) > threshold:
            del_ind.append(i)
        elif peak_pos[i][0] <= 2 :
            del_ind.append(i)
        elif peak_pos[i][0] >= len(img[0])-2 :
            del_ind.append(i)            
        elif peak_pos[i][1] <= 2 :
            del_ind.append(i)
        elif peak_pos[i][1] >= len(img[1])-2 :
            del_ind.append(i)
    delpeakpos_cgan = []
    for i in range(len(del_ind)):
        delpeakpos_cgan.append(list( peak_pos[ del_ind[i] ] ) )
   
    
    peak_pos = np.delete( peak_pos, del_ind, 0 )
    return peak_pos
#%%
def ext_5_img(img,n=5):
    lx = img.shape[0]
    ly = img.shape[1]
    ext = np.zeros([n*lx,n*ly])
    for i in range(lx):
        for j in range(lx):
            ext[i*n:(i+1)*n,j*n:(j+1)*n] = img[i][j]
        
    
    return ext
#%%
def ext_5_pos(Xs,n=5):
    Xs = Xs*n+2
    return Xs
#%%


#%%
#bt = before trimming
toppath = Path("/home/Jeong/work/it/")
bt_rawpath = Path(toppath/"prep")
bt_cganpath = Path(toppath/"prep_cgan")
rawpath = Path(toppath/"t_prep")
rawerrpath = Path(toppath/'err_raw')
cganpath = Path(toppath/"t_prep_cgan")
cganerrpath = Path(toppath/'err_cgan')
# merg7path = Path(toppath/'5m15s_merged_prep')
# err7path = Path(toppath/'err_5m15s_merged_prep')
#%%
raw_plist = list(rawpath.glob("*.fits"))
cgan_plist = list(cganpath.glob("*.fits"))

raw_plist.sort()
cgan_plist.sort()
#%%
# bt_cgan_plist = list(bt_cganpath.glob("*.fits"))

# bt_cgan_plist.sort()
#%%
dat_cgan = CCDData.read(cgan_plist[0]).data.T
#%%
peak_cgan = peak_local_max(np.abs(dat_cgan),min_distance=2,
                           exclude_border=False
                           
                             # ,num_peaks = 1000
                           ,threshold_abs = 9
                           )
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.imshow(dat_cgan.T,
            vmin=-40,
            vmax = 40,
           origin = 'lower')
ax1.set_title("cgan 45s")

# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_cgan = CAp(peak_cgan,r=2)
cir_cgan.plot(ax1,color='r',lw=1,alpha=0.8)

#%%

peak_cgan_ss = sub_NT_1st(peak_cgan, dat_cgan, 60) 

#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.imshow(dat_cgan.T,
            vmin=-40,
            vmax = 40,
           origin = 'lower')
ax1.set_title("cgan 45s")

# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_cgan = CAp(peak_cgan_ss,r=2)
cir_cgan.plot(ax1,color='r',lw=1,alpha=0.8)
#%%
FWHM_ID = 2

x, y = np.mgrid[:256, :256]
#%%
# peak_cgan_gauss = []
FWHM_cgan_gauss = []
# std_cgan_fixed = []

fitter = LevMarLSQFitter()

for peak_pix in (peak_cgan_ss):
    g_init = Gaussian2D(amplitude = np.abs(dat_cgan[peak_pix[0]][peak_pix[1]]), 
                        x_mean = peak_pix[0], y_mean = peak_pix[1],  
                        x_stddev = FWHM_ID * gaussian_fwhm_to_sigma,
                        y_stddev = FWHM_ID * gaussian_fwhm_to_sigma,
                        bounds={'amplitude': (0, 2*dat_cgan[peak_pix[0]][peak_pix[1]]),
                                'x_mean':(peak_pix[0], peak_pix[0]),
                                'y_mean':(peak_pix[1], peak_pix[1]),
                                # 'x_mean':(peak_pix[0]-1*FWHM_ID, peak_pix[0]+1*FWHM_ID),
                                # 'y_mean':(peak_pix[1]-1*FWHM_ID, peak_pix[1]+1*FWHM_ID),
                                'x_stddev':(0, 10*FWHM_ID),
                                'y_stddev':(0, 10*FWHM_ID)})

    fitted = fitter(g_init,x,y,np.abs(dat_cgan))
    # new_pos = [fitted.x_mean.value,fitted.y_mean.value]
    
    # peak_cgan_gauss.append( new_pos )
    FWHM_cgan_gauss.append( [fitted.x_stddev.value/gaussian_fwhm_to_sigma,
                             fitted.y_stddev.value/gaussian_fwhm_to_sigma] )
    
    # if ds(new_pos,peak_pix) >= 1*FWHM_ID :
    #     peak_cgan_gauss.append( list(peak_pix) )
    # else:

    #     peak_cgan_gauss.append( new_pos )    
    
    
    
    
    # g_init2 = Gaussian2D(amplitude = np.abs(dat_cgan[peak_pix[0]][peak_pix[1]]), 
    #                    x_mean = peak_pix[0], y_mean = peak_pix[1],  
    #                    x_stddev = FWHM_ID * gaussian_fwhm_to_sigma,
    #                    y_stddev = FWHM_ID * gaussian_fwhm_to_sigma,
    #                    bounds={'amplitude': (0, 2*dat_cgan[peak_pix[0]][peak_pix[1]]),
    #                            'x_mean':(peak_pix[0], peak_pix[0]),
    #                            'y_mean':(peak_pix[1], peak_pix[1]),
    #                            'x_stddev':(0, 2*FWHM_ID),
    #                            'y_stddev':(0, 2*FWHM_ID)})

    # unfitted = fitter(g_init2,x,y,np.abs(dat_cgan))


    # std_cgan_fixed.append( [unfitted.x_stddev.value,unfitted.y_stddev.value] )
    
    
# peak_gauss = Column(data=peak_gauss, name='pixel_gauss', dtype=float)
# peak_shift = Column(data=peak_gauss - ID_init['pixel_init'],
#                     name='pixel_shift', dtype=float)
# ID_init["pixel_gauss"] = peak_gauss
# ID_init["pixel_shift"] = peak_gauss - ID_init['pixel_init']
# ID_init.sort('wavelength')
# ID_init.pprint()
# # If you want to save:
# # ID_init.write(DATAPATH/'user_input_identify.csv', 
# #              format='ascii.csv', overwrite=True)
#%%
for i in range(0,754,20):
    print(peak_cgan_ss[i],
          # peak_cgan_gauss[i],
          FWHM_cgan_gauss[i])   

#%%
peak_std_5 = []
j=0
for i in range(len(peak_cgan_ss)):
    if (FWHM_cgan_gauss[i][0] > 5/gaussian_fwhm_to_sigma): 
        print(FWHM_cgan_gauss[i])
        peak_std_5.append(peak_cgan_ss[i])
        j += 1
    elif FWHM_cgan_gauss[i][1] > 5/gaussian_fwhm_to_sigma:
        print(FWHM_cgan_gauss[i])
        peak_std_5.append(peak_cgan_ss[i])
        j += 1

print(j)

#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
# gs = gridspec.GridSpec(1,2)
# ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,0])
# ax1.imshow(dat_cgan.T,
#             vmin=-40,
#             vmax = 40,
#            origin = 'lower')
ax2.imshow(dat_cgan.T,
            vmin=-40,
            vmax = 40,
           origin = 'lower')
# ax1.set_title("cgan 45s")

# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_cgan = CAp(peak_cgan_ss,r=1)
cir_cgan.plot(ax2,color='r',lw=1,alpha=0.8)
cir_cgan_fit = CAp(peak_std_5,r=2)
cir_cgan_fit.plot(ax2,color='b',lw=1,alpha=0.8)
#%%
# peak_cgan_gauss2 = peak_cgan_gauss
# std_cgan_gauss2 = std_cgan_gauss
# std_cgan_fixed2 = std_cgan_fixed
#%%
dat_cgan_5 = ext_5_img(dat_cgan)
peak_cgan_ss_5 = ext_5_pos(peak_cgan_ss)

#%%

fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax1.imshow(dat_cgan.T,
            vmin=-40,
            vmax = 40,
            origin = 'lower')
ax2.imshow(dat_cgan_5.T,
            vmin=-40,
            vmax = 40,
           origin = 'lower')
# ax1.set_title("cgan 45s")

# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_cgan = CAp(peak_cgan_ss,r=1)
cir_cgan.plot(ax1,color='r',lw=1,alpha=0.8)
cir_cgan_fit = CAp(peak_cgan_ss_5,r=2)
cir_cgan_fit.plot(ax2,color='b',lw=1,alpha=0.8)
#%%
aperture = CAp(peak_cgan_ss,r=2)
#%%
AP(dat_cgan.T, aperture)
#%%
IN = AP(dat_cgan.T, aperture)
IN['aperture_sum'].info.format = '%.8g'
print(IN)
#%%
mask = (np.abs(dat_cgan) > 9 ).astype(np.int)
#%%
test = mask*dat_cgan
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.imshow(np.abs(test.T),
            # vmin=-40,
            vmax = 40,
            origin = 'lower')
#%%
dat_cgan[mask]
#%%
A = np.ones((3,3))
B = np.ones((3,3))
A[1:3,2]=5
B[1,1:3]=3
#%%
A*B
