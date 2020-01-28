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
from astropy.wcs import WCS
from astropy.time import TimeDelta

from skimage.feature import peak_local_max
from ccdproc import trim_image

from skimage.feature import peak_local_max

from ccdproc import Combiner, combine

import sunpy.map
from sunpy.physics.differential_rotation import diff_rot, solar_rotate_coordinate
# import sunpy.data.sample


from photutils.aperture import CircularAnnulus as CAn
from photutils.aperture import CircularAperture as CAp

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
# #%%
# def flow(x,y1,y2,dt):
#     dt = dt * u.second
#     y1 *= np.pi/180
#     y2 *= np.pi/180
#     lat = np.arcsin( ( np.sin(y1) + np.sin(y2) )/2 )
#     lat *= 180/np.pi*u.deg
#     rot = [diff_rot( dt*u.sec, (y2+y1)/2) ]
    

#     return
# #%%
# def pix2long(x,data = img[find]):
# # pix to wcs
    

#     return


# help(diff_rot)

# diff_rot(25*u.day,0*u.deg)
#%%
# y = 2
# def test(x = y):
#     return(print(x))
# #%%
# test2 = fits.open(bt_raw_plist[0])
# hdu = test2[0]
# wcss = WCS(hdu.header)


# test = CCDData.read(bt_raw_plist[0],wcs = wcss)
# #%%
# test.heade
# #%%
# WCS(test.header)


# #%%
# ccd = trim_image( test,
#                              fits_section = "[1500:2000,1500:2000]")
# tpath = Path(toppath/"t_prep_cgan")
# tpath.mkdir(exist_ok=True)
# ccd.write(tpath/"test.fits",overwrite=True)
# #%%
# test = CCDData.read(tpath/"test.fits")
# #%%
# WCS(test.header)
# #%%
# WCS(ccd.header)

# #%%
# test2 = fits.open(bt_raw_plist[0])
# hdu = test2[0]
# wcss = WCS(hdu.header)

# wcs=WCS(test.header)
# wcs
     
     
#%%
def trim1(img, x = 30, y = 0 , l =256, m = 2, t = 45):
    
    
    if not m in [1,2,3,4]:
        print("m should be one of among 1,2,3,4")
        return
    
    if m == 2:
        x1,x2 = 2048 - 1 - x, 2048 - l - x
        y1,y2 = 2049 + 1 + y, 2049 + l + y

    # elif m = 3:
    # elif m = 4:
    # elif m = 1:

    if type(img) == list:
        for find in len(img):
            name = os.path.basename(img[find])
            x1 = flow(x1,t*find)
            x2 = flow(x2,t*find)
            ccd = trim_image( CCDData.read( img[find] ),
                             fits_section = f"[{x1}:{x2},{y1},{y2}]"

            
                            # fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
                            # fits_section = "[1921:2176, 1921:2176]"
                             )
            
            if "cgan" in str(img[find]):
                tpath = Path(toppath/"t_prep_cgan")
                tpath.mkdir(exist_ok=True)
                ccd.write(tpath/name,overwrite=True)
            else:
                tpath = Path(toppath/"t_prep")
                tpath.mkdir(exist_ok=True)
                ccd.write(tpath/name,overwrite=True)        
        # for fpath in img:
        #     name = os.path.basename(fpath)
        #     x1 = flow(x1,t*i)
        #     x2 = flow(x2,t*i)
        #     ccd = trim_image( CCDData.read( fpath),
        #                      fits_section = f"[{x1}:{x2},{y1},{y2}]"

            
        #                     # fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
        #                     # fits_section = "[1921:2176, 1921:2176]"
        #                      )
            
        #     if "cgan" in str(fpath):
        #         tpath = Path(toppath/"t_prep_cgan")
        #         tpath.mkdir(exist_ok=True)
        #         ccd.write(tpath/name,overwrite=True)
        #     else:
        #         tpath = Path(toppath/"t_prep")
        #         tpath.mkdir(exist_ok=True)
        #         ccd.write(tpath/name,overwrite=True)
    else:
        name = os.path.basename(img)
        ccd = trim_image( CCDData.read( img),
                        fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
                         # fits_section = "[1792:2047, 1792:2047]"
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
# test = CCDData.read(bt_raw_plist[0])    
# dat = test.data.T
# hdr = test.header
# wcs = WCS(hdr)

# wcs
#%%
ccd= CCDData.rea
            
            
#%%
fopen = fits.open(bt_raw_plist[0])

hdu = fopen[0]
wcs = WCS(hdu.header)
dat = hdu.data.T

#%%    
ccd = trim_image( dat,
                             # fits_section = f"[{x1}:{x2},{y1},{y2}]"

            
                            # fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
                            fits_section = "[1921:2176, 1921:2176]"
                             )
            
#%%
fopen = fits.open(bt_raw_plist[0])

hdu = fopen[0]
wcs = WCS(hdu.header)
dat = hdu.data.T
# hdu = fits.open(bt_raw_plist[0])[0]

#%%
fopen.close()
#%%
dat[1826][1578]
#%%
w1 = wcs.wcs_pix2world([2047.5,2047.5], [2048,2048], 0)

#%%
fig = plt.figure()
fig.add_subplot(111,
                # projection=wcs
                )
plt.imshow(dat.T, origin='lower', vmin=-100,vmax=100 ,cmap=plt.cm.viridis)
plt.xlabel('RA')
plt.ylabel('Dec')
#%

#%%

def err1(img):
    for findex in range(len(img) - 1 ):
        
        ccd1 = CCDData.read(img[findex]).data
        ccd2 = CCDData.read(img[findex + 1]).data
        
        ccderr = fits.PrimaryHDU( (ccd1 - ccd2) )
        ccderr = fits.HDUList( [ccderr] )
        
        if 'cgan' in str( img[findex]) :
            ccderr.writeto(cganerrpath/f"err_cgan{findex+1:03}.fits",overwrite=True)
        else:
            ccderr.writeto(rawerrpath/f"err_raw{findex+1:03}.fits",overwrite=True)
        
        
#    return
#%%
#h = j*n*45 // 3600
#m = (j*n*45%3600) // 60
#s = (j*n*45%60)
def merg5_15(img, n = 7):
    
    
    
    merg_path = Path(toppath/f'{(3*n)//4}m{n*45%60:02}s_merged_prep')
    error_path = Path(toppath/f'err_{(3*n)//4}m{n*45%60:02}s_merged_prep')
    merg_path.mkdir(exist_ok=True)
    error_path.mkdir(exist_ok=True)
    for j in range( len(img) // n ):
        merg_list = []
        error_list = []
        for i in range( n ):
            ccd = CCDData.read(img[i+n*j])
            merg_list.append(ccd)
            error_list.append(np.array(ccd))
        merg_list = combine(merg_list, method = 'average')
        error_list = np.std(error_list,0)
        h = j*n//80
        m = (j*n*45%3600) // 60
        s = (j*n*45%60)
        merg_list.header = CCDData.read(raw_plist[n*j + 1]).header
        merg_list.header.add_history(f"{n} images median combined, {(3*n)//4:02}m {n*45%60:02}s image")
        merg_list.write(merg_path/f"hmi_{(3*n)//4:02}m_{n*45%60:02}s_2017_06_14_{16+h:02}_{m:02}_{s:02}_tai_magnetogram.fits",overwrite=True)
        error_list = fits.PrimaryHDU(error_list)
        error_list = fits.HDUList([error_list])
        error_list.writeto(error_path/f"err_hmi_{(3*n)//4:02}m_{n*45%60:02}s_2017_06_14_{16+h:02}_{m:02}_{s:02}_tai_magnetogram.fits",overwrite=True)
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
    
    delpeakpos_cgan = []
    for i in range(len(del_ind)):
        delpeakpos_cgan.append(list( peak_pos[ del_ind[i] ] ) )
    
    
    peak_pos = np.delete( peak_pos, del_ind, 0 )
    return peak_pos


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
# bt_raw_nlist = glob.glob(f"{bt_rawpath}/*.fits")
bt_raw_plist = list(bt_rawpath.glob("*.fits"))

bt_raw_plist.sort()
# bt_cgan_nlist = glob.glob(f"{bt_cganpath}/*.fits")
bt_cgan_plist = list(bt_cganpath.glob("*.fits"))

bt_cgan_plist.sort()

# trim1(bt_raw_plist)
# trim1(bt_cgan_plist)
#%%
raw_nlist = glob.glob(f"{rawpath}/*.fits")
raw_plist = list(rawpath.glob("*.fits"))

cgan_nlist = glob.glob(f"{cganpath}/*.fits")
cgan_plist = list(cganpath.glob("*.fits"))

raw_nlist.sort()
raw_plist.sort()
cgan_nlist.sort()
cgan_plist.sort()
#%%
merg_7_plist = list(merg7path.glob("*.fits"))
err_7_plist = list(err7path.glob("*.fits"))
merg_7_plist.sort()
err_7_plist.sort()
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax1.imshow(dat_7mer.T,
            vmin=-40,
            vmax=40,          
            origin = 'lower')
ax2.imshow(dat_7err.T,
            vmin= 10,
            vmax = 15,
            origin = 'lower')
# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)
