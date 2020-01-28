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

#%%
def trim1(img):
#def trim1(x1=256,x2=256,y1=256,y2=256,img):
#    if (x1 >= x2) or (y1 >= y2):
#        return print("interval error")
    if type(img) == list:
        for fpath in img:
            name = os.path.basename(fpath)
            ccd = trim_image( CCDData.read( fpath),
    #                         fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
                             fits_section = "[1921:2176, 1921:2176]"
                             )
            if "cgan" in str(fpath):
                ccd.write(cganpath/name,overwrite=True)
            else:
                ccd.write(rawpath/name,overwrite=True)
    else:
        name = os.path.basename(img)
        ccd = trim_image( CCDData.read( img),
#                         fits_section = f"[{x1} : {x2}, {y1} : {y2}]"
                         fits_section = "[1921:2176, 1921:2176]"
#                         fits_section = "[1:4096, 1:4096]"
                         )
        if "cgan" in str(img):
            ccd.write(cganpath/name,overwrite=True)
        else:
            ccd.write(rawpath/name,overwrite=True)
        
#    return
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
def merg4_30(img, n = 6):
    
    
    
    merg_path = Path(toppath/f'{(3*n)//4}m{n*45%60:02}s_merged_prep')

    merg_path.mkdir(exist_ok=True)
    
    for j in range( len(img) // n ):
        merg_list = []
        for i in range( n ):
            ccd = CCDData.read(img[i+n*j])
            merg_list.append(ccd)
        
        merg_list = combine(merg_list, method = 'average')

        h = j*n//80
        m = (j*n*45%3600) // 60
        s = (j*n*45%60)
        merg_list.header = CCDData.read(raw_plist[n*j + 1]).header
        merg_list.header.add_history(f"{n} images median combined, {(3*n)//4:02}m {n*45%60:02}s image")
        merg_list.write(merg_path/f"hmi_{(3*n)//4:02}m_{n*45%60:02}s_2017_06_14_{16+h:02}_{m:02}_{s:02}_tai_magnetogram.fits",overwrite=True)

#    return            

#%%
#bt = before trimming
toppath = Path("/home/Jeong/work/it/")
bt_rawpath = Path(toppath/"prep")
bt_cganpath = Path(toppath/"prep_cgan")
rawpath = Path(toppath/"t_prep")
rawerrpath = Path(toppath/'err_raw')
cganpath = Path(toppath/"t_prep_cgan")
cganerrpath = Path(toppath/'err_cgan')
mergedrawpath = Path(toppath/'4m30s_merged_prep')
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
# I did not transposed
# for i in range( len(raw_plist) - 1 ):
#     ccd1 = CCDData.read(raw_plist[i]).data
#     ccd2 = CCDData.read(raw_plist[i+1]).data
    
#     ccderr = fits.PrimaryHDU( (ccd1 - ccd2) )
#     ccderr = fits.HDUList( [ccderr] )
#     ccderr.writeto(rawerrpath/f"err_raw{i+1}.fits",overwrite=True)
#%%

#err1(raw_plist)
#%%
raw_err_plist = list(rawerrpath.glob("*.fits"))
cgan_err_plist = list(cganerrpath.glob("*.fits"))

raw_err_plist.sort()
cgan_err_plist.sort()

#%%
for i in range( len(raw_err_plist) ):
    ccd = CCDData.read(raw_err_plist[i], unit = 'G').data
    print( np.sqrt(np.mean(np.square( ccd ))) )

#%%
for i in range( len(cgan_err_plist) ):
    ccd = CCDData.read(cgan_err_plist[i], unit = 'G').data
    print( np.sqrt(np.mean(np.square( ccd ))) )

#%%
merr_raw = []

for i in range( len(raw_err_plist) ):
    ccd = CCDData.read(raw_err_plist[i], unit = 'G')
    merr_raw.append(ccd)
merr_raw = combine(merr_raw, method = 'median')

merr_raw.header.add_history(f"{len(raw_err_plist)} image(s) median combined error")

merr_raw.write(rawerrpath/"MErr_raw.fits",overwrite=True)
#%%
np.mean(merr_raw.data)    

#%%
merr_cgan = []

for i in range( len(cgan_err_plist) ):
    ccd = CCDData.read(cgan_err_plist[i], unit = 'G')
    merr_cgan.append(ccd)
merr_cgan = combine(merr_cgan, method = 'median')

merr_cgan.header.add_history(f"{len(cgan_err_plist)} image(s) median combined error")

merr_cgan.write(cganerrpath/"MErr_cgan.fits",overwrite=True)
#%%
np.mean(merr_cgan.data)    

#%%
raw_err_plist = list(rawerrpath.glob("*.fits"))
cgan_err_plist = list(cganerrpath.glob("*.fits"))

raw_err_plist.sort()
cgan_err_plist.sort()


#%%

# merged_raw2 = merged_raw
#%%
merged_raw = []

j = 1

for i in range(6):
    ccd = CCDData.read(raw_plist[i+6*j])
    merged_raw.append(ccd)

merged_raw = combine(merged_raw, method = 'average')


merged_raw.header = CCDData.read(raw_plist[6*j + 1]).header
merged_raw.header.add_history("6 images median combined, 4m 30s image")
merged_raw.write(mergedrawpath/"hmi_4m_30s_2017_06_14_16_04_30_tai_magnetogram.fits",overwrite=True)

    # images=[]
    # # Bias Image를 불러 Median Combine
    # for i in range(len(biastab)):
    #     cc = CCDData.read(biastab[i]['FILE'], unit = u.adu)
    #     print(cc.shape)
    #     images.append(cc)
    # mbias = combine(images, method = 'median')
    
    # # 첫 번째 Raw Flat의 Header를 가져와 Master Flat에 넣고 History 추가, 32bit로 바꿔 저장
    # mbias.header = cc.header
    # mbias.header.add_history(f"{len(biastab)} image(s) median combined bias frames")
    # mbias = yfu.CCDData_astype(mbias,dtype = 'float32')
    # mbias.write(prepath/bias_fname,overwrite=True)
    # print(type(mbias))
#%%
merged_raw.header



#%%
ccd1 = CCDData.read(raw_plist[0]).data.T
ccd2 = CCDData.read(raw_plist[1]).data.T
ccderr = ccd1 - ccd2

#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,3)
ax2 = plt.subplot(gs[0,0])
ax1 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])
ax1.imshow(merged_raw.data,
            vmin = -100,
            vmax = 100,
            origin = 'lower')
ax2.imshow(CCDData.read(raw_plist[0]).data,
            vmin = -100,
            vmax = 100,
            origin = 'lower')
ax3.imshow(CCDData.read(cgan_plist[2]).data,
            vmin = -100,
            vmax = 100,
            origin = 'lower')
#%%
np.sqrt(np.mean(np.square( merged_raw.data - merged_raw2.data )))
np.sqrt(np.mean(np.square( CCDData.read(raw_plist[0]).data - CCDData.read(raw_plist[1]).data )))
np.sqrt(np.mean(np.square( CCDData.read(cgan_plist[2]).data - CCDData.read(cgan_plist[3]).data )))

#%%
#Find i range which makes root mean sqaure of merged same with cgan
#%%
# merg4_30(raw_plist,4)
#%%
a = np.asarray([[1,2,3],[6,5,4],[7,8,9]]) 
#%%
a
a**2
