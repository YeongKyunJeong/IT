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
#bt = before trimming
toppath = Path("/home/Jeong/work/it/")
bt_rawpath = Path(toppath/"prep")
bt_cganpath = Path(toppath/"prep_cgan")
rawpath = Path(toppath/"t_prep")
cganpath = Path(toppath/"t_prep_cgan")
#%%
# bt_raw_nlist = glob.glob(f"{bt_rawpath}/*.fits")
# bt_raw_plist = list(bt_rawpath.glob("*.fits"))

# bt_cgan_nlist = glob.glob(f"{bt_cganpath}/*.fits")
# bt_cgan_plist = list(bt_cganpath.glob("*.fits"))

#%% 
#triming and saving
# trim1(bt_cgan_plist)
# #%%
# trim1(bt_raw_plist)
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
ccd_cgan = CCDData.read(cgan_plist[0])
dat_cgan = ccd_cgan.data.T
hdr_cgan = ccd_cgan.header
#%%
ccd_raw = CCDData.read(raw_plist[0])
dat_raw = ccd_raw.data.T
hdr_raw = ccd_raw.header
#%%
raw_plist[0:10]
cgan_plist[0:10]
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,3)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])
ax1.imshow(dat_raw.T,vmin=-100,vmax=100,origin='lower')
ax2.imshow(dat_cgan.T,vmin=-100,vmax=100,origin='lower')
ax3.imshow(dat_raw-dat_cgan,origin='lower')

#%%


peak_raw = peak_local_max(np.abs(dat_raw),min_distance=1
                           ,num_peaks = 1000
                          ,threshold_abs = 12
                          )

peak_cgan = peak_local_max(np.abs(dat_cgan),min_distance=1
                             ,num_peaks = 1000
                           ,threshold_abs = 12
                           # ,indices = False
                           )
# mask_peak_cgan = peak_local_max(np.abs(dat_cgan),min_distance=1
#                             # ,num_peaks=100
#                            ,threshold_abs = 10
#                            ,indices = False
#                            )
#%%

print(len(peak_raw),len(peak_cgan))

#%%
# peak_pos_cgan = dict(x=[],y=[])
# peak_pos_cgan = dict(x=[],y=[])
# for i in peak_cgan:

#     peak_pos_cgan["x"].append(i[0])
#     peak_pos_cgan["y"].append(i[1])
# #%%
# len(peak_pos_cgan["x"])
# len(peak_pos_cgan["y"])
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,3)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])
ax1.imshow(np.abs(dat_raw.T),
           # vmin=-100,
           vmin = 0,
           vmax = 100,
           origin = 'lower')
ax2.imshow(np.abs(dat_cgan.T),
            # vmin=-100,
           vmin = 0,
           vmax = 100,
           # vmin = 0,
           # vmax=105,
           origin = 'lower')
ax3.imshow(dat_raw-dat_cgan)

cir_raw = CAp(peak_raw,r=2)
cir_raw.plot(ax1,color='r',lw=1,alpha=0.8)
cir_cgan = CAp(peak_cgan,r=2)
cir_cgan.plot(ax2,color='r',lw=1,alpha=0.8)

#%%
# peak_cgan[0]
# peak_pos_cgan['x'][0]
# peak_pos_cgan['y'][0]
#%%
# mask_cgan = [np.abs(dat_cgan) > 40]
# #%%
# for i in range(500):
#     if np.abs(dat_cgan[peak_pos_cgan["x"][i],peak_pos_cgan["y"][i]]) < 15 :
#         print(np.abs(dat_cgan[peak_pos_cgan["x"][i],peak_pos_cgan["y"][i]]))
# #%%
# for i in range( len(peak_pos_cgan["x"]) ):
#     if np.abs(dat_cgan[peak_pos_cgan["x"][i],peak_pos_cgan["y"][i]]) > 40:
#         np.delete( peak_cgan, i, 0 )
# #%%
# for i in range(500):
#     if np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1] ]) > 15 :
#         print(np.abs(dat_cgan[peak_pos_cgan["x"][i],peak_pos_cgan["y"][i]]))

# #%%
# len(peak_cgan)

# np.abs(dat_cgan[peak_pos_cgan["x"][54],peak_pos_cgan["y"][54]]) > 40
# peak_cgan = np.delete(peak_cgan,54,0)

#%%

delpeaklist_raw = []

for i in range(len(peak_raw)):
    if np.abs(dat_raw[ peak_raw[i][0], peak_raw[i][1]]) > 40:
        delpeaklist_raw.append(i)

peak_raw = np.delete( peak_raw, delpeaklist_raw, 0 )



delpeaklist_cgan = []

for i in range(len(peak_cgan)):
    if np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1]]) > 40:
        delpeaklist_cgan.append(i)

peak_cgan = np.delete( peak_cgan, delpeaklist_cgan, 0 )


#%%

for i in range(len(peak_cgan)) :
    if np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1] ]) > 40 :
        print(np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1] ]))
        
for i in range(len(peak_raw)) :
    if np.abs(dat_raw[ peak_raw[i][0], peak_raw[i][1] ]) > 40 :
        print(np.abs(dat_raw[ peak_raw[i][0], peak_raw[i][1] ]))
#%%
len(peak_cgan)
len(peak_raw)

        
#%%
# inter_pos_cgan = dict(x=[],y=[])
# inter_pos_cgan = dict(x=[],y=[])
# for i in peak_cgan:

#     inter_pos_cgan["x"].append(i[0])
#     inter_pos_cgan["y"].append(i[1])
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,3)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])
ax1.imshow(np.abs(dat_raw.T),
           # vmin=-100,
           # vmax=100,
            vmin = 0,
            vmax = 40,           
           origin = 'lower')
ax2.imshow(np.abs(dat_cgan.T),
            # vmin=-100,
            # vmax = 100,
            vmin = 0,
            vmax = 40,
           origin = 'lower')
ax3.imshow(dat_raw-dat_cgan)

cir_raw = CAp(peak_raw,r=2)
cir_raw.plot(ax1,color='r',lw=1,alpha=0.8)
cir_cgan = CAp(peak_cgan,r=2)
cir_cgan.plot(ax2,color='r',lw=1,alpha=0.8)