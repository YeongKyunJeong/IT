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
    
    
    
#     for i in range(len(peak_pos)):
#         y = find_dir8( list(peak_pos[i]) , delpeakpos_cgan )
#         if y == "pos6":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]+1 ][ peak_pos[i][1] ] > 0:
#                 del_ind.append(i)
#         elif y == "pos4":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]-1 ][ peak_pos[i][1] ] > 0:
#                 del_ind.append(i)
#         elif y == "pos8":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[peak_pos[i][0] ][ peak_pos[i][1]+1 ] > 0:
#                 del_ind.append(i)
#         elif y == "pos2":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0] ][ peak_pos[i][1]-1 ] > 0:
#                 del_ind.append(i)
#         elif y == "pos9":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]+1 ][ peak_pos[i][1]+1 ] > 0:
#                 del_ind.append(i)
#         elif y == "pos7":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]-1 ][ peak_pos[i][1]+1 ] > 0:
#                 del_ind.append(i)
#         elif y == "pos3":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]+1 ][ peak_pos[i][1]-1 ] > 0:
#                 del_ind.append(i)
#         elif y == "pos1":
#             if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]-1 ][ peak_pos[i][1]-1 ] > 0:
#                 del_ind.append(i)
     
    
    peak_pos = np.delete( peak_pos, del_ind, 0 )
    return peak_pos
# #%%
        
# def sub_NT(peak_pos, img, threshold = 40 ):
#     del_ind = []

#     for i in range(len(peak_pos)):
#         if np.abs(img[ peak_pos[i][0], peak_pos[i][1]]) > threshold:
#             del_ind.append(i)
    
#     delpeakpos_cgan = []
#     for i in range(len(del_ind)):
#         delpeakpos_cgan.append(list( peak_pos[ del_ind[i] ] ) )
    
#     indicator = []
    
#     while indicator != list( np.unique(del_ind) ):
#         for i in range(len(peak_pos)):
#             indicator = list( np.unique(del_ind) )
#             y = find_dir8( list(peak_pos[i]) , delpeakpos_cgan )
            
#             if y == "pos6":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]+1 ][ peak_pos[i][1] ] > 0:
#                     del_ind.append(i)
#             elif y == "pos4":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]-1 ][ peak_pos[i][1] ] > 0:
#                     del_ind.append(i)
#             elif y == "pos8":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[peak_pos[i][0] ][ peak_pos[i][1]+1 ] > 0:
#                     del_ind.append(i)
#             elif y == "pos2":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0] ][ peak_pos[i][1]-1 ] > 0:
#                     del_ind.append(i)
#             elif y == "pos9":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]+1 ][ peak_pos[i][1]+1 ] > 0:
#                     del_ind.append(i)
#             elif y == "pos7":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]-1 ][ peak_pos[i][1]+1 ] > 0:
#                     del_ind.append(i)
#             elif y == "pos3":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]+1 ][ peak_pos[i][1]-1 ] > 0:
#                     del_ind.append(i)
#             elif y == "pos1":
#                 if  img[peak_pos[i][0]][peak_pos[i][1]] * img[ peak_pos[i][0]-1 ][ peak_pos[i][1]-1 ] > 0:
#                     del_ind.append(i)
     


    
#     peak_pos = np.delete( peak_pos, del_ind, 0 )
#     return peak_pos
#%%
# test = np.array([[0 for i in range(7)] for j in range(7)])
# #%%
# test[0][0]=4
# test[3][3]=8
# test[4][0]=4
# test[1][2]=0
# test[4][3]=6
# test[1][1]=5
# test[2][2]=6
# test[5][5]=7
# test
# #%%
# peak_test = peak_local_max(test,min_distance=1,exclude_border=True
#                              # ,num_peaks = 1000
#                            ,threshold_abs = 3
#                            ,indices=False
#                            )
# peak_test
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
merg7path = Path(toppath/'5m15s_merged_prep')
err7path = Path(toppath/'err_5m15s_merged_prep')
# merg8path = Path(toppath/'6m00s_merged_prep')
# merg4path = Path(toppath/'3m00s_merged_prep')
#%%
# bt_raw_nlist = glob.glob(f"{bt_rawpath}/*.fits")
# bt_raw_plist = list(bt_rawpath.glob("*.fits"))

# bt_cgan_nlist = glob.glob(f"{bt_cganpath}/*.fits")
# bt_cgan_plist = list(bt_cganpath.glob("*.fits"))

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
merg5_15(raw_plist)



#%%
# raw_err_plist = list(rawerrpath.glob("*.fits"))
# cgan_err_plist = list(cganerrpath.glob("*.fits"))

# raw_err_plist.sort()
# cgan_err_plist.sort()
#%%



#%%
# merg_8_plist = list(merg8path.glob("*.fits"))
# merg_4_plist = list(merg4path.glob("*.fits"))
merg_7_plist = list(merg7path.glob("*.fits"))
err_7_plist = list(err7path.glob("*.fits"))
# merg_8_plist.sort()
# merg_4_plist.sort()
merg_7_plist.sort()
err_7_plist.sort()
#%%

# dat_raw = CCDData.read(raw_plist[0]).data.T
dat_cgan = CCDData.read(cgan_plist[2]).data.T
dat_7mer = CCDData.read(merg_7_plist[0]).data.T
dat_7err = CCDData.read(err_7_plist[0],unit='G').data.T
# dat_8mer = CCDData.read(merg_8_plist[0]).data.T
# dat_4mer = CCDData.read(merg_4_plist[0]).data.T
#%%
# np.mean(dat_raw)
np.mean(dat_cgan)
np.mean(dat_7mer)
# np.mean(dat_4mer)

# np.max(dat_raw)
np.max(dat_cgan)
np.max(dat_7mer)
# np.max(dat_4mer)
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

#%%
"""err check"""

test1 = CCDData.read(merg_7_plist[0]).data.T
test2 = CCDData.read(merg_7_plist[1]).data.T
err = test1 - test2

#%%


fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,3)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])
ax1.imshow(test1.T,
            vmin=-40,
            vmax=40,          
            origin = 'lower')
ax1.set_title("A = 16:00:00")
ax2.imshow(test2.T,
            vmin= -40,
            vmax = 40,
            origin = 'lower')
ax2.set_title("B = 16:05:15")
ax3.imshow(err.T,
            vmin= -40,
            vmax = 40,
            origin = 'lower')
ax3.set_title("A - B")


#%%
#Error became bigger because the time interval between the images
#%%
# peak_raw = peak_local_max(np.abs(dat_raw),min_distance=1
#                            # ,num_peaks = 1000
#                           ,threshold_abs = 12
#                           )
peak_cgan = peak_local_max(np.abs(dat_cgan),min_distance=2
                             # ,num_peaks = 1000
                           ,threshold_abs = 10
                           )
peak_7mer = peak_local_max(np.abs(dat_7mer),min_distance=2
                              # ,num_peaks = 1000
                            ,threshold_abs = 10
# peak_8mer = peak_local_max(np.abs(dat_8mer),min_distance=1
#                              # ,num_peaks = 1000
#                            ,threshold_abs = 10
                           )
# peak_4mer = peak_local_max(np.abs(dat_4mer),min_distance=1
#                              # ,num_peaks = 1000
#                            ,threshold_abs = 12
#                            )
#%%



# delpeaklist_cgan = []

# for i in range(len(peak_cgan)):
#     if np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1]]) > 40:
#         delpeaklist_cgan.append(i)

# peak_cgan = np.delete( peak_cgan, delpeaklist_cgan, 0 )

# delpeaklist_8mer = []

# for i in range(len(peak_8mer)):
#     if np.abs(dat_8mer[ peak_8mer[i][0], peak_8mer[i][1]]) > 40:
#         delpeaklist_8mer.append(i)

# peak_8mer = np.delete( peak_8mer, delpeaklist_8mer, 0 )

# delpeaklist_4mer = []

# for i in range(len(peak_4mer)):
#     if np.abs(dat_4mer[ peak_4mer[i][0], peak_4mer[i][1]]) > 40:
#         delpeaklist_4mer.append(i)

# peak_4mer = np.delete( peak_4mer, delpeaklist_4mer, 0 )

# delpeaklist_raw = []

# for i in range(len(peak_raw)):
#     if np.abs(dat_raw[ peak_raw[i][0], peak_raw[i][1]]) > 40:
#         delpeaklist_raw.append(i)

# peak_raw = np.delete( peak_raw, delpeaklist_raw, 0 )

#%%

for i in range(len(peak_cgan)) :
    if np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1] ]) > 40 :
        print(np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1] ]))

for i in range(len(peak_7mer)) :
    if np.abs(dat_7mer[ peak_7mer[i][0], peak_7mer[i][1] ]) > 40 :
        print(np.abs(dat_7mer[ peak_7mer[i][0], peak_7mer[i][1] ]))

# for i in range(len(peak_raw)) :
#     if np.abs(dat_raw[ peak_raw[i][0], peak_raw[i][1] ]) > 40 :
#         print(np.abs(dat_raw[ peak_raw[i][0], peak_raw[i][1] ]))

# for i in range(len(peak_4mer)) :
#     if np.abs(dat_4mer[ peak_4mer[i][0], peak_4mer[i][1] ]) > 40 :
#         print(np.abs(dat_4mer[ peak_4mer[i][0], peak_4mer[i][1] ]))
        
# for i in range(len(peak_8mer)) :
#     if np.abs(dat_8mer[ peak_8mer[i][0], peak_8mer[i][1] ]) > 30 :
#         print(np.abs(dat_8mer[ peak_8mer[i][0], peak_8mer[i][1] ]))

#%%
# len(peak_4mer)
len(peak_7mer)
# len(peak_raw)
len(peak_cgan)
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax1.imshow(dat_cgan.T,
            vmin=-40,
            vmax = 40,
           origin = 'lower')
ax1.set_title("cgan 45s")
ax2.imshow(dat_7mer.T,
            vmin=-40,
            vmax = 40,
           origin = 'lower')
ax2.set_title("raw 45*8 6m")
# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_cgan = CAp(peak_cgan,r=2)
cir_cgan.plot(ax1,color='r',lw=1,alpha=0.8)
cir_7mer = CAp(peak_7mer,r=2)
cir_7mer.plot(ax2,color='r',lw=1,alpha=0.8)


#%%
# peak_7mer_sub = sub_NT( peak_7mer, dat_7mer, 40 )
peak_7mer_ss = sub_NT_1st(peak_7mer, dat_7mer, 40)
# peak_8mer_sub = sub_NT( peak_8mer, dat_8mer, 30 )
# peak_cgan_sub = sub_NT( peak_cgan, dat_cgan, 40 )
peak_cgan_ss = sub_NT_1st(peak_cgan, dat_cgan, 40) 
#%%
len(peak_7mer)
len(peak_7mer_ss)
# len(peak_7mer_sub)

len(peak_cgan)
# len(peak_cgan_sub)
len(peak_cgan_ss)
#%%

fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax1.imshow(dat_cgan.T,
            vmin = -40,
            vmax = 40,
           origin = 'lower')
ax2.imshow(dat_7mer.T,
           vmin = -40,
           vmax = 40,
           origin = 'lower')
ax1.set_title(f"cgan 45s NW sub 40 : {len(peak_cgan_ss)}")
ax2.set_title(f"raw 45*8 6m NW sub 40: {len(peak_7mer_ss)}")
# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_7mer = CAp(peak_cgan_ss,r=2)
cir_7mer.plot(ax1,color='r',lw=1,alpha=0.5)
cir_7mer_sub = CAp(peak_7mer_ss,r=2)
cir_7mer_sub.plot(ax2,color='r',lw=1,alpha=0.5)
#%%
len(peak_7mer_sub)
#%%
#error find


#%%

fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax1.imshow(dat_cgan.T,
            vmin = -40,
            vmax = 40,
           origin = 'lower')
ax2.imshow(dat_cgan.T,
           vmin = -40,
           vmax = 40,
           origin = 'lower')
ax1.set_title("cgan 45s")
ax2.set_title("cgan 45s NW sub")
# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_cgan = CAp(peak_cgan,r=2)
cir_cgan.plot(ax1,color='r',lw=1,alpha=0.5)
cir_cgan_sub = CAp(peak_cgan_sub,r=2)
cir_cgan_sub.plot(ax2,color='r',lw=1,alpha=0.5)

#%%
problem = np.array([np.array([128,174]),np.array([127,14])] )
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.imshow(dat_cgan.T,
            vmin = -40,
            vmax = 40,
           origin = 'lower')
ax1.set_title("cgan 45s")
# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_problem = CAp(problem,r=10)
cir_problem.plot(ax1,color='m',lw=3,alpha=0.8)
#%%
peak_7_list = []
dat_7_list = []

for i in range( len(merg_7_plist) ):
    ccd = CCDData.read(merg_7_plist[i]).data.T
    peak = peak_local_max(np.abs(ccd), min_distance = 2,exclude_border=False
                              # ,num_peaks = 1000
                            ,threshold_abs = 10)
    peak = sub_NT_1st(peak,ccd,40)
    peak_7_list.append(peak)
    dat_7_list.append(ccd)
    
print(i,"F")

#%%

    
fig = plt.figure(
    figsize=(24,6)
    )
gs = gridspec.GridSpec(1,4)
ax = []
for i in range(4):
    ax.append(plt.subplot(gs[0,i]))
    ax[i].imshow(dat_7_list[i].T,
            vmin = -40,
            vmax = 40,
           origin = 'lower')
    cir_7mer = CAp(peak_7_list[i],r=2)
    cir_7mer.plot(ax[i],color='r',lw=1,alpha=0.5)
    ax[i].set_title(len(peak_7_list[i]))
    # ax[i].set_xlim(100,150)
    # ax[i].set_ylim(100,150)

    






#%%
#Find the pixels whose values are higher than 40 G 

delpeaklist_cgan = []

for i in range(len(peak_cgan)):
    if np.abs(dat_cgan[ peak_cgan[i][0], peak_cgan[i][1]]) > 30:
        delpeaklist_cgan.append(i)

delpeakpos_cgan = []
for i in range(len(delpeaklist_cgan)):
    delpeakpos_cgan.append(list( peak_cgan[ delpeaklist_cgan[i] ] ) )

#%%

for i in range(len(peak_cgan)):
    y = find_dir8( list(peak_cgan[i]) , delpeakpos_cgan )
    if y == "pos6":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0]+1 ][ peak_cgan[i][1] ] > 0:
            delpeaklist_cgan.append(i)
    elif y == "pos4":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0]-1 ][ peak_cgan[i][1] ] > 0:
            delpeaklist_cgan.append(i)
    elif y == "pos8":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0] ][ peak_cgan[i][1]+1 ] > 0:
            delpeaklist_cgan.append(i)
    elif y == "pos2":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0] ][ peak_cgan[i][1]-1 ] > 0:
            delpeaklist_cgan.append(i)
    elif y == "pos9":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0]+1 ][ peak_cgan[i][1]+1 ] > 0:
            delpeaklist_cgan.append(i)
    elif y == "pos7":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0]-1 ][ peak_cgan[i][1]+1 ] > 0:
            delpeaklist_cgan.append(i)
    elif y == "pos3":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0]+1 ][ peak_cgan[i][1]-1 ] > 0:
            delpeaklist_cgan.append(i)
    elif y == "pos1":
        if  dat_cgan[peak_cgan[i][0]][peak_cgan[i][1]] * dat_cgan[ peak_cgan[i][0]-1 ][ peak_cgan[i][1]-1 ] > 0:
            delpeaklist_cgan.append(i)
 

peak_cgan = np.delete( peak_cgan, delpeaklist_cgan, 0 )

len(np.unique(delpeaklist_cgan))


#%%
dat_cgan_6list = []
peak_cgan_6list = []

for i in range(6):
    dat_cgan_6list.append(CCDData.read(cgan_plist[i]).data.T)
    peak_cgan_6list.append( peak_local_max(np.abs(dat_cgan_6list[i]),min_distance=1
                             # ,num_peaks = 1000
                           ,threshold_abs = 10
                           ) 
                           )
    delpeaklist_cgan_6list = []

    for j in range(len(peak_cgan_6list[i])):
        if np.abs(dat_cgan_6list[i][ peak_cgan_6list[i][j][0],
                                 peak_cgan_6list[i][j][1] ]) > 40:
            delpeaklist_cgan_6list.append(j)
    
    peak_cgan_6list[i] = np.delete( peak_cgan_6list[i], delpeaklist_cgan_6list, 0 )

#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(2,5)
axs = []
axs.append(plt.subplot(gs[0,0]))
axs.append(plt.subplot(gs[0,1]))
axs.append(plt.subplot(gs[0,2]))
axs.append(plt.subplot(gs[1,0]))
axs.append(plt.subplot(gs[1,1]))
axs.append(plt.subplot(gs[1,2]))
ax = plt.subplot(gs[0:2,3:5])

cir_cgan_6list = []
for i in range(6):
    axs[i].imshow(np.abs(dat_cgan_6list[i].T),
                  vmin = 0,
                  vmax = 40,
                  origin = 'lower')
    axs[i].set_xlim(100,150)
    axs[i].set_ylim(100,150)
    cir_cgan_6list.append(CAp(peak_cgan_6list[i],r=1))
    cir_cgan_6list[i].plot(axs[i],color='r',lw=1,alpha=0.8)

    
    
ax.imshow(np.abs(dat_8mer.T),
          vmin = 0,
          vmax = 40,
          origin = 'lower')
ax.set_xlim(100,150)
ax.set_ylim(100,150)

cir_8mer = CAp(peak_8mer,r=2)
cir_8mer.plot(ax,color='r',lw=1,alpha=0.8)

#%%

fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(2,3)
axs = []
axs.append(plt.subplot(gs[0,0]))
axs.append(plt.subplot(gs[0,1]))
axs.append(plt.subplot(gs[0,2]))
axs.append(plt.subplot(gs[1,0]))
axs.append(plt.subplot(gs[1,1]))
axs.append(plt.subplot(gs[1,2]))
# ax = plt.subplot(gs[0:2,3:5])

cir_cgan_6list = []
for i in range(6):
    axs[i].imshow(np.abs(dat_cgan_6list[i].T),
                  # vmin = 0,
                  # vmax = 40,
                  origin = 'lower')
    # axs[i].set_xlim(100,150)
    # axs[i].set_ylim(100,150)
    # cir_cgan_6list.append(CAp(peak_cgan_6list[i],r=1))
    # cir_cgan_6list[i].plot(axs[i],color='r',lw=1,alpha=0.8)

    
    
# ax.imshow(np.abs(dat_8mer.T),
#           vmin = 0,
#           vmax = 40,
#           origin = 'lower')
# ax.set_xlim(100,150)
# ax.set_ylim(100,150)

# cir_8mer = CAp(peak_8mer,r=2)
# cir_8mer.plot(ax,color='r',lw=1,alpha=0.8)

#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(2,2)
# ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
# ax3 = plt.subplot(gs[1,0])
ax4 = plt.subplot(gs[1,1])
# ax1.imshow(dat_raw.T,
#             vmin=-100,
#             vmax=100,          
#            origin = 'lower')
# ax1.set_title("raw 45s")
ax2.imshow(dat_cgan.T,
            vmin=-100,
            vmax = 100,
           origin = 'lower')
ax2.set_title("cgan 45s")
# ax3.imshow(dat_4mer.T,
#             vmin=-100,
#             vmax=100,           
#            origin = 'lower')
# ax3.set_title("raw 45*4 3m")
ax4.imshow(dat_8mer.T,
            vmin=-100,
            vmax = 100,
           origin = 'lower')
ax4.set_title("raw 45*8 6m")

# cir_raw = CAp(peak_raw,r=2)
# cir_raw.plot(ax1,color='r',lw=1,alpha=0.8)
cir_cgan = CAp(peak_cgan,r=2)
cir_cgan.plot(ax2,color='r',lw=1,alpha=0.8)
# cir_4mer = CAp(peak_4mer,r=2)
# cir_4mer.plot(ax3,color='r',lw=1,alpha=0.8)
cir_8mer = CAp(peak_8mer,r=2)
cir_8mer.plot(ax4,color='r',lw=1,alpha=0.8)


