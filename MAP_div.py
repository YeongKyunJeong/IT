import os
import glob

from pathlib import Path
import numpy as np

import matplotlib as mp
from matplotlib import pyplot as plt
from matplotlib import gridspec, rcParams, rc
from matplotlib import colors as mcolors

import astropy.units as u
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.wcs import WCS

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
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

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
        
def sub_NT_1st(peak_pos, img, threshold = 70, lim = 4 ):
    del_ind = []

    for i in range(len(peak_pos)):
        if np.abs(img[ peak_pos[i][0], peak_pos[i][1]]) > threshold:
            del_ind.append(i)
        elif peak_pos[i][0] <= lim :
            del_ind.append(i)
        elif peak_pos[i][0] >= len(img[0])-lim :
            del_ind.append(i)            
        elif peak_pos[i][1] <= lim :
            del_ind.append(i)
        elif peak_pos[i][1] >= len(img[1])-lim :
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
ap = [0,1,2,3,4,5]
r = [1,2,3,4,5,6]
area = [1,9,25,45,69,109]
sec2cm = 1.496e13/3600/180*np.pi
pix2sec = 0.6
pix22cm2 = (pix2sec*sec2cm)**2
#%%
print(f"{pix22cm2:.4}")
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
dat_raw = CCDData.read(raw_plist[0]).data.T
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.imshow(dat_raw.T,
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

# cir_cgan = CAp(peak_cgan_ss,r=1)
# cir_cgan.plot(ax1,color='r',lw=1,alpha=0.8)
# # cir_cgan = CAp(peak_cgan_ss,r=2)
# # cir_cgan.plot(ax1,color='b',lw=1,alpha=0.8)
# # cir_cgan = CAp(peak_cgan_ss,r=3)
# # cir_cgan.plot(ax1,color='y',lw=1,alpha=0.8)
# cir_cgan = CAp(peak_cgan_ss,r=4)
# cir_cgan.plot(ax1,color='g',lw=1,alpha=0.8)


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
ax2.set_title(f"45s Peak Position: {len(peak_cgan_ss)}\nmin distance: 2 pxiels, min = 9, max = 60")

# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

cir_cgan = CAp(peak_cgan_ss,r=1)
cir_cgan.plot(ax2,color='r',lw=1,alpha=0.8)
#%%
aperture = []
for i in range(1,7):
    aperture.append(CAp(peak_cgan_ss,r=i))
#%%
# AP(dat_cgan.T, aperture,method='center')
#%%

IN = AP(dat_cgan.T, aperture,method='center')
IN['aperture_sum_0'].info.format = '%.8g'
IN['aperture_sum_1'].info.format = '%.8g'
IN['aperture_sum_2'].info.format = '%.8g'
IN['aperture_sum_3'].info.format = '%.8g'
IN['aperture_sum_4'].info.format = '%.8g'
IN['aperture_sum_5'].info.format = '%.8g'
#%%
print(IN)
#%%
dslist = []
for i in range(len(peak_cgan_ss)):
    dss = []
    for j in range(len(peak_cgan_ss)):
        if i != j:
            dss.append(ds(peak_cgan_ss[i],peak_cgan_ss[j]))
    dslist.append(np.min(dss))
#%%
len(dslist)
#%%
INtest = IN

for i in range(len(peak_cgan_ss)):    
    if dslist[i] < 5:
        print(i,dslist[i])
#%%
index = []
for i in dslist:
    index.append(int(i/2))
#%%
for i in range(len(index)):
    print(dslist[i],index[i])
    
#%%
for i in range(len(peak_cgan_ss)):
    for j in range(index[i],6):
        INtest[i][f'aperture_sum_{j}'] = 0

#%%
for i in range(0,500,40):
    print(INtest[i],dslist[i],index[i])
#%%
# xxx = np.arange(1,7)

flux_r = []
for i in range(len(peak_cgan_ss)):
    fr = []
    for j in range(0,6):
        fr.append(INtest[i][f'aperture_sum_{j}'])
    flux_r.append(fr)
#%%
ap_list=[]
for i in range(len(peak_cgan_ss)):
    if flux_r[i][0] > 0 :
        if peak_local_max(np.array(flux_r[i])).shape[0] > 0:
           r_1 =  peak_local_max(np.array(flux_r[i]))[0][0]
        else:
           x = np.max(flux_r[i])
           r_1 = np.where( flux_r[i] == x)[0][0]

    else:
        if peak_local_max(-np.array(flux_r[i])).shape[0] > 0:
           r_1 =  peak_local_max(-np.array(flux_r[i]))[0][0]
        else:
           x = np.min(flux_r[i])
           r_1 = np.where( flux_r[i] == x)[0][0]

    if r_1 > 0:
        k = 0
        while k < 0.2  : #?
            k = (np.abs(flux_r[i][r_1])-np.abs(flux_r[i][r_1-1]))/np.abs(flux_r[i][r_1-1])
            if k < 0.2:
                r_1 = r_1 - 1
    ap_list.append(r_1)
#%%

area_list = []
for i in range(len(peak_cgan_ss)):
    if ap_list[i] == 0:
        area_list.append(area[0])
    elif ap_list[i] == 1 :
        area_list.append(area[1])
    elif ap_list[i] == 2 :
        area_list.append(area[2])
    elif ap_list[i] == 3 :
        area_list.append(area[3])
    elif ap_list[i] == 4 :
        area_list.append(area[4])
    elif ap_list[i] == 5 :
        area_list.append(area[5])
    print(area_list[i],ap_list[i])        

area_list = np.array(area_list)*pix22cm2
#%%
ap_list
a=0
b=0
c=0
d=0
e=0
f=0
for i in range(len(peak_cgan_ss)):
    if ap_list[i] == 0 :
        a+=1
    elif ap_list[i] == 1 :
        b+=1
    elif ap_list[i] == 2 :
        c+=1
    elif ap_list[i] == 3 :
        d+=1
    elif ap_list[i] == 4 :
        e+=1
    elif ap_list[i] == 5 :
        f+=1
print(a,b,c,d,e,f, a+b+c+d+e+f)
#%%
flux_tot = []
for i in range(len(peak_cgan_ss)):
    flux_tot.append(flux_r[i][ap_list[i]])

Mx = np.array(flux_tot)*pix22cm2

#%%

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.plot(area_list, np.abs(Mx),marker = 'o',ls = '',alpha = 0.7, color = colors['coral'])
ax1.set_xlabel('Area (cm^2)',fontsize=15)

# ax1.plot(area_list, np.abs(Mx),marker = 'o',ls = '',alpha = 0.7, color = colors['coral'])
# ax1.set_xlabel('Area (pixel)',fontsize=15)

ax1.set_ylabel('Flux Strength (Mx)',fontsize=15)
ax1.set_title("Flux Strength to Area",fontsize = 17)
# ax1.set_xticks([16,17,18])
# ax1.set_yscale('log') 

#%%
for i in range(len(peak_cgan_ss)):    
    if index[i] >= 5:
        print(i,index[i],peak_cgan_ss[i],ap_list[i],flux_tot[i])
        print(flux_r[i])
#%%
peak_r1 = []
peak_r2 = []
peak_r3 = []
peak_r4 = []
peak_r5 = []
peak_r6 = []
for i in range(len(peak_cgan_ss) ):
    if ap_list[i] == 0:
        peak_r1.append(peak_cgan_ss[i])
    elif ap_list[i] == 1:
        peak_r2.append(peak_cgan_ss[i])
    elif ap_list[i] == 2:
        peak_r3.append(peak_cgan_ss[i])
    elif ap_list[i] == 3:
        peak_r4.append(peak_cgan_ss[i])
    elif ap_list[i] == 4:
        peak_r5.append(peak_cgan_ss[i])
    elif ap_list[i] == 5:
        peak_r6.append(peak_cgan_ss[i])
#%%
fig = plt.figure(figsize=(10,8))
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

# ax2.scatter(delete_18_x,delete_18_y,
#             marker = 'x',
#             c = 'm',
#             s = 100,
#             # label = 'dot graph'
#             )

if len(peak_r1) != 0:
    cir_cgan = CAp(peak_r1,r=1)
    cir_cgan.plot(ax2,color='r',lw=1.5,alpha=0.8)
if len(peak_r2) != 0:
    cir_cgan = CAp(peak_r2,r=2)
    cir_cgan.plot(ax2,color='b',lw=1.5,alpha=0.8)
if len(peak_r3) != 0:
    cir_cgan = CAp(peak_r3,r=3)
    cir_cgan.plot(ax2,color='y',lw=1.5,alpha=0.8)
if len(peak_r4) != 0:
    cir_cgan = CAp(peak_r4,r=4)
    cir_cgan.plot(ax2,color='g',lw=1.5,alpha=0.8)
if len(peak_r5) != 0:
    cir_cgan = CAp(peak_r5,r=5)
    cir_cgan.plot(ax2,color='w',lw=1.5,alpha=0.8)
if len(peak_r6) != 0:
    cir_cgan = CAp(peak_r6,r=6)
    cir_cgan.plot(ax2,color='k',lw=1.5,alpha=0.8)

#%%
print(f"{pix22cm2:.4}")

#%%
area_list = []
for i in range(len(peak_cgan_ss)):
    area_list.append(area[ap_list[i]])
#%%
for i in range(len(peak_cgan_ss)):
    print(area_list[i],ap_list[i]+1)    

#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.plot(area_list, np.abs(Mx),marker = '.',ls = '')
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
# np.log10(np.abs(Mx)),

# ax1.set_xscale('log')
ax1.hist(np.log10(np.abs(Mx)),
         log = 'True',
         bins=40,
         range = [16,18],
         rwidth = 0.6,
         color=colors['silver']
         )
ax1.set_title("IN Element Flux Distribution")
ax1.set_xlabel("log Flux (Mx)",fontsize=15)
ax1.set_ylabel("Number", fontsize=15)
ax1.set_yticks([10,20,30,40,50,60],minor = False)
# ax1.grid(axis='y', color='k',alpha=0.25)
# ax1.set_xticks([16,17,18])


#%%
print(f"{10**17.1:.4}")


#%%
# 10*14 Mx interval
num_Mx = np.zeros([41],dtype=int)
for i in range(len(peak_cgan_ss)):
    if np.abs(Mx[i]) < 1e18:
        j = int( (np.log10(np.abs(Mx[i])) - 16 ) / ( 2/40) )
        num_Mx[j] += 1
    else:
        num_Mx[40] += 1
        print(Mx[i])

#%%
np.max(Mx)
np.median(np.abs(flux_tot))
np.mean(np.abs(area_list))
print(f"{np.sum(np.abs(area_list)):.4}")
np.where(np.abs(Mx) > 1e18)
#%%
np.mean(np.abs(Mx))
np.sum(np.abs(Mx))
#%%
# by visual way, exclude [221,28]
delete_18 = []
delete_18_x = []
delete_18_y = []
delete_r5 = []
delete_r5_x = []
delete_r5_y = []

for i in (np.where(np.abs(Mx)>1e18)[0]):
    delete_18.append(peak_cgan_ss[i])
    delete_18_x.append(peak_cgan_ss[i][0])
    delete_18_y.append(peak_cgan_ss[i][1])
# peak_cgan_ss[555]

for i in ( np.where(np.array(ap_list) > 3)[0] ):
    delete_r5.append(peak_cgan_ss[i])
    delete_r5_x.append(peak_cgan_ss[i][0])
    delete_r5_y.append(peak_cgan_ss[i][1])
#%%
fig = plt.figure(figsize=(10,10))
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

ax2.set_title("Aperture of Magnetic Element \n Between 9G and 60G",fontsize = 17)
# ax2.set_title("45s 2017.06.14. 16:00:00",fontsize = 17)
# ax1.set_xlim(100,150)
# ax1.set_ylim(100,150)
# ax2.set_xlim(100,150)
# ax2.set_ylim(100,150)

delete_18_x= [210, 205, 192, 117, 77]
delete_18_y = [87, 14, 245, 121, 5]

ax2.scatter(delete_18_x,delete_18_y,
            marker = 'x',
            c = colors['fuchsia'],
            s = 100,
            label = 'Exluded',
            alpha = 0.8
            )

if len(peak_r1) != 0:
    cir_cgan = CAp(peak_r1,r=1)
    cir_cgan.plot(ax2,color='r',lw=1.5,alpha=0.8)
if len(peak_r2) != 0:
    cir_cgan = CAp(peak_r2,r=2)
    cir_cgan.plot(ax2,color='b',lw=1.5,alpha=0.8)
if len(peak_r3) != 0:
    cir_cgan = CAp(peak_r3,r=3)
    cir_cgan.plot(ax2,color='y',lw=1.5,alpha=0.8)
if len(peak_r4) != 0:
    cir_cgan = CAp(peak_r4,r=4)
    cir_cgan.plot(ax2,color= colors['orangered'],lw=1.5,alpha=0.8)
if len(peak_r5) != 0:
    cir_cgan = CAp(peak_r5,r=5)
    cir_cgan.plot(ax2,color='w',lw=1.5,alpha=0.8)
if len(peak_r6) != 0:
    cir_cgan = CAp(peak_r6,r=6)
    cir_cgan.plot(ax2,color='k',lw=1.5,alpha=0.8)
# ax2.legend()

#%%

# num_Mx = np.delete( peak_pos, del_ind, 0 ) 


# 10*14 Mx interval
num_Mx = np.zeros([41],dtype=int)
for i in range(len(peak_cgan_ss)):
    if np.abs(Mx[i]) < 1e18:
        j = int( (np.log10(np.abs(Mx[i])) - 16 ) / ( 2/40) )
        num_Mx[j] += 1
    else:
        num_Mx[40] += 1
        print(Mx[i])
#%%
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
ax1.plot(num_Mx,
         marker = '.',
         linestyle = "None")  

#%%
np.log(10)

#%%    
np.where(peak_cgan_ss == np.array([142,64]))
peak_cgan_ss[252]
INtest[252]
ap_list[442]
peak_local_max(np.array(flux_r[442]))[0][0]




#%%

test = np.ones([19,19])
cet = [9,9]
aa = CAp(cet,r=5)  
testin = AP(test,
            aa,
            method = 'center')
#%%
testin

#%%
aperture = []
for i in range(1,7):
    aperture.append(CAp(peak_cgan_ss,r=i))
#%%
# AP(dat_cgan.T, aperture,method='center')
#%%

IN = AP(dat_cgan.T, aperture,method='center')
IN['aperture_sum_0'].info.format = '%.8g'
IN['aperture_sum_1'].info.format = '%.8g'
IN['aperture_sum_2'].info.format = '%.8g'
IN['aperture_sum_3'].info.format = '%.8g'
IN['aperture_sum_4'].info.format = '%.8g'
IN['aperture_sum_5'].info.format = '%.8g'