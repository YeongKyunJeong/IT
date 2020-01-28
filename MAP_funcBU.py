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

# def err1(img):
#     for findex in range(len(img) - 1 ):
        
#         ccd1 = CCDData.read(img[findex]).data
#         ccd2 = CCDData.read(img[findex + 1]).data
        
#         ccderr = fits.PrimaryHDU( (ccd1 - ccd2) )
#         ccderr = fits.HDUList( [ccderr] )
        
#         if 'cgan' in str( img[findex]) :
#             ccderr.writeto(cganerrpath/f"err_cgan{findex+1:03}.fits",overwrite=True)
#         else:
#             ccderr.writeto(rawerrpath/f"err_raw{findex+1:03}.fits",overwrite=True)
        
        
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
        merg_list.header = CCDData.read(img[n*j + 1]).header
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
# def div(pos,img):
    
#     ps = np.zeros([img.shape[0],img.shape[1]])

#     for i in range(img.shape[0]):
#         for j in range(img.shape[1]):
#             dst = []
#             for k in range(len(pos)):
#                 dst.append(ds(k,[i,j]))
#             for k in range(len(dst)):
#                 if dst[k] == np.min(dst):
#                     ps[i][j] = k
#             print(i,j)
#     return ps
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
