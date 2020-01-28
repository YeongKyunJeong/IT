import os
import glob

from pathlib import Path
import numpy as np
import matplotlib as mp

from matplotlib import pyplot as plt

import astropy.units as u
from astropy.nddata import CCDData
from astropy.io import fits

from skimage.feature import peak_local_max
from ccdproc import trim_image

#%%