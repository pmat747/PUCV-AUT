# Standard python numerical analysis imports:
import numpy as np
import scipy
from scipy import signal
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math

#########################################################################
mpl.rcParams['text.usetex'] = True

#########################################################################
# Get ET spectrum
ET_D_psd=np.loadtxt('ET_D_sensitivity.txt', skiprows=0)
ET_B_psd=np.loadtxt('ET_B_sensitivity.txt', skiprows=0)
ET_C_psd=np.loadtxt('ET_C_sensitivity.txt', skiprows=0)

# Get CE spectrum
CE1_psd=np.loadtxt('curves_Jan_2020/ce1.txt', skiprows=0)

# Get CE spectrum
CE2_psd=np.loadtxt('curves_Jan_2020/ce2.txt', skiprows=0)

# Get aLIGO spectrum
AL_psd=np.loadtxt('ALIGO_sensitivity.txt', skiprows=7)

# Get ADV spectrum
ADV_psd=np.loadtxt('ADV_sensitivity.txt', skiprows=7)


# interpolation to get regular sampling
deltaf=.5                                # deltaf=0.5 Hz

fNET=ET_D_psd[-1,0];                       # Nyquist freq
freq_ET_psd = np.arange(1., fNET, deltaf)
psd = interpolate.interp1d(ET_D_psd[:,0],ET_D_psd[:,3])
ET_D_psd_res = psd(freq_ET_psd)

fNCE=CE1_psd[-1,0];                       # Nyquist freq
freq_CE1_psd = np.arange(5., fNCE, deltaf)
psd = interpolate.interp1d(CE1_psd[:,0],CE1_psd[:,1])
CE1_psd_res = psd(freq_CE1_psd)

fNCE=CE2_psd[-1,0];                       # Nyquist freq
freq_CE2_psd = np.arange(5., fNCE, deltaf)
psd = interpolate.interp1d(CE2_psd[:,0],CE2_psd[:,1])
CE2_psd_res = psd(freq_CE2_psd)


freq_AL_psd = np.arange(9., 8000., deltaf)
psd = interpolate.interp1d(AL_psd[:,0],AL_psd[:,6])
AL_psd_res = psd(freq_AL_psd)


#########################################################################

# plot the spectra
fmin=10
fmax=3000

fig1 = plt.figure()
plt.loglog(ET_D_psd[:,0],ET_D_psd[:,3],'green',label='ET-D HF+LF')
plt.loglog(ET_B_psd[:,0],ET_B_psd[:,1],'red',label='ET-B')
plt.loglog(ET_C_psd[:,0],ET_C_psd[:,1],'k',label='ET-C')
plt.loglog(CE1_psd[:,0],CE1_psd[:,1],'cyan',label='CE stage1')
plt.loglog(CE2_psd[:,0],CE2_psd[:,1],'m',label='CE stage2')
plt.loglog(AL_psd[:,0],AL_psd[:,6],'b',label='aLIGO design')
#plt.loglog(ADV_psd[:,0],ADV_psd[:,6],'m',label='ADV design')

plt.axis([1, fmax, 1e-25, 4e-20])

plt.grid(True)
plt.ylabel('$\mathrm{ASD [strain/\sqrt{Hz}]}$')
plt.xlabel('$\mathrm{Frequency [Hz]}$')
plt.legend(loc='upper right')
plt.savefig('spectrum.png')
plt.savefig('spectrum.pdf')
plt.savefig('spectrum.eps')
plt.show()




