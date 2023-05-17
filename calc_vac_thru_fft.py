#!/usr/local/bin/python3

### 2-FEB-2023: Code shared by Jeremy Kach
### 20-FEB-2023: Modified by RK to handle: 
#                (a) .mat files imported from MATLAB 
#                (b) 2-dimensional simulation data
#                (c) Averaging over specified particle indices when needed, as opposed to 
#                    averaging over all the particles 


### The unwrapping is also performed within the Python script. Only the raw
### coordinates are imported from MATLAB.

import os
import sys
import numpy as np
import time
from vac_tools import DisplacementTensor, fft, ifft, s1
import h5py

l = len(sys.argv)

if l<8:
    print("\n Correct syntax is [calc_vac_thru_fft.py] [length_of_box] [delta t] [start trim] [end trim] [steps] [path to input file] [path to list of indices]\n")
    exit(0)

else :
    pass

L=float(sys.argv[1])        # length of simulation box 
dt=float(sys.argv[2])       # timestep width used in simulation 
start_trim=int(sys.argv[3]) # line number from which data is used for processing
end_trim=int(sys.argv[4])   # line number upto which data is used for processing
steps=int(sys.argv[5])      # number of timesteps to cover in a single step. Use 1 if each timestep needs to be used
file_path=sys.argv[6]       # path to input file, containing the raw, wrapped positions
index_path=sys.argv[7]      # path to file containing the indices to average over. If this is -1, average over all input

f = h5py.File(file_path,'r')
for key in f:
    data = f.get(key)
    data = np.array(data)

raw_pos = data.T
tot_len = raw_pos.shape[0]
N_particles = raw_pos.shape[1]

if (start_trim<0):
    start_trim=0

if (end_trim>(tot_len-1)):
    end_trim=tot_len-1


if (index_path == '-1'):
    av_ind = np.arange(0,N_particles).reshape(N_particles,1)
else:
    f = h5py.File(index_path,'r')
    for key in f:
        av_ind = f.get(key)
        av_ind = np.array(av_ind)   
        av_ind = av_ind -1 #conversion from MATLAB's 1-based indexing to Python's 0-based indexing

only_val = av_ind[:,0].astype(int) # stores the list of indices to average over


#print("av_ind: ")
#print(av_ind)

##### unwrapping process begins here
unwrap_pos = np.zeros((tot_len, N_particles, 2))
unwrap_pos[0, :, :] = raw_pos[0,:,:]

for i in range(1,tot_len):
    snap = raw_pos[i,:,:]
    snapprev = raw_pos[i-1,:,:]
    dr = snap - snapprev
    dr[dr > L/2] -= L
    dr[dr < -L/2] += L
    unwrap_pos[i, :, :] = unwrap_pos[i-1, :, :] + dr

vel_vec = np.zeros((tot_len-1, N_particles, 2))
for i in range(tot_len-1):
    vel_vec[i,:,:] = (unwrap_pos[i+1,:,:] - unwrap_pos[i,:,:])/dt


### Remember: the unwrapping is performed all timesteps and particles. Appropriate trimming and
##  selection of particles in subsequent steps, as shown below

sel_vel=vel_vec[:,only_val,:]

if (start_trim<0):
    start_trim=0

if (end_trim>(tot_len-2)):
    end_trim=tot_len-2

st = time.time()
vac = DisplacementTensor(sel_vel, dt, start_trim, end_trim, steps)
vac_tens = vac.compute()
et = time.time()

print("VAC Tensor calculated in {} seconds".format(et-st))
#
tval = dt*np.arange(0,tot_len-1)
vac_xx = vac_tens[:, 0, 0]
vac_yy = vac_tens[:, 1, 1]
vac_xy = vac_tens[:, 0, 1]
vac_tot = vac_xx + vac_yy

Nsample = vac_xx.shape[0]
#
#
#opfile=open("msd_xx_from_fft.dat","w")
#for i in range(Nsample):
#    opfile.write("%20.8f\t%20.15f\n"
#    %(tval[i],msd_xx[i]))
#opfile.close()
#
#opfile=open("msd_yy_from_fft.dat","w")
#for i in range(Nsample):
#    opfile.write("%20.8f\t%20.15f\n"
#    %(tval[i],msd_yy[i]))
#opfile.close()
#
#opfile=open("msd_xy_from_fft.dat","w")
#for i in range(Nsample):
#    opfile.write("%20.8f\t%20.15f\n"
#    %(tval[i],msd_xy[i]))
#opfile.close()
#
opfile=open("vac_from_fft.dat","w")
for i in range(Nsample):
    opfile.write("%20.8f\t%20.15f\n"
    %(tval[i],vac_tot[i]))
opfile.close()
#
