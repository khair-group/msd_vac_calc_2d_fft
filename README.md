# msd_vac_calc_2d_fft

Calculation of the mean square displacement and velocity autocorrelation is simplified through a Fourier component analysis.

Code developed by Jeremy Kach.
 
See: Calandrini et al. École thématique de la Société
    Française de la Neutronique, 12:201-232, 2011. URL: https://www.neutron-sciences.org/articles/sfn/pdf/2011/01/sfn201112010.pdf 
for theoretical explanation and relevant formulae.

A general code to implement this was written in 
Python by Jeremy Kach. RK modified the code to handle: 
                (a) .mat files imported from MATLAB 
                (b) 2-dimensional simulation data
                (c) Averaging over specified particle indices when needed, as opposed to 
                    averaging over all the particles

The core-calculations are performed in "displacement.py", while the driver routine "calc_msd_thru_fft.py" performs the unwrapping of the coordinates and passes the appropriate input needed for "displacement.py". A sample usage statement, along with the explanation of the input parameters is given below:

Sytntax: [calc_msd_thru_fft.py] [length_of_box] [delta t] [start trim] [end trim] [steps] [path to input file] [path to list of indices]

[length_of_box]           : length of periodic box (float)
[delta t]                 : timestep width used in simulation
[start trim]              : line number from which data is used for processing, set as -1 if you don't want to trim any lines from the beginning
[end trim]	          : line number upto which data is used for processing, set to a large number if you want to consider all lines till end of file
[steps]                   : number of timesteps to cover in a single step. Use 1 if each timestep needs to be processed
[path to input file]      : path to input file, containing the raw, wrapped positions
[path to list of indices] : path to file containing the indices to average over. If this is -1, average over all input

Sample usage, with the input files included in this folder:

calc_msd_thru_fft.py 10. 1. -1 9999999 1 inp_wrp_pos.mat -1

should produce the following files as output:

msd_tot_from_fft.dat 
msd_xx_from_fft.dat  
msd_xy_from_fft.dat  
msd_yy_from_fft.dat 

Running the MATLAB code "plot_msd.m" produces a plot of the analytical MSD result of a Peruani-Morelli disk, and compares it with the numerically computed result stored in "msd_tot_from_fft.dat"

NOTE: the .mat files have to be saved using the flag '-v7.3', so that they may be interpreted as HDF5 entities. See
https://stackoverflow.com/questions/874461/read-mat-files-in-python/19340117#19340117
for more details.


