#!/usr/local/bin/python3

### 2-FEB-2023: Code shared by Jeremy Kach
### 20-FEB-2023: Modified byy RK to handle 2-dimensional simulation data


import numpy as np
import pyfftw


def fft(x, n, axis):
    a = pyfftw.empty_aligned(x.shape, 'complex128')
    a[:] = x
    fft_object = pyfftw.builders.fft(a, n=n, axis=axis)
    return fft_object()


def ifft(x, axis):
    a = pyfftw.empty_aligned(x.shape, 'complex128')
    a[:] = x
    fft_object = pyfftw.builders.ifft(a, axis=axis)
    return fft_object()


def s1(x, y):
    Nt = x.shape[0]
    dprod = np.multiply(x, y)
    dprod = np.append(dprod, np.zeros((1, dprod.shape[1])), axis=0)
    sumsq = 2*dprod.sum(axis=0)
    S1 = np.zeros((Nt, x.shape[1]))
    for tau in range(Nt):
        sumsq = sumsq - dprod[tau-1] - dprod[Nt-tau]
        S1[tau, :] = sumsq / (Nt - tau)
    return S1


class DisplacementTensor:
    """
    Class containing functions needed to compute and store displacement tensor.
    Method used from Calandrini et al. École thématique de la Société
    Française de la Neutronique, 12:201-232, 2011. Code adapted from
    https://stackoverflow.com/questions/34222272/
    computing-mean-square-displacement-using-python-and-fft
    ...

    Attributes
    ----------
    pos : float
        (N_T,N,2) ndarray of unwrapped particle positions, where N_T is the
        total number of timesteps from the simulation, N is the total number of
        particles, and 2 is for positions in two dimensions
    dt : float
        size of time step
    st : int
        timestep to start computing displacement tensor
    et : int
        last timestep to include in computation
    step_size : int
        Number of timesteps covered in a single step (e.g. if set to 10, will
        take positions every 10 timesteps). Default is 1.
    steps : int
        Number of timesteps used for entire displacement calculation
    rr : float
        (Nt,) ndarray, where Nt is the number of timesteps used in the
        displacement calculation (N_T / step_size). Each component (xx, xy, ...
        zz) is stored individually - e.g. to access the xy displacement, use
        the call DisplacementTensor.xy


    Methods
    -------
    autocorr
        Computes autocorrelation for 6 unique components of the displacement
        tensor
    displacementTensor
        Main method that calls autocorrFFT to calculate
        full displacement tensor. Note that the displacement tensor is
        symmetric (i.e. DisplacementTensor.xy = DisplacementTensor.yx),
        therefore only 3 components are computed (xx, xy, yy)
    """

    def __init__(
        self,
        unwrap_pos,
        dt,
        start_trim,
        end_trim,
        step_size
    ):
        self.dt = dt
        self.st = start_trim
        self.et = end_trim
        self.step_size = step_size
        self.pos = unwrap_pos[start_trim:end_trim:step_size, :, :]

    def autocorr(self):
        """
        Computes autocorrelation function for specified timesteps

        Output
        ------
        S2_AB : float
            (Nt, 3) ndarray of autocorrelation function for the 3 unique tensor
            components in the order (xx, xy, yy)
        """

        Nt = self.pos.shape[0]

        Fx = fft(self.pos[:, :, 0], n=2*Nt, axis=0)
        Fy = fft(self.pos[:, :, 1], n=2*Nt, axis=0)

        PSDxx = 2*Fx.conjugate()*Fx
        PSDxy = Fx.conjugate()*Fy + Fy.conjugate()*Fx
        PSDyy = 2*Fy.conjugate()*Fy

        resxx = ifft(PSDxx, axis=0)
        resxy = ifft(PSDxy, axis=0)
        resyy = ifft(PSDyy, axis=0)

        resxx = (resxx[:Nt]).real
        resxy = (resxy[:Nt]).real
        resyy = (resyy[:Nt]).real

        n = np.arange(1, Nt+1)[::-1]  # divide res(m) by (N-m)

        return resxx/n[:, np.newaxis], resxy/n[:, np.newaxis],\
               resyy/n[:, np.newaxis]

    def displacementTensor(self):
        """
        Parameters
        ----------
        unwrap_pos : float
            Nt x N x 2 ndarray, where Nt is the total number of timesteps and
            N is the number of particles. It is assumed that the positions are
            unwrapped if periodic boundaries are used.
        dt : float
            size of time step
        start_trim : int
            timestep to start computing displacement tensor
        end_trim : int
            last timestep to include in computation
        step_size : int
            Number of timesteps covered in a single step (e.g. if set to 10,
            will take positions every 10 timesteps). Default is 1.

        Returns
        -------
        times : float
            Nt ndarray of times over which MSD is calculated
        disp_tens : float
            2 x 2 ndarray. The symmetric mean square displacement tensor.
            Elements are calculated via the dyadic product of the displacement
            vector and averaged over all particles. Expression is
            < [r(t + \tau) - r(t)][r(t + \tau) - r(t)] >

            The scalar mean square displacement can be computed by taking the
            trace of disp_tens
        """

        N_particles = self.pos.shape[1]
        Nt = self.pos.shape[0]
        self.times = np.arange(self.dt, Nt*self.dt, self.dt)

        # Initiate displacement tensor
        disp_tens = np.zeros((Nt, 2, 2))

        S1xx = s1(self.pos[:, :, 0], self.pos[:, :, 0])
        S1xy = s1(self.pos[:, :, 0], self.pos[:, :, 1])
        S1yy = s1(self.pos[:, :, 1], self.pos[:, :, 1])

        S2xx, S2xy, S2yy = self.autocorr()

        disp_tens[:, 0, 0] = 1/N_particles*(S1xx.sum(axis=1)-S2xx.sum(axis=1))
        disp_tens[:, 0, 1] = 1/N_particles*(S1xy.sum(axis=1)-S2xy.sum(axis=1))
        disp_tens[:, 1, 1] = 1/N_particles*(S1yy.sum(axis=1)-S2yy.sum(axis=1))
        disp_tens[:, 1, 0] = disp_tens[:, 0, 1]

        return disp_tens

    def compute(self):

        return self.displacementTensor()
