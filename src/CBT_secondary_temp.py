# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 15:45:21 2013

@author: Leif Roschier/Aivon Oy


This class CBT_secondary_temp is used as interpolator to find temperature from
conductance minimum.

List of resistances is calculated from temperatures and interpolation function 
is generated. 
"""

from CBT_lib import *
from scipy import interpolate
from numpy import *
import matplotlib.pyplot as plt

class CBT_secondary_temp():
    def __init__(self,C_sigma,R_T,sigma=0.2e9,island_vol=200e-15,
                 T_min=5e-3,T_max=200e-3,T_step=0.5e-3,
                 N_parallel=20, N_series=33,
                 const_P=0.5e-15, excitation=10e-6):
        """
        C_sigma: sum capacitance per island
        R_T: tunnelling resistance
        sigma: electron-phonon coupling
        island_vol: island volume
        T_min, T_max, T_step: params for interpolating function calc
        N_parallel: junctions in parallel
        N_series: junctions in series
        const_P: constant power heating per junction
        excitation: measurement AC voltage, used in numerical derivative
        """
        G=[]
        R=[]
        T=[]
        T_list_1=arange(5.0e-3,20.0e-3,0.5e-3)
        T_list_2=arange(21.0e-3,100.0e-3,1.0e-3)
        T_list_3=arange(100.0e-3,200.0e-3,5e-3)
        T_list_4=arange(200.0e-3,1.0,50e-3)
        T_list_5=arange(1.0,10.0,0.5)
        T_list=concatenate((T_list_1,T_list_2,T_list_3,T_list_4,T_list_5))
        #for t in arange(T_min, T_max+T_step,T_step):
        for t in T_list:
            T.append(t)
            #sigma,N,V,R_T,C_sigma,T_p,island_volume, const_P, eps=1e-9
            this_G = (calc_G(sigma,N_series,1e-9,R_T,C_sigma,t,island_vol,
                                const_P,excitation/N_series))*N_parallel
            #print "T:%g,R.%g"%(t,1.0/this_G)
            G.append(this_G)
            R.append(1.0/this_G)
        #z = polyfit(G, T, 20)
        #p = poly1d(z)
        #GG = p(T)
        if 0:
             plt.figure(8)
             #plt.plot(T,G,'ro',T,GG,'b-')
             plt.plot(T,G,'ro')
        R.reverse()
        T.reverse()
        if 0:
            plt.figure(18)
            plt.plot(T,R,'ro')
        self.T_func = interpolate.interp1d(array(R), array(T),bounds_error=False)
        
if __name__ == '__main__': 
    C_sigma = 231.033e-15
    R_T = 25.4083e3
    sigma=0.2e9
    island_vol=200e-15
    cbt_secondary = CBT_secondary_temp(C_sigma,R_T)
    R_test=48687.271436
    print "temp at %g kOhm = %g mK"%(R_test,cbt_secondary.T_func(R_test)*1000)       
        