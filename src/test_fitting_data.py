# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:23:46 2013

@author: Leif Roschier
"""

"""
Example of fitting data
"""

from CBT_lib import *
from CBT_fit_lib import *

parallel_arrays = 1

# actual data in
meas_V0,meas_R = loadtxt("data/RT_CBT_4.txt", delimiter=None,unpack=True)
#meas_V0,meas_R = loadtxt("data/cbt14.txt", delimiter=None,unpack=True)
meas_G0=1.0/(meas_R*parallel_arrays) #
meas_G0=meas_G0*50.0 # make 100 junctions -> 2 junctions
meas_R = meas_R/50.0 # make 100 junctions -> 2 junctions
meas_V0=meas_V0/50.0 # make 100 junctions -> 2 junctions
# remove values too close to zero, maybe not needed
indices = [x for x, y in enumerate(meas_V0) if (y >1e-9 or y<-1e-9)]

meas_G = meas_G0[indices]
meas_V0 = meas_V0[indices]
meas_R = meas_R[indices]

if False:
    fig1 = figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(meas_V,meas_R,'rx')
    fig1.canvas.draw()
    show()

T_init = 30e-3
island_size_init=400.0
R_tunnel_init=30e3
TEC_init=30e-3

#R,C,T,vol
#bounds = [(1.0,200.0),(1.0,200.0),(5.0,100.0),(0.2,200.0)]
bounds = [(1.0,200.0),(1.0,300.0),(5.0,100.0)]
bounds = None

# actual fitting class generation
fitter=CBT_fitter(filename=None,T_init=T_init,island_size_init=island_size_init,
                  R_tunnel_init=R_tunnel_init,TEC_init=TEC_init,
                  meas_V0=array(meas_V0),meas_R=array(meas_R), bounds=bounds,const_P=500e-18)

fitter.find_offset(5)
fitter.print_offset_curve()
fitter.fit_classic_curve()
fitter.plot_result_initial()
fitter.fit_full_curve()
fitter.print_elapset_time()
fitter.plot_results()

