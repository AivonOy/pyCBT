# -*- coding: utf-8 -*-
"""
Created on Fri Aug 09 12:28:35 2013

@author: Leif
"""


from CBT_lib import *
from CBT_fit_lib import *
from CBT_plot_data_pyx import *

parallel_arrays = 20.0
N_junctions=33.0
excitation = 30e-6/N_junctions


T_init = 10e-3 # K
island_size_init=400.0 #x 1e-15 m3
R_tunnel_init=28.0 # kOhm 
#R_tunnel_init=1.0 # kOhm 
TEC_init=5e-3 # K

#R,C,T,vol
#bounds = [(1.0,200.0),(1.0,200.0),(5.0,100.0),(0.2,200.0)]
#bounds = [(1.0,100.0),(100.0,920.0),(5.0,250.0)]
bounds = [(1.0,100.0),(230.0,240.0),(5.0,250.0)]
#bounds = None
filenames = ["data/test_meas_data.txt"]
fitters=[]
for i,filename in enumerate(filenames):
    fitter=CBT_fitter(filename=filename,T_init=T_init,island_size_init=island_size_init,
                      R_tunnel_init=R_tunnel_init,TEC_init=TEC_init,
                        bounds=bounds,const_P=500e-18,
                        excitation=excitation,parallel_arrays=parallel_arrays,
                        junctions_in_series = N_junctions)
    fitter.find_offset(5)
    fitter.print_offset_curve()
    fitter.fit_classic_curve()
    fitter.plot_result_initial()
    fitter.fit_full_curve()
    fitter.print_elapset_time()
    fitter.plot_nonlin_results()
    fitters.append(fitter)
#fitter.plot_all_results()
cc = CBT_plot_data_pyx(fitters)
cc.plot_data(filename="res_8p5_50_100")


