# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 13:56:30 2013

@author: Leif Roschier/Aivon Oy
"""


"""
This is example script to fit CBT R-V curves at different temperatures
simultaneously in order to find C_sigma and R_T.

Script fits first all separate curves with class CBT_fitter to get initial 
values for final fitting using class CBT_multi_fitter.
"""

from CBT_lib import *
from CBT_fit_lib import *
from CBT_plot_data_pyx import *
from CBT_fit_multi_T import *
import time

start = time.time() # start taking time

parallel_arrays = 20.0  # arrays parallel    
N_junctions=33.0        # junctions in each array
excitation = 30e-6/N_junctions  # excitation of measurement
island_size_init=400.0 #x 1e-15 m3, island size for heating calculations

# initial values for optimization
T_init = 10e-3 # K, initial guess for temperature
R_tunnel_init=28.0 # kOhm 
TEC_init=5e-3 # K

# If bounds given, they are used to limit optimization
#R,C,T
bounds = [(1.0,100.0),(230.0,240.0),(5.0,250.0)] #(min,max)
#bounds = None 

# actual measurement files to be fitted
filenames = [
"data/sample_T2.txt",
"data/sample_T3.txt",
"data/sample_T4.txt",
"data/sample_T5.txt",
"data/sample_T6.txt",
]
"""
filenames = [
"data/sample_T2.txt",
"data/sample_T5.txt",
]
"""

results_txt_file = "mfit_results_41A4_t.txt" # C_sigma and R_T of fit
result_pdf_1="41A4_multifit1.pdf" # conductance curves and fits
result_pdf_2="41A4_multifit2.pdf"
interpolation_filename="41A4_interpolation_table.txt"


# separate fitting of each curve
fitters=[] # list of classes CBT_fitter that also store results
for i,filename in enumerate(filenames):
    fitter=CBT_fitter(filename=filename,T_init=T_init,island_size_init=island_size_init,
                      R_tunnel_init=R_tunnel_init,TEC_init=TEC_init,
                        bounds=bounds,const_P=500e-18,
                        excitation=excitation,parallel_arrays=parallel_arrays,
                        junctions_in_series = N_junctions)
    fitter.find_offset(5) # finds offset with 5 points from minimum each direction
    #fitter.print_offset_curve()
    fitter.fit_classic_curve() # fit analytic curve to get starting point for nonlin fit
    #fitter.plot_result_initial()
    fitter.fit_full_curve() # nonlin fitting
    fitter.print_elapset_time()
    #fitter.plot_nonlin_results()
    fitters.append(fitter) # save results
#fitter.plot_all_results()
    
# fitting all together
multi_fitter = CBT_multi_fitter(fitters,
                                result_text_filename=results_txt_file,
                                interpolation_table_file=interpolation_filename)

# plot results    
cc = CBT_plot_data_pyx(fitters)
cc.plot_multi_data(filename=result_pdf_1)
cc.plot_data(filename=result_pdf_2)
print "it took", time.time() - start, "seconds."