# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 09:51:03 2013

@author: Leif Roschier/Aivon Oy

This class fits single C_sigma, R_T from list of already fitted fitted curves
inside classes CBT_fitter (from file CBT_fit_lib.py).
"""
from scipy import optimize
from CBT_lib import *
from multiprocessing import Pool,Queue,Value, Process

class CBT_multi_fitter():
    def __init__(self,fitters,bounds=None,
                 result_text_filename="mfit_results.txt",
                 multiprocess=False,t_data=None,
                 interpolation_table_file="mfit_interpolation_table.txt"):
        """
        fitters = list of classes CBT_fitter
        """
        self.fitters=fitters
        self.bounds = bounds
        self.multiprocess=multiprocess
        self.fit_curves() # do actual fitting
        self.print_results(result_file=result_text_filename)
        self.built_interpolation_table(t_data=t_data)
        self.save_interp_table_to_file(interpolation_table_file)
        
    def built_interpolation_table(self,t_data=None):
        """
        builds interpolation curve
        """
        G_list=[] # conductances
        R_list=[] # resistances, inverse of conductances
        T_list=[] # list of temperatures
        if t_data==None: # no explicitely given temperatures
            T_list_1=arange(5.0e-3,20.0e-3,0.5e-3)
            T_list_2=arange(21.0e-3,100.0e-3,1.0e-3)
            T_list_3=arange(100.0e-3,200.0e-3,5e-3)
            T_list=concatenate((T_list_1,T_list_2,T_list_3)).tolist()
        else:
            T_list=t_data
        T_list.reverse() # make resistances to appear ascending
        fitter=self.fitters[0] # should be at least one, all have same _multi values
        for t in T_list:
            #sigma,N,V,R_T,C_sigma,T_p,island_volume, const_P, eps=1e-9
            #calc_G(fitter.sigma,fitter.N,V,R_T,C_sigma,T_p,fitter.island_size,
            #                    fitter.const_P,eps=fitter.excitation)
            this_G = (calc_G(fitter.sigma,fitter.junctions_in_series,1e-9,fitter.R_T_multi*1e3,
                             fitter.C_sigma_multi*1e-15,t,fitter.island_size,
                                fitter.const_P,eps=fitter.excitation))*fitter.parallel_arrays
            print "T:%g,R.%g"%(t,1.0/this_G)
            G_list.append(this_G)
            R_list.append(1.0/this_G)
        self.T_list=T_list
        self.G_list=G_list
        self.R_list=R_list
        # interpolating function
        self.T_func = interpolate.interp1d(array(R_list), array(T_list),bounds_error=False)

    def save_interp_table_to_file(self,filename="mfit_interpolation_table.txt"):
        """
        saves interpolation curve to file
        """
        fo = open(filename, "wb")
        for R,T in zip(self.T_list,self.R_list):
             fo.write("%g\t%g\n"%(R,T))            
        fo.close

    def fit_curves(self):
        """
        fits all curves
        """
        # Take average of fit params
        self.T_fit = 0.0
        self.R_T = 0.0
        self.C_sigma = 0.0
        n=0.0
        for fitter in self.fitters:
            self.R_T=self.R_T+fitter.R_T
            self.C_sigma=self.C_sigma+fitter.C_sigma
            n=n+1
        self.R_T=self.R_T/n
        self.C_sigma=self.C_sigma/n
        print "Average C_sigma:%g"%self.C_sigma
        print "Average R_T:%g"%self.R_T
        # actual fitting
        x1=[self.R_T,self.C_sigma] # initial conditions
        for fitter in self.fitters: # add temperatures to initial conditions
            x1.append(fitter.T_fit)
        print "Initial vector"
        print x1
        # check if multiprocessing is used, does not work in Windows
        if self.multiprocess:
            fit_fcn = self.fit_function_multiprocessing
        else:
            fit_fcn = self.fit_function
        if self.bounds == None:
            self.xopt = optimize.fmin_bfgs(fit_fcn,
                                           x1, gtol=1e-4,full_output=1,
                                           disp=1,callback=self.call_func)                
        else:
            "Print optimizing with bounds"
            self.xopt = optimize.fmin_l_bfgs_b(fit_fcn, x1, factr=1e7, approx_grad=True, bounds=self.bounds)

        print "xopt:"
        print self.xopt
        print "==== After multi optimization: ===="
        print "R_T = %g"%(self.xopt[0][0])
        print "C_sigma = %g "%(self.xopt[0][1])
        for idx,T in enumerate(self.xopt[0]):
            if idx>1:
                print "T%i = %g mK"%((idx-1),T)
        # save results
        for idx,fitter in enumerate(self.fitters):
            fitter.R_T_multi = self.xopt[0][0]
            fitter.C_sigma_multi = self.xopt[0][1]
            fitter.T_multi = self.xopt[0][idx+2]


    def fit_function(self,x):
        """
        function to be minimized
        """
        R_T = x[0]*1e3
        C_sigma = x[1]*1e-15
        res=0.0
        for idx,fitter in enumerate(self.fitters): # error for every curve
            T_p=x[2+idx]*1e-3
            G = []
            for V in fitter.meas_V: # error for every V point
                #calc_G(sigma,N,V,R_T,C_sigma,T_p,island_volume, const_P, eps=1e-9):
                G.append(calc_G(fitter.sigma,fitter.N,V,R_T,C_sigma,T_p,fitter.island_size,
                                fitter.const_P,eps=fitter.excitation))
            res = res+sum((array(fitter.meas_G)-array(G))**2)*1e10
            print "Idx:%g,Rt:%g Csigma:%g Tp:%g res:%g"%(idx,R_T,C_sigma,T_p,res)
        return res

    def fit_function_pool(self,x):
        """
        function to be minimized, DOES NOT WORK
        """
        R_T = x[0]*1e3
        C_sigma = x[1]*1e-15
        res=0.0
        pool = Pool(processes=8) 
        for idx,fitter in enumerate(self.fitters): # error for every curve
            T_p=x[2+idx]*1e-3
            G = []
            def G_func_local(V):
                res = calc_G(fitter.sigma,fitter.N,V,R_T,C_sigma,T_p,fitter.island_size,
                                fitter.const_P,eps=fitter.excitation)
                return res
            G = pool.map(G_func_local, fitter.meas_V)
            res = res+sum((array(fitter.meas_G)-array(G))**2)*1e10
            print "Idx:%g,Rt:%g Csigma:%g Tp:%g res:%g"%(idx,R_T,C_sigma,T_p,res)
        return res

    def fit_function_multiprocessing(self,x):
        """
        minimized funtio with multiprocessing
        """
        R_T = x[0]*1e3
        C_sigma = x[1]*1e-15
        results=[]
        procs=[]
        values = []
        for idx,fitter in enumerate(self.fitters): # error for every curve
            T_p=x[2+idx]*1e-3
            results.append(0.0)
            values.append(Value('d', 0.0))
            p = Process(
                target=self.fit_G_worker,
                args=(idx,fitter,R_T,C_sigma,T_p,values[-1]))
            procs.append(p)
            p.start()
        for p in procs:
            p.join(10)
        total_error=0.0
        for error in values:
            total_error = total_error + error.value
            print "error:%g"%error.value
        return total_error

    def fit_G_worker(self,idx,fitter,R_T,C_sigma,T_p,result_error):
        """
        worker for multiprocessing
        """
        G = []
        for V in fitter.meas_V: # error for every V point
        #calc_G(sigma,N,V,R_T,C_sigma,T_p,island_volume, const_P, eps=1e-9):
            G.append(calc_G(fitter.sigma,fitter.N,V,R_T,C_sigma,T_p,fitter.island_size,
                        fitter.const_P,eps=fitter.excitation))
        result_error.value = sum((array(fitter.meas_G)-array(G))**2)*1e10
        print "Idx:%g,Rt:%g Csigma:%g Tp:%g res:%g"%(idx,R_T,C_sigma,T_p,result_error.value)
        return


    def call_func(self,x):
        """
        can be used to print information each fit cycle
        """
        #print x
        pass

    def print_results(self,result_file="dummy_results.txt"):
        """
        Prints fit results, original and multi
        """
        fo = open(result_file, "wb")
        print "===================="
        print "     Final results  "
        for idx,fitter in enumerate(self.fitters):
            print "file:%s"%fitter.filename
            fo.write("file:%s \n"%fitter.filename)
            print "R_T(0):%g R_T(multi):%g"%(fitter.R_T,fitter.R_T_multi)
            fo.write("R_T(0):%g R_T(multi):%g \n"%(fitter.R_T,fitter.R_T_multi))
            print "C_sigma(0):%g C_sigma(multi):%g"%(fitter.C_sigma,fitter.C_sigma_multi)
            fo.write("C_sigma(0):%g C_sigma(multi):%g \n"%(fitter.C_sigma,fitter.C_sigma_multi))
            print "T(0):%g T(multi):%g"%(fitter.T_fit,fitter.T_multi)
            fo.write("T(0):%g T(multi):%g \n"%(fitter.T_fit,fitter.T_multi))

        # Open a file
        fo.close()

        # Close opend file


