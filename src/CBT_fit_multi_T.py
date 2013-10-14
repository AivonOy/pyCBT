# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 09:51:03 2013

@author: Leif Roschier/Aivon Oy

This class fits single C_sigma, R_T from list of already fitted fitted curves
inside classes CBT_fitter (from file CBT_fit_lib.py).
"""
from scipy import optimize
from CBT_lib import *

class CBT_multi_fitter():
    def __init__(self,fitters,bounds=None):
        """
        fitters = list of classes CBT_fitter
        """
        self.fitters=fitters
        self.bounds = bounds
        self.fit_curves() # do actual fitting
        
        
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
        if self.bounds == None:
            self.xopt = optimize.fmin_bfgs(self.fit_function, x1, gtol=1e-4,full_output=1, disp=1,callback=self.call_func)
        else:
            "Print optimizing with bounds"
            self.xopt = optimize.fmin_l_bfgs_b(self.fit_function, x1, factr=1e7, approx_grad=True, bounds=self.bounds)

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
    
    def call_func(self,x): 
        """
        can be used to print information each fit cycle
        """
        #print x
        pass
    