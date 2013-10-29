# -*- coding: utf-8 -*-
"""

@author: Leif Roschier/Aivon Oy

This class fits single C_sigma, R_T from a given curve.
"""
from matplotlib import mlab
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from pylab import figure, show
from scipy import optimize
from scipy import interpolate
from scipy.constants import *
from CBT_lib import *
from copy import *
import time

N=2.0
#n_max=20
class CBT_fitter():
    """ 
    Class for fitting measured CBT data
    """
    def __init__(self,filename="data/cbt14.txt", T_init=50e-3, island_size_init=1.0,
                 R_tunnel_init=30,TEC_init=100e-3, n_max=200, N=2.0,
                 sigma=0.2e9,meas_V0=[],meas_R=[],bounds=None,v_offset=0.0,const_P=1e-18,show_first_fit=True,
                excitation=1e-9,parallel_arrays=20.0, junctions_in_series = 33.0):
        """
            Give data as link to file or with arrays meas_V0 and meas_R with filename=None
            T_init: initial guess for temperature in (K)            
            island size: island size in 1e-15 (m**3)
            R_tunnel_init: tunnelling resistance in (kOhm) 
            TEC_init: charging energy in (K)
            excitation: value for calculating conductance
        """
        self.tic = time.clock()
        self.parallel_arrays = parallel_arrays
        self.junctions_in_series= junctions_in_series
        self.sigma =sigma
        self.n_max=n_max
        self.N=N
        self.v_offset=0.0
        self.const_P=const_P
        self.island_size = island_size_init*1e-15
        self.R_tunnel_init = R_tunnel_init # kOhm
        self.T_init = T_init
        self.TEC_init = TEC_init
        self.bounds = bounds
        self.excitation = excitation
        self.filename=filename
        if not filename==None:
            meas_V0,meas_R = loadtxt(filename, delimiter=None,unpack=True)
        # make chain equal to 2 junctions
        meas_R = meas_R/(junctions_in_series/2.0)
        meas_G0=1.0/(meas_R*parallel_arrays) #
        meas_V0 = meas_V0/(junctions_in_series/2.0)
        # remove voltages around zero
        indices = [x for x, y in enumerate(meas_V0) if (y >1e-8 or y<-1e-8)]
        self.meas_G= meas_G0[indices]
        self.meas_V= meas_V0[indices]
        self.meas_R = meas_R[indices]
        self.v_offset = v_offset = 0.0
        #print "Initial v_offset %g"%v_offset
        #print "Const p:%g"%self.const_P
        self.classic_curve_fitted = False
        self.full_curve_fitted = False        
        self.v_offset =0.0

    def find_offset(self,points=None):
        """
        finds offset voltage
        """
        #find minimum conductance
        G = array(self.meas_G)
        min_G_idx = G.argmin()
        if points==None:
            num_points = G.size/10
        else:
            num_points = points
        min_idx = min_G_idx-num_points
        max_idx = min_G_idx+num_points
        if min_idx<0:
            min_idx=0
        if max_idx>G.size:
            max_idx = G.size
        #print "min_idx:%g"%min_idx
        #print "max_idx:%g"%max_idx
        indices = range(min_idx,max_idx+1)
        self.offset_data={}
        self.offset_data['x'] = x = array(self.meas_V[indices])
        self.offset_data['y'] = y = array(self.meas_G[indices])
        self.offset_data['z'] = z = polyfit(x, y, 2)
        self.offset_data['p'] = p = poly1d(z)
        self.offset_data['y2'] = p(x)
        self.v_offset = -z[1]/(2.0*z[0])
        print "offset:%g"%self.v_offset
        self.meas_V = self.meas_V-self.v_offset

    def print_offset_curve(self):
        """
        plots offset calculation
        """
        plt.figure(8)
        plt.plot( (self.offset_data['x']-self.v_offset)*1e3,self.offset_data['y']*1e6,'ro',
                 self.meas_V*1e3, self.meas_G*1e6,'-x',
                 (self.offset_data['x']-self.v_offset)*1e3,self.offset_data['y2']*1e6,'bo')
        #plt.plot( meas_V,peval( meas_V,plsq[0])/plsq[0][0], meas_V, meas_G/plsq[0][0],'--')
        plt.title('Offset corrected curve')
        plt.legend(['Fit range', 'Measurement', 'Fit'])
        plt.xlabel('voltage (mV)')
        plt.ylabel(r'conductance ($\mu S$)')
        plt.show()

    def print_elapset_time(self):
        toc = time.clock()
        print "It has taken %g seconds = %g minutes"%(toc-self.tic,(toc-self.tic)/60.0)

    def fit_classic_curve(self):
        # initial values for "classic" curve fit                        
        p0=[1.0/(self.R_tunnel_init),self.T_init,self.TEC_init]
        # initial optimization
        self.plsq=plsq = optimize.fmin_bfgs(self.residuals1, p0, args=(self.meas_G,self.meas_V),
                                                gtol=1e-6, full_output=1,maxiter=50,
                                                callback=self.call_func)          
        print "==== After initial optimization: ===="     
        print "R_T = %g kOhm"%(1.0/(plsq[0][0]*self.N))
        print "T = %g mK"%(plsq[0][1]*1000.0)
        print "Ec = %g mK"%(plsq[0][2]*(N-1.0)/self.N*1000)
        print "Csigma = %g fF"%(e**2/(plsq[0][2]*k)*1e15)
        self.classic_curve_fitted = True

    def fit_full_curve(self):
        # main optimization
        if self.classic_curve_fitted==True:
            R_T_init = (1.0/(self.plsq[0][0]*self.N))
            C_sigma_init=e**2/(self.plsq[0][2]*k)*1e15
            T_p_init = self.plsq[0][1]*1e3
        else:
            R_T_init = self.R_tunnel_init
            C_sigma_init = e**2/(self.TEC_init*k)*1e15
            T_p_init = self.T_init*1e3
        island_volume_init = self.island_size
        x1=[R_T_init,C_sigma_init,T_p_init]
        if self.bounds == None:
            self.xopt1 = optimize.fmin_bfgs(self.optimize_1, x1, gtol=1e-3,full_output=1, disp=1,callback=self.call_func)
        else:
            "Print optimizing with bounds"
            self.xopt1 = optimize.fmin_l_bfgs_b(self.optimize_1, x1, factr=1e7, approx_grad=True, bounds=self.bounds)
        toc = time.clock()
        print "==== After main optimization: ===="     
        print "R_T = %g"%(self.xopt1[0][0])
        print "T = %g mK"%(self.xopt1[0][2])
        print "C_sigma = %g "%(self.xopt1[0][1])
        self.full_curve_fitted = True
        self.T_fit = self.xopt1[0][2]
        self.R_T = self.xopt1[0][0]
        self.C_sigma = self.xopt1[0][1]

    def plot_result_initial(self):
        meas_V = self.meas_V.tolist()
        meas_G = self.meas_G.tolist()
        #plt.plot( self.meas_V,self.peval( self.meas_V,self.plsq[0]),'r', self.meas_V, self.meas_G,'--',self.meas_V,self.peval_1( self.meas_V,self.xopt1[0]),'g--')
        plt.figure(9)        
        plt.plot( self.meas_V*1e3,array(self.peval(self.meas_V,self.plsq[0]))*1e6,'r',
                 array(meas_V)*1e3, array(meas_G)*1e6,'-x')    
        plt.title('Fit to analytic curve')
        plt.xlabel('voltage (mV)')
        plt.ylabel(r'conductance ($\mu S$)')
        plt.legend(['Measurement', 'Fit'])
        plt.show()
    
    def plot_all_results(self):
        meas_V = self.meas_V.tolist()
        meas_G = self.meas_G.tolist()
        #plt.plot( self.meas_V,self.peval( self.meas_V,self.plsq[0]),'r', self.meas_V, self.meas_G,'--',self.meas_V,self.peval_1( self.meas_V,self.xopt1[0]),'g--')
        plt.figure(10)        
        plt.plot( self.meas_V,self.peval( self.meas_V,self.plsq[0]),'r', meas_V, meas_G,'-x',meas_V,self.peval_1( meas_V,self.xopt1[0]),'g--')
        #plt.plot( meas_V,peval( meas_V,plsq[0])/plsq[0][0], meas_V, meas_G/plsq[0][0],'--')
        plt.title('Least-squares fit numerical model')
        plt.legend(['Fit1', 'Measurement', 'Fit2'])
        plt.show()
        
    def plot_nonlin_results(self):
        meas_V = self.meas_V.tolist()
        meas_G = self.meas_G.tolist()
        #plt.plot( self.meas_V,self.peval( self.meas_V,self.plsq[0]),'r', self.meas_V, self.meas_G,'--',self.meas_V,self.peval_1( self.meas_V,self.xopt1[0]),'g--')
        plt.figure(11)        
        plt.plot( array(meas_V)*1e3, array(meas_G)*1e6,'-x',
                 array(meas_V)*1e3,array(self.peval_1( meas_V,self.xopt1[0]))*1e6,'g--')
        #plt.plot( meas_V,peval( meas_V,plsq[0])/plsq[0][0], meas_V, meas_G/plsq[0][0],'--')
        plt.title('Nonlinear fit results')
        plt.legend([ 'Measurement', 'nonlin fit'],loc=4)
        plt.xlabel('voltage (mV)')
        plt.ylabel(r'conductance ($\mu S$)')
        plt.show()        
        
    def plot_start(self,x):
        meas_V = self.meas_V.tolist()
        meas_G = self.meas_G.tolist()
        #plt.plot( self.meas_V,self.peval( self.meas_V,self.plsq[0]),'r', self.meas_V, self.meas_G,'--',self.meas_V,self.peval_1( self.meas_V,self.xopt1[0]),'g--')
        plt.figure(12)
        plt.plot( self.meas_V,self.peval( self.meas_V,self.plsq[0]),'r', meas_V, meas_G,'-x',meas_V,self.peval_1( meas_V,x),'g--')
        #plt.plot( meas_V,peval( meas_V,plsq[0])/plsq[0][0], meas_V, meas_G/plsq[0][0],'--')
        plt.title('Least-squares fit numerical model')
        plt.legend(['Fit1', 'Measurement', 'Fit2'],loc=4)
        plt.show()        
    # FIRST OPTIMIZATION
    def g(self,x):
        return (x*sinh(x)-4*sinh(x/2.0)**2)/(8.0*sinh(x/2.0)**4)

    def G_optimize_func(self,x,V):
        #N = 2
        Gt_opt = x[0]
        T_opt = x[1]
        EC_opt = x[2]
        return Gt_opt*(1.0-2*(self.N-1.0)/self.N*EC_opt/T_opt*g(e*V/(self.N*k*T_opt)))
    
    def residuals(self,p, y, V):
         #N = 2
         Gt_opt,T_opt,EC_opt = p
         err = y-Gt_opt*(1.0-2*(self.N-1.0)/self.N*EC_opt/T_opt*self.g(e*V/(self.N*k*T_opt)))
         return err
     
    def residuals1(self,p, meas_G, meas_V):
         #N = 2
        Gt_opt,T_opt,EC_opt = p
        meas_G = meas_G*1000.0; # make into kOhm
        total_error = 0.0
        for idx,V in enumerate(meas_V):
            #print "idx %g"%idx
            total_error=total_error + (meas_G[idx]-Gt_opt*(1.0-2*(self.N-1.0)/self.N*EC_opt/T_opt*self.g(e*V/(self.N*k*T_opt))))**2
        return total_error*1e5

    def peval(self,V, p):
        #N = 2
        Gt_opt,T_opt,EC_opt = p[0],p[1],p[2]
        Gt_opt = Gt_opt/1000.0
        return Gt_opt*(1.0-2*(self.N-1.0)/self.N*EC_opt/T_opt*self.g(e*(V)/(self.N*k*T_opt)))
    
    def call_func(self,x): 
        print x

    # Full optimization

    def fit_func(self,x,R_T,C_sigma,T_p,island_volume): #calc_conductance_curve_full_spline(sigma,N,meas_V,R_T,C_sigma,T_p,island_volume)
        return calc_conductance_curve_full_spline(sigma,self.N,x,R_T*1e3,C_sigma*1e-15,T_p,island_volume*1e-15)

    def optimize_1(self,x):
        R_T = x[0]*1e3
        C_sigma = x[1]*1e-15
        T_p = x[2]*1e-3
        #island_volume = abs(x[3])*1e-15
        G = []
        for V in self.meas_V:
            #calc_G(sigma,N,V,R_T,C_sigma,T_p,island_volume, const_P, eps=1e-9):
            G.append(calc_G(self.sigma,self.N,V,R_T,C_sigma,T_p,self.island_size,
                            self.const_P,eps=self.excitation))
        res = sum((array(self.meas_G)-array(G))**2)*1e10
        print "Rt:%g Csigma:%g Tp:%g res:%g"%(R_T,C_sigma,T_p,res)
        return res

    def peval_1(self,meas_V,x):
        R_T = x[0]*1e3
        C_sigma = x[1]*1e-15
        T_p = x[2]*1e-3
        #island_volume = x[3]*1e-15
        calc_G_few = calc_conductance_curve_basic(self.sigma,self.N,meas_V,R_T,C_sigma,T_p,self.island_size)
        return calc_G_few

    def original_V(self):
        #gives original data
        return self.meas_V.tolist()
        
    def original_G(self):
        #gives original data
        return self.meas_G.tolist()
        
    def nonlinfit_G(self):
        #gives original data
        nonlinfit_G = self.peval_1(self.meas_V,self.xopt1[0])
        return nonlinfit_G
        
    def G_curve(self,R_T0,C_sigma0,T_0):
        R_T = R_T0*1e3
        C_sigma = C_sigma0*1e-15
        T_p = T_0*1e-3
        calc_G = calc_conductance_curve_basic(self.sigma,self.N,self.meas_V,R_T,C_sigma,T_p,self.island_size)
        return calc_G
        
        

if __name__ == '__main__': 
    
    fitter=CBT_fitter()
    fitter.plot_results()


