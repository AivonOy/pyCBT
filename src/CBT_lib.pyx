# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 12:09:04 2011

@author: Leif
"""
# LRR 23.2.2011

from numpy import *
from scipy import *
from scipy.constants import *
from scipy.interpolate import UnivariateSpline
from scipy import interpolate

n_max=150

cdef double fii(double n,double C_sigma):
    return n*e/C_sigma

cdef double DF1plus(double n,double V,double C_sigma):
    return (e/2.0)*(fii(n,C_sigma)+fii(n+1,C_sigma)-(V))

cdef double DF1minus(double n,double V,double C_sigma):
    return (e/2.0)*(V-(fii(n,C_sigma)+fii(n-1,C_sigma)))

cdef double DF2plus(double n,double V,double C_sigma):
    return (e/2.0)*(-V-(fii(n,C_sigma)+fii(n-1,C_sigma)))

cdef double DF2minus(double n,double V,double C_sigma):
    return (e/2.0)*((fii(n,C_sigma)+fii(n+1,C_sigma))+V)

cdef double gamma1plus(double n,double V,double T,double R_T,double C_sigma):
    if (abs(DF1plus(n,V,C_sigma)/(k*T))>1e-12):
        return -1.0/(e**2*R_T)*DF1plus(n,V,C_sigma)/(1.0-exp(DF1plus(n,V,C_sigma)/(k*T)))
    else:
        return -k*T/(e**2*R_T)*1.0

cdef double gamma1minus(double n,double V,double T,double R_T,double C_sigma):
    if (abs(DF1minus(n,V,C_sigma)/(k*T))>1e-12):
        return -1.0/(e**2*R_T)*DF1minus(n,V,C_sigma)/(1.0-exp(DF1minus(n,V,C_sigma)/(k*T)))
    else:
        return -k*T/(e**2*R_T)*1.0
    #return -1.0/(e**2*R_T)*DF1minus(n,V,C_sigma)/(1.0-exp(DF1minus(n,V,C_sigma)/(k*T)))

cdef double gamma2plus(double n,double V,double T,double R_T,double C_sigma):
    if (abs(DF2plus(n,V,C_sigma)/(k*T))>1e-12):
        return -1.0/(e**2*R_T)*DF2plus(n,V,C_sigma)/(1.0-exp(DF2plus(n,V,C_sigma)/(k*T)))
    else:
        return -k*T/(e**2*R_T)*1.0
    #return -1.0/(e**2*R_T)*DF2plus(n,V,C_sigma)/(1.0-exp(DF2plus(n,V,C_sigma)/(k*T)))

cdef double gamma2minus(double n,double V,double T,double R_T,double C_sigma):
    #print abs(DF2minus(n,V,C_sigma)/(k*T))
    if (abs(DF2minus(n,V,C_sigma)/(k*T))>1e-12):
        return -1.0/(e**2*R_T)*DF2minus(n,V,C_sigma)/(1.0-exp(DF2minus(n,V,C_sigma)/(k*T)))
    else:
        return -k*T/(e**2*R_T)*1.0
    #return -1.0/(e**2*R_T)*DF2minus(n,V,C_sigma)/(1.0-exp(DF2minus(n,V,C_sigma)/(k*T)))

cdef double calc_current_1(double V,double T,double R_T,double C_sigma):
    """
    calculates current when voltage V is over 2 (two) junctions
    """
    sigma=array([1.0])
    sigma2 = 1.0/((gamma1plus(-1,V,T,R_T,C_sigma)+gamma2minus(-1,V,T,R_T,C_sigma)+gamma1minus(1,V,T,R_T,C_sigma)+gamma2plus(1,V,T,R_T,C_sigma))/(gamma1plus(0,V,T,R_T,C_sigma)+gamma1minus(0,V,T,R_T,C_sigma)+gamma2plus(0,V,T,R_T,C_sigma)+gamma2minus(0,V,T,R_T,C_sigma)))
    sigma=append(sigma,sigma2)
    #print sigma
    n_real=1
    ref_current=1 # will be calculated below
    stop_calculation = False # stop further calculation
    for n in range(1,n_max+1):
        if stop_calculation==False:
            sigma2 = -((gamma1plus(n-1,V,T,R_T,C_sigma)+gamma2minus(n-1,V,T,R_T,C_sigma))*sigma[-2]-(gamma1plus(n,V,T,R_T,C_sigma)+gamma1minus(n,V,T,R_T,C_sigma)+gamma2plus(n,V,T,R_T,C_sigma)+gamma2minus(n,V,T,R_T,C_sigma))*sigma[-1])/(gamma1minus(n+1,V,T,R_T,C_sigma)+gamma2plus(n+1,V,T,R_T,C_sigma))
            sigma=append(sigma,sigma2)
            n_real=n_real+1
            if n==1:
                ref_curr = sigma[0]*(gamma1plus(0,V,T,R_T,C_sigma)-gamma1minus(0,V,T,R_T,C_sigma))
            if n>1:
                current_curr = sigma[n]*(gamma1plus(n,V,T,R_T,C_sigma)-gamma1minus(n,V,T,R_T,C_sigma))
                print sigma                
                print "N:%g current_curr:%g ref_curr:%g"%(n,current_curr,ref_curr)                
                if abs(current_curr/ref_curr)<1e-10: # criteria for stopping calculation
                    stop_calculation=True
                    print "stopped calculation at n = %g"%n
            if False: # for debug
                print "------N = %i------"%n
                print gamma1plus(n,V,T,R_T,C_sigma)
                print gamma1minus(n,V,T,R_T,C_sigma)
                print gamma2plus(n,V,T,R_T,C_sigma)
                print gamma2minus(n,V,T,R_T,C_sigma)
        sigma_sum = sigma[0]+2*sum(sigma[1:])
    #print sigma_sum
    sigma=sigma/sigma_sum
    #print sigma
    current = sigma[0]*(gamma1plus(0,V,T,R_T,C_sigma)-gamma1minus(0,V,T,R_T,C_sigma))
#    for n in range(1,n_max+1):
    for n in range(1,n_real+1):
        current = current + sigma[n]*(gamma1plus(n,V,T,R_T,C_sigma)-gamma1minus(n,V,T,R_T,C_sigma))
        current = current + sigma[n]*(gamma1plus(-n,V,T,R_T,C_sigma)-gamma1minus(-n,V,T,R_T,C_sigma))
    current = current*e
    return current
    #print current
    #print current*2*R_T
    
# tests
cdef double calc_current(V,T,R_T,C_sigma):
    """
    calculates current when voltage V is over 2 (two) junctions
    """
    sigma=array([1.0])
    #sigma[1] = 1.0/((gamma1plus(-1,V,T,R_T,C_sigma)+gamma2minus(-1,V,T,R_T,C_sigma)+gamma1minus(1,V,T,R_T,C_sigma)+gamma2plus(1,V,T,R_T,C_sigma))/(gamma1plus(0,V,T,R_T,C_sigma)+gamma1minus(0,V,T,R_T,C_sigma)+gamma2plus(0,V,T,R_T,C_sigma)+gamma2minus(0,V,T,R_T,C_sigma)));
    sigma2 = 1.0/((gamma1plus(-1,V,T,R_T,C_sigma)+gamma2minus(-1,V,T,R_T,C_sigma)+gamma1minus(1,V,T,R_T,C_sigma)+gamma2plus(1,V,T,R_T,C_sigma))/(gamma1plus(0,V,T,R_T,C_sigma)+gamma1minus(0,V,T,R_T,C_sigma)+gamma2plus(0,V,T,R_T,C_sigma)+gamma2minus(0,V,T,R_T,C_sigma)))
    sigma=append(sigma,sigma2)
    #print sigma
    #ref_current=1 # will be calculated below
    ref_curr = sigma[0]*(gamma1plus(0,V,T,R_T,C_sigma)-gamma1minus(0,V,T,R_T,C_sigma))
    current = ref_curr + sigma[1]*((gamma1plus(1,V,T,R_T,C_sigma)-gamma1minus(1,V,T,R_T,C_sigma))+(gamma1plus(-1,V,T,R_T,C_sigma)-gamma1minus(-1,V,T,R_T,C_sigma)))
    stop_calculation = False # stop further calculation
    #if False:
    for n in range(2,n_max+1):
        if stop_calculation==False:
            #sigma2 = -((gamma1plus(n-1,V,T,R_T,C_sigma)+gamma2minus(n-1,V,T,R_T,C_sigma))*sigma[-2]-(gamma1plus(n,V,T,R_T,C_sigma)+gamma1minus(n,V,T,R_T,C_sigma)+gamma2plus(n,V,T,R_T,C_sigma)+gamma2minus(n,V,T,R_T,C_sigma))*sigma[-1])/(gamma1minus(n+1,V,T,R_T,C_sigma)+gamma2plus(n+1,V,T,R_T,C_sigma))
            sigma2 = -((gamma1plus(n-2,V,T,R_T,C_sigma)+gamma2minus(n-2,V,T,R_T,C_sigma))*sigma[n-2]-(gamma1plus(n-1,V,T,R_T,C_sigma)+gamma1minus(n-1,V,T,R_T,C_sigma)+gamma2plus(n-1,V,T,R_T,C_sigma)+gamma2minus(n-1,V,T,R_T,C_sigma))*sigma[n-1])/(gamma1minus(n,V,T,R_T,C_sigma)+gamma2plus(n,V,T,R_T,C_sigma))           
            sigma=append(sigma,sigma2)
            current_curr_plus = sigma[n]*(gamma1plus(n,V,T,R_T,C_sigma)-gamma1minus(n,V,T,R_T,C_sigma))
            current_curr_minus = sigma[n]*(gamma1plus(-n,V,T,R_T,C_sigma)-gamma1minus(-n,V,T,R_T,C_sigma))
            if abs((abs(current_curr_plus+current_curr_minus))/ref_curr)<1e-10: # criteria for stopping calculation
                stop_calculation=True
                #print "MAIN: stopped calculation at n = %g"%n
            #print sigma
            #print "ref_curr %g"%ref_curr
            #print "current_curr_plus %g"%current_curr_plus
            #print "current_curr_minus %g"%current_curr_minus
            #print abs((abs(current_curr_plus+current_curr_minus))/ref_curr)
            #print n
            if False: # for debug
                print "------N = %i------"%n
                print gamma1plus(n,V,T,R_T,C_sigma)
                print gamma1minus(n,V,T,R_T,C_sigma)
                print gamma2plus(n,V,T,R_T,C_sigma)
                print gamma2minus(n,V,T,R_T,C_sigma)

            current += current_curr_plus + current_curr_minus
    sigma_sum = sigma[0]+2*sum(sigma[1:])
    #print sigma_sum
    ##sigma=sigma/sigma_sum
    #print sigma
    current = current/sigma_sum*e
    return current


def calc_current_full(sigma,N,V,R_T,C_sigma,T_p,island_volume,const_P=0.1e-15):
    """
    calculates current with given parameters
    sigma: electron phonon coupling constant
    N: junctions
    V: voltage
    R_T: tunnel resistance
    C_sigma: total island capacitance
    T_p: phonon (bath) temperature
    island_volume: island volume for e-p coupling
    """
#<<<<<<< .mine
#    T_e = (4*(V/N)**2/(R_T*sigma*island_volume)+T_p**5.0)**(1.0/5.0)
    #print "V %g"%V
    #print "T_e %g"%T_e
    #print "T_p %g"%T_p
#=======
    #T_e = (4*(V/N)**2/(R_T*sigma*island_volume)+T_p**5.0)**(1.0/5.0)
    P_per_junc =  pow(V/N,2.0)/R_T+const_P
    #T_e = (1*(V/N)**2/(R_T*sigma*island_volume)+T_p**5.0)**(1.0/5.0)
    T_e = (P_per_junc/(sigma*island_volume)+T_p**5.0)**(1.0/5.0)
    if False:
        print "V %g"%V
        print "T_e %g"%T_e
        print "T_p %g"%T_p
#>>>>>>> .r1091
    #print "T_e = %g"%T_e
    current = calc_current(V/(N/2.0),T_e,R_T,C_sigma)
    return current

def calc_temp(sigma,N,V,R_T,C_sigma,T_p,island_volume,const_P=0.1e-15):
    """
    calculates electron temperature
    sigma: electron phonon coupling constant
    N: junctions
    V: voltage
    R_T: tunnel resistance
    C_sigma: total island capacitance
    T_p: phonon (bath) temperature
    island_volume: island volume for e-p coupling
    """
    P_per_junc =  pow(V/N,2.0)/R_T+const_P
    T_e = (P_per_junc/(sigma*island_volume)+T_p**5.0)**(1.0/5.0)
    return T_e    
    
    
def calc_conductance_curve(V_list,T,R_T,C_sigma):
    #test_voltages = arange(-v_max,v_max,v_step)
    test_currents = []
    for V in V_list:
        test_currents.append(calc_current(V,T,R_T,C_sigma))
        #print "V: %g, current %g"%(V,test_currents[-1])

    ## calc conductances manually
    #test_conductances = []
    #for idx,V in enumerate (test_currents[1:-2]):
    #    if idx==0:
    #        print idx
    #    test_conductances.append((test_currents[idx+2]-test_currents[idx])/(2.0*v_step))
    #
    #test_voltages_G = test_voltages[1:-2]

    #
    # SPLINE
    #
    spline = UnivariateSpline(V_list,test_currents,s=0)
    #print "test_conductances"
    #indices = [x for x, y in enumerate(col1) if (y >0.7 or y<-0.7)]
    test_conductances = []
    for v_iter in V_list:
        test_conductances.append(spline.derivatives(v_iter)[1])
    return test_conductances
#print test_conductances
'''
extern double calc_G(double sigma,double N,double V,double R_T,double C_sigma,double T_p,double island_volume, double const_P,double eps) {
    double curr[4];
    double eps2 = eps/2.0; /* then differential is calculated in given range */
    if (eps2==0) { /* 0.1 nanovolts if zero */
        eps2=1e-10; 
        }
    curr[0] = calc_current_full(sigma,N,V+2*eps2,R_T,C_sigma,T_p,island_volume,const_P);
    curr[1] = calc_current_full(sigma,N,V+1*eps2,R_T,C_sigma,T_p,island_volume,const_P);
    curr[2] = calc_current_full(sigma,N,V-1*eps2,R_T,C_sigma,T_p,island_volume,const_P);
    curr[3] = calc_current_full(sigma,N,V-2*eps2,R_T,C_sigma,T_p,island_volume,const_P);
    double R = (-curr[0]+8*curr[1]-8*curr[2]+curr[3])/(12.0*eps2);
    return R;
 '''
def calc_G(sigma,N,V,R_T,C_sigma,T_p,island_volume, const_P, eps=1e-9):
    curr=[0.0,0.0,0.0,0.0]
    eps2 = eps/2.0 # then differential is calculated in given range */
    if (eps2==0): # 0.1 nanovolts if zero */
            eps2=1e-10 
    curr[0] = calc_current_full(sigma,N,V+2.0*eps2,R_T,C_sigma,T_p,island_volume,const_P)
    curr[1] = calc_current_full(sigma,N,V+1.0*eps2,R_T,C_sigma,T_p,island_volume,const_P)
    curr[2] = calc_current_full(sigma,N,V-1.0*eps2,R_T,C_sigma,T_p,island_volume,const_P)
    curr[3] = calc_current_full(sigma,N,V-2.0*eps2,R_T,C_sigma,T_p,island_volume,const_P)
    R = (-curr[0]+8.0*curr[1]-8.0*curr[2]+curr[3])/(12.0*eps2)
    return R

def calc_conductance_curve_full(sigma,N,V_list,R_T,C_sigma,T_p,island_volume):
    """
    Calculates full conductance curve taking into account self-heating
    """
    test_currents = []
    for V in V_list:
        test_currents.append(calc_current_full(sigma,N,V,R_T,C_sigma,T_p,island_volume))
    # SPLINE
    spline = UnivariateSpline(V_list,test_currents,s=0)
    test_conductances = []
    for v_iter in V_list:
        test_conductances.append(spline.derivatives(v_iter)[1])
    return test_conductances

def calc_conductance_curve_basic(sigma,N,V_list,R_T,C_sigma,T_p,island_volume, const_P=1e-18, eps=1e-9):
    """
    Calculates conductance curve without any spline or other tricks
    """
    conductances=[]
    for point in V_list:
        conductances.append(calc_G(sigma,N,point,R_T,C_sigma,T_p,island_volume, const_P, eps))
    return conductances
    

def calc_conductance_curve_few_points(sigma,N,V_list,R_T,C_sigma,T_p,island_volume):
    """
    When only few points extra points are added in order to make spline approximation correct
    """
    extra_points = []
    for point in V_list:
        extra_points.append(0.99001*point)
        extra_points.append(0.95001*point)
        extra_points.append(1.01001*point)
        extra_points.append(1.05001*point)
    
    all_V_points = sort(V_list + extra_points).tolist()
    test_conductances = calc_conductance_curve_full(sigma,N,all_V_points,R_T,C_sigma,T_p,island_volume)
    # get the results in original order
    sparse_conductances=[]
    for point in V_list:
        sparse_conductances.append(test_conductances[all_V_points.index(point)])
    return sparse_conductances
        


def calc_conductance_curve_full_spline(sigma,N,meas_V,R_T,C_sigma,T_p,island_volume):
    v_min = min(meas_V)
    v_max = max(meas_V)
    v_step = (v_max-v_min)/201.0
    test_voltages = arange(v_min,v_max+v_step,v_step)
    calc_G_few = calc_conductance_curve_full(sigma,N,test_voltages,R_T,C_sigma,T_p,island_volume)
    interpolator = interpolate.interp1d(test_voltages, calc_G_few, kind='cubic')
    calc_G = []
    for V in meas_V:
        calc_G.append(interpolator(V))
    return calc_G