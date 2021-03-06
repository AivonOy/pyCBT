# -*- coding: utf-8 -*-
"""
Created on Tue Aug 06 14:03:31 2013

@author: Leif
"""

from pyx import *
from numpy import *

class CBT_plot_data_pyx:
    def __init__(self,fitters):
        self.fitters=fitters    
     
    def plot_data_given_params(self,T,R_T,C_sigma,filename="result_pfit.pdf"):
        """
        Plots data with given params.
        T: array of temperatures
        R_T: scalar tunnelling resistance
        C_sigma: scalar island capacitance
        """
        datas=[]
        v_min=0.0
        v_max=0.0
        for i,fitter in enumerate(self.fitters):
            datas.append({})
            datas[i]['V_data'] = (array(fitter.original_V())*1e6*fitter.junctions_in_series)
            datas[i]['G_orig'] = (array(fitter.original_G())*1e6)
            datas[i]['G_tuned'] = (array(fitter.G_curve(R_T,
                                 C_sigma,T[i]))*1e6)
            datas[i]['T_tuned'] = T[i]
            v_min_0 =datas[i]['V_data'].min()
            v_max_0 =datas[i]['V_data'].max()
            if (v_min_0<v_min): v_min=v_min_0
            if (v_max_0>v_max): v_max=v_max_0    
            
        #plot datas
        g = graph.graphxy(width=8,y=graph.axis.linear(title="$G$ ($\mu$S)"),
                          x=graph.axis.linear(title="$V$($\mu$V)",
                                                min=v_min,max=v_max))
        for data in datas:
            g.plot(graph.data.values(x=data['V_data'],y=data['G_tuned'],
                                     title="$T_{fit}$ = %g"%data['T_tuned']),
                   [graph.style.line(lineattrs=[style.linewidth.thick, style.linestyle.solid, 
                                                color.rgb.red])])
            g.plot(graph.data.values(x=data['V_data'],y=data['G_orig']),
                   [graph.style.symbol(symbol=graph.style.symbol.circle,size=0.05*unit.v_cm)])
            #g.plot(graph.data.values(x=data['V_data'],y=data['G_orig']),
            #       [graph.style.line(lineattrs=[style.linewidth.thin, style.linestyle.dotted, 
            #                                    color.rgb.blue])])
        g.dolayout()

        # or provide one list containing the whole points
        g.writePDFfile(filename)
    
    def plot_data(self,filename='result.pdf'):
        """
        plots fit result
        """
        datas=[]
        v_min=0.0
        v_max=0.0
        for i,fitter in enumerate(self.fitters):
            datas.append({})
            datas[i]['V_data'] = (array(fitter.original_V())*1e6*fitter.junctions_in_series)
            datas[i]['G_orig'] = (array(fitter.original_G())*1e6)
            datas[i]['G_fit'] = (array(fitter.nonlinfit_G())*1e6)
            datas[i]['T_fit'] = fitter.T_fit
            v_min_0 =datas[i]['V_data'].min()
            v_max_0 =datas[i]['V_data'].max()
            if (v_min_0<v_min): v_min=v_min_0
            if (v_max_0>v_max): v_max=v_max_0    

        #plot datas        
        g = graph.graphxy(width=8,y=graph.axis.linear(title="$G$ ($\mu$S)"),
                          x=graph.axis.linear(title="$V$($\mu$V)",
                                                min=v_min,max=v_max),
                                                key=graph.key.key(pos="tr", dist=0.1))
        for data in datas:
            g.plot(graph.data.values(x=data['V_data'],y=data['G_fit'],
                               title="$T_{fit}$ = %3.1f"%data['T_fit']),
                   [graph.style.line(lineattrs=[style.linewidth.thick, style.linestyle.solid, 
                                                color.rgb.red])])
            g.plot(graph.data.values(x=data['V_data'],y=data['G_orig'],
                                     title="meas data"),
                   [graph.style.symbol(symbol=graph.style.symbol.circle,size=0.05*unit.v_cm)])
            #g.plot(graph.data.values(x=data['V_data'],y=data['G_orig']),
            #       [graph.style.line(lineattrs=[style.linewidth.thin, style.linestyle.dotted, 
            #                                    color.rgb.blue])])
        g.dolayout()
        g.writePDFfile(filename)

    def plot_multi_data(self,filename='result.pdf'):
        """
        plots data
        """
        datas=[]
        for i,fitter in enumerate(self.fitters):
            datas.append({})
            datas[i]['V_data'] = (array(fitter.original_V())*1e6*fitter.junctions_in_series)
            datas[i]['G_orig'] = (array(fitter.original_G())*1e6)
            datas[i]['G_fit'] = (array(fitter.G_curve(fitter.R_T_multi,
                                 fitter.C_sigma_multi,fitter.T_multi))*1e6)
            datas[i]['T_fit'] = fitter.T_fit
            v_min =datas[i]['V_data'].min()
            v_max =datas[i]['V_data'].max()
        
        #print datas
        
        g = graph.graphxy(width=8,y=graph.axis.linear(title=r"$G$ ($\mu$S)"),
                          x=graph.axis.linear(title=r"$V$($\mu$V)",
                                                min=v_min,max=v_max))
        for data in datas:
            g.plot(graph.data.values(x=data['V_data'],y=data['G_fit'],
                                     title="$T_{fit}$ = %g"%data['T_fit']),
                   [graph.style.line(lineattrs=[style.linewidth.thick, style.linestyle.solid, 
                                                color.rgb.red])])
            g.plot(graph.data.values(x=data['V_data'],y=data['G_orig']),
                   [graph.style.symbol(symbol=graph.style.symbol.circle,size=0.05*unit.v_cm)])
            #g.plot(graph.data.values(x=data['V_data'],y=data['G_orig']),
            #       [graph.style.line(lineattrs=[style.linewidth.thin, style.linestyle.dotted, 
            #                                    color.rgb.blue])])
        g.dolayout()

        # or provide one list containing the whole points
        g.writePDFfile(filename)
        
    def plot_data_1(self,filename='result1.pdf'):
        
        datas=[]
        for i,fitter in enumerate(self.fitters):
            datas.append({})
            datas[i]['V_data'] = (array(fitter.original_V())*1e6*fitter.junctions_in_series)
            datas[i]['G_orig'] = (array(fitter.original_G())*1e6)
            datas[i]['G_fit'] = (array(fitter.nonlinfit_G())*1e6)
            datas[i]['T_fit'] = fitter.T_fit
            v_min =datas[i]['V_data'].min()
            v_max =datas[i]['V_data'].max()
        
        g = graph.graphxy(width=8, key=graph.key.key())

        As = [0.3, 0.6, 0.9]

        d = [graph.data.join([graph.data.values(x_a=datas[i]['V_data'],y_a=datas[i]['G_fit'],context=dict(A=A)),
                      graph.data.values(x_b=datas[i]['V_data'],y_b=datas[i]['G_orig'])],
                     title=r"$A=%g$" % A)
                     for i, A in enumerate(As)]

        attrs = [color.gradient.RedBlue]

        g.plot(d,
               [graph.style.pos(usenames=dict(x="x_a", y="y_a")),
                graph.style.line(attrs),
                graph.style.pos(usenames=dict(x="x_b", y="y_b")),
                graph.style.symbol(graph.style.symbol.changesquare, symbolattrs=attrs, size=0.1)])
    