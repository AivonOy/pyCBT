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
        
    
    def plot_data(self,filename='result.pdf'):
        #print self.xdata
        # either provide lists of the individual coordinates
        datas=[]
        for i,fitter in enumerate(self.fitters):
            datas.append({})
            datas[i]['V_data'] = (array(fitter.original_V())*1e6*fitter.junctions_in_series)
            datas[i]['G_orig'] = (array(fitter.original_G())*1e6)
            datas[i]['G_fit'] = (array(fitter.nonlinfit_G())*1e6)
            datas[i]['T_fit'] = fitter.T_fit
            v_min =datas[i]['V_data'].min()
            v_max =datas[i]['V_data'].max()
        
        #print datas
        
        g = graph.graphxy(width=8,y=graph.axis.linear(title="$G$ ($\mu$S)"),
                          x=graph.axis.linear(title="$V$($\mu$V)",
                                                min=v_min,max=v_max))
        for data in datas:
            g.plot(graph.data.values(x=data['V_data'],y=data['G_fit'],
                                     title="$T_{fit}$ = %g"%data['T_fit']),
                   [graph.style.line(lineattrs=[style.linewidth.thick, style.linestyle.solid, 
                                                color.rgb.red])])
            g.plot(graph.data.values(x=data['V_data'],y=data['G_orig']),
                   [graph.style.symbol()])
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
    