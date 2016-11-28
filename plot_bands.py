# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 09:11:27 2015

@author: nacho
"""

import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

def plot_bands(t,data,axis,color=(0.,0.,1.),label='somemodel'):
    p25 = stats.scoreatpercentile(data,25,axis=axis)
    p75 = stats.scoreatpercentile(data,75,axis=axis)
    ctrans =(color[0],color[1],color[2],0.1)
#    p10 = stats.scoreatpercentile(data,10,axis=axis)
#    p90 = stats.scoreatpercentile(data,90,axis=axis)
#    plt.fill_between(t,p10,p90,color=ctrans)
    ctrans =(color[0],color[1],color[2],0.2)
    plt.fill_between(t,p25,p75,color=ctrans)
    me = stats.scoreatpercentile(data,50,axis=axis)
    return plt.plot(t,me,lw=2.0,color=color,label=label)