#!/usr/bin/env python3
#######################################################################
#filename: 'meanvariance_anaplot.py'
#Plots progression of mean and variance through multiple compartments.
#
#Copyright 2018 Marvin A. Böttcher
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
########################################################################

import numpy as np
import itertools
import matplotlib.pyplot as plt
import plothelpers
import math

def anamean(inmean,alpha,comps):
    return inmean+comps/(1.-alpha)

def anavar(invar,alpha,comps):
    return invar+comps*alpha/(1.-alpha)**2

def alpha_for_N(comps,Influx=10,Outflux=None):
    if Outflux is None:
        Outflux=Influx*2**comps
    c=(Outflux/Influx)**(1/(comps))
    alpha_sol= (c-2)/(c-1)
    return alpha_sol

def plot(ax,inmean=None,invar=None,maxcomps=10,alpha=None,p=None,d=None,color=None,with_comps_label=None):
    if alpha is None:
        if p is None or d is None:
            print("error: alpha, p or d not specified")
            exit(0)
        else:
            alpha=1+p-d


    if with_comps_label is None:
         label="$\\alpha={0:.2}$".format(alpha)
    elif with_comps_label == True:
        label="$C={0:d}$".format(maxcomps)
    else :
        v1='$C={0:d}$'.format(int(maxcomps))
        v2='$\\alpha={0:.2}$'.format(alpha)
        label='\\begin{tabular}{p{1.cm} p{1.35cm}}'+v1+"&\\hspace{-0.4cm}"+v2+"\\end{tabular}"

    ls=''
    marker='o'
    ms=6.
    mew=0.2

    ##generate plotdata
    comps=np.arange(maxcomps)
    if inmean is not None:
        data=anamean(inmean,alpha,comps)
    elif invar is not None:
        data=anavar(invar,alpha,comps)
    x=np.linspace(0,1,len(comps))

    ##plot
    ax.plot(x,data,linewidth=0.1,
             linestyle=ls,mew=mew,marker=marker,mec="black",
             ms=ms,color=color,label=label)#,where='mid')


def plot_meanvariance(ax,ax2,inmean,invar,compartments):
    color_cycle=plothelpers.create_colorcyle(len(compartments),cmapname="summer_r")

    largestcomps=max(compartments)
    for comps in compartments[::-1]:
        alpha=alpha_for_N(comps,Influx=10,Outflux=10*2**largestcomps)
        col=next(color_cycle)
    
        plot(ax,inmean=inmean,maxcomps=comps,alpha=alpha,color=col,with_comps_label="both")
        plot(ax2,invar=invar,maxcomps=comps,alpha=alpha,color=col,with_comps_label="both")

#start plotting
# plt.xkcd()
plothelpers.latexify(columns=2,fig_height=2)
fig,axes=plt.subplots(1,2)
ax=axes[0]
ax2=axes[1]

#influx distribution
inmean=10
invar=10


comps=[10,15,20]
plot_meanvariance(ax,ax2,inmean,invar,compartments=comps)


ax.set_ylabel("mean $\mu$")
ax.set_xlabel("relative compartment $c/C$")
ax2.set_ylabel("variance $\sigma^2$")
ax2.set_xlabel("relative compartment $c/C$")

letters=itertools.cycle('(a),(b),(c),(d)'.split(','))
for b in [ax,ax2]:
    letter=next(letters)
    b.text(-.17,.94,letter,
    horizontalalignment='left',family="sans-serif",weight="heavy",
    transform=b.transAxes)

# ax2.set_yscale('log')
ax.legend()
plt.tight_layout(pad=0.4)

plt.show()

