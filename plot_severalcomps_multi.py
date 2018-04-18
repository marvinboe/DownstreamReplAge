#!/usr/bin/env python3
#######################################################################
#filename: 'plot_several_comps_multi.py'
#Plots distribution of replicative age after many downstream compartment for
#different influx distributions.
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

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patheffects as pe
import math
import scipy.stats
import scipy.special
import plothelpers
import itertools
import numpy as np
import warnings
warnings.filterwarnings("error")

ymaxlim=0.17 #limit on yaxis (frequency) before data is cutoff.

def conv(input_v,js=None,a=0.5):
    """convolution sum of influx with exponential"""
    if not js:
        js=np.arange(len(input_v))
    retvec=[]
    for j in js:
        retsum=np.nansum([input_v[k]*a**(j-k) for k in np.arange(j+1)])
        retvec.append(retsum)
    return np.array(retvec)

def gen_comp_cycle(influx,compartments=1,alpha=0.5,r=1):
    """create distribution for progression through several compartments"""
    a=alpha
    color_cycle=plothelpers.create_colorcyle(compartments+1)
    if not isinstance(influx, np.ndarray):
        raise ValueError("wrong influx input")

    statslist=[]
    age_dist_buf=influx[:]

    for i in np.arange(compartments):
        output=(1-alpha)*np.roll(age_dist_buf,1) #cells always become older by one
        output[0]=0. #cells always become older by one
        age_dist=conv(output,a=a)/r #convolution sum, calculates age distribution
        age_dist=age_dist[:] #/np.sum(age_dist) no normalization here
        age_dist_buf=np.array(age_dist[:]) #copy old distribution
        ###plot progression
        # color=next(color_cycle)
        # ax.plot(np.arange(len(age_dist)),age_dist,color=color)

    return age_dist_buf

def alpha_for_N(comps,Influx=10,Outflux=1000):
    """calculates alpha based on required outflux for given influx"""
    c=(Outflux/Influx)**(1/float(comps))
    alpha_sol= (c-2)/(c-1)
    return alpha_sol


def plot_full(ax,influx,comps=1,alpha=None,p=None,d=None,color=None,with_comps_label=None,linewidth=None):
    """plot distribution for given parameters"""

    if not isinstance(influx, np.ndarray):
        raise ValueError("wrong influx input")
    if linewidth is None:
        linewidth=1.5

    if alpha is None:
        if p is None or d is None:
            print("error alpha, p or d not specified")
            exit(0)
        else:
            alpha=1+p-d

    dist=gen_comp_cycle(influx,compartments=comps,alpha=alpha)

    # color=next(color_cycle)
    if with_comps_label is None:
         label="$\\alpha={0:.2}$".format(alpha)
    elif with_comps_label == True:
        label="$C={0:d}$".format(comps)
    else :
        v1='$C={0:d}$'.format(int(comps))
        v2='$\\alpha={0:.2}$'.format(alpha)
        label='\\begin{tabular}{p{0.95cm} p{1.35cm}}'+v1+"&\\hspace{-0.4cm}"+v2+"\\end{tabular}"

    x=np.arange(len(dist))
    ax.bar(x,dist,width=0.95,color=color,
            label=label,alpha=0.98,
            lw=0.1)#,linestyle='None',marker='_',ms=6)
    ax.step(x,dist,alpha=.8,color="black",lw=0.2,where='mid')# outline
    if max(dist)>ymaxlim:
        plothelpers.plot_datacap(ax,x[np.argmax(dist)],ymaxlim-0.01,xint=3,yint=0.008)

def plot_several_comps_dist(ax,influx,compartments=None,with_axes='',with_legend=False,distlabel=None):

    if compartments is None:
        compartments=[int(i) for i in np.linspace(1,maxcomps,comps_to_plot)]
    maxcomps=max(compartments)
    color_cycle=plothelpers.create_colorcyle(len(compartments),cmapname="summer_r")
    
    ### plot influx distribution as step function
    linewidth=1.1
    ax.step(x,influx,alpha=.9,label="influx",color="black",lw=linewidth,where='mid')# $b="+str(age_param)+"$")
    if max(influx)>ymaxlim:
        plothelpers.plot_datacap(ax,x[np.argmax(influx)],ymaxlim-0.01,xint=3,yint=0.008)
    
    ### cycle trogh number of compartments and plot distribution
    for comps in compartments[::-1]:
        alpha=alpha_for_N(comps,Influx=10,Outflux=10*2**maxcomps)
        col=next(color_cycle)
    
        plot_full(ax,influx,comps,alpha=alpha,color=col,with_comps_label="both")
    
    if 'x' in with_axes:
        ax.set_xlabel("replicative age")
    else:
        ax.tick_params(labelbottom='off') 
        # ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
    if 'y' in with_axes:
        ax.set_ylabel("frequency")
    else:
        ax.tick_params(labelleft='off') 
        # ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())

    if distlabel is not None:
        ax.text(0.98,0.97,distlabel, horizontalalignment='right',
        verticalalignment='top', transform=ax.transAxes,fontsize=9)
    if with_legend is True:
        leg=ax.legend(loc=0)
        if distlabel is not None:
            # Get the bounding box of the original legend
            bb = leg.get_bbox_to_anchor().inverse_transformed(ax.transAxes)
            bb.y1 += -0.15
            bb.x1 += -0.01
            leg.set_bbox_to_anchor(bb)#, transform = ax.transAxes)
            bb = leg.get_bbox_to_anchor().inverse_transformed(ax.transAxes)




### distributions
xmax=80
x=np.arange(xmax)
age=15
r0=1.
N0=1.

def poisson(x,b):
    if x > 100:
        f1=(b*math.e/x)**x
        return f1*math.exp(-b)/math.sqrt(2*math.pi*x)
    else:
        return b**x/math.gamma(x+1) *math.exp(-b)
poisson=np.vectorize(poisson)
poiss=poisson(x,r0*age/N0)

deltadist=np.zeros(len(x))
deltadist[5]=1

geomdist=scipy.stats.geom(.15,loc=-1)
geompmf=geomdist.pmf(x)
geompmf=geompmf

def werner_sym(xval,p,tstar,N0=N0):
    # r=1
    c=1000
    # tstar=r*p*t/N0+1
    if xval<c:
        retvalue=N0/math.gamma(xval+1) *((1+p)/p)**(xval)*(tstar)**(1./-p)*math.log(tstar)**xval
    elif xval==c:
        retvalue=N0*((1+p)**(xval-1)*(1-(scipy.special.gammaincc(xval,p**(-1)*math.log(tstar)))/math.gamma(xval)))
    else:
        retvalue=0
    return retvalue
werner_sym=np.vectorize(werner_sym)
symprobstemcell=0.1
wernerdist=werner_sym(x,symprobstemcell,r0*age*symprobstemcell/N0+1)
wernerdist=wernerdist/np.sum(wernerdist) #normalize



plothelpers.latexify(columns=2)
fig,ax=plt.subplots(2,2)#,sharey=True)

### plot for four different influx distributions on different axes
comps=[10,15,20]
plot_several_comps_dist(ax[0,0],deltadist,compartments=comps,with_axes='y',distlabel="single age")
plot_several_comps_dist(ax[1,0],geompmf,compartments=comps,with_axes='xy',distlabel="geometric dist.")
plot_several_comps_dist(ax[0,1],poiss,compartments=comps,distlabel="purely asymmetric $p_s=0$")
plot_several_comps_dist(ax[1,1],wernerdist,compartments=comps,with_axes='x',with_legend=True,distlabel="symmetric $p_s={:}$".format(symprobstemcell))


### plot styling
letters=itertools.cycle('(a),(b),(c),(d)'.split(','))
for a in ax:
    for b in a:
        letter=next(letters)
        b.text(.02,.9,letter,
        horizontalalignment='left',family="sans-serif",weight="heavy",
        transform=b.transAxes)

        b.set_ylim(ymax=ymaxlim,ymin=-0.005)

for k, spine in ax[0,0].spines.items():  #ax.spines is a dictionary
    spine.set_zorder(10)

fig.tight_layout(pad=0.4)
fig.subplots_adjust(hspace=0.08,wspace=0.05,right=0.99,left=0.09)
plt.show()



