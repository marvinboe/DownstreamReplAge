#!/usr/bin/env python3
#######################################################################
#filename: 'plot_cml_effect.py'
#Plots expected change of replicative age distribution with CML.
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

# import sympy as sp
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
ymaxlim=0.17
ymaxdata=0.09


def conv(input_v,js=None,a=0.5):
    if not js:
        js=np.arange(len(input_v))
    retvec=[]
    for j in js:
        retsum=np.nansum([input_v[k]*a**(j-k) for k in np.arange(j+1)])
        retvec.append(retsum)
    return np.array(retvec)

def gen_comp_cycle(influx,compartments=1,alpha=0.5,r=1,cmap="viridis"):
    a=alpha
    # color_cycle=plothelpers.create_colorcyle(compartments+1,cmapname=cmap)
    if not isinstance(influx, np.ndarray):
        raise ValueError("wrong influx input")

    statslist=[]
    age_dist_buf=influx[:]

    for i in np.arange(compartments):
        output=(1-alpha)*np.roll(age_dist_buf,1) #cells always become older
        output[0]=0. #cells always become older
        age_dist=conv(output,a=a)/r
        age_dist=age_dist[:]/np.sum(age_dist)
        age_dist_buf=np.array(age_dist[:])
        ###plot progression
        # color=next(color_cycle)
        # ax.plot(np.arange(len(age_dist)),age_dist,color=color)

    return age_dist_buf

def mean(muin,compartments,alpha):
    return muin+compartments*1/(1-alpha)

def variance(varin,compartments,alpha):
    return varin+compartments*alpha/(1-alpha)**2

def alpha_for_N(comps,Influx=10,Outflux=1000):
    c=(Outflux/Influx)**(1/float(comps))
    alpha_sol= (c-2)/(c-1)
    return alpha_sol


def plot_full(ax,influx,comps=1,alpha=None,p=None,d=None,color=None,with_comps_label=None,linewidth=None):
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
    elif with_comps_label=="both":
        v1='$C={0:d}$'.format(int(comps))
        v2='$\\alpha={0:.2}$'.format(alpha)
        label='\\begin{tabular}{p{0.95cm} p{1.35cm}}'+v1+"&\\hspace{-0.4cm}"+v2+"\\end{tabular}"
    else:
        label=with_comps_label

    # ax.step(np.arange(len(dist)),dist,color=color,alpha=1,
    #     path_effects=[pe.Stroke(linewidth=linewidth+0.80, foreground='black',alpha=0.8), pe.Normal()],
    #         label=label,lw=linewidth,where='mid')

    x=np.arange(len(dist))
    ax.bar(x,dist,width=0.95,color=color,
            label=label,alpha=0.98,
            lw=0.1)#,linestyle='None',marker='_',ms=6)
    ax.step(x,dist,alpha=.8,color="black",lw=0.2,where='mid')# $b="+str(age_param)+"$")
    if max(dist)>ymaxlim:
        plot_datacap(ax,x[np.argmax(dist)],ymaxlim-0.01)

    meanvalue=mean(10,comps,alpha)
    meanlinemax=max(dist)+0.010
    ax.plot([meanvalue,meanvalue],[0,meanlinemax],
            color="black",lw=1.0,linestyle='--')
    ax.annotate('\\footnotesize $\\mu={:.3}$'.format(meanvalue), xy=(meanvalue,meanlinemax),
            xytext=(meanvalue,meanlinemax+0.005),ha='center')

    stdvalue=math.sqrt(variance(10,comps,alpha))
    stdliney=dist[int(meanvalue-stdvalue)]
    offset=1.
    if alpha <0.5:
        textpos=(meanvalue-stdvalue-offset,stdliney)
        ha="right"
    else:
        textpos=(meanvalue+stdvalue+offset,stdliney)
        ha="left"

    ax.plot([meanvalue-stdvalue,meanvalue+stdvalue],[stdliney,stdliney],
            color="black",lw=1.0,linestyle='--')
    ax.annotate('\\footnotesize $\\sigma={:.3}$'.format(stdvalue), xy=(meanvalue,stdliney),
            xytext=textpos,ha=ha)


def plot_several_comps_dist(ax,influx,alphas=None,compartments=None,labels=None,with_axes='',with_legend=False,distlabel=None):

    if compartments is None:
        compartments=30

    color_cycle=plothelpers.create_colorcyle(len(alphas),cmapname="Set1_r")
    
    linewidth=1.1
    # ax.step(x,influx,alpha=.9,label="influx",color="black",lw=linewidth,where='mid')# $b="+str(age_param)+"$")
    if max(influx)>ymaxlim:
        plot_datacap(ax,x[np.argmax(influx)],ymaxlim-0.01)
    
    if labels is not None:
        for alpha,label in zip(alphas,labels):
            col=next(color_cycle)
            plot_full(ax,influx,compartments,alpha=alpha,color=col,with_comps_label=label)
    
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
xmax=140
x=np.arange(xmax) #discrete ages
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
poiss=poisson(x,10)
# poiss=poiss/np.sum(poiss)

plothelpers.latexify(columns=1)
fig,ax=plt.subplots(1,1)#,sharey=True)

labels=["healthy $\\alpha=0.3$","CML\\hspace{0.31cm} $\\alpha=0.58$"]
plot_several_comps_dist(ax,poiss,compartments=29,alphas=[0.3,0.58],with_axes='xy',with_legend=True,labels=labels)

ax.set_ylim(ymin=-0.005,ymax=0.11)#ymax=ymaxlim,
ax.set_xlim(xmax=130,xmin=22)

for k, spine in ax.spines.items():  #ax.spines is a dictionary
    spine.set_zorder(10)

fig.tight_layout(pad=0.1)
fig.subplots_adjust(left=0.21)
plt.show()



