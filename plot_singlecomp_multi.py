#!/usr/bin/env python3
#######################################################################
#filename: 'plot_singlecomp_multi.py'
#Plots replicative age distribution in progenitor compartment for
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

alpha_range=[0.8,0.5,0.2]
cont_xmax=30
xvals=np.arange(30) #discrete ages

def conv(input_v,js=None,alpha=0.5):
    """convolution sum of influx with exponential"""
    if not js:
        js=np.arange(len(input_v))
    retvec=[]
    for j in js:
        retsum=alpha**j*np.nansum([input_v[k]*alpha**(-k) for k in np.arange(j+1)])
        retvec.append(retsum)
    return np.array(retvec)

def plot_numdist(ax,indist,cmapname=None,with_axes='',with_legend=False,distlabel=None):

    scaling=(cont_xmax/float(len(indist)))
    xvals=np.arange(len(indist))#*scaling
    colors=plothelpers.create_colorcyle(len(alpha_range),cmapname=cmapname)

    #plot influx
    linewidth=1.5
    ax.step(xvals,indist,color="black",label="influx",where='mid',lw=linewidth)#,linestyle='None',marker='_',ms=6)
    for alpha in alpha_range:
        color=next(colors)
        label="$\\alpha={}$".format(alpha)

        dist=conv(indist,alpha=alpha)
        # ax.step(xvals,dist,color=color,label=label,
        #     path_effects=[pe.Stroke(linewidth=linewidth+0.80, foreground='black',alpha=0.8), pe.Normal()],
        #         where='mid',lw=linewidth)#,linestyle='None',marker='_',ms=6)

        ax.bar(xvals,dist,color=color,label=label,alpha=0.95,
                lw=linewidth)#,linestyle='None',marker='_',ms=6)

    if 'x' in with_axes:
        ax.set_xlabel("replicative age")
    else:
        ax.tick_params(labelbottom='off') 
        # ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
    if 'y' in with_axes:
        ax.set_ylabel("frequency")
        # ax.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=4))
    else:
        ax.tick_params(labelleft='off') 
        # ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
    
    if distlabel is not None:
        ax.text(0.98,0.97,distlabel, horizontalalignment='right',
        verticalalignment='top', transform=ax.transAxes,fontsize=9)
    if with_legend:
        leg=ax.legend(loc=0)
        if distlabel is not None:
            # Get the bounding box of the original legend
            bb = leg.get_bbox_to_anchor().inverse_transformed(ax.transAxes)
            bb.y1 += -0.20
            bb.x1 += -0.10
            leg.set_bbox_to_anchor(bb)#, transform = ax.transAxes)
            bb = leg.get_bbox_to_anchor().inverse_transformed(ax.transAxes)



plothelpers.latexify(columns=2)
fig,ax=plt.subplots(2,2)#,sharey=True)


#### distributions
lambwerner=10.
age=10
r0=1.
N0=1.

def poisson(x,b):
    if x > 100:
        f1=(b*math.e/x)**x
        return f1*math.exp(-b)/math.sqrt(2*math.pi*x)
    else:
        return b**x/math.gamma(x+1) *math.exp(-b)
poisson=np.vectorize(poisson)
poiss=poisson(xvals,age*r0/N0)
poiss=poiss

deltadist=np.zeros(len(xvals))
deltadist[5]=1

geomdist=scipy.stats.geom(.15,loc=-1)
geompmf=geomdist.pmf(xvals)
geompmf=geompmf/np.sum(geompmf)

def werner_sym(xval,p,tstar):
    N0=1
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
symprobstemcell=0.10
wernerdist=werner_sym(xvals,symprobstemcell,1+age*symprobstemcell*r0/N0)
wernerdist=wernerdist/np.sum(wernerdist) #normalize
#### end distributions



### plotting
plot_numdist(ax[0,0],deltadist,with_axes='y',distlabel="single age")
plot_numdist(ax[1,0],geompmf,with_axes='xy',distlabel="geometric dist.")
plot_numdist(ax[0,1],poiss,distlabel="purely asymmetric $p_s=0$")
plot_numdist(ax[1,1],wernerdist,with_axes='x',with_legend=True,distlabel="symmetric $p_s={:}$".format(symprobstemcell))

letters=itertools.cycle('(a),(b),(c),(d)'.split(','))
for a in ax:
    for b in a:
        b.set_ylim(ymin=-0.05,ymax=1.05)
        letter=next(letters)
        b.text(.02,.9,letter,
        horizontalalignment='left',family="sans-serif",weight="heavy",
        transform=b.transAxes)

        # b.set_ylim(ymax=.7)

fig.tight_layout(pad=0.4)
fig.subplots_adjust(hspace=0.08,wspace=0.05,right=0.99,left=0.09)
plt.show()



