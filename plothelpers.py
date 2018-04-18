#######################################################################
#filename: 'plothelpers.py'
#Library with useful functions for plotting.
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

import matplotlib
import matplotlib.pyplot
import itertools
import numpy as np
import math
def latexify(fig=None,fig_width=None, fig_height=None, columns=1):
    """ Sets standard parameters of matplotlib. Call before plotting.
    adapted from https://nipunbatra.github.io/blog/2014/latexify.html

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 2.825 if columns==1 else 5.788 # width in inches
        # fig_width = 3.38 if columns==1 else 7. # width in inches
        # fig_width = 3.176 if columns==1 else 6.491 # width in inches
        # fig_width = 3.39 if columns==1 else 6.9 # width in inches
        # 1 inch= 2.54 cm

    if fig_height is None:
        golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height + 
                "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {#'backend': 'ps',
            'axes.labelsize': 9, # fontsize for x and y labels (was 10)
            'axes.titlesize': 9,
            'font.size': 10, # was 10
            'legend.fontsize': 8, # was 10
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
            'text.usetex': True,
            'figure.figsize': [fig_width,fig_height],
            'font.family': 'sans-serif',
            'font.sans-serif': ['Helvetica'],#['computer modern roman'], #avoid bold axis label
            'text.latex.preamble': [r'\usepackage{helvet}',# set the normal font here
                r'\usepackage[EULERGREEK]{sansmath}',  # load up the sansmath so that math -> helvet
                r'\sansmath'  # <- tricky! -- gotta actually tell tex to use!
                ] 
            }
    if fig:
        print("texify figure dimensions set: ",fig_width,fig_height)
        fig.set_size_inches((fig_width,fig_height),forward=True)
    matplotlib.rcParams.update(params)
    return params


def set_axes_size(width=None,height=None, ax=None):
    """ width, height in inches """

    if width is None: width=2.625 

    if height is None:
        golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
        height = width*golden_mean # height in inches
    if ax is None: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    prevfigsize=ax.figure.get_size_inches()
    prevfigw=prevfigsize[0]
    prevfigh=prevfigsize[1]

    figw=float(width)+(l+1-r)*prevfigw
    figh=float(height)+(b+1-t)*prevfigh
    newl=l*prevfigw/figw
    newr=1-(1-r)*prevfigw/figw
    newb=b*prevfigh/figh
    newt=1-(1-t)*prevfigh/figh

    ax.figure.set_size_inches(figw, figh,forward=True)
    ax.figure.subplots_adjust(left=newl,right=newr,top=newt,bottom=newb)


def create_colorcyle(number,cmap=None,cmapname="viridis"):
    if not cmap:
        cmap = matplotlib.pyplot.get_cmap(cmapname)
    indices = np.linspace(0, cmap.N, number)
    my_colors = itertools.cycle([cmap(int(i)) for i in indices])
    return my_colors


def plot_datacap(ax,x,y,xint=None,yint=None,color="black",lw=0.8,offset=None):
    '''plots two short diagonal lines to denote capping of
       data yaxis. 
       x,y: (center) position
       xint,yint: interval taken up by lines
    '''

    if xint is None:
        xint=1
    if yint is None:
        yint=1

    xint=xint/2.
    yint=yint/2.
    if offset is None:
        offset=yint
    steps=20
    xvals=np.linspace(x-xint,x+xint,steps)
    yvals=np.linspace(y+yint,y-yint,steps)

    ax.plot(xvals,yvals,color=color,lw=lw,zorder=5)
    ax.plot(xvals,yvals+offset,color=color,lw=lw,zorder=5)

    vertx=[xvals[0],xvals[0],xvals[-1],xvals[-1]]
    verty=[yvals[0],yvals[0]+offset,yvals[-1]+offset,yvals[-1]]
    xy=np.vstack([vertx,verty]).T
    # print(xy)
    patch=matplotlib.patches.Polygon(xy,facecolor='white',zorder=4)
    ax.add_patch(patch)
