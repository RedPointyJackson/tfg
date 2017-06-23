#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import scipy.optimize as opt

def cm2inch(value): return value/2.54

import numpy as np
import pandas as pd

lsc = mpl.colors.LinearSegmentedColormap

niceblue = '#64B5CD'
nicered = '#C44E52'
niceyellow = '#CCB974'

plt.style.use('custom')

data = pd.read_csv('measures.csv').groupby(['mc','parameter']).mean().reset_index()
data = data.rename(columns={'mc':'tw', 'parameter':'t', 'value':'C'})

colormap = lsc.from_list('foo',[(0.0,niceblue),
                                (0.5,niceyellow),
                                (1.0,nicered)])

maxt, mint = data['t'].max(), data['t'].min()

def get_color(t):
    τ = (t-mint)/(maxt-mint)
    return colormap(τ)

fig = plt.figure(figsize=(cm2inch(15), cm2inch(8)))
ax = fig.add_subplot(111)




###################################
#
#           Fit
#
###################################


tlist = []
alist = []
blist = []
clist = []

mintw = 2e3
maxtw = 3e8
def fitfunc(x, a, b, c):
    return a+b*x**(-c)


for t,dft in data.groupby('t'):
    goodpoints = [mintw<tw<maxtw for tw in dft['tw']]
    x = dft['tw'][goodpoints]
    y = dft['C'][goodpoints]
    popt, pcov = opt.curve_fit(fitfunc,x,y,
                               p0=(0.7,-1.5,0.25),
                               bounds=[(0.2,-5,0),(0.8,0,0.5)])
    perr = np.sqrt(np.diag(pcov))
    tlist.append(t)
    alist.append((popt[0],perr[0]))
    blist.append((popt[1],perr[1]))
    clist.append((popt[2],perr[2]))


###################################
#
#        General plot
#
###################################


for t,dft in data.groupby('t'):
    ax.plot(dft['tw'], dft['C'],
            linewidth=0.5,
            color=get_color(t))

fullxx = np.logspace(np.log10(data['tw'].min()),
                     np.log10(data['tw'].max()),
                     1000)
xx = [x for x in fullxx if mintw<x<maxtw]


# Overlay a fit for each lowest t
for tidx,t in enumerate(tlist):
    popt = (alist[tidx][0], blist[tidx][0], clist[tidx][0])
    perr = (alist[tidx][1], blist[tidx][1], clist[tidx][1])
    # ax.plot(fullxx, fitfunc(fullxx,*popt), # white background for the line
    #         linewidth=2,
    #         color='white')
    ax.plot(fullxx, fitfunc(fullxx,*popt), '--',
            linewidth=0.5,
            alpha=0.55,
            color=get_color(t))
    ax.plot(xx, fitfunc(xx,*popt),
            linewidth=0.5,
            color=get_color(t))


ax.set_ylabel(r'$[C(t_w,t_w+t)]$')
ax.set_xlabel(r'$t_w$')

ax.set_xscale('log')
ax.set_ylim(0,0.8)
ax.set_xlim(left=1)

ax.grid(c='#E5E5E5')

ax.set_xticks([1e2,1e4,1e6,1e8])

# Create a fake colorbar
def fmt(x, pos):
    return '%d' % (np.floor(x/1e7))

sm = cm.ScalarMappable(cmap=colormap,norm=plt.Normalize(vmin=0,vmax=maxt))
sm.set_array([])


cbaxes = fig.add_axes([0.10, 0.8, 0.3, 0.01])
cbar = plt.colorbar(sm,
                    format=mpl.ticker.FuncFormatter(fmt),
                    orientation='horizontal',
                    cax=cbaxes)
cbar.outline.set_visible(False)


cbar.ax.set_title('t [×10⁷]')
cbar.ax.title.set_size(12)
cbar.ax.title.set_color('#555555')

fig.savefig('fit.pdf')

# Now, the parameters
pfig = plt.figure(figsize=(cm2inch(15), cm2inch(8)))
mini = plt.figure(figsize=(cm2inch(5), cm2inch(5)))

a_ax = pfig.add_subplot(3,1,1)
b_ax = pfig.add_subplot(3,1,2)
c_ax = pfig.add_subplot(3,1,3)
a_ax.set_xscale('log')
b_ax.set_xscale('log')
c_ax.set_xscale('log')
miniax = mini.add_subplot(1,1,1)
miniax.set_xscale('log')


miniax.set_ylabel('a')
miniax.set_xlabel('t')
# miniax.set_ylim(0.6,0.8)
miniax.grid(c='#E5E5E5')

a_ax.set_xticklabels([])
a_ax.set_ylabel('a')
# a_ax.set_ylim(0.6,0.8)
a_ax.grid(c='#E5E5E5')

b_ax.set_xticklabels([])
b_ax.set_ylabel('b')
# b_ax.set_ylim(-2,0)
b_ax.grid(c='#E5E5E5')

c_ax.set_ylabel('c')
c_ax.set_xlabel('t')
# c_ax.set_ylim(0,0.5)
c_ax.grid(c='#E5E5E5')

for i in range(len(tlist)):
    t = tlist[i]
    a,a_err = alist[i]
    b,b_err = blist[i]
    c,c_err = clist[i]
    a_ax.errorbar(t,a,fmt='o-',yerr=a_err, ms=1, capsize= 2, lw=0.5, capthick=0.5, color=niceblue)
    miniax.errorbar(t,a,fmt='o-',yerr=a_err, ms=1, capsize= 2, lw=0.5, capthick=0.5, color=niceblue)
    b_ax.errorbar(t,b,fmt='o-',yerr=b_err, ms=1, capsize= 2, lw=0.5, capthick=0.5, color=niceblue)
    c_ax.errorbar(t,c,fmt='o-',yerr=c_err, ms=1, capsize= 2, lw=0.5, capthick=0.5, color=niceblue)


# fmt = mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False)
# fmt.set_powerlimits((-3,3))
# c_ax.xaxis.set_major_formatter(fmt)
# miniax.xaxis.set_major_formatter(fmt)

pfig.savefig('params.pdf')
mini.savefig('params_minia.pdf')
