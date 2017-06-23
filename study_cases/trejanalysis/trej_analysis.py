#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn.apionly as sns
import scipy.optimize as opt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid.inset_locator import inset_axes

def cm2inch(value):
    return value/2.54


plt.style.use('custom')

import matplotlib as mpl
colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
primary = colors[0]
secondary = colors[1]
tertiary = colors[3]


lsc = mpl.colors.LinearSegmentedColormap

niceblue = '#64B5CD'
nicered = '#C44E52'
niceyellow = '#CCB974'

colormap = lsc.from_list('foo',[(0.0,niceblue),
                                (0.5,niceyellow),
                                (1.0,nicered)])

maxtw1, mintw1 = 1e8, 1e4

def get_color_tw1(tw1):
    lmaxtw1 = np.log10(maxtw1)
    lmintw1 = np.log10(mintw1)
    τ = (np.log10(tw1)-lmintw1)/(lmaxtw1-lmintw1)
    return colormap(τ)



mementodf = pd.read_csv('memento_rejuvenation.csv')
janusdf = pd.read_csv('janus_rejuvenation.csv')
trejdf = pd.concat([mementodf,janusdf], ignore_index=True)
trejdf['Trel'] = (trejdf['T1']-trejdf['T2'])/trejdf['T1']

# tw1=10⁸ can't be compatible with t=10^8. Too near.
ingoodtw1 = trejdf['tw1'] < 9e7
ingoodt = trejdf['t'] < 9e7
trejdf = trejdf[ingoodtw1 | ingoodt]


fig = plt.figure(figsize=(cm2inch(15),cm2inch(20)))

trejdf8  = trejdf[trejdf['L'] == 8]
trejdf12 = trejdf[trejdf['L'] == 12]
trejdf48 = trejdf[trejdf['L'] == 48]
trejdf80 = trejdf[trejdf['L'] == 80]


######################
#                    #
#    Left column     #
#                    #
######################

################## x = t

ax = fig.add_subplot(4,3,1)
ax.set_xscale('log')

for p,pdf in trejdf8.groupby(['tw1','Trel']):
    tw1, Trel = p
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+1.1*pdf['Δχ_σ'] ,pdf['Δχ_μ']-1.1*pdf['Δχ_σ'] , color='w')
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+pdf['Δχ_σ']     ,pdf['Δχ_μ']-pdf['Δχ_σ']     , color=get_color_tw1(tw1))

ax.set_xticks([1e4,1e6,1e8])
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylabel('Δχ')
ax.set_ylim(-0.6,0.1)


ax = fig.add_subplot(4,3,4)
ax.set_xscale('log')

for p,pdf in trejdf12.groupby(['tw1','Trel']):
    tw1, Trel = p
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+1.1*pdf['Δχ_σ'] ,pdf['Δχ_μ']-1.1*pdf['Δχ_σ'] , color='w')
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+pdf['Δχ_σ']     ,pdf['Δχ_μ']-pdf['Δχ_σ']     , color=get_color_tw1(tw1))


ax.set_xticks([1e4,1e6,1e8])
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylabel('Δχ')
ax.set_ylim(-0.4,0.1)

##

ax = fig.add_subplot(4,3,7)
ax.set_xscale('log')

for p,pdf in trejdf48.groupby(['tw1','Trel']):
    tw1, Trel = p
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+1.1*pdf['Δχ_σ'] ,pdf['Δχ_μ']-1.1*pdf['Δχ_σ'] , color='w')
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+pdf['Δχ_σ']     ,pdf['Δχ_μ']-pdf['Δχ_σ']     , color=get_color_tw1(tw1))

ax.set_xticks([1e4,1e6,1e8])
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylabel('Δχ')
ax.set_ylim(-0.3,0.1)

##

ax = fig.add_subplot(4,3,10)
ax.set_xscale('log')

for p,pdf in trejdf80.groupby(['tw1','Trel']):
    tw1, Trel = p
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+1.1*pdf['Δχ_σ'] ,pdf['Δχ_μ']-1.1*pdf['Δχ_σ'] , color='w')
    ax.fill_between(pdf['t'] ,pdf['Δχ_μ']+pdf['Δχ_σ']     ,pdf['Δχ_μ']-pdf['Δχ_σ']     , color=get_color_tw1(tw1))

ax.set_xticks([1e4,1e6,1e8])
ax.set_xlabel('$t_0$')
ax.set_ylabel('Δχ')
ax.set_ylim(-0.3,0.1)

######################
#                    #
#    Middle column   #
#                    #
######################

################## x = t

ax = fig.add_subplot(4,3,2)
ax.set_xscale('log')

for p,pdf in trejdf8.groupby(['tw1']):
    df = pdf.groupby('t').mean().reset_index()
    ax.plot(df['t'], df['Δχ_μ'], color=get_color_tw1(p))

ax.set_xticks([1e4,1e6,1e8])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_ylim(-0.6,0.1)

# Fake legend
xx = np.logspace(np.log10(1e4), np.log10(1e6), 3)
yy = [-0.4,-0.50]
ax.text(xx[0] ,yy[0] ,'    $t_w^1:$' ,backgroundcolor='w' , color='gray')
ax.text(xx[1] ,yy[0] ,'$● 10^4$'     ,backgroundcolor='w' , color=get_color_tw1(1e4))
ax.text(xx[2] ,yy[0] ,'$● 10^5$'     ,backgroundcolor='w' , color=get_color_tw1(1e5))
ax.text(xx[0] ,yy[1] ,'$● 10^6$'     ,backgroundcolor='w' , color=get_color_tw1(1e6))
ax.text(xx[1] ,yy[1] ,'$● 10^7$'     ,backgroundcolor='w' , color=get_color_tw1(1e7))
ax.text(xx[2] ,yy[1] ,'$● 10^8$'     ,backgroundcolor='w' , color=get_color_tw1(1e8))


##

ax = fig.add_subplot(4,3,5)
ax.set_xscale('log')

for p,pdf in trejdf12.groupby(['tw1']):
    df = pdf.groupby('t').mean().reset_index()
    ax.plot(df['t'], df['Δχ_μ'], color=get_color_tw1(p))

ax.set_xticks([1e4,1e6,1e8])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_ylim(-0.4,0.1)

##

ax = fig.add_subplot(4,3,8)
ax.set_xscale('log')

for p,pdf in trejdf48.groupby(['tw1']):
    df = pdf.groupby('t').mean().reset_index()
    ax.plot(df['t'], df['Δχ_μ'], color=get_color_tw1(p))

ax.set_xticks([1e4,1e6,1e8])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_ylim(-0.3,0.1)

##

ax = fig.add_subplot(4,3,11)
ax.set_xscale('log')

for p,pdf in trejdf80.groupby(['tw1']):
    df = pdf.groupby('t').mean().reset_index()
    ax.plot(df['t'], df['Δχ_μ'], color=get_color_tw1(p))

ax.set_xticks([1e4,1e6,1e8])
ax.set_yticklabels([])
ax.set_xlabel('$t_0$')
ax.set_ylabel('')
ax.set_ylim(-0.3,0.1)

##


######################
#                    #
#     Third column   #
#                    #
######################

################## x = (T₂-T₁)/T₁

ax = fig.add_subplot(4,3,3)
ax.set_xscale('linear')

for p,pdf in trejdf8.groupby(['t','tw1']):
    t,tw1 = p
    ax.scatter(pdf['Trel'],pdf['Δχ_μ'], s=3, color=get_color_tw1(tw1))

ax.set_xlabel('')
ax.set_xticks([0.4,0.2])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_ylabel('')
ax.set_ylim(-0.6,0.1)
ax2 = ax.twinx()
ax2.set_ylabel('L=8')
ax2.set_yticks([])

##

ax = fig.add_subplot(4,3,6)
ax.set_xscale('linear')

for p,pdf in trejdf12.groupby(['t','tw1']):
    t,tw1 = p
    ax.scatter(pdf['Trel'],pdf['Δχ_μ'], s=3, color=get_color_tw1(tw1))

ax.set_xlabel('')
ax.set_xticks([0.4,0.2])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_ylabel('')
ax.set_ylim(-0.4,0.1)
ax2 = ax.twinx()
ax2.set_ylabel('L=12')
ax2.set_yticks([])

##

ax = fig.add_subplot(4,3,9)
ax.set_xscale('linear')

for p,pdf in trejdf48.groupby(['t','tw1']):
    t,tw1 = p
    ax.scatter(pdf['Trel'],pdf['Δχ_μ'], s=3, color=get_color_tw1(tw1))

ax.set_xlabel('')
ax.set_xticks([0.4,0.2])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_ylabel('')
ax.set_ylim(-0.3,0.1)
ax2 = ax.twinx()
ax2.set_ylabel('L=48')
ax2.set_yticks([])



##

ax = fig.add_subplot(4,3,12)
ax.set_xscale('linear')

for p,pdf in trejdf80.groupby(['t','tw1']):
    t,tw1 = p
    ax.scatter(pdf['Trel'],pdf['Δχ_μ'], s=3, color=get_color_tw1(tw1))

ax.set_xlabel('$(T_1-T_2)/T_1$')
ax.set_xticks([0.4,0.2])
ax.set_yticklabels([])
ax.set_ylabel('')
ax.set_ylim(-0.3,0.1)
ax2 = ax.twinx()
ax2.set_ylabel('L=80')
ax2.set_yticks([])


fig.savefig('1D_dependence.pdf')
