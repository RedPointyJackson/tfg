#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import pandas as pd

def cm2inch(value): return value/2.54

plt.style.use('custom')


colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
primary = colors[0]

jdf = pd.read_csv('data_jack.csv')
bdf = pd.read_csv('data_boot.csv')

# General

fig = plt.figure(figsize=(cm2inch(12), cm2inch(6)))
ax1 = fig.add_subplot(131)

ax1.errorbar(x=bdf['i'], y=bdf['μ'], yerr=bdf['σ'], fmt='o',
             color=primary, ms=3, lw=1, capsize=3)

ax1.set_xlabel('Bootstrap')
ax1.set_ylabel('μ ± σ')


ax2 = fig.add_subplot(132)

ax2.errorbar(x=jdf['i'], y=jdf['μ'], yerr=jdf['σ'], fmt='o',
             color=primary, ms=3, lw=1, capsize=3)

ax2.set_xlabel('Jacknife')
ax2.set_ylabel('')


ax1.set_xlim(0,15.5)
ax2.set_xlim(0,15.5)

ax1.set_ylim(-2,10)
ax2.set_ylim(-2,10)

ax1.set_xticks([])
ax2.set_xticks([])

ax1.set_yticks([-2,0,1,2,4,6,8,10])
ax2.set_yticks([-2,0,1,2,4,6,8,10])
ax2.set_yticklabels([])

# Comparison

ax3 = fig.add_subplot(133)

xx = jdf['μ']
yy = bdf['μ']

ρ = np.mean( (xx-np.mean(xx)) * (yy-np.mean(yy)) )/np.std(xx)/np.std(yy)
print("ρ =",ρ)

ax3.scatter(xx, yy, alpha=0.3, edgecolor='none')

ax3.set_xlabel('μ (Jacknife)')
ax3.set_ylabel('μ (Bootstrap)')
ax3.set_xlim(0,6)
ax3.set_ylim(0,6)

fig.savefig('general.pdf')
