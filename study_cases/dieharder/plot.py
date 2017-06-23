#!/usr/bin/env python3

def cm2inch(value):
    return value/2.54

import matplotlib.pyplot as plt
import pandas as pd
import seaborn.apionly as sns

plt.style.use('custom')

failcolor = '#C44E52'
passcolor = '#55A868'
warncolor = '#FFA574'

df = pd.read_csv('data.csv')

fig = plt.figure(figsize=(cm2inch(15),cm2inch(6)))

ax = sns.stripplot(x='generator', y='value', data=df, jitter=True)
ax.set_xlabel('Generador')
ax.set_ylabel('p values')

fig.savefig('summary.pdf')
# for gen,gendf in df.groupby('generator'):
#     fig = plt.figure(figsize=(cm2inch(2),cm2inch(2)))
#     ax = fig.add_subplot(1,1,1)
#     ax.set_ylabel('')
#     ax.set_xlabel('')
#     L = len(gendf['value'])
#     vals = list(gendf['value'])
#     for j in range(L):
#         val = vals[j]
#         color = passcolor
#         marker = 'o'
#         ps = 2
#         if val < 0.01: color=failcolor; marker='x'; ps=5
#         if val > 0.997: color=warncolor; marker='x'; ps=5
#         ax.plot(j, val, marker, ms=ps, color=color)
#     ax.set_xticks([])
#     ax.set_yticks([0,0.5,1])
#     ax.set_yticklabels([])
#     fig.savefig('%s.pdf' % gen)
