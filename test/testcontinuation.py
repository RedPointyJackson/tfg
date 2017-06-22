#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn.apionly as sns

def cm2inch(value): return value/2.54

plt.style.use('custom')
dfA = pd.read_csv('/tmp/A.csv')
dfB = pd.read_csv('/tmp/B.csv')
dfB['mc'] += 10
dfC = pd.read_csv('/tmp/C.csv')
dfC['mc'] += 110

dfA['kind'] = "Direct"
dfB['kind'] = "Continuation A"
dfC['kind'] = "Continuation B"

alldf = dfA.append(dfB).append(dfC)

alldfE = alldf[alldf['observable'] == 'Energy']
alldfC = alldf[alldf['observable'] == 'Correlation']

fig = plt.figure( figsize=(cm2inch(10), cm2inch(5)) )
axE = fig.add_subplot(211)
axC = fig.add_subplot(212)

axE.set_xscale('log')
axE.legend()
axE.set_xlabel('mc')
axE.set_ylabel('E')

axC.set_xscale('log')
axC.legend()
axC.set_xlabel('mc')
axC.set_ylabel('β(1-C)')

for k,kdf in alldfE.groupby('kind'):
    axE.plot(kdf['mc'], kdf['value'], label=k, lw=1)

axE.legend()

for k,kdf in alldfC.groupby('kind'):
    axC.scatter(kdf['mc'], kdf['value'], label=k, s=1)

plt.show()
