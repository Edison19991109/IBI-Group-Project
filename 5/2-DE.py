# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 17:30:29 2019

@author: Valerya
"""
import matplotlib.pyplot as plt
from numpy import log

P = {'A':(6.8,5500),'B':(8.7,550),'C':(4.8,7890)}
fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111)
for k in P:
    pH=P[k][0]
    Mw=10-log(P[k][1]/1000)
    ax.plot(pH,Mw,'o',label=k+":pI=%d,Mw=%d"%(P[k][0],P[k][1]))
plt.title('2-DE')
ax.set_xlabel('pH')
ax.set_ylabel('Mw(Da)')
ax.set_xticks(range(4,10))
ax.set_yticks([10-log(10),10-log(5),10-log(2),10-log(1.5),10-log(1),10-log(0.5),10-log(0.2),10-log(0.1)])
ax.set_yticklabels([10000,5000,2000,1500,1000,500,200,100])
ax.set_xlim(4,10)
ax.set_ylim(10-log(10),)
ax.legend(loc='lower right')
plt.show()
