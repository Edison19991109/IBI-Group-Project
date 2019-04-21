# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 17:30:29 2019

@author: Valerya
"""
import matplotlib.pyplot as plt
from numpy import log

P = {'A':(6.8,5500),'B':(8.7,550),'C':(4.8,7890)}
plt.figure(figsize=(10,12))
for k in P:
    pH=P[k][0]
    Mw=10-log(P[k][1]/1000)
    plt.plot(pH,Mw,'o',label=k+":pH=%d,mass=%d"%(P[k][0],P[k][1]))
plt.title('2-DE')
plt.xticks([])
plt.yticks([])
plt.xlim(4,10)
plt.ylim(10-log(10),)
plt.legend(loc='lower right')
plt.show()
