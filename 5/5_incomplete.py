# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 20:20:20 2019

@author: Valerya
"""

#现实中同种氨基酸R基在一条多肽链不同位置pKa会发生变化，此程序无法模拟。当同种氨基酸在一条多肽链中出现次数过多，可能无法计算pI
L = {'A':'E',
     'B':'AAKA',
     'C':'EPHINAQIM',
     'D':'LGSRQKHSLPDLPYDYGAL',
     'E':'INHSIFWTNLSPNGGGEPKGELLEAIKRDFGSFCKFKE',
     'F':'MLCRAACSTGRRLGPVAGAAGSRHKHSLPDLPYDYGALEPHINAQIMQL'} #可以换成input输入多肽

A ={'G':(75,2.34,9.60),'A':(89,2.34,9.69),'P':(115,1.99,10.96),'V':(117,2.32,9.62),
    'L':(131,2.36,9.60),'I':(131,2.36,9.68),'M':(149,2.28,9.21),'F':(165,1.83,9.13),
    'Y':(181,2.20,9.11,10.07,0),'W':(204,2.38,9.39),'S':(105,2.21,9.15),'T':(119,2.11,9.62),
    'C':(121,1.96,10.28,8.18,0),'N':(132,2.02,8.80),'Q':(146,2.17,9.13),'K':(146,2.18,8.95,10.53,1),
    'H':(155,1.82,9.17,6.00,1),'R':(174,2.17,9.04,12.48,1),'D':(133,1.88,9.60,3.65,0),'E':(147,2.19,9.67,4.25,0)}
#Nelson, David L.(David Lee),1942-author.  Lehninger principles of biochemistry /. 7th ed. 
#New York : W. H. Freeman and Company ; Houndmills, Basingstoke : Macmillan Higher Education, c2017.
P = {} #这个字典用来存多肽信息

def count(L,a,b): #数带正/负电氨基酸残基
    x = 0
    for i in L:
        if b == 1:
            if i<a:
                x += 1    
        else:
            if i>a:
                x += 1
    return x
    

for k in L:
    Mw = 0 #分子量
    pI = 0 #等电点
    pKa = [A[L[k][0]][1]] #C端羧基 3HN+/COOH <=A[k][1]=> 3HN+/COO-
    pKb = [A[L[k][len(L[k])-1]][2]] #N端氨基 3HN+/COO- <=A[k][2]=> 2HN/COO-
    for i in range(len(L[k])):
        Mw += A[L[k][i]][0] #加上当前氨基酸分子量
        if len(A[L[k][i]]) == 5: #若当前氨基酸R基可电离
            if A[L[k][i]][4] == 0: #R <=A[k][3]=> R-
                pKa.append(A[L[k][i]][3]) #pKa：pH升高时带负电->不带电
            else: #R+ <=A[k][3]=> R
                pKb.append(A[L[k][i]][3]) #pKa：pH升高时不带电->带正电
    pK = pKa+pKb
    pK.sort()
    ###print(pKa,'\n',pKb,'\n',pK)
    for i in range(len(pK)-1):
        pH = (pK[i]+pK[i+1])/2
        if count(pKa,pH,1) == count(pKb,pH,0) > 0: #还要再讨论：某些多肽链中同种氨基酸多次出现，同时电离导致pH>某点时多肽带负电，<某点时带负电的情况
            pI = float('%.2f' % pH)
            ###print(pK[i],pK[i+1])
            break
        if count(pKa,pH,1) > count(pKb,pH,0):
            pI = pK[i]/2+(pK[i-1]+pK[i+1])/4
            pI = float('%.2f' % pI)
            ###print(pK[i],pK[i+1],pK[i-1])
            break
    P[k] = (pI,Mw-18*(len(L[k])-1))

import pandas as pd
import matplotlib.pyplot as plt
from numpy import log

df = pd.DataFrame(P)
df = df.T
df.columns = ['pI','Mw']
print(df)

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111)
for k in P:
    pI=P[k][0]
    Mw=10-log(P[k][1]/1000) #logMw与距离成线性关系，Mw越大距离越短
    ax.plot(pI,Mw,'o',label=k+':pI='+str(P[k][0])+' Mw='+str(P[k][1]))
plt.title('2-DE')
ax.set_xlabel('pH')
ax.set_ylabel('Mw(Da)')
ax.set_xticks(range(3,112))
ax.set_yticks([10-log(10),10-log(5),10-log(2),10-log(1.5),10-log(1),10-log(0.5),10-log(0.2),10-log(0.1)])
ax.set_yticklabels([10000,5000,2000,1500,1000,500,200,100])
ax.set_xlim(3,11)
ax.set_ylim(10-log(10),)
ax.legend(loc='upper right')
plt.show()
