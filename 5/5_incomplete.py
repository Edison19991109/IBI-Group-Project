# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 20:20:20 2019

@author: Valerya
"""

#现实中同种氨基酸R基在一条多肽链不同位置pKa会发生变化，此程序无法模拟。当同种氨基酸在一条多肽链中出现次数过多，可能无法计算pI
L = {'A':'PDAKVL',
     'B':'EPHINAQIM',
     'C':'LGSRQKHSLPDLPYDYGAL',
     'D':'MLSRAVCGTSRQLAPVLAY',
     'E':'INHSIFWTNLSPNGGGEPKGELLEAIKRDFGSFDKFKE',
     'F':'QLHHSKADAAYVNNLNVTVEAYQEALAKGDVTAQIALQPALKFNGGGV'}

A ={'G':(75,2.34,9.60),'A':(89,2.34,9.69),'P':(115,1.99,10.96),'V':(117,2.32,9.62),
    'L':(131,2.36,9.60),'I':(131,2.36,9.68),'M':(149,2.28,9.21),'F':(165,1.83,9.13),
    'Y':(181,2.20,9.11,10.07,0),'W':(204,2.38,9.39),'S':(105,2.21,9.15),'T':(119,2.11,9.62),
    'C':(121,1.96,10.28,8.18,0),'N':(132,2.02,8.80),'Q':(146,2.17,9.13),'K':(146,2.18,8.95,10.53,1),
    'H':(155,1.82,9.17,6.00,1),'R':(174,2.17,9.04,12.48,1),'D':(133,1.88,9.60,3.65,0),'E':(147,2.19,9.67,4.25,0)}
#Nelson, David L.(David Lee),1942-author.  Lehninger principles of biochemistry /. 7th ed. 
#New York : W. H. Freeman and Company ; Houndmills, Basingstoke : Macmillan Higher Education, c2017.
P = {}

def count(L,a,b):
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
    Mw = 0
    pI = 0
    pKa = [A[L[k][0]][1]]
    pKb = [A[L[k][len(L[k])-1]][2]]
    for i in range(len(L[k])):
        Mw += A[L[k][i]][0]
        if len(A[L[k][i]]) == 5:
            if A[L[k][i]][4] == 0:
                pKa.append(A[L[k][i]][3])
            else:
                pKb.append(A[L[k][i]][3])
    pKa.sort(reverse=True)
    pKb.sort()
    print(pKa,'\n',pKb)
    for i in range(min([len(pKa),len(pKb)])-1):
        pH = (pKa[i]+pKb[i])/2
        if count(pKa,pH,0) == count(pKb,pH,1):
            pI = pH
            print(pKa[i],pKb[i])
            break
    P[k] = (pI,Mw-18*(len(L[k])-1))

import matplotlib.pyplot as plt
from numpy import log

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111)
for k in P:
    pI=P[k][0]
    Mw=10-log(P[k][1]/1000)
    ax.plot(pI,Mw,'o',label=k+":pI=%d,Mw=%d"%(P[k][0],P[k][1]))
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
