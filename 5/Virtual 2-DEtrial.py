# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 20:20:20 2019

@author: Valerya
"""

import re
def Cap(s):
    c = ''
    for char in s:
        if re.match('[A-Za-z]',char):
            c += char.upper()
    return c #去数字空格，全部大写

#现实中同种氨基酸R基在一条多肽链不同位置pKa会发生变化，此程序无法模拟。当同种氨基酸在一条多肽链中出现次数过多，可能无法计算pI
L = {'Oxytocin':'CYIQNCPLG', #催产素
     'Vasopressin':'CYFQNCPRG', #血管加压素
     'InsulinA':'GIVEQCCTSICSLYQLENYCN', #胰岛素A链
     'Neurophysin2':'AVLDLDVRTCLPCGPGGKGRCFGPSICCGDELGCFVGTAEALRCQEENYLPSPCQSGQKPCGSGGRCAAAGICCSPDGCHEDPACDPEAAFSQH', #后叶激素运载蛋白2
     'CytochromeC':'MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE', #细胞色素C
     'RndSeq1':Cap(input('input a sequence of a (30 AA) peptide:\n')), #输入随机序列，转大写
     'RndSeq2':Cap(input('input a sequence of a (50 AA) peptide:\n'))}

A ={'G':(4,75.052,2.34,9.60),'A':(1,89.079,2.34,9.69),'P':(1,115.117,1.99,10.96),'V':(1,117.133,2.32,9.62),
    'L':(1,131.160,2.36,9.60),'I':(1,131.160,2.36,9.68),'M':(1,149.199,2.28,9.21),'F':(1,165.177,1.83,9.13),
    'Y':(4,181.176,2.20,9.11,10.07,0),'W':(1,204.213,2.38,9.39),'S':(4,105.078,2.21,9.15),'T':(4,119.105,2.11,9.62),
    'C':(4,121.145,1.96,10.28,8.18,0),'N':(4,132.104,2.02,8.80),'Q':(4,146.131,2.17,9.13),'K':(2,146.170,2.18,8.95,10.53,1),
    'H':(2,155.141,1.82,9.17,6.00,1),'R':(2,174.188,2.17,9.04,12.48,1),'D':(3,133.089,1.88,9.60,3.65,0),'E':(3,147.116,2.19,9.67,4.25,0)}
#Nelson, David L.(David Lee),1942-author.  Lehninger principles of biochemistry /. 7th ed. 
#New York : W. H. Freeman and Company ; Houndmills, Basingstoke : Macmillan Higher Education, c2017.
P = {} #这个字典用来存多肽信息

def count(L,a,b): #数带正/负电氨基酸残基
    x = 0
    for i in L:
        if b == 1:
            x += 10**a/(10**a + 10**i)
        else:
            x +=10**i/(10**a + 10**i)
    return x

if re.search(r'[BJOUXZ]',L['RndSeq1']) or re.search(r'[BJOUXZ]',L['RndSeq2']):
    print('\nInvalid sequence')
else:
    print('\nAmino acids composition:')
    print('\033[33mNonpolar uncharged(NU), \033[36mPositively charged(+C), \033[32mNegatively charged(-C), \033[37mPolar uncharged(PU)\033[0m R groups') #红，蓝，绿，灰
    for k in L: #for key in L
        Mw = pI = a = b = c = d = 0 #分子量 #等电点
        pKa = [A[L[k][0]][2]] #C端羧基 3HN+/COOH <=A[k][2]=> 3HN+/COO-
        pKb = [A[L[k][len(L[k])-1]][3]] #N端氨基 3HN+/COO- <=A[k][3]=> 2HN/COO-
        print(k+': H-',end='') #不换行
        for i in range(len(L[k])):
            if A[L[k][i]][0] == 1: a += 1; print('\033[33m',end=''); print(L[k][i],end='') #A = {AminoAcid:(class,Mw,C端羧基pk，N端氨基pk，R基pk，R基类型)
            elif A[L[k][i]][0] == 2: b += 1; print('\033[36m',end=''); print(L[k][i],end='')
            elif A[L[k][i]][0] == 3: c += 1;print('\033[32m',end=''); print(L[k][i],end='')
            else: d += 1; print('\033[37m',end=''); print(L[k][i],end='') #class=1红色, class=2蓝色，class=3绿色，class=4灰色
            Mw += A[L[k][i]][1] #加上当前氨基酸分子量
            if len(A[L[k][i]]) == 6: #若当前氨基酸R基可电离
                if A[L[k][i]][5] == 0: #R <=A[k][4]=> R-
                    pKa.append(A[L[k][i]][4]) #pKa：pH升高时不带电->带负电
                else: #R+ <=A[k][3]=> R
                    pKb.append(A[L[k][i]][4]) #pKb：pH升高时带正电->不带电
        print('\033[30m-OH') #换行
        pK = pKa+pKb #pK：所有解离常数
        pK.sort() #从小到大排列
        ###print(pKa,'\n',pKb,'\n',pK) ###for code used for checking
        Charge = {}
        for i in range(len(pK)-1):
            pH = (pK[i]+pK[i+1])/2 #响铃两解离常数中位数， 可能的pI
            pH = float('%.2f' % pH) #2位小数
            Charge[pH] = abs(count(pKb,pH,0) - count(pKa,pH,1)) #还要再讨论：某些多肽链中同种氨基酸多次出现，同时电离导致pH>某点时多肽带负电，<某点时带负电的情况
        pI = min(zip(Charge.values(),Charge.keys()))[1]
        ###print(min(zip(Charge.values(),Charge.keys()))[0])
        Mw = Mw-18.015*(len(L[k])-1) #氨基酸总质量减脱水质量
        P[k] = (len(L[k]),float('%.1f' % Mw),pI,'%.1f%%' % (b*100/len(L[k])),'%.1f%%' % (c*100/len(L[k])),'%.1f%%' % (a*100/len(L[k])),'%.1f%%' % (d*100/len(L[k])))
        #P = {Peptide:(NumberOfResidues, MolecularWeight, IsoelectricPoint, +CPercentage(一位小数), -CPercentage, NUPercentage, PUPercentage)}
    import pandas as pd 
    import matplotlib.pyplot as plt
    from numpy import log

    df = pd.DataFrame(P)
    df = df.T #行列置换
    df.columns = ['Residue','Mw','pI','+C','-C','NU','PU']
    pd.set_option('display.max_columns', None) #print dataframe 不省略中间的列
    print('\nPeptide properties:'); print(df)

    fig = plt.figure(figsize=(10,12)) #图像大小
    ax = fig.add_subplot(111) #子图功能更多
    for k in P:
        pI=P[k][2]
        Mw=-log(P[k][1]/1000) #logMw与距离成线性关系，Mw越大距离越短
        ax.plot(pI,Mw,'o',label=k)
    plt.title('Virtual Two-dimensional Electrophoresis')
    ax.set_xlabel('pH')
    ax.set_ylabel('Mw(Da)')
    ax.set_xticks(range(3,11)) #x轴标刻度
    ax.set_yticks([-log(20),-log(10),-log(5),-log(3),-log(2),-log(1.5),-log(1),-log(0.7)]) #y轴标刻度
    ax.set_yticklabels([20000,10000,5000,3000,2000,1500,1000,700]) #y轴刻度标对应分子量
    ax.set_xlim(3,11) #x轴上界下界
    ax.set_ylim(-log(0.7),-log(20)) #y轴上界下界
    ax.legend(loc='lower right') #图例放在右下角
    plt.show()
