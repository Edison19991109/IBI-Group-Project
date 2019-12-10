#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 18:58:57 2019

@author: cuijiajun
"""
import re
def Cap(s):
    c = ''
    for char in s:
        if re.match('[A-Za-z]',char):
            c += char.upper()
    return c
dna_seq = Cap(input('Input a sequence of DNA (5\'-3\'):\n'))
complement = {'A':'T','G':'C','C':'G','T':'A'}
com_seq = ''
if re.search(r'[^ATGC]',dna_seq):
    print('\nInvalid sequence')
else:
    for i in dna_seq:
        com_seq += complement[i] 
    com_seq = com_seq[::-1]
    print ('\nComplementary DNA strand (5\'-3\') is:\n'+com_seq)