#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 18:45:37 2019

@author: cuijiajun
"""
# a function to calculate the radio of GC according to the input
import re
def Cap(s):
    c = ''
    for char in s:
        if re.match('[A-Za-z]',char):
            c += char.upper()
    return c
seq = Cap(input('Input a sequence of DNA:\n'))
if re.search(r'[^ATCG]',seq):
    print('\nInvalid sequence')
else:
    GC_count = float(seq.count('C') + seq.count('G'))
    GC_ratio = '%.2f%%' % (GC_count / len(seq) * 100)
    print('\nGC content:',GC_ratio)