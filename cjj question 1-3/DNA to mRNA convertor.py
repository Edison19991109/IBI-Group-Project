#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 19:08:39 2019

@author: cuijiajun
"""
import re
def Cap(s):
    c = ''
    for char in s:
        if re.match('[A-Za-z]',char):
            c += char.upper()
    return c
dna_seq = Cap(input('give me a sequence of DNA template strand (5\'-3\'):\n'))
if re.search('[^ATCG]',dna_seq) :
    print ('\nInvalid sequence')
else:
    transcription = {'A':'U','G':'C','C':'G','T':'A'}
    rna_seq = ''
    for i in dna_seq:
        rna_seq += transcription[i]
    rna_seq = rna_seq[::-1]
    print ('\nThe hnRNA (mRNA precursor) sequence (5\'-3\') is:\n'+rna_seq)
    