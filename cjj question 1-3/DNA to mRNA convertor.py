#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 19:08:39 2019

@author: cuijiajun
"""
import re
dna_seq = input('give me a sequence of DNA (please use capital letters): ')
if re.search('a|t|c|g',dna_seq) :
    print ('please use capital letters')
elif re.search('A|T|G|C',dna_seq):
    transcription = {'A':'U','G':'C','C':'G','T':'A'}
    rna_seq = ''
    dna_seq = dna_seq.strip()
    for i in dna_seq:
        rna_seq += transcription[i]
    print ('the mRNA sequence is:',rna_seq)
    