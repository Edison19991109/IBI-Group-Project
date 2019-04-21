#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 18:58:57 2019

@author: cuijiajun
"""
import re
dna_seq = input('give me a sequence of DNA (please use capital letters): ')
if re.search('a|t|c|g',dna_seq) :
    print ('please use capital letters')
elif re.search('A|T|G|C',dna_seq):
    complement = {'A':'T','G':'C','C':'G','T':'A'}
    com_seq = ''
    dna_seq = dna_seq.strip()
    for i in dna_seq:
        com_seq += complement[i]
    com_seq = com_seq[::-1]
    print ('complementary DNA strand (from 5 - 3) is:',com_seq)