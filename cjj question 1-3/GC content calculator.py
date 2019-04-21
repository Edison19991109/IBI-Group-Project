#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 18:45:37 2019

@author: cuijiajun
"""
# a function to calculate the radio of GC according to the input
def seq_GCRatio(sequence):
    GC_count = float(sequence.count('C') + sequence.count('G'))
    seq_length = len(sequence)
    GC_ratio = str(GC_count / seq_length * 100) + '%'
    return GC_ratio
seq = input()
print (seq_GCRatio(seq))