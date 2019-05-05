# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:07:31 2019

@author: Valerya
"""

import re
def Cap(s):
    c = ''
    for char in s:
        if re.match('[A-Za-z]',char):
            c += char.upper()
    return c
RNA = Cap(input('Input a sequence of mRNA (5\'-3\'):\n'))
if re.search(r'[^AUCG]',RNA):
    print('\nInvalid sequence')
else:
    initiator = re.search(r'AUG',RNA)
    if initiator:
        a = initiator.start()
        RNA = RNA[a:]
        codons = []
        match = {"UUU":"F","UUC":"F","UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L","AUU":"I","AUC":"I","AUA":"I","AUG":"M","GUU":"V","GUC":"V","GUA":"V","GUG":"V","UCU":"S","UCC":"S","UCA":"S","UCG":"S","CCU":"P","CCC":"P","CCA":"P","CCG":"P","ACU":"T","ACC":"T","ACA":"T","ACG":"T","GCU":"A","GCC":"A","GCA":"A","GCG":"A","UAU":"Y","UAC":"Y","CAU":"H","CAC":"H","CAA":"Q","CAG":"Q","AAU":"N","AAC":"N","AAA":"K","AAG":"K","GAU":"D","GAC":"D","GAA":"E","GAG":"E","UGU":"C","UGC":"C","UGG":"W","CGU":"R","CGC":"R","CGA":"R","CGG":"R","AGU":"S","AGC":"S","AGA":"R","AGG":"R","GGU":"G","GGC":"G","GGA":"G","GGG":"G"}
        while len(RNA)>2:
            codons.append(RNA[:3])
            RNA = RNA[3:]
        protein = ""
        if "UAA" in codons or "UAG" in codons or "UGA" in codons:
            for i in range(len(codons)-1):
                if codons[i] == "UAA" or codons[i] == "UAG" or codons[i] == "UGA":
                    break
                else:
                    protein += match[codons[i]]
            print('\nPrimary translation product:\nH-'+protein+'-OH')
        else:
            print("\nInvalid sequence: No terminator \nPlease check if it's a sequence of mature mRNA")
    else:
        print("\nInvalid sequence: No initiator \nPlease check if it's a sequence of mature mRNA")