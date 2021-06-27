#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 18:49:22 2021

@author: ericwolos
"""

import uniprotToPaper as up
import numpy as np

uniprotFile = "Genetouniprot.txt"  
OriginalIdentifier = np.genfromtxt(uniprotFile,dtype=str,delimiter='\t',usecols=(0))
FinalGeneNamelist = OriginalIdentifier[1:]

initialGeneNames = np.loadtxt("geneNames.txt",dtype=str) 

up.compareLists(initialGeneNames,FinalGeneNamelist)

