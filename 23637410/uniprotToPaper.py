#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:07:21 2021

@author: ericwolos
"""

import numpy as np

def addLine(list):
    toBeAdded = ""
    for n in list:
        toBeAdded += n+"\t"
    toBeAdded += "\n"
    return toBeAdded

def extractUNIID(uniprotFile):
    geneToUniprot = np.genfromtxt(uniprotFile,dtype=str,delimiter='\t',usecols=(0,1,2,3,4))
    for colNum in geneToUniprot[0,:]:
        if (geneToUniprot[0,1] == "Entry"):
            return (geneToUniprot[1:,1])
    print("Error, no column titled Entry found")    
    
def writeSimilarMolecules(UNIIDlist,fileToWrite):
    
    molecule = ""
    with open("molecule.txt",'r') as fh:
        for ln in fh:
            molecule += ln
    for UNIID in UNIIDlist:        
        cmolecule = molecule.replace("%UNIPROT%","uprot: " + UNIID)
        fileToWrite.write(cmolecule)
    
def addInteraction(fileToWrite,Interactionfile):
    interaction = ""
    with open(Interactionfile,'r') as fh:
        for ln in fh:
            interaction += ln
    fileToWrite.write(interaction)
    
def addHeader(fileToWrite,pmid,source):
    header = "source\t" + source + '\n'+"pmid\t" + pmid + "\n"
    fileToWrite.write(header)
    
def compareLists(largerList,smallerList):
    dictionary = {}
    for element in largerList:
        dictionary[element] = 0
    for element in smallerList:
        if element in dictionary:
            dictionary[element] += 1
        else:
            print(element," Not found")
            
    for key in dictionary:
        if dictionary[key] > 1:
            print(key," occurred ",dictionary[key], "times")
        elif dictionary[key] == 0:
            print(key," was never found")
