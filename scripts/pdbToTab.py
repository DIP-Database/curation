#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 16:12:32 2021

@author: ericwolos
"""


class TabFile:
    
    # Initialize by writing/opening a file
    def __init__(self,filename):
        self.tabfile = open(filename+".txt","w")
        
    # Add molecule to interaction using a template file    
    def addMolecule(self,uniprotID, templateFile):
        self.tabfile.write('\n')
        molecule = ""
        template = open(templateFile,'r')
        for ln in template:
            molecule += ln
        cmolecule = molecule.replace("%UNIPROT%","uprot:" + uniprotID)
        self.tabfile.write(cmolecule)
        
        
    # Add interaction using a template file
    def addInteraction(self,templateFile):
        self.tabfile.write('\n')
        interaction = ""
        template = open(templateFile,'r')
        for ln in template:
            interaction += ln
        self.tabfile.write(interaction)
        
    def addHeader(self,doi,source,pmid):
        header = "source\t" + source + '\n'+"doi\t" + doi + "\n" + "pmid\t" + pmid + "\n"
        self.tabfile.write(header)

    
