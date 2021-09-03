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
    def addMolecule(self,uniprotID, taxid, templateFile):
        self.tabfile.write('\n')
        molecule = ""
        template = open(templateFile,'r')
        for ln in template:
            molecule += ln
        cmolecule = molecule.replace("%UNIPROT%","uprot:" + uniprotID).replace("%TAXID%",taxid)
        self.tabfile.write(cmolecule)
        
        
    # Add interaction using a template file
    def addInteraction(self,templateFile,isPhysical, pdb, taxid):
        self.tabfile.write('\n')
        interaction = ""
        
        interactionType = "MI:0407(direct)"
        if (isPhysical):
            interactionType = "MI:0915(physical)"
        
        template = open(templateFile,'r')
        for ln in template:
            interaction += ln
        cinteraction = interaction.replace("%INTERACTION%",interactionType).replace("%PDB%",pdb).replace("%TAXID%",taxid)
        self.tabfile.write(cinteraction)
        
    def addHeader(self,templateFile,doi,source,pmid,pdb):
        header = ""
    
        template = open(templateFile,'r')
        for ln in template:
            header += ln
        cheader = header.replace("%DOI%",doi).replace("%SOURCE%",source).replace("%PMID%",pmid).replace("%PDB%",pdb)
        self.tabfile.write(cheader)


    
