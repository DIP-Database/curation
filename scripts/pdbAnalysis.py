#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 19:07:10 2021

@author: ericwolos
"""


import pdbRelated as pb
import argparse
from lxml import etree
import os
import pdbToTab as tb

myparser = argparse.ArgumentParser( description='PDB Analyzer' )


myparser.add_argument( '--pdb', '-s',  dest="pdb", type=str, required=True,
                     help='Four character pdb')

myparser.add_argument( '--pmid', '-p',  dest="pmid", type=str, required=True,
                     help='Eight digit PubMed Identification Code')

myparser.add_argument( '--source', '-src',  dest="source", type=str, required=False,
                      default="DIP",
                     help='Curation Source (default is DIP)')

args = myparser.parse_args()
pdbname = args.pdb.lower()
source = args.source
if args.pmid.isnumeric() and len(args.pmid) == 8:
    pmid = args.pmid
else:
    print("Invalid pmid.")

multimerUrl = "https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?" + pdbname



parentdir = os.getcwd() + "/"
path = parentdir + "cache/" + pdbname[:2] + "/" + pdbname + "/"
try:
    os.makedirs(path, exist_ok = True)
except OSError:
    print("Directory '%s' can not be created" % path)
    

multimerPath = path + pdbname + ".mlt.xml"
multimerFile = pb.getPage(multimerUrl, multimerPath)
stoichDict = pb.scrapeComposition(multimerFile)

interfacePath = path + pdbname + ".xml"
interface = pb.findInteractions(pdbname,interfacePath)
chainlist = interface[0]
interfacePairs = interface[1]



cifPath = path + pdbname + ".cif"
cifFile = pb.getCif(pdbname,cifPath)

chainDict = pb.scrapeCIF(cifFile,path)

listOfSU = []
for chain in chainlist:
    if chain in chainDict:
        uniprotID = chainDict[chain][0]
        if len(uniprotID) != 6:
            print (chain," is likely an antibody")
        else:
            stoich = 0
            if chain in stoichDict:
                stoich = stoichDict[chain]
            seq = chainDict[chain][0]
            newSU = pb.subunit(chain, stoich, chainDict[chain][0], seq)
            listOfSU.append(newSU)


newpdb = pb.pdb(pdbname, listOfSU)

tabfile = tb.TabFile(pmid)

doi = "testdoi"
tabfile.addHeader(doi, source, pmid)
tabfile.addInteraction("interaction.txt")

SUnames = []

for SU in listOfSU:
    SUnames.append(SU.identifier)
    tabfile.addMolecule(SU.uniprotID, "molecule.txt")
    
for pair in interfacePairs:
    if pair[0] in SUnames and pair[1] in SUnames and pair[0] != pair[1]:
        print(pair[0], " interacts with ", pair[1])
        
        tabfile.addInteraction("interaction.txt")
        
        chain1 = chainDict[pair[0]][0]
        chain2 = chainDict[pair[1]][0]
        
        tabfile.addMolecule(chain1, "molecule.txt")
        tabfile.addMolecule(chain2, "molecule.txt")
        


