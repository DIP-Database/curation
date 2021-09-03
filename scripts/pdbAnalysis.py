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

interfacePath = path + pdbname + ".int.xml"
interface = pb.findInteractions(pdbname,interfacePath)
chainlist = interface[0]
interfacePairs = interface[1]



cifPath = path + pdbname + ".cif"
cifFile = pb.getCif(pdbname,cifPath)
# Copy of cif File necessary because parser module changes the object so it cannot be parsed twice
cifFile2 = pb.getCif(pdbname,cifPath)
cifFile3 = pb.getCif(pdbname,cifPath)


chainDict = pb.scrapeCIF(cifFile,cifFile2,path)


listOfSU = []
listOfTXID = []
for chain in chainlist:
    if chain in chainDict:
        uniprotID = chainDict[chain]["likelyUniprotID"]
        
        stoich = 0
        if chain in stoichDict:
            stoich = stoichDict[chain]
        SUtaxid = chainDict[chain]["taxid"]    
        newSU = pb.subunit(chain, stoich, uniprotID, chainDict[chain]["sequence"], SUtaxid)
        listOfSU.append(newSU)
        listOfTXID.append(SUtaxid)
        


newpdb = pb.pdb(pdbname, listOfSU)

tabfile = tb.TabFile(pmid)

headerInfo = pb.getHeaderInfo(cifFile3)
pmid = headerInfo[0][0]
doi = headerInfo[1][0]
tabfile.addHeader("header.txt",doi, source, pmid, pdbname)

overallTaxid = -1
if (len(set(listOfTXID)) == 1):
    overallTaxid = listOfTXID[0]
    
tabfile.addInteraction("interaction.txt",True,pdbname,overallTaxid)

SUnames = []

for SU in listOfSU:
    SUnames.append(SU.identifier)
    tabfile.addMolecule(SU.uniprotID, SU.taxid, "molecule.txt")
    
for pair in interfacePairs:
    if pair[0] in SUnames and pair[1] in SUnames and pair[0] != pair[1]:
        
        chain1 = chainDict[pair[0]]
        chain2 = chainDict[pair[1]]
        
        chain1TXID = chain1["taxid"]
        chain2TXID = chain2["taxid"]
        
        # -1 represent in vitro
        taxid = "-1"
        if (chain1TXID == chain2TXID):
            taxid = chain1TXID
        
        tabfile.addInteraction("interactionInterfaces.txt",False,pdbname,taxid)
        
        chain1UP = chain1["likelyUniprotID"]
        chain1NCBI = chain1["taxid"]
        chain2UP = chain2["likelyUniprotID"]
        chain2NCBI = chain2["taxid"]
        
        tabfile.addMolecule(chain1UP, chain1NCBI, "molecule.txt")
        tabfile.addMolecule(chain2UP, chain1NCBI, "molecule.txt")
        


