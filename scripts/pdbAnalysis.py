#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 19:07:10 2021

@author: ericwolos
"""


import pdbRelated as pb
import argparse
import os
import pdbToTab as tb

myparser = argparse.ArgumentParser( description='PDB Analyzer' )


myparser.add_argument( '--pdb', '-s',  dest="pdb", type=str, required=True,
                     help='Four character pdb')

myparser.add_argument( '--source', '-src',  dest="source", type=str, required=False,
                      default="DIP",
                     help='Curation Source (default is DIP)')

args = myparser.parse_args()
pdbname = args.pdb.lower()
source = args.source

multimerUrl = "https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?" + pdbname


parentdir = os.getcwd() + "/"

    
pathPISA = parentdir + "pisa" + "/" + pdbname[:1] + "/" + pdbname[:2] + "/" + pdbname + "/"
pb.makePath(pathPISA)

multimerPath = (pathPISA + pdbname + ".mlt.xml")
multimerFile = pb.getPage(multimerUrl, multimerPath)
stoichDict = pb.scrapeComposition(multimerFile)



interfacePath = (pathPISA + pdbname + ".int.xml")
interface = pb.findInteractions(pdbname,interfacePath)
chainlist = interface[0]
interfacePairs = interface[1]


pathPDB = parentdir + "pdb" + "/" + pdbname[:1] + "/" + pdbname[:2] + "/" + pdbname + "/"
pb.makePath(pathPDB)
cifPath = pathPDB + pdbname + ".cif"

cifFile = pb.getCif(pdbname,cifPath)
# Copy of cif File necessary because parser module changes the object so it cannot be parsed more than once
cifFile2 = pb.getCif(pdbname,cifPath)
cifFile3 = pb.getCif(pdbname,cifPath)


chainDict = pb.scrapeCIF(cifFile,cifFile2,pathPDB)

listOfSU = []
listOfTXID = []
SUnames = []


for chain in chainlist:
    if chain in chainDict:
        uniprotID = chainDict[chain]["likelyUniprotID"]
        
        stoich = 0
        if chain in stoichDict:
            stoich = stoichDict[chain]
        SUtaxid = chainDict[chain]["taxid"]
        SUrefid = chainDict[chain]["refid"]
        newSU = pb.subunit(chain, stoich, uniprotID, chainDict[chain]["sequence"], SUtaxid, SUrefid)
        listOfSU.append(newSU)
        listOfTXID.append(SUtaxid)
        SUnames.append(chain)
        


newpdb = pb.pdb(pdbname, listOfSU)



headerInfo = pb.getHeaderInfo(cifFile3)
pmid = headerInfo[0][0]
doi = headerInfo[1][0]
pathOUTPUT = parentdir + "output" + "/" + pdbname[:1] + "/" + pdbname[:2] + "/" + pdbname + "/"
pb.makePath(pathOUTPUT)
outputPath = (pathOUTPUT + pdbname + "Tab").replace("$REPLACE","output")

tabfile = tb.TabFile(outputPath)
tabfile.addHeader("header.txt",doi, source, pmid, pdbname)

overallTaxid = -1
if (len(set(listOfTXID)) == 1):
    overallTaxid = listOfTXID[0]
    
tabfile.addInteraction("interaction.txt",True,pdbname,overallTaxid)


refIDs = {}
for SU in listOfSU:
    if SU.refid not in refIDs:
        refIDs[SU.refid] = 1
    else:
        refIDs[SU.refid] += 1

for refID in refIDs:
    for SU in listOfSU:
        if refID == SU.refid:

            tabfile.addMolecule(SU.uniprotID, SU.taxid, refIDs[refID], "molecule.txt", chainDict[SU.identifier]["isfmRange"])
            break
    

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
        chain1BR = chain1["isfmRange"]
        
        chain2UP = chain2["likelyUniprotID"]
        chain2NCBI = chain2["taxid"]
        chain2BR = chain2["isfmRange"]
        
        tabfile.addMolecule(chain1UP, chain1NCBI, 0, "molecule.txt",chain1BR)
        tabfile.addMolecule(chain2UP, chain2NCBI, 0, "molecule.txt", chain2BR)
        


