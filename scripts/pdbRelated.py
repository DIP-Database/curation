#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 21:23:53 2021

@author: ericwolos
"""
import re
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio import pairwise2
import urllib
from lxml import etree
import os.path
import requests
from Bio import SeqIO

class subunit:
    def __init__(self,identifier,stoichiometry,uniprotID,seq):
        self.identifier = identifier
        self.stoichiometry = stoichiometry
        self.uniprotID = uniprotID
        self.seq = seq


class pdb:
    def __init__(self, identifier, subunits):
        self.identifier = identifier
        self.subunits = subunits
        
def scrapeComposition(composition,chainDictionary):
    listofSU = []
    
    p = re.compile('([^\[\]]+)\[(\d+)\]')
    parsed = p.findall(composition)
    
    for element in parsed:
        identifier = element[0]
        uniprotID = chainDictionary[identifier][0]
        seq = chainDictionary[identifier][1]
        listofSU.append(subunit(identifier,element[1],uniprotID,seq))
    
    return listofSU

def printPDBInfo(newPDB):
    print(newPDB.identifier)
    for subunit in range(0,len(newPDB.subunits)):
        curSubUnit = newPDB.subunits[subunit]
        return (curSubUnit.identifier + " (" + curSubUnit.uniprotID + ")" 
                + " occurs " + str(curSubUnit.stoichiometry) + " times " +
                "with sequence: " + curSubUnit.seq)
    
def getIsoform(uniprotID, dirPath):
    url = "https://www.uniprot.org/uniprot/" + uniprotID + ".fasta"
    filePath = dirPath + uniprotID + ".fasta"
    if (os.path.exists(filePath)):
        file = open(filePath,"r")
    else:
        result = requests.get(url)
        with open(filePath, 'wb') as file:
            file.write(result.content)
        file = open(filePath,"r")
    return file

def getIsoformNames(uniprotID,dirPath):
    url = "https://www.uniprot.org/uniprot/" + uniprotID + ".xml"
    file = getPageXml(url, dirPath + uniprotID + ".xml")
    parsedurl = etree.parse(file)
    isoforms = parsedurl.xpath('/u:uniprot/u:entry/u:comment[@type="alternative products"]/u:isoform/u:id/text()',namespaces = {"u":"http://uniprot.org/uniprot"})
    return isoforms

def scrapeCIF(cifFile, dirPath):
    
    mmcifDictionary = MMCIF2Dict(cifFile)
    mmcifUniprotID = mmcifDictionary["_struct_ref_seq.pdbx_db_accession"]
    mmcifChainID = mmcifDictionary["_struct_ref_seq.pdbx_strand_id"]
    
    mmcifSeqUniprotID = mmcifDictionary["_struct_ref.pdbx_db_accession"]
    mmcifChainSeq = mmcifDictionary["_struct_ref.pdbx_seq_one_letter_code"]
    
    chainDict = {}
    for proteinNum in range(0,len(mmcifUniprotID)):
        
        curChain = mmcifChainID[proteinNum]
        seq = ""
        curUniprotID = mmcifUniprotID[proteinNum]
    
        
        # Get the FASTA sequence from cif by searching with the Uniprot ID
        
        for proteinNum in range(0,len(mmcifChainSeq)):
            if mmcifSeqUniprotID[proteinNum] == curUniprotID:
                seq = re.sub(r'[\W_]+', '', mmcifChainSeq[proteinNum])
            
                
            
        isfmNames = getIsoformNames(curUniprotID, dirPath)
        
        #Generate a dictionary of format:
        # isfmMasterDict[isoform uniprot ID] = sequence of isoform
        isfmMasterDict = {}
        likelyID = curUniprotID
        if (len(isfmNames) > 0):
            for isfm in isfmNames:
                isfmFile = getIsoform(isfm, dirPath)
                isfmDict = SeqIO.to_dict(SeqIO.parse(isfmFile, "fasta"))
                for record in isfmDict:
                    isfmMasterDict[record] = isfmDict[record].seq
                
                
            # rank this sequence against possible isoforms on uniprot
            likelyID = rankAlignments(seq, isfmMasterDict, curUniprotID)
        
        chainDict[curChain] = [likelyID,seq]
    
    
    return chainDict


def getPageXml(url, filePath):
    
    if (os.path.exists(filePath)):
        file = open(filePath,"r")
    else:
        parsedurl = etree.parse(urllib.request.urlopen(url))
        result = etree.tostring(parsedurl, method="html")
        file = open(filePath, 'w')
        file.write(result.decode("utf-8"))
        file = open(filePath,"r")
    return file


def getCif(pdb, filePath):
    url = "https://files.rcsb.org/download/" + pdb + ".cif"
    if (os.path.exists(filePath)):
        file = open(filePath,"r")
    else:
        result = requests.get(url)
        with open(filePath, 'wb') as file:
            file.write(result.content)
        file = open(filePath,"r")
    return file
    




def rankAlignments(sequence,listOfPossibleIsoforms, uniprotID):
    maxScore = -1
    maxScoreName = ""
    for isfmName in listOfPossibleIsoforms:
        isoform = listOfPossibleIsoforms[isfmName]#.seq
        alignmentScore = pairwise2.align.globalms(sequence, isoform, 2, -1, -.5, -.5,score_only=True)
        if (alignmentScore > maxScore):
            maxScore = alignmentScore
            maxScoreName = isfmName
    if (maxScoreName != ""):
        print ("found likely isoform for ", uniprotID, " with score ", maxScore)
        p = re.compile('\w+\-\d')
        parsed = p.findall(maxScoreName)
        if (len(parsed) > 0):
            return parsed[0]
        else:
            return maxScoreName
    else:
        print ("unable to find likely isoform for ", uniprotID, "with score ", maxScore)
        return uniprotID
    
    
    
    
    
    
    

        
        
