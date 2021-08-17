#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 21:23:53 2021

@author: ericwolos
"""
import re
from Bio import pairwise2
import urllib
from lxml import etree
import os.path
import requests
from Bio import SeqIO
from Bio.Align import substitution_matrices


# For debugging
import pprint
pp = pprint.PrettyPrinter(indent=4) 

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
    file = getPage(url, dirPath + uniprotID + ".xml")
    parsedurl = etree.parse(file)
    isoforms = parsedurl.xpath('/u:uniprot/u:entry/u:comment[@type="alternative products"]/u:isoform/u:id/text()',namespaces = {"u":"http://uniprot.org/uniprot"})
    return isoforms

def scrapeCIF(cifFile, dirPath):

    
    chainDict = {}
    for record in SeqIO.parse(cifFile,"cif-seqres"):
        if (len(record.dbxrefs) > 0):
            
            curUniprotID = record.dbxrefs[0][4:]
            curChain = re.sub('\w+\:', '', record.id)
            seq = record.seq
            
            isfmNames = getIsoformNames(curUniprotID, dirPath)
            
            #Generate a dictionary of format:
            # isfmMasterDict[isoform uniprot ID] = sequence of isoform
            
            isfmMasterDict = {}
            likelyID = (curUniprotID,False)
            if (len(isfmNames) > 0):
                for isfm in isfmNames:
                    isfmFile = getIsoform(isfm, dirPath)
                    isfmDict = SeqIO.to_dict(SeqIO.parse(isfmFile, "fasta"))
                    for record in isfmDict:
                        isfmMasterDict[isfm] = isfmDict[record].seq
                    
                    
                # rank this sequence against possible isoforms on uniprot
                likelyID = rankAlignments(seq, isfmMasterDict, curUniprotID)
            
            likelyIsfmName = likelyID[0]
            isfmRange = likelyID[1]
            chainDict[curChain] = [likelyIsfmName,seq, isfmRange]     
        
    
    return chainDict


def getPage(url, filePath):
    
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
    blosum = substitution_matrices.load("BLOSUM62")
    maxScore = -1
    maxScoreName = ""
    scores = []
    isBR = False
    for isfmName in listOfPossibleIsoforms:
        isoform = listOfPossibleIsoforms[isfmName] #isoform sequence
        alignmentScore = pairwise2.align.globalms(sequence, isoform, 2, -1 , -2, -0.1, score_only=True)
        scores.append(alignmentScore)
        if (alignmentScore > maxScore):
            maxScore = alignmentScore
            maxScoreName = isfmName
            
    
    # Verify that any likely isoform was identified
    if (maxScoreName == ""):
        print ("unable to find likely isoform for ", uniprotID)
        return (uniprotID, isBR)
    
    # Check if mutiple isoforms get the same score
    duplicates = False
    if (scores.count(maxScore) > 1):
        duplicates = True
        
    # Detect range by finding first identical residue in both first and reversed sequences
    bestAlignment = pairwise2.align.globaldx(sequence, listOfPossibleIsoforms[maxScoreName] , blosum)
    seqA = bestAlignment[0][0]
    seqB = bestAlignment[0][1]
    beginIndex = findFirstIdenticalResidue(seqA, seqB)
    # reverses sequences
    reverseSeqA = seqA[::-1]
    reverseSeqB = seqB[::-1]
    endIndex = len(seqA) - findFirstIdenticalResidue(reverseSeqA, reverseSeqB) + 1

    print(seqA)
    print(seqB)

    dashes = seqA[beginIndex-1:endIndex]
    print(dashes)
    dashNum = dashes.count('-')
    
    alignRange = (beginIndex - dashNum, endIndex - dashNum)
        
    
    if (duplicates == False):
        print ("found likely isoform for ", uniprotID, " with score ", maxScore, " and range", alignRange)
        p = re.compile('\w+\-\d')
        parsed = p.findall(maxScoreName)
        if (len(parsed) > 0):
            return (parsed[0],alignRange)
        else:
            return maxScoreName
    else:
        print ("multiple likely isoforms for ", uniprotID, "with score ", maxScore, " and range", alignRange)
        return (uniprotID, alignRange)

def findFirstIdenticalResidue(seqA,seqB):
    for residueNum in range(0,len(seqA)):
        if seqA[residueNum] == seqB[residueNum]:
            beginIndex = residueNum+1
            return beginIndex
    
def findInteractions(pdb, filePath):
    interfaces = "https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?" + pdb.lower()
    interfaceFile = getPage(interfaces, filePath)
    parsedurl = etree.parse(interfaceFile)
    chain1 = parsedurl.xpath('/pisa_interfaces/pdb_entry/interface/h-bonds/bond/chain-1/text()')
    chain2 = parsedurl.xpath('/pisa_interfaces/pdb_entry/interface/h-bonds/bond/chain-2/text()')
    
    numchains = len(set(chain1+chain2))
    
    interactorPairs = []
    for index in range(0,len(chain1)):
        if chain1[index] < chain2[index]:
            interactorPairs.append((chain1[index],chain2[index]))
        else:
            interactorPairs.append((chain2[index],chain1[index]))
    result = []
    [result.append(x) for x in interactorPairs if x not in result]

    # Use of a graph data structure to store the interactions
    # Each chain is a key to a dictionary whose value is a list of connected chains
    graph = {}
    for pair in result:
        if pair[0] != pair[1]:
            if not pair[0] in graph:
                graph[pair[0]] = list(pair[1])
            elif graph[pair[0]].count(pair[1]) < 1:
                graph[pair[0]].append(pair[1])
                
    isDirect = True
    for key in graph:
        if len(graph[key]) != numchains - 1:
            isDirect = False
            break
        
    if isDirect:
        print ("direct association")
    else:
        print ("physical association")
        
    pp.pprint(graph)
    
    return graph
    

    
    
    
    
    
    

        
        
