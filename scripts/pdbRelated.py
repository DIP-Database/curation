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
from Bio.PDB import MMCIF2Dict as mm
from Bio.Align import substitution_matrices



# For debugging
import pprint
pp = pprint.PrettyPrinter(indent=4) 

class subunit:
    def __init__(self,identifier,stoichiometry,uniprotID,seq,taxid):
        self.identifier = identifier
        self.stoichiometry = stoichiometry
        self.uniprotID = uniprotID
        self.seq = seq
        self.taxid = taxid


class pdb:
    def __init__(self, identifier, subunits):
        self.identifier = identifier
        self.subunits = subunits
        
    def printPDBInfo(self):
        print(self.identifier)
        for subunit in range(0,len(self.subunits)):
            curSubUnit = self.subunits[subunit]
            print (curSubUnit.identifier + " (" + curSubUnit.uniprotID + ")" 
                   + " occurs " + str(curSubUnit.stoichiometry) + " times " +
                   "with sequence: " + curSubUnit.seq)
        

# Analyzes composition of pdb via a PISA multimers file
def scrapeComposition(multimerFile):
    stoichDict = {}
    multparsedurl = etree.parse(multimerFile)
    
    # Verify that a valid entry is present
    if (multparsedurl.xpath('/pisa_multimers/pdb_entry/status/text()')[0] != "Ok"):
        return stoichDict
    
    composition = multparsedurl.xpath('/pisa_multimers/pdb_entry/asm_set/assembly/composition/text()')[0]
    
    # Use regular expressions to find the stoichiometry of each chain
    p = re.compile('([^\[\]]+)\[(\d+)\]')
    parsed = p.findall(composition)
    
    # Fill up dictionary with chain and stoichiometry
    for element in parsed:
        identifier = element[0]
        stoichDict[identifier] = element[1]
        
    return stoichDict


# Gets fasta files from uniprot
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


# Scrapes isoform names based off of a uniprot ID
def getIsoformNames(uniprotID,dirPath):
    url = "https://www.uniprot.org/uniprot/" + uniprotID + ".xml"
    file = getPage(url, dirPath + uniprotID + ".xml")
    parsedurl = etree.parse(file)
    isoforms = parsedurl.xpath('/u:uniprot/u:entry/u:comment[@type="alternative products"]/u:isoform/u:id/text()',namespaces = {"u":"http://uniprot.org/uniprot"})
    return isoforms



def scrapeCIF(cifFile, cifFileCopy,dirPath):
    
    chainDict = {}
    
    parsedCif = SeqIO.parse(cifFile,"cif-seqres")

    for record in parsedCif:
        if (len(record.dbxrefs) > 0):
            print(record.description)
            curUniprotID = record.dbxrefs[0][4:]
            curChain = re.sub('\w+\:', '', record.id)
            seq = record.seq
            
            # Unsed to seperate out invalid uniprotIDs
            if (len(curUniprotID) != 6 and len(curUniprotID) != 10):
                continue
            
            
            isfmNames = getIsoformNames(curUniprotID, dirPath)
            
            # Generate a dictionary of format:
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
            chainInfo = {}
            chainInfo["sequence"] = seq
            chainInfo["likelyUniprotID"] = likelyIsfmName
            chainInfo["isfmRange"] = isfmRange
            chainDict[curChain] = chainInfo     
        
    # Now, to determine which taxonomic ID corresponds to each subunit. This feature is not supported by Bio.SeqIO so lower level access with MMCIF2Dict is necessary
    mmcifDict = mm.MMCIF2Dict(cifFileCopy)

    ncbiDict = dict(zip(mmcifDict["_entity_src_gen.entity_id"],mmcifDict["_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id"]))
    refDict = dict(zip(mmcifDict["_struct_ref_seq.pdbx_strand_id"],mmcifDict["_struct_ref_seq.ref_id"]))
    
    for key in chainDict:
        chainDict[key]["taxid"] = ncbiDict[refDict[key]]

    pp.pprint(chainDict)
    return chainDict

def getHeaderInfo(cifFile):
    mmcifDict = mm.MMCIF2Dict(cifFile)
    pmid = mmcifDict["_citation.pdbx_database_id_PubMed"]
    doi = mmcifDict["_citation.pdbx_database_id_DOI"]
    return (pmid,doi)

# The following two functions are designed to check if the file exists locally, and if not, get them from the internet
def getPage(url, filePath):
    
    if (os.path.exists(filePath)):
        file = open(filePath,"r")
    else:
        print(url)
        print(filePath)
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
    # First, align:
        
    bestAlignment = pairwise2.align.globalms(sequence, listOfPossibleIsoforms[maxScoreName] , 2,-1,-2,-0.1)
    
    # Both alignment sequences:
    seqA = bestAlignment[0][0]
    seqB = bestAlignment[0][1]
    beginIndex = findFirstIdenticalResidue(seqA, seqB)
    
    # reverse sequences to find the first identical residue at the end of the range:
    reverseSeqA = seqA[::-1]
    reverseSeqB = seqB[::-1]
    endIndex = len(seqA) - findFirstIdenticalResidue(reverseSeqA, reverseSeqB) + 1

    # interpret with the alignment output to find the range:
    dashes = seqA[beginIndex-1:endIndex]

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
    
    # conversion of chain lists to one single list of pairs of interacting residues
    interactorPairs = []
    for index in range(0,len(chain1)):
        interactorPairs.append((chain1[index],chain2[index]))
        
    # removal of duplicates
    result = []
    [result.append(x) for x in interactorPairs if x not in result]
    
    #pp.pprint(result)
    
    chainlist = []
    for pair in result:
        if pair[0] not in chainlist:
            chainlist.append(pair[0])
        if pair[1] not in chainlist:
            chainlist.append(pair[1])
    
    
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
        
    return (chainlist,result)  
  
    
    # The graph structure is useful for visualizing what subunits interact with each other
    #pp.pprint(graph)
    
    #return graph
    

    
    
    
    
    
    

        
        
