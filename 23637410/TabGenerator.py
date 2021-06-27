#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:11:55 2021

@author: ericwolos
"""

import uniprotToPaper as up
import argparse

myparser = argparse.ArgumentParser( description='Tab delimited txt file Constructor\nEnsure that the proper molecule.txt and interaction.txt files are included' )

myparser.add_argument( '--pmid', '-p',  dest="pmid", type=str, required=True,
                     help='Eight digit PubMed Identification Code')

myparser.add_argument( '--source', '-s',  dest="source", type=str, required=False,
                      default="DIP",
                     help='Curation Source (default is DIP)')


args = myparser.parse_args()

pmid = "1"
source = "DIP"
if args.pmid.isnumeric() and len(args.pmid) == 8:
    pmid = args.pmid
else:
    print("Invalid pmid.")


newPaper = open(pmid+".txt","w")

up.addHeader(newPaper,args.source,pmid)

up.addInteraction(newPaper,"interaction.txt")


uniprotFile = "Genetouniprot.txt"  
UNIIDlist = up.extractUNIID(uniprotFile)


up.writeSimilarMolecules(UNIIDlist,newPaper)

newPaper.close()