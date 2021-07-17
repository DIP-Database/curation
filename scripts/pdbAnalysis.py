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


myparser = argparse.ArgumentParser( description='PDB Analyzer' )


myparser.add_argument( '--pdb', '-s',  dest="pdb", type=str, required=True,
                     help='Four character pdb')


args = myparser.parse_args()
pdbname = args.pdb.lower()
multimerUrl = "https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?" + pdbname



parentdir = os.getcwd() + "/"
path = parentdir + "cache/" + pdbname[0] + "/" + pdbname + "/"
try:
    os.makedirs(path, exist_ok = True)
except OSError:
    print("Directory '%s' can not be created" % path)
    

multimerPath = path + pdbname + ".html"

multimerFile = pb.getPageXml(multimerUrl, multimerPath)


multparsedurl = etree.parse(multimerFile)
composition = multparsedurl.xpath('/pisa_multimers/pdb_entry/asm_set/assembly/composition/text()')[0]


cifPath = path + pdbname + ".cif"
cifFile = pb.getCif(pdbname,cifPath)

chainDict = pb.scrapeCIF(cifFile,path)

newPDB = pb.pdb(pdbname,pb.scrapeComposition(composition,chainDict))
print(pb.printPDBInfo(newPDB))


