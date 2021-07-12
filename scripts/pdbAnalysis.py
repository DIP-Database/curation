#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 19:07:10 2021

@author: ericwolos
"""

import urllib
import pdbRelated as pb
import argparse
from lxml import etree
import os.path
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
    
# Multimer
multimerFile = path + pdbname + ".html"

if (os.path.exists(multimerFile)):
    file = open(multimerFile,"r")
else:
    parsedurl = etree.parse(urllib.request.urlopen(multimerUrl))
    result = etree.tostring(parsedurl, method="html")
    file = open(multimerFile, 'w')
    file.write(result.decode("utf-8"))
    file = open(multimerFile,"r")
    

parsedurl = etree.parse(file)
composition = parsedurl.xpath('/pisa_multimers/pdb_entry/asm_set/assembly/composition/text()')[0]
# End multimer

cifFile = pdbname + ".cif"

chainDict = pb.scrapeCIF(cifFile)

newPDB = pb.pdb(pdbname,pb.scrapeComposition(composition,chainDict))
print(pb.printPDBInfo(newPDB))


