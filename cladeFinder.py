# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:23:13 2018

@author: hunte
"""

"""
python cladeFinder.py YFull_YTree_v6.02__20180402.tree.json pos negatives output

create negatives
bcftools filter -O z -i '(GT=="1/1" && AA==ALT) || (GT=="0/0" && AA=REF)' chrY_cleaned_1_hg38.vcf | bcftools query -f '%ID,' > negatives

create positives
bcftools query -f '%ID,' chrY_derived_1_hg38.vcf.gz > pos

"""


import json
import sys


if len(sys.argv) > 1:
    treeFile = sys.argv[1]
    tabixFilePath = sys.argv[2]
    outputFile = sys.argv[3]
    maxThreads = int(sys.argv[4])

import time
start = time.time()

import tabix

#TODO get unique column values tabix query?

tb = tabix.open(tabixFilePath)
idresults = tb.querys("id:1-9999999")
ids = []
for theid in idresults:
    ids.append(theid[1])


def getPositivesAndNegatives(sampleId, tb):
    posResults = tb.querys("pos:" + sampleId + "-" + sampleId)
    positives = []
    for p in posResults:
        positives.append(p[3])
    positives = set(positives)
    negResults = tb.querys("neg:" + sampleId + "-" + sampleId)
    negatives = []
    for n in negResults:
        if len(n) < 4:
            print('negative results too few rows',n)
        else:
            negatives.append(n[3])
    negatives = set(negatives)
    return list(positives) + list(negatives)


allStart = time.time()



import requests

def findClade(snps):
    url = "https://cladefinder.yseq.net/json.php"
    PARAMS = {"input": ", ".join(snps),"json":""}
    r = requests.post(url=url,data=PARAMS)
    response = json.loads(r.text)
    if "clade" in response:
        return response["clade"]
    else:
        return "?"

cladeMap = {}
for theid in ids:
    cladeMap[theid] = findClade(getPositivesAndNegatives(theid, tb))

allEnd = time.time()
print('clade finder executed on ' + str(len(ids)) + ' ' + str(round((allEnd-allStart)/60,2)) + " minutes\n")
print(str(round((allEnd-allStart)/len(ids),2)) + " seconds per sample\n")


def writeSampleCladesFile(cladeMap):
    w = open(outputFile,"w+")
    for sample in cladeMap:
        w.write(",".join([sample,cladeMap[sample]]) + "\n")


writeSampleCladesFile(cladeMap)

