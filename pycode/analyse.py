import ipdb
import math
import pdb
import re
# import Image
import os

import numpy as np
import vigra

import basicOperations


## A window centered at GT coordinate
winSize = 231
winSizeR = winSize/2
winSizeR2 = winSize/2+1
cut_orig = True
cut_result = True

## directories

path_result = "/home/seawave/work/output/test1b"


results = basicOperations.getFiles(path_result,'txt')


mitosis = []
for result in results:

    textTemp = open(result,'r')
    lines =  textTemp.readlines()

    if len(lines)==0:
        continue

    i = 0
    for line in lines:
        tempfeat = []
        words = line.split(' ')
        if int(words[-1]) == 0:
            continue
        filename = os.path.split(result)[-1]
        tempfeat.append(filename)
        for x in words:
            tempfeat.append(x)
        mitosis.append(tempfeat)

    textTemp.close()
 
tempfile = open("mitosis_candi.txt", "w")
for x in mitosis:
    for i in range(len(x)-1):
        tempfile.write(x[i] + " ")
    tempfile.write(x[-1])

tempfile.close()


ipdb.set_trace()
