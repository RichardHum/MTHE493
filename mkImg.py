#!/usr/bin/python2

import cv2
import sys
import csv
import numpy as np

with open(sys.argv[1],'r') as f:
    r = csv.reader(f,delimiter=',')
    img = []
    for row in r:
        img.append(row)

for row in img:
    del row[len(row)-1]

img = np.array(img)

cv2.imwrite(sys.argv[2],img.astype(int))
