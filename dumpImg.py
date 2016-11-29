#!/usr/bin/python2

import cv2
import sys
import csv

img = cv2.imread(sys.argv[1],0)

with open(sys.argv[2],'wb') as f:
    w = csv.writer(f)
    w.writerows(img)

