# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:20:59 2019

@author: Deza
"""
##############################################################################################################
def number_detection(string,L):
    if len(string) > L[0]:
        newstr= ""
        for i  in range(len(string)-L[1]):
            if string[L[0]+i].isdigit() or string[L[0]+i] == ".":
                newstr = newstr + string[L[0]+i]
        return float(newstr)
    else:
        return 0
def mean(L):
    s = 0
    for i in range(len(L)):
        s += L[i]
    return s/len(L)
##############################################################################################################
mydir = "D:/Cours/L3 Chimie/L3 PCM Stage fin annee/ccs/"
import os
import matplotlib.pyplot as plt
#import numpy as np
import time
startTime = time.time()
ldir = os.listdir(mydir)

#Filtering usefull file
lfil = []
for i in ldir:
    if ".out" in i:
        lfil.append(i)

#Big table containing n array
tabletxt = []
titre = []
for i in range(len(lfil)):
    newtxt = open(mydir+lfil[i],'r')
    tabletxt.append(newtxt.readlines())
    newtxt.close()
    titre.append("Dimère "+lfil[i][11:-6]+" cluster "+lfil[i][14:-4])
print("Opening done")


#table of values
color=["limegreen","cyan","orangered","darkorchid","firebrick","navy","dimgray"]
table = []
for i in range(len(tabletxt)):
    L = []
    for j in range(len(tabletxt[i])):
        L.append(number_detection(tabletxt[i][j],[39,46]))
    table.append(L)
#table.append(table[0]+table[1]+table[2])
#titre.append("t2 cluster 1&2")

#plotting
for i in range(len(table)):
    plt.hist(table[i], 20, facecolor = color[i],edgecolor = 'black')
    plt.show()
    print(mean(table[i]))
print('\nWork Done in: '+str(time.time() - startTime)+" s")