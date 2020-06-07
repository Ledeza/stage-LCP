# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 12:56:36 2019

@author: Deza
"""

mydir = "D:/Cours/L3 Chimie/L3 PCM Stage fin annee/dynamique t1 t2/"
import os
import matplotlib.pyplot as plt
import numpy as np
import time
startTime = time.time()
ldir = os.listdir(mydir)

#Filtering usefull file
lfil = []
for i in ldir:
    if ".dat" in i:
        lfil.append(i)

#Big table containing n array
table = []
titre = []
for i in range(len(lfil)):
    table.append(np.loadtxt(mydir+lfil[i])[:,1:])
    titre.append(lfil[i][5:-4])

#Plotting
col = ["black","green","red","blue"]
for i in range(len(table)):
    plt.subplot(2,2,i+1)    
    for j in range(1):    
        plt.scatter(table[i][j*40000:(j+1)*40000,0],table[i][j*40000:(j+1)*40000,1],s=0.00015,color=col[j])
    plt.title(titre[i])
    plt.xlabel("Angle \u03B81")
    plt.ylabel("Angle \u03B82")
    plt.xlim(-80,80)
plt.tight_layout()
plt.savefig(mydir+"quadriplot.png")
plt.show()
print('\nWork Done in: '+str(time.time() - startTime)+" s")