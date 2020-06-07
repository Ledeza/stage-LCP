# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 14:53:49 2019

@author: Deza
"""
############################################################################################
def nb_count(L):
    n = 0
    for i in range(len(L)):
        if L[i]>499:
            n += 1
    return n
############################################################################################
mydir = "D:/Cours/L3 Chimie/L3 PCM Stage fin annee/ccs/"
import os
import matplotlib.pyplot as plt
import numpy as np
import time
startTime = time.time()
ldir = os.listdir(mydir)

#Filtering usefull file
lfil = []
for i in ldir:
    if ".xvg" in i:
        lfil.append(i)

#Big table containing n array
table = []
titre = []
for i in range(len(lfil)):
    table.append(np.loadtxt(mydir+lfil[i]),)
    titre.append(lfil[i][11:-4])

#Explode and labels:
lexp = []
lname = []
for i in range(len(table)):
    N = nb_count(list(table[i][:,1]))
    L = []
    L2 =[]
    for j in range(len(table[i][:,1])):
        if j < N:
            L.append(0.03*(N-j))
            L2.append("Cluster "+str(j+1))
        else:
            L.append(0)
            L2.append("")
    lexp.append(L)
    lname.append(L2)

#Pie
for i in range(len(table)):
    fig1, ax1 = plt.subplots()
    ax1.pie(table[i][:,1], explode=lexp[i], labels=lname[i],shadow=True, startangle=90)
    ax1.axis('equal')
    plt.title(titre[i])
    plt.tight_layout()
    plt.show()
    #fig1.savefig(mydir+titre[i]+".png")
print('\nWork Done in: '+str(time.time() - startTime)+" s")
#, autopct='%1.1f%%'