#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 10:37:18 2019

@author: deza
"""
#DEF##########################################################################
def reduction(chain):
    position =[]
    for i in range(len(chain)):
        if chain[i] != " ":
            position.append(i)
    chain2 = ""
    for i in range(len(position)):
        if position[i]+1 in position:
            chain2 += chain[position[i]]
        else:
            chain2 += chain[position[i]]+" "
    return chain2
def split_int(chain):
    L = chain.split(" ")
    L = L[:-1]
    for i in range(len(L)):
        if L[i][-1].isdigit():
            L[i]=float(L[i])
    return L
def dihed_4cos_4par(x,k1,k2,k3,k4,delta1,delta2,delta3,delta4,C):
    return k1*(1+np.cos(np.radians(1*x-delta1)))+k2*(1+np.cos(np.radians(2*x-delta2)))+k3*(1+np.cos(np.radians(3*x-delta3)))+ k4*(1+np.cos(np.radians(4*x-delta4)))+C
def moindre_carre(L1,L2):
    if len(L1) == len(L2):
        mc = 0
        for i in range(len(L1)):
            mc += (L1[i]-L2[i])**2
        return mc
    else:
        return 9999999
#START########################################################################
#mydir = '/userTMP/stage6/deza/fitting/'
mydir = 'D:/Cours/L3 Chimie/L3 PCM Stage fin annee/resultats qm mm/'
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import time
startTime = time.time()
ldir = os.listdir(mydir)

#Puting apart .dat files and recovering names
ltrie,lqm = [],[]
for i in range(0,len(ldir)):
    if 'dat' in ldir[i] and 'delta' not in ldir[i]:
        ltrie.append(ldir[i])
        if '_qm' in ldir[i]:
            lqm.append(ldir[i])
ldir = ltrie

#Check
print("You have "+str(int(len(ldir)/2))+" conditions:")
for i in range(len(lqm)):
    print("- "+lqm[i][:-7])

#Putting list together
lcond = []
for i in range(len(lqm)):
    L = [lqm[i]]
    for j in range(len(ldir)):
        if (lqm[i][:-7]+"_mm") in ldir[j]:
            L.append(ldir[j])
    lcond.append(L)
print("\nPairing Done")

#table
megatable = []
for i in range(0,len(lcond)):
    a,b,table = [],[],[]
    for j in range(0,len(lcond[i])):
        if j == 0: #Meaning qm work
            for line in open(mydir+lcond[i][j]):
                a.append(split_int(reduction(line)))
            a = a[:-1]
            table.append([elem[0] for elem in a])
            table.append([elem[1] for elem in a])
        else:
            for line in open(mydir+lcond[i][j]):
                line = reduction(line)
                b.append(split_int(reduction(line)))
            table.append([elem[2] for elem in b])

    megatable.append(table)
print("Megatable Done")

#Calculs
for i in range(0,len(megatable)):
    megatable[i].append([])
    for j in range(0,len(megatable[i][1])):
        megatable[i][-1].append(megatable[i][1][j]*627.509-megatable[i][2][j])
    megatable[i].append([])
    for j in range(0,len(megatable[i][1])):
        megatable[i][-1].append(megatable[i][-2][j]-min(megatable[i][-2]))
print("Calculation Done")

#Plot
lfigname = ["Energie \u03B81 avec 2 angles calibrés","Energie \u03B81 avec \u03B82 calibré","Energie QM \u03B81 avec \u03B82 fixé et calibré","Energie QM \u03B81 avec \u03B82 fixé","Energie QM \u03B81","Energie QM \u03B82 calibré avec \u03B81 fixé","Energie QM \u03B82 avec \u03B81 fixé","Energie QM \u03B82"]
col = ['blue','red','purple','grey','cyan','green','orange','black','pink''yellow']
for i in range(0,len(megatable)):
    plt.plot(megatable[i][0],megatable[i][-1], 'ro', markersize=3, color=col[i])
    plt.title(lfigname[i])
    plt.xlabel("Angle")
    plt.ylabel("Energie (kcal/mol)")
    plt.show()    
for i in range(0,len(megatable)):
    plt.plot(megatable[i][0],megatable[i][-1], 'ro', markersize=3, color=col[i])
plt.title("Energie QM")
plt.xlabel("Angle")
plt.ylabel("Energie (kcal/mol)")
plt.show()
print("Plot Qm done")

#fitting
lfigname = ["\u0394Energie \u03B81 avec 2 angles calibrés","\u0394Energie \u03B81 avec \u03B82 calibré","\u0394Energie \u03B81 avec \u03B82 fixé et calibré","\u0394Energie \u03B81 avec \u03B82 fixé","\u0394Energie \u03B81","\u0394Energie \u03B82 calibré avec \u03B81 fixé","\u0394Energie \u03B82 avec \u03B81 fixé","\u0394Energie \u03B82"]
X = np.linspace(-180,180,1001)
for i in range (0,len(lqm)):
    poptk, pcov = curve_fit(dihed_4cos_4par,megatable[i][0],megatable[i][-1],bounds=([0,0,0,0,-180,-180,-180,-180,-10],[5,5,5,5,180,180,180,180,10]),method='trf')
    Y = dihed_4cos_4par(X,*poptk)
    mc = moindre_carre(megatable[i][-1],dihed_4cos_4par(np.asarray(megatable[i][0]),*poptk))
    somme = sum(megatable[i][-1])
    print("\n--------------------------------------------------")
    print(lqm[i][:-7])
    plt.plot(X,Y)
    plt.plot(megatable[i][0],megatable[i][-1],'ro', markersize=4)
    plt.title(lfigname[i])
    plt.xlabel("Angle")
    plt.ylabel("Energie (kcal/mol)")
    plt.savefig(lfigname[i])
    plt.show()
    print('Parameters are: ')
    for j in range(0,int((len(poptk)-1)/2)):
        print('   -k'+str(j+1)+': '+str(poptk[j])+' delta'+str(j+1)+': '+str(poptk[j+4]))
    print("Shift is :",poptk[-1])
    print('Least square represent '+str(mc**0.5/somme*100)+'% of the total sum of point')

#Saving
for i in range(0,len(lqm)):
    newtxt = open(mydir+"delta_"+lqm[i][:-7]+".dat",'w')
    for j in range(0,len(megatable[i][0])):
        newtxt.write("%f"%megatable[i][0][j]+" "+"%f"%megatable[i][-1][j]+"\n")
    newtxt.close()
print('\nWork Done in: '+str(time.time() - startTime)+" s")

#Figure prod
##Fig plt QM for Theta 1 and 2
#lfigname = ["Energie \u03B81 avec 2 angles calibrés","Energie \u03B81 avec \u03B82 calibré","Energie QM \u03B81 avec \u03B82 fixé et calibré","Energie QM \u03B81 avec \u03B82 fixé","Energie QM \u03B81","Energie QM \u03B82 calibré avec \u03B81 fixé","Energie QM \u03B82 avec \u03B81 fixé","Energie QM \u03B82"]
#plt.subplot(211)
#plt.plot(megatable[4][0],megatable[4][1], 'ro', markersize=3, color="green")
#plt.title(lfigname[4])
#plt.xlabel("Angle")
#plt.ylabel("Energie (Hartree)")
#plt.subplot(212)
#plt.plot(megatable[7][0],megatable[7][1], 'ro', markersize=3, color="blue")
#plt.title(lfigname[7])
#plt.xlabel("Angle")
#plt.ylabel("Energie (Hartree)")
#plt.tight_layout()
##Fig plt QM for theta 1 and 2 fixed
#lfigname = ["Energie \u03B81 avec 2 angles calibrés","Enerxgie \u03B81 avec \u03B82 calibré","Energie QM \u03B81 avec \u03B82 fixé et calibré","Energie QM \u03B81 avec \u03B82 fixé","Energie QM \u03B81","Energie QM \u03B82 calibré avec \u03B81 fixé","Energie QM \u03B82 avec \u03B81 fixé","Energie QM \u03B82"]
#plt.subplot(211)
#plt.plot(megatable[3][0],megatable[3][1], 'ro', markersize=3, color="green")
#plt.title(lfigname[3])
#plt.xlabel("Angle")
#plt.ylabel("Energie (Hartree)")
#plt.subplot(212)
#plt.plot(megatable[6][0],megatable[6][1], 'ro', markersize=3, color="blue")
#plt.title(lfigname[6])
#plt.xlabel("Angle")
#plt.ylabel("Energie (Hartree)")
#plt.tight_layout()
#plt.savefig("biplotQM_fixed.png")
#Showing energy difference

