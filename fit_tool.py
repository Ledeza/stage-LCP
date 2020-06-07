#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:18:24 2019

@author: deza
"""
###################################################################################
def dihed_4cos_2par(x,k1,k2,k3,k4,delta1,delta2,delta3,delta4,C):
    return k1*(1+np.cos(np.radians(1*x-delta1)))+k2*(1+np.cos(np.radians(2*x-delta2)))+k3*(1+np.cos(np.radians(3*x-delta3)))+ k4*(1+np.cos(np.radians(4*x-delta4)))+C
def dihed_3cos_2par(x,k1,k2,k3,delta1,delta2,delta3,C):
    return k1*(1+np.cos(np.radians(1*x-delta1)))+k2*(1+np.cos(np.radians(2*x-delta2)))+k3*(1+np.cos(np.radians(3*x-delta3)))+C
def dihed_2cos_2par(x,k1,k2,delta1,delta2,C):
    return k1*(1+np.cos(np.radians(2*x-delta1)))+k2*(1+np.cos(np.radians(4*x-delta2)))+C
def dihed_1cos_2par(x,k1,delta1,C):
    return k1*(1+np.cos(np.radians(2*x+delta1)))+C
def moindre_carre(L1,L2):
    if len(L1) == len(L2):
        mc = 0
        for i in range(len(L1)):
            mc += (L1[i]-L2[i])**2
        return mc
    else: 
        return 9999999
def optimise(funct1,array,n):
    lbd, hbd =[], [] #lower bounds, higher bounds
    for i in range(n):
        lbd.append(0)
        hbd.append(5)
    for i in range(n):
        lbd.append(-180)
        hbd.append(180)
    lbd.append(0)
    hbd.append(10)
    X = np.linspace(-180,180,1001)
    poptk, pcov = curve_fit(funct1,array[:,0],array[:,1],bounds=(lbd,hbd),method='trf')
    Y1 = funct1(X,*poptk)
    plt.plot(X,Y1)
    plt.plot(array[:,0],array[:,1],'ro')
    plt.show()
    print('Parameters are: ')
    for i in range(0,int((len(poptk)-1)/2)):
        print('   -k'+str(i+1)+': '+str(poptk[i])+' delta'+str(i+1)+': '+str(poptk[i+n]))
    print(poptk[-1])
    mc = moindre_carre(array[:,1],funct1(array[:,0],*poptk))
    somme = sum(array[:,1])
    print('Least square represent '+str(mc**0.5/somme*100)+'% of the total sum of point')
    for i in range(n):
        delta = input("Fix the value of delta"+str(i+1)+" please: ")
        lbd[n+i],hbd[n+i] = delta-0.0000001,delta
    poptk, pcov = curve_fit(funct1,array[:,0],array[:,1],bounds=(lbd,hbd),method='trf')
    Y2 = funct1(X,*poptk)
    plt.plot(X,Y2)
    plt.plot(array[:,0],array[:,1],'ro')
    plt.show()
    print('Parameters are: ')
    for i in range(0,(len(poptk)-1)/2):
        print('   -k'+str(i+1)+': '+str(poptk[i])+' delta'+str(i+1)+': '+str(poptk[i+n]))
    mc = moindre_carre(array[:,1],funct1(array[:,0],*poptk))
    somme = sum(array[:,1])
    print('Least square represent '+str(mc**0.5/somme*100)+'% of the total sum of point')
###################################################################################
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#Diff_data = np.loadtxt('/userTMP/stage6/deza/fitting/delta_o1_cal.dat')
Diff_data = np.loadtxt("D:/Cours/L3 Chimie/L3 PCM Stage fin annee/resultats qm mm/delta_o2_o1fixed_cal.dat")
flag =0
while flag == 0:
    ncos = int(input("How many cosinus are you keeping form plot_fit: "))
    if ncos == 1:
        optimise(dihed_1cos_2par,Diff_data,ncos)
        flag = 1
    elif ncos == 2:
        optimise(dihed_2cos_2par,Diff_data,ncos)
        flag = 1
    elif ncos == 3:
        optimise(dihed_3cos_2par,Diff_data,ncos)
        flag = 1
    elif ncos == 4:
        optimise(dihed_4cos_2par,Diff_data,ncos)
        flag = 1
    else:
        print("Number of cosinus is invalid, please make sure ncos is in [1,4] range")
print("\nDone")