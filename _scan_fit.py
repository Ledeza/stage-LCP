#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 12:27:16 2019

@author: thirion
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

QM_data = np.loadtxt('/userTMP/stage6/deza/fitting/')
Diff_data = np.loadtxt('/userTMP/stage6/deza/fitting/')

def dihed_k(x,k):
    return k*(1+np.cos)

def dihed_kbis(x,k):
    return 2*(1+np.cos(np.radians(1*x-180)))+2*(1+np.cos(np.radians(1*(x-180)-180)))

def dihed_kter(x,k2,k4,k6):
    return 2*k2*(1+np.cos(np.radians(2*x-180)))+2*k6*(1+np.cos(np.radians(6*x-0)))

def dihed_kk(x,k1,k2):
    return 1.5*(1+np.cos(np.radians(3*x-180))) + 1*(1+np.cos(np.radians(1*x-170))) + 1*(1+np.cos(np.radians(1*x-190)))

def dihed_n(x,n):
    return 5.5*(1+np.cos(np.radians(n*x-180)))

def dihed_nk(x,k,n):
    return k*(1+np.cos(np.radians(n*x-180)))

def dihed_deltak(x,k1,delta1,k2,delta2):
    return 1.5*(1+np.cos(np.radians(3*x-delta1)))+k2*(1+np.cos(np.radians(1*x-delta2)))+k2*(1+np.cos(np.radians(1*x-(360-delta2))))

def dihed_full(x,k1,k2,k3,n1,n2,n3,delta1,delta2,delta3):
    return k1*(1+np.cos(np.radians(n1*x-delta1)))+k2*(1+np.cos(np.radians(n2*x-delta2)))+k3*(1+np.cos(np.radians(n3*x-delta3)))
def dihed_3cos_3par(x,k1,k2,k3,delta1,delta2,delta3,C):
    return k1*(1+np.cos(np.radians(1*x-delta1)))+k2*(1+np.cos(np.radians(2*x-delta2)))+k3*(1+np.cos(np.radians(3*x-delta3)))+C

def dihed_4cos_4par(x,k1,k2,k3,k4,delta1,delta2,delta3,delta4,C):
    return k1*(1+np.cos(np.radians(1*x-delta1)))+k2*(1+np.cos(np.radians(2*x-delta2)))+k3*(1+np.cos(np.radians(3*x-delta3)))+ k4*(1+np.cos(np.radians(4*x-delta4)))+C

def dihed_4cos_2par(x,k1,k2,k3,k4,C):
    return k1*(1+np.cos(np.radians(1*x-27)))+k2*(1+np.cos(np.radians(2*x+54)))+k3*(1+np.cos(np.radians(3*x+62)))+ k4*(1+np.cos(np.radians(4*x-88)))+C

def dihed_4cos_par(x,n4,k4,delta):
    return 0.186*(1+np.cos(np.radians(2.474*x-180)))+3.293*(1+np.cos(np.radians(1.998*x-0)))+3.645*(1+np.cos(np.radians(0.123*x-60)))+k4*(1+np.cos(np.radians(n4*x-delta)))

def dihed_3cos_par(x,n3,k3,delta,C):
    return 2.38*(1+np.cos(np.radians(2.0*x+50)))+3.293*(1+np.cos(np.radians(2.0*x-0)))+k3*(1+np.cos(np.radians(n3*x-delta)))+C

def dihed_2cos(x,k1,k2,n1,n2,delta1,delta2,C):
    return k1*(1+np.cos(np.radians(n1*x-delta1)))+k2*(1+np.cos(np.radians(n2*x-delta2)))+C

def dihed_2cos_par(x,k1,k2,delta1,delta2,C):
    return k1*(1+np.cos(np.radians(2*x-delta1)))+k2*(1+np.cos(np.radians(1*x-delta2)))+C

def dihed_2cos_2par(x,k1,k2,C):
    return k1*(1+np.cos(np.radians(2*x+50)))+k2*(1+np.cos(np.radians(1*x-27)))+C

def dihed_1cos(x,k1,C):
    return k1*(1+np.cos(np.radians(2.0*x+50)))+C

X = np.linspace(-180,180,1001)   
poptk, pcov = curve_fit(dihed_4cos_4par,QM_data[0:71],Diff_data[0:71],bounds=([-5,-5,-5,-5,-180,-180,-180,-180,-10],[5,5,5,5,180,180,180,180,10]),method='trf')
Y1 = dihed_4cos_4par(X,*poptk)
moindre_carre = 0
for i in range(len(QM_data)):
    moindre_carre = moindre_carre + (Diff_data[i] - dihed_4cos_4par(QM_data[i],*poptk))**2

poptk, pcov = curve_fit(dihed_4cos_2par,QM_data[0:71],Diff_data[0:71],bounds=([-5,-5,-5,-5,-10],[5,5,5,5,10]),method='trf')
Y1 = dihed_4cos_2par(X,*poptk)
print(poptk,"4cos")
#poptk, pcov = curve_fit(dihed_4cos_par,QM_data[0:71],Diff_data[0:71],bounds=([0,0,0],[5,5,180]),method='trf')
#poptk, pcov = curve_fit(dihed_3cos_par,QM_data[0:71],Diff_data[0:71],bounds=([0,0,0],[5,5,180]),method='trf')
#poptk, pcov = curve_fit(dihed_2cos_par,QM_data[0:71],Diff_data[0:71],bounds=([-5,-5,-180,-180,-10],[5,5,180,180,10]),method='trf')
#poptk, pcov = curve_fit(dihed_2cos,QM_data[0:71],Diff_data[0:71],bounds=([0,0,0,0,-180,-180,-10],[20,20,5,5,180,180,10]),method='trf')
#poptk, pcov = curve_fit(dihed_2cos_2par,QM_data[0:71],Diff_data[0:71],bounds=([0,0,-10],[5,5,10]),method='trf')
#print(poptk,"stop")
#print(pcov,"stop")
#print(np.sqrt(np.diag(pcov)))
#print(*poptk)
#Y = dihed_4cos_4par(X,*poptk)
#Y = dihed_4cos_par(X,*poptk)
#Y = dihed_3cos_par(X,*poptk)
poptk, pcov = curve_fit(dihed_3cos_3par,QM_data[0:71],Diff_data[0:71],bounds=([-5,-5,-5,-180,-180,-180,-10],[5,5,5,180,180,180,10]),method='trf')
Y = dihed_3cos_3par(X,*poptk)
moindre_carre = 0
for i in range(len(QM_data)):
    moindre_carre = moindre_carre + (Diff_data[i] - dihed_3cos_3par(QM_data[i],*poptk))**2
#print(moindre_carre,"3cos")
#Y = dihed_2cos(X,*poptk)
poptk, pcov = curve_fit(dihed_2cos_par,QM_data[0:71],Diff_data[0:71],bounds=([-5,-5,-180,-180,-10],[5,5,180,180,10]),method='trf')
Y3 = dihed_2cos_par(X,*poptk)
moindre_carre = 0
for i in range(len(QM_data)):
    moindre_carre = moindre_carre + (Diff_data[i] - dihed_2cos_par(QM_data[i],*poptk))**2
#print(moindre_carre,"2cos")
#print(poptk)
#Y = dihed_2cos_par(X,*poptk)
#Y = dihed_2cos_2par(X,*poptk)
poptk, pcov = curve_fit(dihed_1cos,QM_data[0:71],Diff_data[0:71],bounds=([0,-1],[5,5]),method='trf')
Y2 = dihed_1cos(X,*poptk)
moindre_carre = 0
for i in range(len(QM_data)):
    moindre_carre = moindre_carre + (Diff_data[i] - dihed_1cos(QM_data[i],*poptk))**2
#print(moindre_carre,"1cos")
#print(pcov,"stop")
plt.plot(QM_data[0:71],Diff_data[0:71],'ro',)
#plt.plot(QM_data[2:71],Diff_data_nosym[2:71],'yo',)
#plt.plot(X,dihed_full(X,*poptkk),'b-')
#plt.plot(X,dihed_kbis(X,*poptkk),'g-')
#plt.plot(X,Y,'-',color="black")
plt.plot(X,Y1,'-',color="blue")
#plt.plot(X,Y2,'-',color="green")
#plt.plot(X,Y3,'-',color="purple")
#print(pcov)