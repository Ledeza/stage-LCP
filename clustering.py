#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 16:16:41 2019

@author: stage6
"""
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
array = np.loadtxt("D:/Cours/L3 Chimie/L3 PCM Stage fin annee/dynamique t1 t2/geom_dityr_3c_t1_300K.dat")



#array = array[:,1:]
#plt.scatter(array[:,0],array[:,1],s=0.0001)
#plt.show()
#array = array - np.mean(array,0)
#array2 = array/np.std(array,0)
#cov = array2.T.dot(array2)/(np.shape(array2)[0]-1)
#diagonale, passage = np.linalg.eig(cov)
#print(diagonale)
#array2 = np.dot(array2, passage)
#array = np.dot(array, passage)
#c = ['r','b','purple','g']
#plt.scatter(array[:,0],array[:,1],s=0.0001)
#plt.show()
#plt.scatter(array2[:,0],array2[:,1],s=0.0001)
#plt.show()
#for i in range(4):
#    plt.scatter(array[40000*i:40000*(i+1),0],array[40000*i:40000*(i+1),3], c=c[i],s=0.0001)
#    plt.show()
#plt.scatter(array[:,0],array[:,2],s=0.0001)
#plt.show()


#for i in range(3):
#    ax = plt.axes(projection='3d')
#    ax.scatter(array[:,0],array[:,1],array[:,i+2],s=0.001)
#    ax.show()