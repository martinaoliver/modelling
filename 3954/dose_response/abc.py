# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:10:40 2019

@author: Marti
"""
import random 
import statistics
from abc_functions import gillespie
from abc_functions import celldivision_gillespie


time_storage={}
mRNA_storage={}
mRNAdivision_storage = {}

def simdiff(mRNA):
    mean=statistics.mean(mRNA)#[400:len(mRNA)])
    variance = statistics.variance(mRNA)#[400:len(mRNA)])
    simdiff_value = (10-mean)**2+(10-variance)**2
    return simdiff_value  

def abc():
    for n in range (1,1000):
        i1=random.uniform(0,1)
        i2=random.uniform(0,1)
        if i1>i2:
            time, mRNA= gillespie(0.001,1000)
            simdiff_value1 = simdiff(mRNA)
            return (time,mRNA, simdiff_value1)

        if i2>i1:
            time, mRNA= celldivision_gillespie(0.01,10000)
#            simdiff_value2 = simdiff(mRNA)
#            print (simdiff_value2)
        

x,y,z=abc()


  