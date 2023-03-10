# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:47:11 2019

@author: Marti
"""

import random
import statistics

def gillespie (timestep_scale, final_time):
    ko = 0.2#random.uniform(0,0.2)
    k1 = 0.01#random.uniform(0,0.02)
    t=0
    m=0
    time_list = [0]
    mRNA_list = [0]
    while t < final_time:
        t1=random.randint(0,final_time+1)*timestep_scale
        t2=random.randint(0,final_time+1)*timestep_scale
        
        if t1>t2:
            ao = ko
            m+=ao
            t+=t1
        if t1<t2:
            if m>0:
                a1 = k1*m
                m-=a1
                t+= t2
        mRNA_list.append(m)
        time_list.append(t)
    return (time_list, mRNA_list)
    

def celldivision_gillespie (timestep_scale, final_time):
    ko = random.uniform(0,0.2)
    k1 = random.uniform(0,0.02)
    t=0
    m=0
    time_list = [0]
    mRNA_list = [0]
    division = 0
    while t < final_time:
        t1=random.randint(0,final_time+1)*timestep_scale
        t2=random.randint(0,final_time+1)*timestep_scale
        
        if t1>t2:
            ao = ko
            m+=ao
            t+=t1
        if t1<t2:
            if m>0:
                a1 = k1*m
                m-=a1
                t+= t2
        mRNA_list.append(m)
        time_list.append(t)
        if t >= 1200*division:
            m/=2
            division+=1
    return (time_list, mRNA_list)

def fano(dictionary, dataset, model):
    simulation_number = 5
    time_range = 10000
    if model == celldivision_gillespie:
        time_storage['t_%s' % n], mRNA_storage['m_%s' % n], mRNAdivision_storage['md_%s' % n] = celldivision_gillespie(0.001,time_range)    
    if model == gillespie:
        for n in range (1,simulation_number+1):
            time_storage['t_%s' % n], mRNA_storage['m_%s' % n] = gillespie(0.001,time_range)    
            
    mean_l = []
    variance_l = []
    for n in range (1,simulation_number+1):
        mean = statistics.mean(dictionary[dataset + '_%s' %n])
        variance = statistics.variance(dictionary[dataset + '_%s' %n])
        mean_l.append(mean)
        variance_l.append(variance)
    mean_simulations = statistics.mean(mean_l)
    variance_simulations = statistics.mean(variance_l)
    return (mean_simulations, variance_simulations)
