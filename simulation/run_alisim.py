import csv 
import os
from ete3 import Tree
import numpy as np
from scipy.stats import powerlognorm,beta,invweibull,alpha,genlogistic,exponweib,exponnorm,norm, gengamma,dweibull,laplace,loglaplace,loggamma,genexpon,exponpow,maxwell,powernorm,dgamma,chi,weibull_min,invgamma,lognorm,gennorm,invgauss,weibull_max
#from scipy import stats
from multiprocessing import Pool
from functools import partial
import random
import pandas as pd
import ast

tuple_list = []

with open('modelpara.csv') as b:
    for line in b.readlines():
        paras = line.split(',')
        tuple_list.append(paras)
        
tuple_list.pop(0)
for i in range(len(tuple_list)):
    tuple_list[i].append(i)

def run_alisim(line):
    file_name, classes, rates, length, ntaxa = line[0:5]
    classes = int(classes)
    rates = int(rates)
    length = int(length)
    ntaxa = int(ntaxa)
    seed_num = line[-1]
    
    # gtr parameters
    if classes == 1:
        gtr_para = 'GTR{'+line[10]+'/'+line[11] +'/'+line[12]+'/'+line[13]+'/'+line[14]+'}+F{'+line[15]+'/'+line[16]+'/'+line[17]+'/'+line[18]+'}'
        model_para = ' -m ' + gtr_para
    else:        
    # class weights
        model_para = ' -m "MIX{'
        for i in range(classes):
            if i == 0:
                gtr_para = 'GTR{'+line[10]+'/'+line[11] +'/'+line[12]+'/'+line[13]+'/'+line[14]+'}+F{'+line[15]+'/'+line[16]+'/'+line[17]+'/'+line[18]+'}:1:' + line[9]         
                model_para = model_para + gtr_para
            else:
                gtr_para = 'GTR{'+line[i*11+10]+'/'+line[i*11+11] +'/'+line[i*11+12]+'/'+line[i*11+13]+'/'+line[i*11+14]+'}+F{'+line[i*11+15]+'/'+line[i*11+16]+'/'+line[i*11+17]+'/'+line[i*11+18]+'}:1:' + line[i*11+9]          
                model_para = model_para + ',' + gtr_para
        model_para = model_para + '}'   
    
    # rate parameters
    if rates == 0:
        rate_para = ''
    elif rates == 2:
        rate_para = '+G{' + line[7] + '}'
        
    if classes == 1:
        model_para = model_para +rate_para
    else:
        model_para = model_para + rate_para + '"'
        
    # alisim
    alisim_cmd = 'iqtree2 --alisim ' + file_name + model_para + ' --length ' + str(length) + ' -seed ' + str(seed_num) + ' -af fasta -t RANDOM{yh/' + str(ntaxa) + '} -redo'
    os.system(alisim_cmd)
    
partial_running = partial(run_alisim)
with Pool(20) as p:
    p.map(partial_running, tuple_list)