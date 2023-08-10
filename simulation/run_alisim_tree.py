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

branch_list = []

with open('branch_length.csv') as b:
    for line in b.readlines():
        paras = line.split(',')
        branch_list.append(paras)
        
tuple_list.pop(0)
branch_list.pop(0)
for i in range(len(tuple_list)):
    tuple_list[i] = tuple_list[i] + branch_list[i] + [i]

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
    #alisim_cmd = 'iqtree2 --alisim ' + file_name + model_para + ' --length ' + str(length) + ' -seed ' + str(seed_num) + ' -af fasta -t RANDOM{yh/' + str(ntaxa) + '} -redo'
    #os.system(alisim_cmd)
    #edit tree length
    treefile0 = file_name + '.treefile'
    tree0 = open(treefile0,'r').read()
    t=Tree(tree0)
    ext_n = 0
    int_n = 0
    ext_branches = line[65:65+ntaxa]
    int_branches = line[66+ntaxa:66+ntaxa+ntaxa-3]
    for n in t.traverse():
        if not n.is_root():
            if n.is_leaf():
                n.dist = ext_branches[ext_n]
                ext_n += 1
            else:
                n.dist = int_branches[int_n]
                int_n += 1
    tree1 = t.write(format=1)
    treefile1 = file_name + '.new.treefile'
    with open(treefile1, 'w+') as result:
        result.write(tree1 + '\n')
    # alisim again with new branch lengths
    alisim_cmd2 = 'iqtree2 --alisim ' + file_name + model_para + ' --length ' + str(length) + ' -seed ' + str(seed_num) + ' -af fasta -t ' + treefile1 + ' -redo'
    os.system(alisim_cmd2) 
    
partial_running = partial(run_alisim)
with Pool(20) as p:
    p.map(partial_running, tuple_list)