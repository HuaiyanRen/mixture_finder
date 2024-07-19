import numpy as np
import os
from multiprocessing import Pool
from functools import partial


classes = [1,2,3,4,5] 
rates = [2] # 0: +E, 1: +I, 2: +I+G
length = [1000,2000,5000,10000]
ntaxa = [100]
replicates = list(np.arange(0,20,1))

classes_list = [m for m in classes for i in range(len(rates)*len(length)*len(ntaxa)*len(replicates))]
rates_list = [m for m in rates for i in range(len(length)*len(ntaxa)*len(replicates))]*len(classes)
length_list = [m for m in length for i in range(len(ntaxa)*len(replicates))]*len(classes)*len(rates)
ntaxa_list = [m for m in ntaxa for i in range(len(replicates))]*len(classes)*len(rates)*len(length)
replicates_list = replicates*len(classes)*len(rates)*len(length)*len(ntaxa)
    
tuple_list = [0]*len(classes)*len(rates)*len(length)*len(ntaxa)*len(replicates)

for i in range(len(tuple_list)):
    tuple_list[i] = classes_list[i], rates_list[i], length_list[i], ntaxa_list[i], replicates_list[i]
    

def running_tuple(tuple_list):
    classes, rates, length, ntaxa, replicates = tuple_list
    #file_name
    file_name = 'c' + str(classes) + '_r' + str(rates) + '_l' + str(length) + '_t' + str(ntaxa) + '_rep' + str(replicates)
    tree0 = file_name + '.treefile'

    if os.path.isfile('all/'+ file_name + '.iqtree'): 
        if not os.path.isfile('all/'+ file_name + '_rf.rfdist'): 
            treem = 'all/'+ file_name + '.treefile'
            cmdm = '/scratch/dx61/hr8997/software/iqtree-2.3.5.onnxupdate-Linux-intel/bin/iqtree2 -rf ' + treem + ' ' + tree0 + ' -pre all/' + file_name + '_rf -redo -nt 1'
            os.system(cmdm)
        
    if os.path.isfile('one/'+ file_name + '.iqtree'):
        if not os.path.isfile('one/'+ file_name + '_rf.rfdist'):
            tree1 = 'one/'+ file_name + '.treefile'
            cmd1 = '/scratch/dx61/hr8997/software/iqtree-2.3.5.onnxupdate-Linux-intel/bin/iqtree2 -rf ' + tree1 + ' ' + tree0 + ' -pre one/' + file_name + '_rf -redo -nt 1'
            os.system(cmd1)
        
    if os.path.isfile('gtr/'+ file_name + '.iqtree'):
        if not os.path.isfile('gtr/'+ file_name + '_rf.rfdist'):
            treeg = 'gtr/'+ file_name + '.treefile'
            cmdg = '/scratch/dx61/hr8997/software/iqtree-2.3.5.onnxupdate-Linux-intel/bin/iqtree2 -rf ' + treeg + ' ' + tree0 + ' -pre gtr/' + file_name + '_rf -redo -nt 1'
            os.system(cmdg)
    
    
partial_running = partial(running_tuple)
with Pool(10) as p:
    score_list = p.map(partial_running, tuple_list)