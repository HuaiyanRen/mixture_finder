from multiprocessing import Pool
from functools import partial
import subprocess
import numpy as np
import csv 


with open('result_pythia.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes',  'ntaxa', 'sites', 'difficulty'])
                       
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
    #run
    cmd = 'pythia --msa ' + file_name + '.fa --raxmlng /mnt/data/dayhoff/home/u7151703/software/raxml-ng'
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    with open(file_name + '_pythia.txt', 'w+') as f:
        f.write(result.stdout)
        f.write(result.stderr)


partial_running = partial(running_tuple)
with Pool(10) as p:
    p.map(partial_running, tuple_list)


for paras in tuple_list:
    classes, rates, length, ntaxa, replicates = paras
    
    # dataset_name
    file_name = 'c' + str(classes) + '_r' + str(rates) + '_l' + str(length) + '_t' + str(ntaxa) + '_rep' + str(replicates)
    pythia_file = file_name + '_pythia.txt'
    
    with open(pythia_file) as b:
        for line in b.readlines():
            if 'The predicted difficulty for MSA' in line:
                diff = float(line.split()[-1])
    
    result_row = [file_name, classes, ntaxa, length, diff]
    
    with open('result_pythia.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(result_row)   

