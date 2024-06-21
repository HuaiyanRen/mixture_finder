import csv 
import numpy as np
import os
from datetime import datetime, timedelta

with open('result_time.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes',  'ntaxa', 'sites', 'invariable', 'rate', 'tree_length',
                        'optimal1', 'invar1', 'time1', 'mem1', 
                        'optimal2', 'invar2', 'time2', 'mem2',
                        'optimal3', 'invar3', 'time3', 'mem3'])
                       
classes = [5] 
rates = [2] # 0: +E, 1: +I, 2: +I+G
length = [10000]
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

def time_str_to_seconds(time_str):
    parts = time_str.split(':')
    if len(parts) == 3:
        hours = int(parts[0])
        minutes = int(parts[1])
        seconds = int(parts[2])
        total = hours*3600 + minutes*60 + seconds
    elif len(parts) == 2:
        minutes = int(parts[0])
        seconds = float(parts[1])
        total = minutes*60 + seconds
    return total

def get_para(line):
    gtr_all = line.split()[-1]
    gtr_list = gtr_all.split(',')
    
    if '+FO{' in gtr_all:
        A = float(gtr_list[-4].split('{')[-1])
        C = float(gtr_list[-3])
        G = float(gtr_list[-2])
        T = float(gtr_list[-1].split('}')[0])
    else:
        A,C,G,T = 0.25,0.25,0.25,0.25

    if 'JC' in line or 'F81' in line:
        AC,AG,AT,CG,CT = 1,1,1,1,1
    elif 'K2P' in line or 'HKY' in line or 'K80' in line:
        AC,AT,CG = 1,1,1
        AG = float(gtr_list[0].split('{')[1].split('}')[0])
        CT = AG
    elif 'TN' in line:
        AC,AT,CG = 1,1,1
        AG = float(gtr_list[0].split('{')[1])
        CT = float(gtr_list[1].split('}')[0])
    elif 'K3P' in line or 'K81' in line:
        AC = 1
        AG = float(gtr_list[0].split('{')[1])
        AT = float(gtr_list[1].split('}')[0])
        CG = AT
        CT = AG
    elif 'TPM2' in line:
        CG = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1].split('}')[0])
        AT = AC
        CT = AG
    elif 'TPM3' in line:  
        AT = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1].split('}')[0])
        CG = AC
        CT = AG
    elif 'TIM2' in line:
        CG = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        AT = AC
        CT = float(gtr_list[2].split('}')[0])
    elif 'TIM3' in line:
        AT = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        CG = AC
        CT = float(gtr_list[2].split('}')[0])
    elif 'TIM' in line:
        AC = 1
        AG = float(gtr_list[0].split('{')[1])
        AT = float(gtr_list[1])
        CG = AT
        CT = float(gtr_list[2].split('}')[0])
    elif 'TVM' in line:
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        AT = float(gtr_list[2])
        CG = float(gtr_list[3].split('}')[0])
        CT = AG
    elif 'SYM' in line or 'GTR' in line:
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        AT = float(gtr_list[2])
        CG = float(gtr_list[3])
        CT = float(gtr_list[4].split('}')[0])
    return AC,AG,AT,CG,CT,A,C,G,T     
    
for paras in tuple_list:
    classes, rates, length, ntaxa, replicates = paras
    
    # dataset_name
    file_name = 'c' + str(classes) + '_r' + str(rates) + '_l' + str(length) + '_t' + str(ntaxa) + '_rep' + str(replicates)
    simu_file = file_name + '.treefile.log'
    invariable = 0
    with open(simu_file) as b:
        for line in b.readlines():
            if 'Proportion of invariable sites:' in line:
                invariable = float(line.split()[-1])
                
    true_treefile = file_name + '.treefile'
    true_tree = open(true_treefile,'r').read()
    tts = true_tree.split(':')
    true_tree_length = 0
    for i in range(1,len(tts)):
        true_tree_length = true_tree_length + float(tts[i].split(',')[0].split(')')[0])
        
    result_row = [file_name, classes, ntaxa, length, invariable, 2, true_tree_length]
    
    if os.path.isfile('allto/'+ file_name + '.iqtree'):
        invar = 0 
        with open('allto/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'Best-fit model according to BIC:' in line:
                    optimal = line.count(',') + 1
                    model_re = line.split()[-1]
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
        
        #time and mem
        t_list = []
        m_list = []
        if os.path.isfile('allto/'+ file_name + '_time.txt'):
            with open('allto/'+ file_name + '_time.txt') as b:
                for line in b.readlines():
                    if 'Elapsed (wall clock) time (h:mm:ss or m:ss):' in line:
                        t = line.split()[-1]
                        t_list.append(t)
                    if 'Maximum resident set size (kbytes):' in line:
                        m = float(line.split()[-1])/1024
                        m_list.append(m)
            
        if len(t_list) == 1:
            timem = time_str_to_seconds(t_list[0]) + 48*3600
        else:
            timem = time_str_to_seconds(t_list[0])
        memm = max(m_list)
                
            
        result_row = result_row + [optimal, str(invar), timem, memm]
    else:
        result_row = result_row + ['','','','']     

        
    if os.path.isfile('gtrto/'+ file_name + '.iqtree'):
        invar = 0 
        with open('gtrto/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'Best-fit model according to BIC:' in line:
                    optimal = line.count(',') + 1
                    model_re = line.split()[-1]
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
       
        #time and mem
        tg_list = []
        mg_list = []
        if os.path.isfile('gtrto/'+ file_name + '_time.txt'):
            with open('gtrto/'+ file_name + '_time.txt') as b:
                for line in b.readlines():
                    if 'Elapsed (wall clock) time (h:mm:ss or m:ss):' in line:
                        tg = line.split()[-1]
                        tg_list.append(tg)
                    if 'Maximum resident set size (kbytes):' in line:
                        mg = float(line.split()[-1])/1024
                        mg_list.append(mg)
        

        timeg = time_str_to_seconds(tg_list[0])
        memg = max(mg_list)
                    
        result_row = result_row + [optimal, str(invar), timeg, memg]
    else:
        #print(paras, 'method' + str(i))
        result_row = result_row + ['','','',''] 
        
        
    if os.path.isfile('onet/'+ file_name + '.iqtree'):
        invar = 0 
        with open('onet/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
       
        #time and mem
        tg_list = []
        mg_list = []
        if os.path.isfile('onet/'+ file_name + '_time.txt'):
            with open('onet/'+ file_name + '_time.txt') as b:
                for line in b.readlines():
                    if 'Elapsed (wall clock) time (h:mm:ss or m:ss):' in line:
                        tg = line.split()[-1]
                        tg_list.append(tg)
                    if 'Maximum resident set size (kbytes):' in line:
                        mg = float(line.split()[-1])/1024
                        mg_list.append(mg)
        

        timeg = time_str_to_seconds(tg_list[0])
        memg = max(mg_list)
                    
        result_row = result_row + [1, str(invar), timeg, memg]
    else:
        #print(paras, 'method' + str(i))
        result_row = result_row + ['','','','']     

    
    with open('result_time.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(result_row)         
    