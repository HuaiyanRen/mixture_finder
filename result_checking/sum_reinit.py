import csv 
import numpy as np
import os

with open('result_init.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes',  'ntaxa', 'sites', 'invariable', 
                        'optimala', 'invara', 
                        'reinit_a2', 'JC_a2', 'reinit_a3', 'JC_a3', 'reinit_a4', 'JC_a4',
                        'reinit_a5', 'JC_a5', 'reinit_a6', 'JC_a6', 'reinit_a7', 'JC_a7', 
                        'optimalg', 'invarg',
                        'reinit_g2', 'JC_g2', 'reinit_g3', 'JC_g3', 'reinit_g4', 'JC_g4',
                        'reinit_g5', 'JC_g5', 'reinit_g6', 'JC_g6', 'reinit_g7', 'JC_g7', 
                        ])
                       
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
    simu_file = file_name + '.new.treefile.log'
    invariable = 0
    with open(simu_file) as b:
        for line in b.readlines():
            if 'Proportion of invariable sites:' in line:
                invariable = float(line.split()[-1])
                
    result_row = [file_name, classes, ntaxa, length, invariable]
    
    if os.path.isfile('all/'+ file_name + '.iqtree'):
        invar = 0 
        with open('all/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'odel of substitution:' in line:
                    optimal = line.count(',') + 1
                    model_re = line.split()[-1]
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
    
        reinit_list = [0]*8
        jc_list = ['']*8
        
                
        with open('all/' + file_name + '.log') as b:
            for line in b.readlines():
                if 'with initial weight:' in line:
                    n = line.count(',')
                    reinit_list[n-1] += 1
                if 'est-fit model: M' in line:
                    n = line.count(',')
                    if ',JC}' in line:
                        jc_list[n-1] = True
                    else:
                        jc_list[n-1] = False
            
        result_row =  result_row + [optimal, str(invar),
                                        reinit_list[0], jc_list[0], reinit_list[1], jc_list[1], reinit_list[2], jc_list[2], 
                                        reinit_list[3], jc_list[3], reinit_list[4], jc_list[4], reinit_list[5], jc_list[5]]
    else:
        result_row = result_row + ['','',
                                       '','','','','','',
                                       '','','','','','']
    
    if os.path.isfile('gtr/'+ file_name + '.iqtree'):
        invar = 0 
        with open('gtr/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'odel of substitution:' in line:
                    optimal = line.count(',') + 1
                    model_re = line.split()[-1]
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
    
        reinit_list = [0]*8
        jc_list = ['']*8
        
                
        with open('gtr/' + file_name + '.log') as b:
            for line in b.readlines():
                if 'with initial weight:' in line:
                    n = line.count(',')
                    reinit_list[n-1] += 1
                if 'est-fit model: M' in line:
                    n = line.count(',')
                    if ',JC}' in line:
                        jc_list[n-1] = True
                    else:
                        jc_list[n-1] = False
    
        result_row =  result_row + [optimal, str(invar),
                                    reinit_list[0], jc_list[0], reinit_list[1], jc_list[1], reinit_list[2], jc_list[2], 
                                    reinit_list[3], jc_list[3], reinit_list[4], jc_list[4], reinit_list[5], jc_list[5]]
    else:
        result_row = result_row + ['','',
                                   '','','','','','',
                                   '','','','','','']
    
    with open('result_init.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(result_row)  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    