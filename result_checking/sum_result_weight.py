import numpy as np
import csv 
from scipy.optimize import linear_sum_assignment

with open('result_weight.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes',  'ntaxa', 'sites', 'invariable',
                        't_ali',
                        'e_ali', 'rmse', 'rmse_matrix'])
                       
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


def file_name_row(file_name):
    with open('modelpara.csv', mode='r', newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            if row and row[0] == file_name:
                return row
    return None 

def normal_q(r,f):
    q= [[0,        r[0]*f[1],r[1]*f[2],r[2]*f[3]],
        [r[0]*f[0],0        ,r[3]*f[2],r[4]*f[3]],
        [r[1]*f[0],r[3]*f[1],0        ,1*f[3]],
        [r[2]*f[0],r[4]*f[1],1*f[2]   ,0]]

    q[0][0] = -sum(q[0])
    q[1][1] = -sum(q[1])
    q[2][2] = -sum(q[2])
    q[3][3] = -sum(q[3])

    dia_sum = - (f[0]*q[0][0] + f[1]*q[1][1] + f[2]*q[2][2] + f[3]*q[3][3])
    for i in range(4):
        for j in range(4):
            q[i][j] = q[i][j]/dia_sum
        
    return q

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
    
    result_row = [file_name, classes, ntaxa, length, invariable]            
    
    if int(classes) == 1:
        tw_list = [1]
        ew_list = [1]   
        result_row = result_row + ['', '', '', '']
        
    else:
        line = file_name_row(file_name)
        tq_list = []
        tw_list = []
        for i in range(0, int(classes)):
            s = [float(line[i*11+10]), float(line[i*11+11]), float(line[i*11+12]), float(line[i*11+13]), float(line[i*11+14])]
            f = [float(line[i*11+15]), float(line[i*11+16]), float(line[i*11+17]), float(line[i*11+18])]
            q = normal_q(s, f)
            tq_list.append([q[0][1], q[0][2], q[0][3], q[1][2], q[1][3], q[2][3]])
            tw_list.append(float(line[i*11+9]))
            
            
        with open('cor/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'odel of substitution:' in line:
                    est_q_list = line.split('{')[-1].split('}')[0].split(',')
                    for q in range(int(classes)):
                        if '+F' in est_q_list[q]:
                            est_q_list[q] = est_q_list[q].replace("+FO", "")
            
        eq_list = []                    
        ew_list = []
        with open('cor/'+ file_name + '.iqtree') as b:
            for line in b.readlines():
                for i in range(1,int(classes)+1):
                    if str(i) + '  '+ est_q_list[i-1] in line:
                        ew_list.append(float(line.split()[-2]))
                        q = get_para(line)[0:5]
                        f = get_para(line)[5:9]
                        q = normal_q(s, f)
                        eq_list.append([q[0][1], q[0][2], q[0][3], q[1][2], q[1][3], q[2][3]])
        ew_list = [n/sum(ew_list) for n in ew_list]
    
        rmse_matrix = []
        for i in range(int(classes)):
            row_list = []
            for j in range(int(classes)):
                rmse = np.sqrt(np.mean((np.array(tq_list[i]) - np.array(eq_list[j])) ** 2))
                row_list.append(rmse)
            rmse_matrix.append(row_list)
        
        t_ali, e_ali = linear_sum_assignment(rmse_matrix)
        sort_tw_list = [tw_list[i] for i in t_ali]
        sort_ew_list = [ew_list[i] for i in e_ali]
        w_rmse = np.sqrt(np.mean((np.array(sort_tw_list) - np.array(sort_ew_list)) ** 2))

        result_row = result_row + [t_ali, e_ali, w_rmse, rmse_matrix]
    
    
    with open('result_weight.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(result_row)         
    
