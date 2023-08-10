import csv 
import numpy as np
from scipy.stats import expon


with open('modelpara.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes', 'rates', 'length', 'taxa', 'replicates',
                        'invar', 'gamma',
                        'q1','weight1','AC1', 'AG1', 'AT1', 'CG1', 'CT1', 'A1', 'C1', 'G1', 'T1',
                        'q2','weight2','AC2', 'AG2', 'AT2', 'CG2', 'CT2', 'A2', 'C2', 'G2', 'T2',
                        'q3','weight3','AC3', 'AG3', 'AT3', 'CG3', 'CT3', 'A3', 'C3', 'G3', 'T3',
                        'q4','weight4','AC4', 'AG4', 'AT4', 'CG4', 'CT4', 'A4', 'C4', 'G4', 'T4',
                        'q5','weight5','AC5', 'AG5', 'AT5', 'CG5', 'CT5', 'A5', 'C5', 'G5', 'T5',
                        'end'])
    
    
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
    
r_list = list(np.arange(5,50,1))  
#f_list = list(np.arange(1,10,1))    
#instance = np.random.RandomState(seed = 20230803) 
np.random.seed(20230803)





for r in range(len(tuple_list)):
    classes, rates, length, ntaxa, replicates = tuple_list[r]
    # dataset_name
    file_name = 'c' + str(classes) + '_r' + str(rates) + '_l' + str(length) + '_t' + str(ntaxa) + '_rep' + str(replicates)
    para_row = [file_name,classes,rates,length,ntaxa,replicates]+['']*58
    # gtr parameters
    if classes == 1:
        ran_para = np.random.choice(r_list,6)
        f_list = np.random.uniform(1, 10, 4)
        ran_f = list(f_list/sum(f_list))
        AC,AG,AT,CG,CT = ran_para[0]/ran_para[5], ran_para[1]/ran_para[5], ran_para[2]/ran_para[5], ran_para[3]/ran_para[5], ran_para[4]/ran_para[5]
        A,C,G,T = ran_f[0]/sum(ran_f), ran_f[1]/sum(ran_f), ran_f[2]/sum(ran_f), ran_f[3]/sum(ran_f)
        para_row[8:19] = 'GTR+F',1,AC,AG,AT,CG,CT,A,C,G,T
    else:        
        # class weights
        weight_list = np.random.uniform(0.1, 1, classes)
        weights = list(weight_list/sum(weight_list))
        # q matrices
        for i in range(classes):
            ran_para = np.random.choice(r_list,6)
            ran_f = np.random.choice(f_list,4)
            AC,AG,AT,CG,CT,A,C,G,T = ran_para[0]/ran_para[5], ran_para[1]/ran_para[5], ran_para[2]/ran_para[5], ran_para[3]/ran_para[5], ran_para[4]/ran_para[5], ran_f[0]/sum(ran_f), ran_f[1]/sum(ran_f), ran_f[2]/sum(ran_f), ran_f[3]/sum(ran_f)    
            para_row[i*11+8:i*11+19] = 'GTR+F',weights[i],AC,AG,AT,CG,CT,A,C,G,T
        
    # rate parameters
    if  rates == 2:
        gamma = expon.rvs(scale=1.0).round(5)
        para_row[6:8] = 0, gamma

    #record model parameters
    with open('modelpara.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(para_row) 
    