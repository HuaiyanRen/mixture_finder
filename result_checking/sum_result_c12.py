import numpy as np
import csv 
from scipy.stats import chi2

with open('result_c12.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes',  'ntaxa', 'sites', 'invariable',
                        'df1', 'llh1', 'df2', 'llh2', 'p_value', 'op1', 'op2', 'op3', 'op4'])
                       
classes = [1] 
rates = [2] # 0: +E, 1: +I, 2: +I+G
length = [1000, 10000] #[1000,2000,5000,10000]
ntaxa = [100]
replicates = list(np.arange(0,250,1))

classes_list = [m for m in classes for i in range(len(rates)*len(length)*len(ntaxa)*len(replicates))]
rates_list = [m for m in rates for i in range(len(length)*len(ntaxa)*len(replicates))]*len(classes)
length_list = [m for m in length for i in range(len(ntaxa)*len(replicates))]*len(classes)*len(rates)
ntaxa_list = [m for m in ntaxa for i in range(len(replicates))]*len(classes)*len(rates)*len(length)
replicates_list = replicates*len(classes)*len(rates)*len(length)*len(ntaxa)
    
tuple_list = [0]*len(classes)*len(rates)*len(length)*len(ntaxa)*len(replicates)

for i in range(len(tuple_list)):
    tuple_list[i] = classes_list[i], rates_list[i], length_list[i], ntaxa_list[i], replicates_list[i]
    
    
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
    
    df_list = []
    llh_list = []
    
    with open('two/' + file_name + '.log') as b:
        for line in b.readlines():
            if '; df:' in line:
                df = float(line.split(';')[1].split()[-1])
                llh = float(line.split(';')[2].split()[-1])
                df_list.append(df)
                llh_list.append(llh)
                
    p_value = chi2.cdf(2*(llh_list[1] - llh_list[0]), df_list[1] - df_list[0])
    
    op1 = 1
    op2 = 1
    op3 = 1
    op4 = 1
    if p_value > 0.9:
        op1= 2
    if p_value > 0.95:
        op2= 2
    if p_value > 0.99:
        op3= 2
    if p_value > 0.999:
        op4= 2
    
    
    result_row = result_row + [df_list[0], llh_list[0], df_list[1], llh_list[1], 1-p_value, op1, op2, op3, op4]
    
    with open('result_c12.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(result_row)    
    
