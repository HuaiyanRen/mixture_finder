import numpy as np
import pandas as pd
import csv 
import os


with open('result_weight.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes',  'ntaxa', 'sites', 'invariable',
                        'tw1', 'tw2', 'tw3', 'tw4', 'tw5',
                        'ew1', 'ew2', 'ew3', 'ew4', 'ew5', 
                        'sort_tw', 'sort_ew', 'rmse'])
                       
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
    df = pd.read_csv('modelpara.csv')    
    matching_rows = df[df.iloc[:, 0] == file_name]
    name_index = matching_rows.index[0]
    return name_index

def extract_q(row_index, col_names):
    df = pd.read_csv('modelpara.csv')
    
    row_data = df.loc[row_index, col_names]
    return row_data 

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
                

    tw_list = []
    row_index = file_name_row(file_name)
    tw1 = extract_q(row_index, 'weight1')
    tw_list.append(tw1)
    tw2 = extract_q(row_index, 'weight2')
    if isinstance(tw2, float) and str(tw2) == 'nan':
        tw_list.append('')
    else:
        tw_list.append(tw2)
    tw3 = extract_q(row_index, 'weight3')
    if isinstance(tw3, float) and str(tw3) == 'nan':
        tw_list.append('')
    else:
        tw_list.append(tw3)
    tw4 = extract_q(row_index, 'weight4')
    if isinstance(tw4, float) and str(tw4) == 'nan':
        tw_list.append('')
    else:
        tw_list.append(tw4)
    tw5 = extract_q(row_index, 'weight5')
    if isinstance(tw5, float) and str(tw5) == 'nan':
        tw_list.append('')
    else:
        tw_list.append(tw5)
      
    
    result_row = [file_name, classes, ntaxa, length, invariable] + tw_list
    
    if os.path.isfile('cor/'+ file_name + '.iqtree'):
        with open('cor/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'odel of substitution:' in line:
                    if int(classes) > 1:
                        est_q_list = line.split('{')[-1].split('}')[0].split(',')
                        for q in range(int(classes)):
                            if '+F' in est_q_list[q]:
                                est_q_list[q] = est_q_list[q].replace("+FO", "")
    
    ew_list = ['','','','','']
    if int(classes) == 1:
        ew_list[0] = 1
    else:
        with open('cor/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                for i in range(1,int(classes)+1):
                    if str(i) + '  '+ est_q_list[i-1] in line:
                        ew_list[i-1] = float(line.split()[-2])

    result_row = result_row + ew_list

    tw_list = [i for i in tw_list if i != '']
    ew_list = [i for i in ew_list if i != '']
    tw_list = sorted(tw_list, reverse=True)
    ew_list = sorted(ew_list, reverse=True)
    
    rmse = np.sqrt(np.mean((np.array(tw_list) - np.array(ew_list)) ** 2))

    result_row = result_row + [tw_list, ew_list, rmse]
    
    
    with open('result_weight.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(result_row)         
    
