from multiprocessing import Pool
from functools import partial
import subprocess
import numpy as np
import pandas as pd
import argparse
#import ast

def file_name_row(file_name):
    df = pd.read_csv('modelpara.csv')    
    matching_rows = df[df.iloc[:, 0] == file_name]
    name_index = matching_rows.index[0]
    return name_index

def extract_q(row_index, col_names):
    df = pd.read_csv('modelpara.csv')
    
    row_data = df.loc[row_index, col_names]
    return row_data 

def running_iqtree(classes, rates, length, ntaxa, rep1, rep2, method, pool_num):
    #classes = list(ast.literal_eval(classes))
    #rates = list(ast.literal_eval(rates))
    #length = list(ast.literal_eval(length))
    #ntaxa = list(ast.literal_eval(ntaxa))    
    rep1 = int(rep1)
    rep2 = int(rep2)
    #method = int(method)
    pool_num = int(pool_num)
    replicates = list(np.arange(rep1,rep2,1))

    classes_list = [m for m in classes for i in range(len(rates)*len(length)*len(ntaxa)*len(replicates)*len(method))]
    rates_list = [m for m in rates for i in range(len(length)*len(ntaxa)*len(replicates)*len(method))]*len(classes)
    length_list = [m for m in length for i in range(len(ntaxa)*len(replicates)*len(method))]*len(classes)*len(rates)
    ntaxa_list = [m for m in ntaxa for i in range(len(replicates)*len(method))]*len(classes)*len(rates)*len(length)
    replicates_list = [m for m in replicates for i in range(len(method))]*len(classes)*len(rates)*len(length)*len(ntaxa)
    method_list = method*len(classes)*len(rates)*len(length)*len(ntaxa)*len(replicates)
    
    tuple_list = [0]*len(classes)*len(rates)*len(length)*len(ntaxa)*len(replicates)*len(method)

    for i in range(len(tuple_list)):
        tuple_list[i] = classes_list[i], rates_list[i], length_list[i], ntaxa_list[i], replicates_list[i], method_list[i] 
    partial_running = partial(running_tuple)
    with Pool(pool_num) as p:
        p.map(partial_running, tuple_list)
    
def running_tuple(tuple_list):
    classes, rates, length, ntaxa, replicates, method = tuple_list
    #file_name
    file_name = 'c' + str(classes) + '_r' + str(rates) + '_l' + str(length) + '_t' + str(ntaxa) + '_rep' + str(replicates)
    #model
    if int(classes) == 1:
        row_index = file_name_row(file_name)
        q1 = extract_q(row_index, 'q1')
        q1 = q1.replace('+F', '+FO')
        model_cmd = q1 + '+I+G'
    elif int(classes) > 1:
        row_index = file_name_row(file_name)
        q_list = []
        q1 = extract_q(row_index, 'q1')
        q1 = q1.replace('+F', '+FO')
        q_list.append(q1)
        q2 = extract_q(row_index, 'q2')
        if not isinstance(q2, float):
            q2 = q2.replace('+F', '+FO')
            q_list.append(q2)
            q3 = extract_q(row_index, 'q3')
            if not isinstance(q3, float):
                q3 = q3.replace('+F', '+FO')
                q_list.append(q3)
                q4 = extract_q(row_index, 'q4')
                if not isinstance(q4, float):
                    q4 = q4.replace('+F', '+FO')
                    q_list.append(q4)
                    q5 = extract_q(row_index, 'q5')
                    if not isinstance(q5, float):
                        q5 = q5.replace('+F', '+FO')
                        q_list.append(q5)
        model_cmd = 'MIX"{'
        for i in range(len(q_list)):
            model_cmd = model_cmd + 'GTR+FO,'
        model_cmd = model_cmd[:-1]
        model_cmd = model_cmd + '}+I+G"'    
    #run
    cmd = '/usr/bin/time -v /mnt/data/dayhoff/home/u7151703/software/iqtree-2.2.6.mix-Linux/bin/iqtree2 -s ' + file_name + '.fa -m ' + model_cmd + ' -pre cor/'+ file_name + ' -nt 1'
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    with open('cor/' + file_name + '_time.txt', 'w') as f:
        #f.write(result.stdout)
        f.write(result.stderr)
        
    
# running
parser = argparse.ArgumentParser(description='')
parser.add_argument('--classes', '-c', help='', nargs='+',
                    required=True)
parser.add_argument('--rates', '-r', help='', nargs='+',
                    required=True)
parser.add_argument('--length', '-l', help='', nargs='+',
                    required=True)
parser.add_argument('--ntaxa', '-t', help='', nargs='+',
                    required=True)
parser.add_argument('--rep1', '-r1', help='', 
                    default = 0)
parser.add_argument('--rep2', '-r2', help='', 
                    required=True)
parser.add_argument('--method', '-m', help='', nargs='+',
                    required=True)
parser.add_argument('--pool_num', '-p', help='', 
                    required=True)
args = parser.parse_args()

if __name__ == '__main__':
    try:
        running_iqtree(args.classes, args.rates, args.length, args.ntaxa, args.rep1, args.rep2, args.method, args.pool_num)
    except Exception as e:
        print(e)

#classes = [1,2,3,4,5] 
#rates = [0,1,2] # 0: +E, 1: +I, 2: +I+G
#length =list(np.arange(1000,10001,1000))
#ntaxa = [100]
#replicates = list(np.arange(0,10,1))