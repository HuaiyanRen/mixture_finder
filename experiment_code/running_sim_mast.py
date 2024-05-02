from multiprocessing import Pool
from functools import partial
import subprocess
import numpy as np
import argparse
#import ast

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
    #run
    cmd = '/usr/bin/time -v /data/huaiyan/software/iqtree-2.2.6.mix-Linux/bin/iqtree2 -s ' + file_name + '.fa -m MIX+MFP -lrt 0 -merit BIC -mrate E,I,G,I+G,R,I+R -pre allt/'+ file_name + ' -nt 1'
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    with open('allt/' + file_name + '_time.txt', 'w') as f:
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