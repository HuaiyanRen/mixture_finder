from multiprocessing import Pool
from functools import partial
import os
import argparse
import numpy as np

type_dict = {0:'JC+FO',
             1:'F81+FO',
             2:'K80+FO',
             3:'HKY+FO',
             4:'TN+FO',
             5:'TNe+FO',
             6:'K81+FO',
             7:'K81u+FO',
             8:'TPM2+FO',
             9:'TPM2u+FO',
             10:'TPM3+FO',
             11:'TPM3u+FO',
             12:'TIM+FO',
             13:'TIMe+FO',
             14:'TIM2+FO',
             15:'TIM2e+FO',
             16:'TIM3+FO',
             17:'TIM3e+FO',
             18:'TVM+FO',
             19:'TVMe+FO',
             20:'SYM+FO',
             21:'GTR+FO'}


def optimize_classes(file_position, iqtree_position, file_name, score_type, repeats, pool_num, type_dict):
    '''This is an algorithm the find the classes of mixture model of lowest score.
    The funtion returns the optimal numbers of classes and prints the full model.
    Module 'os' and 'numpy' are required. '''
    
    '''file_position: cmd statement that open (should include 'cd' comment) the folder including sequence file. Use '\'.
       
       iqtree_position: the folder including iqtree files (iqtree2, iqtree2-click, libiomp5md.dll). Use '/'.
       
       file_name: the name of sequence file.
       
       score_type: AIC, AICc, BIC.
       
       repeats: the times of fitting the model in each fitting
       
       pool_num: the maximum of multiprocessing pool.'''

# create a folder for results
    os.system('mkdir ' + file_name + '.results')
    os.system('mkdir ' + file_name + '.results/logs')
    with open(file_name + '.results/result.txt', 'w+') as result:
        result.write('-step modeltype lowest_score score_lists \n \n')

# initial variables
    repeats = int(repeats)
    pool_num = int(pool_num)
    score_type = '(' + score_type + ')'    
    
    last_score = 9999999999 # the score of last step, used for comparing
    classes = 1 # start with 1 class
    new_score = 0 # the lowest score generated in each step
    model_num_seq = [] # record submodels in ecah class
    
# adding class 1
    sub_min_list = repeatitions_a(repeats, type_dict, pool_num, file_name, classes, score_type, model_num_seq) # repeatly run iqtree and list the scores    
    score_min = np.min(sub_min_list) # get the lowest score among different submodels

    
    # print and record results in adding class 1
    model_num = sub_min_list.index(score_min) # get the submodel type of lowest score    
    model_num_seq.append(model_num) # record the submodel type
    new_score = score_min
    with open(file_name + '.results/result.txt', 'a+') as result:   
        print(str(classes)+' class,','the lowest score:'+ str(new_score), print_model_names(model_num_seq), '\n', file = result) # list the score of each step
    
# adding one more class and loop
    while new_score < last_score: # if getting lower new score when adding classes, keep running
        
        last_score = new_score # record the score of last step
        classes += 1 # add new class

        sub_min_list = repeatitions_a(repeats, type_dict, pool_num, file_name, classes, score_type, model_num_seq) # repeatly run iqtree and list the scores    
        score_min = np.min(sub_min_list) # get the lowest score among different submodels

        # print and record results in adding one more class
        model_num = sub_min_list.index(score_min) # get the submodel type of lowest score    
        model_num_seq.append(model_num) # record the submodel type
        new_score = score_min
        with open(file_name + '.results/result.txt', 'a+') as result:
            print(str(classes)+' classes,','the lowest score:'+ str(new_score), print_model_names(model_num_seq), '\n', file = result) # list the score of each step
        
    # once add a new class, try other sub models in previous classes
        if new_score < last_score: # criteria of succussfully adding a new class  
            model_num_seq_last = list(model_num_seq) # record the last full model
            #model_num_seq_traceback = list(model_num_seq) 
            pos = 0 # position of the changing class 
            
            # check class 1
            sub_min_list = repeatitions_r(repeats, type_dict, pool_num, file_name, classes, score_type, model_num_seq, pos) # repeatly run iqtree and list the scores    
            score_min = np.min(sub_min_list) # get the lowest score among different submodels
            
            # record results in changing submodels in class 1
            model_num = sub_min_list.index(score_min) # get the submodel type of lowest score           
            if score_min < new_score: # once find better score, select the model                        
                new_score = score_min
                model_num_seq[pos] = model_num
            
            # print the change
            with open(file_name + '.results/result.txt', 'a+') as result:   
                print('traceback result for class ' + str(pos+1) +': ', print_model_names(model_num_seq), 'with the score ' + str(new_score), '\n', file = result)
            
        # if class 1 is changed, check classes after          
            while model_num_seq[pos] != model_num_seq_last[pos]:                
                pos += 1
                if pos < len(model_num_seq) - 1:
                
                    sub_min_list = repeatitions_r(repeats, type_dict, pool_num, file_name, classes, score_type, model_num_seq, pos) # repeatly run iqtree and list the scores    
                    score_min = np.min(sub_min_list) # get the lowest score among different submodels
            
                    # record results in changing submodels in class 1
                    model_num = sub_min_list.index(score_min) # get the submodel type of lowest score           
                    if score_min < new_score: # once find better score, select the model                        
                        new_score = score_min
                        model_num_seq[pos] = model_num

                    # print the change
                    with open(file_name + '.results/result.txt', 'a+') as result:   
                        print('traceback result for class ' + str(pos+1) +': ', print_model_names(model_num_seq), 'with the score ' + str(new_score), '\n', file = result)
                    
# once the score increasing, when adding a new class, return the nubmer classes with lowest score.
    model_num_seq.pop() # delete the added new class with increasing score
    # the final result
    with open(file_name + '.results/result.txt', 'a+') as result: 
        print('the optimal number of classes is:'+ str(classes-1), print_model_names(model_num_seq), file = result)
    
    return 'ALGORITHM FINISH'


# sub_functions

def model_select(model_num_seq):
    '''the parameter is a list of number and return right cmd statement for fitting mixture models'''
    model_type_seq = []
    for n in range(len(model_num_seq)):
        model_type_seq.append(type_dict[model_num_seq[n]])
    return ' -m "MIX{' + ','.join(model_type_seq) + '}"'

def print_model_names(model_num_seq):
    '''the parameter is a list of number the function returns corresponding models'''
    model_type_seq = []
    for n in range(len(model_num_seq)):
        model_type_seq.append(type_dict[model_num_seq[n]])
    return ','.join(model_type_seq) 


def iqtree_running_a(type_order, file_name, classes, score_type, model_num_seq):    
    
    model_type, order_num = type_order
    log_file = file_name + '.results/logs/' + str(classes) + 'classes.' + type_dict[model_type] + str(order_num)
    model_num_seq_loop = model_num_seq + [model_type]
    
    cmd = 'iqtree2 -s ' + file_name + model_select(model_num_seq_loop) + ' -pre ' + log_file + ' -redo -nt 1'
    os.system(cmd)
    
    tree_file = log_file + '.iqtree'
    with open(tree_file) as iqtree_result:
        for line in iqtree_result.readlines():
            if score_type in line:
                list_line = line.split()
                score = list_line[-1]
    return float(score)

def iqtree_running_r(type_order, file_name, classes, score_type, model_num_seq, pos):    
    
    model_type, order_num = type_order
    log_file = file_name + '.results/logs/' + str(classes) + 'classes.' + type_dict[model_type] + str(order_num)
    model_num_seq_loop = list(model_num_seq)
    model_num_seq_loop[pos] = model_type
    
    cmd = 'iqtree2 -s ' + file_name + model_select(model_num_seq_loop) + ' -pre ' + log_file + ' -redo -nt 1'
    os.system(cmd)
    
    tree_file = log_file + '.iqtree'
    with open(tree_file) as iqtree_result:
        for line in iqtree_result.readlines():
            if score_type in line:
                list_line = line.split()
                score = list_line[-1]
    return float(score)   

def repeatitions_a(repeats, type_dict, pool_num, file_name, classes, score_type, model_num_seq):
    
    type_list = list(type_dict)
    type_repeats = [m for m in type_list for i in range(repeats) ]   
    order_num = list(np.arange(1, repeats + 1, 1))*len(type_list)
    
    type_order = [0]*len(type_repeats)
    for i in range(len(type_repeats)):
        type_order[i] = type_repeats[i], order_num[i]

    partial_iqtree_running = partial(iqtree_running_a, 
                                     file_name = file_name, 
                                     classes = classes, 
                                     score_type = score_type,
                                     model_num_seq = model_num_seq)
    with Pool(pool_num) as p:
        score_list = p.map(partial_iqtree_running, type_order)

    sub_min_list = []
    for i in range(len(type_list)):
        sub_score_list = score_list[i*repeats: (i+1)*repeats]
        sub_min = np.min(sub_score_list)
        sub_min_list.append(sub_min)
        with open(file_name + '.results/result.txt', 'a+') as result:
            print('-class:'+ str(classes), type_dict[i], sub_min, sub_score_list, file = result) # list all the scores generated
    return sub_min_list

def repeatitions_r(repeats, type_dict, pool_num, file_name, classes, score_type, model_num_seq, pos):
    
    type_list = list(type_dict)
    type_repeats = [m for m in type_list for i in range(repeats) ]   
    order_num = list(np.arange(1, repeats + 1, 1))*len(type_list)
    type_order = [0]*len(type_repeats)
    for i in range(len(type_repeats)):
        type_order[i] = type_repeats[i], order_num[i]

    partial_iqtree_running = partial(iqtree_running_r, 
                                     file_name = file_name, 
                                     classes = classes, 
                                     score_type = score_type,
                                     model_num_seq = model_num_seq,
                                     pos = pos)
    with Pool(pool_num) as p:
        score_list = p.map(partial_iqtree_running, type_order)
        
    sub_min_list = []
    for i in range(len(type_list)):
        sub_score_list = score_list[i*repeats: (i+1)*repeats]
        sub_min = np.min(sub_score_list)
        sub_min_list.append(sub_min)
        with open(file_name + '.results/result.txt', 'a+') as result:
            print('-traceback class:'+ str(pos+1), type_dict[i], sub_min, sub_score_list, file = result) # list all the scores generated
    return sub_min_list   


# running
parser = argparse.ArgumentParser(description='find the optimal number of classes in mixture models')
parser.add_argument('--file_position', '-fp', help='cmd statement that open (should include \'cd\' comment) the folder including sequence file.  Use \'\ \'. \n Used in windows',
                    default='')
parser.add_argument('--iqtree_position', '-tp', help='the position of the folder including iqtree files (iqtree2, iqtree2-click, libiomp5md.dll). Use \'/\'. \n Used in windows',
                    default='')
parser.add_argument('--file_name', '-n', help='the name of sequence file', 
                    required=True)
parser.add_argument('--score_type', '-s', help='one of AIC, AICc, BIC', 
                    required=True)
parser.add_argument('--repeats', '-r', help='the times of fitting the model in each fitting', 
                    default=5)
parser.add_argument('--pool_num', '-p', help='the maximum of multiprocessing pool.', 
                    default=10)
args = parser.parse_args()

if __name__ == '__main__':
    try:
        optimize_classes(args.file_position, args.iqtree_position, args.file_name, args.score_type, args.repeats, args.pool_num, type_dict)
    except Exception as e:
        print(e)


# here are the variables that I test my functions.

#file_position ='d: && cd D:\masterresearch\congnato'
#iqtree_position = 'D:/masterresearch/iqtree2'
#file_name = 'alignment_COI_3rdpos-out.nex'
#score_type = '(BIC)'
#repeats = 10