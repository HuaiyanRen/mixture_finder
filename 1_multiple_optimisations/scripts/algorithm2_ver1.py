import os
import numpy as np

# here are the variables that I test my functions.

file_position ='d: && cd D:\masterresearch\congnato'
iqtree_position = 'D:/masterresearch/iqtree2'
file_name = 'subalignment.nex'
score_type = '(BIC)'
repeats = 10

def optimize_classes(file_position, iqtree_position, file_name, score_type, repeats):
    '''This is an algorithm the find the classes of mixture model of lowest score.
    The funtion returns the optimal numbers of classes and prints all the score have been generated.'''
    
    '''file_position: cmd statement that open (should include 'cd' comment) the folder including sequence file. Use '\'.
       
       iqtree_position: the folder including iqtree files (iqtree2, iqtree2-click, libiomp5md.dll). Use '/'.
       
       file_name:the name of sequence file.
       
       score_type: '(AIC)', '(AICc)', '(BIC)'
       
       repeats: the times of fitting the model in each numbers of classes'''
# initial variables       
    new_score = 0
    score_list = []
    score_list_loop = []

# run classes from 1 to 10
    for classes in range(1,11):

        cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(classes) # build the cmd to run iqtree fitting mixture model

        tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

        score_list_loop = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores

        new_score = np.min(score_list_loop) # get the lowest score of 10 repeatitions
    
        score_list.append(new_score) # list the lowsest score of each class
    
        print(str(classes)+' class,','the lowest score:'+ str(new_score), score_list_loop) # list the score of each step
    
# find the number of classes with lowest score
    score_min = np.min(score_list)
    classes = score_list.index(score_min) + 1

    print('the optimal number of classes is:'+ str(classes))

# optimaze the model in fixed number of classes  
    type_dict = {0:'GTR+FO', 1:'JC+FO', 2: 'HKY+FO'}
    model_type_list = [0] * classes
    
    # change the fisrt class
    for i in range(1,len(type_dict)): # try other sub models except GTR
        model_type_list_loop = list(model_type_list)
        model_type_list_loop[0] = i

        cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select_type(model_type_list) # build the cmd to run iqtree fitting mixture model
        
        tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

        score_list_loop = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores
        
        new_score = np.min(score_list_loop)
        
        print(print_model_names(model_type_list_loop), 'the lowest score:'+ str(new_score), score_list_loop) # list the score of each step
        
        if new_score < score_min:
            score_min = new_score
            model_type_list[0] = i
            
    # if model optimazed in last step, keep on change
    class_position = 0
    while model_type_list[class_position] != 0 and class_position < len(model_type_list) - 1:
        class_position += 1
        
        for i in range(1,len(type_dict)): # try other sub models except GTR
            model_type_list_loop = list(model_type_list)
            model_type_list_loop[class_position] = i

            cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select_type(model_type_list) # build the cmd to run iqtree fitting mixture model
        
            tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

            score_list_loop = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores
        
            new_score = np.min(score_list_loop)
        
            print(print_model_names(model_type_list), 'the lowest score:'+ str(new_score), score_list_loop) # list the score of each step
        
            if new_score < score_min:
                score_min = new_score
                model_type_list[class_position] = i
                
# return the nubmer classes with lowest score.
    return 'the optimal number of classes is:'+ str(classes), 'the lowest score is:'+ str(score_min), 'the full model is:' + print_model_names(model_type_list)


# sub_functions

def find_score(filename,score_type):
    '''return AIC or BIC score from iqtree file'''
    with open(filename) as iqtree_result:
        for line in iqtree_result.readlines():
            if score_type in line:
                list_line = line.split()
                score = list_line[-1]
    return score

def model_select(classes):
    ''' return right cmd statement for fitting mixture models of correct numbers classes'''
    if classes == 1:
        return ' -m GTR+FO -pre 1classes -redo'
    else:
        repeats = ''
        for i in range(classes-1):
            repeats = repeats + ',GTR+FO'
        return ' -m "MIX{GTR+FO' + repeats + '}" -pre ' + str(classes) + 'classes -redo'
    
def tree_file_name(file_position,classes):
    ''' return right cmd statement for renaming the iqtree files with their numbers of classes'''
    return file_position.split()[-1] + '/' + str(classes) +'classes.iqtree'

def iqtree_repeatitions(cmd, tree_file, score_type, repeats):
    '''repeat the iqtree commend ten times and return a list of ten scores'''
    score_list = []

    for i in range(repeats):

        os.system(cmd)
    
        new_score = float(find_score(tree_file, score_type))
    
        score_list.append(new_score)
        
    return score_list  

def model_select_type(model_num_seq):
    '''the parameter is a list of number and return right cmd statement for fitting mixture models'''
    model_type_seq = []
    for n in range(len(model_num_seq)):
        model_type_seq.append(model_type(model_num_seq[n]))
    return ' -m "MIX{' + ','.join(model_type_seq) + '}" -pre ' + str(len(model_num_seq)) + 'classes -redo'

def print_model_names(model_num_seq):
    '''the parameter is a list of number the function returns corresponding models'''
    model_type_seq = []
    for n in range(len(model_num_seq)):
        model_type_seq.append(model_type(model_num_seq[n]))
    return ','.join(model_type_seq) 
        
def model_type(numbers):
    '''This function receives a number and returns corresponding model, based on type_dict'''
    type_dict = {0:'GTR+FO', 1:'JC+FO', 2: 'HKY+FO'}
    return type_dict[numbers]