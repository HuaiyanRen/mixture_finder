import os
import numpy as np

# here are the variables that I test my functions.

# file_position ='d: && cd D:\masterresearch\congnato'
# iqtree_position = 'D:/masterresearch/iqtree2'
# file_name = 'alignment.nex'
# score_type = '(BIC)'
# repeats = 10

def optimize_classes(file_position, iqtree_position, file_name, score_type, repeats):
    '''This is an algorithm the find the classes of mixture model of lowest score.
    The funtion returns the optimal numbers of classes and prints all the score have been generated.'''
    
    '''file_position: cmd statement that open (should include 'cd' comment) the folder including sequence file. Use '\'.
       
       iqtree_position: the folder including iqtree files (iqtree2, iqtree2-click, libiomp5md.dll). Use '/'.
       
       file_name:the name of sequence file.
       
       score_type: '(AIC)', '(AICc)', '(BIC)'
       
       repeats: the times of fitting the model in each numbers of classes'''
# initial variables       
    the_score = 9999999999
    classes = 1
    new_score = 0
    score_list = []

# run class 1

    cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(classes) # build the cmd to run iqtree fitting mixture model

    tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

    score_list = iqtree_ten_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores

    new_score = np.min(score_list) # get the lowest score of 10 repeatitions
    
    print(classes, new_score, score_list) # list the score of each step
    
# add classes and loop
    while new_score < the_score: # if getting lower new score when adding classes, keep running
        
        the_score = new_score # accept the lower score
        classes += 1 # add new class
    
        cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(classes) # build the cmd to run iqtree fitting mixture model
        
        tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step
        
        score_list = iqtree_ten_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores
    
        new_score = np.min(score_list) # get the lowest score of 10 repeatitions
        
        print(classes, new_score, score_list) # list the score of each step
    
# once the score increasing, return the nubmer classes with lowest score.
    return classes-1


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
        return ' -m GTR+FO -pre 1class -redo'
    else:
        repeats = ''
        for i in range(classes-1):
            repeats = repeats + ',GTR+FO'
        return ' -m "MIX{GTR+FO' + repeats + '}" -pre ' + str(classes) + 'classes -redo'
    
def tree_file_name(file_position,classes):
    ''' return right cmd statement for renaming the iqtree files with their numbers of classes'''
    if classes == 1:
        return file_position.split()[-1] + '/' + str(classes) +'class.iqtree'
    else:
        return file_position.split()[-1] + '/' + str(classes) +'classes.iqtree'

def iqtree_ten_repeatitions(cmd, tree_file, score_type, repeats):
    '''repeat the iqtree commend ten times and return a list of ten scores'''
    score_list = []

    for i in range(repeats):

        os.system(cmd)
    
        new_score = float(find_score(tree_file, score_type))
    
        #print(new_score)
    
        score_list.append(new_score)
        
    return score_list  






