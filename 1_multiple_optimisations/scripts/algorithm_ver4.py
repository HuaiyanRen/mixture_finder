import os
import numpy as np

# here are the variables that I test my functions.

#file_position ='d: && cd D:\masterresearch\congnato'
#iqtree_position = 'D:/masterresearch/iqtree2'
#file_name = 'subalignment.nex'
#score_type = '(BIC)'
#repeats = 10

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
    model_num_seq = []
    type_dict = {0:'JC+FO', 1:'HKY+FO', 2: 'GTR+FO'} # this should be same with sub model 'model type'. further model types could be added
    
# run class 1
    score_list_loop = []

    for i in range(len(type_dict)): # run all model types in class 1
        
        model_num_seq_loop = model_num_seq + [i]
        
        cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model
        
        tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

        score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores
        
        min_of_list = np.min(score_list)
        
        print('class:'+ str(classes), type_dict[i], min_of_list, score_list) # list all the scores generated

        score_list_loop.append(min_of_list) # get the lowest score of 10 repeatitions
        
    score_min = np.min(score_list_loop) # get the lowest score of different models
    
    model_num = score_list_loop.index(score_min) # get the model type of lowest model
    
    model_num_seq.append(model_num) # record the model
    
    new_score = score_min
    
    print(str(classes)+' class,','the lowest score:'+ str(new_score), print_model_names(model_num_seq)) # list the score of each step
    
# add classes and loop
    while new_score < the_score: # if getting lower new score when adding classes, keep running
        
        the_score = new_score # accept the lower score
        classes += 1 # add new class
        score_list_loop = []
        
        for i in range(len(type_dict)): # run all model types in added class
        
            model_num_seq_loop = model_num_seq + [i]
        
            cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model
        
            tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

            score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores
            
            min_of_list = np.min(score_list)
            
            print('class:'+str(classes), type_dict[i], min_of_list, score_list) # list all the scores generated

            score_list_loop.append(min_of_list) # get the lowest score of 10 repeatitions
    
        score_min = np.min(score_list_loop) # get the lowest score of different models
    
        model_num = score_list_loop.index(score_min) # get the model type of lowest model
        
        model_num_seq.append(model_num) # record the model
    
        new_score = score_min
    
        print(str(classes)+' classes,','the lowest score:'+ str(new_score), print_model_names(model_num_seq)) # list the score of each step
        
# once add a new class, try other sub models in previous classes
        if new_score < the_score: 
        
            model_num_seq_origin = list(model_num_seq) # record the model before checking
        
# check class 1           
            for i in range(len(type_dict)): # run all models in the class
            
                model_num_seq_change = list(model_num_seq)
                model_num_seq_change[0] = i # start with class 1
                        
                cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_change) # build the cmd to run iqtree fitting mixture model
        
                tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

                score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores
                    
                min_of_list = np.min(score_list)

                print('traceback', print_model_names(model_num_seq_change), min_of_list, score_list)
                        
                if min_of_list < new_score: # once find better score, select the model
                        
                    new_score = min_of_list
                    model_num_seq = list(model_num_seq_change)
                        
                    print('better score:' + str(new_score),print_model_names(model_num_seq))
            
# if class 1 is changed, check classes after
            pos = 0
           
            while model_num_seq[pos] != model_num_seq_origin[pos]:
                
                pos += 1
                if pos < len(model_num_seq) - 1:
                
                    for i in range(len(type_dict)): # run all models in classes
                              
                        model_num_seq_change = list(model_num_seq)
                        model_num_seq_change[pos] = i
                        
                        cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_change) # build the cmd to run iqtree fitting mixture model
        
                        tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step

                        score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # run iqtree 10 times and list the scores
                    
                        min_of_list = np.min(score_list)

                        print('traceback', print_model_names(model_num_seq_change), min_of_list, score_list)
                        
                        if min_of_list < new_score: # once find better score, select the model
                        
                            new_score = min_of_list
                            model_num_seq = list(model_num_seq_change)
                        
                            print('better score:' + str(new_score),print_model_names(model_num_seq))
        
# once the score increasing, return the nubmer classes with lowest score.
    model_num_seq.pop()

    return 'the optimal number of classes is:'+ str(classes-1), print_model_names(model_num_seq)

# sub_functions

def find_score(filename,score_type):
    '''return AIC or BIC score from iqtree file'''
    with open(filename) as iqtree_result:
        for line in iqtree_result.readlines():
            if score_type in line:
                list_line = line.split()
                score = list_line[-1]
    return score

def model_select(model_num_seq):
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
    type_dict = {0:'JC+FO', 1:'HKY+FO', 2: 'GTR+FO'}
    return type_dict[numbers]

def tree_file_name(file_position,classes):
    ''' return right cmd statement for renaming the iqtree files with their numbers of classes'''
    return file_position.split()[-1] + '/' + str(classes) +'classes.iqtree'

def iqtree_repeatitions(cmd, tree_file, score_type, repeats):
    '''repeat the iqtree commend ten times and return a list of ten scores'''
    score_list = []

    for i in range(repeats):

        os.system(cmd)
    
        new_score = float(find_score(tree_file, score_type))
    
        #print(new_score)
    
        score_list.append(new_score)
        
    return score_list 