import os
import numpy as np

def optimize_classes(file_position, iqtree_position, file_name, score_type, repeats):
    '''This is an algorithm the find the classes of mixture model of lowest score.
    The funtion returns the optimal numbers of classes and prints the full model.
    Module 'os' and 'numpy' are required. '''
    
    '''file_position: cmd statement that open (should include 'cd' comment) the folder including sequence file. Use '\'.
       
       iqtree_position: the folder including iqtree files (iqtree2, iqtree2-click, libiomp5md.dll). Use '/'.
       
       file_name: the name of sequence file.
       
       score_type: '(AIC)', '(AICc)', '(BIC)'
       
       repeats: the times of fitting the model in each fitting'''
# create a folder for results
    os.system('mkdir results')
    with open('results/result.txt', 'w+') as result:
        result.write('-step modeltype lowest_score score_lists \n \n')

# initial variables       
    last_score = 9999999999 # the score of last step, used for comparing
    classes = 1 # start with 1 class
    new_score = 0 # the lowest score generated in each step
    model_num_seq = [] # record submodels in ecah class
    type_dict =  {0:'JC+FO', 1:'HKY+FO', 2: 'GTR+FO'} # this should be same with the dict in sub function 'model_type()'. further model types could be added
    
# adding class 1
    score_list_loop = [] # record the lowest scores of all considered submodels

    # run all model types in class 1  
    for i in range(len(type_dict)): 
        
        model_num_seq_loop = model_num_seq + [i] # add type of submodels
        
        # run iqtree software
        #cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model
        cmd = 'iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model
        tree_file = tree_file_name(file_position, classes) # assign the iqtree file name generated in last step
        score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # repeatly run iqtree and list the scores
        min_of_list = np.min(score_list) # get the lowest score in repeatitions
        
        # print and record results in repeatition of running a submodel
        with open('results/result.txt', 'a+') as result:
            print('-class:'+ str(classes), print_model_names(model_num_seq_loop), min_of_list, score_list, file = result) # list all the scores generated
        score_list_loop.append(min_of_list) # record the lowest score

    # get the lowest score among different submodels        
    score_min = np.min(score_list_loop)   
    model_num = score_list_loop.index(score_min) # get the submodel type of lowest score    
    model_num_seq.append(model_num) # record the submodel type
    
    # print and record results in adding class 1
    new_score = score_min
    with open('results/result.txt', 'a+') as result:   
        print(str(classes)+' class,','the lowest score:'+ str(new_score), print_model_names(model_num_seq), '\n', file = result) # list the score of each step
    
# adding one more class and loop
    while new_score < last_score: # if getting lower new score when adding classes, keep running
        
        last_score = new_score # record the score of last step
        classes += 1 # add new class
        score_list_loop = [] # record the lowest scores of all considered submodels
 
        # run all types of submodel in added class
        for i in range(len(type_dict)):
        
            model_num_seq_loop = model_num_seq + [i] # add type of submodels

            # run iqtree software
            #cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model        
            cmd = 'iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model
            tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step
            score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # repeatly run iqtree and list the scores            
            min_of_list = np.min(score_list) # get the lowest score in repeatitions

            # print and record results in repeatition of running a submodel
            with open('results/result.txt', 'a+') as result:
                print('-class:'+str(classes), print_model_names(model_num_seq_loop), min_of_list, score_list, file = result) # list all the scores generated
            score_list_loop.append(min_of_list) # record the lowest score

        # get the lowest score among different submodels
        score_min = np.min(score_list_loop)    
        model_num = score_list_loop.index(score_min) # get the model type of lowest model        
        model_num_seq.append(model_num) # record the submodel type

        # print and record results in adding one more class
        new_score = score_min
        with open('results/result.txt', 'a+') as result:
            print(str(classes)+' classes,','the lowest score:'+ str(new_score), print_model_names(model_num_seq), '\n', file = result) # list the score of each step
        
    # once add a new class, try other sub models in previous classes
        if new_score < last_score: # criteria of succussfully adding a new class  
            model_num_seq_traceback = list(model_num_seq) # record candidated models
            pos = 0 # position of the changing class 
            
            # check class 1           
            for i in range(len(type_dict)): # run all models in the class

                # change types of submodel            
                model_num_seq_change = list(model_num_seq_traceback)
                model_num_seq_change[pos] = i

                # run iqtree software        
                #cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_change) # build the cmd to run iqtree fitting mixture model        
                cmd = 'iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model
                tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step
                score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # repeatly run iqtree and list the scores                    
                min_of_list = np.min(score_list) # get the lowest score in repeatitions

                # print results in repeatition of running a submodel
                with open('results/result.txt', 'a+') as result:
                    print('-traceback', print_model_names(model_num_seq_change), min_of_list, score_list, file = result)

                # compare the score of changing a submodel to last lowest score        
                if min_of_list < new_score: # once find better score, select the model                        
                    new_score = min_of_list
                    model_num_seq_traceback = list(model_num_seq_change)
                
            model_num_seq = list(model_num_seq_traceback)
            # print the change
            with open('results/result.txt', 'a+') as result:   
                print('traceback result for class ' + str(pos+1) +': ', print_model_names(model_num_seq_traceback), 'with the score ' + str(new_score), '\n', file = result)
            
        # if class 1 is changed, check classes after          
            while model_num_seq_traceback[pos] != model_num_seq[pos]:                
                pos += 1
                if pos < len(model_num_seq) - 1:
                
                    for i in range(len(type_dict)): # run all types of submodels in class

                        # change types of submodel                              
                        model_num_seq_change = list(model_num_seq_traceback)
                        model_num_seq_change[pos] = i

                        # run iqtree software
                        #cmd = file_position + ' && ' + iqtree_position + '/iqtree2 -s ' + file_name + model_select(model_num_seq_change) # build the cmd to run iqtree fitting mixture model        
                        cmd = 'iqtree2 -s ' + file_name + model_select(model_num_seq_loop) # build the cmd to run iqtree fitting mixture model
                        tree_file = tree_file_name(file_position,classes) # assign the iqtree file name generated in last step
                        score_list = iqtree_repeatitions(cmd, tree_file, score_type, repeats) # repeatly run iqtree and list the scores                    
                        min_of_list = np.min(score_list) # get the lowest score in repeatitions

                        # print results in repeatition of running a submodel
                        print('-traceback', print_model_names(model_num_seq_change), min_of_list, score_list)

                        # compare the score of changing a submodel to last lowest score                        
                        if min_of_list < new_score: # once find better score, select the model                        
                            new_score = min_of_list
                            model_num_seq_traceback = list(model_num_seq_change)
                            # print the change                        
                            with open('results/result.txt', 'a+') as result:   
                                print('better score:' + str(new_score), print_model_names(model_num_seq_traceback), '\n', file = result)
                
                    model_num_seq = list(model_num_seq_traceback)
                    # print the change
                    with open('results/result.txt', 'a+') as result:   
                        print('traceback result for class ' + str(pos+1) +': ', print_model_names(model_num_seq_traceback), 'with the score ' + str(new_score), '\n', file = result)
            
        
# once the score increasing, when adding a new class, return the nubmer classes with lowest score.
    model_num_seq.pop() # delete the added new class with increasing score
    # the final result
    with open('results/result.txt', 'a+') as result: 
        print('the optimal number of classes is:'+ str(classes-1), print_model_names(model_num_seq), file = result)
    
    return 'ALGORITHM FINISH'


# sub_functions

def model_select(model_num_seq):
    '''the parameter is a list of number and return right cmd statement for fitting mixture models'''
    model_type_seq = []
    for n in range(len(model_num_seq)):
        model_type_seq.append(model_type(model_num_seq[n]))
    return ' -m "MIX{' + ','.join(model_type_seq) + '}" -pre ' + str(len(model_num_seq)) + 'classes -redo'

def model_type(numbers):
    '''This function receives a number and returns corresponding model, based on type_dict'''
    type_dict =  {0:'JC+FO', 1:'HKY+FO', 2: 'GTR+FO'} # all GTR submodels can be add in. this dict should be same with the dict in main function
    return type_dict[numbers] 
def tree_file_name(file_position,classes):
    ''' return right cmd statement for renaming the iqtree files with their numbers of classes'''
    #return file_position.split()[-1] + '/' + str(classes) +'classes.iqtree'
    return str(classes) +'classes.iqtree'

def print_model_names(model_num_seq):
    '''the parameter is a list of number the function returns corresponding models'''
    model_type_seq = []
    for n in range(len(model_num_seq)):
        model_type_seq.append(model_type(model_num_seq[n]))
    return ','.join(model_type_seq)         

def iqtree_repeatitions(cmd, tree_file, score_type, repeats):
    '''repeat the iqtree commend ten times and return a list of ten scores'''
    score_list = []
    for i in range(repeats):
        print(cmd)
        os.system(cmd)    
        new_score = float(find_score(tree_file, score_type))    
        score_list.append(new_score)        
    return score_list 

def find_score(filename,score_type):
    '''return score from iqtree file'''
    with open(filename) as iqtree_result:
        for line in iqtree_result.readlines():
            if score_type in line:
                list_line = line.split()
                score = list_line[-1]
    return score


# here are the variables that I test my functions. Users should assign the variables and run the function

#file_position ='d: && cd D:\masterresearch\congnato'
file_position ='' # linux users don't need change this variable
#iqtree_position = 'D:/masterresearch/iqtree2'
iqtree_position = '' # linux users don't need change this variable
file_name = 'alignment.nex'
score_type = '(BIC)'
repeats = 10

# run
optimize_classes(file_position, iqtree_position, file_name, score_type, repeats)