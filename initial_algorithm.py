import os


# here are the variables that I test my functions.

# file_position = 'd: && cd D:\masterresearch\congnato'
# file_name ='alignment.nex'
# tree_file_position = 'D:/masterresearch/congnato'
# score_type = '(BIC)'


def optimize_classes(file_position, file_name, tree_file_position, score_type):
    '''This is an algorithm the find the classes of mixture model of lowest score.
    The funtion returns the count of classes and the lowest score.'''
    
    '''file_position: cmd statement that open the folder including BOTH sequence file and iqtree file.
       Here is a problem that if files for running iqtree are not in the same folder, python cannot run.
       
       file_name: name of sequence file.
       
       tree_file_position: it is the same folder as sequence file, but should be input with '/' in windows.
       
       score_type: '(AIC)', '(AICc)', '(BIC)' '''
# initial variables       
    the_score = 9999999999
    classes = 1
    new_score = 0

# run class 1

    cmd = file_position + '&& iqtree2 -s ' + file_name + model_select(classes) # build the cmd to run iqtree fitting mixture model

    os.system(cmd) # run
    
    tree_file = tree_file_name(tree_file_position,classes) # assign the iqtree file name generated in last step

    new_score = float(find_score(tree_file, score_type)) # assign the score of mixture model including class 1

# add classes and loop
    while new_score < the_score: # if getting lower new score when adding classes, keep running
        
        print(classes, new_score) # list the score of each step
        
        the_score = new_score # accept the lower score
        classes += 1 # add new class
    
        cmd = file_position + '&& iqtree2 -s ' + file_name + model_select(classes) # build the cmd to run iqtree fitting mixture model
    
        os.system(cmd) # run 
        
        tree_file = tree_file_name(tree_file_position,classes) # assign the iqtree file name generated in last step
    
        new_score = float(find_score(tree_file, score_type)) # assign the score of mixture model including multiple classes
    
# once the score increasing, return the nubmer classes with lowest score.
    return classes - 1, the_score


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
        return ' -m GTR+FO -pre 1class'
    else:
        repeats = ''
        for i in range(classes-1):
            repeats = repeats + ',GTR+FO'
        return ' -m "MIX{GTR+FO' + repeats + '}" -pre ' + str(classes) + 'classes'
    
def tree_file_name(tree_file_position,classes):
    ''' return right cmd statement for renaming the iqtree files with their numbers of classes'''
    if classes == 1:
        return tree_file_position + '/' + str(classes) +'class.iqtree'
    else:
        return tree_file_position + '/' + str(classes) +'classes.iqtree'

    






