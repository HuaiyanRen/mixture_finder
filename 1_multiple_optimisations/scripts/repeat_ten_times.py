import os

def model_select(classes):
    ''' return right cmd statement for fitting mixture models of correct numbers classes'''
    if classes == 1:
        return ' -m GTR+FO -pre 1class -redo'
    else:
        repeats = ''
        for i in range(classes-1):
            repeats = repeats + ',GTR+FO'
        return ' -m "MIX{GTR+FO' + repeats + '}" -pre ' + str(classes) + 'classes -redo'

def find_score(filename,score_type):
    '''return AIC or BIC score from iqtree file'''
    with open(filename) as iqtree_result:
        for line in iqtree_result.readlines():
            if score_type in line:
                list_line = line.split()
                score = list_line[-1]
    return score

def tree_file_name(file_position,classes):
    ''' return right cmd statement for renaming the iqtree files with their numbers of classes'''
    if classes == 1:
        return file_position.split()[-1] + '/' + str(classes) +'class.iqtree'
    else:
        return file_position.split()[-1] + '/' + str(classes) +'classes.iqtree'




for classes in range(3,10):
    
    cmd = 'd: && cd D:\masterresearch\congnato && D:/masterresearch/iqtree2/iqtree2 -s alignment.nex' + model_select(classes)
    
    tree_file = tree_file_name('d: && cd D:\masterresearch\congnato',classes)
    
    for i in range(10):
        
        score_list = []

        os.system(cmd)
    
        new_score = float(find_score(tree_file, 'BIC'))
       
        score_list.append(new_score)
        
        print(classes, new_score)