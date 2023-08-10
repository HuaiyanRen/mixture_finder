import os
import argparse
#import numpy as np
from scipy.stats import chi2

type_dict_simple = {0:'JC',
                    1:'F81',
                    2:'K2P',
                    3:'HKY',
                    4:'TN',
                    5:'TNe',
                    6:'K3P',
                    7:'K3Pu',
                    8:'TPM2',
                    9:'TPM2u',
                    10:'TPM3',
                    11:'TPM3u',
                    12:'TIM',
                    13:'TIMe',
                    14:'TIM2',
                    15:'TIM2e',
                    16:'TIM3',
                    17:'TIM3e',
                    18:'TVM',
                    19:'TVMe',
                    20:'SYM',
                    21:'GTR'}

type_dict = {0:'JC+FQ',
             1:'F81+FO',
             2:'K2P+FQ',
             3:'HKY+FO',
             4:'TN+FO',
             5:'TNe+FQ',
             6:'K3P+FQ',
             7:'K3Pu+FO',
             8:'TPM2+FQ',
             9:'TPM2u+FO',
             10:'TPM3+FQ',
             11:'TPM3u+FO',
             12:'TIM+FO',
             13:'TIMe+FQ',
             14:'TIM2+FO',
             15:'TIM2e+FQ',
             16:'TIM3+FO',
             17:'TIM3e+FQ',
             18:'TVM+FO',
             19:'TVMe+FQ',
             20:'SYM+FQ',
             21:'GTR+FO'}

rate_dict = {0:'',
             1:'+I',
             2:'+G',
             3:'+I+G',
             4:'+R2',
             5:'+R3',
             6:'+R4',
             7:'+R5',
             8:'+R6',
             9:'+R7',
             10:'+R8',
             11:'+I+R2',
             12:'+I+R3',
             13:'+I+R4',
             14:'+I+R5',
             15:'+I+R6',
             16:'+I+R7',
             17:'+I+R8'}


def mixfinder(iqtree_loc, file_name, method, candi_q, nt, pre, bb):
    
    # create a folder for results
    folder_name = file_name + '.' + pre
    if not os.path.isfile(folder_name):
        os.system('mkdir ' + folder_name)
    if not os.path.isfile(folder_name + '/logs'):    
        os.system('mkdir ' + folder_name + '/logs')
    if not os.path.isfile(folder_name + '/result.txt'):    
        with open(folder_name + '/result.txt', 'w+') as result:
            result.write('mixfinder log \n')
        
    # initialization
    classes = 1
    if_adding = True
    model_num_seq = [0]
    if candi_q == 'mf':
        mset = ''
        type_list = list(type_dict)
    else:
        mset = ' -mset '+ type_dict_simple[int(candi_q)]
        type_list = list([int(candi_q)])
        
    cmd = iqtree_loc + 'iqtree2 -s ' + file_name + ' -m MF'+mset+' -mfreq FO -mrate E,I,G,I+G,R,I+R -pre ' + folder_name + '/logs/c1 -nt ' +nt
    os.system(cmd)
    
    tree_file = folder_name + '/logs/c1.iqtree'

    with open(tree_file) as iqtree_result:
        for line in iqtree_result.readlines():
            if 'Model of substitution:' in line:
                MF_result = line.split()[-1]
            if 'Log-likelihood of the tree:' in line:
                llh_old = float(line.split()[4])
            if 'Number of free parameters (#branches + #model parameters):' in line:
                df_old = float(line.split()[-1])
            if 'Bayesian information criterion (BIC) score:'  in line:
                bic_old = float(line.split()[-1])

    for i in range(len(type_dict_simple)):
        if type_dict_simple[i] in MF_result:
            model_num_seq[0] = i
    for i in range(len(rate_dict)):    
        if rate_dict[i] in MF_result:
            rate_num = i

    tree1 = folder_name + '/logs/c1.treefile' 
    with open(folder_name + '/result.txt', 'a+') as result:   
        print('c1: '+ MF_result + ', df: '+str(df_old) +', LnL: '+str(llh_old) + ', BIC: '+str(bic_old) + '\n', file = result)
        
    # adding classes
    while if_adding:
        classes = classes+1
        score_list = []
        llh_list = []
        
        for i in range(len(type_list)):
            log_file = folder_name + '/logs/c' + str(classes) + type_dict_simple[type_list[i]]
            model_num_seq_loop = model_num_seq + [type_list[i]]
            
            cmd = iqtree_loc + 'iqtree2 -s ' + file_name + model_select(model_num_seq_loop) + rate_dict[rate_num] + ' -te ' + tree1 + ' -pre ' + log_file + ' -nt ' +nt
            os.system(cmd)
            
            iqtree_file = log_file + '.iqtree'
            with open(iqtree_file) as iqtree_result:
                for line in iqtree_result.readlines():
                    if '(BIC)' in line:
                        list_line = line.split()
                        score = float(list_line[-1])
                    if 'Log-likelihood of the tree:' in line:
                        llh = float(line.split()[4])
            score_list.append(score)
            llh_list.append(llh)
        
        sorted_score_list = sorted(zip(score_list,type_list,llh_list))
        sorted_model = [x[1] for x in sorted_score_list]
        sorted_score = [x[0] for x in sorted_score_list] 
        sorted_llh = [x[2] for x in sorted_score_list]
                
        for i in range(len(sorted_score)):
            with open(folder_name + '/result.txt', 'a+') as result:
                print('-MF'+ str(classes) +': '+ type_dict_simple[sorted_model[i]]+' LnL: '+str(sorted_llh[i]) + ' BIC: '+str(sorted_score[i]), file = result) 
        
        with open( folder_name + '/logs/c' + str(classes) + type_dict_simple[sorted_model[0]] + '.iqtree') as iqtree_result:
            for line in iqtree_result.readlines():
                if 'Model of substitution:' in line:
                    MF_result = line.split()[-1]
                if 'Log-likelihood of the tree:' in line:
                    llh_new = float(line.split()[4])
                if 'Number of free parameters (#branches + #model parameters):' in line:
                    df_new = float(line.split()[-1])
                if 'Bayesian information criterion (BIC) score:'  in line:
                    bic_new = float(line.split()[-1])
        
        with open(folder_name + '/result.txt', 'a+') as result:   
            print('c'+str(classes)+': '+ print_model_names(model_num_seq+[sorted_model[0]]) + ', df: '+str(df_new) +', LnL: '+str(llh_new) + ', BIC: '+str(bic_new) + '\n', file = result)        

        if int(method) == 2:
            p_value = chi2.cdf(2*(llh_new - llh_old), df_new - df_old)
            if p_value > 0.95:
                if_adding = True
                #with open(folder_name + '.results/result.txt', 'a+') as result:
                #    print('2*('+ str(llh_new)+' - '+str(llh_old)+'), '+str(df_new)+' - '+str(df_old)+') p:' + str(p_value),file = result)
            else:
                if_adding = False
        elif int(method) == 3:
            if bic_new < bic_old:
                if_adding = True
            else:
                if_adding = False
        if if_adding:
            llh_old = llh_new 
            df_old = df_new
            bic_old = bic_new
            model_num_seq.append(sorted_model[0])

    #re-estimating
    if classes > 2:
        para_m = ['']*(classes-1)
        para_i = ''
        para_g = ''
        para_r = ''
    
        para_file = folder_name + '/logs/c' + str(classes-1) +  type_dict_simple[model_num_seq[-1]] + '.iqtree'
        with open(para_file) as b:
            for line in b.readlines():
                for i in range(classes-1):
                    if str(i+1) + '  ' + type_dict_simple[model_num_seq[i]] in line:
                        if model_num_seq[i] == 0:
                            para_m[i] = line.split()[-1] + '+FQ'
                        else:
                            para_m[i] = line.split()[-1].replace(',','/') 
                if 'Proportion of invariable sites:' in line:
                    para_i = '+I{' + line.split()[-1] + '}'
                if 'Gamma shape alpha:' in line:
                    para_g = '+G4{' + line.split()[-1] + '}'
                if 'Site proportion and rates:' in line:
                    line_list = line.split()
                    para_r_list = []
                    para_r_num = 0
                    for j in range(4,len(line_list)):
                        para_r_list.append(line_list[j].split('(')[-1].split(')')[0])
                        para_r_num = para_r_num + 1
                    para_r = '+R' + str(para_r_num) + '{' + ",".join(para_r_list) + '}'
    
        model_para = '"MIX{' + para_m[0] 
        for i in range(1, len(para_m)):
            model_para = model_para +',' + para_m[i]
        model_para = model_para + '}'
        if len(para_i) > 0:
            model_para = model_para + para_i
        if len(para_g) > 0:
            model_para = model_para + para_g
        if len(para_r) > 0:
            model_para = model_para + para_r
        model_para = model_para + '"'
        
        if bb == '':
            bb_str = ''
        else:
            bb_str = ' -bb '+ bb
        cmd = iqtree_loc + 'iqtree2 -s ' + file_name + ' -m ' + model_para + ' -pre ' + folder_name + '/logs/reest' + bb_str + ' -nt ' + nt
        os.system(cmd)
            
        with open(folder_name + '/result.txt', 'a+') as result:
            print('optimal is ' + str(classes-1) + ' classes model: ' + model_para, file = result)
    else:
        with open(folder_name + '/result.txt', 'a+') as result:
            print('optimal is ' + str(classes-1) + ' classes model', file = result)
    

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
    return '{' + ', '.join(model_type_seq) + '}'

parser = argparse.ArgumentParser(description='')
parser.add_argument('--iqtree_loc', '-loc', help='',
                    default = '')
parser.add_argument('--file_name', '-s', help='',
                    required=True)
parser.add_argument('--method', '-med', help='',
                    required=True)
parser.add_argument('--candi_q', '-m', help='',
                    default = 'mf')
parser.add_argument('--nt', '-nt', help='',
                    default = '1')
parser.add_argument('--pre', '-pre', help='',
                    default = 'results')
parser.add_argument('--bb', '-bb', help='',
                    default = '')
args = parser.parse_args()


if __name__ == '__main__':
    try:
        mixfinder(args.iqtree_loc,args.file_name, args.method, args.candi_q, args.nt, args.pre, args.bb)
    except Exception as e:
        print(e)









               
