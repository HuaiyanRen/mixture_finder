import csv 

import numpy as np
from scipy.stats import powerlognorm,beta,invweibull,alpha,genlogistic,exponweib,exponnorm,norm, gengamma,dweibull,laplace,loglaplace,loggamma,genexpon,exponpow,maxwell,powernorm,dgamma,chi,weibull_min,invgamma,lognorm,gennorm,invgauss,weibull_max,genpareto
#from scipy import stats

import pandas as pd
import ast

with open('modelpara.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes', 'rates', 'length', 'taxa', 'replicates',
                        'invar', 'gamma',
                        'q1','weight1','AC1', 'AG1', 'AT1', 'CG1', 'CT1', 'A1', 'C1', 'G1', 'T1',
                        'q2','weight2','AC2', 'AG2', 'AT2', 'CG2', 'CT2', 'A2', 'C2', 'G2', 'T2',
                        'q3','weight3','AC3', 'AG3', 'AT3', 'CG3', 'CT3', 'A3', 'C3', 'G3', 'T3',
                        'q4','weight4','AC4', 'AG4', 'AT4', 'CG4', 'CT4', 'A4', 'C4', 'G4', 'T4',
                        'q5','weight5','AC5', 'AG5', 'AT5', 'CG5', 'CT5', 'A5', 'C5', 'G5', 'T5',
                        'end'])

 
    
classes = [1] 
rates = [2] # 0: +E, 1: +I, 2: +I+G
length = [1000,2000,5000,10000]
ntaxa = [100]
replicates = list(np.arange(0,250,1))

classes_list = [m for m in classes for i in range(len(rates)*len(length)*len(ntaxa)*len(replicates))]
rates_list = [m for m in rates for i in range(len(length)*len(ntaxa)*len(replicates))]*len(classes)
length_list = [m for m in length for i in range(len(ntaxa)*len(replicates))]*len(classes)*len(rates)
ntaxa_list = [m for m in ntaxa for i in range(len(replicates))]*len(classes)*len(rates)*len(length)
replicates_list = replicates*len(classes)*len(rates)*len(length)*len(ntaxa)
    
tuple_list = [0]*len(classes)*len(rates)*len(length)*len(ntaxa)*len(replicates)

for i in range(len(tuple_list)):
    tuple_list[i] = classes_list[i], rates_list[i], length_list[i], ntaxa_list[i], replicates_list[i]

with open('branch_length.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name']+['ext' + str(i) for i in range(1, ntaxa[0]+1)] + ['check'] + ['int' + str(i) for i in range(1, ntaxa[0]-2)] + ['end'])   
     
#instance = np.random.RandomState(seed = 20230803) 
np.random.seed(20230803)

modeldist = pd.read_csv('model_para_dst.csv', index_col=None)

def random_q(classes):
    q_list = list(np.arange(0,18,1))
    times_list = [20,2531,893,11,1265,766,567,1370,548,1687,390,615,1109,1092,2939,2360,1873,2288]
    prob_list = [i/sum(times_list) for i in times_list]
    return np.random.choice(q_list, p = prob_list, size=classes)


def q_parameter(q_n):
    AC = dist(modeldist.AC_dst_name[q_n], modeldist.AC_dst_paras[q_n])
    AG = dist(modeldist.AG_dst_name[q_n], modeldist.AG_dst_paras[q_n])
    AT = dist(modeldist.AT_dst_name[q_n], modeldist.AT_dst_paras[q_n], a_c = AC)
    CG = dist(modeldist.CG_dst_name[q_n], modeldist.CG_dst_paras[q_n], a_c = AC, a_t =  AT)
    CT = dist(modeldist.CT_dst_name[q_n], modeldist.CT_dst_paras[q_n], a_g = AG)
    Ap = dist(modeldist.FA_dst_name[q_n], modeldist.FA_dst_paras[q_n])
    Cp = dist(modeldist.FC_dst_name[q_n], modeldist.FC_dst_paras[q_n])
    Gp = dist(modeldist.FG_dst_name[q_n], modeldist.FG_dst_paras[q_n])
    Tp = dist(modeldist.FT_dst_name[q_n], modeldist.FT_dst_paras[q_n])
    Fsum = Ap + Cp + Gp + Tp
    A, C, G, T = Ap/Fsum, Cp/Fsum, Gp/Fsum, Tp/Fsum
    return AC,AG,AT,CG,CT,A,C,G,T 

q_name ={0:'F81+F',
         1:'GTR+F',
         2:'HKY+F',
         3:'JC',
         4:'K2P',
         5:'K3P',
         6:'K3Pu+F',
         7:'SYM',
         8:'TIM2e',
         9:'TIM3e',
         10:'TIMe',
         11:'TNe',
         12:'TPM2',
         13:'TPM2u+F',
         14:'TPM3',
         15:'TPM3u+F',
         16:'TVM+F',
         17:'TVMe'}

def dist(name, paras, a_c = None, a_t = None, a_g = None):
    if name == '1':
        return 1
    elif name == '0.25':
        return 0.25
    elif name == 'AC':
        return a_c
    elif name == 'AT':
        return a_t
    elif name == 'AG':
        return a_g
    elif name == 'invweibull':
        return abs(invweibull.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'loglaplace':
        return abs(loglaplace.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'loggamma':
        return abs(loggamma.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'exponnorm':
        return abs(exponnorm.rvs(*ast.literal_eval(paras)).round(5))    
    elif name == 'dgamma':
        return abs(dgamma.rvs(*ast.literal_eval(paras)).round(5)) 
    elif name == 'dweibull':
        return abs(dweibull.rvs(*ast.literal_eval(paras)).round(5)) 
    elif name == 'alpha':
        return abs(alpha.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'norm':
        return abs(norm.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'gengamma':
        return abs(gengamma.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'laplace':
        return abs(laplace.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'chi':
        return abs(chi.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'beta':
        return abs(beta.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'weibull_min':
        return abs(weibull_min.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'invgamma':
        return abs(invgamma.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'genlogistic':
        return abs(genlogistic.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'lognorm':
        return abs(lognorm.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'exponweib':
        return abs(exponweib.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'maxwell':
        return abs(maxwell.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'genlogistic':
        return abs(genlogistic.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'exponpow':
        return abs(exponpow.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'gennorm':
        return abs(gennorm.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'genexpon':
        return abs(genexpon.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'powernorm':
        return abs(powernorm.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'powerlognorm':
        return abs(powerlognorm.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'invgauss':
        return abs(invgauss.rvs(*ast.literal_eval(paras)).round(5))
    elif name == 'weibull_max':
        return abs(weibull_max.rvs(*ast.literal_eval(paras)).round(5))
    else:
        print(name)

for r in range(len(tuple_list)):
    classes, rates, length, ntaxa, replicates = tuple_list[r]
    # dataset_name
    file_name = 'c' + str(classes) + '_r' + str(rates) + '_l' + str(length) + '_t' + str(ntaxa) + '_rep' + str(replicates)
    para_row = [file_name,classes,rates,length,ntaxa,replicates]+['']*58
    # gtr parameters
    if classes == 1:
        q_list = random_q(classes)
        AC,AG,AT,CG,CT,A,C,G,T = q_parameter(q_list[0])
        para_row[8:19] =  q_name[q_list[0]],1,AC,AG,AT,CG,CT,A,C,G,T
    else:        
        # class weights
        weight_list = np.random.uniform(0.1, 1, classes)
        weights = list(weight_list/sum(weight_list))
        # q matrices
        for i in range(classes):
            q_list = random_q(classes)
            AC,AG,AT,CG,CT,A,C,G,T = q_parameter(q_list[0])
            para_row[i*11+8:i*11+19] = q_name[q_list[i]],weights[i],AC,AG,AT,CG,CT,A,C,G,T
        
    # rate parameters
    if rates == 0:
        rate_para = ''
    elif rates == 1:
        invar = genpareto.rvs(*(-0.425, -1.346e-10, 0.345)).round(5)
        while invar >=1 or invar < 0:            
            invar = genpareto.rvs(*(-0.425, -1.346e-10, 0.345)).round(5)
        rate_para = '+I{' + str(invar) + '}'
        para_row[6] = invar
    elif rates == 2:
        invar = genpareto.rvs(*(-0.425, -1.346e-10, 0.345)).round(5)     
        while invar >=1 or invar < 0:            
            invar = genpareto.rvs(*(-0.425, -1.346e-10, 0.345)).round(5)
        gamma = invweibull.rvs(*(5.576, -1.346, 2.408)).round(5)
        rate_para = '+I{' + str(invar) + '}+G{' + str(gamma) + '}'
        para_row[6:8] = invar, gamma
        
    #record model parameters
    with open('modelpara.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(para_row) 

    # branch
    int_n = ntaxa - 3
    ext_n = ntaxa 
    
    ext_bran = abs(alpha.rvs(*(3.091e-11, -0.006, 0.012), size=ext_n).round(5))
    int_bran = abs(alpha.rvs(*(1.136e-09, -0.003, 0.007), size=int_n).round(5))
    branch_row = [file_name]+ list(ext_bran) +['check']+ list(int_bran) + ['']
    
    #record branch lengths
    with open('branch_length.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(branch_row) 
