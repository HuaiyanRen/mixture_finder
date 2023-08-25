import csv 
import numpy as np
import os

with open('result.csv','w+',newline='') as csvf:
    csv_write = csv.writer(csvf)
    csv_write.writerow(['name', 'classes',  'ntaxa', 'sites', 'invariable', 'rate', 'tree_length', 'min_w',
                        'optimal1', 'invar1', 'rate_o1',  'rate_re1', 'change1', 'tree_length1', 'nrf1', 'llh1', 'bic1', 'time1', 'ise_q1', 'ise_f1',
                        'optimal2', 'invar2', 'rate_o2',  'rate_re2', 'change2', 'tree_length2', 'nrf2', 'llh2', 'bic2', 'time2', 'ise_q2', 'ise_f2',
                        'tree_length0','nrf0', 'llh0', 'bic0', 'time0', 'ise_q0', 'ise_f0'])
                       
classes = [1,2,3,4,5] 
rates = [2] # 0: +E, 1: +I, 2: +I+G
length = [1000,2000,5000,10000]
ntaxa = [100]
replicates = list(np.arange(0,20,1))

classes_list = [m for m in classes for i in range(len(rates)*len(length)*len(ntaxa)*len(replicates))]
rates_list = [m for m in rates for i in range(len(length)*len(ntaxa)*len(replicates))]*len(classes)
length_list = [m for m in length for i in range(len(ntaxa)*len(replicates))]*len(classes)*len(rates)
ntaxa_list = [m for m in ntaxa for i in range(len(replicates))]*len(classes)*len(rates)*len(length)
replicates_list = replicates*len(classes)*len(rates)*len(length)*len(ntaxa)
    
tuple_list = [0]*len(classes)*len(rates)*len(length)*len(ntaxa)*len(replicates)

for i in range(len(tuple_list)):
    tuple_list[i] = classes_list[i], rates_list[i], length_list[i], ntaxa_list[i], replicates_list[i]

def get_para(line):
    gtr_all = line.split()[-1]
    gtr_list = gtr_all.split(',')
    
    if '+FO{' in gtr_all:
        A = float(gtr_list[-4].split('{')[-1])
        C = float(gtr_list[-3])
        G = float(gtr_list[-2])
        T = float(gtr_list[-1].split('}')[0])
    else:
        A,C,G,T = 0.25,0.25,0.25,0.25

    if 'JC' in line or 'F81' in line:
        AC,AG,AT,CG,CT = 1,1,1,1,1
    elif 'K2P' in line or 'HKY' in line or 'K80' in line:
        AC,AT,CG = 1,1,1
        AG = float(gtr_list[0].split('{')[1].split('}')[0])
        CT = AG
    elif 'TN' in line:
        AC,AT,CG = 1,1,1
        AG = float(gtr_list[0].split('{')[1])
        CT = float(gtr_list[1].split('}')[0])
    elif 'K3P' in line or 'K81' in line:
        AC = 1
        AG = float(gtr_list[0].split('{')[1])
        AT = float(gtr_list[1].split('}')[0])
        CG = AT
        CT = AG
    elif 'TPM2' in line:
        CG = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1].split('}')[0])
        AT = AC
        CT = AG
    elif 'TPM3' in line:  
        AT = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1].split('}')[0])
        CG = AC
        CT = AG
    elif 'TIM2' in line:
        CG = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        AT = AC
        CT = float(gtr_list[2].split('}')[0])
    elif 'TIM3' in line:
        AT = 1
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        CG = AC
        CT = float(gtr_list[2].split('}')[0])
    elif 'TIM' in line:
        AC = 1
        AG = float(gtr_list[0].split('{')[1])
        AT = float(gtr_list[1])
        CG = AT
        CT = float(gtr_list[2].split('}')[0])
    elif 'TVM' in line:
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        AT = float(gtr_list[2])
        CG = float(gtr_list[3].split('}')[0])
        CT = AG
    elif 'SYM' in line or 'GTR' in line:
        AC = float(gtr_list[0].split('{')[1])
        AG = float(gtr_list[1])
        AT = float(gtr_list[2])
        CG = float(gtr_list[3])
        CT = float(gtr_list[4].split('}')[0])
    return AC,AG,AT,CG,CT,A,C,G,T     
    
for paras in tuple_list:
    classes, rates, length, ntaxa, replicates = paras
    
    # dataset_name
    file_name = 'c' + str(classes) + '_r' + str(rates) + '_l' + str(length) + '_t' + str(ntaxa) + '_rep' + str(replicates)
    simu_file = file_name + '.treefile.log'
    invariable = 0
    with open(simu_file) as b:
        for line in b.readlines():
            if 'Proportion of invariable sites:' in line:
                invariable = float(line.split()[-1])
                
    true_treefile = file_name + '.treefile'
    true_tree = open(true_treefile,'r').read()
    tts = true_tree.split(':')
    true_tree_length = 0
    for i in range(1,len(tts)):
        true_tree_length = true_tree_length + float(tts[i].split(',')[0].split(')')[0])
    
    #ise
    true_q = [[0,0,0,0,0],
             [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
             [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
             [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
    true_f = [[0,0,0,0],
             [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
             [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
             [0,0,0,0],[0,0,0,0],[0,0,0,0]]
    true_w = [0]*int(classes)
    true_n = [[0,0,0,0,0],
             [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
             [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
             [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
    
    if classes == 1:
        with open(simu_file) as b:
            for line in b.readlines():
                if 'A-C:' in line:
                    true_q[1][0] = float(line.split()[-1])
                if 'A-G:' in line:
                    true_q[1][1] = float(line.split()[-1])
                if 'A-T:' in line:
                    true_q[1][2] = float(line.split()[-1])
                if 'C-G:' in line:
                    true_q[1][3] = float(line.split()[-1])
                if 'C-T:' in line:
                    true_q[1][4] = float(line.split()[-1])
                if 'pi(A)' in line:
                    true_f[1][0] = float(line.split()[-1])
                if 'pi(C)' in line:
                    true_f[1][1] = float(line.split()[-1])
                if 'pi(G)' in line:
                    true_f[1][2] = float(line.split()[-1])
                if 'pi(T)' in line:
                    true_f[1][3] = float(line.split()[-1])
        true_w[0] = 1 - invariable
                    
    elif classes > 1:
        with open(simu_file) as b:
            for line in b.readlines():
                for i in range(1,classes+1):
                    if str(i) + '  GTR' in line:
                        true_w[i-1] = float(line.split()[-2])
                        gtr_all = line.split()[-1]
                        gtr_list = gtr_all.split(',')
                        true_q[i][0] = float(gtr_list[0].split('{')[1])
                        true_q[i][1] = float(gtr_list[1])
                        true_q[i][2] = float(gtr_list[2])
                        true_q[i][3] = float(gtr_list[3])
                        true_q[i][4] = float(gtr_list[4].split('}')[0])
                        true_f[i][0] = float(gtr_list[4].split('{')[1])
                        true_f[i][1] = float(gtr_list[5])
                        true_f[i][2] = float(gtr_list[6])
                        true_f[i][3] = float(gtr_list[7].split('}')[0])

    for j in range(1, int(classes)+1):
        true_deno = sum(true_q[j]) + 1
        true_n[j][0] = true_q[j][0]/true_deno
        true_n[j][1] = true_q[j][1]/true_deno
        true_n[j][2] = true_q[j][2]/true_deno
        true_n[j][3] = true_q[j][3]/true_deno
        true_n[j][4] = true_q[j][4]/true_deno
    
    min_w = min(true_w) 
    
    result_row = [file_name, classes, ntaxa, length, invariable, 2, true_tree_length, min_w]

    
    if os.path.isfile('all/'+ file_name + '.iqtree'):
        invar = 0 
        with open('all/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'Best-fit model according to BIC:' in line:
                    optimal = line.count(',') + 1
                    model_re = line.split()[-1]
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
                if '(sum of branch lengths)' in line:
                    tree_length = line.split()[-1]
                if 'Log-likelihood of the tree:' in line:
                    llhm = float(line.split()[4])
                if 'Bayesian information criterion (BIC) score:'  in line:
                    bicm = float(line.split()[-1])
           
        with open('all/' + file_name + '.log') as b:
            for line in b.readlines():
                if 'chosen according to BIC' in line:
                    model_o = line.split()[2]
                    break
                
        if os.path.isfile('all/'+ file_name + '_rf.rfdist'):
            with open('all/' + file_name + '_rf.rfdist') as b:
                for line in b.readlines():
                    if 'Tree0' in line:
                        nrfm = float(line.split()[-1])/(2*ntaxa-6)
        else:
            nrfm = ''
                                        
        if '+R' in model_re:
            for k in range(2,11):
                if 'R'+str(k) in model_re:
                    rate_re = k+1
        elif '+G' in model_re:
            rate_re = 2
        elif '+I' in model_re:
            rate_re = 1    
        else:
            rate_re = 0
                
        if '+R' in model_o:
            for k in range(2,11):
                if 'R'+str(k) in model_o:
                    rate_o = k+1
        elif '+G' in model_o:
            rate_o = 2
        elif '+I' in model_o:
            rate_o = 1    
        else:
            rate_o = 0    
            
            
        change = 'equal'
        if rate_re < rate_o:
            change = 'small'
        elif rate_re > rate_o:
            change = 'large'
        
        #time
        m_time = [] 
        with open('all/' + file_name + '.log') as b:
            for line in b.readlines():
                if 'Time for fast ML tree search:' in line:
                    m_time.append(float(line.split()[6]))
                if 'CPU time for ModelFinder:' in line:
                    m_time.append(float(line.split()[4]))
                if 'CPU time used for tree search:' in line:
                    m_time.append(float(line.split()[6]))
        
        #ise
        with open('all/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'odel of substitution:' in line:
                    if optimal > 1:
                        est_q_list = line.split('{')[-1].split('}')[0].split(',')
                        for q in range(optimal):
                            if '+F' in est_q_list[q]:
                                est_q_list[q] = est_q_list[q].replace("+FO", "")

        est_q = [[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
        est_f = [[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0]]
        est_w = [0]*optimal

        est_n = [[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]

        if optimal == 1:
            est_f[1][0:4] = 0.25,0.25,0.25,0.25
            with open('all/' + file_name + '.iqtree') as b:
                for line in b.readlines():
                    if 'A-C:' in line:
                        est_q[1][0] = float(line.split()[-1])
                    if 'A-G:' in line:
                        est_q[1][1] = float(line.split()[-1])
                    if 'A-T:' in line:
                        est_q[1][2] = float(line.split()[-1])
                    if 'C-G:' in line:
                        est_q[1][3] = float(line.split()[-1])
                    if 'C-T:' in line:
                        est_q[1][4] = float(line.split()[-1])
                    if 'pi(A)' in line:
                        est_f[1][0] = float(line.split()[-1])
                    if 'pi(C)' in line:
                        est_f[1][1] = float(line.split()[-1])
                    if 'pi(G)' in line:
                        est_f[1][2] = float(line.split()[-1])
                    if 'pi(T)' in line:
                        est_f[1][3] = float(line.split()[-1])
            est_w[0] = 1 - invar                     
        else:
            with open('all/' + file_name + '.iqtree') as b:
                for line in b.readlines():
                    for i in range(1,optimal+1):
                        if str(i) + '  '+ est_q_list[i-1] in line:
                            est_w[i-1] = float(line.split()[-2])
                            est_q[i][0:5] = get_para(line)[0:5]
                            est_f[i][0:4] = get_para(line)[5:9]
            est_w = [n/sum(est_w) for n in est_w]
                    
        for j in range(1, optimal+1):
            est_deno = sum(est_q[j]) + 1
            est_n[j][0] = est_q[j][0]/est_deno
            est_n[j][1] = est_q[j][1]/est_deno
            est_n[j][2] = est_q[j][2]/est_deno
            est_n[j][3] = est_q[j][3]/est_deno
            est_n[j][4] = est_q[j][4]/est_deno                   
        addendq_1 = 0
        addendq_2 = 0
        addendq_3 = 0
        addendf_1 = 0
        addendf_2 = 0
        addendf_3 = 0
        for c1 in range(optimal):
            for c2 in range(optimal):
                    
                factorq = est_w[c1]*est_w[c2]
                for qi in range(5):
                    factorq = factorq*(1-max([est_n[c1+1][qi],est_n[c2+1][qi]]))
                addendq_1 = addendq_1 + factorq
                    
                factorf = est_w[c1]*est_w[c2]
                for fi in range(4):
                    factorf = factorf*(1-max([est_f[c1+1][fi],est_f[c2+1][fi]]))
                addendf_1 = addendf_1 + factorf
                   
        for c1 in range(optimal):
            for c2 in range(int(classes)):
                
                factorq = est_w[c1]*true_w[c2]
                for qi in range(5):
                    factorq = factorq*(1-max([est_n[c1+1][qi],true_n[c2+1][qi]]))             
                addendq_2 = addendq_2 + factorq
            
                factorf = est_w[c1]*true_w[c2]
                for fi in range(4):                    
                    factorf = factorf*(1-max([est_f[c1+1][fi],true_f[c2+1][fi]]))
                addendf_2 = addendf_2 + factorf
                    
        for c1 in range(int(classes)):
            for c2 in range(int(classes)):
                    
                factorq = true_w[c1]*true_w[c2]
                for qi in range(5):
                    factorq = factorq*(1-max([true_n[c1+1][qi],true_n[c2+1][qi]]))
                addendq_3 = addendq_3 + factorq
            
                factorf = true_w[c1]*true_w[c2]
                for fi in range(4):
                    factorf = factorf*(1-max([true_f[c1+1][fi],true_f[c2+1][fi]]))
                addendf_3 = addendf_3 + factorf             
        ise_q = addendq_1 - 2*addendq_2 + addendq_3
        ise_f = addendf_1 - 2*addendf_2 + addendf_3          
        
        result_row = result_row + [optimal, str(invar), str(rate_o), str(rate_re), change, tree_length, nrfm, llhm, bicm, sum(m_time), str(ise_q), str(ise_f)]
    else:
        #print(paras, 'method' + str(i))
        result_row = result_row + ['','','','','','','','','','','','']
        
    if os.path.isfile('gtr/'+ file_name + '.iqtree'):
        invar = 0 
        with open('gtr/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'Best-fit model according to BIC:' in line:
                    optimal = line.count(',') + 1
                    model_re = line.split()[-1]
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
                if '(sum of branch lengths)' in line:
                    tree_length = line.split()[-1]
                if 'Log-likelihood of the tree:' in line:
                    llhg = float(line.split()[4])
                if 'Bayesian information criterion (BIC) score:'  in line:
                    bicg = float(line.split()[-1])
                    
        with open('gtr/' + file_name + '.log') as b:
            for line in b.readlines():
                if 'chosen according to BIC' in line:
                    model_o = line.split()[2]
                    break
        
        if os.path.isfile('gtr/'+ file_name + '_rf.rfdist'):
            with open('gtr/' + file_name + '_rf.rfdist') as b:
                for line in b.readlines():
                    if 'Tree0' in line:
                        nrfg = float(line.split()[-1])/(2*ntaxa-6)
        else:
            nrfg = ''
                                        
        if '+R' in model_re:
            for k in range(2,11):
                if 'R'+str(k) in model_re:
                    rate_re = k+1
        elif '+G' in model_re:
            rate_re = 2
        elif '+I' in model_re:
            rate_re = 1    
        else:
            rate_re = 0
                
        if '+R' in model_o:
            for k in range(2,11):
                if 'R'+str(k) in model_o:
                    rate_o = k+1
        elif '+G' in model_o:
            rate_o = 2
        elif '+I' in model_o:
            rate_o = 1    
        else:
            rate_o = 0    
            
            
        change = 'equal'
        if rate_re < rate_o:
            change = 'small'
        elif rate_re > rate_o:
            change = 'large'
        
        #time
        g_time = [] 
        with open('gtr/' + file_name + '.log') as b:
            for line in b.readlines():
                if 'Time for fast ML tree search:' in line:
                    g_time.append(float(line.split()[6]))
                if 'CPU time for ModelFinder:' in line:
                    g_time.append(float(line.split()[4]))
                if 'CPU time used for tree search:' in line:
                    g_time.append(float(line.split()[6]))
                    
        #ise
        with open('gtr/'+ file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'odel of substitution:' in line:
                    if optimal > 1:
                        est_q_list = line.split('{')[-1].split('}')[0].split(',')
                        for q in range(optimal):
                            if '+F' in est_q_list[q]:
                                est_q_list[q] = est_q_list[q].replace("+FO", "")

        est_q = [[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
        est_f = [[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0]]
        est_w = [0]*optimal

        est_n = [[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]

        if optimal == 1:
            est_f[1][0:4] = 0.25,0.25,0.25,0.25
            with open('gtr/'+ file_name + '.iqtree') as b:
                for line in b.readlines():
                    if 'A-C:' in line:
                        est_q[1][0] = float(line.split()[-1])
                    if 'A-G:' in line:
                        est_q[1][1] = float(line.split()[-1])
                    if 'A-T:' in line:
                        est_q[1][2] = float(line.split()[-1])
                    if 'C-G:' in line:
                        est_q[1][3] = float(line.split()[-1])
                    if 'C-T:' in line:
                        est_q[1][4] = float(line.split()[-1])
                    if 'pi(A)' in line:
                        est_f[1][0] = float(line.split()[-1])
                    if 'pi(C)' in line:
                        est_f[1][1] = float(line.split()[-1])
                    if 'pi(G)' in line:
                        est_f[1][2] = float(line.split()[-1])
                    if 'pi(T)' in line:
                        est_f[1][3] = float(line.split()[-1])
            est_w[0] = 1 - invar  
                                
        else:
            with open('gtr/'+ file_name + '.iqtree') as b:
                for line in b.readlines():
                    for i in range(1,optimal+1):
                        if str(i) + '  '+ est_q_list[i-1] in line:
                            est_w[i-1] = float(line.split()[-2])
                            est_q[i][0:5] = get_para(line)[0:5]
                            est_f[i][0:4] = get_para(line)[5:9]
            est_w = [n/sum(est_w) for n in est_w]
                    
        for j in range(1, optimal+1):
            est_deno = sum(est_q[j]) + 1
            est_n[j][0] = est_q[j][0]/est_deno
            est_n[j][1] = est_q[j][1]/est_deno
            est_n[j][2] = est_q[j][2]/est_deno
            est_n[j][3] = est_q[j][3]/est_deno
            est_n[j][4] = est_q[j][4]/est_deno
                            
        addendq_1 = 0
        addendq_2 = 0
        addendq_3 = 0
        addendf_1 = 0
        addendf_2 = 0
        addendf_3 = 0
        for c1 in range(optimal):
            for c2 in range(optimal):
                    
                factorq = est_w[c1]*est_w[c2]
                for qi in range(5):
                    factorq = factorq*(1-max([est_n[c1+1][qi],est_n[c2+1][qi]]))
                addendq_1 = addendq_1 + factorq
                    
                factorf = est_w[c1]*est_w[c2]
                for fi in range(4):
                    factorf = factorf*(1-max([est_f[c1+1][fi],est_f[c2+1][fi]]))
                addendf_1 = addendf_1 + factorf
                   
        for c1 in range(optimal):
            for c2 in range(int(classes)):
                
                factorq = est_w[c1]*true_w[c2]
                #print('- ',factorq)
                for qi in range(5):
                    factorq = factorq*(1-max([est_n[c1+1][qi],true_n[c2+1][qi]]))
                    #print(factorq, max([est_q[c1+1][qi],true_q[c2+1][qi]]),' from ',[est_q[c1+1][qi],true_q[c2+1][qi]])                
                addendq_2 = addendq_2 + factorq
            
                factorf = est_w[c1]*true_w[c2]
                for fi in range(4):                    
                    factorf = factorf*(1-max([est_f[c1+1][fi],true_f[c2+1][fi]]))
                addendf_2 = addendf_2 + factorf
                    
        for c1 in range(int(classes)):
            for c2 in range(int(classes)):
                    
                factorq = true_w[c1]*true_w[c2]
                for qi in range(5):
                    factorq = factorq*(1-max([true_n[c1+1][qi],true_n[c2+1][qi]]))
                addendq_3 = addendq_3 + factorq
            
                factorf = true_w[c1]*true_w[c2]
                for fi in range(4):
                    factorf = factorf*(1-max([true_f[c1+1][fi],true_f[c2+1][fi]]))
                addendf_3 = addendf_3 + factorf
                     
        ise_q = addendq_1 - 2*addendq_2 + addendq_3
        ise_f = addendf_1 - 2*addendf_2 + addendf_3
            
        result_row = result_row + [optimal, str(invar), str(rate_o), str(rate_re), change, tree_length, nrfg, llhg, bicg, sum(g_time),str(ise_q),str(ise_f)]
    else:
        #print(paras, 'method' + str(i))
        result_row = result_row + ['','','','','','','','','','','',''] 
        
    if os.path.isfile('one/'+ file_name + '.iqtree'):
        invar = 0
        with open('one/' + file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'Proportion of invariable sites:' in line:
                    invar = float(line.split()[-1])
                if '(sum of branch lengths)' in line:
                    tree_length = line.split()[-1]
                if 'Log-likelihood of the tree:' in line:
                    llh = float(line.split()[4])
                if 'Bayesian information criterion (BIC) score:'  in line:
                    bic = float(line.split()[-1])
                    
        if os.path.isfile('one/'+ file_name + '_rf.rfdist'):
            with open('one/' + file_name + '_rf.rfdist') as b:
                for line in b.readlines():
                    if 'Tree0' in line:
                        nrf1 = float(line.split()[-1])/(2*ntaxa-6)   
        else:
            nrf1 = ''
        
        #time
        o_time = []
        with open('one/' + file_name + '.log') as b:
            for line in b.readlines():
                if 'Time for fast ML tree search:' in line:
                    o_time.append(float(line.split()[6]))
                if 'CPU time for ModelFinder:' in line:
                    o_time.append(float(line.split()[4]))
                if 'CPU time used for tree search:' in line:
                    o_time.append(float(line.split()[6]))
        
        #ise
        optimal = 1
        est_q = [[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
        est_f = [[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                 [0,0,0,0],[0,0,0,0],[0,0,0,0]]
        est_w = [0]*optimal

        est_n = [[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
        
        est_f[1][0:4] = 0.25,0.25,0.25,0.25
        with open('one/'+ file_name + '.iqtree') as b:
            for line in b.readlines():
                if 'A-C:' in line:
                    est_q[1][0] = float(line.split()[-1])
                if 'A-G:' in line:
                    est_q[1][1] = float(line.split()[-1])
                if 'A-T:' in line:
                    est_q[1][2] = float(line.split()[-1])
                if 'C-G:' in line:
                    est_q[1][3] = float(line.split()[-1])
                if 'C-T:' in line:
                    est_q[1][4] = float(line.split()[-1])
                if 'pi(A)' in line:
                    est_f[1][0] = float(line.split()[-1])
                if 'pi(C)' in line:
                    est_f[1][1] = float(line.split()[-1])
                if 'pi(G)' in line:
                    est_f[1][2] = float(line.split()[-1])
                if 'pi(T)' in line:
                    est_f[1][3] = float(line.split()[-1])
        est_w[0] = 1 - invar  
                
                    
        for j in range(1, 2):
            est_deno = sum(est_q[j]) + 1
            est_n[j][0] = est_q[j][0]/est_deno
            est_n[j][1] = est_q[j][1]/est_deno
            est_n[j][2] = est_q[j][2]/est_deno
            est_n[j][3] = est_q[j][3]/est_deno
            est_n[j][4] = est_q[j][4]/est_deno
                            
        addendq_1 = 0
        addendq_2 = 0
        addendq_3 = 0
        addendf_1 = 0
        addendf_2 = 0
        addendf_3 = 0
        for c1 in range(optimal):
            for c2 in range(optimal):
                    
                factorq = est_w[c1]*est_w[c2]
                for qi in range(5):
                    factorq = factorq*(1-max([est_n[c1+1][qi],est_n[c2+1][qi]]))
                addendq_1 = addendq_1 + factorq
                    
                factorf = est_w[c1]*est_w[c2]
                for fi in range(4):
                    factorf = factorf*(1-max([est_f[c1+1][fi],est_f[c2+1][fi]]))
                addendf_1 = addendf_1 + factorf
                   
        for c1 in range(optimal):
            for c2 in range(int(classes)):
                
                factorq = est_w[c1]*true_w[c2]
                #print('- ',factorq)
                for qi in range(5):
                    factorq = factorq*(1-max([est_n[c1+1][qi],true_n[c2+1][qi]]))
                    #print(factorq, max([est_q[c1+1][qi],true_q[c2+1][qi]]),' from ',[est_q[c1+1][qi],true_q[c2+1][qi]])                
                addendq_2 = addendq_2 + factorq
            
                factorf = est_w[c1]*true_w[c2]
                for fi in range(4):                    
                    factorf = factorf*(1-max([est_f[c1+1][fi],true_f[c2+1][fi]]))
                addendf_2 = addendf_2 + factorf
                    
        for c1 in range(int(classes)):
            for c2 in range(int(classes)):
                    
                factorq = true_w[c1]*true_w[c2]
                for qi in range(5):
                    factorq = factorq*(1-max([true_n[c1+1][qi],true_n[c2+1][qi]]))
                addendq_3 = addendq_3 + factorq
            
                factorf = true_w[c1]*true_w[c2]
                for fi in range(4):
                    factorf = factorf*(1-max([true_f[c1+1][fi],true_f[c2+1][fi]]))
                addendf_3 = addendf_3 + factorf
                     
        ise_q = addendq_1 - 2*addendq_2 + addendq_3
        ise_f = addendf_1 - 2*addendf_2 + addendf_3
            
        result_row = result_row + [tree_length,nrf1, llh, bic, sum(o_time),str(ise_q),str(ise_f)]
    else:
        result_row = result_row + ['','','','','','','']
    
    with open('result.csv','a+',newline='') as csvf:
        csv_write = csv.writer(csvf)
        csv_write.writerow(result_row)         
    