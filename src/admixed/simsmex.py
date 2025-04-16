import numpy as np
import msprime
import random
from numpy import linalg as LNG 
from random import randint, randrange
import pandas as pd


import HMMmex2 as hmm
N=5

#несколько вспомогательных функций
def connected(m):
    for i in range(len(m)-1):
        if m[i][1] == m[i+1][0]:
            return True
    return False
        
def remove_one(m):
    mas = m
    while connected(mas) == True:
        for i in range(len(mas)-1):
            if mas[i][1] == mas[i+1][0]:
                mas[i][1] = mas[i+1][1]
                mas.pop(i+1)
                break
    return mas


#Вход: ts, название популяции, индивид(которого мы препарируем), время предка
def get_migrating_tracts_ind(ts, pop_name, ind, T_anc):
    
    pop=-1
    for i in ts.populations():
        if i.metadata['name']==pop_name:
            pop=i.id
    
    mig = ts.tables.migrations
    migration_int = []

    for tree in ts.trees():  #перебираем все деревья. Как известно, каждому дереву отвечает участок днк  
        anc_node = ind #выбираем мексиканца
        while tree.time( tree.parent(anc_node) ) <= T_anc : #идем в прошлое до вершины anc_node по предкам нашего мексиканца, пока не наткнемся на миграцию 
            anc_node = tree.parent(anc_node)
        migs = np.where(mig.node == anc_node)[0] #выбирем все строки, соответствующие заданному узлу

        #идем по таблице миграций с anc_node и проверяем, чтобы миграции попадали в тот самый участок днк
        for i in migs:

            stroka = mig[i]
            if stroka.time == T_anc and stroka.dest == pop and tree.interval.left >= stroka.left and tree.interval.right <= stroka.right:
                migration_int.append([tree.interval.left, tree.interval.right])

    migration_int2 = []
    for i in range(len(migration_int)):
        if migration_int[i][0] != migration_int[i][1]:
            migration_int2.append(migration_int[i])
    migration_int = migration_int2
    

    
    mi = remove_one(migration_int)
    mi.sort()  

    return mi  




# Creates the sequence of observations.
def createSeqObsS(ts,cut,ind,pop) -> list:
    tables = ts.dump_tables()
    nodes = tables.nodes
    seq = np.zeros(int(ts.sequence_length/cut),dtype=int)  #The list which will contain the result of our function.
    pop_id = [p.id for p in ts.populations() if p.metadata['name']==pop][0] 
    start=-1
    end=-1
    for i in range(len(nodes)):
        if nodes[i].population==pop_id and start ==-1:
            start = i
        if nodes[i].population!=pop_id and start !=-1:
            end = i-1
            break
    
    for v in ts.variants():
        i=int(v.site.position/cut)
        c = False
        x=0
        while x < end: 
            if(nodes[x].individual==ind): # Check is this position corresponds to individual 'ind'
                b=v.genotypes[x]
                
                x=start
            elif(nodes[x].population==pop_id): # Check is this position corresponds to an individual in 'pop'
                
                if(v.genotypes[x]==b): 
                    c = True # The mutation found in 'ind' is also present in the population 'pop' 
                    x = end # We can skip all the remaining individuals
            x=x+1                     
        if not c : # If c is False it means that the mutation in 'ind' has not been found in the population 'pop' 
            seq[i]=seq[i]+1
    return seq

#функция, наблюдения для неандертальцев в зависимости от их количества
def createSeqObs_neand(ts,cut, ind, n_neanderthals, n_mexicans, n, n_neand ):
    #если в массиве есть повторяющиеся значения, то возвращает True, иначе False
    def repetitions(arr):

        for elem in arr:
            if arr.count(elem) > 1:
                return True
        return False
    inx = [randrange(3*n+n_mexicans, 3*n+n_mexicans+n_neand) for i in range(n_neanderthals)]
    
    while  repetitions(inx) is True:
        inx = [randrange(3*n+n_mexicans, 3*n+n_mexicans+n_neand) for i in range(n_neanderthals)]
        
    inx.sort()
#    print(inx)    
    seq = np.zeros(int(ts.sequence_length/cut),dtype=int)  #The list which will contain the result of our function.
    
    for v in ts.variants():
        i=int(v.site.position/cut)
        
        flag = False
        for j in inx:
            if abs(v.genotypes[j]-v.genotypes[ind])==0:
                flag = True
                break
                
            
        if flag == False:
            seq[i] += 1
              
    return seq

#функция, наблюдения
def createSeqObs_main(ts,cut,pop, ind,  n_mexicans, n, n_ref_pop ):
    
    inx = []
    if pop == 'EU':
        inx=[i for i in range(n_mexicans, n_mexicans+n_ref_pop)]
    if pop == 'NA':
        inx = [i for i in range(n_mexicans+n, n+n_mexicans +n_ref_pop)]
    if pop == 'AF':
        inx = [i for i in range(n_mexicans+2*n, 2*n+n_mexicans+n_ref_pop)]        

    seq = np.zeros(int(ts.sequence_length/cut),dtype=int)    
    for v in ts.variants():
        i=int(v.site.position/cut)
        
        flag = False
        for j in inx:
            if abs(v.genotypes[j]-v.genotypes[ind])==0 and v.genotypes[ind]==1:
                flag = True
                break
                
            
        if flag == False:
            seq[i] += 1
              
    return seq


#функция, наблюдения
def createSeqObs_main2(ts,cut, ind,  n_mexicans, n, n_ref_pop, n_neand ):
    c=0


    inx0 = [i for i in range(n_mexicans, n_mexicans+n_ref_pop)]
    inx1 = [i for i in range(n_mexicans+n, n+n_mexicans +n_ref_pop)]
    inx2 = [i for i in range(n_mexicans+2*n, 2*n+n_mexicans+n_ref_pop)] 
    inx3 = [i for i in range(3*n+n_mexicans, 3*n+n_mexicans+n_neand)]

    seq = np.zeros((int(ts.sequence_length/cut), 4), dtype=int )    
    for v in ts.variants():
        i=int(v.site.position/cut)
        
        flag = [False, False, False, False]
        for j in inx0:
            if abs(v.genotypes[j]-v.genotypes[ind])==0:
                flag[0] = True
                break
                
        for j in inx1:
            if abs(v.genotypes[j]-v.genotypes[ind])==0:
                flag[1] = True
                break
  
        for j in inx2:
            if abs(v.genotypes[j]-v.genotypes[ind])==0:
                flag[2] = True
                break
                
        for j in inx3:
            if abs(v.genotypes[j]-v.genotypes[ind])==0:
                flag[3] = True
                break
    
            
            
        for j in range(0,len(flag)):
            if flag[j]==False:
                seq[i][j]+=1
                
#        if flag[0]==False and flag[1]==False:
#            seq[i][0]=seq[i][0]-1
#            seq[i][1]=seq[i][1]-1
#            c+=1            
#    print(c)          
    return seq


def color_states(cut,tractsND, tractsEU, tractsAS, tractsAF, seq_length):
    
    
    # лежит ли точка p в интервале i?
    def point_in_interval(p, interval):
        if p >= interval[0] and p <= interval[1]:
            return True
        else:
            return False
    
    state_seq = np.ones(int(seq_length/cut))
    state_seq = -state_seq
    for t in range(0, int(seq_length), cut):



        f=False
        for j in tractsAF:
            if point_in_interval(t, j)==True:
                f= True
        if f==True:
            state_seq[int(t/cut)]=4

        f=False
        for j in tractsEU:
            if point_in_interval(t, j)==True:
                f= True
        if f==True :
            state_seq[int(t/cut)]=0



        f=False
        for j in tractsAS:
            if point_in_interval(t, j)==True:
                f= True
        if f==True :
            state_seq[int(t/cut)]=2

        f=False
        for j in tractsND:
            if point_in_interval(t, j)==True:
                f= True
        if f==True and state_seq[int(t/cut)]==0:
            state_seq[int(t/cut)]=1
        if f==True and state_seq[int(t/cut)]==2:
            state_seq[int(t/cut)]=3
    return state_seq




def intersections(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]
        ri += 1
    return ranges


def len_tracts(set_intervals):
    if len(set_intervals)==0:
        return 0
    else:
        s=0
        for j in range(len(set_intervals)):
            s+= set_intervals[j][1]-set_intervals[j][0]+1
        return s
    
def confusion_mtrx(real, res_HMM):
    N=5
    conf_matrix = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            conf_matrix[i,j] =int(len_tracts(intersections(real[i], res_HMM[j])))
            

    return conf_matrix

def classification_rpt(conf_matrix):
    clas_report = {}
    for i in range(N):
        dd={}
        dd['precision'] = round(conf_matrix[i,i]/sum(conf_matrix[:,i]),7)
        dd['recall'] = round(conf_matrix[i,i]/sum(conf_matrix[i,:]),7)
        dd['f1-score'] = round(2*dd['recall']*dd['precision']/(dd['recall']+dd['precision']), 7)
        clas_report[str(i)] = dd
    return clas_report


def colored_intervals(tracts_nd, tracts_eu, tracts_na, tracts_af):
    col_tracts=[]
    
    modern_eu = []


    for tr in tracts_eu:
        tr = [tr]
        simple_int = intersections(tr, tracts_nd)
        if simple_int != []:
            if simple_int[0][0]> tr[0][0]:
                modern_eu.append([tr[0][0], simple_int[0][0]])
            else: 
                modern_eu.append([tr[0][0], simple_int[0][1] ])
            for j in range(len(simple_int)-1):
                modern_eu.append([simple_int[j][1], simple_int[j+1][0]])
            if simple_int[-1][1]< tr[0][1]:
                modern_eu.append([ simple_int[-1][1], tr[0][1]])
        else:
            modern_eu.append(tr[0])
            
    modern_na = []
    for tr in tracts_na:
        tr = [tr]
        simple_int = intersections(tr, tracts_nd)
        if simple_int != []:
            if simple_int[0][0]> tr[0][0]:
                modern_na.append([tr[0][0], simple_int[0][0]])
            else: 
                modern_na.append([tr[0][0], simple_int[0][1] ])
            for j in range(len(simple_int)-1):
                modern_na.append([simple_int[j][1], simple_int[j+1][0]])
            if simple_int[-1][1]< tr[0][1]:
                modern_na.append([ simple_int[-1][1], tr[0][1]])
        else:
            modern_na.append(tr[0])
            
    return [modern_eu, intersections(tracts_nd, tracts_eu), modern_na, intersections(tracts_nd, tracts_na), tracts_af]





def colored_intervals2(tracts_nd, tracts_eu, tracts_na, tracts_af):
    col_tracts=[]
    
    modern_eu = []


    for tr in tracts_eu:
        tr = [tr]
        simple_int = intersections(tr, tracts_nd)
        if simple_int != []:
            if simple_int[0][0]> tr[0][0]:
                modern_eu.append([tr[0][0], simple_int[0][0]])
            for j in range(len(simple_int)-1):
                modern_eu.append([simple_int[j][1], simple_int[j+1][0]])
            if simple_int[-1][1]< tr[0][1]:
                modern_eu.append([ simple_int[-1][1], tr[0][1]])
        else:
            modern_eu.append(tr[0])
            
    modern_na = []
    for tr in tracts_na:
        tr = [tr]
        simple_int = intersections(tr, tracts_nd)
        if simple_int != []:
            if simple_int[0][0]> tr[0][0]:
                modern_na.append([tr[0][0], simple_int[0][0]])
#            else: 
 #               modern_na.append([ simple_int[0][1], simple_int[1][0] ])
            for j in range(len(simple_int)-1):
                modern_na.append([simple_int[j][1], simple_int[j+1][0]])
            if simple_int[-1][1]< tr[0][1]:
                modern_na.append([ simple_int[-1][1], tr[0][1]])
        else:
            modern_na.append(tr[0])
            
    return [modern_eu, intersections(tracts_nd, tracts_eu), modern_na, intersections(tracts_nd, tracts_na), tracts_af]


def createDataFrameMex(seq_length,rr,mu,cut,Lmbd_opt, n_st, seq, n_neanderthal,  
                       n_ref_pop,  nd_tracts,  eu_tracts, na_tracts, af_tracts, idx_mas):
    
    

    df= pd.DataFrame(columns=['State', 'Value', 'Score', 'n_eu',
                                       'n_neand', 'L',  'n_ref_pop'])    
    d =  mu * cut  

    P=[0.4, 0.05, 0.4, 0.05, 0.1]
    iii=0

    for idx in idx_mas:



        lmbd_opt=Lmbd_opt[iii]
        
        a = hmm.initA(lmbd_opt[5]/d, lmbd_opt[6]/d, rr, cut, lmbd_opt[7],  lmbd_opt[8],  lmbd_opt[9],  lmbd_opt[10])
        b=hmm.initB(mu, cut, lmbd_opt[0:5], n_st) 
        tracts_HMM =  hmm.get_HMM_tracts(hmm.viterbi(seq [idx], P, a, b))
        for k in range(N):
            for j in range(len(tracts_HMM[k])):
                tracts_HMM[k][j][0]= cut * tracts_HMM[k][j][0]
                tracts_HMM[k][j][1]= cut * (tracts_HMM[k][j][1]+1)-1

        iii+=1
 
        real_tracts = colored_intervals2(nd_tracts[idx], eu_tracts[idx], na_tracts[idx], af_tracts[idx])

        cl_report = classification_rpt(confusion_mtrx(real_tracts, tracts_HMM))
     
        print(tracts_HMM[4])

            
        for j in range(N):
            df.loc[len(df.index)] = [j, cl_report[str(j)]['precision'], 'precision',idx, n_neanderthal, cut,
                                         n_ref_pop]
            df.loc[len(df.index)] = [j, cl_report[str(j)]['recall'], 'recall',idx, n_neanderthal, cut, 
                                         n_ref_pop]

            
    return df

def createDataFrameMex_main(tracts_HMM, cut, n_neanderthal,  
                       n_ref_pop,  nd_tracts,  eu_tracts, na_tracts, af_tracts, idx_mas):

    df= pd.DataFrame(columns=['State', 'Value', 'Score', 'n_eu',
                                       'n_neand', 'L',  'n_ref_pop'])    


    for idx in idx_mas:
        real_tracts = colored_intervals2(nd_tracts[idx], eu_tracts[idx], na_tracts[idx], af_tracts[idx])

        cl_report = classification_rpt(confusion_mtrx(real_tracts, [tracts_HMM[0][idx], tracts_HMM[1][idx], 
                                                                    tracts_HMM[2][idx], tracts_HMM[3][idx],
                                                                    tracts_HMM[4][idx]]))
            
        for j in range(N):
            df.loc[len(df.index)] = [j, cl_report[str(j)]['precision'], 'precision',idx, n_neanderthal, cut,
                                         n_ref_pop]
            df.loc[len(df.index)] = [j, cl_report[str(j)]['recall'], 'recall',idx, n_neanderthal, cut, 
                                         n_ref_pop]

            
    return df

