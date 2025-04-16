import argparse
import argparse
import numpy as np
import sys
import useful as usfl
import HMM
import EM 


import pathlib




parser = argparse.ArgumentParser(description='DAIseg') 
parser.add_argument('--bed', type=str, help='Region bed file')
parser.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')
parser.add_argument('--EM_steps', type=str, help='number of EMsteps')
parser.add_argument('--EM_samples', type=int, help='Number of samples used in training')
parser.add_argument('--HMM_par', type= str, help='File with parameters')
parser.add_argument('--o', type= str, help = 'Name of output file' )
parser.add_argument('--EM_est', type= str, help = 'Make estimation of the all parameters or only coalescent times' )
parser.add_argument('--prepared_file', type=str, help='dlkfjgljk')
parser.add_argument('--arch_cover', type=str)
parser.add_argument('--obs_samples', type=str, help='File with samples names')
parser.add_argument('--decoding', type=str, help='Viterbi or aposteriory decoding')
parser.add_argument('--cut_off', type=float, help='Decoding cut off')

args = parser.parse_args()






dict_all = usfl.main_read2(args.prepared_file)
domain = usfl.read_bed(args.bed)
seq_start, seq_end = domain[0][0], domain[-1][1]
N = 2 # number of hidden states
GEN_time, MU, RR, L, Lambda_0 = usfl.read_par_HMM(args.HMM_par)
n_eu = len(dict_all[list(dict_all.keys())[0]]['Obs'])


#read the windows archaic covering 
with open(args.arch_cover,'r') as f:
    cover=f.readlines()
cover = [float(cover[i].replace('\n','')) for i in range(len(cover))]



gaps_numbers, seq_start_mas, seq_end_mas, len_mas=[], [], [], []
for i in range(len(domain)):
    if (domain[i][0] // L)*L + (domain[0][0]% L) >= domain[i][0]:
        seq_start_mas.append((domain[i][0] // L)*L + (domain[0][0]% L))
    else:
        seq_start_mas.append((domain[i][0] // L)*L + (domain[0][0]% L)+L)
    if (domain[i][1]//L)*L-1+(domain[0][0]% L) <= domain[i][1]:
        seq_end_mas.append((domain[i][1]//L)*L-1+(domain[0][0]% L))
    else:
        seq_end_mas.append((domain[i][1]//L)*L-1+(domain[0][0]% L)-L)        
        
    len_mas.append(int((domain[i][1]-domain[i][0])/1000))

    if i!=len(domain)-1:
        gaps_numbers.append([int((domain[i][1]-domain[0][0])/1000),int((domain[i+1][0]-domain[0][0])/1000)] )







#creating observations sequence (in gaps is zero)
SEQ=[]
N_ST=[]
for ind in range(n_eu):
    sq=np.vstack([usfl.make_obs_ref(dict_all, domain, ind, L,  'Outgroup'), usfl.make_obs_ref(dict_all, domain, ind, L,  'Archaic')])
    sq=sq.transpose()
    n_st = sq.max()+1
    SEQ.append(sq)
    N_ST.append(n_st)
SEQ=np.array(SEQ)
N_st=SEQ.max()+1


#split observations by windows, removing gaps
SEQ_mas=[]
arch_cover=[]
for i in range(len(len_mas)):
    p1=int((seq_start_mas[i]-seq_start_mas[0])/1000)
    p2=int((seq_end_mas[i]-seq_start_mas[0])/1000)
    SEQ_mas.append(SEQ[:,p1:(p2+1)])
    arch_cover.append(cover[p1:(p2+1)])



def run_daiseg(lmbd_opt,seq, n_st, idx, start, ar_cover):
    d = MU * L
    A = HMM.initA(L,RR, lmbd_opt[4]/d, lmbd_opt[3])
    

    
    B_our_mas = np.array([HMM.initB_arch_cover(MU,L, lmbd_opt, n_st, 0.1+i*0.1) for i in range(10)])
    B_Skov = HMM.initBwN(L, lmbd_opt[0:3], n_st)
    P=[0.97, 0.03]

    tracts_HMM =  HMM.get_HMM_tracts(HMM.viterbi_modified(seq [idx], P, A, B_our_mas, B_Skov, ar_cover))

    for k in range(N):
       for j in range(len(tracts_HMM[k])):
           tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]+start
           tracts_HMM[k][j][1]= L * (tracts_HMM[k][j][1]+1)+start-1

    return tracts_HMM

def run_daiseg_all(lmbd_0):
    tracts_HMM_mas=[]

    
    for idx in range(0, n_eu):    
        tracts_HMM=[[],[]]
        for i in range(len(SEQ_mas)):
            tr=run_daiseg(lmbd_0, SEQ_mas[i], N_st, idx, seq_start_mas[i], arch_cover[i])
            for j in range(N):   
               for k in tr[j]:             
                   tracts_HMM[j].append( k )
 

        tracts_HMM_mas.append([tracts_HMM[j] for j in range(N)])
    return tracts_HMM_mas



cut_off=args.cut_off
def run_daiseg_posterior(lmbd_opt,seq, n_st, idx, start, ar_cover, gaps_numbers, cut_off):
    d = MU * L
    A = HMM.initA(L,RR, lmbd_opt[4]/d, lmbd_opt[3])
    

    
    B_our_mas = np.array([HMM.initB_arch_cover(MU,L, lmbd_opt, n_st, 0.1+i*0.1) for i in range(10)])
    B_Skov = HMM.initBwN(L, lmbd_opt[0:3], n_st)
    P=[0.97, 0.03]

    tracts_HMM =  HMM.get_HMM_tracts(HMM.posterior(seq [idx], P, A, B_our_mas, B_Skov, ar_cover, gaps_numbers, cut_off))

    for k in range(N):
       for j in range(len(tracts_HMM[k])):
           tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]+start
           tracts_HMM[k][j][1]= L * (tracts_HMM[k][j][1]+1)+start-1

    return tracts_HMM


def run_daiseg_all_posterior(lmbd_0, cut_off, gaps_numbers):
    tracts_HMM_mas=[]

    
    for idx in range(0, n_eu):    
        tracts_HMM=[[],[]]
        for i in range(len(SEQ_mas)):
            tr=run_daiseg_posterior(lmbd_0, SEQ_mas[i], N_st, idx, seq_start_mas[i], arch_cover[i], gaps_numbers, cut_off)
            for j in range(N):   
               for k in tr[j]:             
                   tracts_HMM[j].append( k )
 

        tracts_HMM_mas.append([tracts_HMM[j] for j in range(N)])
    return tracts_HMM_mas

def EM_gaps(seq, lambda_0, n_st, cover):
    return EM.EM_algorithm_gaps(P, seq, n_st, MU, RR, lambda_0, epsilon, L, int(args.EM_steps), gaps_numbers, cover )

P=[0.95, 0.05]

epsilon = 1e-8


if args.EM=='no': 
    if args.decoding=='viterbi':
        Tracts_HMM_mas = run_daiseg_all(Lambda_0)
    else:
        Tracts_HMM_mas=run_daiseg_all_posterior(Lambda_0, cut_off, gaps_numbers)

if args.EM=='yes': 


    Lambda_opt = EM_gaps(np.array(SEQ[0 : (args.EM_samples+1)]), Lambda_0, N_st, cover)    
    if args.decoding=='viterbi':
        Tracts_HMM_mas = run_daiseg_all(Lambda_opt)
    else:
        Tracts_HMM_mas=run_daiseg_all_posterior(Lambda_opt, cut_off, gaps_numbers)



with open(args.obs_samples,'r') as f:
    names=f.readlines()
names=[str(names[i].replace('\n','')) for i in range(len(names))]

#write into file arg.o results #Sample #Haplotype_number #Archaic tracts
with open(args.o+'.archaic.txt', "w") as f:
   for i in range(len(Tracts_HMM_mas)):
       if i % 2 ==0:
           f.write(names[int(i // 2)]+'\t0\t'+str(Tracts_HMM_mas[i][1])+'\n')
       else:
           f.write(names[int(i // 2)]+'\t1\t'+str(Tracts_HMM_mas[i][1])+'\n')      


with open(args.o+'.modern.txt', "w") as f:
   for i in range(len(Tracts_HMM_mas)):
       if i % 2 ==0:
           f.write(names[int(i // 2)]+'\t0\t'+str(Tracts_HMM_mas[i][0])+'\n')
       else:
           f.write(names[int(i // 2)]+'\t1\t'+str(Tracts_HMM_mas[i][0])+'\n')  


