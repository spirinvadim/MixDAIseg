import numpy as np
import math

N=5

cover_cut=0.25
step=0.05
#Ti: Introgression of Nd
#Tmex: Time in Mexico
def initA(Ti, Tmex,r,L,a,b,c1,c2) -> np.array:
    A = np.zeros((5,5))
    
   
    d=1-a-b
    A[0][1]=Ti*r*L*a*c1
    A[0][2]=Tmex*r*L*b*(1-c2)
    A[0][3]=Tmex*r*L*b*c2
    A[0][4]=Tmex*r*L*d
    A[0][0]=1-A[0][1]-A[0][2]-A[0][3]-A[0][4]
 
    A[1][0]=Ti*r*L*a*(1-c1)
    A[1][2]=Tmex*r*L*b*(1-c2)
    A[1][3]=Tmex*r*L*b*c2
    A[1][4]=Tmex*r*L*d
    A[1][1]=1-A[1][0]-A[1][2]-A[1][3]-A[1][4]
    
    A[2][0]=Tmex*r*L*a*(1-c1)
    A[2][1]=Tmex*r*L*a*c1
    A[2][3]=Ti*r*L*b*c2
    A[2][4]=Tmex*r*L*d
    A[2][2]=1-A[2][0]-A[2][1]-A[2][3]-A[2][4]
    
    A[3][0]=Tmex*r*L*a*(1-c1)
    A[3][1]=Tmex*r*L*a*c1
    A[3][2]=Ti*r*L*b*(1-c2)
    A[3][4]=Tmex*r*L*d
    A[3][3]=1-A[3][0]-A[3][1]-A[3][2]-A[3][4]

    A[4][0]=Tmex*r*L*a*(1-c1)
    A[4][1]=Tmex*r*L*a*c1
    A[4][2]=Tmex*r*L*b*(1-c2)
    A[4][3]=Tmex*r*L*b*c2
    A[4][4]=1-A[4][0]-A[4][1]-A[4][2]-A[4][3]
    
    return A

#Ti: Introgression of Nd
#Tea: Split between Asia and Europe
#Tmex: Time of migration to Mexixo
#Taf: Time out of Africa
#Tn: Time of Split between Nd and Sapiens

def initB( lmbd, n_st) -> np.array: 
    
    B = np.empty(shape=(5,n_st[1],n_st[2],n_st[3],n_st[4]))
    meani = lmbd[0]
    meann = lmbd[1]
    meanaf = lmbd[2]
    meanea =lmbd[3]
    meanmex =lmbd[4]



    
    Pi = np.empty(n_st[3])
    Pea = np.empty(max(n_st[0], n_st[1]))
    Pmex = np.empty(max([n_st[0], n_st[1], n_st[2]]))
    Paf=np.empty(max([n_st[0], n_st[1], n_st[2]]))
    Pn=np.empty(max([n_st[2], n_st[3]]))

    
    Pi[0]=np.exp(-meani)
    Pea[0]=np.exp(-meanea)
    Pmex[0]=np.exp(-meanmex)
    Paf[0]=np.exp(-meanaf)
    Pn[0]=np.exp(-meann)
    
    sumi=0
    sumea=0
    summex=0
    sumaf=0
    sumn=0

    for i in range(1, n_st[3]):
        Pi[i]=Pi[i-1]*meani/i
        sumi=sumi+Pi[i]

    for i in range(1, max([n_st[0], n_st[1], n_st[2]])):
        Pmex[i]=Pmex[i-1]*meanmex/i
        Paf[i]=Paf[i-1]*meanaf/i
        summex=summex+Pmex[i]
        sumaf=sumaf+Paf[i]

    for i in range(1, max(n_st[0], n_st[1])):
        Pea[i]=Pea[i-1]*meanea/i
        sumea=sumea+Pea[i]
        
    for i in range(1, max([n_st[2], n_st[3]])):
        Pn[i]=Pn[i-1]*meann/i
        sumn=sumn+Pn[i]

    


    Pea[0]=1-sumea
    Pi[0]=1-sumi
    Pmex[0]=1-summex
    Paf[0]=1-sumaf
    Pn[0]=1-sumn
    
    for i in range(n_st[0]): 
        for j in range(n_st[1]):
            for k in range(n_st[2]):
                for t in range(n_st[3]):
                    B[0][i][j][k][t]=Pmex[i]*Pea[j]*Paf[k]*Pn[t]
                    B[1][i][j][k][t]=Pmex[i]*Pea[j]*Pn[k]*Pi[t]
                    B[2][i][j][k][t]=Pea[i]*Pmex[j]*Paf[k]*Pn[t]
                    B[3][i][j][k][t]=Pea[i]*Pmex[j]*Pn[k]*Pi[t]
                    B[4][i][j][k][t]=Paf[i]*Paf[j]*Pmex[k]*Pn[t]    


    return B




def initB_arch_cover( lmbd, n_st, cover):
    
    B = np.empty(shape=(N,n_st[0],n_st[1], n_st[2], n_st[3]))
    
    
    meani, meann, meanaf, meanea, meanmex, meani2, meann2 = lmbd[0], lmbd[1], lmbd[2], lmbd[3] , lmbd[4], lmbd[0]*cover, lmbd[1]*cover
    
    Pi, Pi2 = np.empty(n_st[3]), np.empty(n_st[3])
    Pea = np.empty(max(n_st[0], n_st[1]))
    Pmex = np.empty(max([n_st[0], n_st[1], n_st[2]]))
    Paf=np.empty(max([n_st[0], n_st[1], n_st[2]]))
    Pn, Pn2=np.empty(max([n_st[2], n_st[3]])), np.empty(max([n_st[2], n_st[3]]))        
    
    

    Pi[0], Pn[0], Paf[0], Pea[0], Pmex[0], Pi2[0], Pn2[0] = np.exp(-meani), np.exp(-meann), np.exp(-meanaf), np.exp(-meanea), np.exp(-meanmex), np.exp(-meani2), np.exp(-meann2)
    sumi, sumaf, sumn, sumi2, sumn2, sumea, summex = [0 for i in range(7)]


    for i in range(1, n_st[3]):
        Pi[i]=Pi[i-1]*meani/i
        sumi=sumi+Pi[i]
        Pi2[i]=Pi2[i-1]*meani2/i
        sumi2=sumi2+Pi2[i]

    for i in range(1, max([n_st[0], n_st[1], n_st[2]])):
        Pmex[i]=Pmex[i-1]*meanmex/i
        Paf[i]=Paf[i-1]*meanaf/i
        summex=summex+Pmex[i]
        sumaf=sumaf+Paf[i]

    for i in range(1, max(n_st[0], n_st[1])):
        Pea[i]=Pea[i-1]*meanea/i
        sumea=sumea+Pea[i]
        
    for i in range(1, max([n_st[2], n_st[3]])):
        Pn[i]=Pn[i-1]*meann/i
        sumn=sumn+Pn[i]
        Pn2[i]=Pn2[i-1]*meann2/i
        sumn2=sumn2+Pn2[i]


        

    Pi[0], Pn[0], Paf[0], Pea[0], Pmex[0], Pi2[0], Pn2[0]  = 1-sumi, 1-sumn, 1-sumaf, 1-sumea, 1-summex,  1-sumi2, 1-sumn2
    
    for i in range(n_st[0]): 
        for j in range(n_st[1]):
            for k in range(n_st[2]):
                for t in range(n_st[3]):
                    B[0][i][j][k][t]=Pmex[i]*Pea[j]*Paf[k]*Pn2[t]
                    B[1][i][j][k][t]=Pmex[i]*Pea[j]*Pn[k]*Pi2[t]
                    B[2][i][j][k][t]=Pea[i]*Pmex[j]*Paf[k]*Pn2[t]
                    B[3][i][j][k][t]=Pea[i]*Pmex[j]*Pn[k]*Pi2[t]
                    B[4][i][j][k][t]=Paf[i]*Paf[j]*Pmex[k]*Pn2[t]    
                  
    return B



def initBwN( lmbd, n_st) -> np.array: 
    
    B = np.empty(shape=(N,n_st[0],n_st[1],n_st[2]))

    meann = lmbd[1]
    meanaf = lmbd[2]
    meanea =lmbd[3]
    meanmex =lmbd[4]
    
    Pea = np.empty(max(n_st[0], n_st[1]))
    Pmex = np.empty(max([n_st[0], n_st[1], n_st[2]]))
    Paf=np.empty(max([n_st[0], n_st[1], n_st[2]]))
    Pn=np.empty(max([n_st[2], n_st[3]]))



    

    Pea[0]=np.exp(-meanea)
    Pmex[0]=np.exp(-meanmex)
    Paf[0]=np.exp(-meanaf)
    Pn[0]=np.exp(-meann)
    
    sumi=0
    sumea=0
    summex=0
    sumaf=0
    sumn=0


    for i in range(1, max([n_st[0], n_st[1], n_st[2]])):
        Pmex[i]=Pmex[i-1]*meanmex/i
        Paf[i]=Paf[i-1]*meanaf/i
        summex=summex+Pmex[i]
        sumaf=sumaf+Paf[i]

    for i in range(1, max(n_st[0], n_st[1])):
        Pea[i]=Pea[i-1]*meanea/i
        sumea=sumea+Pea[i]
        
    for i in range(1, max([n_st[2], n_st[3]])):
        Pn[i]=Pn[i-1]*meann/i
        sumn=sumn+Pn[i]

    Pea[0]=1-sumea
    Pmex[0]=1-summex
    Paf[0]=1-sumaf
    Pn[0]=1-sumn
    
    for i in range(n_st[0]): 
        for j in range(n_st[1]):
            for k in range(n_st[2]):                
                B[0][i][j][k]=Pmex[i]*Pea[j]*Paf[k]
                B[1][i][j][k]=Pmex[i]*Pea[j]*Pn[k]
                B[2][i][j][k]=Pea[i]*Pmex[j]*Paf[k]
                B[3][i][j][k]=Pea[i]*Pmex[j]*Pn[k]
                B[4][i][j][k]=Paf[i]*Paf[j]*Pmex[k]

    return B


def viterbi(V, initial_distribution, a, b):
    
    T = len(V)
    M = a.shape[0]
 
    omega = np.zeros((T, M))
    omega[0, :] = np.log(initial_distribution * b[:, V[0][0],V[0][1],V[0][2],V[0][3]])
 
    prev = np.zeros((T - 1, M))
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(a[:, j]) + np.log(b[j, V[t][0], V[t][1], V[t][2], V[t][3]])
 
            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)
 
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
 
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
 
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
 
    # Convert numeric values to actual hidden states
 
    result = []
    for s in S:
        if s == 0:
            result.append(0)
        elif s == 1:
            result.append(1)
        elif s == 2:
            result.append(2)
        elif s == 3:
            result.append(3)
        elif s == 4:
            result.append(4)

    return result


def get_HMM_tracts(seq):
    migrating_tracts = []
    for i in range(N):
        migrating_tracts.append([])
    start=0
    for i in range(1,len(seq)):
        if seq[i]!=seq[i-1]:
            migrating_tracts[seq[i-1]].append([start,i-1])
            start=i
    migrating_tracts[seq[len(seq)-1]].append([start,len(seq)-1])
    return migrating_tracts






# define transition matrix


def initA2(Ti, Tmex, r, L, a, b, c1, c2):
    tr_mat = np.zeros((5, 5))
    
    
    f=1-a-b
    
    et_a = math.exp(-r*L*Ti)
    et_m = math.exp(-r*L*Tmex)
    
    tr_mat[0, 0] = et_a * et_m + ((1 - et_a) * et_m+(1 - et_m) * a) * (1-c1)
    tr_mat[0, 1] = ((1 - et_a) * et_m + (1 - et_m) * a) * c1
    tr_mat[0, 2] = (1 - et_m) * b * (1-c2)
    tr_mat[0, 3] = (1 - et_m) * b * c2
    tr_mat[0, 4] = (1 - et_m) * f
    
    tr_mat[1, 0] = ((1 - et_a) * et_m + (1 - et_m) * a) * (1-c1)
    tr_mat[1, 1] = et_a * et_m + ((1 - et_a) * et_m + (1 - et_m) * a) * c1
    tr_mat[1, 3] = (1 - et_m) * b * c2
    tr_mat[1, 2] = (1 - et_m) * b * c2
    tr_mat[1, 4] = (1 - et_m) * f 
    
    
    tr_mat[2, 0] = (1 - et_m) * a * (1-c1)
    tr_mat[2, 1] = (1 - et_m) * a * c1
    tr_mat[2, 2] = et_a * et_m + ((1-et_a) * et_m +(1 - et_m) * b) * (1-c2)
    tr_mat[2, 3] = ((1 - et_a) * et_m + (1 - et_m) * b) *c2
    tr_mat[2, 4] = (1 - et_m) * f
    

    tr_mat[3, 0] = (1 - et_m) * a * (1-c1)
    tr_mat[3, 1] = (1 - et_m) * a * c1
    tr_mat[3, 3] = et_a * et_m + ((1-et_a) * et_m+(1 - et_m) * b) * c2
    tr_mat[3, 2] = ((1 - et_m) * et_m + (1 - et_m) * b) * (1-c2)
    tr_mat[3, 4] = (1 - et_m) * f
    
    tr_mat[4, 4] = et_m + (1-et_m) * f
    tr_mat[4, 0] = (1 - et_m) * a * (1-c1)
    tr_mat[4, 1] = (1 - et_m) * a * c2
    tr_mat[4, 3] = (1 - et_m) * b * c2
    tr_mat[4, 2] = (1 - et_m) * b * (1-c2)
    
    return tr_mat

def viterbi_modified(V, initial_distribution, a, b_our_mas, b_Skov, archaic_cover):
    
    T = len(V)
    M = a.shape[0]
    
 
    omega = np.zeros((T, M))
   

    if archaic_cover[0]>=cover_cut:
        omega[0, :] = np.log(initial_distribution * b_our_mas[int((archaic_cover[0]-cover_cut)/step),:, V[0][0],V[0][1], V[0][2], V[0][3]])

    else:
        
        omega[0, :] = np.log(initial_distribution * b_Skov[:, V[0][0],V[0][1], V[0][2]])        
 
    prev = np.zeros((T - 1, M))
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            if archaic_cover[t]>=cover_cut:
                probability = omega[t - 1] + np.log(a[:, j]) + np.log(b_our_mas[int((archaic_cover[t]-cover_cut)/step),j, V[t][0], V[t][1], V[t][2], V[t][3]])
                
            else:
                probability = omega[t - 1] + np.log(a[:, j]) + np.log(b_Skov[j, V[t][0], V[t][1], V[t][2]])
                
               
 
            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)
 
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
 
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
 
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
 
    # Convert numeric values to actual hidden states
 
    result = []
    for s in S:
        if s == 0:
            result.append(0)
        elif s == 1:
            result.append(1)
        elif s == 2:
            result.append(2)
        elif s == 3:
            result.append(3)
        elif s == 4:
            result.append(4)             

    return result
    



