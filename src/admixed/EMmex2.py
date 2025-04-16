import HMMmex2 as HMM
import numpy as np
from numpy import linalg as LNG 
import math
import math
from numpy import linalg as LNG 

import useful2 as usfl


N=5
# each window is covered by cover[t] parameter which in [-0.001...0.999]. Cover_cut remove windows with low covering and transforms method to Skov. 
cover_cut = 0.25 
step = 0.05

#transitions and emissions in case of gaps in observable chromosome (centromeres/telomeres), not archaic gaps. 
a_gaps=np.identity(N)
b_gaps=np.ones((N,1))    



def alpha_scaled_opt_gaps(a,b, b_mas, o, p, gaps, cover):   
    c = np.zeros(len(o)) #scaling factors, которые как раз позволяют не обнулиться
    
    alpha = np.zeros((N, len(o)))

    if cover[0]>=cover_cut:
        alpha[:,0] = b_mas[int((cover[0]-cover_cut)/step),: , o[0][0],o[0][1], o[0][2], o[0][3]]*p
    else:
        alpha[:, 0] = b[:, o[0][0],o[0][1],o[0][2]] * p
    
    
    c[0] = 1 / sum(alpha[:, 0])
    alpha[:, 0] = alpha[:, 0] / sum(alpha[:, 0])
    for t in range(1, len(o)):   
        if usfl.point_in_set(t, gaps)==True:

            for i in range(0, N):
                
                alpha[i, t] = np.dot(alpha[:, t-1],a_gaps[:,i]) * b_gaps[i][0]     

        else:
            if cover[t] >= cover_cut:
                for i in range(0, N):  
                    alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b_mas[int((cover[t]-cover_cut)/step), i, o[t][0],o[t][1], o[t][2], o[t][3]]   

    
            else:
        
                for i in range(0, N):            
                    alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b[i, o[t][0],o[t][1],o[t][2]] 
            
        c[t] = 1 / sum(alpha[:,t]) #сохраняем множители      

            
        alpha[:, t] = alpha[:, t] / sum(alpha[:,t])     
        


        
    return alpha, c

# Backward procedure. Scaled case.
def beta_scaled_opt_gaps(a,b, b_mas, o, scaling_factors, gaps, cover):
    
    beta = np.zeros((N, len(o)))
    
    length = len(o)
    beta[:, length - 1] = np.ones(N)*scaling_factors[length-1] 
    

    for t in range(len(o)-2,-1,-1):        
        if usfl.point_in_set(t, gaps)==True:
            for i in range(0, N):             
                for l in range(0, N):
                    beta[i, t] += a_gaps[i, l] * b_gaps[l][0] * beta[l, t+1]       
        else:

            if cover[t+1] >= cover_cut:
                for i in range(0, N):             
                    for l in range(0, N):
                        beta[i, t] += a[i, l] * b_mas[int((cover[t+1]-cover_cut)/step),l, o[t+1][0],o[t+1][1], o[t+1][2],o[t+1][3]] * beta[l, t+1]
                
            
            else:    
                for i in range(0, N):             
                    for l in range(0, N):
                        beta[i, t] += a[i, l] * b[l, o[t+1][0],o[t+1][1],o[t+1][2]] * beta[l, t+1]
                
        beta[:, t] = beta[:, t] * scaling_factors[t]

    return beta 

# gamma matrix
def def_gamma_gaps(alpha, beta, gaps):
        
    gamma = np.zeros((N,len(alpha[0])))
    for m in range(0, len(alpha[0])):
        denom = sum(alpha[:, m]*beta[:,m])
        
        for i in range(0,N):
            gamma[i, m] = (alpha[i, m] * beta[i, m]) / denom
#            if gamma[i,m]==0:
#                print('oh mo', alpha[i, m], beta[i, m], i, m)
#                break
    return gamma


# ksi[i, j, t]
def def_ksi_gaps( a, b, o, alpha, beta, gaps):
    
    M = len(o)
    ksi = np.zeros((N, N, M-1))
    
    for t in range(0, M-1):
        if usfl.point_in_set(t, gaps)==True:
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a_gaps[i, j] * b_gaps[j][0] * beta[j, t+1]        
            for i in range(0, N):
                for j in range(0, N):       
                    ksi[i, j, t] = (alpha[i, t]*a_gaps[i, j]*b_gaps[j,0]* beta[j, t+1]) / denom
        else:
        
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a[i, j] * b[j, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[j, t+1]
                
        
            for i in range(0, N):
                for j in range(0, N):
                    ksi[i, j, t] = (alpha[i, t]*a[i, j]*b[j, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[j, t+1]) / denom
    
    return ksi









def new_lambda_mex_common_gaps(o_mas, gamma_mas, gaps):

    s=0    
    for j in gaps:
        s+=j[1]-j[0]+1
        
    lmbd=0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i]
        for t in range(1, len(o), 1):
            if usfl.point_in_set(t, gaps)==False:
                lmbd += o[t, 0] * ( gamma[0, t] + gamma[1, t]) + o[t, 1] * (gamma[2, t] + gamma[3, t]) + o[t, 2]  * gamma[4, t]     
        
    return lmbd/((len(o)-s-1)*len(o_mas))




def new_lambda_af_common_gaps(o_mas, gamma_mas, gaps):
    nom, denom = 0, 0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i] 
        for t in range(1, len(o), 1): 
            if usfl.point_in_set(t, gaps)==False:  
                nom += o[t, 2] * ( gamma[0, t] + gamma[2, t]) + (o[t, 0] + o[t, 1]) * gamma[4, t]
                denom += gamma[0, t] + gamma[2, t]+ 2 * gamma[4, t]    

    
    return nom/ denom

def new_lambda_ea_common_gaps(o_mas, gamma_mas, gaps):
    nom, denom = 0, 0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i]
        for t in range(1, len(o), 1):  
            if usfl.point_in_set(t, gaps)==False:
                nom += o[t, 1] * ( gamma[0, t] + gamma[1, t]) + o[t, 0]  * (gamma[2, t]+gamma[3,t])
                denom += gamma[0, t] + gamma[1, t] +  gamma[2, t] + gamma[3, t]

    return nom/ denom
 
 
 
 
 

def new_lambda_n_common_gaps(o_mas, gamma_mas, gaps, cover):

    s=0    
    for j in gaps:
    
        s+=j[1]-j[0]+1

    lmbd = 0
    denom=0
 
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i]

        for t in range(1, len(o), 1):
            if usfl.point_in_set(t, gaps)==False:
                if cover[t] >= cover_cut:
                    
                    lmbd +=  o[t, 3] * ( gamma[0, t] + gamma[2, t]+gamma[4, t]) + o[t, 2] * (gamma[1, t] + gamma[3, t]) 
                    denom += (int(cover[t]/step)*step) * ( gamma[0, t] + gamma[2, t]+gamma[4, t]) + (gamma[1, t] + gamma[3, t]) 
                else:
                    lmbd += o[t, 2] * (gamma[1, t] + gamma[3, t]) 
                    denom += (gamma[1, t] + gamma[3, t])
  
   
    return lmbd/ denom

def new_lambda_i_common_gaps(o_mas, gamma_mas, gaps, cover):
    nom=0
    denom=0
    
    for i in range(len(o_mas)):
        o=o_mas[i]
        g=gamma_mas[i]   
        
        
        
       
        for t in range(1, len(o), 1):   
                    

                
            if usfl.point_in_set(t, gaps)==False:                
                if cover[t] >= cover_cut:
                

                    nom += o[t, 3] * (g[1, t] + g[3, t])
                    
                    
                    denom += (int(cover[t]/step)*step) * (g[1, t] + g[3, t])
                  
    if denom==0: # only in case of no archaic covering
        return 0
                
    return nom / denom
    
    
    
    
    
    
def E_step_gaps(cut,  p, O, n_st, mu,rr, lambd_old, gaps, cover, matrix_type):

    d=mu*cut

    if matrix_type=='simple':

        a = HMM.initA(lambd_old[5]/d, lambd_old[6]/d, rr, cut, lambd_old[7],  lambd_old[8],  lambd_old[9],  lambd_old[10])    
    else:
        a = HMM.initA2(lambd_old[5]/d, lambd_old[6]/d, rr, cut, lambd_old[7],  lambd_old[8],  lambd_old[9],  lambd_old[10])  
    
    b_our_mas = np.array([HMM.initB_arch_cover( lambd_old, n_st, cover_cut+i*step) for i in range(int((1-cover_cut)/step)+1)])

    
    b_Skov = HMM.initBwN(lambd_old[0:5], n_st)

    GAMMA=[]
    for o in O:
        alpha, sc_factors = alpha_scaled_opt_gaps(a,b_Skov, b_our_mas, o, p, gaps, cover)
        beta = beta_scaled_opt_gaps(a, b_Skov, b_our_mas, o, sc_factors, gaps, cover)    
        gamma = def_gamma_gaps(alpha, beta, gaps)    


        GAMMA.append(gamma)       

        
        
        
    

    return  [new_lambda_i_common_gaps(O, GAMMA, gaps, cover),new_lambda_n_common_gaps(O, GAMMA, gaps, cover), 
                new_lambda_af_common_gaps(O, GAMMA, gaps), 
                new_lambda_ea_common_gaps(O, GAMMA,gaps),
                new_lambda_mex_common_gaps(O, GAMMA, gaps),
                lambd_old[5], lambd_old[6],
                lambd_old[7], lambd_old[8], lambd_old[9], lambd_old[10]] 
                
                
                
def EM_common_gaps(p, o_mas, n_states, mut_rate, rr, lambda_0, epsilon, cut,  em_steps, gaps, cover, matrix_type):
    d=mut_rate * cut
    lmbd = np.array(lambda_0)


    for i in range(em_steps):
        


        lmbd_new = np.array(E_step_gaps(cut, p, o_mas, n_states, mut_rate, rr,lmbd,  gaps, cover, matrix_type))
        print(LNG.norm(lmbd_new-lmbd))
        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new

    print('Number EM steps',i)
    return lmbd
