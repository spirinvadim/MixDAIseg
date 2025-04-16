import numpy as np

from functools import reduce
from operator import xor
from itertools import chain


def point_in_set(p, m_a):
    if m_a==[[]]:
        return False
    def point_in(p, a):
        if p>=a[0] and p<=a[1]:
            return True
        else:
            return False
    f=False
    for j in m_a:
        f=point_in(p,j)
        if f==True:
            return f
    return f

# two set of intervals intersection 
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

def read_bed(f):


    with open(f,'r') as f:
        domain=f.readlines()
    domain=[domain[i].replace('\n','').split('\t')[1:3] for i in range(len(domain))]    

    for i in  range(len(domain)):
        domain[i][0]=int(domain[i][0])
        domain[i][1] = int(domain[i][1])
    return domain
    
    
def read_aa_file(f_aa): #read file with ancestral alleles
    with open(f_aa,'r') as f:
        l=f.readlines()
    dct={}
    for  i in l:    
        m=i.replace('\n','').split('\t')
        dct[int(m[1])]=m[2]
    return dct




def main_read(file, f_aa): # return the dictionary with REF,ALT, Outgroup and Archaic allels. And Ancestral allles (if possible)
    dct_AA=read_aa_file(f_aa)
    with open(file,'r') as f :
        lns=f.readlines()
        
    nested_dict = {}
        
    for i in range(1,len(lns)):
        m=lns[i].replace(' \n','').split('\t') 
        
        eu = [int(m[3][k])  for k in range(len(m[3])) if m[3][k]!='.']      
        na = [int(m[4][k])  for k in range(len(m[4])) if m[4][k] != '.']    
        af=[int(m[5][k])  for k in range(len(m[5])) if m[5][k]!='.']  
        arch=[int(m[6][k])  for k in range(len(m[6])) if m[6][k]!='.']  
        o=m[7].split(' ')
        o = [int(o[j]) for j in range(len(o))]
   
        nested_dict[int(m[0])]={'REF': m[1], 'ALT':m[2], 'EU': eu, 'NA':na, 'AF': af, 'Archaic': arch, 'Obs': o  }

    for j in dct_AA.keys():
        AA=dct_AA[j]
        if nested_dict.get(j) is not None:
            st=nested_dict[j]['REF']+nested_dict[j]['ALT'].replace(',','')
        
        
            for i in range(len(st)):
                if st[i]==AA:
                    k=i        
        
            nested_dict[j]['AA']=k
    for i, j in nested_dict.items():
        if j.get('AA') is None:
            nested_dict[i]['AA']='.'

    return nested_dict

def main_read2(file2): # return the dictionary with REF,ALT, Outgroup and Archaic allels. And Ancestral allles (if possible)
    
    with open(file2,'r') as f :
        lns=f.readlines()
    
    nested_dict = {}
        
    for i in range(1,len(lns)):
        
        m=lns[i].replace(' \n','').split('\t') 

        eu=m[4].replace(',','')               
        if eu!='':
            eu = [int(eu[k])  for k in range(len(eu))]  
        else:
            eu = []
            
            
        na=m[5].replace(',','')             
                  
        if na!='':
            na = [int(na[k])  for k in range(len(na))]  
        else:
            na = []
            
        af=m[6].replace(',','')             
                  
        if af!='':
            af = [int(af[k])  for k in range(len(af))]  
        else:
            af = []
        
        arch= m[7].replace('.','')
        if len(arch)!=0:
            arch=[int(k) for k in arch.split(',')]
        else:
            arch=[]

        o=m[8].replace('\n','').split(' ')
        o=[int(k) for k in o ]
        


        if m[3]=='.': 
            s=''
            nested_dict[int(m[0])]={'REF': m[1], 'ALT':m[2], 'EU': eu, 'NA':na, 'AF':af, 'Archaic': arch, 'Obs': o  }
            
        else:
            s=m[3]   
            nested_dict[int(m[0])]={'REF': m[1], 'ALT':m[2], 'EU': eu, 'NA':na, 'AF':af, 'Archaic': arch, 'Obs': o, 'AA': s  }


    return nested_dict
    
    
def read_par_HMM(file):
    with open(file,'r') as f:

        GEN_time = float(f.readline())
        MU = float(f.readline())
        RR = float(f.readline())
        L = int(f.readline())

        Lambda_0=np.zeros(11)
        Lambda_0[1] = float(f.readline())/GEN_time*MU*L
        Lambda_0[2] = float(f.readline())/GEN_time*MU*L
        Lambda_0[0] = float(f.readline())/GEN_time*MU*L
        Lambda_0[3] = float(f.readline())/GEN_time*MU* L
        Lambda_0[4] = float(f.readline())/GEN_time*MU* L
        Lambda_0[5] = float(f.readline())/GEN_time*MU* L
        Lambda_0[6] = float(f.readline())/GEN_time*MU* L
        Lambda_0[9] = float(f.readline())
        Lambda_0[10] = Lambda_0[9]
        Lambda_0[7] = float(f.readline())
        Lambda_0[8] = float(f.readline())

    return GEN_time, MU, RR, L, Lambda_0

def make_obs_ref(nested_dict, domain, ind, L,  ref):


        
    start, end = domain[0][0], domain[-1][1]
    if float(end-start) % L != 0:
        T = int((end-start)/ L) + 1
    obs_ref = np.zeros(T, dtype=int)

    
    for k in nested_dict.keys():
        j=int((k-start)/L) #window-positions
        
        if nested_dict[k].get('AA') is None:
        
            if nested_dict[k]['Obs'][ind] not in nested_dict[k][ref] and len(nested_dict[k][ref])!=0:
                obs_ref[j]+=1

                    
        else:
            if nested_dict[k]['Obs'][ind] not in nested_dict[k][ref] and len(nested_dict[k][ref])!=0 and nested_dict[k]['Obs'][ind]!=nested_dict[k]['AA']:
                
                obs_ref[j]+=1
                


    return np.array(obs_ref)



def make_obs_indep(nested_dict, domain, ind, L):
    start, end = domain[0][0], domain[-1][1]
    if float(end-start) % L != 0:
        T = int((end-start)/ L) + 1
    


    
    for k in nested_dict.keys():
        j=int((k-start)/L) #window-positions
        
        ref_n=0
        s=[0 for i in range(4)]
        
        for ref in ['EU', 'NA', 'AF', 'Archaic']:
        
            if nested_dict[k].get('AA') is None:        
                if nested_dict[k]['Obs'][ind] not in nested_dict[k][ref] and len(nested_dict[k][ref])!=0:
                    s[ref_n] += 1
            else:
                if nested_dict[k]['Obs'][ind] not in nested_dict[k][ref] and len(nested_dict[k][ref])!=0 and nested_dict[k]['Obs'][ind]!=nested_dict[k]['AA']:                
                    s[ref_n] += 1
            ref_n+=1
            
        o_mas=[ [np.zeros(T, dtype=int) for i in range(4)] for i in range(5)]
        
        #Modern EU
        
        if [s[0], s[1], s[2], s[3]] == [1,1,1,1]:     
            o_mas[0][0][j] += 1
        if [s[0], s[1], s[2], s[3]] == [0,1,1,1]:     
            o_mas[0][1][j] += 1           
        if [s[0], s[1], s[2], s[3]] == [0,0,1,1]:     
            o_mas[0][2][j] += 1     
        if [s[0], s[1], s[2], s[3]] == [0,0,0,1]:     
            o_mas[0][3][j] += 1  

        #Archaic EU
        if [s[0], s[1], s[2], s[3]] == [1,1,1,1]:     
            o_mas[1][0][j] += 1
        if [s[0], s[1], s[2], s[3]] == [0,1,1,1]:     
            o_mas[1][1][j] += 1           
        if [s[0], s[1], s[2], s[3]] == [0,0,1,1]:     
            o_mas[1][2][j] += 1     
        if [s[0], s[1], s[2], s[3]] == [0,0,1,0]:     
            o_mas[1][3][j] += 1    


             
  
        if [s[0], s[1], s[2], s[3]] == [1,0,1,1]:     
            o_mas[2][0][j] += 1
        if [s[0], s[1], s[2], s[3]] == [1,1,1,1]:     
            o_mas[2][1][j] += 1           
        if [s[0], s[1], s[2], s[3]] == [0,0,1,1]:     
            o_mas[2][2][j] += 1     
        if [s[0], s[1], s[2], s[3]] == [0,0,0,1]:     
            o_mas[2][3][j] += 1   
            
            
        if [s[0], s[1], s[2], s[3]] == [1,0,1,1]:     
            o_mas[3][0][j] += 1
        if [s[0], s[1], s[2], s[3]] == [1,1,1,1]:     
            o_mas[3][1][j] += 1           
        if [s[0], s[1], s[2], s[3]] == [0,0,1,1]:     
            o_mas[3][2][j] += 1     
        if [s[0], s[1], s[2], s[3]] == [0,0,1,0]:     
            o_mas[3][3][j] += 1           
        

        if [s[0], s[1], s[2], s[3]] == [1,1,0,1]:     
            o_mas[4][0][j] += 1     
        if [s[0], s[1], s[2], s[3]] == [1,1,1,1]:     
            o_mas[4][2][j] += 1     
        if [s[0], s[1], s[2], s[3]] == [0,0,0,1]:     
            o_mas[4][3][j] += 1   
            
    return o_mas

def read_out(file):
    with open(file, 'r') as f:
        l=f.readlines()

    l=[l[i].split('\t') for i in range(len(l))]
    l=[l[i][2] for i in range(len(l))]

    nd_HMM_tracts=[]
    for j in range(len(l)):
        m1=[]
        l2=l[j].replace('[[','').replace(']]','').replace('\n','').split('], [')
        for j1 in range(len(l2)):
            if len(l2[j1])>2:
                m2=l2[j1].split(', ')
                m2=[int(m2[j2]) for j2 in range(len(m2))]
            else:
                m2=[]
            m1.append(m2)
        nd_HMM_tracts.append(m1)
    return nd_HMM_tracts


def daiseg_mex_tracts(prefix, n_mexicans):
    tr_hmm_eu, tr_hmm_nd_eu, tr_hmm_na, tr_hmm_nd_na, tr_hmm_af= [],[], [], [], []
    for idx in range(n_mexicans):
        tr_hmm_eu.append(read_out(str(prefix)+'.modern.eu.txt')[idx])
        tr_hmm_nd_eu.append(read_out(str(prefix)+'.neand.eu.txt')[idx])
        tr_hmm_na.append(read_out(str(prefix)+'.modern.na.txt')[idx])
        tr_hmm_nd_na.append(read_out(str(prefix)+'.neand.na.txt')[idx])
        tr_hmm_af.append(read_out(str(prefix)+'.af.txt')[idx])
    

    return tr_hmm_eu, tr_hmm_nd_eu, tr_hmm_na, tr_hmm_nd_na, tr_hmm_af
