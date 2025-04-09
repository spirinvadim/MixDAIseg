import pysam
from collections import defaultdict 
import numpy as np


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


def neand_al_freq(fsamples, fvcf, ftracts, CHR):
    with open(fsamples,'r') as f:
        ids=f.readlines()
    ids=[ids[i].replace('\n','') for i in range(len(ids))]

    with open(ftracts,'r') as f:
        l=f.readlines()

    l=[l[i].split('\t')[2] for i in range(len(l))]
    l=[l[i].replace('\n','').replace('[','').replace(']','').split(',') for i in range(len(l))]
    l_mas=[]
    for j in l:
        if j!=[]:
        
            m=[[int(j[2*i]), int(j[2*i+1])] for i in range(int(len(j)/2))]       
        else:
            m=[[]]
        l_mas.append(m)
    
    dict_intervals={}
    k=0
    for j in range(len(ids)):
        dict_intervals[ids[j]]=[l_mas[k],l_mas[k+1]]
        
        k+=2

    panel = pysam.VariantFile(fvcf)
    
    samples = list(panel.header.samples)

    

    samples_dict = dict()
    for i in range(len(samples)):
        if samples[i] in ids:
            samples_dict[samples[i]] = i

#    dict_intervals2={}
#    for key, value in samples_dict.items():
#        dict_intervals2[value]=dict_intervals[key]


#    d_freq = defaultdict(lambda: 0) 
#    d_counts=defaultdict(lambda: 0)
    
#    freq_dict={}
#    for rec in panel.fetch(CHR):
        
#        gt=np.concatenate([[rec.samples.values()[idx]['GT'][0],rec.samples.values()[idx]['GT'][1]] for idx in samples_dict.values()])
#        in_id=np.concatenate([[point_in_set(rec.pos,dict_intervals2[idx][0]), point_in_set(rec.pos,dict_intervals2[idx][1])] 
#                              for idx in samples_dict.values()])
#        fr=0
#        n_fr=0
    
#        for i,j in zip(gt,in_id):
#            if j==True and i==1:
#                fr+=1
#                n_fr+=1
#            if j==True and i==0:
#                n_fr+=1
    
#        if n_fr!=0:
#            freq_dict[rec.pos]=[float(fr/n_fr), n_fr]

#    return freq_dict, samples_dict, dict_intervals
    return samples_dict, dict_intervals