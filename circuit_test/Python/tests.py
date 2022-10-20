'''Just an area to test out different functions if it works or not'''
'''import numpy as np

def mat_ext(M):
    a = np.shape(M)
    M1=np.zeros((a[0]+1,a[1]+1),M.dtype)
    M1[:-1,:-1]=M
    return M1

def Vs_assigner(M,x):
    Value = np.array([[x]])
    a = np.concatenate((M, Value), axis=0)
    return a

def add_values(M,node1,node2):
    # Create the column for voltage stamp
    size_a = np.shape(M)
    va= np.zeros((size_a[0],1))
    va[node1-1,0] = 1
    va[node2-1,0] = -1
    print(va)
    
    # Create the row for voltage stamp
    zero_ext = np.zeros((1,1))
    ha = va[..., None]
    haz = np.concatenate((ha,zero_ext),axis=None)
    
    #Contatenate both for the new matrix
    M1 = np.hstack((M,va))
    M2 = np.vstack((M1,haz))
    
    return M2

M=np.arange(16).reshape((4,4))    

M1=np.arange(4).reshape((4,1))    

print(Vs_assigner(M1,20))

print(M)
#a = mat_ext(M)
a = add_values(M,2,3)
b = add_values(a,1,2)

print(a)
print(b)'''

import numpy as np
def F(i):
    x = np.zeros((len(i),1))
    return x
x = [5]
lol = [
    1,
    1,
    1,
    1,
    1
]

print(x)