import numpy as np
from sympy import * 
import seaborn as sns
import pandas as pd

def r_beta(length, beta, element_number):
    k = length * beta
    r = np.zeros((4,1))
    
    if element_number == 1: #bound 1-2
        r[0] = k
        r[1] = k
        
    elif element_number == 2: #bound 2-3
        r[1] = k
        r[2] = k
        
    elif element_number == 3: #bound 3-4
        r[2] = k
        r[3] = k
        
    elif element_number == 4: #bound 4-1
        r[3] = k
        r[0] = k
    
    else:
        print("R_beta Error")
    
    return r
        
        
def k_alfa(length, constant, element_number):
    a1 = (-2)*length*constant/3
    a2 = (-length)*constant/3
    k_alfa = np.zeros((4,4))

    if element_number == 1: #bound 1-2
        k_alfa[0,0] = a1
        k_alfa[0,1] = a2
        k_alfa[1,0] = a2
        k_alfa[1,1] = a1
        
    elif element_number == 2: #bound 2-3
        k_alfa[1,1] = a1
        k_alfa[1,2] = a2
        k_alfa[2,1] = a2
        k_alfa[2,2] = a1
        
    elif element_number == 3: #bound 3-4
        k_alfa[2,2] = a1
        k_alfa[2,3] = a2
        k_alfa[3,2] = a2
        k_alfa[3,3] = a1
        
    elif element_number == 4: #bound 4-1
        k_alfa[0,0] = a1
        k_alfa[0,3] = a2
        k_alfa[3,0] = a2
        k_alfa[3,3] = a1
        
    else:
        print(element_number)
        print("k_alfa error")
        
        
    return k_alfa

def r_q(a,b,q):
    r = np.ones((4,1))
    r_q = a*b*q*r
    
    return r_q


def k_k(a,b,k_x,k_y):
    k_k = np.array([ 
        (k_y*a/(3*b)+k_x*b/(3*a),  k_y*a/(6*b)-k_x*b/(3*a), -k_y*a/(6*b)-k_x*b/(6*a),  k_x*b/(6*a)-k_y*a/(3*b)),
        (k_y*a/(6*b)-k_x*b/(3*a),  k_y*a/(3*b)+k_x*b/(3*a),   k_x*b/(6*a)-k_y*a/(3*b), -k_x*b/(6*a)-k_y*a/(6*b)),
       (-k_y*a/(6*b)-k_x*b/(6*a), -k_y*a/(3*b)+k_x*b/(6*a),   k_y*a/(3*b)+k_x*b/(3*a), -k_x*b/(3*a)+k_y*a/(6*b)),
       (-k_y*a/(3*b)+k_x*b/(6*a), -k_y*a/(6*b)-k_x*b/(6*a),   k_y*a/(6*b)-k_x*b/(3*a),  k_x*b/(3*a)+k_y*a/(3*b))
      ])
    return k_k

def nodes_temperature(nodes_q, elements_q, A, r_beta_sum, k_k_sum):
    T = np.zeros((nodes_q,1))
    matrix_global = np.zeros((nodes_q,nodes_q))    
    A = A - 1 
    i = 0
    j = 0
    for i in range(elements_q):
        for j in range(4):
            T[int(A[i,j])] = T[int(A[i,j])] + r_beta_sum[i,j]
               
    n = 0
    i = 0
    j = 0
    for n in range(elements_q):
        for i in range(4):
            for j in range(4):
                matrix_global[int(A[n,i]),int(A[n,j])] = matrix_global[int(A[n,i]),int(A[n,j])] + k_k_sum[i,j,n]
         
    return np.dot(np.linalg.inv(matrix_global),T)

def elements_temperature(a, b,A, nodes_temperature):       
    A = A - 1
    y = symbols('y')
    x = symbols('x')    
    s = x - a
    t = y - b
    
    N1 = ((a-s)*(b-t))/(4*a*b)
    N2 = ((a+s)*(b-t))/(4*a*b)
    N3 = ((a+s)*(b+t))/(4*a*b)
    N4 = ((a-s)*(b+t))/(4*a*b)
    
    N_t = np.array([N1, N2, N3, N4])        
    elements_temperature = zeros(elements_q,1)
    i = 0

    for i in range(elements_q):
        Temperature = [nodes_temperature[int(A[i,0])], nodes_temperature[int(A[i,1])], nodes_temperature[int(A[i,2])], nodes_temperature[int(A[i,3])]]
        q = np.transpose(N_t)*np.transpose(Temperature)
        elements_temperature[i] = sum(sum(q))
        
    return elements_temperature

    

if __name__ == '__main__':
    #Data set
    h = 55
    T_inf = 20
    k_x = 45
    k_y = 45
    n_2 = 0
    q = 0
    beta_1 = 5000
    beta_2 = h * T_inf
    alfa = -h
    a = 0.01/2
    b = a
    
    #Number of nodes and elements
    elements_q = 10
    nodes_q = 22
    
    #Target matrixes 
    A = np.zeros((elements_q,4))
    r_beta_sum = np.zeros((elements_q,4))
    k_k_sum = np.zeros((4,4,elements_q))
   
    #Elements
    
    elements_data = np.array([[0,0,alfa,alfa,alfa,0,beta_2,beta_2,beta_2,[5,6,2,1]],
                                [1,0,alfa,alfa,alfa,0,beta_2,beta_2,beta_2,[7,8,4,3]],
                                [2,0,alfa,alfa,0,0,beta_2,0,beta_2,[9, 10, 6, 5]],
                                [3,0,alfa,0,alfa,0,beta_2,0,beta_2,[14, 15, 8, 7]],
                                [4,0,0,0,alfa,beta_1,0,0,beta_2,[16, 17, 10, 9]],
                                [5,0,0,alfa,0,beta_1,0,beta_2,0,[17, 18, 11, 10]],
                                [6,0,0,alfa,0,beta_1,0,beta_2,0,[18, 19, 12, 11]],
                                [7,0,0,alfa,0,beta_1,0,beta_2,0,[19, 20, 13, 12]],
                                [8,0,0,alfa,0,beta_1,0,beta_2,0,[20, 21, 14, 13]],
                                [9,0,alfa,0,0,beta_1,beta_2,0,0,[21, 22, 15, 14]]])
    
    for row in range(len(elements_data)):
        k_k_11 = k_k(a,b,k_x,k_y)
        
        if elements_data[row,1] == 0:
            k_alfa_112 = 0
        else:
            k_alfa_112 = k_alfa(2*b,elements_data[row,1],1) #Brzeg 1-2 wezly 5 6
        if elements_data[row,2] == 0:
            k_alfa_123 = 0
        else:
            k_alfa_123 = k_alfa(2*a,elements_data[row,2],2)   #Brzeg 2-3 wezly 6 2
        
        if elements_data[row,3] == 0:
            k_alfa_134 = 0
        else:
            k_alfa_134 = k_alfa(2*b,elements_data[row,3],3) #Brzeg 3-4 wezly 2 1
        
        if elements_data[row,4] == 0:
            k_alfa_141 = 0
        else:
            k_alfa_141 = k_alfa(2*a,elements_data[row,4],4)   #Brzeg 4-1 wezly 1 5    
               
        k_k_1 = k_k_11 + k_alfa_112 + k_alfa_123 + k_alfa_134 + k_alfa_141
        
        if elements_data[row,5] == 0:
            r_beta_112 = 0
        else:
            r_beta_112 = r_beta(2*b,elements_data[row,5],1) # r_q(2*a,2*b,q) #Brzeg 1-2 wezly 5 6
        
        if elements_data[row,6] == 0:
            r_beta_123 = 0
        else:
            r_beta_123 = r_beta(2*a,elements_data[row,6],2) #Brzeg 2-3 wezly 6 2
        
        if elements_data[row,7] == 0:
            r_beta_134 = 0
        else:
            r_beta_134 = r_beta(2*b,elements_data[row,7],3) #Brzeg 3-4 wezly 2 1
        
        if elements_data[row,8] == 0:
            r_beta_141 = 0
        else:
            r_beta_141 = r_beta(2*a,elements_data[row,8],4) #Brzeg 4-1 wezly 1 5
        
        r_beta_1 = r_beta_112 + r_beta_123 + r_beta_134 + r_beta_141
        
        A[row,:] = elements_data[row,9]
        r_beta_sum[row,:] = np.transpose(r_beta_1)    
        k_k_sum[:,:,row] = k_k_1
        
        
        
        '''
        k_alfa_112 = k_alfa(2*b,elements_data[row,1],1) #Brzeg 1-2 wezly 5 6
        k_alfa_123 = k_alfa(2*a,elements_data[row,2],2)   #Brzeg 2-3 wezly 6 2
        k_alfa_134 = k_alfa(2*b,elements_data[row,3],3) #Brzeg 3-4 wezly 2 1
        k_alfa_141 = k_alfa(2*a,elements_data[row,4],4)   #Brzeg 4-1 wezly 1 5    
        
        k_k_1 = k_k_11 + k_alfa_112 + k_alfa_123 + k_alfa_134 + k_alfa_141
        
        r_beta_112 = r_beta(2*b,elements_data[row,5],1) # r_q(2*a,2*b,q) #Brzeg 1-2 wezly 5 6
        r_beta_123 = r_beta(2*a,elements_data[row,6],2) #Brzeg 2-3 wezly 6 2
        r_beta_134 = r_beta(2*b,elements_data[row,7],3) #Brzeg 3-4 wezly 2 1
        r_beta_141 = r_beta(2*a,elements_data[row,8],4) #Brzeg 4-1 wezly 1 5
        
        r_beta_1 = r_beta_112 + r_beta_123 + r_beta_134 + r_beta_141
        
        
        A[row,:] = elements_data[row,9]
        r_beta_sum[row,:] = np.transpose(r_beta_1)    
        k_k_sum[:,:,row] = k_k_1
        
       '''
        
    '''
    
    #Element 1 - wezly 5 6 2 1
    k_k_11 = k_k(a,b,k_x,k_y)
    
    k_alfa_123 = k_alfa(2*a,alfa,2)   #Brzeg 2-3 wezly 6 2
    k_alfa_134 = k_alfa(2*b,alfa,3) #Brzeg 3-4 wezly 2 1
    k_alfa_141 = k_alfa(2*a,alfa,4)   #Brzeg 4-1 wezly 1 5    
    k_k_1 = k_k_11 + k_alfa_123 + k_alfa_134 + k_alfa_141
    
    r_beta_112 = 0; # r_q(2*a,2*b,q)   #Brzeg 1-2 wezly 5 6
    r_beta_123 = r_beta(2*a,beta_2,2) #Brzeg 2-3 wezly 6 2
    r_beta_134 = r_beta(2*b,beta_2,3) #Brzeg 3-4 wezly 2 1
    r_beta_141 = r_beta(2*a,beta_2,4) #Brzeg 4-1 wezly 1 5
    
    r_beta_1 = r_beta_112 + r_beta_123 + r_beta_134 + r_beta_141
    
    
    A[0::] = [5,6,2,1]
    r_beta_sum[0::] = np.transpose(r_beta_1)    
    k_k_sum[:,:,0] = k_k_1
    
    print(A)
    print(r_beta_sum)
    print(k_k_sum)

    #Element 2 - wezly 7 8 4 3
    k_k_22 = k_k(a,b,k_x,k_y)
    k_alfa_223 = k_alfa(2*a,alfa,2)   #Brzeg 2-3 wezly 8 4
    k_alfa_234 = k_alfa(2*b,alfa,3)   #Brzeg 3-4 wezly 4 3
    k_alfa_241 = k_alfa(2*a,alfa,4)   #Brzeg 4-1 wezly 3 7
    
    k_k_2 = k_k_22 + k_alfa_223 + k_alfa_234 + k_alfa_241;
    
    r_beta_212 = 0; # r_q(2*a,2*b,q)   %Brzeg 1-2 wezly 7 8
    r_beta_223 = r_beta(2*a,beta_2,2) #Brzeg 2-3 wezly 8 4 
    r_beta_234 = r_beta(2*b,beta_2,3) #Brzeg 3-4 wezly 4 3
    r_beta_241 = r_beta(2*a,beta_2,4) #Brzeg 4-1 wezly 3 7
    
    r_beta_2 = r_beta_212 + r_beta_223 + r_beta_234 + r_beta_241;
    
    A[1::] = [7, 8, 4, 3]
    r_beta_sum[1::] = np.transpose(r_beta_2)   
    k_k_sum[:,:,1] = k_k_2
    
    
    # Element 3 - wezly 9 10 6 5 
    k_k_33 = k_k(a,b,k_x,k_y)
    k_alfa_323 = k_alfa(2*a,alfa,2)   #Brzeg 2-3 wezly 10 6 
    k_alfa_341 = k_alfa(2*a,alfa,4)   #Brzeg 4-1 wezly 5 9
    
    k_k_3 = k_k_33 + k_alfa_323 + k_alfa_341
    
    r_beta_312 = 0 # r_q(2*a,2*b,q)   #Brzeg 1-2 wezly 9 10
    r_beta_323 = r_beta(2*a,beta_2,2) #Brzeg 2-3 wezly 10 6 
    r_beta_334 = 0 # r_q(2*a,2*b,q)   #Brzeg 3-4 wezly 6 5
    r_beta_341 = r_beta(2*a,beta_2,4) #Brzeg 4-1 wezly 5 9
    
    r_beta_3 = r_beta_312 + r_beta_323 + r_beta_334 + r_beta_341;
    
    A[2::] = [9, 10, 6, 5]
    r_beta_sum[2::] = np.transpose(r_beta_3)
    k_k_sum[:,:,2] = k_k_3
    
    
    # Element 4 - wezly 14 15 8 7
    k_k_44 = k_k(a,b,k_x,k_y)
    k_alfa_423 = k_alfa(2*a,alfa,2)   #Brzeg 2-3 wezly 15 8 
    k_alfa_441 = k_alfa(2*a,alfa,4)   #Brzeg 4-1 wezly 7 14
    
    k_k_4 = k_k_44 + k_alfa_423 + k_alfa_441
    
    r_beta_412 = 0 # r_q(2*a,2*b,q)   #Brzeg 1-2 wezly 14 15
    r_beta_423 = r_beta(2*a,beta_2,2) #Brzeg 2-3 wezly 15 8
    r_beta_434 = 0 # r_q(2*a,2*b,q)   #Brzeg 3-4 wezly 8 7
    r_beta_441 = r_beta(2*a,beta_2,4) #Brzeg 4-1 wezly 7 14
    
    r_beta_4 = r_beta_412 + r_beta_423 + r_beta_434 + r_beta_441
    
    A[3::] = [14, 15, 8, 7]
    r_beta_sum[3::] = np.transpose(r_beta_4)
    k_k_sum[:,:,3] = k_k_4
    
    # Element 5 - wezly 16 17 10 9
    k_k_55 = k_k(a,b,k_x,k_y)
    k_alfa_541 = k_alfa(2*a,alfa,4) #Brzeg 4-1 wezly 9 16
    
    k_k_5 = k_k_55 + k_alfa_541;
    
    r_beta_512 = r_beta(2*b,beta_1,1) #Brzeg 1-2 wezly 16 17
    r_beta_523 = 0 # r_q(2*a,2*b,q)   #Brzeg 2-3 wezly 17 10
    r_beta_534 = 0 # r_q(2*a,2*b,q)   #Brzeg 3-4 wezly 10 9
    r_beta_541 = r_beta(2*a,beta_2,4) #Brzeg 4-1 wezly 9 16
    
    r_beta_5 = r_beta_512 + r_beta_523 + r_beta_534 + r_beta_541
    
    A[4::] = [16, 17, 10, 9]
    r_beta_sum[4::] = np.transpose(r_beta_5)
    k_k_sum[:,:,4] = k_k_5; 

    # Element 6 - wezly 17 18 11 10
    k_k_66 = k_k(a,b,k_x,k_y);
    k_alfa_634 = k_alfa(2*b,alfa,3); #Brzeg 3-4 wezly 11 10
    
    k_k_6 = k_k_66 + k_alfa_634;
    
    r_beta_612 = r_beta(2*b,beta_1,1) #Brzeg 1-2 wezly 17 18
    r_beta_623 = 0 # r_q(2*a,2*b,q)   #Brzeg 2-3 wezly 18 11
    r_beta_634 = r_beta(2*b,beta_2,3) #Brzeg 3-4 wezly 11 10
    r_beta_641 = 0 # r_q(2*a,2*b,q)   #Brzeg 4-1 wezly 10 17
    
    r_beta_6 = r_beta_612 + r_beta_623 + r_beta_634 + r_beta_641
    
    A[5::]  = [17, 18, 11, 10]
    r_beta_sum[5::]  = np.transpose(r_beta_6)
    k_k_sum[:,:,5] = k_k_6
    
    
    # Element 7 - wezly 18 19 12 11
    k_k_77 = k_k(a,b,k_x,k_y)
    k_alfa_734 = k_alfa(2*b,alfa,3); #Brzeg 3-4 wezly 12 11
    
    k_k_7 = k_k_77 + k_alfa_734
    
    r_beta_712 = r_beta(2*b,beta_1,1) #Brzeg 1-2 wezly 18 19
    r_beta_723 = 0 # r_q(2*a,2*b,q)   #Brzeg 2-3 wezly 19 12
    r_beta_734 = r_beta(2*b,beta_2,3) #Brzeg 3-4 wezly 12 11
    r_beta_741 = 0 # r_q(2*a,2*b,q)   #Brzeg 4-1 wezly 11 18
    
    r_beta_7 = r_beta_712 + r_beta_723 + r_beta_734 + r_beta_741
    
    A[6::] = [18, 19, 12, 11]
    r_beta_sum[6::] = np.transpose(r_beta_7)
    k_k_sum[:,:,6] = k_k_7
    
    
    #Element 8 - wezly 19 20 13 12
    k_k_88 = k_k(a,b,k_x,k_y)
    k_alfa_834 = k_alfa(2*b,alfa,3) #Brzeg 3-4 wezly 13 12
    
    k_k_8 = k_k_88+ k_alfa_834
    
    r_beta_812 = r_beta(2*b,beta_1,1) #Brzeg 1-2 wezly 19 20
    r_beta_823 = 0 # r_q(2*a,2*b,q)   #Brzeg 2-3 wezly 20 13
    r_beta_834 = r_beta(2*b,beta_2,3) #Brzeg 3-4 wezly 13 12
    r_beta_841 = 0 # r_q(2*a,2*b,q)   #Brzeg 4-1 wezly 12 19
    
    r_beta_8 = r_beta_812 + r_beta_823 + r_beta_834 + r_beta_841
    
    A[7::] = [19, 20, 13, 12];
    r_beta_sum[7::] = np.transpose(r_beta_8)
    k_k_sum[:,:,7] = k_k_8; 
    
    
    # Element 9 - wezly 20 21 14 13
    k_k_99 = k_k(a,b,k_x,k_y)
    k_alfa_934 = k_alfa(2*b,alfa,3) #Brzeg 3-4 wezly 14 13
    
    k_k_9 = k_k_99+ k_alfa_934
    
    r_beta_912 = r_beta(2*b,beta_1,1) #Brzeg 1-2 wezly 20 21
    r_beta_923 = 0 # r_q(2*a,2*b,q)   #Brzeg 2-3 wezly 21 14
    r_beta_934 = r_beta(2*b,beta_2,3) #Brzeg 3-4 wezly 14 13
    r_beta_941 = 0 # r_q(2*a,2*b,q)   #Brzeg 4-1 wezly 13 20
    
    r_beta_9 = r_beta_912 + r_beta_923 + r_beta_934 + r_beta_941
    
    A[8::] = [20, 21, 14, 13]
    r_beta_sum[8::] = np.transpose(r_beta_9)
    k_k_sum[:,:,8] = k_k_9
    
    
    # Element 10 - wezly 21 22 15 14
    k_k_1010 = k_k(a,b,k_x,k_y)
    k_alfa_1023 = k_alfa(2*a,alfa,2) #Brzeg 2-3 wezly 22 15
    
    k_k_10 = k_k_1010 + k_alfa_1023
    
    r_beta_1012 = r_beta(2*b,beta_1,1) #Brzeg 1-2 wezly 21 22
    r_beta_1023 = r_beta(2*a,beta_2,2) #Brzeg 2-3 wezly 22 15
    r_beta_1034 = 0 # r_q(2*a,2*b,q)   #Brzeg 3-4 wezly 15 14
    r_beta_1041 = 0 # r_q(2*a,2*b,q)   #Brzeg 4-1 wezly 14 21
    
    r_beta_10 = r_beta_1012 + r_beta_1023 + r_beta_1034 + r_beta_1041
    
    A[9::] = [21, 22, 15, 14];
    r_beta_sum[9::] = np.transpose(r_beta_10)
    k_k_sum[:,:,9] = k_k_10;
    
    '''
    
    #Temperatures
    nodes_temperature = nodes_temperature(nodes_q, elements_q, A, r_beta_sum, k_k_sum) 
    elements_temperature = elements_temperature(a, b, A, nodes_temperature)
    
    #Visualisation - Seaborn Heatmap chart    
    n = 1
    nodes = np.zeros((n,n,elements_q))
    x1 = 0.01
    print(x1)
    y = symbols('y')
    x = symbols('x')

    max_val = elements_temperature[0].subs([(x,0),(y,0)])
    min_val = max_val
    
    
    for k in range(elements_q):
        for i in range(n):
            for j in range(n):
                nodes[n-i-1,j,k]=elements_temperature[k].subs([(x,x1),(y,x1)])
                t=nodes[n-i-1,j,k]
                if t > max_val:
                    max_val = t
                if t < min_val:
                    min_val = t
    
    print(nodes)
    
    heatmap=np.zeros((3*n,6*n))
    heatmap[0,0] = nodes[:,:,0]
    heatmap[0,5] = nodes[:,:,1]
    heatmap[1,0] = nodes[:,:,2]
    heatmap[1,5] = nodes[:,:,3]
    heatmap[2,0:6] = nodes[:,:,4:10]
     
    ax = sns.heatmap(heatmap, mask = heatmap < 0.01, xticklabels = [], yticklabels = [] ,cmap = 'Blues', vmin=50, vmax=85,linewidths=.5, annot = heatmap, square = True, linecolor='black', cbar_kws={'label': 'Temperatura [°C] '}).set_title('Rozkład 2D temperatury w elemencie')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    