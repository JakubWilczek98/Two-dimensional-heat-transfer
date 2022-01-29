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
        # print('r_beta bound 1-2')
        # print(r)
        
    elif element_number == 2: #bound 2-3
        r[1] = k
        r[2] = k
        # print('r_beta bound 2-3')
        # print(r)
        
    elif element_number == 3: #bound 3-4
        r[2] = k
        r[3] = k
        # print('r_beta bound 3-4')
        # print(r)
        
    elif element_number == 4: #bound 4-1
        r[3] = k
        r[0] = k
        # print('r_beta bound 4-1')
        # print(r_beta,2)
    
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
        # print('k_alfa bound 1-2')
        # print(np.round(k_alfa,2))
        
    elif element_number == 2: #bound 2-3
        k_alfa[1,1] = a1
        k_alfa[1,2] = a2
        k_alfa[2,1] = a2
        k_alfa[2,2] = a1
        # print('k_alfa bound 2-3')
        # print(np.round(k_alfa,2))
        
    elif element_number == 3: #bound 3-4
        k_alfa[2,2] = a1
        k_alfa[2,3] = a2
        k_alfa[3,2] = a2
        k_alfa[3,3] = a1
        # print('k_alfa bound 3-4')
        # print(np.round(k_alfa,2))
        
    elif element_number == 4: #bound 4-1
        k_alfa[0,0] = a1
        k_alfa[0,3] = a2
        k_alfa[3,0] = a2
        k_alfa[3,3] = a1
        # print('k_alfa bound 4-1')
        # print(np.round(k_alfa,2))
        
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
    
    matrix_g_latex = pd.DataFrame(np.round(matrix_global,2)).to_latex()
    print(matrix_g_latex)
    
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
    h = 85
    T_inf = 28
    k_x = 55
    k_y = 55
    nodes = 0
    q = 0
    beta_1 = 6000
    beta_2 = h * T_inf
    alfa = -h
    a = 0.02/2
    b = a
    
    #Number of nodes and elements
    elements_q = 40
    nodes_q = 63

    #Target matrixes 
    A = np.zeros((elements_q,4))
    r_beta_sum = np.zeros((elements_q,4))
    k_k_sum = np.zeros((4,4,elements_q))
   
    #Boundary conditions for elements
    #Thickened mesh
    
    
    elements_data = np.array([[0,0,0,alfa,alfa,0,0,beta_2,beta_2,[7,8,2,1]],
                                [1,0,alfa,alfa,0,0,beta_2,beta_2,0,[8,9,3,2]],
                                [2,0,0,alfa,alfa,0,0,beta_2,beta_2,[10,11,5,4]],
                                [3,0,alfa,alfa,0,0,beta_2,beta_2,0,[11,12,6,5]],
                                [4,0,0,0,alfa,0,0,0,beta_2,[13,14,8,7]],
                                [5,0,alfa,0,0,0,beta_2,0,0,[14,15,9,8]],
                                [6,0,0,0,alfa,0,0,0,beta_2,[16,17,11,10]],
                                [7,0,alfa,0,0,0,beta_2,0,0,[17,18,12,11]],
                                [8,0,0,0,alfa,0,0,0,beta_2,[19,20,14,13]],
                                [9,0,alfa,0,0,0,beta_2,0,0,[20,21,15,14]],
                                [10,0,0,0,alfa,0,0,0,beta_2,[22,23,17,16]],
                                [11,0,alfa,0,0,0,beta_2,0,0,[23,24,18,17]],
                                [12,0,0,0,alfa,0,0,0,beta_2,[25,26,20,19]],
                                [13,0,alfa,0,0,0,beta_2,0,0,[26,27,21,20]],
                                [14,0,0,0,alfa,0,0,0,beta_2,[35,36,23,22]],
                                [15,0,alfa,0,0,0,beta_2,0,0,[36,37,24,23]],
                                [16,0,0,0,alfa,0,0,0,beta_2,[38,39,26,25]],
                                [17,0,0,0,0,0,0,0,0,[39,40,27,26]],
                                [18,0,0,alfa,0,0,0,beta_2,0,[40,41,28,27]],
                                [19,0,0,alfa,0,0,0,beta_2,0,[41,42,29,28]],
                                [20,0,0,alfa,0,0,0,beta_2,0,[42,43,30,29]],
                                [21,0,0,alfa,0,0,0,beta_2,0,[43,44,31,30]],
                                [22,0,0,alfa,0,0,0,beta_2,0,[44,45,32,31]],
                                [23,0,0,alfa,0,0,0,beta_2,0,[45,46,33,32]],
                                [24,0,0,alfa,0,0,0,beta_2,0,[46,47,34,33]],
                                [25,0,0,alfa,0,0,0,beta_2,0,[47,48,35,34]],
                                [26,0,0,0,0,0,0,0,0,[48,49,36,35]],
                                [27,0,alfa,0,0,0,beta_2,0,0,[49,50,37,36]],
                                [28,0,0,0,alfa,beta_1,0,0,beta_2,[51,52,39,28]],
                                [29,0,0,0,0,beta_1,0,0,0,[52,53,40,39]],
                                [30,0,0,0,0,beta_1,0,0,0,[53,54,41,40]],
                                [31,0,0,0,0,beta_1,0,0,0,[54,55,42,41]],
                                [32,0,0,0,0,beta_1,0,0,0,[55,56,43,42]],
                                [33,0,0,0,0,beta_1,0,0,0,[56,57,44,43]],
                                [34,0,0,0,0,beta_1,0,0,0,[57,58,45,44]],
                                [35,0,0,0,0,beta_1,0,0,0,[58,59,46,45]],
                                [36,0,0,0,0,beta_1,0,0,0,[59,60,47,46]],
                                [37,0,0,0,0,beta_1,0,0,0,[60,61,48,47]],
                                [38,0,0,0,0,beta_1,0,0,0,[61,62,49,48]],
                                [39,0,alfa,0,0,beta_1,beta_2,0,0,[62,63,50,49]]])
    '''
    #Basic mesh
    elements_data = np.array([[0,0,alfa,alfa,alfa,0,beta_2,beta_2,beta_2,[5,6,2,1]],
                                [1,0,alfa,alfa,alfa,0,beta_2,beta_2,beta_2,[7,8,4,3]],
                                [2,0,alfa,0,alfa,0,beta_2,0,beta_2,[9, 10, 6, 5]],
                                [3,0,alfa,0,alfa,0,beta_2,0,beta_2,[14, 15, 8, 7]],
                                [4,0,0,0,alfa,beta_1,0,0,beta_2,[16, 17, 10, 9]],
                                [5,0,0,alfa,0,beta_1,0,beta_2,0,[17, 18, 11, 10]],
                                [6,0,0,alfa,0,beta_1,0,beta_2,0,[18, 19, 12, 11]],
                                [7,0,0,alfa,0,beta_1,0,beta_2,0,[19, 20, 13, 12]],
                                [8,0,0,alfa,0,beta_1,0,beta_2,0,[20, 21, 14, 13]],
                                [9,0,alfa,0,0,beta_1,beta_2,0,0,[21, 22, 15, 14]]])
    ''' 
    
        
        
        
        
        
    for row in range(len(elements_data)):
        k_k_11 = k_k(a,b,k_x,k_y)
        
        if elements_data[row,1] == 0:
            k_alfa_112 = 0
        else:
            k_alfa_112 = k_alfa(2*b,elements_data[row,1],1) #Bound 1-2 nodes 5 6
        if elements_data[row,2] == 0:
            k_alfa_123 = 0
        else:
            k_alfa_123 = k_alfa(2*a,elements_data[row,2],2)   #Bound 2-3 nodes 6 2
        
        if elements_data[row,3] == 0:
            k_alfa_134 = 0
        else:
            k_alfa_134 = k_alfa(2*b,elements_data[row,3],3) #Bound 3-4 nodes 2 1
        
        if elements_data[row,4] == 0:
            k_alfa_141 = 0
        else:
            k_alfa_141 = k_alfa(2*a,elements_data[row,4],4)   #Bound 4-1 nodes 1 5    
               
        k_k_1 = k_k_11 + k_alfa_112 + k_alfa_123 + k_alfa_134 + k_alfa_141
        
        # print(np.round(k_k_1,2))
        
        if elements_data[row,5] == 0:
            r_beta_112 = 0
        else:
            r_beta_112 = r_beta(2*b,elements_data[row,5],1) # r_q(2*a,2*b,q) #Bound 1-2 nodes 5 6
        
        if elements_data[row,6] == 0:
            r_beta_123 = 0
        else:
            r_beta_123 = r_beta(2*a,elements_data[row,6],2) #Bound 2-3 nodes 6 2
        
        if elements_data[row,7] == 0:
            r_beta_134 = 0
        else:
            r_beta_134 = r_beta(2*b,elements_data[row,7],3) #Bound 3-4 nodes 2 1
        
        if elements_data[row,8] == 0:
            r_beta_141 = 0
        else:
            r_beta_141 = r_beta(2*a,elements_data[row,8],4) #Bound 4-1 nodes 1 5
        
        r_beta_1 = r_beta_112 + r_beta_123 + r_beta_134 + r_beta_141
        
        # print(np.round(r_beta_1,2))
        
        A[row,:] = elements_data[row,9]
        r_beta_sum[row,:] = np.transpose(r_beta_1)    
        k_k_sum[:,:,row] = k_k_1
        
    
    #Temperatures
    nodes_temperature = nodes_temperature(nodes_q, elements_q, A, r_beta_sum, k_k_sum) 
    elements_temperature = elements_temperature(a, b, A, nodes_temperature)
    
    #print(np.round(nodes_temperature,2))
    #print(elements_temperature)
    
    #Visualisation - Seaborn Heatmap chart    
    n = 1
    nodes = np.zeros((n,n,elements_q))
    x1 = 0.02
    y = symbols('y')
    x = symbols('x')

    max_val = elements_temperature[0].subs([(x,0),(y,0)])
    min_val = max_val
    
    #print(elements_temperature)
    
    for k in range(elements_q):
        for i in range(n):
            for j in range(n):
                nodes[n-i-1,j,k]=elements_temperature[k].subs([(x,x1),(y,x1)])
                t=nodes[n-i-1,j,k]
                if t > max_val:
                    max_val = t
                if t < min_val:
                    min_val = t
    
    #print(nodes)
    
    
    
    heatmap=np.zeros((6*n,12*n))
    heatmap[0,0:2] = nodes[:,:,0:2]
    heatmap[0,10:12] = nodes[:,:,2:4]
    heatmap[1,0:2] = nodes[:,:,4:6]
    heatmap[1,10:12] = nodes[:,:,6:8]
    heatmap[2,0:2] = nodes[:,:,8:10]
    heatmap[2,10:12] = nodes[:,:,10:12]
    heatmap[3,0:2] = nodes[:,:,12:14]
    heatmap[3,10:12] = nodes[:,:,14:16]
    heatmap[4,0:12] = nodes[:,:,16:28]
    heatmap[5,0:12] = nodes[:,:,28:40]
    
     
    ax = sns.heatmap(heatmap, mask = heatmap < 0.01, xticklabels = [], yticklabels = [] ,cmap = 'Blues', vmin=44, vmax=75,linewidths=.5, annot = heatmap, square = True, linecolor='black', cbar_kws={'label': 'Temperatura [°C] '}).set_title('Rozkład 2D temperatury w elemencie')
    '''
    
    
    
    heatmap=np.zeros((3*n,6*n))
    heatmap[0,0] = nodes[:,:,0]
    heatmap[0,5] = nodes[:,:,1]
    heatmap[1,0] = nodes[:,:,2]
    heatmap[1,5] = nodes[:,:,3]
    heatmap[2,0:6] = nodes[:,:,4:10]
     
    ax = sns.heatmap(heatmap, mask = heatmap < 0.01, xticklabels = [], yticklabels = [] ,cmap = 'Blues', vmin=40, vmax=65,linewidths=.5, annot = heatmap, square = True, linecolor='black', cbar_kws={'label': 'Temperatura [°C] '}).set_title('Rozkład 2D temperatury w elemencie')
    
    '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    