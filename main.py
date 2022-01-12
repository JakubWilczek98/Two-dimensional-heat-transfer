import numpy as np
import matplotlib.pyplot as plt
from sympy import * 

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
        
    elif element_number == 3: #bound 2-3
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
    
    elements_q = 10
    nodes_q = 22
    
    A = np.zeros((elements_q,4))
    r_beta_sum = np.zeros((elements_q,4))
    k_k_sum = np.zeros((4,4,elements_q))
   

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
    k_k_sum[:,:,0] = k_k_1 #Tu moze byc błąd, wyżej sprawdzone!!!!!
    
    
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
    
    # Temperatura w poszczegolnych wezlach
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
         
    results = np.dot(np.linalg.inv(matrix_global),T)
    
    #Obliczanie temperatury w poszczegolnych elementach 
    
    y = symbols('y')
    x = symbols('x')
    
    s = x - a
    t = y - b
    
    N1 = ((a-s)*(b-t))/(4*a*b)

    N2 = ((a+s)*(b-t))/(4*a*b)
    N3 = ((a+s)*(b+t))/(4*a*b)
    N4 = ((a-s)*(b+t))/(4*a*b)
    
    N_t = np.array([N1, N2, N3, N4])
    
    
    Eq = zeros(elements_q,1)
    i = 0

    for i in range(elements_q):
        Temperature = [results[int(A[i,0])], results[int(A[i,1])], results[int(A[i,2])], results[int(A[i,3])]]
        q = np.transpose(N_t)*np.transpose(Temperature)
        Eq[i] = sum(sum(q))
        
    print(Eq)
    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    