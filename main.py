import numpy as np
import matplotlib.pyplot as plt

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


    
    
    
    
    