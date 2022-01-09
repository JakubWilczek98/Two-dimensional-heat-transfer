import numpy as np
import matplotlib.pyplot as plt

def r_beta(length, beta, element_number):
    k = length * beta
    r = np.zeros((4,1))
    
    if element_number == "1": #bound 1-2
        r[0] = k
        r[1] = k
        
    if element_number == "2": #bound 2-3
        r[1] = k
        r[2] = k
        
    if element_number == "3": #bound 3-4
        r[2] = k
        r[3] = k
        
    if element_number == "4": #bound 4-1
        r[3] = k
        r[1] = k
    
    else:
        print("R_beta Error")
    
    return r
        
        
def k_alfa(length, constatn, element_number):
    a1 = (-2)*length*constant/3
    a2 = (-length)*constant/3
    k_alfa = np.zeros((4,4))

    if element_number == "1": #bound 1-2
        k_alfa[0,0]
        k_alfa[0,1]
        k_alfa[1,0]
        k_alfa[1,1]
        
    if element_number == "2": #bound 2-3
        k_alfa[1,1]
        k_alfa[1,2]
        k_alfa[2,1]
        k_alfa[2,2]
        
    if element_number == "3": #bound 2-3
        k_alfa[2,2]
        k_alfa[2,3]
        k_alfa[3,2]
        k_alfa[3,3]
        
    if element_number == "4": #bound 4-1
        k_alfa[0,0]
        k_alfa[0,3]
        k_alfa[3,0]
        k_alfa[3,3]
        
    return k_alfa

def r_q(a,b,q):
    r = np.ones((4,1))
    r_q = a*b*q*r
    
    return r_q



    