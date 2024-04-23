from .PhotonicCristal import *
import numpy as np
import matplotlib.pyplot as plt

def dispersion(v:float,theta:float,A:Layer,B:Layer) -> float:
    
    def MAT() -> float:
        
        def M_t() -> float:
            
            def M_a():
                
                def na_2():
                    return A.e * A.u
                
                def Q_a():
                    return (np.pi / 300) * v * (A.a/2) * np.sqrt(na_2() - np.sin(theta)**2)
                
                return np.array([[np.cos(1*Q_a()), A.u*A.a/(2 * Q_a()) * np.sin(1 * Q_a()) ],
                  [-2 * Q_a()/(A.u * A.a) * np.sin(Q_a()) , np.cos(1*Q_a())]])
            
            def M_b():
                
                def nb_2():
                    return B.e * B.u
                
                def Q_b():
                    return (np.pi / 300) * v * (B.a/2) * np.sqrt(nb_2() - np.sin(theta)**2)
                
                return np.array([[np.cos(1*Q_b()), A.u*A.a/(2 * Q_b()) * np.sin(1 * Q_b()) ],
                  [-2 * Q_b()/(B.u * B.a) * np.sin(Q_b()) , np.cos(1*Q_b())]])
            
            return np.dot(M_b() , M_a())

        return np.trace(M_t())

    if ((1/2)*MAT()) < -1 :
        return 1
    if  ((1/2)*MAT()) > 1 :
        return -1
    else:
        return (1/np.pi) * ((np.arccos((1/2)* MAT()))).real

        

def transmision(v:float,theta:float,A:Layer,B:Layer):

    def E():
        return np.cos(theta)
    
    def M_a():
                
        def na_2():
            return A.e * A.u
        
        def Q_a():
            return (np.pi / 300) * v * (A.a/2) * np.sqrt(na_2() - np.sin(theta)**2)
        
        return np.array([[np.cos(1*Q_a()), A.u*A.a/(2 * Q_a()) * np.sin(1 * Q_a()) ],
            [-2 * Q_a()/(A.u * A.a) * np.sin(Q_a()) , np.cos(1*Q_a())]])
            
    def M_b():
        
        def nb_2():
            return B.e * B.u
        
        def Q_b():
            return (np.pi / 300) * v * (B.a/2) * np.sqrt(nb_2() - np.sin(theta)**2)
        
        return np.array([[np.cos(1*Q_b()), B.u*B.a/(2 * Q_b()) * np.sin(1 * Q_b()) ],
            [-2 * Q_b()/(B.u * B.a) * np.sin(Q_b()) , np.cos(1*Q_b())]])

    
    MP = np.dot(M_b(),M_a())
    MT = np.linalg.matrix_power(MP,50)
    Tras = np.abs( (2 * E()) /
                E() * ((MT[0,0] + MT[1,1]) - MT[1,0] - E()**2 * MT[0,1]))**2
    Refl = np.abs((E() * ((-MT[1,1] + MT[0,0]) + MT[1,0])) - (E()**2 * MT[0,1]) / 
                (E() * (MT[1,1] + MT[0,0]) - MT[1,0])  - (E()**2 * MT[0,1]))**2 
    

    return Tras,Refl
    
      
    

def refractions(cristal:Cristal1D):
    '''
    Función para crear el arreglo de índices de refracción
    usamos las siguientes variables:
    Argumentos
    num = numero de repeticiones
    ni = indice de refracción del incidente
    nf = indice de refracción de salida
    na = indice de refracción de la capa a
    nb = indice de refracción de la capa b
    '''
    n = [cristal.ni]   
    for i in range(int(cristal.num/2)):
        n.append(cristal.A.e)
        n.append(cristal.B.e)
    n.append(cristal.nf)
    return np.array(n)

def sizes(cristal:Cristal1D,refraction_list:np.array):
    '''
    función para crear el arreglo de los tamaños de las capas
    Argumentos
    size_a = tamaño de la capa A
    size_b = tamaño de la capa B
    refraction_list = la lista con los índices de refracción, para usar como referencia
    '''
    d = refraction_list[1:-1]
    l = len(d)
    #nas = [size_a for i in range(l)]
    #nbs = [size_b for i in range(l)]
    n = []
    for i in range(int(cristal.num/2)):
        n.append(cristal.A.a)
        n.append(cristal.B.a)

    return np.array(n)

def positions(sizes:np.array):
    '''
    Función para crear el arreglo de las posiciones de cada capa
    Argumentos
    sizes = arreglo de los tamaños de las capas
    '''
    pos = [0]
    for i in range(len(sizes)):
        pos.append(pos[i] + sizes[i])
    return np.array(pos)


def angles(cristal:Cristal1D, refractions:np.array):
    '''
    función para el cálculo de los ángulos basado en los índices de refracción
    Argumentos:
    a_0 = ángulo inicial 
    refractions = arreglo que contiene los índices de refracción al interior del cristal
    '''
    a0 = cristal.a_0
    Ang = np.zeros(cristal.num+2)
    Ang[0] = a0
    for i in range(1,len(Ang)):
        Ang[i] = np.arcsin(refractions[i-1] / refractions[i] * np.sin(Ang[i-1]))

    
    return Ang[1:]

def k(Lambda:float,Angles:np.array,refractions:np.array):
    tabla = np.ones((len(Angles)))
    for i in range(len(tabla)):
        tabla[i] = tabla[i] * 2 * np.pi/Lambda * np.cos(Angles[i]) * refractions[i]        
    return tabla[:-1] 

def wave_vector_i(Lambda, a0, ni):
    return 2 * (np.pi/Lambda) * np.cos(a0) * ni

def wave_vector_f(Lambda, Angles, nf ):
    return 2 * (np.pi/Lambda) * np.cos(Angles[-1]) * nf

def dinamic_matrix(cristal:Cristal1D,refractions:np.array, angles:np.array):
    '''
    función para calcular la matriz dinámica de las capas internas
    refractions = lista que contiene los indices de refracción de las capas
    angles = lista que contiene los angulos de refracción en las capas
    num = número de capas
    '''
    md = []
    for i in range(cristal.num):
        inner_list = [[1, 1], [-refractions[i] * np.cos(angles[i]), refractions[i] * np.cos(angles[i])]]
        md.append(inner_list)
    return np.array(md)    

def dinamic_matrix_i(cristal:Cristal1D):
    mdi = [[1,1],[-cristal.ni*np.cos(cristal.a_0), cristal.ni*np.cos(cristal.a_0)]]
    return np.array(mdi)

def dinamic_matrix_f(cristal:Cristal1D,angles):
    mdf = [[1,1],[-cristal.nf* np.cos(angles[-1]), cristal.nf*np.cos(angles[-1])]]
    return np.array(mdf)


def propagation_matrix(w,angles:np.array,sizes:np.array,refractions:np.array,cristal:Cristal1D):
    mp = []

    for i in range(cristal.num):
        matrix = np.array([
            [np.exp(-1j * k(w,angles,refractions)[i] * sizes[i]), 0],
            [0, np.exp(1j * k(w,angles,refractions)[i] * sizes[i])]
        ])
        mp.append(matrix)

    return np.array(mp)

def pr_A(x,cristal:Cristal1D):
    k_x = k(x)
    return np.exp(-1j * k_x[0] * cristal.A.a) * (np.cos(k_x[1] * cristal.B.a) - 1j / 2 * np.sin(k_x[1] * cristal.B.a) * (k_x[1] / k_x[0] + k_x[0] / k_x[1]))

def pr_B(x,cristal:Cristal1D):
    k_x = k(x)
    return np.exp(1j * k_x[0] * cristal.A.a) * (-1j / 2 * np.sin(k_x[1] * cristal.B.a) * (k_x[1] / k_x[0] - k_x[0] / k_x[1]))

def pr_Cs(x,cristal:Cristal1D):
    k_x = k(x)
    return np.exp(-1j * k_x[0] * cristal.A.a) * (1j / 2 * np.sin(k_x[1] * cristal.B.a) * (k_x[1] / k_x[0] - k_x[0] / k_x[1]))

def pr_Ds(x,cristal:Cristal1D):
    k_x = k(x)
    return np.exp(1j * k_x[0] * cristal.A.a) * (np.cos(k_x[1] * cristal.B.a) + 1j / 2 * np.sin(k_x[1] * cristal.B.a) * (k_x[1] / k_x[0] + k_x[0] / k_x[1]))

def wave(x):
    return np.arccos(1 / 2 * (pr_A(x) + pr_Ds(x))) 

def U(x, m):
    return np.sin((m + 1) * wave(x)) / np.sin(wave(x))

def Ma(x, m):
    U_prev_1 = U(x, m - 1)
    U_prev_2 = U(x, m - 2)
    return np.array([[pr_A(x) * U_prev_1 - U_prev_2, pr_B(x) * U_prev_1],
                     [pr_Cs(x) * U_prev_1, pr_Ds(x) * U_prev_1 - U_prev_2]])

def mas(x,cristal:Cristal1D):
    return Ma(x, cristal.num_rep)

def Rss(x):
    ma_x = mas(x)
    return np.abs(ma_x[1, 0])**2 / np.abs(ma_x[0, 0])**2

def Tss(x):
    return 1 / np.abs(mas(x)[0, 0])**2

def masUC(x):
    return Ma(x, 1)

def X(x):
    return (1 / np.pi) * np.real(np.arccos((1 / 2) * np.trace(masUC(x))))