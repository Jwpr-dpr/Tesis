from PhotonicCristal import *

def planeWaves1D(cristal:Cristal1D,numG:int):
    '''
    Argumentos:
    cristal(Cristal): Objeto de tipo cristal que representa el cristal fotónico 1D
    numG (float): Número de ondas planas a utilizar

    Salidas:
    plot: Estructura de bandas del cristal 
    '''

    l1 = cristal.A.a
    l2 = cristal.B.a

    eps1 = cristal.A.e
    eps2 = cristal.B.e

    a = l1+l2
    numG = 50

    G  = np.arange(-numG* 2 * np.pi/a, numG * 2* np.pi/a, 2* np.pi/a)
    G1 = np.arange(-numG* 2 * np.pi/a, numG * 2* np.pi/a, 2* np.pi/a)
    chi = np.zeros((len(G),len(G1)),dtype = 'complex_')

    for i in range(len(G)):
        for J in range(len(G1)):
            if(G[i] - G1[J]) == 0:
                chi[i,J] = 1/(l1+l2) * (1/eps1 * l1 + 1/eps2*l2)
            else:
                chi[i,J] = 1j /(l1+l2)/(G[i]-G1[J]) * (1/eps1 * np.exp(-1j * (G[i]-G1[J]) * l1) - 1) + 1/eps2 * (np.exp(-1j * (G[i]-G1[J]) * (l1+l2)) - np.exp(-1j * (G[i]-G1[J])*l1))
    
    M = np.zeros((len(G),len(G1)),dtype = 'complex_')
    ks = np.linspace(-2*np.pi/a, 2*np.pi/a,100)
    Ks = list()
    dispe = list()

    for k in range(len(ks)):
        for i in range(len(G)):
            for j in range(len(G1)):
                M[j,i] = chi[j,i] * (ks[k] + G1[j]) * (ks[k] + G[i])
        V = np.linalg.eig(M)[0]
        dispe.append(np.sqrt(np.sort(np.abs(V))) * a/2/np.pi)
        Ks.append(ks[k] * a/np.pi)

    dispersion = []
    points = []
    for i in range(len(Ks)):
        dispersion.append([Ks[i]] * len(dispe[i]))
        points.append(dispe[i])
        
    return dispersion, points
        