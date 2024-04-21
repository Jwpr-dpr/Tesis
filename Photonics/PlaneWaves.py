from PhotonicCristal import *
from scipy.linalg import eigh


def planeWaves1D_dispersion(cristal:Cristal1D,numG:int):
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
    
    G  = np.arange(-numG* 2 * np.pi/a, numG * 2* np.pi/a, 2* np.pi/a)
    G1 = np.arange(-numG* 2 * np.pi/a, numG * 2* np.pi/a, 2* np.pi/a)
    chi = np.zeros((len(G),len(G1)),dtype = 'complex_')

    for i in range(len(G)):
        for J in range(len(G1)):
            if(G[i] - G1[J]) == 0:
                chi[i,J] = 1/(l1+l2) * (1/eps1 * l1 + 1/eps2*l2)
            else:
                chi[i,J] = 1j /(l1+l2)/(G[i]-G1[J]) * ((1/eps1 * (np.exp(-1j * (G[i]-G1[J]) * l1) - 1)) + (1/eps2 * (np.exp(-1j * (G[i]-G1[J]) * (l1+l2)) - np.exp(-1j * (G[i]-G1[J])*l1))))

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

def planeWaves1D_field(cristal:Cristal2D,numG):
    
    l1 = cristal.A.a
    l2 = cristal.B.a

    eps1 = cristal.A.e
    eps2 = cristal.B.e

    a = l1+l2

    G  = np.arange(-numG* 2 * np.pi/a, numG * 2* np.pi/a, 2* np.pi/a)
    G1 = np.arange(-numG* 2 * np.pi/a, numG * 2* np.pi/a, 2* np.pi/a)
    chi = np.zeros((len(G),len(G1)),dtype = 'complex_')

    for i in range(len(G)):
        for J in range(len(G1)):
            if(G[i] - G1[J]) == 0:
                chi[i,J] = 1/(l1+l2) * (1/eps1 * l1 + 1/eps2*l2)
            else:
                chi[i,J] = 1j /(l1+l2)/(G[i]-G1[J]) * ((1/eps1 * (np.exp(-1j * (G[i]-G1[J]) * l1) - 1)) + (1/eps2 * (np.exp(-1j * (G[i]-G1[J]) * (l1+l2)) - np.exp(-1j * (G[i]-G1[J])*l1))))

    
    def calculate_field(x, a, numG, k, Vs, in1):
        H = 0
        countt = 0
        for G in np.arange(-numG*2*np.pi/a, numG*2*np.pi/a, 2*np.pi/a):
            countt += 1
            H += Vs[countt,in1]*np.exp(1j*(k+G)*x)
        return H

    
    k = np.pi/(15*a)
    in1 = 20

    countG1 = 1
    countG = 1
    M1 = np.zeros((2*numG+1, 2*numG+1), dtype=complex)

    for G in np.arange(-numG*2*np.pi/a, numG*2*np.pi/a, 2*np.pi/a):
        for G1 in np.arange(-numG*2*np.pi/a, numG*2*np.pi/a, 2*np.pi/a):
            M1[countG1-1,countG-1] = chi[countG1-1,countG-1]*(k+G1)*(k+G)
            countG += 1
        countG1 += 1
        countG = 1

    countG1 = 1

    eigvals, eigvecs = eigh(M1)
    ind = np.argsort(eigvals)
    Ds = eigvals[ind]
    Vs = eigvecs[:,ind]

    freq = np.sqrt(abs(Ds[in1-1]))*a/2/np.pi

    div = 15
    div1 = 1500
    x = np.linspace(0, a, div1+1)
    FieldE = np.abs(calculate_field(x, a, numG, k, Vs, in1))**2

    alt = np.max(np.abs(FieldE))
    x1 = np.linspace(0, alt, div+1)
    y1 = l1*np.ones(div+1)

    plt.figure(2)
    plt.plot(y1, x1, linewidth=5)
    plt.plot(x, FieldE, linewidth=2)
    plt.xlabel('x(m)', fontsize=20)
    plt.ylabel('E (V/m)', fontsize=20)
    #plt.ylim([0, alt])
    plt.xlim([0, a])
    plt.show()

    return
        