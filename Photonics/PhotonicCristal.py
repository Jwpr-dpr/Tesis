import numpy as np
import matplotlib.pyplot as plt

class Layer:
    '''Clase que representa una capa en un cristal fotónico'''
    def __init__(self,a:float,e:float,u:float):
        '''
        Constructor de la clase Layer.

        Args:
            e (float): Valor de permitividad.
            u (float): Valor de permeabilidad.
            a (float): Ancho de la capa.
        '''
        self.a = a
        self.e = e
        self.u = u
        return

    def ni_2(self):
        return self.e * self.u
    
    
    def M_i(self,v,theta):
        def ni_2():
            return self.e * self.u

        def Q_i():
            return (np.pi/300) * v * (self.a/2) * np.sqrt(ni_2() - np.sin(theta)**2)       
        
        m = np.array([[np.cos(1*Q_i()), self.u_a*self.a/(2 * Q_i()) * np.sin(1 * Q_i()) ],
                    [-2 * Q_i()/(self.u_a * self.a) * np.sin(Q_i()) , np.cos(1*Q_i())]])
        return m

    def M_t(v,theta,M_a,M_b):
        return np.dot(M_b(v,theta) , M_a(v,theta))

class Cylinder:
    '''
    Clase que representa un cilindro sólido.
    '''
    def __init__(self,r:float,eps:float,u:float) -> None:
        '''
        Args:
        r(float): radio del cilindro
        eps(float): permitividad del cilindro
        u(float): permeabilidad del cilindro
        '''
        self.r = r
        self.eps = eps
        self.u = u
        pass

class Cube:
    '''
    Clase que representa un cubo sólido
    '''
    def __init__(self,x:float,y:float,eps:float) -> None:
        '''
        Args:
        x(float): Tamaño en x del cubo
        y(float): Tamaño en y del cubo
        eps(float): permitividad del cubo 
        '''
        self.x = x
        self.y = y
        self.eps = eps
        pass

class Cristal1D:
    '''
    Clase que representa un cristal fotónico 1D.
    La celda unitaria está compuesta por dos capas con distinto índice de refracción.
    '''
    def __init__(self,num_rep:int,ni:float,nf:float,a_0:float,A:Layer,B:Layer) -> None:
        '''
        Inicializa la celda unitaria y las repéticiones que tendrá
        Args:
        num_rep(int): Número de repeticiones periódicas de la celda unitaria
        ni(float): 
        nf(float):
        a_0(float):
        A(layer): Objeto de tipo Layer, es la primera capa de la celda unitaria
        B(layer): Objeto de tipo Layer, es la segunda capa de la celda unitaria
        '''
        self.num_rep = num_rep
        self.ni = ni
        self.nf = nf
        self.a_0 = a_0
        self.A = A
        self.B = B
        self.num = num_rep 
        pass
    def graficar_cristal(self):
    # Definir el tamaño del plot
        plt.figure(figsize=(8, 6))
        
        # Repetir la creación de los rectángulos Num veces a lo largo del eje x
        for i in range(self.num):
            # Calcular las coordenadas y dimensiones de los rectángulos
            rect_A = plt.Rectangle((i * (self.A.a + self.B.a), 0), self.A.a, self.A.a + self.B.a, color='blue')
            rect_B = plt.Rectangle(((i+1) * (self.A.a + self.B.a) - self.B.a, 0), self.B.a, self.A.a + self.B.a, color='green')
            
            # Agregar los rectángulos al plot
            plt.gca().add_patch(rect_A)
            plt.gca().add_patch(rect_B)
        
        # Establecer límites y etiquetas de los ejes
        plt.xlim(0, self.num * (self.A.a + self.B.a))
        plt.ylim(0, self.A.a + self.B.a)
        plt.xlabel('Repeticiones')
        plt.ylabel('Espesor')
        plt.title('Cristal Fotónico 1D')
        
        # Mostrar el plot
        
        plt.show()

class Cristal2D:
    '''
    clase que representa un cristal fotónico 2D.
    La celda unitaria consiste de una barra cilíndrica rodeada de otro material.
    '''
    def __init__(self,num_rep:int,bar:Cylinder,cube:Cube) -> None:
        '''
        Inicializa la celda unitaria y las repéticiones en cada dirección.
        Args:
        num_rep(int): Número de repeticiones periódicas de la celda unitaria.
        cube(Cube): Cubo de un material, base de la celda unitaria 
        bar(Cylinder): barra de un material que va inscrita dentro del cubo en la celda unitaria.
        '''
        self.num_rep = num_rep
        self.bar = bar
        self.cube = cube
        pass
