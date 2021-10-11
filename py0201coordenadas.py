# -*- coding: utf-8 -*-
"""
Autores:
Alejandro E. Martínez Castro, R. Gallego, E. Puertas, F. Ávila. 
amcastro@ugr.es
Licencia Creative-Commons CC BY-NC-ND 3.0 ES
"""
#===============================================================
#Librerías y funciones necesarias
#===============================================================
import numpy as np

#===============================================================
#Clases y funciones para matriz de rigidez A3D
#Programación orientada a objetos
#===============================================================

#---------------------------------------------------------------
#Definición de la clase "Punto". Coordenadas (x, y, z) de un punto
#    Observación: Se inicializa con las coordenadas (0,0,0)
#    Ejemplo: p1 = Punto() generará un punto p1 de coordenadas (0,0,0)
#             p1 = Punto(1,2,3) generará un punto de coordenadsa (1,2,3)
#    Uso de la función coords()
#          Muestra por pantalla el punto, en formato (x,y,z)
#    Ejemplo: p1 = Punto(1,2,3)
#              p1.coords()
#             (Se muestra en pantalla (1,2,3) )
#---------------------------------------------------------------

class Punto:
    """ Clase para representar los puntos en 3D, coordenadas x, y, z """
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y 
        self.z = z
    def coords(self):
        return "({0}, {1},{2})".format(self.x, self.y, self.z)

#---------------------------------------------------------------
#Definición de la función "longitud" 
#     Esta función calcula la distancia entre dos puntos
#     Notas : += es una función interesante, heredada de C
#          Ejemplo: a = a + b se escribe a += b
#     Esta función trabaja con p1,p2. Es una función "amiga" de la clase Punto
#--------------------------------------------------------------
def longitud(p1, p2):
    """ Devuelve la longitud (módulo) del vector p1-p2 """
    Longitud = 0
    Longitud += (p2.x - p1.x ) ** 2
    Longitud += (p2.y - p1.y ) ** 2
    Longitud += (p2.z - p1.z ) ** 2
    Longitud = np.sqrt(Longitud)
    return Longitud

#---------------------------------------------------------------
# Definición de la clase A3D, para barra articulada 3D
#---------------------------------------------------------------
class A3D:
    """ Clase para representar la barra articulada 3D """
    def __init__(self, p1 = Punto(), p2 = Punto(), EA = 0):
        self.p1 = p1
        self.p2 = p2
        self.ea = EA
        
#Función para obtener la matriz de rigidez A3D en locales (2x2)
    def mrig_local(self):
        L = longitud(self.p1,self.p2)
        K_L = np.array([[self.ea / L, -self.ea / L],
                     [-self.ea / L, self.ea / L]])
        return K_L

#Función para obtener la matriz de cambio de base de Local a Global
#Para matrices A3D
    def m_lg(self):
        L = longitud (self.p1, self.p2)   #Se usa la función definida antes "longitud"     
        #Definición del vector unitario de dirección p1 apuntando a p2
        vx = (self.p2.x - self.p1.x) / L 
        vy = (self.p2.y - self.p1.y) / L
        vz = (self.p2.z - self.p1.z) / L
        #Montaje de la matriz de cambio de base de Local a Global. 
        #Recuerde: Por columnas, se indican los vectores del sistema Local, referenciados en el sistema Global.         
        M_LG = np.array([[vx, 0],
                         [vy, 0],
                         [vz, 0],
                         [0,  vx],
                         [0,  vy],
                         [0,  vz]])
        return M_LG

#Función para obtener la matriz de rigidez A3D en globales (6x6)
    def mrig_global(self): 
        
        M_LG = self.m_lg() # Se llama a la función previamente creada
        K_L = self.mrig_local() # Se llama a otra función previamente creada
        #Producto matricial para construir la matriz de rigidez en globales.         
        K_G = np.dot(M_LG, np.dot(K_L , M_LG.T ))  
        
        return K_G

#Función para escribir por pantalla la matriz de rigidez A3D
    def imprime_mrig(self):
        np.set_printoptions(precision=3)
        print(self.mrig_global())


#===============================================================
# Datos de entrada, usamos la clase "Punto" creada         
#===============================================================
Punto1 = Punto(10,10,0)    
Punto2 = Punto(20,30,60)
Elast = 210E6 # kN/m2
Area = 2E-4 # 2 cm2, pasado a m2

EA = Elast * Area

#===============================================================
# Creación de un objeto de tipo A3D, a través de la clase A3D
#===============================================================
barraA = A3D(Punto1, Punto2, EA)

#===============================================================
# Generación de las matrices de rigidez y de cambio
#===============================================================
# Usaremos las funciones del objeto "barraA" de clase A3D. 
# Para usar una función de la clase A3D, se escribe nombredeobjeto.función(argumentos)

matrig = barraA.mrig_global()      # matriz de rigidez local
mat_lg = barraA.m_lg()             # matriz cambio local-global
matrig_local = barraA.mrig_local() # matriz de rigidez global

#===============================================================
# Impresión de resultados
#===============================================================
print("Matriz de rigidez en coord. locales")
print(matrig_local)
print("")

print("Matriz de cambio de base Local - Global")
print(mat_lg)
print("")

print("Matriz de rigidez en coord. globales")
barraA.imprime_mrig()
