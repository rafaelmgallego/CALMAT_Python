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

import numpy as np              #Para importar la librería numpy
import matplotlib.pyplot as plt #Para dibujar

#===============================================================
#Función para crear la matriz de rigidez de un elemento de 6 gdl
#Pórtico 2D (P2D)
#===============================================================
def k_p2d(ea, ei, lon):
    k_local = np.zeros([6,6])

    #Se rellenan primero los términos de la triangular superior
    #En Python las matrices y vectores se indexan desde el 0,1,2,3...
    #Una matriz 6x6 tendrá índices desde el 0,1,2,3,4,5.
    k_local[0,:] = np.array([ ea / lon, 0., 0., -ea / lon, 0., 0.])
    k_local[1,1:] = np.array([12 * ei / lon**3, 6* ei / lon**2, 0,
                     -12 * ei / lon**3, 6 * ei / lon**2])
    k_local[2,2:] = np.array([4 * ei / lon, 0, -6 * ei / lon**2, 2 * ei / lon ])
    k_local[3,3:] = np.array([ea / lon, 0, 0])
    k_local[4,4:] = np.array([12 * ei / lon**3, -6 * ei / lon**2])
    k_local[5,5]  = 4 * ei / lon

    #El resto de la matriz se construye por simetría.
    #Se suma k_local + traspuesta (k_local)
    #Se resta la diagonal de k_local para que no aparezca 2 veces
    #Para escribir la traspuesta de k_local, simplemente se escribe k_local.T
    k_local = (k_local + k_local.T) -  np.diag ( np.diag(k_local) )

    #Salida de la función:
    return k_local

#===============================================================
#Datos de entrada
#===============================================================
#Sistema de unidades: [kN], [m]
#Datos de un IPN-200

Elast = 210E6       #210 GPa = 210E6 kN/m2
Area  =  33.5E-4    #33.5 cm^2 = 33.5E-4 m^2
Inercia  = 2140E-8  #2140 cm^4 = 2140E-8 m^4
Longitud = 2.8      #Longitud en metros

#===============================================================
#Datos de desplazamientos conocidos en los extremos
#===============================================================
#Datos del nudo inicial i
ui_x = 0
ui_y = 0
Theta_i = 0

#Datos del nudo final j
uj_x = 1E-3
uj_y = 1E-3
Theta_j = 1E-3

#===============================================================
#Matriz de rigidez
#===============================================================

EA = Elast * Area
EI = Elast * Inercia

Mrig = k_p2d(EA, EI, Longitud)

print("")
print("Matriz de Rigidez")
print(Mrig)

#Para que la salida por pantalla sea más clara,
#fijamos la salida con 3 cifras decimales
np.set_printoptions(precision=3)
print("")
print("Matriz de Rigidez")
print(Mrig)

#===============================================================
#Vector desplazamientos
#===============================================================
u = np.array([ui_x,ui_y,Theta_i,uj_x,uj_y,Theta_j])

#===============================================================
#Vector de fuerzas
#===============================================================
#Producto Matriz de rigidez * Vector de desplazamientos
#(observe que no es necesario escribir u como matriz-columna):
P_local = np.dot(Mrig,u)

print("")
print("Vector de fuerzas")
print(P_local)

#Asignación de las componentes del vector de fuerza:
Pi_x = P_local[0]
Pi_y = P_local[1]
M_i  = P_local[2]
Pj_x = P_local[3]
Pj_y = P_local[4]
Mj   = P_local[5]

#Salida por pantalla:
print("")
print("Pi_x = ", P_local[0])
print("Pi_y = ", P_local[1])
print("M_i  = ", P_local[2])
print("Pj_x = ", P_local[3])
print("Pj_y = ", P_local[4])
print("Mj   = ", P_local[5])

#===============================================================
#Estudio paramétrico
#===============================================================
Longitudes = np.linspace(0.5, 5, 50) #Vector desde 0.5 hasta 5, con 50 puntos intermedios
VecPi_y = []

for Longitud in Longitudes:
    Mrigidez = k_p2d(EA, EI, Longitud)
    P = np.dot(Mrigidez,u)
    VecPi_y = np.append(VecPi_y,P[1]) #La posición [1] es la segunda, pues la primera es P[0]

print("")
print("Longitud,     Pi_y")
print("==================" )
for i in range(len(Longitudes)):
    print(Longitudes[i], ",", VecPi_y[i])


line, = plt.plot(Longitudes, VecPi_y, label = "P_{iy}", linewidth = 2.0)


plt.xlabel("Longitud de viga (m)")   #Etiqueta del eje x
plt.ylabel("kN")                     #Etiqueta del eje y

plt.legend(loc=0)  #Mostrar leyenda en posición óptima (loc=0)
plt.show()         #Mostrar figura final
