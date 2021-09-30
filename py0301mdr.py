
from __future__ import division
import numpy as np

#===============================================================
# FUNCIONES PARA P2D
#===============================================================
#---------------------------------------------------------------
# Función K_P2D: Matriz de rigidez pórtico bidimensional P2D
#---------------------------------------------------------------
def K_P2D(EA, EI, L, Angulo):
## Función que obtiene la matriz de rigidez 6x6, en coordenadas globales, del P2D
    ## EA = Producto EA de la barra
    ## EI = Producto EI de la barra
    ## L = Longitud de la barra
    ## Angulo = angulo de orientación de la barra, en grados sexagesimales
     Angulo = Angulo * np.pi / 180.
     
# Matriz de cambio de base de local a global  M_{L,G}:   
     M_LG = np.array([[np.cos(Angulo), -np.sin(Angulo), 0],
                      [np.sin(Angulo), np.cos(Angulo), 0],
                      [0,             0,             1]])

# Se expande a una 6x6. Se construye por bloques, el bloque fuera de la diagonal es una matriz 3x3 de ceros. Esto se genera con np.zeros((3,3))
     M_LG = np.vstack([np.hstack([M_LG, np.zeros((3,3))]),
                       np.hstack([np.zeros((3,3)), M_LG])] )
         
# Primero se inicia en cero: 
     K_L = np.zeros((6,6))

# Se rellenan los términos de la diagonal superior:
     K_L[0,:] = np.array([ EA / L, 0., 0., - EA / L, 0., 0.])
     K_L[1,1:] = np.array([12 * EI / L**3, 6* EI / L**2, 0, 
                     -12 * EI / L**3, 6 * EI / L**2])
     K_L[2,2:] = np.array([4 * EI / L, 0, -6 * EI / L**2, 2 * EI / L ])
     K_L[3,3:] = np.array([EA / L, 0, 0])
     K_L[4,4:] = np.array([12 * EI / L**3, -6 * EI / L**2])
     K_L[5,5]  = 4 * EI / L
     
# Se construye el resto de la matriz por simetría. Se suma K_local + traspuesta (K_local), en numpy la traspuesta se obtiene con .T
# Se resta la diagonal de k_local para que no aparezca 2 veces
     K_L = (K_L + K_L.T) -  np.diag ( np.diag(K_L) )
     
# Cambio de base a coordenadas globales:  K_G = M(L,G) * K_L * M(G,L) 
# El producto se realiza sabiendo que el producto matricial A*B = np.dot(A,B)
     K_G = np.dot(M_LG, np.dot(K_L , M_LG.T ) )      
    
# Salida de la función: la submatriz cambiada a globales
     return K_G


#---------------------------------------------------------------
# Función para definir una submatriz de un P2D
#---------------------------------------------------------------
def K_P2D_sub(ea, ei, lon, angle, icol, jcol):
## Esta función obtiene la submatriz 3x3 de k_lineal de un P2D cambiada a globales
    ## ea = Producto EA de la barra
    ## ei = Producto EI de la barra
    ## lon = Longitud de la barra
    ## angle = angulo de orientación de la barra, en grados sexagesimales
    ## icol = posición (puede ser 1 o 2) 
    ## jcol = posición (puede ser 1 o 2)
     if (icol not in ([1,2]) or jcol not in ([1,2]) ):
         print("Valores de icol, jcol fuera de rango")
         return
        
# Puesto que Python indexa los vectores y matrices comenzando con el cero, indicaremos las posiciones restando 1 a las posiciones que se numeran comenzando desde el 1     
     icol = icol -1
     jcol = jcol -1 
     
# Llamada a la función K_P2D para generar la matriz 6x6:
     K_globales = K_P2D(ea, ei, lon, angle)

# Extracción del bloque:
     K_Bloque = K_globales[icol * 3 : icol * 3 + 3 , jcol * 3 : jcol * 3 + 3]
     
     return K_Bloque


#===============================================================
# RESOLUCIÓN DEL PROBLEMA
# Uso de funciones para definir submatrices en P2D
#===============================================================

#---------------------------------------------------------------
# Definición de las barras
#---------------------------------------------------------------
Elast = 20E6 # 20 GPa = 20E+6 MPa, común para las dos barras

# Barra a
Aa = 0.1**2 # Sección cuadrada de 0.1 m de lado
Ia = 1./12 * 0.1**4 # Inercia en m^4
La = 5 # Longitud en metros
alfa_a = 30 # Ángulo de la barra: 30 grados

EA_a = Elast * Aa
EI_a = Elast * Ia


# Barra b
Ab = 0.2**2 # Sección cuadrada, de 0.2 m de lado
Ib = 1./12 * 0.2**4 # Inercia en m^4
Lb = 7 # Longitud en metros
alfa_b = 0 #Angulo de la barra b

EA_b = Elast * Ab
EI_b = Elast * Ib

#---------------------------------------------------------------
# Submatrices para el nudo 22
#---------------------------------------------------------------
# Para la barra a, el nudo 2 es el nudo final, por lo que se debe extraer la submatriz del bloque "2,2":
k22_a = K_P2D_sub(EA_a, EI_a, La, alfa_a, 2, 2) 

# Para la barra b, el nudo 2 es el nudo inicial, por lo que se debe extraer la submatriz del bloque "1,1":
k22_b = K_P2D_sub(EA_b, EI_b, Lb, alfa_b, 1, 1) 

# Mostramos estas submatrices
print("K22_a =") 
print(k22_a)
print("")
print("K22_b =")
print(k22_b)

#---------------------------------------------------------------
# Ensamblaje: K22_a + K22_b
#---------------------------------------------------------------
k22 = k22_a + k22_b
print("")
print("K22_a + k22_b =")
print(k22)

#---------------------------------------------------------------
# Vector de fuerzas en el nodo 2
#---------------------------------------------------------------
fuerzas = [100, -200, 300]

#---------------------------------------------------------------
# Resolución del sistema de ecuaciones
# Obtención de desplazamientos
#---------------------------------------------------------------
u_2 = np.linalg.solve(k22, fuerzas)
print("")
print("Solución del sistema ")
print("u_2x =", u_2[0], "m")
print("u_2y =", u_2[1], "m")
print("theta_2 =", u_2[2], "rad")