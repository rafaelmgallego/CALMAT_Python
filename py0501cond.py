# -*- coding: utf-8 -*-
"""
Autores:
Alejandro E. Martínez Castro, R. Gallego, E. Puertas, F. Ávila. 
amcastro@ugr.es
Licencia Creative-Commons CC BY-NC-ND 3.0 ES
"""

import numpy as np              
import matplotlib.pyplot as plt

#===============================================================
# DEFINICIÓN DE CLASES Y FUNCIONES 
# Matrices de Rigidez P2D y A2D en globales
# Completa y con submatrices
# Programación Orientada a Objetos
#===============================================================

#---------------------------------------------------------------
# Definición de la clase "Punto". Coordenadas (x, y) de un punto
#    Observación: Se inicializa con las coordenadas (0,0)
#    Se define bidimensional, pues se usará para el P2D
#---------------------------------------------------------------
class Punto2D:
    """ Clase para representar los puntos en 3D, coordenadas x, y, z """
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y 
    def coords(self):
        return "({0}, {1})".format(self.x, self.y)

#---------------------------------------------------------------
# Definición de la función "longitud" 
# Calcula la distancia entre dos puntos
# Trabaja con p1,p2. Es una función "amiga" de la clase Punto 
#---------------------------------------------------------------
def longitud(p1, p2):
    """ Devuelve la longitud (módulo) del vector p1-p2 """
    Longitud = 0
    Longitud += (p2.x - p1.x ) ** 2
    Longitud += (p2.y - p1.y ) ** 2
    Longitud = np.sqrt(Longitud)
    return Longitud

#---------------------------------------------------------------
# Definición de la función "alfa"
# Devuelve el ángulo que define el vector P1 -> P2
# La función "arctan2" devuelve el ángulo buscado
#---------------------------------------------------------------
def alfa(p1,p2):
    """Angulo P1 -> P2"""
    vy = p2.y -p1.y
    vx = p2.x -p1.x
    angulo = np.arctan2(vy,vx)
    return angulo

#---------------------------------------------------------------
# Definición de la clase "P2D", para barra P2D
#---------------------------------------------------------------
# Inicialización: A partir de dos puntos, y los parámetros EA, EI
# Funciones de la clase:  
#     K_local: Para obtener la matriz de rigidez en locales (6x6)
#     M_lg: Matriz de cambio de base de Locales a Globales (6x6)
#     K_global: Para obtener la matriz de rigidez en globales (6x6)
#     K_global_caja(i,j): Para obtener la submatriz ij en globales (3x3)

class P2D:
    """ Clase para representar la barra P2D """
    def __init__(self, p1 = Punto2D(), p2 = Punto2D(), ea = 0, ei = 0):
        self.p1 = p1
        self.p2 = p2
        self.ea = ea
        self.ei = ei

    def K_local(self):
        L = longitud(self.p1,self.p2)
        EA = self.ea
        EI = self.ei        
       
        K_L = np.zeros((6,6))

        K_L[0,:] = np.array([ EA / L, 0., 0., - EA / L, 0., 0.])
        K_L[1,1:] = np.array([12 * EI / L**3, 6* EI / L**2, 0, 
                     -12 * EI / L**3, 6 * EI / L**2])
        K_L[2,2:] = np.array([4 * EI / L, 0, -6 * EI / L**2, 2 * EI / L ])
        K_L[3,3:] = np.array([EA / L, 0, 0])
        K_L[4,4:] = np.array([12 * EI / L**3, -6 * EI / L**2])
        K_L[5,5]  = 4 * EI / L
        
        K_L = (K_L + K_L.T) -  np.diag ( np.diag(K_L) )
        
        return K_L

    def M_lg(self):
        Angulo = alfa (self.p1, self.p2)   
        
        M_LG = np.array([[np.cos(Angulo), -np.sin(Angulo), 0],
                      [np.sin(Angulo), np.cos(Angulo), 0],
                      [0,             0,             1]])
        M_LG = np.vstack([np.hstack([M_LG, np.zeros((3,3))]),
                       np.hstack([np.zeros((3,3)), M_LG])] )

        return M_LG

    def K_global(self): 
        
        M_LG = self.M_lg() 
        K_L = self.K_local() 
        
        K_G = np.dot(M_LG, np.dot(K_L , M_LG.T ))  
        
        return K_G
        
# Función para extraer una submatriz
    def K_global_caja(self,icol=1,jcol=1):
        if (icol not in ([1,2]) or jcol not in ([1,2]) ):
            print("Valores de icol, jcol fuera de rango")
            return
        icol = icol -1
        jcol = jcol -1 
        kg  = self.K_global()

        K_sub = kg[icol * 3:icol * 3 + 3 , jcol * 3 : jcol * 3 + 3]
     
        return K_sub


# Matriz de rigidez con giro condensado en globales
#     La operación de condensación se realiza en esta función

    def K_Condensa_Giro(self,nodo=2) : 
        # Si nodo = 1, giro condensado en nodo inicial
        # Si nodo = 2, giro condensado en nodo final
        indice = nodo*3-1 # Indice local del GDL en giro, en vector de Python
            
        kg = self.K_global()
        
#     k11,k12,k21,k22 son matrices-bloque necesarias para condensar el giro
#      k11: submatriz 5x5. A partir de kg, eliminando F-C del GDL condensado
#      k12: columna 5x1. En la posición del GDL condensado
#      k21: fila 1x5. En la posición del GDL condensado. 
#      k22: componente de la diagonal en el GDL condensado   
        k11 = np.delete(kg,indice,axis=0)
        k11 = np.delete(k11,indice,axis=1)
      # Se extrae la columna k12 en formato matriz 
      # Observe cómo se extraen en Python columnas, respetando el formato columna     
        k12 = kg.T[[indice]].T    
        k12 = np.delete(k12,indice,axis=0) #se elimina la fila "indice"
      # Se extrae la fila "indice" en formato matriz
        k21 = kg[[indice]] # Directamente extrae la fila "indice" en forma de matriz
        k21 = np.delete(k21,indice,axis=1) #se elimina la columa "indice"     
        k22 = kg[indice,indice]
        
        return k11 - k12.dot(k21) / k22 # k12.dot(k21)=np.dot(k12,k21)


#---------------------------------------------------------------
# Definición de la clase "A2D", para barra A2D
#---------------------------------------------------------------
# Inicialización: A partir de dos puntos, y el parámetro ea
# Funciones de la clase:  
#     K_local : Para obtener la matriz de rigidez en locales (2x2)
#     M_lg : Matriz de cambio de base de Locales a Globales (4x2)
#     K_global: Para obtener la matriz de rigidez en globales (6x6)
#     K_global_caja(i,j): Para obtener la submatriz en ij en globales (2x2)

     
class A2D:
    """ Clase para representar la barra A2D """
    def __init__(self, p1 = Punto2D(), p2 = Punto2D(), ea = 0):
        self.p1 = p1
        self.p2 = p2
        self.ea = ea

    def K_local(self):
        L = longitud(self.p1,self.p2)
        EA = self.ea     
        K_L = np.array([[EA / L, -EA / L],
                     [-EA / L, EA / L]])
        return K_L

    def M_lg(self):
        Angulo = alfa (self.p1, self.p2)   
        M_LG = np.array([[np.cos(Angulo), 0],
               [np.sin(Angulo), 0],
               [0, np.cos(Angulo)],
               [0, np.sin(Angulo)]])
        return M_LG

    def K_global(self): 
        M_LG = self.M_lg() 
        K_L = self.K_local() 
        K_G = np.dot(M_LG, np.dot(K_L , M_LG.T ))  
        return K_G
        
# Función para extraer una submatriz
    def K_global_caja(self,icol=1,jcol=1):
        if (icol not in ([1,2]) or jcol not in ([1,2]) ):
            print("Valores de icol, jcol fuera de rango")
            return
        icol = icol -1
        jcol = jcol -1 
        kg  = self.K_global()
        K_sub = kg[icol * 2:icol * 2 + 2 , jcol * 2 : jcol * 2 + 2]
        return K_sub


#===============================================================
# RESOLUCIÓN DEL PROBLEMA
# Cálculo de estructura P2D - A2D con condensación de gdl
#===============================================================
#---------------------------------------------------------------
# Definición de los parámetros geométricos
#---------------------------------------------------------------
L1 = 30 #m
L2 = 45 #m
h1 = 20 #m
d1 = 5  #m
h2 = 20 #m
d2 = 5  #m

# Definición de los puntos a partir de los parámetros geométricos

# Puntos 1,2,3,4, que definen el tablero
p1 = Punto2D(0,0)
p2 = Punto2D(L1,0)
p3 = Punto2D(L1+L2,0)
p4 = Punto2D(L1+L2+L1,0)

# Puntos 5,6, apoyos. Observe el uso de coordenadas relativas a los puntos 1 a 4
p5 = Punto2D(p2.x-d1, p2.y-h1)
p6 = Punto2D(p3.x+d2, p3.y-h2)

# Verificación de los puntos. Salida por pantalla
listapuntos = [p1,p2,p3,p4,p5,p6]

print("Listado de puntos")
npunto = 1
for punto in listapuntos:
    print("Punto", npunto, "Coordenadas", punto.coords())
    npunto +=1
print("")  
  
#---------------------------------------------------------------
# Definición de los parámetros mecánicos de las barras
#---------------------------------------------------------------
Elast = 35E6 #Módulo de elasticidad, kN/m2
At = 3.5     #Área sección tablero, m2
It = 0.12    #Inercia sección tablero, m4

Ap1 = 0.2    #Área pila 1, m2
Ip1 = 5.3E-3 #Inercia pila 1, m4

Ap2 = 0.2      #Área pila 2, m2
Ip2 = 5.3E-3   #Inercia pila 2, m4

# Parámetros de las barras:
EA_tablero = Elast * At
EI_tablero = Elast * It

EA_p1 = Elast * Ap1
EI_p1 = Elast * Ip1

EA_p2 = Elast * Ap2
EI_p2 = Elast * Ip2

#---------------------------------------------------------------
# Definición de las barras
#---------------------------------------------------------------
barra_a = P2D(p1,p2,EA_tablero,EI_tablero)
barra_b = P2D(p2,p3,EA_tablero,EI_tablero)
barra_c = P2D(p3,p4,EA_tablero,EI_tablero)
barra_d = P2D(p5,p2,EA_p1,EI_p1)
barra_e = P2D(p6,p3,EA_p2,EI_p2)

#---------------------------------------------------------------
# Dibujo aproximado de la estructura
#---------------------------------------------------------------
viga_x = [p1.x, p2.x, p3.x, p4.x]
viga_y = [p1.y, p2.y, p3.y, p4.y]
line, = plt.plot(viga_x, viga_y, label = "Viga", linewidth = 2.0)
pilar1_x = [p5.x, p2.x]
pilar1_y = [p5.y, p2.y]
line, = plt.plot(pilar1_x, pilar1_y, label = "Pilar1", linewidth = 2.0)
pilar2_x = [p6.x, p3.x]
pilar2_y = [p6.y, p3.y]
line, = plt.plot(pilar2_x, pilar2_y, label = "Pilar2", linewidth = 2.0)

plt.axes().set_aspect('equal','datalim')
plt.show()

#---------------------------------------------------------------
# Definición de las submatrices de rigidez necesarias
#---------------------------------------------------------------
# Barra a: Condensación del giro en 1
#   Observe la definición de la función K_Condensa_Giro.
#   El índice 1 hace referencia a que se condensa el giro del nodo inicial. Si hubiese sido 2 se condensa el del extremo. Por defecto condensa en 2, por eso en este caso hay que indicarlo. 

ka = barra_a.K_Condensa_Giro(1)
np.set_printoptions(precision=3)
print("")
print("Matriz de rigidez de la barra a, giro en 1 condensado")
print(ka)

# Extracción de la caja k22_a: Atención a los índices en Python
k22_a = ka[2:5,2:5]
# También se podía haber hecho directamente:
# k22_a = barra_a.K_Condensa_Giro(1)[2:5,2:5]


# Barra b: No se condensa nada
# Con la función específica de la clase se hace así: 
k22_b_clase = barra_b.K_global_caja(1,1)
k23_b_clase = barra_b.K_global_caja(1,2)
k33_b_clase = barra_b.K_global_caja(2,2)

# Se han designado "clase" para indicar que se han obtenido con la función K_global_caja()
# Sin la función de la clase, se hace así: 
k22_b = barra_b.K_global()[0:3,0:3]
k23_b = barra_b.K_global()[0:3,3:6]
k33_b = barra_b.K_global()[3:6,3:6]
 
# Observación: La operación de selección es inclusiva en el primer índice y exclusiva en el segundo 0:3 significa que tome los valores 0,1,2, y no el 3
# Si se comparan K22_b con K22_b_clase y siguientes, se observa que son iguales.


#  Barra c: Condensación del giro en el nodo final
kc = barra_c.K_Condensa_Giro(2) # No es necesario indicar (2) pues está por defecto
k33_c = kc[0:3,0:3]
k34_c = kc[0:3,3:5] #Caja de dimensiones 3x2
k44_c = kc[3:5,3:5] #Caja de dimensiones 2x2


# Barra d: Condensación del giro en el nodo final

kd = barra_d.K_Condensa_Giro() # Observe que por defecto toma el giro condensado en 2

# Se expande a una 3x3
k22_d = np.zeros([3,3]) 
k22_d[0:2,0:2] = kd[3:5,3:5]


# Barra e: Condensación del giro en el extremo
ke = barra_e.K_Condensa_Giro() # Observe que por defecto toma el giro condensado en 2
# Se expande a una 3x3
k33_e = np.zeros([3,3]) 
k33_e[0:2,0:2] = ke[3:5,3:5]


# Suma de matrices en nudos 2 y 3
k22 = k22_a + k22_b + k22_d
k33 = k33_b + k33_c + k33_e

#---------------------------------------------------------------
# Ensamblaje de la matriz de rigidez
#---------------------------------------------------------------

krigidez = np.zeros([8,8])

krigidez[0:3,0:3] = k22

krigidez[0:3,3:6] = k23_b
krigidez[3:6,0:3] = k23_b.T # Completamos con la traspuesta

krigidez[3:6,3:6] = k33

krigidez[3:6,6:8] = k34_c
krigidez[6:8,3:6] = k34_c.T #Completamos con la traspuesta

krigidez[6:8,6:8] = k44_c

#---------------------------------------------------------------
#  Vector de fuerzas de empotramiento
#---------------------------------------------------------------
q = 100 # Carga distribuida: 10 kN/m

La = longitud(p1,p2)
F2_a = np.array([0, 5./8 *q *La,-q*La**2./8.])

Lb = longitud(p2,p3)
F2_b = np.array([0, q * Lb / 2., q * Lb**2 / 12.])
F3_b = np.array([0, q * Lb / 2., -q * Lb**2 / 12.])

Lc = longitud(p3,p4)
F3_c = np.array([0, 5./8 *q*Lc, q * Lc**2/8.])
F4_c = np.array([0, 3./8.*q*Lc])

F2 = F2_a + F2_b
F3 = F3_b + F3_c

Fempotram = np.zeros(8)
Fempotram[0:3] = F2
Fempotram[3:6] = F3
Fempotram[6:8] = F4_c

#---------------------------------------------------------------
# Eliminación de filas y columnas para imponer u4y = 0
# Como Python indexa los vectores desde 0, este gdl es el 7
#---------------------------------------------------------------

krigidez_final = np.delete(krigidez,7,axis = 0)
krigidez_final = np.delete(krigidez_final,7,axis = 1)
Fuerzas = -1 * Fempotram # Cambio de signo porque el sistema es K u = F - Femp
Fuerzas = np.delete(Fuerzas,7,axis=0)

#---------------------------------------------------------------
## Resolución del sistema de ecuaciones
#---------------------------------------------------------------
u = np.linalg.solve(krigidez_final,Fuerzas)
print("")
print("Solución del sistema ")
print("u_2x =",   u[0], "m")
print("u_2y =",   u[1], "m")
print("Theta_2 =",u[2], "rad")
print("")
print("u_3x =",   u[3], "m")
print("u_3y =",   u[4], "m")
print("Theta_3 =",   u[5], "rad")
print("")
print("u_4x =",      u[6], "m")
