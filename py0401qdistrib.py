from __future__ import division 
import numpy as np              

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
# Cálculo de estructura P2D - A2D
#===============================================================
#---------------------------------------------------------------
# Definición de los nudos
#---------------------------------------------------------------
p1 = Punto2D(0,0)
p2 = Punto2D(5,7)
p3 = Punto2D(11,7)
p4 = Punto2D(11,0)

#---------------------------------------------------------------
# Definición de las barras
#---------------------------------------------------------------
# Barra a
E_a = 20E6  #20 kN/m^2
Aa = 0.1**2 # Sección cuadrada de 0.1 m de lado
Ia = 1./12 * 0.1**4

EA_a = E_a * Aa
EI_a = E_a * Ia

barra_a = P2D(p1,p2,EA_a, EI_a)

# Barra b

E_b = 25E6    # 25 kN/m^2
Ab = 0.2**2   # Sección cuadrada de 0.2 m de lado
Ib = 1./12 * 0.2**4

EA_b = E_b * Ab
EI_b = E_b * Ib

barra_b = P2D(p2,p3,EA_b, EI_b)

# Barra c 
E_c = 200E6    # 200 kN/m^2
Ac = np.pi * 0.05**2 / 4. # Barra circular de 5 cm de diámetro

EA_c = E_c * Ac

barra_c = A2D(p2,p4,EA_c)


#---------------------------------------------------------------
# Determinación de las submatrices
#---------------------------------------------------------------
# De forma simbólica, es posible eliminar el bloque correspondiente al nudo 3, pues está empotrado. 
# También es posible eliminar el bloque del nudo 4, pues están impedidos los desplazamientos en los 2 gdl que tiene la barra A2D.
# Del sistema restante, almacenaremos en memoria las submatrices necesarias.

# Barra a
# La barra "a" conecta los nudos 1 y 2, por tanto las posiciones 11, 12 y 22, se encuentran en las posiciones 11, 12, y 22 de la matriz de rigidez en globales
k11_a = barra_a.K_global_caja(1,1)
k12_a = barra_a.K_global_caja(1,2)
k22_a = barra_a.K_global_caja(2,2)


# Barra b
# La barra "b" conecta los nudos 2 y 3. La submatriz k22_b es la única necesaria, y estará en la posición local 11 de la matriz de rigidez de la barra, en globales. 
k22_b = barra_b.K_global_caja(1,1)


# Barra c
# La barra "c" conecta los nudos 2 y 4. Se necesita la submatriz k22_c, que se encuentran en la posición 11 de la matriz de rigidez en globales
k22_c = barra_c.K_global_caja(1,1)

# La barra "c" es un tirante, y necesita expandirse a una 3x3 con ceros para poderla combinar con las matrices de P2D
# k22_c debe expandirse a una matriz de dimensión 3x3, con la última fila y columna de ceros.
k22_c_expandida = np.zeros([3,3]) # Primero se crea una matriz con ceros
k22_c_expandida[0:2,0:2] = k22_c # Luego se copian las posiciones de k22_c


# Suma de submatrices. Sólo lo necesita la posición k22_{a+b+c}
k22_abc = k22_a + k22_b + k22_c_expandida

#---------------------------------------------------------------
# Ensamblaje
#---------------------------------------------------------------
# La matriz de rigidez completa tiene dimensiones 6x6:
krigidez = np.zeros([6,6])

krigidez[0:3,0:3] = k11_a
krigidez[0:3,3:6] = k12_a
krigidez[3:6,0:3] = k12_a.T # Matriz transpuesta de k12_a
krigidez[3:6,3:6] = k22_abc 

print("")
print("Matriz de rigidez")
print(krigidez)

# En el nudo 1 hay un apoyo deslizante, que impide los desplazamientos en los gdl 2 y 3, por lo que es necesario eliminar las filas y columnas 2 y 3 de la matriz de rigidez
# Para Python, estas filas son las 1 y 2, ya que Python comienza indexando en 0
# Esto puede hacerse con la operación de numpy "delete"

posiciones = np.array([2,3]) # Posiciones indexando desde 1
posiciones -= 1 # Posiciones de Python, indexando desde 0
# Esta operación equivale a posiciones = posiciones - [1,1]
# Primero las filas (axis = 0):
krigidez_final = np.delete(krigidez,posiciones,axis = 0)
# Luego, las columnas(axis=1):
krigidez_final = np.delete(krigidez_final,posiciones,axis = 1) 

# Mostramos la nueva matriz de rigidez:
print("")
print("Matriz de rigidez tras eliminar filas-columnas 2 y 3")
print(krigidez_final)


#---------------------------------------------------------------
# Vector de fuerzas de empotramiento
#---------------------------------------------------------------
# Aparecen fuerzas de empotramiento en los nudos 2 y 3 debido a la carga distribuida en la barra b.
# Sólo interviene la fuerza de empotramiento del nudo 2

# Calculamos el vector de fuerzas de empotramiento del nudo 2, debido a la carga q:
q = 100  # kN/m
Lb = longitud(p2,p3)
p2 = np.array([0, q * Lb / 2, q * Lb**2 / 12.]) 

vec_fuerzas = np.zeros(6)
vec_fuerzas[3:6] = p2
vec_fuerzas = np.delete(vec_fuerzas,posiciones)

# Cambio de signo para tener en cuenta que el sistema es K * u = F - F_empotramiento:
vec_fuerzas = -vec_fuerzas 

#---------------------------------------------------------------
## Resolución del sistema de ecuaciones
#---------------------------------------------------------------
u = np.linalg.solve(krigidez_final,vec_fuerzas)
print("")
print("Solución del sistema. Desplazamientos: ")
print("u_1x =", u[0], "m")
print("u_2x =", u[1], "m")
print("u_2y =", u[2], "m")
print("Theta_2 =", u[3], "rad")