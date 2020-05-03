# Laminado de espesor t constituido por 4 laminas 30/90/90/-30

import numpy as np 

# DATOS 
X_t = 2225      # Resistencia a traccion MPA
X_c =1300       # Resistencia a compresion MPA
Y_t=61          # resistencia a traccion en direccion transversal MPA
Y_c =245        # resistencia a compresion en direccion transversa MPA
S =108.85       # resistencia a coratadura MPA
alpha_1 =-1E-6  #coef dilatacion en la direccion fibra
alpha_2=26E-6   #coef dilatacion en la direccion transversal
T_c =180        #temperatura de curado
T_uso =23       #temperatura de uso

#matrices Q 
#(1) Lamina de 30 

Q_1 =np.array([[90.8,27.69,45.56],
              [27.69,19.12,16.51],
              [45.56,16.51,30.41]])
#(2) y (3) laminas de 90
Q_2 =np.array([[8.44,2.53,0.],
              [2.53,151.8,0.],
              [0.,0.,5.25]])
Q_3 =Q_2 

#(4) lamina de -30
Q_4 =np.array([[90.8,27.69,-45.56],
              [27.69,19.12,-16.51],
              [-45.56,-16.51,30.41]])   
#Determinacion de las matrices A, B y D
Q =np.zeros((4,3,3))   #reunimos todas las Q en una unica
Q[0,:,:]=Q_1
Q[1,:,:]=Q_2
Q[2,:,:]=Q_3
Q[3,:,:]=Q_4

#calculo de A
z_v =np.linspace(-1/2,1/2,5,axis=0)    #vector de zs
dz =np.diff(z_v)  #vector de diferencias de z




def calculo_suma(Q,vector):
    '''calcula la suma de Q multiplicada por el vector'''
    matriz =np.zeros((3,3))
    for i in range(4):
        matriz += Q[i,:,:]*vector[i]
    return matriz

A = calculo_suma(Q,dz) #por t

#calculo de B
dz2 =np.diff(z_v**2)    #vector de diferencias al cuadrado.

B =1/2*calculo_suma(Q,dz2) # por t**2

#calculo de D

dz3 =np.diff(z_v**3)     #vector de diferencias al cubo

D =1/3*calculo_suma(Q,dz3)  #por t**3
                  
# Calculo de las deformacion asociadas a la temperatura 

def rotar_vector(v,ang):
    '''Devuelve la matriz M rotada con angulo ang en grados,
    transforma un vector en ejes de ortotropia a ejes del laminado'''
    ang =np.pi/180*ang  #transformacion a radianes
    M =np.array([[np.cos(ang)**2,np.sin(ang)**2,-2*np.sin(ang)*np.cos(ang)],
                [np.sin(ang)**2,np.cos(ang)**2,2*np.sin(ang)*np.cos(ang)],
                [np.sin(ang)*np.sin(ang),-np.sin(ang)*np.cos(ang),np.cos(ang)**2-np.sin(ang)**2]])
    v_rot =M @ v
    return v_rot 

# calculo de los coeficientes en ejes del laminado

angulos =np.array([30,90,90,-30])

alphas =np.zeros((3,4))
alpha_v =np.array([alpha_1,alpha_2,0])

for i in range(4):
    alphas[:,i]=rotar_vector(alpha_v,angulos[i])
print(alphas)

Delta_T =T_uso-T_c 
N =np.zeros((1,3))

for i in range(4):
    N +=Q[i,:,:]  @ alphas[:,i] *dz[i]*Delta_T
    
print(N)







