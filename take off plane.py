import matplotlib.pyplot as plt
from numpy import *
import math
from fractions import Fraction
import numpy as np
pi = math.pi

global rho_0, g, dt, G, M_Terre, P_ext, T_air, M_molaire, coef, S_m

rho_0 = 1.2
g = 9.81
t = 0
dt = 0.1
R_Terre = 6378137
G = 6.674*10**(-11)
M_Terre = 5.972*10**24
P_ext = 101325
T_air = 300
R = 8.314
M_molaire = 28.965/1000
coef = 1

#S_m = 2.1*g/(rho_0*v_x*v_x)
#https://fr.wikipedia.org/wiki/Robin_DR-400 
#cz = S_m*650/14.8
v_z = 0
v_2 =v_z
z = 0.00001
z_2 = z
x = 0
T = []
S_m = 14.2/650*0.5

def vole(v_x,z_c,x_c,t_c):
	
	t = 0
	i = 1
	z = 0.00001
	z_2 = z
	v_z = 0
	v_2 = v_z
	x = 0
	v_c = sqrt(2*g*exp(M_molaire*g*z_c/(R*T_air))/(rho_0*S_m))
	
	#constante d automatisme
	K = 1/10  # K<2 sinon diverge
	Td = 60

	
	vitesse = [0]
	hauteur = [0]
	T = [0]
	speed = [0]
	altitude= [0]
	Pression = [P_ext]
	Pression_test = [P_ext]
	Li_Er = [0]
	v_X = [v_x]
	X = [0]
	acc_x = [0]
	acc_z = [0]
	
	while z > 0 and (x < x_c or t > t_c):
		
		T.append(t)
		
		
		#valeurs atmosphere
		P = P_ext*exp(-M_molaire*g*z/(R*T_air))
		rho = rho_0*P/P_ext
		Pression.append(P)
		
		
		#composante veeticale
		a_z = v_x*v_x*rho*S_m/2-g
		v_z = v_z + a_z*dt
		z = z + v_z*dt +a_z*dt*dt/2
		acc_z.append(a_z)
		vitesse.append(v_z)
		hauteur.append(z)
		
		
		#calcule de correction
		erreur = v_c - v_x - Td*(-v_z*v_x/2*(-M_molaire*g/(R*T_air)))
		Li_Er.append(erreur)
		
		
		#composante horizontale
		a_x = K*erreur #-v_x*v_x*rho*4/(2*650)
		v_x = v_x + a_x*dt#+ K*erreur #+ a_x*dt
		x = x + v_x*dt #+ a_x*dt*dt/2
		acc_x.append(a_x)
		v_X.append(v_x)
		X.append(x)
		
		
		t = t + dt
	
	K = 1/50
	Td = 100
	
	while z > 5 and t < 4*t_c:
		
		T.append(t)
		
		#valeurs atmosphere
		P = P_ext*exp(-M_molaire*g*z/(R*T_air))
		rho = rho_0*P/P_ext
		Pression.append(P)
		
		#composante veeticale
		a_z = v_x*v_x*rho*S_m/2-g
		v_z = v_z + a_z*dt
		z = z + v_z*dt +a_z*dt*dt/2
		acc_z.append(a_z)
		vitesse.append(v_z)
		hauteur.append(z)
		
		#calcule de correction
		erreur = sqrt(2*g/(rho_0*S_m)) - v_x - Td*(-v_z*v_x/2*(-M_molaire*g/(R*T_air)))
		Li_Er.append(erreur)
		
		
		#composante horizontale
		a_x = + K*erreur #-v_x*v_x*rho*4/(2*650)
		v_x = v_x + a_x*dt
		x = x + v_x*dt #+ a_x*dt*dt/2
		acc_x.append(a_x)
		v_X.append(v_x)
		X.append(x)
		
		
		t = t + dt
		
	
	return T,vitesse,hauteur,Pression,Li_Er,X,v_X,acc_x,acc_z

v_x = 40
z_c = 1000 #round(-R*T_air/(M_molaire*g)*log(2*g/(S_m*rho_0*v_x*v_x)))
x_c = 20*1000
t_c = 500


print("L'objectif est d'atteindre",z_c,'m',"en partant à",v_x*3.6,'km/h')
print("la charge alaire est de",round(S_m,5),"m²/kg")
print(sqrt(2*g/(rho_0*S_m)))



liste = vole(v_x,z_c,x_c,t_c)
T = liste[0]
v_z = liste[1]
hauteur = liste[2]
Pression = liste[3]
Li_Er = liste[4]
X = liste[5]
v_x = liste[6]
a_x = liste[7]
a_z = liste[8]

l_g = [g]*len(X)

#liste_2 = normal(v_x)
#altitude = liste_2[4]

j = 5
plt.figure(1,figsize=(13,10))

plt.subplot(3,1,1)
plt.title("altitude")
plt.plot(X,hauteur,"-g")
plt.ylabel('m')
plt.subplot(3,1,2)
plt.title("vitesse verticale")
plt.plot(T,v_z,"-g")
plt.ylabel('m/s')
plt.subplot(3,1,3)
plt.title("acceleration verticale")
plt.plot(T[:j],a_z[:j],"-g",label="acceleration")
plt.plot(T[:j],l_g[:j],"-r",label="attraction gravitationnelle")
plt.legend(loc="lower right")
plt.ylabel('m/s²')


#plt.subplot(4,1,4)
#plt.title("correction de la vitesse horizontale")
#plt.plot(T,Li_Er)
#plt.xlabel('temps')
#plt.ylabel('')

plt.figure(2,figsize=(13,10))
#plt.subplot(3,1,1)
#plt.title("position")
#plt.plot(T,X,"-g")
#plt.ylabel('m')
plt.subplot(2,1,1)
plt.title("vitesse horizontale")
plt.plot(T,v_x,"-g")
plt.ylabel('m/s')
plt.subplot(2,1,2)
plt.title("acceleration horizontale")
plt.plot(T,a_x,"-g")
plt.ylabel('m/s²')
#plt.subplot(2,1,2)
#plt.title("vitesse verticale")
#plt.plot(T,speed,"-g",label="avec correction")
#plt.legend(loc="lower right")
#plt.ylabel('m/s')

#plt.subplot(2,1,1)
#plt.title("vitesse horizontale")
#plt.plot(T,v_z,"-g",label="avec correction")
#plt.legend(loc="lower right")
#plt.ylabel('m/s')

#plt.figure(3,figsize=(13,10))

plt.show()
