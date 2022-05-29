import matplotlib.pyplot as plt
from numpy import *
import math
from fractions import Fraction
import numpy as np
pi = math.pi

global rho_0, g, dt, G, M_Terre, P_ext, T_air, M_molaire, coef, S_m

rho_0 = 1.2
g = 9.81
v_x = 45
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

S_m = 2.1*g/(rho_0*v_x*v_x)
#https://fr.wikipedia.org/wiki/Robin_DR-400 
cz = S_m*650/14.8
v_z = 0
v_2 =v_z
z = 0.00001
z_2 = z
x = 0
T = []



print(round(-R*T_air/(M_molaire*g)*log(2*g/(S_m*rho_0*v_x*v_x))))


print("L'altitude pour une vitesse horizontal de",v_x,'m/s est de',round(-R*T_air/(M_molaire*g)*log(2*g/(S_m*rho_0*v_x*v_x)),0),'m')
print("L'avion va donc à",v_x*3.6,"km/h, la charge alaire est de",round(S_m,5)*1000,"m²/gramme")

def normale(v_x):
	
	t = 0
	z = 0.00001
	z_2 = z
	v_z = 0
	v_2 = v_z
	x = 0
	
	vitesse = [0]
	hauteur = [0]
	temps = [0]
	speed = [0]
	altitude= [0]
	Pression = [P_ext]
	Pression_test = [P_ext]
	Li_Er = [0]
	X = [v_x]
	A_x = [0]
	
	while t < 1000 :#and z > 0 and z_2 > 0:
		
		a_z = v_x*v_x*rho_0*S_m/2-g
		v_z = v_z + a_z*dt
		z = z + v_z*dt +a_z*dt*dt/2
		
		T.append(t)
		vitesse.append(v_z)
		hauteur.append(z)
		
		
		
		#g = -G*M_Terre/(R_Terre+z)**2
		P = P_ext*exp(-M_molaire*g*z_2/(R*T_air))
		rho = rho_0*P/P_ext
		
		
		#test pour viabilité du DL de exp(z)
		P_test = P_ext*(1+(-M_molaire*g*z_2/(R*T_air)))
		erreur = (P-P_test)/P_ext
		Li_Er.append(erreur)
		rho_test = rho_0*P_test/P_ext
		
		a_2 = v_x*v_x*rho*S_m/2-g
		v_2 = v_2+a_2*dt
		z_2 = z_2 + v_2*dt + a_2*dt*dt/2
		
		speed.append(v_2)
		altitude.append(z_2)
		Pression.append(P)
		Pression_test.append(P_test)
		
		a_x = 0.0001*sqrt(2*g*rho_0/S_m)*M_molaire*g*v_2*exp(M_molaire*g*z_2/(2*R*T_air))
		v_x = v_x 
		x = x + v_x*dt #+ a_x*dt*dt/2
		A_x.append(a_x)
		X.append(v_x)
		
		t = t + dt
	return T,vitesse,speed,hauteur,altitude,Pression,Pression_test,Li_Er,A_x,X

def correction(v_x):
	
	t = 0
	i = 1
	z = 0.00001
	z_2 = z
	v_z = 0
	v_2 = v_z
	x = 0
	z_c = round(-R*T_air/(M_molaire*g)*log(2*g/(S_m*rho_0*v_x*v_x)),0)
	v_c = v_x
	
	#constante d automatisme
	K = 1  # K<2 sinon diverge
	Td = 20

	
	vitesse = [0]
	hauteur = [0]
	T = [0]
	speed = [0]
	altitude= [0]
	Pression = [P_ext]
	Pression_test = [P_ext]
	Li_Er = [0]
	v_X = [v_x]
	A_x = [0]
	
	while t < 1000 and z > 0 :
		
		T.append(t)
		vitesse.append(v_z)
		hauteur.append(z)
		
		#g = -G*M_Terre/(R_Terre+z)**2
		P = P_ext*exp(-M_molaire*g*z/(R*T_air))
		rho = rho_0*P/P_ext
		
		a_z = v_x*v_x*rho*S_m/2-g
		v_z = v_z + a_z*dt
		z = z + v_z*dt +a_z*dt*dt/2
		
		#test pour viabilité du DL de exp(z)
		#P_test = P_ext*(1+(-M_molaire*g*z_2/(R*T_air)))
		erreur = v_c - v_x - Td*(-v_z*v_x/2*(-M_molaire*g/(R*T_air)))
		Li_Er.append(erreur)
		#rho_test = rho_0*P_test/P_ext
		
		#a_2 = v_x*v_x*rho*S_m/2-g
		#v_2 = v_2+a_2*dt
		#z_2 = z_2 + v_2*dt + a_2*dt*dt/2
		
		#speed.append(v_2)
		#altitude.append(z_2)
		Pression.append(P)
		
		a_x = 0
		v_x = v_x + K*erreur
		x = x + v_x*dt #+ a_x*dt*dt/2
		A_x.append(x)
		v_X.append(v_x)
		
		t = t + dt
		i=+ 1
	return T,vitesse,hauteur,Pression,Li_Er,A_x,v_X

liste_2 = normale(v_x)
altitude = liste_2[4]
vitesse = liste_2[2]
v_2 = liste_2[9]

liste = correction(v_x)
T = liste[0]
speed = liste[1]
v_z = liste[6]
hauteur = liste[2]
Pression_test = liste[3]
Li_Er = liste[4]
X = liste[5]

#liste_2 = normal(v_x)
#altitude = liste_2[4]


plt.figure(1,figsize=(15,10))

plt.subplot(2,1,1)
plt.title("altitude sans correction de v_x")
plt.plot(X,altitude)

plt.subplot(2,1,2)
plt.title("altitude avec correction de v_x")
plt.plot(X,hauteur)

#plt.subplot(4,1,4)
#plt.title("correction de la vitesse horizontale")
#plt.plot(T,Li_Er)
#plt.xlabel('temps')
#plt.ylabel('')

plt.figure(2,figsize=(15,10))

plt.subplot(3,1,2)
plt.title("vitesse verticale")
plt.plot(T,speed,"-g",label="avec correction")
plt.plot(T,vitesse,"-r",label="sans correction")
plt.legend(loc="lower right")
plt.ylabel('m/s')

plt.subplot(3,1,1)
plt.title("vitesse horizontale")
plt.plot(T,v_z,"-g",label="avec correction")
plt.plot(T,v_2,"-r",label="sans correction")
plt.legend(loc="lower right")
plt.ylabel('m/s')

plt.subplot(3,1,3)
plt.title("valeur de la correction")
plt.plot(T,Li_Er)
plt.xlabel('temps')

plt.show()
