import matplotlib.pyplot as plt
from numpy import *
import math
from fractions import Fraction
import numpy as np
pi = math.pi

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
v_z = 0
v_2 =v_z
z = 0.00001
z_2 = z
x = 0
T = []
vitesse = []
hauteur = []
temps = []
speed = []
altitude= []
Pression = []
Pression_test = []
Li_Er = []
X = []
A_x = []

print("l'altitude pour une vitesse horizontal de",v_x,'m/s est de',round(-R*T_air/(M_molaire*g)*log(2*g/(S_m*rho_0*v_x*v_x)),0),'m')


while t < 1000 and z > 0 and z_2 > 0:
	
	a_z = v_x*v_x*rho_0*S_m/2-g
	v_z = v_z + a_z*dt
	z = z + v_z*dt +a_z*dt*dt/2
	
	T.append(t)
	vitesse.append(v_z)
	hauteur.append(z)
	
	
	
	#g = -G*M_Terre/(R_Terre+z)**2
	P = P_ext*exp(-M_molaire*g*z_2/(R*T_air))
	P_test = P_ext*(1+(-M_molaire*g*z_2/(R*T_air)))
	rho = rho_0*P/P_ext
	
	erreur = (P-P_test)/P_ext
	
	Li_Er.append(erreur)
	
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




plt.subplot(4,1,1)
plt.title("altitude avec rho constant")
plt.plot(T,hauteur)

plt.subplot(4,1,2)
plt.title("altitude avec rho non constant")
plt.plot(T,altitude)

plt.subplot(4,1,3)
plt.title("Pression test")
plt.plot(T,Pression_test)

plt.subplot(4,1,4)
plt.title("Pression")
plt.plot(T,Li_Er)

plt.show()
