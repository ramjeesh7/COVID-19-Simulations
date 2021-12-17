print(f"COVID-19 model by Ramjee Sharma!")
print(f"No part of these codes will be reproduced and used for any purpose without the auther's permission.")
    
#Import Python libraries
import numpy as np
import matplotlib.pyplot as plt
import time
    
#3d 4th Order Runge-Kutta Function
def rk4_3d(x,y,z,dt):
    k1_1=dt*f(x,y,z)
    k2_1=dt*g(x,y,z) 
    k3_1=dt*h(x,y,z) 
        
    k1_2=dt*f(x+0.5*k1_1,y+0.5*k1_1,z+0.5*k1_1)
    k2_2=dt*g(x+0.5*dt*k2_1,y+0.5*k2_1,z+0.5*k2_1) 
    k3_2=dt*h(x+0.5*k3_1,y+0.5*k3_1,z+0.5*k3_1) 
        
    k1_3=dt*f(x+0.5*k1_2,y+0.5*k1_2,z+0.5*k1_2)
    k2_3=dt*g(x+0.5*k2_2,y+0.5*k2_2,z+0.5*k2_2) 
    k3_3=dt*h(x+0.5*k3_2,y+0.5*k3_2,z+0.5*k3_2) 
        
    k1_4=dt*f(x+k1_3,y+k1_3,z+k1_3) 
    k2_4=dt*g(x+k2_3,y+k2_3,z+k2_3) 
    k3_4=dt*h(x+k3_3,y+k3_3,z+k3_3)
        
    x=x+(1/6)*(k1_1+2*k1_2+2*k1_3+k1_4)
    y=y+(1/6)*(k2_1+2*k2_2+2*k2_3+k2_4)
    z=z+(1/6)*(k3_1+2*k3_2+2*k3_3+k3_4)
    return x,y,z
    
def f(S,I,R):#Right hand side of dS/dt
    return -beta*I*S
def g(S,I,R):#Right hand side of dI/dt
    return beta*I*S-kappa*I
def h(S,I,R):#Right hand side of dR/dt
    return kappa*I
  
I0=0.01 # Initial number(Percentage) of infected people
S0=0.99 # Initial number of potential for infection people 
R0=0 #Initial number of recovered people
beta=3.2 # Number of contact per per unit time*Infection probability
kappa=0.32 #Recovery Rate 1/(total recovery time)
  
nSteps=1000
#Initialization
S=[]# Number of suspectible people
I=[] # Number of infected people
R=[] # Number of recovered people
dt=0.01
#Time integration using 4th order RK method
S.append(S0)
I.append(I0)
R.append(R0)
for i in range (nSteps):
    S0,I0,R0=rk4_3d(S0,I0,R0,dt)
    S.append(S0)
    I.append(I0)
    R.append(R0)
    
#Plotting the results
t=range(0,nSteps+1)
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(t, S, label='Susceptible')
ax.plot(t, I, label='Infected')
ax.plot(t, R, label='Recovered')
plt.title('COVID-19  model-Copyright @ Ramjee Sharma')
ax.legend()
plt.show()
#These codes are written by Ramjee Sharma. No parts of these codes will be
#reproduced for commercial or research use without the author's permission.
#All right reserved.
print("Maximum Infection = ",100*np.round(np.max(I),4),"%")