# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 18:28:40 2026

@author: Zhiyi
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d


# Constants
kB = 1.380649e-23       # kB
m = 40 * 1.66054e-27    # Ion(40Ca+) mass
hbar = 1.054e-34

# Trap frequencies
omega = np.array([2*np.pi*1e6, 2*np.pi*0.3e6, 2*np.pi*0.5e6])

# Laser parameters
Gamma = 2*np.pi*22e6     # natural linewidth for 397nm transition
k = 2*np.pi / 397e-9     # wave number
s0 = 0.3                 # saturation parameter
Delta = -Gamma / 2       # detuning that leads to doppler limit
T_D = hbar * Gamma / (2 * kB) # Doppler limit
Rotation_Steps = 10

# List for temperature history lists
T_history = []  # holds T_history_1, T_History_2...
Average_Temp = []

for i in range(Rotation_Steps):
    T_history_i = []  # temperature history for rotation angle (90/Rotation_Steps)*i
    T_history.append(T_history_i)

for i in range(Rotation_Steps):
    # Laser wavevector
    theta = np.deg2rad(45) # laser angle fixed at 45 degrees in x-z plane
    Rotation_Angle = i*(90/Rotation_Steps)
    phi = np.deg2rad(Rotation_Angle)  # angle with x-axis in x-y plane
    
    kvec = k * np.array([
        np.sin(theta)*np.cos(phi),
        np.sin(theta)*np.sin(phi),
        np.cos(theta)
    ])
    
    
    def scattering_rate(delta):
        return (Gamma/2) * s0 / (1 + s0 + (2*delta/Gamma)**2)
    
    R0 = scattering_rate(Delta)
    
    # dR/dDelta
    dR = -(4 * s0 * Delta / Gamma**2) \
         / (1 + s0 + (2*Delta/Gamma)**2)**2
    
    alpha = -np.outer(kvec, kvec) * dR * hbar
    
    # Calculating the diffusion tensor
    
    k_mag = np.linalg.norm(kvec)
    I = np.eye(3)

    # D = 0.5 * hbar**2 * R0 * (np.outer(kvec, kvec) + (k_mag**1 / 3) * I)
    
    # Time step
    dt = 1e-9
    Nsteps = 4000000
    
    #Initial conditions
    T0 = 0.1
    sigma_v = np.sqrt(kB*T0/m)
    
    r = np.zeros(3)
    v = [sigma_v,sigma_v,sigma_v]
    #v = np.random.normal(0, sigma_v, 3)
    
    # Numerical integration
    T_history[i] = []
    
    for j in range(Nsteps):
        # Trap force
        Ftrap = -m * omega**2 * r
        
        # Testing if anharmonicity removes dark modes
        # Ftrap = -m*omega**2*r - lambda*r**3
    
        # Calculating velocity dependent scattering rate
        delta_eff = Delta - np.dot(kvec, v)
        R = scattering_rate(delta_eff)
        
        # Friction
        Fcool = hbar * kvec * R
    
        # Noise
        # D = 0.5 * hbar**2 * R * (np.outer(kvec, kvec) + (k_mag**1 / 3) * I)
        D = (1/6) * (hbar*k)**2 * R * I
        Sigma = 2 * D * dt / m**2
        eta = np.random.multivariate_normal(np.zeros(3), Sigma)
        
        # Velocity update
        v += (Ftrap + Fcool)/m * dt + eta
    
        # Position update
        r += v * dt
    
        if j % 100 == 0:
            T_inst = m*np.dot(v, v)/(3*kB)
            T_history[i].append(T_inst)
    '''
    # Averaging out secular motion
    def T_average(temp, N=100):
        temp = np.array(temp)
        kernel = np.ones(N) / N
        return np.convolve(temp, kernel, mode='same')
    '''
    def T_average(temp, N=100):
        return uniform_filter1d(temp, size=N, mode='reflect')
    
    Average_Temp.append(T_average(T_history[i]))

    
t = np.arange(len(T_history[0])) * dt * 100
plt.axhline(T_D, linestyle='--', color='k', label='Doppler limit')
plt.semilogy(t*1e3, T_history[0])
plt.xlabel("Time [ms]")
plt.ylabel("Temperature [K]")
plt.show()


#%%
for i in range(Rotation_Steps):
    plt.semilogy(t*1e3, Average_Temp[i], linewidth=0.5, label=f'Rotation Angle = {10*i}°')
    
plt.axhline(T_D, linestyle='--', color='k', label='Doppler limit')
plt.xlabel("Time [ms]")
plt.ylabel("Average Temperature in 10^-4s [K]")
#plt.legend()
#plt.xlim(0.01,0.5)
#plt.ylim(0.0005,0.1)
plt.show()


#%%
def average_temperature(time, temperature, t_start, t_end):
    mask = (time >= t_start) & (time <= t_end)
    if not np.any(mask):
        raise ValueError("No data points found in the given time interval.")

    return temperature[mask].mean()

Final_Temperature_list = []
for i in range(Rotation_Steps):
    Final_Temperature = average_temperature(t*1e3, Average_Temp[i], 2.0, 3.5)
    print(f'Final Temperature for Rotation Angle {10*i}°')
    print(Final_Temperature)
    Final_Temperature_list.append(Final_Temperature)

print(T_D)

#%%
Rotation_Angle = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
plt.scatter(Rotation_Angle, Final_Temperature_list)
plt.axhline(T_D, linestyle='--', color='k', label='Doppler limit')
plt.xlabel("Rotation Angle [°]")
plt.ylabel("Average Temperature from 2.0 to 3.5 ms [K]")
plt.show()
