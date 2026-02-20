# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 19:12:34 2026

@author: Zhiyi

This is a library for simulating laser cooling phenomena.

To use the library:
    1) Create objects: the potential, the laser and the ion with relevant parameters
    2) Create object: cooling(takes in the objects in 1) as arguments)
    3) Integration(takes in cooling as argument)

Planned updates:
    1) Making diffusion tensor anistropic
    2) Making the integration tool more generic(can be used for other cooling models)
"""

import numpy as np


# Constants
kB = 1.380649e-23                 # Boltzmann's constant
hbar = 1.054e-34                  # Reduced Planck's constant
atomic_mass_unit = 1.66054e-27    # Atomic mass unit
I = np.eye(3)


class harmonic_potential:
    
    def __init__(self, freq_x, freq_y, freq_z):
        
        self.freq_x = freq_x
        self.freq_y = freq_y
        self.freq_z = freq_z
    
    def generate_potential(self, x_range, y_range, z_range, step):
        
        return potential_data
        
    def rotate_Hessian(self, ):
        
        return potential_data_rotated
        
        
class laser:
    
    def _init_(self, laser_wavelength, laser_linewidth, saturation_parameter):
        """
        laser_wavelength: laser wavelength in nm
        laser_linewidth: laser linewidth in nm
        saturation_parameter: laser saturation parameter, between 0 and 1
        """
        
        self.laser_wavelength = laser_wavelength
        self.laser_linewidth = laser_linewidth
        self.saturation_parameter = saturation_parameter
        
    def wavevector(self, laser_angles):
        """
        laser_angles: should be in the form of [theta, phi], where theta and phi are angles in spherical coordinates
        """
        
        theta = np.deg2rad(laser_angles[0])
        phi = np.deg2rad(laser_angles[1])
        wavenumber = 2 * np.pi / self.laser_wavelength
        
        kvec = wavenumber * np.array([
        np.sin(theta)*np.cos(phi),
        np.sin(theta)*np.sin(phi),
        np.cos(theta)
    ])
        return kvec
    
    def wavenumber(self):
        
        k = 2*np.pi / self.laser_wavelength
        
        return k
 
        
class ion:
    
    def _init_(self, mass, transition_wavelength, natural_linewidth, initial_position, initial_temperature):
        self.mass = mass
        self.transition_wavelength = transition_wavelength
        self.natural_linewidth = natural_linewidth
        self.initial_position = initial_position
        self.initial_temperature = initial_temperature
        
    def mass_kg(self):
        mass_kg = self.mass * atomic_mass_unit
        
        return mass_kg
    
    # Initial condition
    def initial_velocity(self):
        sigma_v = np.sqrt(kB * self.initial_temperature/self.mass_kg())
        v = [sigma_v,sigma_v,sigma_v]
        
        return v



# Calculating necessary parameters and forces in Doppler cooling
class Doppler_forces:
    def __init__(self, potential, laser, ion):
        self. laser = laser
        self. ion = ion
        
        
    # Calculate the scattering rate
    def scattering_rate(self):
        
        # Laser detuning
        delta = self.ion.transition_wavelength - self.laser.wavelength
        
        # Natural linewidth of the atomic transition
        gamma = np.deg2rad(self.ion.natural_linewidth)           
        
        # Detuning that leads to the Doppler limit
        delta_Doppler = -gamma/2            
        
        # Saturation parameter                   
        s0 = laser.saturation_parameter                        
        
        if delta != delta_Doppler:
            print('WARNING: The detuning will not lead to Doppler limit')
        
        return (gamma/2) * s0 / (1 + s0 + (2 * delta/gamma)**2)
    
    # Calculate the Doppler limit
    def cooling_limit(self):
        
        # Natural linewidth of the atomic transition
        gamma = np.deg2rad(self.ion.natural_linewidth)
        
        # Calculate the Doppler limit
        TD = hbar * gamma / (2 * kB)
        
        return TD
    
    # Friction coefficient if approximate the absorption recoil as friction
    def friction_coefficient(self):
        
        # Laser detuning
        delta = self.ion.transition_wavelength - self.laser.wavelength
        
        # Natural linewidth of the atomic transition
        gamma = np.deg2rad(self.ion.natural_linewidth)
        
        # Saturation parameter                   
        s0 = self.laser.saturation_parameter       
        
        alpha = -(4 * s0 * delta / gamma**2) / (1 + s0 + (2*delta/gamma)**2)**2
        return alpha
    
    def diffusion_tensor(self):
        D = (1/6) * (hbar*self.laser.wavenumber)**2 * self.scattering_rate() * I
        
        return D


    
def integration(Doppler_forces, dt, Nsteps):
    position_history = []
    
    for i in range(Nsteps):
        
    
        
    



        
        

