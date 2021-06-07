#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FDTD1 is a class for one-dimensional finite difference time domain
calculations using 
 - boundaries: ABC or reflecting 
 - waveforms: continuous wave or Gaussian
 - sources: soft or hard
"""
from enum import Enum
from math import exp, sin
import numpy as np
from scipy.constants import speed_of_light, epsilon_0, mu_0

class Source(Enum):
    HARD = 0
    SOFT = 1

class Wave(Enum):
    SINE = 0
    GAUSSIAN = 1
    
class Field(Enum):
    ELECTRIC = 0
    MAGNETIC = 1
    
class Boundary(Enum):
    BARE=0
    MUR=1

class FDTD1:
    
    def __init__(self,
                 dx,
                 N,
                 source_position,
                 Courant_factor = 0.5, 
                 source_field = Field.ELECTRIC,
                 source_type = Source.SOFT,
                 source_wave = Wave.GAUSSIAN,
                 boundary_type = Boundary.MUR,
                 Z = (1./(epsilon_0 * speed_of_light))):
        
        if boundary_type == Boundary.MUR and Courant_factor != 0.5:
            raise ValueError("For Mur boundary, courant must be set to default 0.5, not %f"%Courant_factor)
        self.boundary_type = boundary_type
        self.courant = Courant_factor
        self.init_boundaries()
        
        self.N = N
        self.init_fields(N)
        self.dx = dx
        self.dt = self.get_dt(self.courant, self.dx)
        self.time_step = 0 #current time step
        self.Z = Z
        
        self.source_field = source_field
        
        if (source_position < 1) or (source_position >= N):
            raise ValueError("source position %d is outside the valid positions of [1:(N-1)] i.e. [1:%d]"%(source_position,N-1))
        
        self.source_position = source_position
        if not (source_wave == Wave.GAUSSIAN or source_wave == Wave.SINE):
            raise TypeError("Unknown wave type: Source wave should be Wave.GAUSSIAN or Wave.SINE")
        
        self.source_type = source_type
        self.source_wave = source_wave
        self.source_value = 0
        self.init_source()
        self.init_update_fields()
        
                
    def hard_electric(self, n, E):
        self.Ez[n] = E
        
    def soft_electric(self, n, E):
        self.Ez[n] = self.Ez[n] + E

    def init_gaussian(self, temporal_delay=30, pulse_width=10):
        """ 
        See p116 https://eecs.wsu.edu/~schneidj/ufdtd/chap5.pdf 
        We can neglect the dt in the gaussian because it is present
        in the top line and bottom line, so cancels.
        """  
        self.gaussian_delay = temporal_delay
        self.gaussian_width = pulse_width
        self.gaussian_value = 0
            
    def init_sine(self, omega, magnitude = 1):
        self.sine_omega = omega
        self.sine_magnitude = magnitude
        self.sine_value = 0
    
    def init_source(self):
        #Provide usable default sources, for convenience
        
        if self.source_wave == Wave.GAUSSIAN:
            self.init_gaussian()
            
        if self.source_wave == Wave.SINE:    
            self.init_sine(0.3/self.dt)
        
    def update_gaussian(self):
       """
       Using the current time step, calculate the value of the
       Gaussian source, and store it for a source update routine to use
       """
       arg = ((self.time_step-self.gaussian_delay)/self.gaussian_width)**2
       self.gaussian_value = exp(-arg)
   
    def update_sine(self):
        """
        We assume that since 
        t = m *dt 
        then the source is
        sin(omega * m * dt) 
        """
        sine = sin(self.sine_omega * self.dt * self.time_step)
        self.sine_value = self.sine_magnitude * sine
        
     
    def update_source(self):    
        """
        Update self.source_value with current source waveform value
        """
        if self.source_wave == Wave.GAUSSIAN:
            self.update_gaussian()
            self.source_value = self.gaussian_value
              
        if self.source_wave == Wave.SINE:    
            self.update_sine() 
            self.source_value = self.sine_value
            
        if self.source_field == Field.ELECTRIC:
            
            if self.source_type == Source.HARD:
                
                self.Ez[self.source_position] = self.source_value
            
            if self.source_type == Source.SOFT:
                
                self.Ez[self.source_position] = self.Ez[self.source_position] + self.source_value
            
                       
        if self.source_field== Field.MAGNETIC:
            
            if self.source_type == Source.HARD:
                
                self.Hy[self.source_position] = self.source_value
            
            if self.source_type == Source.SOFT:
                
                self.Hy[self.source_position] = self.Hy[self.source_position] + self.source_value
    
    def init_fields(self, N):
        """
        Ez[0] is to the left of Hy[0]
        i.e. Ez[0]       Ez[1]        Ez[2]
                   Hy[0]       Hy[1]        ...
        """
        self.Ez = np.array(np.zeros(N), dtype=float)
        self.Hy = np.array(np.zeros(N-1), dtype=float) 
           
    def init_update_fields(self):
        """
        See p3 in
        https://my.ece.utah.edu/~ece6340/LECTURES/lecture%2014/FDTD.pdf
        For Taflove's normalisation of the fields
        """
        self.field_normalisation = 1/(mu_0 * epsilon_0)**0.5 * self.dt / self.dx
        
    def update_fields(self):
        
        cc = self.field_normalisation
        
        for n in range(1,self.N-1):
            self.Ez[n] = self.Ez[n] + cc * (self.Hy[n-1]-self.Hy[n])
            
        self.update_source()   
        
        for n in range(0,self.N-1):
            self.Hy[n] = self.Hy[n] + cc * (self.Ez[n]-self.Ez[n+1])       

        
    def init_boundaries(self):
        if self.boundary_type == Boundary.BARE:
            return #nothing to do
        if self.boundary_type == Boundary.MUR:
            self.MurE0previous = 0
            self.MurENprevious = 0
            
    def update_boundaries(self):
        if self.boundary_type == Boundary.BARE:
            return #nothing to do
        if self.boundary_type == Boundary.MUR:
            self.Ez[0] = self.MurE0previous
            self.MurE0previous = self.Ez[1]
            self.Ez[-1] = self.MurENprevious
            self.MurENprevious = self.Ez[-2]
   
    def get_dt(self, courant, dx):
        return courant * dx / speed_of_light
    
    
    def iterate(self):
        self.time_step = self.time_step + 1
        self.update_fields()
        self.update_boundaries()
        
    
if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
     
    demo = FDTD1(0.1,50,25, source_wave = Wave.GAUSSIAN, source_type = Source.HARD, boundary_type = Boundary.MUR) 
    source = []
    for n in range(100):
        demo.time_step = demo.time_step + 1
        demo.update_source()
        source.append(demo.source_value)
    plt.figure()    
    plt.plot(source)    
    plt.xlabel('time step')     
    plt.ylabel('E-field amplitude (V)')
    plt.title('Gaussian source')
    plt.savefig('../img/gaussian_source.png', dpi=100)
        
    demo = FDTD1(0.1,50,25, source_wave = Wave.SINE, source_type = Source.HARD, boundary_type = Boundary.MUR) 
    source = []
    for n in range(100):
        demo.time_step = demo.time_step + 1
        demo.update_source()
        source.append(demo.source_value)
    plt.figure()    
    plt.plot(source)    
    plt.xlabel('time step')     
    plt.ylabel('E-field amplitude (V)')
    plt.title('Sine wave source')
    plt.savefig('../img/sine_source.png', dpi=100)
        
    demo = FDTD1(0.1,100,50, source_wave = Wave.GAUSSIAN, source_type = Source.SOFT, boundary_type = Boundary.BARE) 
    source = []
    plt.figure()
    offset = 0
    for n in range(400):
        demo.iterate()
        if n%10==0:
            plt.plot(demo.Ez + offset)
            offset = offset + 1
   
    plt.xlabel('Position (1/dx)')     
    plt.ylabel('E-field amplitude (V)')
    plt.title('Soft Gaussian source in bare bounded domain')
    fig = plt.gcf()
    fig.set_size_inches(5, 10)
    fig.savefig('../img/bare_gaussian_soft.png', dpi=100)
    
    demo = FDTD1(0.1,100,50, source_wave = Wave.GAUSSIAN, source_type = Source.HARD, boundary_type = Boundary.BARE) 
    source = []
    plt.figure()
    offset = 0
    for n in range(400):
        demo.iterate()
        if n%10==0:
            plt.plot(demo.Ez + offset)
            offset = offset + 1
   
    plt.xlabel('Position (1/dx)')     
    plt.ylabel('E-field amplitude (V)')
    plt.title('Hard Gaussian source in bare bounded domain')
    fig = plt.gcf()
    fig.set_size_inches(5, 10)
    fig.savefig('../img/bare_gaussian_hard.png', dpi=100)   

    demo = FDTD1(0.1,100,50, source_wave = Wave.GAUSSIAN, source_type = Source.SOFT, boundary_type = Boundary.MUR) 
    source = []
    plt.figure()
    offset = 0
    for n in range(400):
        demo.iterate()
        if n%10==0:
            plt.plot(demo.Ez + offset)
            offset = offset + 1
   
    plt.xlabel('Position (1/dx)')     
    plt.ylabel('E-field amplitude (V)')
    plt.title('Soft Gaussian source in Mur bounded domain')   
    fig = plt.gcf()
    fig.set_size_inches(5, 10)
    fig.savefig('../img/mur_gaussian_soft.png', dpi=100)