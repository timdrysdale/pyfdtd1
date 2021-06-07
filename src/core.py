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
    GAUSSIAN = 0
    
class Field(Enum):
    ELECTRIC = 0
    MAGNETIC = 0
    
class Boundary(Enum):
    MUR=0
    BARE=0

class SourcePositionNotValid(Exception):
    """
    Exception raised when the source is specified in a position that is not allowed,
    e.g. at a boundary or outside the domain
    
    Attributes:
        source_position -- input source position which caused the error
        domain_size -- size of the domain 
    """    
    def __init__(self, source_position, domain_size):
        self.source_position = source_position
        self.domain_size = domain_size
        self.message = "source position %d is outside the valid positions of [1:(N-1)] i.e. [1:%d]"%(source_position, domain_size-1)
        super().__init__(self.message)
        
    def __str__(self):
        return f'{self.message}'
    
class InvalidSourceWave(Exception):
    """
    Exception raised when the source  wave type is not known
   
    """    
   
    def __init__(self):
        self.message = "Unknown wave type: Source wave should be Wave.GAUSSIAN or Wave.SINE"
        super().__init__(self.message)
        
    def __str__(self):
        return f'{self.message}'    
    

class FDTD1:
    
    def __init__(self,
                 dx,
                 N,
                 source_position,
                 courant = 0.5, 
                 source_field = Field.ELECTRIC,
                 source_type = Source.SOFT,
                 source_wave = Wave.GAUSSIAN,
                 boundary_type = Boundary.MUR,
                 Z = (1./(epsilon_0 * speed_of_light))):
        self.boundary_type = boundary_type
        self.courant = courant
        self.init_fields(N)
        self.dx = dx
        self.dt = self.get_dt(self.courant, self.dx)
        self.time_step = 0 #current time step
        self.Z = Z
        
        self.source_field = source_field
        
        if (source_position < 1) or (source_position >= N):
            raise SourcePositionNotValid(source_position, N)
            raise ValueError("source position %d is outside the valid positions of [1:(N-1)] i.e. [1:%d]"%(source_position,N-1))
        
        self.source_position = source_position
        if not (source_type == Wave.GAUSSIAN or source_type == Wave.SINE):
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
            
    def init_sine(self, omega = 20, magnitude = 1):
        self.sine_omega = omega
        self.sine_magnitude = magnitude
        self.sine_value = 0
    
    def init_source(self):
        #Provide usable default sources, for convenience
        
        if self.source_wave == Wave.GAUSSIAN:
            self.init_gaussian()
            
        if self.source_wave == Wave.SINE:    
            self.init_sine()
        
    def update_gaussian(self):
       """
       Using the current time step, calculate the value of the
       Gaussian source, and store it for a source update routine to use
       """
       arg = ((self.time_step-self.gaussian_delay)/self.gaussian_width)**2
       self.gaussian_value = exp(arg)
   
    def update_sine(self):
        """
        We assume that since 
        t = m *dt 
        then the source is
        sin(omega * m * dt) 
        """
        sine = sin(self.omega * self.dt * self.time_step)
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
            
        if self.source_field == Field.Electric:
            
            if self.source_type == Source.HARD:
                
                self.Ez[self.source_position] = self.source_value
            
            if self.source_type == Source.SOFT:
                
                self.Ez[self.source_position] = self.Ez[self.source_position] + self.source_value
            
                       
        if self.source_field== Field.Magnetic:
            
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
        self.Hy = np.array(np.zeros(N), dtype=float) 
           
    def init_update_fields(self):
        """
        See p3 in
        https://my.ece.utah.edu/~ece6340/LECTURES/lecture%2014/FDTD.pdf
        For Taflove's normalisation of the fields
        """
        self.field_normalisation = 1/(mu_0 * epsilon_0)**0.5 * self.dt / self.dx
        
    def update_fields(self):
        h_gradient = self.Hy[0:-2] - self.Hy[1:-1]
        self.Ez[1:] = self.Ez[1:]  + self.field_normalisation * h_gradient
        e_gradient = self.Ez[0:-2] - self.Ez[1:-1]
        self.Hy[0:-1] = self.Hy[0:-1] + self.field_normalisation * e_gradient
        
    
    def update_boundaries(self):
        if self.boundary_type == Boundary.BARE:
            return #nothing to do
        if self.boundary_type == Boundary.MUR:
            print("boundary not implemented yet")
            pass
            
      
   
    def get_dt(self, courant, dx):
        return courant * dx / speed_of_light
    
    