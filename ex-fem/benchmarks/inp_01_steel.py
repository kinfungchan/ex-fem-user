import numpy as np
from .baseInput import BaseInput
from .bar_mesh import generate_mesh
from database.boundaryConditions import  VelBoundaryConditions as vbc
from database.boundaryConditions import  AccelBoundaryConditions as abc

##Start of Input File 
class SteelBarInput(BaseInput):
    '''
    Input Constructor 
        Geometry: 1x1 Element Cross Section from 0-250 X-Direction
        Material: Steel
        Boundary Conditions: Half Sine Velocity & Spatial Coupled Accelerations
        Time Integration: Co=0.5 (Coupling with 02 for m=2) 
    '''
    def __init__(self):
        super().__init__()

        # Material properties
        self.young = 207 # GPa
        self.density = 7.83e-6 # kg/mm^3
        self.poisson = 0.0 # Poisson's Ratio
        self.num_elems = 250 # Number of elements
        self.num_nodes = 502 # Number of nodes

        # Mesh generation
        self.coordinates, self.connectivity = generate_mesh(self.num_elems, 0.0, 1.0, 1) # 1 Row Bar of 250 elements, 1.0 each length, 1.0 width
        
        # Input Boundary Conditions
        def vel(t): return vbc.velbc(t, self.max_length, self.young, self.density) 
        velBCs = vbc(np.array([1, 252]), np.array([[vel, None], [vel, None]])) # Indexes Node and [Velocities, DOF]
        # Boundary Conditions
        self.v_bc = velBCs
        self.a_bc = abc(np.array([]).reshape((0, 0)), np.array([]).reshape((0, 0)))
        self.f_bc = None
        self.s_bc = None

        # Time variables
        self.tfinal = 0.0002 * 500
        self.Co = 0.5

        ## Additional Parameters
        self.max_length = 125.0 # Maximum length of the bar

        # Validate input
        self.validate_input()