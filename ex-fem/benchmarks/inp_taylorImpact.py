import numpy as np
from .baseInput import BaseInput
from database.boundaryConditions import VelBoundaryConditions as vbc
from database.boundaryConditions import AccelBoundaryConditions as abc

class CopperImpactInput(BaseInput):
    '''
    Input Constructor for Taylor Impact of an OFHC Copper Rod
        Geometry: Dynamically sized quad elements (Half-symmetry model)
        Material: OFHC Copper with Linear Isotropic Hardening
        Boundary Conditions: Rigid frictionless wall at Y=0, Symmetry at X=0
    '''
    def __init__(self, nx=5, ny=50): 
        super().__init__()

        # Material properties (OFHC Copper) - Consistent unit system: mm, ton, s, N, MPa
        self.young = 117.0e3 
        self.poisson = 0.35 
        self.density = 8.93e-9 
        
        # Mesh parameters
        lx = 3.2
        ly = 32.4
        
        self.num_elems = nx * ny
        self.num_nodes = (nx + 1) * (ny + 1)
        
        # Plasticity configurations (Linear hardening: Y = 400 + 100*epv)
        self.is_plastic = True
        self.plastic_curve = np.array([
            [0.0, 400.0],
            [10.0, 400.0 + 100.0 * 10.0] 
        ])

        # Mesh generation
        xs = np.linspace(0, lx, nx + 1)
        ys = np.linspace(0, ly, ny + 1)
        
        coo = []
        for y in ys:
            for x in xs:
                coo.append([x, y])
        self.coordinates = np.array(coo)

        conn = []
        for j in range(ny):
            for i in range(nx):
                n1 = j * (nx + 1) + i + 1
                n2 = n1 + 1
                n3 = n2 + (nx + 1)
                n4 = n1 + (nx + 1)
                conn.append([n1, n2, n3, n4])
        self.connectivity = np.array(conn)
        
        # Boundary Conditions
        bc_indexes = []
        bc_funcs = []
        
        def fix_vel(t): return 0.0
        
        # Apply symmetry on X=0 (left) and rigid wall on Y=0 (bottom)
        for j in range(ny + 1):
            for i in range(nx + 1):
                node_id = j * (nx + 1) + i + 1
                
                if j == 0 and i == 0:
                    # Bottom-Left corner: Fix X and Y
                    bc_indexes.append(node_id)
                    bc_funcs.append([fix_vel, fix_vel])
                elif j == 0:
                    # Bottom edge: Free X, Fix Y
                    bc_indexes.append(node_id)
                    bc_funcs.append([None, fix_vel])
                elif i == 0:
                    # Left edge: Fix X, Free Y
                    bc_indexes.append(node_id)
                    bc_funcs.append([fix_vel, None])
                    
        self.v_bc = vbc(np.array(bc_indexes), np.array(bc_funcs))
        self.u_bc = None
        self.a_bc = abc(np.array([]).reshape((0, 0)), np.array([]).reshape((0, 0)))
        self.f_bc = None

        # Time variables
        self.tfinal = 80.0e-6 
        self.Co = 0.05

        ## Additional Parameters
        self.is_axisymmetric = True

        self.validate_input()