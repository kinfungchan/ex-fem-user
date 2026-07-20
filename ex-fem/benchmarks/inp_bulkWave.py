import numpy as np
from .baseInput import BaseInput
from .bar_mesh import generate_bar_mesh
from database.boundaryConditions import  DisplBoundaryConditions as ubc
from database.boundaryConditions import  VelBoundaryConditions as vbc

class Bulk_Wave_Input_Monolithic(BaseInput):
    '''
    Input Constructor 
        Geometry: 1x1 Coarse Square Domain for Monolithic Subdomain of Chiappa Benchmark
        Material: Steel
        Boundary Conditions: Prescribed Velocities about Centre Square
        Time Integration: Co=0.1
    '''
    def __init__(self):
        super().__init__()

        # Material properties
        self.young = 209e09 # GPa
        self.density = 7800 # kg/mm^3
        self.poisson = 1.0/3 # Poisson's Ratio

        self.num_elems = 1600 # Number of elements
        self.num_nodes = 1681 # Number of nodes
        # Mesh generation
        num_elems_x = 40
        num_elems_y = 40
        min_x = 0.0
        min_y = 0.0
        length_x = 1.0/num_elems_x # Element Length in X - 0.025
        length_y = 1.0/num_elems_y # Element Length in Y - 0.025
        self.coordinates, self.connectivity = generate_bar_mesh(num_elems_x, num_elems_y, min_x, min_y, length_x, length_y, direction='Y')

        # Input Boundary Conditions 
        # parameters from the paper (bulk case, Sec.3.1)
        self.F1 = 1   # x velocity inside D [m/s]
        self.F2 = 0.0   # y velocity inside D [m/s]
        h3 = 0.45
        h1 = 0.45
        h4 = 0.1
        h2 = 0.1

        coords = self.coordinates
        x = coords[:,0]
        y = coords[:,1]

        # subdomain D (rectangular): [h3, h3+h4] x [h1, h1+h2]
        inside_D_mask = (x >= h3) & (x <= (h3 + h4)) & (y >= h1) & (y <= (h1 + h2))
        self.inside_indices = np.where(inside_D_mask)[0]

        c = ((self.young*(1-self.poisson))/(self.density*(1+self.poisson)*(1-2*self.poisson)))**0.5
        L = min(length_x, length_y, (length_x**2 + length_y**2)**0.5)
        CFL = 0.1
        dt_crit = CFL * (L / c)
        def vel(t, factor, velic_fn):
            raw = velic_fn(t, dt_crit+1e-12)
            if raw is None:
                return None
            return factor * raw

        vel_callables = []
        for _ in self.inside_indices:
            vx_fn = lambda t, F1=self.F1, velic_fn=vbc.velic: vel(t, F1, velic_fn)
            vy_fn = lambda t, F2=self.F2, velic_fn=vbc.velic: vel(t, F2, velic_fn)
            vel_callables.append([vx_fn, vy_fn])
        velBCs_array = np.array(vel_callables, dtype=object)
        velICs = vbc(self.inside_indices + 1, velBCs_array)
        self.v_bc = velICs

        tol = 1e-8
        left_mask   = np.isclose(coords[:, 0], 0.0, atol=tol)    # x == 0
        right_mask  = np.isclose(coords[:, 0], 1.0, atol=tol)    # x == 1
        lower_mask  = np.isclose(coords[:, 1], 0.0, atol=tol)    # y == 0
        upper_mask  = np.isclose(coords[:, 1], 1.0, atol=tol)    # y == 1

        left_edge   = np.where(left_mask)[0] + 1
        right_edge  = np.where(right_mask)[0] + 1
        lower_edge  = np.where(lower_mask)[0] + 1
        upper_edge  = np.where(upper_mask)[0] + 1

        all_edge_indices_1based = np.unique(
            np.concatenate([left_edge, right_edge, lower_edge, upper_edge])
        )

        def displ(t):
            return ubc.displbcConstant(t, 0.0)

        N_constrained = len(all_edge_indices_1based)
        displ_callables = np.empty((N_constrained, 2), dtype=object)

        for i, node1 in enumerate(all_edge_indices_1based):
            node0 = node1 - 1 

            ux_callable = None
            uy_callable = None

            # if node on left or right edge -> constrain ux
            if left_mask[node0] or right_mask[node0]:
                ux_callable = displ

            # if node on lower or upper edge -> constrain uy
            if lower_mask[node0] or upper_mask[node0]:
                uy_callable = displ

            displ_callables[i, 0] = ux_callable
            displ_callables[i, 1] = uy_callable
        displBCs = ubc(all_edge_indices_1based, displ_callables)

        # Boundary Conditions
        self.u_bc = displBCs
        self.v_bc = velICs
        self.a_bc = None
        self.f_bc = None

        # Time variables
        self.tfinal = 0.00009
        self.Co = CFL

        ## Additional Parameters
        closest_idx = np.argmin(np.linalg.norm(self.coordinates - np.array([0.35,0.35]), axis=1))
        self.probe_idx = closest_idx

        # Validate input
        self.validate_input()