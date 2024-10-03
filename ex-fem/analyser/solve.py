import numpy as np
from database.boundaryConditions import  VelBoundaryConditions as vbc
from database.boundaryConditions import  AccelBoundaryConditions as abc
from database.boundaryConditions import  ForceBoundaryConditions as fbc
from database.boundaryConditions import  SupportBoundaryConditions as sbc

class SubdomainSolution:

    def __init__(self, input):
        self.input = input

        ## Material Constants
        self.E = self.input.young
        self.rho = self.input.density
        self.poisson = self.input.poisson

        ## Geometry
        self.n_elem = self.input.num_elems
        self.n_nodes = self.input.num_nodes
        self.coo = self.input.coordinates
        self.conn = self.input.connectivity
        self.clen = np.zeros(self.n_elem)
        for i in range(self.n_elem):
            dx = self.coo[self.conn[i][1]-1][0] - self.coo[self.conn[i][0]-1][0]
            dy = self.coo[self.conn[i][1]-1][1] - self.coo[self.conn[i][0]-1][1]
            self.clen[i] = np.sqrt(dx**2 + dy**2)
        self.thickness = 0.5 
        self.vol = np.zeros(self.n_elem)
        for i in range(self.n_elem):
            x1, y1 = self.coo[self.conn[i][0]-1]
            x3, y3 = self.coo[self.conn[i][2]-1]
            self.vol[i] = abs((x3 - x1) * (y3 - y1) * self.thickness)
        
        ## Nodal Quantities
        # Kinematic
        self.a = np.zeros((self.n_nodes, 2), dtype=float) #[self.n_nodes][dof] 
        self.v = np.zeros((self.n_nodes, 2), dtype=float)
        self.u = np.zeros((self.n_nodes, 2), dtype=float)
        # Kinetic
        self.f_int = np.zeros((self.n_nodes, 2), dtype=float)
        self.mass = np.zeros((self.n_nodes, 2), dtype=float)
        self.f_ext = np.zeros((self.n_nodes, 2), dtype=float)

        ## Element Quantities
        # Velocity Gradients
        self.lxx = np.zeros(self.n_elem)
        self.lxy = np.zeros(self.n_elem)
        self.lyx = np.zeros(self.n_elem) 
        self.lyy= np.zeros(self.n_elem)       
        # Spin Tensor
        self.wxx = np.zeros(self.n_elem) 
        self.wyy = np.zeros(self.n_elem) 
        self.wzz = np.zeros(self.n_elem)
        self.r1 = np.zeros(self.n_elem)
        self.r2 = np.zeros(self.n_elem)
        self.r3 = np.zeros(self.n_elem)
        # Rate of Deformation
        self.dxx = np.zeros(self.n_elem)
        self.dyy = np.zeros(self.n_elem)
        self.dxy = np.zeros(self.n_elem)
        # Volumetric Strain
        self.dvol = np.zeros(self.n_elem)
        # Stress Measures
        self.sxx = np.zeros(self.n_elem)
        self.syy = np.zeros(self.n_elem)
        self.sxy = np.zeros(self.n_elem)
        self.szz = np.zeros(self.n_elem)
        self.res_sxx = np.zeros(self.n_elem)
        self.res_syy = np.zeros(self.n_elem)
        self.res_sxy = np.zeros(self.n_elem)
        self.res_szz = np.zeros(self.n_elem)
        # Wave Speeds
        self.waves = np.sqrt(self.E / self.rho)
        # Gauss Point Location(s) adn Quantities
        self.gp = np.array([0.0,0.0]) #single gauss point
        self.J1 = np.zeros(self.n_elem) # dx/dxi
        self.J2 = np.zeros(self.n_elem) # dy/dxi
        self.J3 = np.zeros(self.n_elem) # dx/deta
        self.J4 = np.zeros(self.n_elem) # dy/deta
        self.detJ = np.zeros(self.n_elem)

        # Boundary Conditions
        self.v_bc = self.input.v_bc
        self.a_bc = self.input.a_bc
        self.f_bc = self.input.f_bc
        self.s_bc = self.input.s_bc

        ## Time
        self.n = 0
        self.t = 0.0
        self.tfinal = self.input.tfinal
        self.Co = self.input.Co
        self.dt = self.Co * (np.min(self.clen) / self.waves)

        ## Energies
        self.kinetic_energy = []
        self.internal_energy = []
        self.tot_energy = []
        self.timestamps = []

        ## Global Shape Function Derivatives
        self.dN1dx = np.zeros((self.n_elem), dtype=float)
        self.dN1dy = np.zeros((self.n_elem), dtype=float)
        self.dN2dx = np.zeros((self.n_elem), dtype=float)
        self.dN2dy = np.zeros((self.n_elem), dtype=float)
        self.dN3dx = np.zeros((self.n_elem), dtype=float)
        self.dN3dy = np.zeros((self.n_elem), dtype=float)
        self.dN4dx = np.zeros((self.n_elem), dtype=float)
        self.dN4dy = np.zeros((self.n_elem), dtype=float)

        ## Shape Functions 
        self.N1 = 1/4 * (1 - self.gp[0]) * (1 - self.gp[1])
        self.N2 = 1/4 * (1 + self.gp[0]) * (1 - self.gp[1])
        self.N3 = 1/4 * (1 + self.gp[0]) * (1 + self.gp[1])
        self.N4 = 1/4 * (1 - self.gp[0]) * (1 + self.gp[1])
        ## Local Shape Function Derivatives
        self.dN1dxi = 1/4 * -(1 - self.gp[1])
        self.dN2dxi = 1/4 * (1 - self.gp[1])
        self.dN3dxi = 1/4 * (1 + self.gp[1])
        self.dN4dxi = 1/4 * -(1 + self.gp[1])
        self.dN1deta = 1/4 * -(1 - self.gp[0])
        self.dN2deta = 1/4 * -(1 + self.gp[0])
        self.dN3deta = 1/4 * (1 + self.gp[0])
        self.dN4deta = 1/4 * (1 - self.gp[0])
        
    def el_geom(self):
        # Calculate Jacobian of each Element
        for e in range(self.n_elem):
            x1 = self.coo[self.conn[e][0]-1][0]
            x2 = self.coo[self.conn[e][1]-1][0]
            x3 = self.coo[self.conn[e][2]-1][0]
            x4 = self.coo[self.conn[e][3]-1][0]
            y1 = self.coo[self.conn[e][0]-1][1]
            y2 = self.coo[self.conn[e][1]-1][1]
            y3 = self.coo[self.conn[e][2]-1][1]
            y4 = self.coo[self.conn[e][3]-1][1]
            self.J1[e] = x1 * self.dN1dxi + x2 * self.dN2dxi + x3 * self.dN3dxi + x4 * self.dN4dxi
            self.J2[e] = y1 * self.dN1dxi + y2 * self.dN2dxi + y3 * self.dN3dxi + y4 * self.dN4dxi
            self.J3[e] = x1 * self.dN1deta + x2 * self.dN2deta + x3 * self.dN3deta + x4 * self.dN4deta
            self.J4[e] = y1 * self.dN1deta + y2 * self.dN2deta + y3 * self.dN3deta + y4 * self.dN4deta

        # Calculate the determinant of the Jacobian
        for e in range(self.n_elem):
            self.detJ[e] = self.J1[e] * self.J4[e] - self.J2[e] * self.J3[e]

        # Calculate the inverse of the Jacobian
        invJ1 = np.zeros(self.n_elem)
        invJ2 = np.zeros(self.n_elem)
        invJ3 = np.zeros(self.n_elem)
        invJ4 = np.zeros(self.n_elem)
        for e in range(self.n_elem):
            invJ1[e] = (1 / self.detJ[e]) * self.J4[e]
            invJ2[e] = (1 / self.detJ[e]) * -self.J2[e]
            invJ3[e] = (1 / self.detJ[e]) * -self.J3[e]
            invJ4[e] = (1 / self.detJ[e]) * self.J1[e]

        # Calculate the derivatives of the shape functions with respect to the Global Coordinates
        for e in range(self.n_elem):
            self.dN1dx[e] = invJ1[e] * self.dN1dxi + invJ2[e] * self.dN1deta
            self.dN1dy[e] = invJ3[e] * self.dN1dxi + invJ4[e] * self.dN1deta
            self.dN2dx[e] = invJ1[e] * self.dN2dxi + invJ2[e] * self.dN2deta
            self.dN2dy[e] = invJ3[e] * self.dN2dxi + invJ4[e] * self.dN2deta
            self.dN3dx[e] = invJ1[e] * self.dN3dxi + invJ2[e] * self.dN3deta
            self.dN3dy[e] = invJ3[e] * self.dN3dxi + invJ4[e] * self.dN3deta
            self.dN4dx[e] = invJ1[e] * self.dN4dxi + invJ2[e] * self.dN4deta
            self.dN4dy[e] = invJ3[e] * self.dN4dxi + invJ4[e] * self.dN4deta

    def el_rate(self):
        for e in range(self.n_elem):
            node1 = self.conn[e][0] - 1 #node1, node2, node3, node4 = [i - 1 for i in self.conn[e]] 
            node2 = self.conn[e][1] - 1
            node3 = self.conn[e][2] - 1
            node4 = self.conn[e][3] - 1

            # Calculate velocity gradients (L tensor)
            self.lxx[e] = self.v[node1][0] * self.dN1dx[e] + self.v[node2][0] * self.dN2dx[e] + self.v[node3][0] * self.dN3dx[e] + self.v[node4][0] * self.dN4dx[e]  # ∂u/∂x
            self.lxy[e] = self.v[node1][0] * self.dN1dy[e] + self.v[node2][0] * self.dN2dy[e] + self.v[node3][0] * self.dN3dy[e] + self.v[node4][0] * self.dN4dy[e]  # ∂u/∂y
            self.lyx[e] = self.v[node1][1] * self.dN1dx[e] + self.v[node2][1] * self.dN2dx[e] + self.v[node3][1] * self.dN3dx[e] + self.v[node4][1] * self.dN4dx[e]  # ∂v/∂x
            self.lyy[e] = self.v[node1][1] * self.dN1dy[e] + self.v[node2][1] * self.dN2dy[e] + self.v[node3][1] * self.dN3dy[e] + self.v[node4][1] * self.dN4dy[e]  # ∂v/∂y

            # Calculate rate of deformation tensor (D tensor)
            self.dxx[e] = self.lxx[e]   # ε̇_xx
            self.dyy[e] = self.lyy[e]   # ε̇_yy
            self.dxy[e] = 0.5 * (self.lxy[e] + self.lyx[e])  # ε̇_xy

            # Calculate spin tensor (W tensor)
            self.wxx[e] = 0.0  # ω_xx
            self.wyy[e] = 0.0  # ω_yy
            self.wzz[e] = 0.5 * (self.lyx[e] - self.lxy[e])  

            # Calculate volumetric strain rate (optional, but useful)
            self.dvol[e] = self.dxx[e] + self.dyy[e]

    def matstatupd(self):
        for e in range(self.n_elem): # Calculate Spin Increment
            self.wxx[e] = self.wxx[e] * self.dt
            self.wyy[e] = self.wyy[e] * self.dt
            self.wzz[e] = self.wzz[e] * self.dt

        for e in range(self.n_elem): # Rotate Stress (Sn+1 = Sn + Sn*W - W*Sn)
            self.r1[e] = 2.0 * self.res_sxy[e] * self.wzz[e]
            self.r2[e] = 0.0 # 2.0 * res_sxz * wyy
            self.r3[e] = 0.0 # 2.0 * res_syz * wxx

        for e in range(self.n_elem):
            # Rotate stress (Sn+1 = Sn + Sn*W - W*Sn)
            self.sxx[e] = self.res_sxx[e] - self.r1[e] + self.r2[e]
            self.syy[e] = self.res_syy[e] + self.r1[e] - self.r3[e]
            self.szz[e] = self.res_szz[e] - self.r2[e] + self.r3[e]
            self.sxy[e] = self.res_sxy[e] + self.wzz[e] * (self.res_sxx[e] - self.res_syy[e]) # + (wyy * res_syz) + (wxx * res_sxz)
            
        Gdt = self.E / (2.0 * (1.0 + self.poisson)) * self.dt
        C1dt = self.E * (1.0 - self.poisson) / ((1.0 + self.poisson) * (1.0 - 2.0 * self.poisson)) * self.dt
        C2dt = C1dt * self.poisson / (1.0 - self.poisson)

        for e in range(self.n_elem): # Update Stress (Sn+1 = Sn+1 + C * dε)
            self.res_sxx[e] = self.sxx[e] + (C1dt * self.dxx[e] + C2dt * self.dyy[e])
            self.res_syy[e] = self.syy[e] + (C1dt * self.dyy[e] + C2dt * self.dxx[e])
            self.res_szz[e] = self.szz[e] + (C2dt * (self.dxx[e] + self.dyy[e]))
            self.res_sxy[e] = self.sxy[e] + (Gdt * self.dxy[e])

    def assmb_internal(self):
        # Scale stresses by volume (integration weight)
        sxx = np.zeros(self.n_elem)
        syy = np.zeros(self.n_elem)
        sxy = np.zeros(self.n_elem)
        for e in range(self.n_elem):
            weight = 2.0  # Adjust this weight according to your integration scheme
            sxx[e] = self.res_sxx[e] * (weight * self.detJ[e])
            syy[e] = self.res_syy[e] * (weight * self.detJ[e])
            sxy[e] = self.res_sxy[e] * (weight * self.detJ[e])

        # Compute internal forces due to stresses
        for e in range(self.n_elem):
            # Get the indices of the nodes for the current element
            node1, node2, node3, node4 = [self.conn[e][i] - 1 for i in range(4)]

            # Contribution to internal force from element stress and shape function derivatives
            self.f_int[node1][0] += sxx[e] * self.dN1dx[e] + sxy[e] * self.dN1dy[e]
            self.f_int[node1][1] += sxy[e] * self.dN1dx[e] + syy[e] * self.dN1dy[e]

            self.f_int[node2][0] += sxx[e] * self.dN2dx[e] + sxy[e] * self.dN2dy[e]
            self.f_int[node2][1] += sxy[e] * self.dN2dx[e] + syy[e] * self.dN2dy[e]

            self.f_int[node3][0] += sxx[e] * self.dN3dx[e] + sxy[e] * self.dN3dy[e]
            self.f_int[node3][1] += sxy[e] * self.dN3dx[e] + syy[e] * self.dN3dy[e]

            self.f_int[node4][0] += sxx[e] * self.dN4dx[e] + sxy[e] * self.dN4dy[e]
            self.f_int[node4][1] += sxy[e] * self.dN4dx[e] + syy[e] * self.dN4dy[e]

    def assmb_mass(self):
        for e in range(self.n_elem):
            element_mass = self.vol[e] * self.rho / 4.0  # Divide by 4 for 4-node elements
            nodes = self.conn[e] - 1
            for node in nodes:
                self.mass[node][0] += element_mass
                self.mass[node][1] += element_mass

    def el_state_upd(self):
        self.el_geom()
        self.el_rate()
        self.matstatupd()
        self.assmb_internal()
        # Assume mass-conserving problem
        if self.t == 0.0:
            self.assmb_mass()

    def assmb_vbcs(self, t):
        if self.v_bc:
            for index, velocities in zip(self.v_bc.indexes, self.v_bc.velocities):
                for dof, velocity in enumerate(velocities):
                    if velocity is not None or 0:
                        self.v[index - 1][dof] = velocity(t)

    def assmb_abcs(self):
        if self.a_bc:
            for index, accelerations in zip(self.a_bc.indexes, self.a_bc.accelerations):
                for dof, acceleration in enumerate(accelerations):
                    if acceleration is not None or 0:
                        self.a[index - 1][dof] = acceleration
  
    def solveq(self):
        self.a = (self.f_ext - self.f_int) / self.mass
        self.assmb_abcs()
        if self.n == 0:
            self.v += 0.5 * self.a * self.dt
        else:
            self.v += self.a * self.dt
        self.assmb_vbcs(self.t + 0.5 * self.dt)
        self.u += self.v * self.dt
        self.n += 1
        self.t += self.dt
        self.f_int.fill(0)
        

