import numpy as np

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
        self.thickness = 1.0
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
        self.u_inc = np.zeros((self.n_nodes, 2), dtype=float) # displacement increment
        self.a_tilda = np.zeros((self.n_nodes, 2), dtype=float) # unconstrained accel
        # Kinetic
        self.f_int = np.zeros((self.n_nodes, 2), dtype=float)
        self.mass = np.zeros((self.n_nodes, 2), dtype=float)
        self.f_ext = np.zeros((self.n_nodes, 2), dtype=float)
        # Previous
        self.a_prev = np.zeros((self.n_nodes, 2), dtype=float)
        self.v_prev = np.zeros((self.n_nodes, 2), dtype=float)
        self.u_prev = np.zeros((self.n_nodes, 2), dtype=float)
        self.f_int_prev = np.zeros((self.n_nodes, 2), dtype=float)
        self.f_ext_prev = np.zeros((self.n_nodes, 2), dtype=float)

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
        self.waves = np.sqrt(self.E * (1.0 - self.poisson)
             / (self.rho * (1.0 + self.poisson) * (1.0 - 2.0*self.poisson)))
        # Gauss Point Location(s) and Quantities
        self.gp = np.array([0.0,0.0]) #single gauss point
        self.J1 = np.zeros(self.n_elem) # dx/dxi
        self.J2 = np.zeros(self.n_elem) # dy/dxi
        self.J3 = np.zeros(self.n_elem) # dx/deta
        self.J4 = np.zeros(self.n_elem) # dy/deta
        self.detJ = np.zeros(self.n_elem)

        # Boundary Conditions
        self.v_bc = self.input.v_bc
        self.a_bc = self.input.a_bc

        ## Time
        self.n = 0
        self.t = 0.0
        self.tfinal = self.input.tfinal
        self.Co = self.input.Co
        self.dt = self.Co * (np.min(self.clen) / self.waves)
        self.t_prev = 0.0
        self.dt_prev = 0.0

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

            # Calculate volumetric strain rate
            self.dvol[e] = self.dxx[e] + self.dyy[e]

    def matstatupd(self):
        # Plane-strain rate-form
        use_strain_rates = True
        dt = self.dt
        E = self.E
        nu = self.poisson
        mu = E / (2.0 * (1.0 + nu))                         # shear modulus
        lam = (nu * E) / ((1.0 + nu) * (1.0 - 2.0 * nu))    # Lame's first parameter
        if use_strain_rates:
            lam_dt = lam * dt
            mu_dt = mu * dt
            lp2mu_dt = (lam + 2.0 * mu) * dt
        else:
            lam_dt = lam
            mu_dt = mu
            lp2mu_dt = (lam + 2.0 * mu)

        for e in range(self.n_elem):
            wzz_inc = self.wzz[e] * dt if use_strain_rates else self.wzz[e]
            r1 = 2.0 * self.res_sxy[e] * wzz_inc

            # rotated stresses (before elastic increment)
            sxx_rot = self.res_sxx[e] - r1
            syy_rot = self.res_syy[e] + r1
            szz_rot = self.res_szz[e]
            sxy_rot = self.res_sxy[e] + wzz_inc * (self.res_sxx[e] - self.res_syy[e])

            de_xx = self.dxx[e]
            de_yy = self.dyy[e]
            de_xy = self.dxy[e]   # tensor shear rate: eps_dot_xy

            ds_xx = lp2mu_dt * de_xx + lam_dt * de_yy
            ds_yy = lam_dt * de_xx + lp2mu_dt * de_yy
            ds_zz = lam_dt * (de_xx + de_yy)
            ds_xy = 2.0 * mu_dt * de_xy

            self.res_sxx[e] = sxx_rot + ds_xx
            self.res_syy[e] = syy_rot + ds_yy
            self.res_szz[e] = szz_rot + ds_zz
            self.res_sxy[e] = sxy_rot + ds_xy

    def assmb_internal(self):
        # Scale stresses by volume (integration weight)
        sxx = np.zeros(self.n_elem)
        syy = np.zeros(self.n_elem)
        sxy = np.zeros(self.n_elem)
        for e in range(self.n_elem):
            weight = 4.0
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
            total_mass = self.vol[e] * self.rho
            nodal_mass = total_mass / 4.0
            nodes = self.conn[e] - 1
            for node in nodes:
                self.mass[node][0] += nodal_mass
                self.mass[node][1] += nodal_mass

    def el_state_upd(self):
        self.el_geom()
        self.el_rate()
        self.matstatupd()
        self.assmb_internal()
        # Assume mass-conserving problem
        if self.t == 0.0:
            self.assmb_mass()

    def assmb_vbcs(self, t):
        if not self.v_bc:
            return
        for index, velocities in zip(self.v_bc.indexes, self.v_bc.velocities):
            for dof, velocity_callable in enumerate(velocities):
                if velocity_callable is None:
                    continue
                val = velocity_callable(t)   # numeric or None
                if val is not None:
                    self.v[index - 1][dof] = val

    def assmb_abcs(self):
        if self.a_bc:
            for index, accelerations in zip(self.a_bc.indexes, self.a_bc.accelerations):
                for dof, acceleration in enumerate(accelerations):
                    if acceleration is not None or 0:
                        self.a[index - 1][dof] = acceleration

    def save_prev(self):
        self.a_prev, self.v_prev, self.u_prev = np.copy(self.a), np.copy(self.v), np.copy(self.u)
        self.f_int_prev, self.f_ext_prev = np.copy(self.f_int), np.copy(self.f_ext_prev)
        self.t_prev, self.dt_prev = np.copy(self.t), np.copy(self.dt_prev)
  
    def solveq(self):
        self.save_prev()
        self.a = (self.f_ext - self.f_int) / self.mass
        self.a_tilda = np.copy(self.a)
        self.assmb_abcs()
        if self.n == 0:
            self.v += 0.5 * self.a * self.dt
        else:
            self.v += self.a * self.dt
        self.assmb_vbcs(self.t + 0.5 * self.dt)
        self.u_inc = self.v * self.dt
        self.u += self.v * self.dt
        self.coo = self.coo + self.u_inc
        self.n += 1
        self.t += self.dt
        self.f_int.fill(0)