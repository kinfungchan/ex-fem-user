import numpy as np

class SubdomainSolution:

    def __init__(self, input):
        self.input = input

        ## Material Constants
        self.E = self.input.young
        self.rho = self.input.density
        self.poisson = self.input.poisson
        self.C_q = getattr(self.input, 'C_q', 1.5)
        self.C_l = getattr(self.input, 'C_l', 0.06)
        # Formulation Flags
        self.is_axisymmetric = getattr(self.input, 'is_axisymmetric', False)
        # Plasticity Properties
        self.is_plastic = getattr(self.input, 'is_plastic', False)
        if self.is_plastic:
            self.plastic_curve = self.input.plastic_curve
            self.epv = np.zeros(self.input.num_elems)

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
        self.lyy = np.zeros(self.n_elem)       
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
        self.dzz = np.zeros(self.n_elem) 
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
        self.q_bv = np.zeros(self.n_elem)
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
        self.u_bc = self.input.u_bc
        self.v_bc = self.input.v_bc
        self.a_bc = self.input.a_bc
        self.f_bc = self.input.f_bc

        ## Time
        self.n = 0
        self.t = 0.0
        self.tfinal = self.input.tfinal
        self.Co = self.input.Co
        self.dt = self.Co * (np.min(self.clen) / self.waves)
        self.t_prev = 0.0
        self.dt_prev = 0.0

        ## Energies
        self.kinetic_energy = 0.0
        self.internal_energy = 0.0   
        self.external_energy = 0.0    

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
            node1 = self.conn[e][0] - 1
            node2 = self.conn[e][1] - 1
            node3 = self.conn[e][2] - 1
            node4 = self.conn[e][3] - 1

            # Calculate velocity gradients (L tensor)
            self.lxx[e] = self.v[node1][0] * self.dN1dx[e] + self.v[node2][0] * self.dN2dx[e] + self.v[node3][0] * self.dN3dx[e] + self.v[node4][0] * self.dN4dx[e]
            self.lxy[e] = self.v[node1][0] * self.dN1dy[e] + self.v[node2][0] * self.dN2dy[e] + self.v[node3][0] * self.dN3dy[e] + self.v[node4][0] * self.dN4dy[e]
            self.lyx[e] = self.v[node1][1] * self.dN1dx[e] + self.v[node2][1] * self.dN2dx[e] + self.v[node3][1] * self.dN3dx[e] + self.v[node4][1] * self.dN4dx[e]
            self.lyy[e] = self.v[node1][1] * self.dN1dy[e] + self.v[node2][1] * self.dN2dy[e] + self.v[node3][1] * self.dN3dy[e] + self.v[node4][1] * self.dN4dy[e]

            # Calculate rate of deformation tensor (D tensor)
            self.dxx[e] = self.lxx[e]
            self.dyy[e] = self.lyy[e]
            self.dxy[e] = 0.5 * (self.lxy[e] + self.lyx[e])

            # Calculate spin tensor (W tensor)
            self.wxx[e] = 0.0
            self.wyy[e] = 0.0
            self.wzz[e] = 0.5 * (self.lyx[e] - self.lxy[e])

            # Hoop strain rate for axisymmetric formulation
            if self.is_axisymmetric:
                x1 = self.coo[node1][0]
                x2 = self.coo[node2][0]
                x3 = self.coo[node3][0]
                x4 = self.coo[node4][0]
                r_c = 0.25 * (x1 + x2 + x3 + x4)
                v_xc = 0.25 * (self.v[node1][0] + self.v[node2][0] + self.v[node3][0] + self.v[node4][0])
                
                if r_c > 1e-8:
                    self.dzz[e] = v_xc / r_c
                else:
                    self.dzz[e] = self.dxx[e]
            else:
                self.dzz[e] = 0.0

    def calc_bulk_viscosity(self):
        self.q_bv.fill(0.0)
        for e in range(self.n_elem):
            tr_D = self.dxx[e] + self.dyy[e] + self.dzz[e]
            
            if tr_D < 0.0:
                L_e = np.sqrt(4.0 * self.detJ[e])
                self.q_bv[e] = self.rho * L_e * (self.C_q * L_e * (tr_D ** 2) - self.C_l * self.waves * tr_D)

    def matstatupd(self):
        dt = self.dt
        E = self.E
        nu = self.poisson
        mu = E / (2.0 * (1.0 + nu))
        lam = (nu * E) / ((1.0 + nu) * (1.0 - 2.0 * nu))
        
        lam_dt = lam * dt
        mu_dt = mu * dt

        for e in range(self.n_elem):
            wzz_inc = self.wzz[e] * dt 
            r1 = 2.0 * self.res_sxy[e] * wzz_inc

            sxx_rot = self.res_sxx[e] - r1
            syy_rot = self.res_syy[e] + r1
            szz_rot = self.res_szz[e]      
            sxy_rot = self.res_sxy[e] + wzz_inc * (self.res_sxx[e] - self.res_syy[e])

            de_xx = self.dxx[e]
            de_yy = self.dyy[e]
            de_zz = self.dzz[e]
            de_xy = self.dxy[e]

            tr_de = de_xx + de_yy + de_zz
            ds_xx = 2.0 * mu_dt * de_xx + lam_dt * tr_de
            ds_yy = 2.0 * mu_dt * de_yy + lam_dt * tr_de
            ds_zz = 2.0 * mu_dt * de_zz + lam_dt * tr_de
            ds_xy = 2.0 * mu_dt * de_xy

            self.res_sxx[e] = sxx_rot + ds_xx
            self.res_syy[e] = syy_rot + ds_yy
            self.res_szz[e] = szz_rot + ds_zz
            self.res_sxy[e] = sxy_rot + ds_xy

    def matstatupd_elasto_plastic(self):
        dt = self.dt
        E = self.E
        nu = self.poisson
        mu = E / (2.0 * (1.0 + nu))
        lam = (nu * E) / ((1.0 + nu) * (1.0 - 2.0 * nu))
        
        lam_dt = lam * dt
        mu_dt = mu * dt

        def get_Y_and_H(ep_val, curve):
            if ep_val <= curve[0, 0]:
                H = (curve[1, 1] - curve[0, 1]) / (curve[1, 0] - curve[0, 0]) if len(curve) > 1 else 0.0
                return curve[0, 1], H
            for i in range(len(curve) - 1):
                if curve[i, 0] <= ep_val <= curve[i+1, 0]:
                    H = (curve[i+1, 1] - curve[i, 1]) / (curve[i+1, 0] - curve[i, 0])
                    Y = curve[i, 1] + H * (ep_val - curve[i, 0])
                    return Y, H
            return curve[-1, 1], 0.0

        for e in range(self.n_elem):
            wzz_inc = self.wzz[e] * dt
            r1 = 2.0 * self.res_sxy[e] * wzz_inc
            
            sxx_rot = self.res_sxx[e] - r1
            syy_rot = self.res_syy[e] + r1
            szz_rot = self.res_szz[e]
            sxy_rot = self.res_sxy[e] + wzz_inc * (self.res_sxx[e] - self.res_syy[e])

            de_xx = self.dxx[e]
            de_yy = self.dyy[e]
            de_zz = self.dzz[e] 
            de_xy = self.dxy[e]

            tr_de = de_xx + de_yy + de_zz
            ds_xx = 2.0 * mu_dt * de_xx + lam_dt * tr_de
            ds_yy = 2.0 * mu_dt * de_yy + lam_dt * tr_de
            ds_zz = 2.0 * mu_dt * de_zz + lam_dt * tr_de
            ds_xy = 2.0 * mu_dt * de_xy

            sxx_tr = sxx_rot + ds_xx
            syy_tr = syy_rot + ds_yy
            szz_tr = szz_rot + ds_zz
            sxy_tr = sxy_rot + ds_xy

            p_tr = (sxx_tr + syy_tr + szz_tr) / 3.0
            Sxx = sxx_tr - p_tr
            Syy = syy_tr - p_tr
            Szz = szz_tr - p_tr
            Sxy = sxy_tr

            q_tr = np.sqrt(1.5 * (Sxx**2 + Syy**2 + Szz**2 + 2.0 * Sxy**2))
            
            epv_curr = self.epv[e]
            Y, H = get_Y_and_H(epv_curr, self.plastic_curve)

            # Yield Condition
            if q_tr <= Y + 1e-8:
                # Elastic step
                self.res_sxx[e] = sxx_tr
                self.res_syy[e] = syy_tr
                self.res_szz[e] = szz_tr
                self.res_sxy[e] = sxy_tr
            else:
                # Plastic step: Radial Return Newton-Raphson Solver
                dg = 0.0
                for _ in range(20):
                    Y_temp, H_temp = get_Y_and_H(epv_curr + dg, self.plastic_curve)
                    residual = q_tr - 3.0 * mu * dg - Y_temp
                    if abs(residual) < 1e-6 * Y_temp:
                        break
                    derivative = -3.0 * mu - H_temp
                    dg -= residual / derivative
                    
                self.epv[e] += dg
                Y_final, _ = get_Y_and_H(self.epv[e], self.plastic_curve)
                
                factor = Y_final / q_tr
                Sxx *= factor
                Syy *= factor
                Szz *= factor
                Sxy *= factor
                
                self.res_sxx[e] = Sxx + p_tr
                self.res_syy[e] = Syy + p_tr
                self.res_szz[e] = Szz + p_tr
                self.res_sxy[e] = Sxy

    def assmb_internal(self):
        for e in range(self.n_elem):
            # Get the indices of the nodes for the current element
            node1, node2, node3, node4 = [self.conn[e][i] - 1 for i in range(4)]
            area = 4.0 * self.detJ[e]

            if self.is_axisymmetric:
                x1 = self.coo[node1][0]
                x2 = self.coo[node2][0]
                x3 = self.coo[node3][0]
                x4 = self.coo[node4][0]
                r_c = 0.25 * (x1 + x2 + x3 + x4)
                vol_weight = 2.0 * np.pi * r_c * area
                
                hoop_stress = self.res_szz[e] - self.q_bv[e]
                hoop_force = 2.0 * np.pi * hoop_stress * area * 0.25
            else:
                vol_weight = area * self.thickness
                hoop_force = 0.0
                
            sxx = (self.res_sxx[e] - self.q_bv[e]) * vol_weight
            syy = (self.res_syy[e] - self.q_bv[e]) * vol_weight
            sxy = self.res_sxy[e] * vol_weight

            # Contribution to internal force from element stress and shape function derivatives
            self.f_int[node1][0] += sxx * self.dN1dx[e] + sxy * self.dN1dy[e] + hoop_force
            self.f_int[node1][1] += sxy * self.dN1dx[e] + syy * self.dN1dy[e]

            self.f_int[node2][0] += sxx * self.dN2dx[e] + sxy * self.dN2dy[e] + hoop_force
            self.f_int[node2][1] += sxy * self.dN2dx[e] + syy * self.dN2dy[e]

            self.f_int[node3][0] += sxx * self.dN3dx[e] + sxy * self.dN3dy[e] + hoop_force
            self.f_int[node3][1] += sxy * self.dN3dx[e] + syy * self.dN3dy[e]

            self.f_int[node4][0] += sxx * self.dN4dx[e] + sxy * self.dN4dy[e] + hoop_force
            self.f_int[node4][1] += sxy * self.dN4dx[e] + syy * self.dN4dy[e]

    def assmb_mass(self):
        for e in range(self.n_elem):
            node1, node2, node3, node4 = [self.conn[e][i] - 1 for i in range(4)]
            area = 4.0 * self.detJ[e]

            if self.is_axisymmetric:
                x1 = self.coo[node1][0]
                x2 = self.coo[node2][0]
                x3 = self.coo[node3][0]
                x4 = self.coo[node4][0]
                r_c = 0.25 * (x1 + x2 + x3 + x4)
                vol = 2.0 * np.pi * r_c * area
            else:
                vol = area * self.thickness
                
            total_mass = vol * self.rho
            nodal_mass = total_mass / 4.0
            
            nodes = [node1, node2, node3, node4]
            for node in nodes:
                self.mass[node][0] += nodal_mass
                self.mass[node][1] += nodal_mass

    def el_state_upd(self):
        self.el_geom()
        self.el_rate()
        self.calc_bulk_viscosity()
        if self.is_plastic:
            self.matstatupd_elasto_plastic()
        else: # elastic
            self.matstatupd()
        self.assmb_internal()
        if self.t == 0.0:
            self.assmb_mass()

    def assmb_ubcs(self, t):
        if self.u_bc:
            for index, displacements in zip(self.u_bc.indexes, self.u_bc.displacements):
                for dof, displacement in enumerate(displacements):
                    if displacement is not None or 0:
                        self.u[index - 1][dof] = displacement(t)
                        self.u_inc[index - 1][dof] = self.u[index - 1][dof] - self.u_prev[index - 1][dof]
                        self.v[index - 1][dof] = (self.u[index - 1][dof] - self.u_prev[index - 1][dof]) / self.dt
                        self.a[index - 1][dof] = (self.v[index - 1][dof] - self.v_prev[index - 1][dof]) / self.dt

    def assmb_vbcs(self, t):
        if not self.v_bc:
            return
        for index, velocities in zip(self.v_bc.indexes, self.v_bc.velocities):
            for dof, velocity_callable in enumerate(velocities):
                if velocity_callable is None:
                    continue
                val = velocity_callable(t) 
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
        self.assmb_ubcs(self.t + self.dt)
        self.coo = self.coo + self.u_inc
        self.n += 1
        self.t += self.dt
        self.f_int.fill(0)