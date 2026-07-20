import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
import imageio 
import os

class Output:
    '''
    Output class to store results are each time step
    '''
    def __init__(self, t, coordinates, accel, vel, displ, sxx, syy, sxy):

        self.t = t
        self.coordinates = coordinates
        self.accel = accel
        self.vel = vel
        self.displ = displ
        self.sxx = sxx
        self.syy = syy
        self.sxy = sxy

class History:
    '''
    Output class to store results for ALL time steps
    '''
    def __init__(self, coordinates, n_nodes, n_elem):
        self.t = np.array([0.0])
        self.coordinates = np.array([coordinates], dtype=float)

        # Initialize 3D arrays for acceleration, velocity, and displacement
        self.accel = np.zeros((1, n_nodes, 2), dtype=float)
        self.vel = np.zeros((1, n_nodes, 2), dtype=float)
        self.displ = np.zeros((1, n_nodes, 2), dtype=float)

        # Initialize 2D arrays for stress components
        self.sxx = np.zeros((1, n_elem), dtype=float)
        self.syy = np.zeros((1, n_elem), dtype=float)
        self.sxy = np.zeros((1, n_elem), dtype=float)

    def append_timestep(self, t, coordinates, accel, vel, displ, sxx, syy, sxy):
        self.t = np.append(self.t, t)
        self.coordinates = np.vstack((self.coordinates, [coordinates]))
        
        # Append new time step data to 3D arrays
        self.accel = np.vstack((self.accel, [accel]))
        self.vel = np.vstack((self.vel, [vel]))
        self.displ = np.vstack((self.displ, [displ]))

        # Append new time step data to 2D arrays
        self.sxx = np.vstack((self.sxx, [sxx]))
        self.syy = np.vstack((self.syy, [syy]))
        self.sxy = np.vstack((self.sxy, [sxy]))

    def get_node_dof_data(self, node, dof):
        return self.accel[:, node, dof], self.vel[:, node, dof], self.displ[:, node, dof]

    def plot_node_dof(self, node, dof):

        accel_data, vel_data, displ_data = self.get_node_dof_data(node, dof)
        
        plt.figure(figsize=(12, 8))
        
        plt.subplot(3, 1, 1)
        plt.plot(self.t, accel_data)
        plt.title(f'Node {node}, DOF {dof} Acceleration vs Time')
        plt.xlabel('Time')
        plt.ylabel('Acceleration')
        
        plt.subplot(3, 1, 2)
        plt.plot(self.t, vel_data)
        plt.title(f'Node {node}, DOF {dof} Velocity vs Time')
        plt.xlabel('Time')
        plt.ylabel('Velocity')

        plt.subplot(3, 1, 3)
        plt.plot(self.t, displ_data)
        plt.title(f'Node {node}, DOF {dof} Displacement vs Time')
        plt.xlabel('Time')
        plt.ylabel('Displacement')

        plt.tight_layout()
        plt.show()

    def get_element_data(self, element):
        return self.sxx[:, element], self.syy[:, element], self.sxy[:, element]
    
    def plot_element_data(self, element):

        sxx_data, syy_data, sxy_data = self.get_element_data(element)

        plt.figure(figsize=(12, 8))
        
        plt.subplot(3, 1, 1)
        plt.plot(self.t, sxx_data)
        plt.title(f'Element {element} Sxx vs Time')
        plt.xlabel('Time')
        plt.ylabel('Sxx')
        
        plt.subplot(3, 1, 2)
        plt.plot(self.t, syy_data)
        plt.title(f'Element {element} Syy vs Time')
        plt.xlabel('Time')
        plt.ylabel('Syy')

        plt.subplot(3, 1, 3)
        plt.plot(self.t, sxy_data)
        plt.title(f'Element {element} Sxy vs Time')
        plt.xlabel('Time')
        plt.ylabel('Sxy')

        plt.tight_layout()
        plt.show()

    def plot_mesh_data(self, time_step, variables, dof, coords_list, conn_list, var_name,
                       save_path=None, clim=None, zoom=1.0):
        if time_step >= len(self.t):
            raise IndexError(f"time_step {time_step} is outside stored history")

        coords = np.asarray(coords_list[0])
        conn = np.asarray(conn_list[0], dtype=int) - 1
        variable = np.asarray(variables[0])

        if conn.shape[1] != 4:
            raise ValueError("plot_mesh_data expects 4-node quad elements")
        if coords.ndim != 3 or coords.shape[2] != 2:
            raise ValueError("plot_mesh_data expects coordinate history with shape (time, node, 2)")
        if variable.ndim != 3 or variable.shape[1] != coords.shape[1]:
            raise ValueError("plot_mesh_data expects nodal field history with shape (time, node, dof)")

        xy = coords[time_step]
        values = variable[time_step, :, dof]

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.style.use('ggplot')

        triangles = np.vstack((conn[:, [0, 1, 2]], conn[:, [0, 2, 3]]))
        triangulation = mtri.Triangulation(xy[:, 0], xy[:, 1], triangles)
        mappable = ax.tripcolor(
            triangulation,
            values,
            shading='gouraud',
            cmap='RdYlGn_r',
            vmin=None if clim is None else clim[0],
            vmax=None if clim is None else clim[1],
        )

        edges = []
        for element in conn:
            nodes = xy[element]
            closed = np.vstack((nodes, nodes[0]))
            edges.extend(zip(closed[:-1], closed[1:]))
        ax.add_collection(LineCollection(edges, colors='white', linewidths=0.1, alpha=0.45))

        colorbar = fig.colorbar(mappable, ax=ax)
        colorbar.set_label(var_name)

        ax.set_aspect('equal', adjustable='box')
        ax.set_title(f"{var_name} - Time step {time_step} - Time: {self.t[time_step]:.4e}")
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        xmin, ymin = np.min(xy, axis=0)
        xmax, ymax = np.max(xy, axis=0)
        xmid = 0.5 * (xmin + xmax)
        ymid = 0.5 * (ymin + ymax)
        width = xmax - xmin
        height = ymax - ymin
        if zoom and zoom > 0:
            width /= zoom
            height /= zoom
        ax.set_xlim(xmid - 0.5 * width, xmid + 0.5 * width)
        ax.set_ylim(ymid - 0.5 * height, ymid + 0.5 * height)

        if save_path:
            fig.savefig(save_path, dpi=200, bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

class Animation:
    '''
    Output class to animate results through time
    '''
    def __init__(self, History: History, directory):
        self.P = History
        self.directory = directory
        # Ensure the directory exists
        os.makedirs(self.directory, exist_ok=True)

        self.filenames_accel = []
        self.filenames_vel = []
        self.filenames_disp = []
        self.filenames_stress = []
        self.filenames_bv = []

        self.filenames_mesh_accel = []
        self.filenames_mesh_vel = []
        self.filenames_mesh_disp = []
        self.filenames_mesh_stress = []
        self.filenames_mesh_accel_pv = []
        self.filenames_mesh_vel_pv = []
        self.filenames_mesh_displ_pv = []
        self.filenames_mesh_stress_pv = []
     
    def save_single_plot(self, n_plots, x, y, title, xlabel, ylabel, filenames, n, t):
        filename = f'{self.directory}/FEM1D_{title}{n}.png'
        filenames.append(filename)
        plt.figure(figsize=(10, 6))
        plt.style.use('ggplot')
        # for i in range(n_plots):
        #         plt.plot(x[i], y[i])
        plt.plot(x,y)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.ylim(-0.002, 0.012)
        plt.legend([f"{t}"])
        plt.savefig(filename)
        plt.close()

    def save_mesh_plot(self, coordinates, connectivity, nodal_field, dof, title, ylabel, filenames, n, t):
        filename = f'{self.directory}/FEM2D_{title}{n}.png'
        filenames.append(filename)
        plt.figure(figsize=(12, 6))
        plt.style.use('ggplot')
        plt.tripcolor(coordinates[:, 0], coordinates[:, 1], nodal_field[:, dof], shading='gouraud')
        plt.colorbar(label=ylabel)
        plt.xlabel('Position')
        plt.ylabel(ylabel)
        plt.title(title)
        ## Plot element boundaries
        # for element in connectivity:
        #     plt.plot(coordinates[element-1, 0], coordinates[element-1, 1], 'k-')
        plt.legend([f"{t}"])
        plt.savefig(filename)
        plt.close()

    def save_timesteps(self, start_time, end_time, step, variables, dof, coords_list, conn_list,
                       var_name, filenames, *, clim=None, zoom=1.0, prefix="frame"):
        for time_step in range(start_time, end_time + 1, step):
            save_path = os.path.join(self.directory, f"{prefix}_{time_step:04d}.png")
            try:
                self.P.plot_mesh_data(
                    time_step,
                    variables,
                    dof,
                    coords_list,
                    conn_list,
                    var_name,
                    save_path=save_path,
                    clim=clim,
                    zoom=zoom,
                )
            except IndexError:
                break
            filenames.append(save_path)
        return filenames
            
    def create_gif(self, gif_name, filenames):
        gif_path = os.path.join(self.directory, gif_name)
        with imageio.get_writer(gif_path, mode='I') as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
        for filename in set(filenames):
            os.remove(filename)