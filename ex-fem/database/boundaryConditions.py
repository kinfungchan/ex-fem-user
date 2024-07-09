import numpy as np
"""
This module is to apply Boundary conditions to the domain
"""

class VelBoundaryConditions:
    def __init__(self, indexes: np.array, velocities: np.array):
        self.indexes = indexes
        self.velocities = velocities

    def velbc(t, L, E, rho):
        sinePeriod = (L / 2) * np.sqrt(rho / E)
        freq = 1 / sinePeriod
        if t >= sinePeriod * 0.5:
            return 0
        else:
            return 0.01 * np.sin(2 * np.pi * freq * t)

class AccelBoundaryConditions:
    def __init__(self, indexes: np.array, accelerations: np.array):
        self.indexes = indexes
        self.accelerations = accelerations

class ForceBoundaryConditions:
    def __init__(self, indexes: np.array, forces: np.array):
        self.indexes = indexes
        self.forces = forces

class SupportBoundaryConditions:
    def __init__(self, indexes: np.array, supports: np.array):
        self.indexes = indexes
        self.supports = supports