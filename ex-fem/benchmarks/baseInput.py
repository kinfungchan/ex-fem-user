"""
This class contains the parent class for the input file
"""

##Start of Input File 
class BaseInput:
    def __init__(self):
        # Material properties
        self.young = None  # Young's modulus
        self.density = None  # Density
        self.poisson = None  # Poisson's ratio
        self.num_elems = None  # Number of elements
        self.num_nodes = None  # Number of nodes

        # Mesh generation
        self.coordinates, self.connectivity = None, None # Coordinates and connectivity

        # Boundary conditions 
        self.v_bc = None # Velocity boundary conditions
        self.a_bc = None # Acceleration boundary conditions
        self.f_bc = None # Force boundary conditions
        self.s_bc = None # Support boundary conditions

        # Time variables
        self.tfinal = None  # Final time
        self.Co = None # Courant number

    def validate_input(self):
        """Validate the input parameters."""
        required_params = ["young", "density", "poisson", "num_elems", "num_nodes",
                           "coordinates", "connectivity", 
                           "Co", "tfinal"]
        for param in required_params:
            if getattr(self, param) is None:
                raise ValueError(f"Missing required parameter: {param}")