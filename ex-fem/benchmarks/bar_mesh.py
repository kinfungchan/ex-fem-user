import numpy as np

def generate_mesh(num_elems, min_x, length_x, width_y):
    # Number of nodes in the x-direction and y-direction
    num_nodes_x = num_elems + 1

    # Generate coordinates
    x_coords = np.linspace(min_x, min_x + length_x * num_elems, num_nodes_x)
    y_coords = np.array([0.0, width_y])

    # Create the coordinates array
    coordinates = np.array([[x, y] for y in y_coords for x in x_coords])

    # Generate connectivity
    connectivity = []
    for i in range(num_elems):
        n1 = i + 1
        n2 = n1 + 1
        n3 = n2 + num_nodes_x
        n4 = n1 + num_nodes_x
        connectivity.append([n1, n2, n3, n4])

    connectivity = np.array(connectivity)
    
    return coordinates, connectivity

def generate_bar_mesh(num_elems_x, num_elems_y, min_x, min_y, length_x, length_y, direction):
    # Number of nodes in the x-direction and y-direction
    num_nodes_x = num_elems_x + 1
    num_nodes_y = num_elems_y + 1

    # Generate coordinates
    x_coords = np.linspace(min_x, min_x + length_x * num_elems_x, num_nodes_x)
    y_coords = np.linspace(min_y, min_y + length_y * num_elems_y, num_nodes_y)

    # Create the coordinates array
    coordinates = np.array([[x, y] for y in y_coords for x in x_coords])

    # Generate connectivity
    connectivity = []
    for j in range(num_elems_y):
        for i in range(num_elems_x):
            n1 = i + j * num_nodes_x + 1
            n2 = n1 + 1
            n3 = n2 + num_nodes_x
            n4 = n1 + num_nodes_x
            connectivity.append([n1, n2, n3, n4])

    connectivity = np.array(connectivity)

    print("Coordinates:")
    print(coordinates)
    print("\nConnectivity:")
    print(connectivity)
    
    return coordinates, connectivity