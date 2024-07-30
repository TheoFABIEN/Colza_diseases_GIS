def mesh_creator(x1, x2, y1, y2, mesh_size, spacing = 0):
    """
    Compute the geopgraphic coordinates of all the squares composing the mesh,
    expanding around an initial square. 
    Args:
    x1 (float): x coordinate of the left side of the initial selected area.
    x2 (float): x coordinate of the right side of the initial selected area.
    y1 (float): y coordinate of the bottom side of the initial selected area.
    y2 (float): y coordinate of the top side of the initial selected area.
    mesh_size (int): The side length of the mesh. Needs to be an odd number, 
                     otherwise the function raises a value error. 
    spacing (float): The spacing between adjacent squares. Default is 0.
    Returns:
    surrounding_areas (list): List of tuples containing the (x1, x2, y1, y2)
                              coordinates of all the squares in the mesh.

    """

    if mesh_size % 2 == 0:
        raise ValueError('The size of the mesh should be an odd number.')

    # Calculate the dimentions of the squares
    square_width = x2 - x1
    square_height = y2 - y1

    half_size = mesh_size // 2

    all_areas_coordinates = []

    # Loop through the range to get all squares in the mesh:
    for i in range(-half_size, half_size + 1):
        for j in range(-half_size, half_size + 1):
            adj_x1 = x1 + i * (square_width + spacing)
            adj_x2 = x2 + i * (square_width + spacing)
            adj_y1 = y1 + j * (square_height + spacing)
            adj_y2 = y2 + j * (square_height + spacing)
            all_areas_coordinates.append((adj_x1, adj_x2, adj_y1, adj_y2))

    return all_areas_coordinates 
