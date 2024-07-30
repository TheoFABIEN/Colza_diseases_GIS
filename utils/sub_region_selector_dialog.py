import os
import shutil
import sys
import tkinter
from tkinter import filedialog


def draw_mesh(n):
    """
    Draws an ASCII representation of a square mesh.
    Args:
    n (int): number of squares per side.
    """
    # Define the basic elements
    horizontal_line = "+---" * n + "+"
    vertical_line = "|   " * n + "|"

    # Generate the full mesh
    for i in range(n):
        print(horizontal_line)
        print(vertical_line)
    print(horizontal_line)  # Print the last horizontal line


def dialog():
    """
    Prompt the user to select a directory containing shapefiles, choose a method
    for selecting the area of interest, choose how many batches are needed and
    choose if the newly created shapefiles should be saved on the computer.
    
    Returns:
    path (str): Path to the selected directory.
    coords (list of int): Coordinates of the area of interst, only if the user
                          did not select it on the interactive map.
    geo_selection_method (int): The method chosen for selecting the area of 
                               interest (1 for coordinate-based selection, 2 for
                               manual selection on the interactive map).
    num_batches (int): Size of the mesh containing areas of interest. If 1, only
                       1 square will be used.
    save_shapefile (bool): Whether to save the outputs as shapefiles.
    """
    
    # Open a file dialog to choose the folder containing the shapefiles:
    print(
            "\n Select the folder path containing shapefiles for a region \n"
    )
    tkinter.Tk().withdraw()
    path = filedialog.askdirectory()
    
    #if os.path.exists(path + '/UPDATED_FILES'):
    #    shutil.rmtree(path + '/UPDATED_FILES')
    #os.mkdir(path + '/UPDATED_FILES')
    
    geo_selection_method = input(
            "\n Select the area of interest : \n \n \
            From coordinates (LAMB93)      ->  1 \n \n \
            Manually on the map            ->  2 (default)\n \n" 
    )
    if geo_selection_method == "":
        geo_selection_method = 0
    geo_selection_method = int(geo_selection_method)
    
    x1, x2, y1, y2 = 0, 0, 0, 0
    if geo_selection_method == 1:
        print('Enter the coordinates of the upper left corner of the area: ')
        x1 = int(input('x : '))
        y1 = int(input('y : '))
        print("Enter the length of one side of the square (m): ")  
        l = int(input('Length: '))
        x2 = x1 + l
        y2 = y1 - l
       

    print("Mesh dimension (number of square lengths per side): \n \
    (Enter 1 to have only one area of interest)")
    mesh_size = int(input())
    
    print("Mesh shape: ")
    draw_mesh(mesh_size)

    reg_name = input('\nEnter region name (optional): ')
    if reg_name != '':
        reg_name += '_'

    return path, [x1, x2, y1, y2], geo_selection_method, mesh_size, reg_name
    
    
