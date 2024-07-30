import tkinter as tk 
from tkinter import filedialog
import os

def select_folder():
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(master = root)
    root.destroy()
    return folder_path 
    

def find_shapefiles(directory):
    shapefiles = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.shp'):
                shapefiles.append(root)
    return shapefiles
