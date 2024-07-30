import os 
import shutil 
import sys
import tkinter
from tkinter import filedialog
import warnings

import contextily as cx
import geopandas as gpd 
import matplotlib.pyplot as plt 
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.widgets import RectangleSelector
import matplotlib.widgets as widgets
import numpy as np
import pandas as pd
from pyogrio import read_dataframe, write_dataframe
from shapely.geometry import Polygon

import pyfftw 
import scipy.fftpack
from skimage.draw import polygon2mask
import xarray as xr

import utils


warnings.filterwarnings('ignore')

def prYellow(skk): print('\033[93m {}\033[00m' .format(skk)) #print in yellow

# Ask the user how the script should be used
path, [x1, x2, y1, y2], geo_selection_method, mesh_size, reg_name = (
    utils.dialog()
)

# Open the shapefiles
shapefiles = []
for r, d, f in os.walk(path):
    for file in f:
        if file.endswith('PARCELLES_GRAPHIQUES.shp'):
            shapefiles.append(os.path.join(r, file))
if shapefiles == []:
    prYellow("No shapefile found \n Terminating program...")
    sys.exit()

# Open the .csv that will be used to change crop labels
conversion_codes_cultures = pd.read_csv(
        os.getcwd() + '/Conversion_Codes_Cultures.csv'
).loc[:, ['CodeAgg', 'CODE_CULTURE']]


###############  FUNCTIONS FOR INTERACTIVE SELECTOR IN REGION MAP  ##########


def select_callback(eclick, erelease):
    """
    Callback for line selection
    *eclick* and *erelease* are the press and release events
    """
    global x1, y1, x2, y2 #set coordinates as global to use them later 
    x1, y1 = eclick.xdata, eclick.ydata 
    x2, y2 = erelease.xdata, erelease.ydata
    print(f"({x1:3.2f}, {y1:3.2f}) --> ({x2:3.2f}, {y2:3.2f})")

def toggleselector(event):
    name = type(selector).__name__
    print(f'{name} activated')
    selector.set_active(True)



###############  MAIN LOOP  #################################


if reg_name != '':
    updated_folder = path + f'/UPDATED_FILES_{reg_name[:-1]}'
else:
    updated_folder = path + '/UPDATED_FILES'
if os.path.exists(updated_folder):
    shutil.rmtree(updated_folder)
os.mkdir(updated_folder)

for i, f in enumerate(shapefiles):

    shapefile = read_dataframe(f)
    shapefile = shapefile.rename(columns = {'CODE_CULTU': 'CODE_CULTURE'})

    # Manual selection of the region of interest on the interactive map:
    if (geo_selection_method == 2 or geo_selection_method == 0) and i == 0:  

        fig, axs = plt.subplots()

        print('\n Shapefile loaded, building interactive map... \n')
        shapefile.plot(ax = axs)
        selector = RectangleSelector(
                axs,
                select_callback,
                useblit = True,
                button = [1, 3],
                minspanx = 0,
                minspany = 0,
                spancoords = 'data', 
                interactive = True
        )

        selector.add_state('square')

        cx.add_basemap(axs, crs = shapefile.crs.to_string())

        fig.canvas.mpl_connect('key_press_event', toggleselector)
        plt.title('Close this window after making a selection.')
        plt.show()

    #From this region of interest, compute coordinates for the whole mesh
    all_areas_coordinates = utils.mesh_creator(
        min(x1, x2), max(x1, x2), min(y1, y2), max(y1, y2),
        mesh_size,
        spacing = 0
    )

    for area_idx, area_coord in enumerate(all_areas_coordinates):
        
        area_x1, area_x2, area_y1, area_y2 = area_coord

        # Create new polygon matching the area of interest
        box_coords = (
            (min(area_x1, area_x2), min(area_y1, area_y2)),
            (min(area_x1, area_x2), max(area_y1, area_y2)),
            (max(area_x1, area_x2), max(area_y1, area_y2)),
            (max(area_x1, area_x2), min(area_y1, area_y2))
        )
        box_polygon_shape = Polygon(box_coords)
        box_polygon = gpd.GeoDataFrame(
            index = [0],
            crs = 'epsg:2154',
            geometry = [box_polygon_shape]
        )

        shapefile_area = gpd.GeoDataFrame.copy(shapefile)

        # Add this polygon to the shapefile
        shapefile_area = gpd.pd.concat([box_polygon, shapefile_area])

        # Only keep crops that fit in the big polygon
        shapefile_area = shapefile_area[
            shapefile_area.within(
                box_polygon.geometry.iloc[0]
            )
        ]

        # Create CodeAgg column
        CodeAgg_allrows = pd.merge(
            conversion_codes_cultures,
            shapefile_area,
            how = 'right',
            left_on = 'CODE_CULTURE',
            right_on = 'CODE_CULTURE'
        )['CodeAgg']

        CodeAgg_allrows = list(CodeAgg_allrows)
        shapefile_area.insert(2, "CodeAgg", CodeAgg_allrows)
       
        if not os.path.exists(updated_folder + f'/area_{area_idx + 1}'):
            os.mkdir(updated_folder + f'/area_{area_idx + 1}')
        os.mkdir(
            updated_folder + f'/area_{area_idx + 1}' + 
            f'/{reg_name}area{area_idx + 1}_year{i + 1}'
        )

        name = f'year_{i + 1}.shp'
        shapefile_area.to_file(
                updated_folder + f'/area_{area_idx + 1}' + 
                f'/{reg_name}area{area_idx + 1}_year{i + 1}/' +
                name,
                engine = 'pyogrio'
        )

        # Match colors with corresponding types of crops
        colors = [
                '#FFF305' if x == 'Col'        # Colza
                else '#E2FDD3' if x == 'Ann'   # Annual
                else '#005C00' if x == 'Per'   # Permanent
                else '#EDC9FB' if x == 'Mou'   # Mustard 
                else '#A7835A'#'red'                     # Unknown
                for x in shapefile_area['CodeAgg']
        ]

        #Create and save fig
        fig = shapefile_area.plot(
            column = shapefile_area['CodeAgg'], color = colors
        )
        fig.add_artist(ScaleBar(1))  
        fig.figure.savefig(
                updated_folder + f'/area_{area_idx + 1}/' + 
                reg_name + f'area{area_idx + 1}_year{i + 1}.png',
                bbox_inches = 'tight'
        )

        print(f'Zone {area_idx + 1} completed.', end = '\r')

    print(f'Task completed for file {i + 1} of {len(shapefiles)}.')


##########  EXPORT METRICS  ##############


areas_folders = sorted(os.listdir(updated_folder))
areas_folders = [
    file for file in areas_folders 
    if os.path.isdir(updated_folder + '/' + file)
] 

metrics_dataframes_list = []

for i, area in enumerate(areas_folders):

    metrics_dict = {} # Stores metrics for each year

    files = sorted(os.listdir(updated_folder + '/' + area))
    files = [
        item for item in files
        if os.path.isdir(updated_folder+ '/' + area + '/' + item)
    ]
    n_years = len(files)






    # Create list to store areas of contact zones
    boundary_areas_list = [None]



    contiguity_metrics_list = []



    for year, file in enumerate(files):

        shapefile = read_dataframe(
            f"{updated_folder}/{area}/{file}"
        ) 
        
        # Consider Mustard crops as colza crops 
        # (necessary for calculate_metrics() function to work)
        shapefile['CodeAgg'] = shapefile['CodeAgg'].replace('Mou', 'Col')
        metrics = utils.calculate_metrics(shapefile, round_val = 2)
        metrics_dict.update({f"year_{year + 1}": metrics})
        
        
        # If year > 0, we can compare with the previous year
        if year > 0:
            shapefile_prev_year = read_dataframe(
                f"{updated_folder}/{area}/{files[year - 1]}"
            )
            
            shapefile_prev_year['CodeAgg'] = shapefile_prev_year['CodeAgg'].replace(
                'Mou', 'Col'
            )
            
            # Get contact zones between colza crops of the 2 years 
            contiguity_metrics = utils.contiguity_metrics(
                shapefile_prev_year[shapefile_prev_year['CodeAgg'] == 'Col'],
                shapefile[shapefile['CodeAgg'] == 'Col'],
                buffer_distance = 5,
                round_val = 2,
                crs = shapefile.crs
            )

            contiguity_metrics_list.append(contiguity_metrics)

            boundary_area = round(contiguity_metrics['boundary_area'], 5)
            boundary_areas_list.append(boundary_area)
        
    contiguity_df = pd.DataFrame(contiguity_metrics_list)
    contiguity_df = contiguity_df.transpose()

    metrics_df = pd.DataFrame(metrics_dict)

    # Remove None for plotting
    boundary_areas_list = [0 if x == None else x for x in boundary_areas_list]


    # Add column in contiguity_df to match metrics_df:
    contiguity_df.insert(0, 'Year 0', 0)
    contiguity_df.columns = list(metrics_df.columns)

    metrics_df = pd.concat([metrics_df, contiguity_df], axis = 0)


    # Add the area of contact zones to metrics_df
    #boundary_areas_list = [0 if x == None else x for x in boundary_areas_list]
    #metrics_df.loc['Area_boundary_zones_with_prev_year'] = boundary_areas_list


    
    



    # Save the dataframe 
    metrics_df.to_csv(
        f"{updated_folder}/metrics_{reg_name}area{i + 1}",
        index = True
    )

    # Store it somewhere for later
    metrics_dataframes_list.append(metrics_df)

    ######### EXPORT PLOTS #################

    #plot_folder = f"{updated_folder}/plots"
    #if os.path.exists(plot_folder):
    #    shutil.rmtree(plot_folder)
    #os.mkdir(plot_folder)
    
    plot_folder = updated_folder + f'/area_{i + 1}/plots' 
    if os.path.exists(plot_folder):
        shutil.rmtree(plot_folder)
    os.mkdir(plot_folder)

    years = np.arange(1, len(files) + 1)

    # Number of colza fields plot 

    num_colza_fields = list(metrics_df.iloc[0])
    num_colza_fields = [int(x) for x in num_colza_fields]

    fig, ax = plt.subplots()
    ax.bar(years, num_colza_fields)
    fig.savefig(f"{plot_folder}/num_colza_fields_area{i + 1}.png")
    
    # Proportion of fields that contain colza
    
    prop_colza_fields = list(metrics_df.iloc[1])
    prop_colza_fields = [float(x) for x in prop_colza_fields]
    fig1, ax1 = plt.subplots()
    ax1.bar(years, prop_colza_fields)
    fig1.savefig(f"{plot_folder}/prop_fields_with_colza{i + 1}.png")

    # Surfaces plot 

    fig2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    plt.subplots_adjust(hspace = .5, wspace = .3)

    surf_colza = [float(x) for x in list(metrics_df.iloc[2])]
    surf_ann = [float(x) for x in list(metrics_df.iloc[3])]
    surf_per = [float(x) for x in list(metrics_df.iloc[4])]
    surf_unknown = [float(x) for x in list(metrics_df.iloc[5])]

    ax1.bar(years, surf_colza)
    ax1.set_title('Surface of colza')
    ax2.bar(years, surf_ann)
    ax2.set_title('Surface of annual')
    ax3.bar(years, surf_per)
    ax3.set_title('Surface of permanent')
    ax4.bar(years, surf_unknown)
    ax4.set_title('Surface of unknown')

    fig2.savefig(f"{plot_folder}/surfaces_area{i + 1}.png")

    # Area of contact zones plot

    fig, ax = plt.subplots()
    ax.bar(years, boundary_areas_list)
    fig.savefig(f"{plot_folder}/contact_zones_area{i + 1}.png")


###### Make general plots (across areas) ######

region_plot_folder = updated_folder + '/region_plots' 
if os.path.exists(region_plot_folder):
    shutil.rmtree(region_plot_folder)
os.mkdir(region_plot_folder)

# Number of colza fields 

colza_fields_general = [
    list(df.iloc[0].astype('int')) for df in metrics_dataframes_list
]
fig, ax = plt.subplots()
ax.boxplot(colza_fields_general)
ax.set_title('Number of colza fields depending on the area')
fig.savefig(region_plot_folder + '/number_colza_fields.png')

# Proportion of fields that contain colza

colza_prop_general = [
    list(df.iloc[1].astype('float')) for df in metrics_dataframes_list
]
fig, ax = plt.subplots()
ax.boxplot(colza_prop_general)
ax.set_title('Proportion of annual fields that contain colza')
fig.savefig(region_plot_folder + '/proportion_fields_with_colza.png')

# Total surfaces 

colza_surface_general = [
    list(df.iloc[2].astype('float')) for df in metrics_dataframes_list
]
ann_surface_general = [
    list(df.iloc[3].astype('float')) for df in metrics_dataframes_list
]
per_surface_general = [
    list(df.iloc[4].astype('float')) for df in metrics_dataframes_list
]
unknown_surface_general = [
    list(df.iloc[5].astype('float')) for df in metrics_dataframes_list
]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
plt.subplots_adjust(hspace = .5, wspace = .3)

ax1.boxplot(colza_surface_general)
ax1.set_title('Colza surface')
ax2.boxplot(ann_surface_general)
ax2.set_title('Annual surface')
ax3.boxplot(per_surface_general)
ax3.set_title('Permanent surface')
ax4.boxplot(unknown_surface_general)
ax4.set_title('Unknown surface')

fig.savefig(region_plot_folder + '/surfaces.png')

# Areas of boundary zones

boundary_surf_general = [
    list(df.iloc[6].astype('float')) for df in metrics_dataframes_list
]

fig, ax = plt.subplots()

ax.boxplot(boundary_surf_general)
ax.set_title('Boundary surface')

fig.savefig(region_plot_folder + '/boundary_surfaces.png')
