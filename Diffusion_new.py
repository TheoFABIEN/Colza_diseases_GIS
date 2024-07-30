import os 
import random
import shutil
import warnings

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from parallelbar import progress_starmap
from parallelbar.tools import cpu_bench
from rasterstats import zonal_stats
import scipy.fftpack
import scipy.special
from pyogrio import read_dataframe
from scipy.ndimage import zoom
from tqdm import tqdm
import xarray as xr

from multiprocessing import Pool

import utils

try:
    import pyfftw
    scipy.fftpack = pyfftw.interfaces.scipy_fftpack
    print('🚀 pyfftw library found! Jumping to lightspeed...')
except ImportError:
    print('pyfftw not found. Using fftconvolve() from scipy.signal instead')


warnings.filterwarnings("ignore")




def place_colza_for_year2(gdf_year2, prop_random_fields, seed):
    """
    Prepares the GeoDataFrame for the next year by converting mustard fields to 
    colza (if prop_random_fields is not provided) or selecting a random 
    set of fields that will be considered colza fields for the simulation. By
    setting the random seed here, we ensure that with the same input parameters
    we will always get the same field configuration, allowing us to test 
    multiple kernels in indentical conditions. 

    Parameters:
    gdf_year2 (geopandas.GeoDataFrame): GeoDataFrame for the year n+1.
    prop_random_fields (float, optional): Proportion of random fields to be 
    converted to colza. If None, only does mustard to colza conversion.
    seed (float): The random seed to use.

    Returns:
    geopandas.GeoDataFrame: Modified GeoDataFrame for the next year with 
    updated 'CodeAgg' values to reflect colza field conversions.
    """
    res_df = gdf_year2.copy()
    if not prop_random_fields:
        res_df.loc[res_df['CodeAgg'] == 'Mou', 'CodeAgg'] = 'Col'
    else:
        mask = res_df['CodeAgg'].isin(['Ann', 'Col', 'Mou'])
        filtered_df = res_df[mask]
        try:
            random.seed(seed)  
        except NameError:
            'Seed was not defined'
        sel_indices = random.sample(
            sorted(filtered_df.index),
            k = int(len(filtered_df) * prop_random_fields)
        )
        res_df['CodeAgg'] = np.where(
            res_df.index.isin(sel_indices),
            'Col', 'not_colza'
        )
    return res_df 



def complete_analysis(emission_value, mu, beta2, prop_random_fields = None, 
                      seed = None):


    list_of_metrics = []    # Storage for all shapefiles 


    for shp_path in tqdm(shapefiles, total = len(shapefiles)): 
    
        gdf = read_dataframe(shp_path)
        
        # Convert mustard fields to colza fields:
        gdf.loc[gdf['CodeAgg'] == 'Mou', 'CodeAgg'] == 'Col'

        # Fields borders ('BOR') are considered annual. Let's correct this:
        gdf.loc[gdf['CODE_CULTU'] == 'BOR', 'CodeAgg'] = 'Per'
    
        colza_in_gdf = not gdf[gdf['CodeAgg'] == 'Col'].empty
        last_year = 'year8' in shp_path

        # Compute metrics relative to the content of the GeoDataFrame 
        gdf_metrics = utils.calculate_metrics(gdf, round_val = 5)
        gdf_metrics = {
            'year_n_'+k: v for k, v in gdf_metrics.items()
        }
        
        if colza_in_gdf:
    
            # Rasterize geodataframe
            out_raster = utils.CustomRaster(
                gdf,
                emission = emission_value,
                type_of_plant = 'Col',
                resolution = (- RESOLUTION, RESOLUTION)
            )
            out_grid = out_raster.geocube
    
            # Create convolution kernel
            resized_area_length = max(out_grid.sizes['x'], out_grid.sizes['y'])
            x = y = np.linspace(
                -2 * resized_area_length,
                2 * resized_area_length,
                num = 2 * resized_area_length + 1
            )
            z = utils.dispersion_kernel(x, y, mu = mu, beta2 = beta2)
            
            # Convolve using the diffusion kernel
            def apply_convolution(array):
                return scipy.signal.fftconvolve(array, z, mode = 'same')
            convolved_raster = out_grid.apply(apply_convolution)
    
            # Total number of emitted spores
            total_emitted = convolved_raster.NumValue.sum().item()
    
            if not last_year:
    
                #### Find and prepare GeoDataFrame for next year ####
    
                num_part = shp_path.rstrip('0123456789')
                num = int(shp_path[len(num_part):]) + 1
                shp_year2_path = num_part + str(num)
    
                gdf_year2_raw = read_dataframe(shp_year2_path) 

                # Same correction as for gdf:
                gdf_year2_raw.loc[
                    gdf_year2_raw['CODE_CULTU'] == 'BOR', 'CodeAgg'
                ] = 'Per'

                # Compute metrics relative to the content of gdf_year2
                gdf_yr2_metrics = utils.calculate_metrics(
                    gdf_year2_raw, 
                    round_val = 5
                )
                gdf_yr2_metrics = {
                    'year_n+1_'+k: v for k, v in gdf_yr2_metrics.items()
                }

                if analysis == '2':
                    prop_random_fields = int(gdf_yr2_metrics[
                        'year_n+1_Number_colza_fields'
                    ]) / (len(gdf_year2_raw)-1) # -1 because first row is domain



                for seed in seedlist:

                    gdf_year2 = place_colza_for_year2(
                        gdf_year2_raw, prop_random_fields, seed
                    )

                    colza_in_gdf_year2 = not gdf_year2.loc[
                        gdf_year2['CodeAgg'] == 'Col'
                    ].empty

                    base_metrics = {
                        "gdf_year_n": os.path.basename(shp_path),
                        "gdf_year_n+1": os.path.basename(shp_year2_path),
                        "seed": seed,
                        "emission value": emission_value,
                        "mu": mu,
                        "beta2": beta2,
                        "proportion_random_colza_fields": prop_random_fields,
                        "total_spores_emitted": total_emitted,
                    }
               
    
                    # Compute contiguity metrics:

                    contiguity_metrics = utils.contiguity_metrics(
                        gdf[gdf['CodeAgg'] == 'Col'],
                        gdf_year2[gdf_year2['CodeAgg'] == 'Col'],
                        buffer_distance = 5,    # You can experiment with this
                        round_val = 5,
                        crs = gdf.crs
                    )

                    if colza_in_gdf_year2:
                       
                        """
                        For zonal_stats to work we need to pass the numpy array
                        from the raster and the affine transform mapping its 
                        coordinates to the crs:
                        """
                        raster_numvalue = convolved_raster.NumValue
                        raster_array = raster_numvalue.values
                        affine_transform = raster_numvalue.rio.transform()
                        crs = raster_numvalue.rio.crs

                        stats = zonal_stats(
                            gdf_year2.loc[gdf_year2['CodeAgg'] == 'Col'],
                            raster_array,
                            affine = affine_transform,
                            crs = crs,
                            stats = ['sum', 'mean', 'count', 'std'],
                            all_touched = all_touched
                        )
                        
                        # Rescale 'mean' and 'std' values to match resolution:
                        stats = utils.rescale_stats(stats, RESOLUTION)
    
                        # Total number and proportion of intercepted spores
                        total_intercepted = sum([dico['sum'] for dico in stats])
                        prop_intercepted = total_intercepted / total_emitted

                        # Compute global stats (across polygons)
                        single_geom = gdf_year2[
                            gdf_year2['CodeAgg'] == 'Col'
                        ].unary_union
                        single_geom_gdf = gpd.GeoDataFrame(
                            geometry = [single_geom],
                            crs = crs
                        )

                        stats_global = zonal_stats(
                            single_geom_gdf,
                            raster_array,
                            affine = affine_transform,
                            crs = crs,
                            stats = ['mean', 'std'],
                            all_touched = all_touched
                        ) 

                        stats_global = utils.rescale_stats(stats_global, RESOLUTION)

                        mean_intercepted_per_pixel = stats_global[0]['mean']

                        # std of number of intercepted spores per pixels
                        std_intercepted_per_pixel = stats_global[0]['std']

                        current_metrics = {
                            **base_metrics,
                            "total_spores_intercepted": total_intercepted,
                            "prop_spores_intercepted": prop_intercepted,
                            "mean_intercepted_spores_per_pixels": mean_intercepted_per_pixel,
                            "std_intercepted_spores_per_pixel": std_intercepted_per_pixel 
                        }
                        current_metrics = {
                            **current_metrics,
                            **contiguity_metrics,
                            **gdf_metrics,
                            **gdf_yr2_metrics
                        }
                        list_of_metrics.append(current_metrics)

                    else:       # Case where gdf_year2 contains no colza
                    
                        current_metrics = {
                            **base_metrics,
                            "total_spores_intercepted": 0,
                            "prop_spores_intercepted": 0,
                            "mean_intercepted_spores_per_pixels": 0,
                            "std_intercepted_spores_per_pixel": 0 
                        }
                        current_metrics = {
                            **current_metrics,
                            **contiguity_metrics,
                            **gdf_metrics,
                            **gdf_yr2_metrics
                        }
                        list_of_metrics.append(current_metrics)
            

            else:    # Shapefile is last of the sequence

                current_metrics = {
                    "gdf_year_n": os.path.basename(shp_path),
                    "gdf_year_n+1": None, 
                    "seed": None,
                    "emission value": emission_value,
                    "mu": mu,
                    "beta2": beta2,
                    "proportion_random_colza_fields": prop_random_fields, 
                    "total_spores_emitted": total_emitted,
                    "total_spores_intercepted": None,
                    "prop_spores_intercepted": None,
                    "mean_intercepted_spores_per_pixels": None,
                    "std_intercepted_spores_per_pixel": None
                }
                contiguity_metrics = {
                    'boundary_area': None, 'jaccard_index': None,
                    'mean_distance_centroids': None, 
                    'std_distance_centroids': None,
                    'mean_distance_edges': None, 'std_distance_edges': None,
                    'mean_minimal_distances': None,
                    'std_minimal_distances': None
                }
                gdf_yr2_metrics = {
                    "year_n+1_Number_colza_fields": None,
                    "year_n+1_Proportion_of_fields_that_contain_colza": None,
                    "year_n+1_Total_surface_colza": None,
                    "year_n+1_Surface_Ann": None,
                    "year_n+1_Surface_Per": None,
                    "year_n+1_Surface_unknown": None
                }
                current_metrics = {
                    **current_metrics,
                    **contiguity_metrics,
                    **gdf_metrics,
                    **gdf_yr2_metrics
                }
                list_of_metrics.append(current_metrics)


    # Export results
    metrics_df = pd.DataFrame(list_of_metrics)
    return metrics_df 

                 
    
    
if __name__ == '__main__': 

#==============================================================================
#                           GENERAL PARAMETERS  
#==============================================================================
       
    RESOLUTION = 10  # Downgrade value i.e. length of a pixel in the raster (m)
    

#==============================================================================
#                           DIFFUSION PARAMETERS 
#==============================================================================

    analysis = input('Choose the type of analysis you want:\n\n \
                     1 -> Analysis on a given set of parameters\n \
                     2 -> Analysis on a given set of paremeters + random fields\n \
                     3 -> Analysis of sensitivity with latin hypercube\n')
   
    print(f"analysis: {analysis}")
    
    NUM_SEEDS = 5    # For analysis 2 or 3

#_____Sensitivity analysis_____________________________________________________
    if analysis == '3':           

        seedlist = list(range(NUM_SEEDS)) 

        parameters_ranges = np.array([
            [1, 200],           # Range for emission
            [10, 1000],         # Range for mu 
            [.2, 1.8],          # Range for Beta 2
            [.1, .5]            # Proportion of colza fields sampled randomly
        ])
        SAMPLE_SIZE = 12 
        parameters = utils.scaled_latin_hypercube(
            parameters_ranges,
            sample_size = SAMPLE_SIZE 
        )

#______Test on a given set of parameters_______________________________________
    elif analysis == '1':
        parameters = [
            [50, 250, .2]
        ]
        seedlist = [1]     # Or your favourite number

#_____Test on a given set of parameters + colza fields shuffled randomly_______
    else:
        parameters = [
            [50, 250, .2]
        ]
        seedlist = list(range(NUM_SEEDS)) 
    


#==============================================================================
#               PROGRAM EXECUTION (no parameters to adjust here)
#==============================================================================

    all_touched = False
    if RESOLUTION > 1:
        all_touched = True
    
    # Folder selection
    selected_folder_path = utils.select_folder()
    results_path = f'{selected_folder_path}/0RESULTS_DIFFUSION'
    if os.path.exists(results_path):
        shutil.rmtree(results_path)
    os.mkdir(results_path)
   
    # Parallel execution of the function
    shapefiles = utils.find_shapefiles(selected_folder_path)
    shapefiles = [file for file in shapefiles if "UPDATED_FILES_" in file]
    #results = progress_starmap(complete_analysis, parameters)
    with Pool() as p:
        results = p.starmap(complete_analysis, parameters)
    results_concat = pd.concat(results)
    results_concat = results_concat.sort_values(
        by = ['gdf_year_n', 'emission value', 'mu'],
        inplace = False
    )
    results_concat.to_csv(f'{results_path}/metrics_df.csv')


