import os 
import tkinter as tk
from tkinter import filedialog

import folium
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
from pyogrio import read_dataframe
import scipy.fftpack 
from scipy.ndimage import zoom
from skimage.draw import polygon2mask
import streamlit as st
from streamlit_folium import st_folium
import xarray as xr 
from rasterstats import zonal_stats

import utils 


#  fftw is a nightmare to install on windows, so we make it optional:
try:
    import pyfftw
    scipy.fftpack = pyfftw.interfaces.scipy_fftpack
except ImportError:
    print('pyfftw not found. Using fftconvolve() from scipy.signal instead.')

    
st.set_page_config(
    page_title = "Colza Disease Risk Tool",
    layout = "wide",
    initial_sidebar_state="expanded",
    page_icon = ":seedling:"
)

map_tab, parameters_tab, metrics_tab = st.tabs(
    [r":world_map: $\textsf{\large ~ Interactive map}$", 
     r":gear: $\textsf{\large ~ Diffusion parameters}$",
     r":bar_chart: $\textsf{\large ~ Metrics}$"]
)

with parameters_tab:   
    _, main_column, _ = st.columns([.5, 1, 1])
    with main_column: 
        st.write('#')
        mu = st.number_input(
            r'$\textsf{Mean dispersion distance (mu)}$',
            value = 250 
        )  
        st.write('##')
        beta2 = st.number_input(
            r'$\textsf{Shape parameter (beta 2)}$',
            value = .2
        )
        st.write('##')
        intensity_value = st.number_input(
            r'$\textsf{Initial intensity of infection}$',
            value = 50
        )


tk.Tk().withdraw()  # Hide tkinter root window

# Initialize session_state for folder_path
if 'folder_path' not in st.session_state:
    st.session_state['folder_path'] = []


with st.sidebar:
    
    _, title_column, _ = st.columns([1, 20, 1])  # Center title
    with title_column:
        st.title(':orange[**Colza Disease Risk Tool**]')

    st.markdown('***')   # Add space

    sel_fold_button = st.button('Select shapefiles sequence folder')

    if sel_fold_button:
        st.session_state['folder_path'] = utils.select_folder()
        st.write('Selected folder :  \n ', st.session_state['folder_path'])

    if not st.session_state['folder_path']:
        st.stop()
        
    sel_fold_path = st.session_state['folder_path']

    # Get all files and folders located in sel_fold_path, sort them
    files = sorted(os.listdir(sel_fold_path))
    # Exclude elements that are not folders (mainly to remove .png images) 
    files = [
        item for item in files
        if os.path.isdir(os.path.join(sel_fold_path, item)) and item != 'plots'
    ]
    n_years = len(files)
    year_str_seq = [i+1 for i in range(n_years)]
    year_str_seq.pop(0)  # Remove 1st element as we can't get risk for 1st year
    year_str_seq.append(f'{max(year_str_seq) + 1} (no data for this year)')
    st.write(str(n_years), 'years sequence found.')
    
    selected_year = st.selectbox(
    '**Show risk map for year:**',
    year_str_seq
    )
   
    st.markdown('###')
    st.markdown('***')

    with st.expander('### Tutorial'):
        st.write(utils.tutorial)

    with st.expander('### Longer tutorial'):
        st.write(utils.longer_tutorial)


#### Read the GeoDataFrame for the selected year ####


quality_reduction = 2   # Length of a pixel (m)


all_touched = False
if quality_reduction > 2:
    all_touched = True



predictive = False
if type(selected_year) is str:
    predictive = True
    selected_year = int(selected_year[0])
shp_path = sel_fold_path + '/' + files[selected_year-2]
gdf = read_dataframe(shp_path)

# Convert mustard fields to colza fields:
gdf.loc[gdf['CodeAgg'] == 'Mou', 'CodeAgg'] = 'Col'

# Convert it to raster
out_raster = utils.CustomRaster(
    gdf, 
    emission = 50,
    type_of_plant = 'Col', 
    resolution = (-quality_reduction, quality_reduction)
)


#### Open and prepare the GeoDataFrame for upcoming year ####

if not predictive:
    shp_year2_path = sel_fold_path + '/' + files[selected_year-1]
    gdf_year2 = read_dataframe(shp_year2_path)
if predictive:
    gdf_year2 = gdf 

# Convert mustard fields to colza fields
gdf_year2.loc[gdf_year2['CodeAgg'] == 'Mou', 'CodeAgg'] = 'Col'

gdf_year2['Domain_polygon'] = [
    0 if gdf_year2['ID_PARCEL'][i] is None
    else 1
    for i in range(len(gdf_year2))
]


#### Verify that the GeoDataFrames contain colza crops ####

colza_in_gdf = not gdf[gdf['CodeAgg'] == 'Col'].empty
colza_in_gdf_year2 = not gdf_year2[
    gdf_year2['CodeAgg'] == 'Col'
].empty


if colza_in_gdf:

    #### Create kernel, convolve and resize risk map ####
    
    out_grid = out_raster.geocube
    resized_area_length = max(out_grid.sizes['x'], out_grid.sizes['y'])
    x = y = np.linspace(
        -2*resized_area_length,
        2*resized_area_length,
        num = 2*resized_area_length + 1
    )
    z = utils.dispersion_kernel(x, y, mu = mu, beta2 = beta2)

    def apply_convolution(array):
        return scipy.signal.fftconvolve(array, z, mode = 'same')

    convolved_raster = out_grid.apply(apply_convolution)

    array = np.array(convolved_raster.NumValue)

    _="""
    To use the raster in zonal_stats() we need to extract the numpy array and
    the affine transform mapping the array dimensions to the crs:
    """
    raster_numvalue = convolved_raster.NumValue
    raster_array = raster_numvalue.values
    affine_transform = raster_numvalue.rio.transform()
    crs = raster_numvalue.rio.crs
    
    my_bar = st.progress(0)
    random_loading_text = np.random.choice(
        np.asarray(utils.loading_messages)
    )
    stats = []
    for i in range(len(gdf_year2) - 1):

        element = gdf_year2[
            gdf_year2['Domain_polygon'] != 0
        ].iloc[i:i+1]

        stats.append(
            zonal_stats(
                element, 
                raster_array,
                affine = affine_transform,
                crs = crs,
                stats = ['sum', 'mean', 'count', 'std'],
                all_touched = all_touched
            )[0]
        )

        my_bar.progress(
            int(i*100/(len(gdf_year2) - 1)), 
            text = random_loading_text)
        

    my_bar.empty()

    # Rescale mean and std:
    stats = utils.rescale_stats(stats, quality_reduction) 

    risk_list = [dico['mean'] for dico in stats]

else:
    risk_list = [0] * (len(gdf_year2) - 1)


#### Add info in gdf_year2 and plot ####

risk_list.insert(0, None)  # Add row for 1st polygon (region of interest) 

gdf_year2['risk'] = risk_list 

gdf_year2['Risk of infection:'] = [
    f'{i:.3g}' if i is not None 
    else None for i in risk_list
]


# Create a transformed version for better color visualization
if colza_in_gdf:
    risk_boxcox = scipy.stats.boxcox(risk_list[1:])[0].tolist()
    risk_boxcox.insert(0, None)
    gdf_year2['risk_values_for_visualization'] = risk_boxcox
else:
    gdf_year2['risk_values_for_visualization'] = risk_list

with map_tab:
    map_column, _, colorbar_column = st.columns([.6,.01,.25])

    with map_column:
        m = gdf_year2[gdf_year2['Domain_polygon'] != 0].explore(
            #column = 'risk',
            column = 'risk_values_for_visualization',
            popup = ['ID_PARCEL', 'Risk of infection:'],
            tooltip = 'Risk of infection:',  
            name = 'Risks for this year'
        )

        if not gdf[gdf['CodeAgg'] == 'Col'].empty:
            gdf[gdf['CodeAgg'] == 'Col'].explore(
                m = m,
                color = 'red',
                name = 'Where colza was last year',
                popup = ['ID_PARCEL'],
                tooltip = 'ID_PARCEL'
            )
        else:
            st.write(f'No colza crops were found at year {selected_year - 1} !')

        if not predictive:
            
            if not gdf_year2[gdf_year2['CodeAgg'] == 'Col'].empty:
                gdf_year2[gdf_year2['CodeAgg'] == 'Col'].explore(
                    m = m,
                    color = 'brown',
                    name = 'Colza fields this year',
                    #popup = ['ID_PARCEL'],
                    tooltip = ['ID_PARCEL', 'Risk of infection:'],
                    show = False  # Hidden by default
                )

            # Plot contact zones btwn colza of year 1 and colza of year 2
            boundary_gdf = utils.contact_zones_highlighter(
                gdf[gdf['CodeAgg'] == 'Col'],
                gdf_year2[gdf_year2['CodeAgg'] == 'Col'],
                buffer_distance = 5,  # Arbitrary
                crs = gdf.crs
            )
            boundary_gdf.explore(
                m = m,
                color = 'orange',
                name = 'Contact zones',
                show = False,
                tooltip = False,
                popup = False
            )

        if not colza_in_gdf_year2:
            st.write(f'No colza crops were found at year {selected_year} !')

        folium.LayerControl().add_to(m)

        # Change display with css input, to change font size on map:
        custom_css = """       
        <style>
        .leaflet-tooltip {
            font-size: 12pt !important;
        }
        .leaflet-popup-content {
            font-size: 12pt !important;
        }
        .leaflet-control-container  .legend {
            display: none !important;
        }
        .leaflet-control-container .leaflet-control-scale {
            font-size: 14pt !important;
        }
        </style>
        """
        # Get the map HTML
        map_html = m.get_root().render()
        
        # Inject the custom CSS into the map HTML
        map_html = map_html.replace('</head>', f'{custom_css}</head>')
        
        # Render the modified map in Streamlit
        st.components.v1.html(
            map_html, 
            height = 600
        )


    with colorbar_column:
        
        if colza_in_gdf:

            # Set figure and font size + normalize color values
            px = 1/plt.rcParams['figure.dpi']
            plt.rcParams.update({'font.size':6})
            fig, ax = plt.subplots(figsize = (20*px, 200*px))
            
            #max_value = max(risk_list[1:])
            max_value = max(gdf_year2['risk'][1:])
            min_value = min(gdf_year2['risk'][1:])
            norm = mpl.colors.Normalize(vmin = min_value, vmax = max_value)

            # Create colorbar and labels
            cbar = ax.figure.colorbar(          
                mpl.cm.ScalarMappable(norm=norm, cmap = 'viridis'),
                ax = ax,
                pad = .05,
                fraction = 1,
                #ticks = [min_value, (min_value + max_value)/2, max_value]
                ticks = [min_value, max_value]
            )

            cbar.ax.set_yticklabels([
                #f"{min_value:.2g} - Let's hit that dirt and plant some seeds !",
                f"{min_value:.2g} - Low risk: colza is happy and healthy !",
                #f"{(min_value - max_value)/2:.2g} - You're taking a risky path, partner !",
                f"{max_value:.2g} - Planting here ain't a good idea at all !"
            ])
            #cbar.ax.tick_params(labelsize = 10)
            
            ax.axis('off')
            st.pyplot(fig, use_container_width = True)

            st.write('Note: \nAs the risk repartition is generally highly \
            skewed towards 0, the color palette was adjusted using a Box-Cox\
            transform to enhance visualization. Thus the color in the middle\
            of the colorbar does not mean a "middle" risk.')


with metrics_tab:

    if not predictive:
        # Create table data
        metrics_dict = {
            f'Year {selected_year - 1}': utils.calculate_metrics(
                gdf, 2
            ),
            f'Year {selected_year}': utils.calculate_metrics(
                gdf_year2, 2
            )
        }

        metrics_df = pd.DataFrame(
            data = metrics_dict,
        )
        metrics_df.index = [
            "Number of colza fields",
            "Proportion of fields that contain colza",
            "Total surface of colza (ha)",
            "Surface of annual crops (ha)",
            "Surface of perennial crops (ha)",
            "Surface of unknown (ha)"
        ]

        # Mean closest distance between polygons in year n+1 and polygons in 
        # year n: 
        if colza_in_gdf and colza_in_gdf_year2:
            closest_distances = utils.closest_distances(
                gdf[gdf['CodeAgg'] == 'Col'],
                gdf_year2[gdf_year2['CodeAgg'] == 'Col']
            )
            mean_closest_distances = np.mean(closest_distances)
        else:
            closest_distances = np.array([])
        

        # Display
        st.write('#')
        st.write('### Recap metrics')
        _, main_column, _ = st.columns([.1, 1, .5])

        with main_column:
            st.table(metrics_df)

            boundary_area = round(
                boundary_gdf.area.sum() / 10000, # m² converted to ha
                2
            )
            if colza_in_gdf and colza_in_gdf_year2:
                st.write(
                    'Sum of areas of contact zones (see on map):' +
                    '$~~~~~~~$' +
                    str(boundary_area) + 
                    '$~$ha'
                )
        
        if closest_distances.size > 1:
            st.write(
                f'### Proximity of polygon centroids for year {selected_year} \
                with polygon centroids for year {selected_year - 1}'
            )

            _, main_column, _ = st.columns([.1, 1, .5])
            with main_column:
                fig, ax = plt.subplots(
                    2, 1, gridspec_kw = {'height_ratios': [9, .7]},
                )
                ax[0].hist(closest_distances, color = '#78ABA8')
                ax[0].spines['top'].set_visible(False)
                ax[0].spines['right'].set_visible(False)
                ax[0].set_ylabel('Number of polygons')
                
                _, xmax = ax[0].get_xlim()
                _, ymax = ax[0].get_ylim()
                ax[0].text(
                    xmax*.75,
                    ymax*.9,
                    f'Mean value: {round(np.mean(closest_distances), 2)}',
                    fontsize = 'large',
                    bbox = dict(facecolor = 'white', alpha = .5)
                )

                ax[1].boxplot(closest_distances, vert = False, showfliers = False)
                ax[1].spines['left'].set_visible(False)
                ax[1].spines['top'].set_visible(False)
                ax[1].spines['right'].set_visible(False)
                ax[1].set_yticklabels([])
                ax[1].set_yticks([])
                ax[1].set_xlim(ax[0].get_xlim())
                ax[1].set_xlabel('Distance with closest polygon centroid (m)')
                st.pyplot(fig)

        elif closest_distances.size == 1:
            with main_column:
                st.write(f'Distance between the centroid of the colza crop at\
                year {selected_year} and the centroid of the closest colza\
                crop at year {selected_year - 1}: $~~~~~~~$\
                {round(closest_distances[0], 2)} m')

    
    # Semi-variogram computation

    if colza_in_gdf:
        from scipy.ndimage import zoom
        convolved_array = np.array(convolved_raster.NumValue)
        area_size = gdf.total_bounds[2] - gdf.total_bounds[0]
        zoom_factor = area_size / convolved_array.shape
        convolved_array_resized = zoom(
            convolved_array, zoom_factor, order = 0, mode = 'nearest'
        )
        bin_centers, semivariance_values = utils.semivariogram(
            convolved_array_resized, 
            n_samples = 1000,
            max_distance = len(convolved_array_resized) / 2
        )

        # Display
        st.write(f'### Semi-variogram of the risk map for year {selected_year}')

        _, semi_variogram_column, _ = st.columns([.2, 1, .6])
        with semi_variogram_column:

            fig, ax = plt.subplots()
            ax.scatter(
                bin_centers, 
                semivariance_values, 
                s = 2, 
                color = '#78ABA8'
            )
            ax.set_xlabel('Distance (pixels)')
            ax.set_ylabel('Semivariance')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            st.pyplot(fig)


        # Distribution of total risks and areas after random sampling

        num_iter = 200
        num_samples = (gdf_year2['CodeAgg'] == 'Col').sum()

        risks, areas = utils.random_colza_selection_risk_distribution(
            gdf_year2[(gdf_year2['CodeAgg'] == 'Col') 
                             | (gdf_year2['CodeAgg'] == 'Ann')],
            risk_map = convolved_array,
            crs = crs,
            transform = affine_transform,
            num_samples = num_samples,
            num_iter = num_iter
        )

        if not predictive:

            obs_risk = zonal_stats(
                gdf_year2.loc[gdf_year2['CodeAgg'] == 'Col'],
                raster_array,
                affine = affine_transform,
                crs = crs,
                stats = ['sum']
            ) 
            obs_risk = sum([dico['sum'] for dico in obs_risk])
            ax1 = risks.axes[0]
            ax1.axvline(
                obs_risk, 
                color = '#EF9C66', 
                linestyle = 'dashed', 
                linewidth = 1
            )
            
            obs_area = gdf_year2.loc[
                gdf_year2['CodeAgg'] == 'Col', 'SURF_PARC'
            ].sum()
            ax2 = areas.axes[0]
            ax2.axvline(
                obs_area, 
                color = '#EF9C66', 
                linestyle = 'dashed', 
                linewidth = 1
            )

        # Display
        st.write(f'### Distribution of total risks and areas for year \
        {selected_year} after sampling {num_samples} colza crops randomly\
         {num_iter} times')
        _, risks_col, areas_col, _ = st.columns([.1, 1, 1, .3])
        with risks_col:
            st.pyplot(risks, use_container_width = True)
        with areas_col:
            st.pyplot(areas, use_container_width = True)
        if not predictive:
            st.write('The (beautiful) vertical orange bar shows the observed values with\
            available data for this year.')

          
        
