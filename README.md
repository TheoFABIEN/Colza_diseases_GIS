Summary of the steps finalized in this project:
1) Extracting sub_regions from the RPG dataset and creating relevant categories for polygons (Col, Ann, Per).
2) Using a convolution operation on the raster representing an area at year n, build a heatmap representing the risk of infection for year n+1.
3) Perform sensitivity analysis of these prediction.
4) Extract various metrics (areas characteristics, contingency zones between year n and year n+1, model outputs...)

### DATASET
The dataset used in this project can be found [here](https://geoservices.ign.fr/rpg).

We gathered shapefiles ranging from 2015 to 2021, for the following regions: Bretagne, Centre Val-de-Loire, Franche-Comté, Grand Est, Hauts de France, Ile de France, Normandie, Nouvelle Aquitaine, Occitanie.

### STREAMLIT APP
A streamlit app was built to provide a CLI to interact with, when studying the impact of colza fields placement at year n on risks at year n+1.
In the project repository, use the following command to launch it:
```
python3 -m streamlit-run risk-maps.py
```
### FILE STRUCTURE
utils: contains all utilitary files, each containing custom functions and classes used in the project.

RPG_Sub_Region_Selector.py: See step 1) of the project.

risk_maps.py: see streamlit app.

metrics_between_regions.py: computes large-scale metrics and corresponding boxplots. This is used to compare tendencies between different studied regions.

Diffusion.py: Apply convolution to each shapefile of the dataset and compute all the metrics relative to the shapefile, its equivalent for the upcoming year, their contingency and the risks of infection. 
Can optionally perform sensitivity analysis and random shuffling of the colza fields at year n+1. At the end, everything is put in a big .csv file which can be found in the same location as the dataset. 

tests.ipynb: Notebook displaying basic operations and results used repeatedly in this project.

Conversion_Codes_Cultures.csv: data table used for creating the "CodeAgg" column in GeoDataFrames. The goal of this is to group the crops under different types (Colza (crops of interest), Annual crops, Perennial crops). You will see this table called in the RPG_Sub_Region_Selector file.

### REQUIREMENTS
Requirements are not known to be version-sensitive for now. They can be installed normally with pip or conda.

+ contextily
+ folium
+ geocube
+ geopandas
+ matplotlib
+ matplotlib_scalebar
+ numpy
+ pandas
+ multiprocessing
+ pyogrio
+ rasterstats
+ scipy
+ Shapely
+ scikit-image
+ streamlit
+ streamlit_folium
+ tqdm
+ xarray

