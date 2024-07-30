import geopandas as gpd
from rasterstats import zonal_stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist, pdist, squareform
from shapely.geometry import Polygon, mapping
from skimage.draw import polygon2mask


def calculate_metrics(gdf, round_val):
    """
    Calculates metrics related to field types in a shapefile.

    Parameters:
    gdf (GeoDataFrame): GeoDataFrame containing shapefile data.
    round_val (int): Number of decimal places for rounding metrics.

    Returns:
    dict: Dictionary of calculated metrics.
    """
    # Total surface area for each field type
    surf_colza = round(gdf[gdf['CodeAgg'] == 'Col']['SURF_PARC'].sum(), round_val)
    surf_ann = round(gdf[gdf['CodeAgg'] == 'Ann']['SURF_PARC'].sum(), round_val)
    surf_per = round(gdf[gdf['CodeAgg'] == 'Per']['SURF_PARC'].sum(), round_val)

    # Number of colza fields
    num_colza_fields = len(gdf[gdf['CodeAgg'] == 'Col'])
    
    # Proportion of colza fields relative to annual fields
    num_ann_fields = len(gdf[gdf['CodeAgg'] == 'Ann'])
    prop_colza_fields = round(
        num_colza_fields / num_ann_fields if num_ann_fields > 0 else 0, 3
    )

    # Total area of the shapefile in hectares
    total_area = gdf.iloc[0].geometry.area / 10000  # Convert m² to ha
    known_surface = surf_colza + surf_ann + surf_per
    surf_unknown = round(total_area - known_surface, round_val)
    
    return {
        "Number_colza_fields": str(num_colza_fields),
        "Proportion_of_fields_that_contain_colza": str(prop_colza_fields),
        "Total_surface_colza": str(surf_colza),
        "Surface_Ann": str(surf_ann),
        "Surface_Per": str(surf_per),
        "Surface_unknown": str(surf_unknown)
    }



def semivariogram(risk_map, n_samples = 500, max_distance = None):
    """
    Calculate the semivariogram of the given risk map.

    Parameters:
    risk_map (numpy array of shape (m, n)): The risk map.
    n_samples (int, optional): The number of random samples to extract from the 
    risk map. Default is 500.
    max_distance (float, optional): The maximum distance to consider when 
    calculating the semivariogram. Default is None.

    Returns:
    bin_centers (numpy array of shape (k,)): The centers of the distance bins.
    semivariance_values (numpy array of shape (k,)): The average semivariance 
    for each distance bin.
    """
    # Extract random coordinates and values
    coordinates = np.random.randint(0, len(risk_map), (n_samples, 2))
    values = np.asarray(
        [risk_map[x, y] for x, y in coordinates]
    )

    # Get the pairwise distances between coordinates
    distances = pdist(coordinates)
    semivariances = pdist(values.reshape(-1, 1), metric = 'sqeuclidean') / 2

    # Convert to a square form matrix
    distance_matrix = squareform(distances)
    semivariance_matrix = squareform(semivariances)


    # Get unique distance bins
    if max_distance is None:
        max_distance = np.max(distances)
    bins = np.linspace(0, max_distance, num = 50)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # Calculate average semivariance for each bin 
    semivariance_values = np.zeros_like(bin_centers)
    for i in range(len(bins) - 1):
        mask = (distances >= bins[i]) & (distances < bins[i + 1])
        semivariance_values[i] = np.mean(semivariances[mask])

    return bin_centers, semivariance_values



def closest_distances(geodata1, geodata2):
    """
    Calculate the Euclidean distance between each centroid in shapefile 2
    and the closest centroid from shapefile 1.
    
    Parameters:
    geodata1 (GeoDataFrame)
        The first shapefile (shapefile 1) containing polygons.
    geodata2 (GeoDataFrame)
        The second shapefile (shapefile 2) containing polygons.
        
    Returns:
    min_distances (numpy array)
        An array of the minimum Euclidean distances between each centroid in 
        shapefile 2 and the closest centroid in shapefile 1.
    """
    # Extract centroids
    centroids1 = np.array(
        [(geom.centroid.x, geom.centroid.y) 
         for geom in geodata1.geometry]
    )
    centroids2 = np.array(
        [(geom.centroid.x, geom.centroid.y) 
         for geom in geodata2.geometry]
    )
    
    # Calculate pairwise distances
    distances = cdist(centroids2, centroids1)
    
    # Find the minimum distance for each centroid in shapefile 2
    min_distances = np.min(distances, axis = 1)
    
    return min_distances
    
    
    
    
# Next 2 functions will be used for analysis of contact zones between places 
# where colza was at year n and places where colza was at year n+1
    
    
    
def lines_to_polygons(lines):
    """
    Converts closed LineStrings or MultiLineStrings to Polygon objects.
    This function will only be used in contact_zones_highlighter()

    Parameters:
    lines (list): A list of shapely LineString or MultiLineString objects.

    Returns:
    polygons (list): A list of shapely Polygon objects converted from closed 
    LineStrings or MultiLineStrings.
    """
    polygons = []
    for line in lines:
        if line.geom_type == "LineString" and line.is_closed:  # Check if the LineString is a closed loop
            polygons.append(Polygon(line))
        elif line.geom_type == 'MultiLineString':
            for geom in line.geoms:
                if geom.is_closed:
                    polygons.append(Polygon(geom))
    return polygons


def contact_zones_highlighter(shapefile_1, shapefile_2, 
                               buffer_distance = 1, crs = None):                                              
    """
    Finds and returns the shared boundaries between two shapefiles as a 
    GeoDataFrame of Polygons. It does so by buffering the geometries in the 
    input shapefiles by a given distance, and then find intersections. The 
    intersection are then converted to polygons.

    Parameters:
    shapefile_1 (GeoDataFrame): The first input GeoDataFrame.
    shapefile_2 (GeoDataFrame): The second input GeoDataFrame.
    buffer_distance (float, optional): The distance to buffer the geometries. 
                                       Default is 1.
    crs (str or int, optional): The Coordinate Reference System to use for the 
                                output GeoDataFrame. Default is None.

    Returns:
    polygons_gdf (GeoDataFrame): A GeoDataFrame containing the shared boundaries 
                                 as Polygon objects.
    """      
    shapefile_1 = shapefile_1.copy()
    shapefile_2 = shapefile_2.copy()                                                  
    shapefile_1['geometry'] = shapefile_1.buffer(buffer_distance)
    shapefile_2['geometry'] = shapefile_2.buffer(buffer_distance)
    
    shared_boundaries = []
    
    for poly1 in shapefile_1['geometry']:
        for poly2 in shapefile_2['geometry']:
            intersection = poly1.intersection(poly2)
            if not intersection.is_empty:
                # Extract the shared boundary (line segments)
                shared_boundary = intersection.boundary
                shared_boundaries.append(shared_boundary)
    
    # Create a GeoDataFrame to store the shared boundaries
    shared_gdf = gpd.GeoDataFrame(geometry=shared_boundaries)
    #shared_gdf['length'] = shared_gdf.length

    # Create a gdf containing polygons instead of linestrings
    lines = shared_gdf.geometry
    polygons = lines_to_polygons(lines)
    polygons_gdf = gpd.GeoDataFrame(geometry = polygons, crs = crs)

    return polygons_gdf



def contiguity_metrics(
    shapefile_1, shapefile_2, buffer_distance = 1, round_val = 2, crs = None):
    """
    Calculate various contiguity metrics between two shapefiles containing 
    polygons.

    Parameters:
    shapefile_1 (GeoDataFrame): First shapefile.
    shapefile_2 (GeoDataFrame): Second shapefile.
    buffer_distance (float, optional, default 1): Buffer distance used to 
        highlight contact zones between two shapefiles.
    crs (dict or str, optional): Coordinate reference system to be used for the
        contact zones. If None, the crs of the input shapefiles is used.

    Returns:
    metrics (dict): Dictionary containing the following metrics:
        - 'boundary_area' (float): The total area of contact zones between the 
            two shapefiles (ha).
        - 'jaccard_index' (float): Jaccard index of the two shapefiles.
        - 'mean_distance_centroids' (float): Mean Euclidean distance between 
            the centroids of the polygons in the two shapefiles.
        - 'std_distance_centroids' (float): Variance of the Euclidean 
            distances between the centroids of the polygons in the two 
            shapefiles.
        - 'mean_distance_edges' (float): Mean Euclidean distance between the 
            edges of the polygons in the two shapefiles.
        - 'std_distance_edges' (float): Variance of the Euclidean 
            distances between the edges of the polygons in the two shapefiles.
        - 'mean_minimal_distances' (float): Mean of the minimal Euclidean 
            distances from each centroid in shapefile_2 to the closest centroid 
            in shapefile_1.
        - 'std_minimal_distances' (float): Variance of the minimal 
            Euclidean distances from each centroid in shapefile_2 to the 
            closest centroid in shapefile_1.
            
    """
    if len(shapefile_1) != 0 and len(shapefile_2) != 0:

        # Area of contact zones
        boundary_gdf = contact_zones_highlighter(
            shapefile_1,
            shapefile_2, 
            buffer_distance = buffer_distance,
            crs = crs
        )
        boundary_area = boundary_gdf.area.sum() / 10000 # m² converted to ha
        boundary_area = round(boundary_area, round_val)

        # Jaccard index
        union1 = shapefile_1.unary_union
        union2 = shapefile_2.unary_union
        intersection = union1.intersection(union2)
        union = union1.union(union2)
        jaccard = round(intersection.area / union.area, round_val)

        # Distance between centroids (mean and std)
        centroids1 = shapefile_1.geometry.centroid
        centroids2 = shapefile_2.geometry.centroid
        coords1 = np.array([[point.x, point.y] for point in centroids1])
        coords2 = np.array([[point.x, point.y] for point in centroids2])
        distances_centroids = cdist(coords1, coords2)
        mean_distance_centroids = round(np.mean(distances_centroids), round_val)
        std_distance_centroids = round(np.std(distances_centroids), round_val)

        # Distance between edges (mean and std)
        distances_edges = [
            poly1.distance(poly2) 
            for poly1 in shapefile_1.geometry
            for poly2 in shapefile_2.geometry
        ]
        mean_distance_edges = round(np.mean(distances_edges), round_val)
        std_distance_edges = round(np.std(distances_edges), round_val)

        # Closest distance n n+1 (mean and std)
        minimal_distances = closest_distances(shapefile_1, shapefile_2)
        mean_minimal_distances = round(np.mean(minimal_distances), round_val)
        std_minimal_distances = round(np.std(minimal_distances), round_val)

    else:
        boundary_area = 0
        jaccard = None
        mean_distance_centroids, std_distance_centroids = [None]*2
        mean_distance_edges, std_distance_edges = [None]*2
        mean_minimal_distances, std_minimal_distances = [None]*2


    metrics = {
        'boundary_area': boundary_area,
        'jaccard_index': jaccard,
        'mean_distance_centroids': mean_distance_centroids,
        'std_distance_centroids': std_distance_centroids, 
        'mean_distance_edges': mean_distance_edges,
        'std_distance_edges': std_distance_edges,
        'mean_minimal_distances': mean_minimal_distances,
        'std_minimal_distances': std_minimal_distances
    }

    return metrics




def random_colza_selection_risk_distribution(
    gdf, risk_map, transform, crs, num_samples, num_iter):
    """
    Gives distribution of risks of infection by over random selections of a 
    given number of polygons that may contain colza. 

    Parameters:
    gdf (geopandas.GeoDataFrame): GeoDataFrame containing polygons representing 
                                  colza crops.
    risk_map (numpy.ndarray): Risks of infection in the study area, computed
                              from data of the previous year.
    transform (affine.Affine): Transform to map the risk map coordinates to the
                               CRS of the GeoDataFrame.
    crs (rasterio.crs.CRS): CRS of the GeoDataFrame.
    num_samples (int): Number of polygons to select randomly for each iteration.
    num_iter (int): Number of iterations to run the simulation.

    Returns:
    matplotlib.figure.Figure: Distribution of total risks after a large number
                               of random draws. 
    """
    if risk_map is None:
        raise ValueError('A risk map must be provided.')
    
    risk_stats = zonal_stats(
        gdf, risk_map, affine = transform, crs = crs, stats = ['sum']
    )
    gdf['risk'] = [list(i.values())[0] for i in risk_stats]
    
    total_risks, total_areas = [], []
    for _ in range(num_iter):
        sample = gdf.sample(num_samples, replace = False)
        total_risks.append(sample['risk'].sum())
        total_areas.append(sample['SURF_PARC'].sum())
        
    fig1, ax1 = plt.subplots(figsize = (4,3))
    ax1.hist(
        total_risks,
        bins = 30,
        color = '#78ABA8',
    )
    ax1.set_xlabel('Total Risk'),
    ax1.set_ylabel('Frequency')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    fig2, ax2 = plt.subplots(figsize = (4,3))
    ax2.hist(
        total_areas,
        bins = 30,
        color = '#78ABA8'
    )
    
    ax2.set_xlabel('Total area (ha)')
    ax1.set_ylabel('Frequency')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    fig1.tight_layout()
    fig2.tight_layout()
    
    return fig1, fig2
    
    
    
def rescale_stats(stats, rescale_factor):
    """
    Rescales the 'mean' and 'std' values in a list of dictionaries by the 
    square of the given rescale factor. This is to adjust the statistical values 
    to the resolution of the spatial data.
    
    Parameters:
    stats (list of dict): A list of dictionaries containing statistical values
                          including the 'mean' and 'std'.
    rescale_factor (float): The factor by which to rescale the values.
    
    Returns:
    list of dict: The input list of dictionaries with 'mean' and 'std' values 
                  rescaled.
    """
    for item in stats:
        item['mean'] /= rescale_factor**2
        item['std'] /= rescale_factor**2
    
    return stats



    
    
    
    
    


