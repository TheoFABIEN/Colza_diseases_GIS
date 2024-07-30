import numpy as np 
from geocube.api.core import make_geocube 

class CustomRaster():
    """
    A class to rasterize and extract a simple numpy array from it.

    Args:
        shapefile (gpd.GeoDataFrame): GeoDataFrame containing vector data.
        type_of_plant (string) : Type of plant we want to focus on when we 
        rasterize. Default is None.

    Attributes:
        geocube (xarray.Dataset) : Object containing raster data.
        resolution (tuple) : spatial resolution of the raster (includes the direction
        as indicated by a negative value).
    """
    def __init__(self, shapefile, emission = 1, type_of_plant = None, resolution = None):
        
        # Adjust emission for the chosen resolution:
        adj_emission = emission*abs(resolution[0]**2)
        
        if type_of_plant:
            shapefile['NumValue'] = np.where(
                shapefile['CodeAgg'] == type_of_plant, adj_emission, 0
            )
        else:
            shapefile['NumValue'] = adj_emission
            
        
        if not resolution:
            #Add condition to reduce quality if shapefile is too big:
            if abs(shapefile.total_bounds[0] - shapefile.total_bounds[2]) > 2000:
                resolution = (-5, 5)
            else:
                resolution = (-1, 1)

        self.geocube = make_geocube(
            vector_data = shapefile,
            measurements = ['NumValue'],
            resolution = resolution,
            fill = 0
            )
        
        self.resolution = resolution
        self.emission_value = adj_emission
        

    def as_np_array(self, resize = 0):
        """
        Converts the xarray raster to a numpy array.

        Args:
            resize (int) : resize the square array to the specified value.
            emission (int) : value of emission per m².

        Returns:
            numpy.ndarray : Array representing the raster.
        """
        xar_grid = self.geocube.NumValue
        if resize != 0:
        	resized = xar_grid.rio.reproject(
            	xar_grid.rio.crs,
            	shape = (resize, resize)
        	)
        else :
            resized = xar_grid

        np_raster = resized.to_numpy()
        """
        Polygon values need to be multiplied by square of resolution.
        For example, if resolution = (-2, 2), one pixel is 2m², so one pixel 
        carries values for 2² = 4 pixels of the default map. 
        """
        poly_values = self.emission_value
        print(poly_values)
        np_raster[np_raster == 1] = poly_values

        return np_raster

        
