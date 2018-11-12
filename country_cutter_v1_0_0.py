from osgeo import ogr, gdal, osr
import os
import subprocess
import pcraster as pcr

""" 
A tool for creating region maps from global datasets
Capable of clipping maps to different levels of adminitrative areas
Can reproject, merge tiles, manual clipping, rasterizing shapes
Works with raster datasets and shape datasets.


Created with Python 3.6, gdal 2.2.2, PCRaster 4.2.0

"""


#-----------------------------------------------------------------------------------------------------------------------------------
#Script Set Up (inputs):
#-----------------------------------------------------------------------------------------------------------------------------------
#Case study name, for GADM area selection and for name giving of files (e.g. 'Ethiopia', 'Colombia', 'Poland', 'Thailand')
CASE_STUDY_NAME = 'Spain'
#location of the administrative areas dataset (GADM)
GADM_dataset_location = os.path.join(r'\\ad.geo.uu.nl','FG', 'PLUC/PLUC-OW/datasets/Countries/GADM/gadm34_levels_shp/gadm34_0.shp')
#level of administrative area (e.g. Country level = 0)
GADM_level = 0
#output projection. Origin can be customized with function define_new_projection()
# e.g. target_projection = "+proj=laea +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +units=m +no_defs"
target_projection = "+proj=laea +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +units=m +no_defs"
#resolution of the end product. Must be in same unit as defined in target_projection
target_resolution = 1000.0
#-----------------------------------------------------------------------------------------------------------------------------------


def runcommand(cmd):    
    """run command and show the return message
    """
    print('running command: ' + cmd)
    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               shell=True
                               )
    output, errors = process.communicate()
    if errors == '':
        print(output)
    else:
        print("Error : command '"+ str(cmd)+"'\n"+ str(errors))
    return

def convert_file_to_PCRmap(filename, datatype=""):
    """
    Convert from .tif to .map
    optional: set datatype, if gdal-pcraster does not recognise .tif data type
    datatype = "scalar" or "nominal"
    """
    if datatype!="":
        datatype={"scalar":"-ot Float32 ","nominal":"-ot Int32 "}[datatype]
    outputfilename = filename.split(".tif")[0] + ".map"
    cmd = "gdal_translate -of PCRaster " + datatype + filename + " " + outputfilename
    runcommand(cmd)
    print(filename + " converted to "+ outputfilename)
    return outputfilename

def create_shapefile_from_GADM(GADM_dataset_location, GADM_input_level):
    """Creates a shapefile of an area from the GADM dataset. This can be e.g. a country or a province.
        Requires the GADM dataset of the administrative level you require (level 0 = country level) 
        Parameters:
            GADM_dataset_location: Location of the GADM dataset
            GADM_input_level: Administative level
    """
    #Open GADM data and extract area shape
    GADM_dataset = ogr.Open(GADM_dataset_location)
    shape_data = GADM_dataset.GetLayer(0)
    area_shape = None
    for feature in shape_data:
        if feature['NAME_{}'.format(GADM_input_level)] == CASE_STUDY_NAME:
            area_shape = feature.geometry()
            break
    if GADM_dataset == None:
        raise RuntimeError("Cannot find the GADM dataset.")
    if area_shape == None:
        raise RuntimeError("Cannot find area. Check the area name. Example:'Thailand'")
    output_filename = os.path.join(os_scratch_folder, CASE_STUDY_NAME + "_wgs.shp")
    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    Datasource = shpdriver.CreateDataSource(output_filename)
    SpatialRef = shape_data.GetSpatialRef()
    outLayer = Datasource.CreateLayer(output_filename, srs = SpatialRef, geom_type=ogr.wkbPolygon)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(area_shape)
    outLayer.CreateFeature(outFeature)
    print("Created shapefile {}".format(output_filename))
    return output_filename

def get_bounds_from_shape(shape_filename):
    shapefile = ogr.Open(shape_filename)
    shape_data = shapefile.GetLayer(0)
    area_shape = shape_data.GetFeature(0)  
    area_shape_geometry = area_shape.geometry() 
    area_shape_geometry.ExportToJson()
    bounds = area_shape_geometry.GetEnvelope()
    return bounds

def get_map_middle_coordinates(shape_filename):
    """
    Get coordinates of the middle of the map for the custom projection.
    """
    bounds = get_bounds_from_shape(shape_filename)
    middlePointLon = (abs(bounds[0]) - abs(bounds[1])) / 2 + bounds[1]
    middlePointLat = (abs(bounds[3]) - abs(bounds[2])) / 2 + bounds[3]
    return middlePointLon, middlePointLat

def crop_shapefile(shape_filename, xmin, xmax, ymin, ymax):
    """Perform additional cropping on the shapefile to e.g. remove islands
    """
    output_shape_filename = os.path.join(os_scratch_folder, "{}_wgs_cropped.shp".format(CASE_STUDY_NAME))
    runcommand("ogr2ogr -clipsrc {} {} {} {} {} {}".format(xmin, ymin, xmax, ymax, output_shape_filename, shape_filename))
    print("Performing additional cropping..")
    return output_shape_filename

def extract_projection_from_shape(shape_filename):
    shapefile = ogr.Open(shape_filename)
    shape_data = shapefile.GetLayer(0)
    SpatialRef = shape_data.GetSpatialRef()
    proj4string = SpatialRef.ExportToProj4()
    return proj4string

def define_new_projection(input_projection_proj4, centerpoint_x, centerpoint_y):
    """
    Defines a new projection based on an unspecific input projection and a specified centerpoint
    """
    proj4str_parts = input_projection_proj4.split(" ")
    for i, str_part in enumerate(proj4str_parts):
        if proj4str_parts[i].startswith("+lat"):
            proj4str_parts[i] = "+lat_0={}".format(centerpoint_y)
        if proj4str_parts[i].startswith("+lon"):
            proj4str_parts[i] = "+lon_0={}".format(centerpoint_x)
    proj4str_replaced = " ".join(proj4str_parts)
    print("Defined custom projection: {}".format(proj4str_replaced))    
    return proj4str_replaced  

def reproject_shapefile(input_shape_filename, new_projection_proj4):
    """
    Reproject shapefile to specified projection defined in a proj4 string.
    Defaults to config variable target_projection
    """
    output_shape_filename= os.path.join(os_scratch_folder, CASE_STUDY_NAME + "_laea.shp")
    #Open input shapefile and get geometry and projection
    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(output_shape_filename):
        shpdriver.DeleteDataSource(output_shape_filename)
    input_shapefile = ogr.Open(input_shape_filename)
    shape_data = input_shapefile.GetLayer(0)
    area_shape = shape_data.GetFeature(0) 
    area_shape_geometry = area_shape.geometry()
    inSpatialRef = shape_data.GetSpatialRef()
    #Apply new projection
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromProj4(new_projection_proj4)
    coordinate_system = osr.CoordinateTransformation(inSpatialRef, outSpatialRef) 
    area_shape_geometry.Transform(coordinate_system)
    #Save new shapefile
    Datasource = shpdriver.CreateDataSource(output_shape_filename)
    outLayer = Datasource.CreateLayer(output_shape_filename, srs = outSpatialRef, geom_type=ogr.wkbPolygon)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(area_shape_geometry)
    outLayer.CreateFeature(outFeature)
    return output_shape_filename

def cut_raster_dataset_to_area(dataset_name, dataset_filename, shape_filename, resolution, target_projection, resampling_method = "near", nodata_value = -32768):
    """Cuts a raster dataset with a shapefile
    Parameters:
        dataset_name: Name of (global) raster dataset to cut
        dataset_filename: Filename of (global) raster dataset to cut
        shape_filename: Filename of shapefile to cut with
        resolution: Target resolution
        target_projection: Projection to apply, defined in a proj4 string
        resampling_method: Name of resampling method. Choose between: near, bilinear, cubic, cubicspline, lanczos, average, mode, max, min, med, q1, q3
        nodata_value: Value that represents missing data
    """
    if extract_projection_from_shape(shape_filename).split("+")[1].split("=")[1] != target_projection.split("+")[1].split("=")[1]:
        raise RuntimeError("Input shape of {} is not in the same projection as the target! First use reproject_shapefile.".format(dataset_name))
    output_raster_filename = os.path.join(os_scratch_folder, dataset_name + '_laea.tif')
    bounds = get_bounds_from_shape(shape_filename)
    bounds_gdalwarp_format = ("{} "*4).format(bounds[0], bounds[2], bounds[1], bounds[3])
    runcommand('gdalwarp -overwrite -te {} -tr {} {} -t_srs "{}" -r {} {} {}'.format(bounds_gdalwarp_format, resolution, resolution, target_projection, resampling_method, dataset_filename, output_raster_filename))
    reprojected_dataset_filename = output_raster_filename
    output_raster_filename = os.path.join(os_output_folder, "{}{}{}".format(dataset_name, CASE_STUDY_NAME, '_laea.tif'))
    runcommand('gdalwarp -overwrite -te {} -tr {} {} -cutline {} -crop_to_cutline -dstnodata "{}" {} {}'.format(bounds_gdalwarp_format, resolution, resolution, shape_filename, nodata_value, reprojected_dataset_filename, output_raster_filename))
    return output_raster_filename

def PCR_add_legend(PCR_map_filename, legend_filename):
    runcommand('legend --clone ' + PCR_map_filename + ' -f ' + legend_filename + ' ' + PCR_map_filename)
    return

def read_and_rasterize_GRIP(GRIP_dataset_filename, area_shape_WGS_filename, reference_image):
    """"Reads the GRIP dataset and outputs rasters of all road types
    Parameters:
        GRIP_dataset_filename: location of the GRIP geodatabase
        area_shape_WGS_filename: location of area shape in WGS84 projection
        reference_image: Image that is taken as reference to position the raster"""
    #open geodatabase
    driver = ogr.GetDriverByName("OpenFileGDB")
    ds = driver.Open(GRIP_dataset_filename)
    layer = ds.GetLayer(0)
    #open shapefile and get geometry
    area_shape = ogr.Open(area_shape_WGS_filename)
    layer2 = area_shape.GetLayer(0)
    feature = layer2.GetFeature(0)
    geometry = feature.GetGeometryRef()
    #Create a class map for each road type
    GRIP_road_types = {1:"Highways", 2:"Primary_roads", 3:"Secondary_roads", 4:"Tertiary_roads", 5:"Local_roads"}
    for roadtype in range(1,6):
        #filter extent
        selection = "{}='{}'".format("gp_rtp", roadtype)
        layer.SetAttributeFilter(selection)
        layer.SetSpatialFilter(geometry)
        #create new shapefile with altered extent
        output_filename= os.path.join(os_scratch_folder, "roads_" + GRIP_road_types[roadtype] + "_" + CASE_STUDY_NAME +".shp")
        shpdriver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(output_filename):
            shpdriver.DeleteDataSource(output_filename)
        Datasource = shpdriver.CreateDataSource(output_filename)
        print("copying layer..")
        outLayer=Datasource.CopyLayer(layer, "roads_" + GRIP_road_types[roadtype])
        del Datasource, outLayer
        
        #Rasterize shape
        InputVector = os.path.join(os_scratch_folder,"roads_" + GRIP_road_types[roadtype] + "_" + CASE_STUDY_NAME +".shp")
        OutputImage = os.path.join(os_scratch_folder, "roads_" + GRIP_road_types[roadtype] + "_" + CASE_STUDY_NAME +".tif")
    
        gdalformat = 'GTiff'
        datatype = gdal.GDT_Byte
        burnVal = roadtype #value for the output image pixels
    
        # Get projection info from reference image
        Image = gdal.Open(reference_image, gdal.GA_ReadOnly)

        # Open Shapefile
        Shapefile = ogr.Open(InputVector)
        Shapefile_layer = Shapefile.GetLayer()

    
        # Rasterise
        print("Rasterising shapefile...")
        Output = gdal.GetDriverByName(gdalformat).Create(OutputImage, Image.RasterXSize, Image.RasterYSize, 1, datatype, options=['COMPRESS=DEFLATE'])
        Output.SetProjection(Image.GetProjectionRef())
        Output.SetGeoTransform(Image.GetGeoTransform()) 
    
        # Write data to band 1
        Band = Output.GetRasterBand(1)
        Band.SetNoDataValue(0)
        gdal.RasterizeLayer(Output, [1], Shapefile_layer, burn_values=[burnVal])
    
        # Close datasets
        Band = None
        Output = None
        Image = None
        Shapefile = None
    
        # Build image overviews
        subprocess.call("gdaladdo --config COMPRESS_OVERVIEW DEFLATE "+OutputImage+" 2 4 8 16 32 64", shell=True)
        print("Done.")

#-----------------------------------------------------------------------------------------
# Map finalizing functions. user defined functions to adapt the output to the exact needs
#-----------------------------------------------------------------------------------------

def create_pcr_clone_and_nullmask(land_use_map_filename):
    land_use_map = pcr.readmap(land_use_map_filename)
    clone_map = pcr.boolean(1)
    pcr.report(clone_map, os.path.join(os_output_folder, "clone_" + CASE_STUDY_NAME + ".map"))
    nullmask_map = pcr.ifthen(pcr.boolean(land_use_map), pcr.boolean(0))
    pcr.report(nullmask_map, os.path.join(os_output_folder, "nullMask_" + CASE_STUDY_NAME + ".map"))
    return clone_map, nullmask_map

def create_city_map(ESA_CCI_land_use_map_filename, city_class_id):
    #Create cities map from land use map
    land_use_map = pcr.readmap(ESA_CCI_land_use_map_filename)
    cities = pcr.ifthenelse(land_use_map==city_class_id, pcr.boolean(1), pcr.boolean(0))
    pcr.report(cities, os.path.join(os_output_folder, "cities_" + CASE_STUDY_NAME + ".map"))
    return

def pcr_cover_nodata(map_filename, nullmask_map):
    #Cover nodata
    input_map = pcr.readmap(map_filename)
    map_no_data_covered = pcr.cover(input_map, pcr.scalar(nullmask_map))
    pcr.report(map_no_data_covered, map_filename)
    return

def finalize_prot_area_maps(protected_areas_map, nullmask_map):
    #Finalizing protected areas
    protected_areas_map_country = pcr.cover(protected_areas_map, 0)
    protected_areas_map_country_cut = pcr.ifthen(nullmask_map==0, protected_areas_map_country)
    protected_areas_map_country_bool = pcr.ifthenelse(pcr.scalar(protected_areas_map_country_cut)>0, pcr.boolean(1), pcr.boolean(0))
    pcr.report(protected_areas_map_country_bool, os.path.join(os_output_folder, "protected_areas_map_" + CASE_STUDY_NAME + "_bool.map"))


#Creating and cleaning directories
if not os.path.isdir("scratch"):
    os.mkdir("scratch")
if not os.path.isdir("output"):
    os.mkdir("output")
os_scratch_folder = os.path.join(os.getcwd(), "scratch")
os_output_folder = os.path.join(os.getcwd(), "output")
for file in os.listdir(os_scratch_folder):
    os.remove(os.path.join(os_scratch_folder, file))

#Create area shape from GADM
shape_country_filename = create_shapefile_from_GADM(GADM_dataset_location, GADM_level)
#Additional clipping: take e.g. islands off
if CASE_STUDY_NAME == 'Spain':
    shape_country_filename = crop_shapefile(shape_country_filename, -13, 180, -180, 180)
elif CASE_STUDY_NAME == 'Colombia':
    shape_country_filename = crop_shapefile(shape_country_filename, -80, 180, -180, 14)
#Reproject shapefile to target projection
centerpoint_x, centerpoint_y = get_map_middle_coordinates(shape_country_filename)
custom_projection = define_new_projection(target_projection, centerpoint_x, centerpoint_y) #combine with get map middle?
shapefile_country_laea = reproject_shapefile(shape_country_filename, custom_projection)


#----------------------------------------
# Datasets (PLUC)
#----------------------------------------
print("Cutting datasets..")

#Create land use map from ESA CCI 2015
dataset_name = "land_use_ESA_CCI_2015_"
dataset_filename = r"\\ad.geo.uu.nl\GIS\Data\world\ESA_LandCoverCCI\v2.07\ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif"
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection)
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "nominal")
legend_filename = r"C:\Users\3402142\Desktop\PLUC-OW\script\data-for-practical\legendLU.txt"
PCR_add_legend(cut_dataset_PCR_filename, legend_filename)

#create clone, nullmask and city maps
ESA_CCI_land_use_map_filename = os.path.join(os_output_folder, 'land_use_ESA_CCI_2015_' + CASE_STUDY_NAME + '_laea.map')
clone_map, nullmask_map = create_pcr_clone_and_nullmask(ESA_CCI_land_use_map_filename)
create_city_map(ESA_CCI_land_use_map_filename, 190)

#Prepare DEM from ETOPO1
dataset_name = "dem_"
dataset_filename = r"\\ad.geo.uu.nl\GIS\Data\ETOPO\ETOPO1\ETOPO1_Bed_g.tif"
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection)
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "scalar")

#Prepare FAO suitability maps (cereals, alfalfa) (bilinear interpolation)
dataset_name = "suitability_cereals_"
dataset_filename = r"\\ad.geo.uu.nl\FG\PLUC\PLUC-OW\datasets\Suitability\res03_crav6190h_sxhr_cer.tif"
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection, resampling_method="bilinear")
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "scalar")
dataset_name = "suitability_alfalfa_"
dataset_filename = r"\\ad.geo.uu.nl\FG\PLUC\PLUC-OW\datasets\Suitability\res03_crav6190l_sxlr_alf.tif"
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection, resampling_method="bilinear")
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "scalar")

#Prepare population density map
dataset_name = "pop_density_"
dataset_filename = r'\\ad.geo.uu.nl\FG\PLUC\PLUC-OW\datasets\pop-dens\GHS_POP_GPW42015_GLOBE_R2015A_54009_1k_v1_0_GTIFF.tif'
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection)
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "scalar")

#Prepare cropland/pasture maps from Raman Kutty
dataset_name = "cropland_area_Ramankutty_"
dataset_filename = r'\\ad.geo.uu.nl\FG\PLUC\PLUC-OW\datasets\Other\CroplandPastureArea2000_Geotiff\cropland2000_area.tif'
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection, resampling_method="bilinear")
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "scalar")
pcr_cover_nodata(cut_dataset_PCR_filename, nullmask_map)

dataset_name = "pasture_area_Ramankutty_"
dataset_filename = r'\\ad.geo.uu.nl\FG\PLUC\PLUC-OW\datasets\Other\CroplandPastureArea2000_Geotiff\pasture2000_area.tif'
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection, resampling_method="bilinear")
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "scalar")
pcr_cover_nodata(cut_dataset_PCR_filename, nullmask_map)

#Preparing protected areas
dataset_name = "protected_area_"
dataset_filename = r'\\ad.geo.uu.nl\FG\PLUC\PLUC-OW\datasets\NoGo\lrprtpa5min_package\lr_prt_pa_5min.tif'
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection)
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "nominal")
finalize_prot_area_maps(cut_dataset_PCR_filename, nullmask_map)

#Prepare road class map
reference_image = os.path.join(os_output_folder, "dem_{}_laea.tif".format(CASE_STUDY_NAME))        
read_and_rasterize_GRIP(r"\\ad.geo.uu.nl\GIS\Data\world\GRIP\GRIP4_global_vector_fgdb\GRIP4_GlobalRoads.gdb", shape_country_filename, reference_image)
GRIP_road_types = {1:"Highways", 2:"Primary_roads", 3:"Secondary_roads", 4:"Tertiary_roads", 5:"Local_roads"}   
filenames = [os.path.join(os_scratch_folder,"roads_" + GRIP_road_types[i] + "_" + CASE_STUDY_NAME +".tif") for i in range(1,6)]
#merge maps in 1
runcommand(r'python C:\Users\3402142\AppData\Local\Continuum\Anaconda2\envs\python3\Scripts\gdal_merge.py -o {} -ot Int32 -n 0 {} {} {} {} {}'.format(os.path.join(os_scratch_folder, 'merged_roadmap.tif'), filenames[4], filenames[3], filenames[2], filenames[1], filenames[0])) 
output_filename = os.path.join(os_scratch_folder, "road_class_map_" + CASE_STUDY_NAME + ".map")
runcommand("gdal_translate -OF PCRaster -ot Int32 " + os.path.join(os_scratch_folder, 'merged_roadmap.tif') + " " + output_filename)
dataset_name = "roads_class_map_"
dataset_filename = os.path.join(os_scratch_folder, "road_class_map_" + CASE_STUDY_NAME + ".map")
cut_dataset_filename = cut_raster_dataset_to_area(dataset_name, dataset_filename, shapefile_country_laea, target_resolution, custom_projection)
cut_dataset_PCR_filename = convert_file_to_PCRmap(cut_dataset_filename, "nominal")




print("All maps cut for " + CASE_STUDY_NAME + ". Output in: " + os_output_folder)
