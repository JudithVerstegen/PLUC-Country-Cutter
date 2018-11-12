Readme Country cutter 1.0.0

Requirements:
   Software:
      Python 3.6
      gdal 2.2.2
      PCRaster 4.2.0
   Datasets:
      GADM 3.4 (https://gadm.org/data.html)
      Global datasets to cut (for PLUC dataset information see "PLUC_datasets_info.docx")

A tool for preparing PCRaster maps of an area of interest (e.g. province/country), using data from a collection of global maps.
The tool can cut out an area of interest, based on administrative areas, provided by the GADM dataset.
Input datasets must be rasters including projection information (e.g. .tif).
Support for shapes as input is not generic and is coupled with the processing of the GRIP dataset (roads map).
The script is devided into different sections, each dealing with a specific target.

The sections in which the user can give input are:
#---------------
# Script Set Up
#---------------
   In this section the user can define the name of the area of interest, the location of the GADM dataset,
   the GADM level of administrative area, the target resolution and the target projection.
   Parameters:
      - GADM_level is the administrative level of the input GADM dataset (e.g. Country level = 0)
      - CASE_STUDY_NAME should be the name of the administrative area (e.g. "Spain" in the country level)
      - target_projection should be in the form of a Proj4 string. The target projection can be customized
        with the function define_new_projection, by supplying a new origin, for example, the centerpoint 
        of the area of interest. (coordinate tuple given by the function get_map_middle_coordinates)
      - The target resolution should be given in the unit given in the target_projection proj4 string 

#-----------
# Datasets
#-----------
   Main part of the script. The user can define which data to cut and take other actions, such as adding a legend, 
   or performing some post-processing defined by the user in the map finalizing functions.
   Prefilled are the datasets used in PLUC (see "PLUC_datasets_info.docx")
   Parameters:
      dataset_name is the name of the dataset. this is used for the output filenames.
      dataset_filename should indicate the location of this dataset on the disk
      cut_dataset_filename is the return from the function cut_raster_dataset_to_area() and points to the file created
      cut_dataset_PCR_filename is the output of the convert_file_to_PCRmap() function and points to the file created
   Other usable functions:
      create_pcr_clone_and_nullmask() - creates clone.map and nullmask.map based on cut map
      PCR_add_legend() - adds a legend to a PCR map
      pcr_cover_nodata() - covers the no_data on a PCR map
      
#--------------------------
# Map finalizing functions
#--------------------------
   In this section the user can fine tune the output, or define other functions to adapt the output maps / create new maps based on the output maps 
   The following functions are used in the generating of the PLUC maps:
   - create_pcr_clone_and_nullmask() - creates clone.map and nullmask.map based on cut map, which are needed for the other functions
   - create_city_map() - creates a city map based on the land cover map
   - pcr_cover_nodata() - covers the nodata of a PCRaster map with the value 0
   - finalize_prot_area_maps() - performs a few actions on the protected areas map in order to get it ready as input for PLUC


System requirements:
   Minimal requirements unknown. Run successfully on 3.4GHz, 8GB RAM, 10GB free disk space.







