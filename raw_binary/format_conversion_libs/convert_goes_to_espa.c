/*****************************************************************************
FILE: convert_goes_to_espa.c
  
PURPOSE: Contains functions for reading GOES-R ABI netCDF products and writing
to ESPA raw binary format.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
  1. The XML metadata format written via this library follows the ESPA internal
     metadata format found in ESPA Raw Binary Format vX.Y.doc.  The schema for
     the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/schema/espa_internal_metadata_vX_Y.xsd.
*****************************************************************************/
#include <unistd.h>
#include <ctype.h>
#include "convert_goes_to_espa.h"


/******************************************************************************
MODULE:  doy_to_month_day

PURPOSE: Convert the DOY to month and day.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error converting from DOY to month and day
SUCCESS         Successfully converted from the DOY to the month and day

NOTES:
******************************************************************************/
int doy_to_month_day
(
    int year,            /* I: year of the DOY to be converted */
    int doy,             /* I: DOY to be converted */
    int *month,          /* O: month of the DOY */
    int *day             /* O: day of the DOY */
)
{
    char FUNC_NAME[] = "doy_to_month_day";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    bool leap;                /* is this a leap year? */
    int i;                    /* looping variable */
    int nday_lp[12] = {31, 29, 31, 30,  31,  30,  31,  31,  30,  31,  30,  31};
        /* number of days in each month (for leap year) */
    int idoy_lp[12] = { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336};
        /* starting DOY for each month (for leap year) */
    int nday[12] = {31, 28, 31, 30,  31,  30,  31,  31,  30,  31,  30,  31};
        /* number of days in each month (with Feb being a leap year) */
    int idoy[12] = { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};
        /* starting DOY for each month */

    /* Is this a leap year? */
    leap = (bool) (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));

    /* Determine which month the DOY falls in */
    *month = 0;
    if (leap)
    {  /* leap year -- start with February */
        for (i = 1; i < 12; i++)
        {
            if (idoy_lp[i] > doy)
            {
                *month = i;
                *day = doy - idoy_lp[i-1] + 1;
                break;
            }
        }

        /* if the month isn't set, then it's a December scene */
        if (*month == 0)
        {
            *month = 12;
            *day = doy - idoy_lp[11] + 1;
        }
    }
    else
    {  /* non leap year -- start with February */
        for (i = 1; i < 12; i++)
        {
            if (idoy[i] > doy)
            {
                *month = i;
                *day = doy - idoy[i-1] + 1;
                break;
            }
        }

        /* if the month isn't set, then it's a December scene */
        if (*month == 0)
        {
            *month = 12;
            *day = doy - idoy[11] + 1;
        }
    }

    /* Validate the month and day */
    if (*month < 1 || *month > 12)
    {
        sprintf (errmsg, "Invalid month: %d\n", *month);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (leap)
    {  /* leap year */
        if (*day < 1 || *day > nday_lp[(*month)-1])
        {
            sprintf (errmsg, "Invalid day: %d-%d-%d\n", year, *month, *day);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }
    else
    {  /* non leap year */
        if (*day < 1 || *day > nday[(*month)-1])
        {
            sprintf (errmsg, "Invalid day: %d-%d-%d\n", year, *month, *day);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Successful conversion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  read_geospatial_lat_lon

PURPOSE: Reads the geospatial lat/long bounding coordinates from the
geospatial_lat_lon_extent container.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reading the bounding coords from the netCDF file
SUCCESS         Successfully obtained the bounding coords from the netCDF file

NOTES:
******************************************************************************/
int read_geospatial_lat_lon
(
    int ncid,                 /* I: netCDF file ID */
    int primary_index,        /* I: index of the primary variable */
    double *bound_coords      /* O: bounding coordinates */
)
{
    char FUNC_NAME[] = "read_geospatial_lat_lon";  /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int status;              /* return value */

    /* Read the geospatial_westbound_longitude, geospatial_northbound_latitude,
       geospatial_eastbound_longitude, and geospatial_southbound_latitude
       attributes from the primary variable */
    /* Westbound extent */
    if ((status = nc_get_att_double (ncid, primary_index,
         "geospatial_westbound_longitude", &bound_coords[ESPA_WEST])))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the "
            "geospatial_westbound_longitude attribute value from the "
            "primary variable geospatial_lat_lon_extent.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Eastbound extent */
    if ((status = nc_get_att_double (ncid, primary_index,
         "geospatial_eastbound_longitude", &bound_coords[ESPA_EAST])))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the "
            "geospatial_eastbound_longitude attribute value from the "
            "primary variable geospatial_lat_lon_extent.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Northbound extent */
    if ((status = nc_get_att_double (ncid, primary_index,
         "geospatial_northbound_latitude", &bound_coords[ESPA_NORTH])))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the "
            "geospatial_northbound_latitude attribute value from the "
            "primary variable geospatial_lat_lon_extent.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Southbound extent */
    if ((status = nc_get_att_double (ncid, primary_index,
         "geospatial_southbound_latitude", &bound_coords[ESPA_SOUTH])))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the "
            "geospatial_southbound_latitude attribute value from the "
            "primary variable geospatial_lat_lon_extent.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful conversion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  read_projection_info

PURPOSE: Reads the projection information from the goes_imager_projection
variable.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reading the projection info from the netCDF file
SUCCESS         Successfully obtained the projection info from the netCDF file

NOTES:
******************************************************************************/
int read_projection_info
(
    int ncid,                    /* I: netCDF file ID */
    int primary_index,           /* I: index of the primary variable */
    Espa_proj_meta_t *proj_info  /* O: projection information */
)
{
    char FUNC_NAME[] = "read_projection_info";  /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char tmpstr[STR_SIZE];   /* temporary string for reading attributes */
    char attname[STR_SIZE];  /* attribute name to process */
    int status;              /* return value */
    size_t attlen = 0;       /* length of the string attribute */

    /* Read the grid_mapping_name, semi_major_axis, semi_minor_axis,
       longitude_of_projection_origin, and perspective_point_height attributes
       from the primary variable */

    /* Grid mapping name */
    strcpy (attname, "grid_mapping_name");
    nc_inq_attlen (ncid, primary_index, attname, &attlen);
    if ((status = nc_get_att_text (ncid, primary_index, attname, tmpstr)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the %s attribute value from the "
            "primary variable goes_imager_projection.", attname);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    tmpstr[attlen] = '\0';

    if (strcmp (tmpstr, "geostationary"))
    {
        sprintf (errmsg, "Invalid projection type: %s.  GOES-R ABI data is "
            "expected to be in the geostationary projection.", tmpstr);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    proj_info->proj_type = ESPA_GEOSTATIONARY;

    /* Semi-major axis */
    if ((status = nc_get_att_double (ncid, primary_index, "semi_major_axis",
         &proj_info->semi_major_axis)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the semi_major_axis attribute "
            "value from the primary variable goes_imager_projection.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Semi-minor axis */
    if ((status = nc_get_att_double (ncid, primary_index, "semi_minor_axis",
         &proj_info->semi_minor_axis)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the semi_minor_axis attribute "
            "value from the primary variable goes_imager_projection.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Central meridian */
    if ((status = nc_get_att_double (ncid, primary_index,
         "longitude_of_projection_origin", &proj_info->central_meridian)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the "
            "longitude_of_projection_origin attribute value from the primary "
            "variable goes_imager_projection.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Satellite height */
    if ((status = nc_get_att_double (ncid, primary_index,
         "perspective_point_height", &proj_info->satellite_height)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain the perspective_point_height "
            "attribute value from the primary variable "
            "goes_imager_projection.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Check the semi-major and semi-minor to make sure they match up with
       what is expected for GRS80 */
    if (GCTP_GRS80_SEMI_MAJOR - proj_info->semi_major_axis > 0.0001)
    {
        sprintf (errmsg, "Semi-major axis is unexpected value of %f. Value "
            "was expected to be %f, which is the GRS80 semi-major axis.",
            GCTP_GRS80_SEMI_MAJOR, proj_info->semi_major_axis);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (GCTP_GRS80_SEMI_MINOR - proj_info->semi_minor_axis > 0.0001)
    {
        sprintf (errmsg, "Semi-minor axis is unexpected value of %f. Value "
            "was expected to be %f, which is the GRS80 semi-minor axis.",
            GCTP_GRS80_SEMI_MINOR, proj_info->semi_minor_axis);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Store the remaining projection info.  The data isn't in any particular
       datum.  And the false easting and northing are 0.0. */
    proj_info->datum_type = ESPA_NODATUM;
    proj_info->false_easting = 0.0;
    proj_info->false_northing = 0.0;

    /* Successful conversion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  read_image_bounds

PURPOSE: Reads the x/y image space bounding coordinates from the
[x|y_image_bounds] variable

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reading the image bounds from the netCDF file
SUCCESS         Successfully obtained the image bounds from the netCDF file

NOTES:
******************************************************************************/
int read_image_bounds
(
    int ncid,                 /* I: netCDF file ID */
    int primary_index,        /* I: index of the primary variable */
    float *image_bounds       /* O: image bounds from specified variable */
)
{
    char FUNC_NAME[] = "read_image_bounds";  /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int status;              /* return value */
    size_t start[1];         /* starting location for reading data */
    size_t count[1];         /* how many values within each dimension will be
                                read */

    /* Read the array of two floating point values for the min/max bounds*/
    /* Initialize start variable to start reading at 0 index */
    start[0] = 0;    /* x */

    /* Read the entire variable */
    count[0] = 2;    /* x */

    /* Read the data */
    if ((status = nc_get_vara_float (ncid, primary_index, start, count,
         image_bounds)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Error reading image bounds from primary variable");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful conversion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  read_ncdf_metadata

PURPOSE: Reads the product version, projection information, bounding
coordinates, and lat/long extents from the netCDF file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reading the product version from the netCDF file
SUCCESS         Successfully obtained the product version from the netCDF file

NOTES:
******************************************************************************/
int read_ncdf_metadata
(
    int ncid,               /* I: netCDF file ID */
    char *product_version,  /* O: product version */
    double *bound_coords,   /* O: bounding coordinates */
    Espa_proj_meta_t *proj_info,  /* O: projection information */
    float *x_image_bounds,  /* O: x-direction image bounds */
    float *y_image_bounds   /* O: y-direction image bounds */
)
{
    char FUNC_NAME[] = "read_ncdf_metadata";  /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variable */
    int status;              /* return value */
    int ndims;               /* number of input dimensions in netCDF file */
    int nvars;               /* number of input variables in netCDF file */
    int ngatts;              /* number of global attributes in netCDF file */
    int unlimdimid;          /* ID of the unlimited dimension */
    int primary_index;       /* index of the primary variable */
    size_t attlen = 0;       /* length of the string attribute */
    char attname[STR_SIZE];  /* attribute name to process */
    char in_varnames[MAX_GRIDDED_NVARS][NC_MAX_NAME+1]; /* var names as read */
    int in_var_ndims[MAX_GRIDDED_NVARS];  /* num dims for each var as read */
    int in_var_dimids[MAX_GRIDDED_NVARS][MAX_GRIDDED_NDIMS];
                             /* array for the dimension IDs as read */
    int in_var_natts[MAX_GRIDDED_NVARS];  /* number of var attributes as read */
    nc_type in_data_type[MAX_GRIDDED_NVARS];
                             /* data type for each variable as read */

    /* Determine how many netCDF variables, dimensions, and global attributes
       are in the file; also the dimension id of the unlimited dimension, if
       there is one. */
    if ((status = nc_inq (ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Error inquiring about the variables, dimensions, "
            "global attributes, etc. for the netCDF file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
//printf ("Number of variables: %d\n", nvars);
//printf ("Number of dims: %d\n", ndims);
//printf ("MAX_GRIDDED_NDIMS: %d\n", MAX_GRIDDED_NDIMS);
//printf ("MAX_GRIDDED_NVARS: %d\n", MAX_GRIDDED_NVARS);

    /* Get information about the variables in the file */
    primary_index = -1;     /* initialize index to invalid value */
    for (i = 0; i < nvars; i++)
    {
        if ((status = nc_inq_var (ncid, i, in_varnames[i],
             &in_data_type[i], &in_var_ndims[i], in_var_dimids[i],
             &in_var_natts[i])))
        {
            nc_strerror (status);
            sprintf (errmsg, "Error inquiring about variable %d", i);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
//printf ("Variable %d/%d: %s\n", i, nvars, in_varnames[i]);
//printf ("  ndims - %d\n", in_var_ndims[i]);
//if (in_var_ndims[i] == 2) printf ("  dims - %d %d\n", in_var_dimids[i][0], in_var_dimids[i][1]);
//printf ("  natts - %d\n", in_var_natts[i]);

        /* If this variable is algorithm_product_version_container then
           process it */
        if (!strcmp (in_varnames[i], "algorithm_product_version_container"))
        {
//printf ("  **algorithm_product_version_container found\n");
            primary_index = i;

            /* Read the product_version attribute */
            strcpy (attname, "product_version");
            nc_inq_attlen (ncid, primary_index, attname, &attlen);
            if ((status = nc_get_att_text (ncid, primary_index, attname,
                 product_version)))
            {
                nc_strerror (status);
                sprintf (errmsg, "Not able to obtain the %s attribute value "
                    "from the primary variable "
                    "algorithm_product_version_container.", attname);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }
        product_version[attlen] = '\0';

        /* If this variable is geospatial_lat_lon_extent then process it */
        if (!strcmp (in_varnames[i], "geospatial_lat_lon_extent"))
        {
//printf ("  **geospatial_lat_lon_extent found\n");
            primary_index = i;

            /* Read the bounding coordinates from the
               geospatial_lat_lon_extent */
            status = read_geospatial_lat_lon (ncid, primary_index,
                bound_coords);
            if (status != SUCCESS)
            {
                sprintf (errmsg, "Reading the bounding coordinates.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* If this variable is goes_imager_projection then process it */
        if (!strcmp (in_varnames[i], "goes_imager_projection"))
        {
//printf ("  **goes_imager_projection found\n");
            primary_index = i;

            /* Read the projection information, which will be used for the
               band data */
            status = read_projection_info (ncid, primary_index, proj_info);
            if (status != SUCCESS)
            {
                sprintf (errmsg, "Reading the projection information");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* If this variable is x_image_bounds then process it */
        if (!strcmp (in_varnames[i], "x_image_bounds"))
        {
//printf ("  **x_image_bounds variable found\n");
            primary_index = i;

            /* Read x image bounds for east/west projection extents.  These
               represent the outer extents of the image and will need to be
               adjusted to represent the center of the pixel. */
            status = read_image_bounds (ncid, primary_index, x_image_bounds);
            if (status != SUCCESS)
            {
                sprintf (errmsg, "Reading x_image_bounds east/west coords.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
printf ("DEBUG: x_image_bounds: %f, %f\n", x_image_bounds[0], x_image_bounds[1]);
        }

        /* If this variable is y_image_bounds then process it */
        if (!strcmp (in_varnames[i], "y_image_bounds"))
        {
//printf ("  **y_image_bounds variable found\n");
            primary_index = i;

            /* Read y image bounds for north/south projection extents.  These
               represent the outer extents of the image and will need to be
               adjusted to represent the center of the pixel. */
            status = read_image_bounds (ncid, primary_index, y_image_bounds);
            if (status != SUCCESS)
            {
                sprintf (errmsg, "Reading y_image_bounds north/south coords.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
printf ("DEBUG: y_image_bounds: %f, %f\n", y_image_bounds[0], y_image_bounds[1]);
        }
    }  /* for i in nvars */

    /* Successful conversion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  read_goes_netcdf

PURPOSE: Read the metadata from the GOES-R netCDF file and populate the ESPA
internal metadata structure

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reading the GOES-R ABI file
SUCCESS         Successfully populated the ESPA metadata structure

NOTES:
******************************************************************************/
int read_goes_netcdf
(
    char *ncdf_filename,               /* I: netCDF filename */
    int ncid,                          /* I: netCDF file ID */
    Espa_ncdf_var_attr_t ncdf_attr[2], /* I: file attributes for the netCDF
                                          primary variables
                                          ([0]-CMI, [1]-DQF) */
    char *xml_filename,                /* XML filename */
    Espa_internal_meta_t *metadata     /* I/O: input metadata structure to be
                                          populated from the GOES-R ABI file */
)
{
    char FUNC_NAME[] = "read_goes_netcdf";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    char xml_basename[STR_SIZE]; /* filename without path and extension */
    char product_version[STR_SIZE]; /* product version */
    char production_date[STR_SIZE]; /* production date of the data */
    char project[STR_SIZE];   /* project global attribute */
    char scene_id[STR_SIZE];  /* scene_id global attribute */
    char spatial_resolution[STR_SIZE]; /* spatial_resolution global attribute */
    char bandname[STR_SIZE];  /* GOES band name for the current file */
    char tmpstr[STR_SIZE];    /* temporary string for date/time info */
    char attname[STR_SIZE];   /* attribute name to process */
    char goesstr[3];          /* string to hold the GOES instrument number */
    char yearstr[5];          /* string to hold acquisition/production year */
    char doystr[4];           /* string to hold acquisition/production DOY */
    char timestr[8];          /* string to hold acquisition/production time */
    char *cptr = NULL;        /* character pointer for strings */
    int i;                    /* looping variable */
    int status;               /* return status */
    int count;                /* number of chars copied in snprintf */
    int acq_doy;              /* acquisition DOY */
    int acq_year;             /* acquisition year */
    int acq_month;            /* acquisition month */
    int acq_day;              /* acquisition day */
    int prod_doy;             /* production DOY */
    int prod_year;            /* production year */
    int prod_month;           /* production month */
    int prod_day;             /* production day */
    int prod_hour;            /* production hour */
    int prod_min;             /* production minute */
    int prod_sec;             /* production seconds */
    size_t attlen = 0;        /* length of the string attribute */
    float x_image_bounds[2];  /* west/east image coordinates for the x dim */
    float y_image_bounds[2];  /* north/south image coordinates for the y dim */

    Espa_global_meta_t *gmeta = &metadata->global;  /* pointer to the global
                                                       metadata structure */
    Espa_band_meta_t *bmeta;      /* pointer to the array of bands metadata */

    /* Strip the extension off the XML file to get the xml_basename */
    strcpy (xml_basename, xml_filename);
    cptr = strrchr (xml_basename, '.');
    if (cptr != NULL)
    {
        /* File extension found */
        *cptr = '\0';
    }
    else
    {
        sprintf (errmsg, "Unexpected XML filename: %s. Filenames are expected "
            "to have a file extension.", xml_filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Use the netCDF filename to determine the data_provider, satellite,
       instrument, and product.
       Example - OR_ABI-L2-CMIPC-M3C02_G16_s20171721702192_e20171721704565_c20171721705067.nc */
    strcpy (gmeta->data_provider, "NOAA NESDIS");
    strcpy (gmeta->instrument, "GOES-R Advanced Baseline Imager");

    if (strncpy (goesstr, &ncdf_filename[23], 2) == NULL)
    {
        sprintf (errmsg, "Error pulling the GOES number from the netCDF "
            "filename: %s", ncdf_filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    goesstr[2] = '\0';
    sprintf (gmeta->satellite, "GOES-%s", goesstr);

    /* Use the netCDF filename to determine the acquisition date as yyyyddd.
       We will use the start of acquisition, which is available with the '_s'
       section of the filename.  '_e' is the ending time. */
    if (strncpy (yearstr, &ncdf_filename[27], 4) == NULL)
    {
        sprintf (errmsg, "Error pulling the acquisition year from the netCDF "
            "filename: %s", ncdf_filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    yearstr[4] = '\0';
    acq_year = atoi (yearstr);

    if (strncpy (doystr, &ncdf_filename[31], 3) == NULL)
    {
        sprintf (errmsg, "Error pulling the acquisition DOY from the base "
            "filename: %s", ncdf_filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    doystr[3] = '\0';
    acq_doy = atoi (doystr);

    /* Year and DOY need to be converted to yyyy-mm-dd */
    if (doy_to_month_day (acq_year, acq_doy, &acq_month, &acq_day) != SUCCESS)
    {
        sprintf (errmsg, "Error converting %d-%d to yyyy-mm-dd", acq_year,
            acq_doy);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    count = snprintf (gmeta->acquisition_date, sizeof (gmeta->acquisition_date),
        "%04d-%02d-%02d", acq_year, acq_month, acq_day);
    if (count < 0 || count >= sizeof (gmeta->acquisition_date))
    {
        sprintf (errmsg, "Overflow of gmeta->acquisition_date string");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Use the netCDF filename to determine the production date/time, which
       comes with the '_c' section of the filename. */
    if (strncpy (yearstr, &ncdf_filename[59], 4) == NULL)
    {
        sprintf (errmsg, "Error pulling the production year from the netCDF "
            "filename: %s", ncdf_filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    yearstr[4] = '\0';
    prod_year = atoi (yearstr);

    if (strncpy (doystr, &ncdf_filename[63], 3) == NULL)
    {
        sprintf (errmsg, "Error pulling the production DOY from the base "
            "filename: %s", ncdf_filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    doystr[3] = '\0';
    prod_doy = atoi (doystr);

    if (strncpy (timestr, &ncdf_filename[66], 7) == NULL)
    {
        sprintf (errmsg, "Error pulling the production date/time from the "
            "netCDF filename");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    timestr[7] = '\0';

    /* Year and DOY need to be converted to yyyy-mm-dd */
    if (doy_to_month_day (prod_year, prod_doy, &prod_month, &prod_day) !=
        SUCCESS)
    {
        sprintf (errmsg, "Error converting %d-%d to yyyy-mm-dd", prod_year,
            prod_doy);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Production time string needs to be in hours, minutes, seconds.  The
       seconds in the netCDF filename are three characters (xyz) representing
       the floating point seconds with one decimal place (xy.z).  This needs
       to be truncated to just two digits in the output production time
       string.  NOTE: Rounding would be more accurate, but would also require
       possibly changing the minutes, hours, and date if the time were to round
       upward.  For simplicity, we will just truncate. */
    strncpy (tmpstr, &timestr[0], 2);
    prod_hour = atoi (tmpstr);

    strncpy (tmpstr, &timestr[2], 2);
    prod_min = atoi (tmpstr);

    strncpy (tmpstr, &timestr[4], 2);
    prod_sec = atoi (tmpstr);

    /* Store the production date/time as YYYY-MM-DDThh:mm:ssZ */
    count = snprintf (production_date, sizeof (production_date),
        "%04d-%02d-%02dT%02d:%02d:%02dZ", prod_year, prod_month, prod_day,
        prod_hour, prod_min, prod_sec);
    if (count < 0 || count >= sizeof (production_date))
    {
        sprintf (errmsg, "Overflow of production_date string");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Use the netCDF filename to determine which band this file is */
    if (!strncmp (&ncdf_filename[18], "C02", 3))
        strcpy (bandname, "b2");
    else if (!strncmp (&ncdf_filename[18], "C03", 3))
        strcpy (bandname, "b3");
    else
    {
        sprintf (errmsg, "Unexpected GOES-R ABI filename.  It is expected that "
            "this is either channel 2 (red) or channel 3 (NIR) data.  %s",
            ncdf_filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the various global and variable attributes from the netCDF file */
    status = read_ncdf_metadata (ncid, product_version, gmeta->bounding_coords,
        &gmeta->proj_info, x_image_bounds, y_image_bounds);
    if (status != SUCCESS)
    {
        sprintf (errmsg, "Reading the global and variable netCDF attributes.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the project global attribute */
    strcpy (attname, "project");
    nc_inq_attlen (ncid, NC_GLOBAL, attname, &attlen);
    if ((status = nc_get_att_text (ncid, NC_GLOBAL, attname, project)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain %s global attribute.", attname);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    project[attlen] = '\0';

    /* Read the scene_id global attribute */
    strcpy (attname, "scene_id");
    nc_inq_attlen (ncid, NC_GLOBAL, attname, &attlen);
    if ((status = nc_get_att_text (ncid, NC_GLOBAL, attname, scene_id)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain %s global attribute.", attname);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    scene_id[attlen] = '\0';

    /* Switch 'Full Disk' to uppercase.  Check for Full Disk or CONUS.  Those
       are the only two supported scene types. */
    if (!strcmp (scene_id, "Full Disk"))
        strcpy (scene_id, "goes_full_disk");
    else if (!strcmp (scene_id, "CONUS"))
        strcpy (scene_id, "goes_conus");
    else
    {
        sprintf (errmsg, "Unexpected scene_id type %s.  Only Full Disk and "
            "CONUS are expected.", scene_id);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the spatial_resolution global attribute */
    strcpy (attname, "spatial_resolution");
    nc_inq_attlen (ncid, NC_GLOBAL, attname, &attlen);
    if ((status = nc_get_att_text (ncid, NC_GLOBAL, attname,
         spatial_resolution)))
    {
        nc_strerror (status);
        sprintf (errmsg, "Not able to obtain %s global attribute.", attname);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    spatial_resolution[attlen] = '\0';

    /* Adjust the UL and LR outer extent coordinates by the pixel size to
       represent the center of the pixel.  Use the pixel size from the CMI
       variable. */
    x_image_bounds[0] += ncdf_attr[GOES_CMI].pixel_size[0] * 0.5;
    x_image_bounds[1] -= ncdf_attr[GOES_CMI].pixel_size[0] * 0.5;
    y_image_bounds[0] -= ncdf_attr[GOES_CMI].pixel_size[1] * 0.5;
    y_image_bounds[1] += ncdf_attr[GOES_CMI].pixel_size[1] * 0.5;
printf ("DEBUG: adjusted x_image_bounds: %f, %f\n", x_image_bounds[0], x_image_bounds[1]);
printf ("DEBUG: adjusted y_image_bounds: %f, %f\n", y_image_bounds[0], y_image_bounds[1]);

    /* Projection information is in radians */
    strcpy (gmeta->proj_info.units, "radians");
    gmeta->proj_info.ul_corner[0] = x_image_bounds[0];
    gmeta->proj_info.ul_corner[1] = y_image_bounds[0];
    gmeta->proj_info.lr_corner[0] = x_image_bounds[1];
    gmeta->proj_info.lr_corner[1] = y_image_bounds[1];
    strcpy (gmeta->proj_info.grid_origin, "CENTER");
    gmeta->orientation_angle = 0.0;

    /* Allocate bands for the XML structure */
    metadata->nbands = NGOES_BANDS;
    if (allocate_band_metadata (metadata, metadata->nbands) != SUCCESS)
    {   /* Error messages already printed */
        return (ERROR);
    }
    bmeta = metadata->band;

    /* Loop through the bands (CMI, DQF), fill in the information from above */
    for (i = 0; i < metadata->nbands; i++)
    {
        /* Fill in the band level information already obtained. Use the global
           scene_id attribute for the product type. Use GOESABI for the short
           name.  Use 'image' for the CMI category and 'qa' for the DQF
           category. */
        strcpy (bmeta[i].product, scene_id);
        strcpy (bmeta[i].short_name, "GOESABI");
        bmeta[i].nsamps = ncdf_attr[i].nsamps;
        bmeta[i].nlines = ncdf_attr[i].nlines;

        /* Setup the band name. If this is the DQF band then append _dqf to the
           band name.  Use 'image' for the CMI category and 'qa' for the DQF
           category.  The DQF band needs to contain the DQF extension. */
        if (i == GOES_DQF)
        {
            strcpy (bmeta[i].category, "qa");

            count = snprintf (bmeta[i].name, sizeof (bmeta[i].name), "%s_dqf",
                bandname);
            if (count < 0 || count >= sizeof (bmeta[i].name))
            {
                sprintf (errmsg, "Overflow of bmeta[].name string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            count = snprintf (bmeta[i].long_name,
                sizeof (bmeta[i].long_name), "band %c data quality flag",
                bandname[1]);
            if (count < 0 || count >= sizeof (bmeta[i].long_name))
            {
                sprintf (errmsg, "Overflow of bmeta[].long_name string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            count = snprintf (bmeta[i].data_units, sizeof (bmeta[i].data_units),
                "quality/feature classification");
            if (count < 0 || count >= sizeof (bmeta[i].data_units))
            {
                sprintf (errmsg, "Overflow of bmeta[].data_units string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Set up the 4 classes for the data quality flag */
            bmeta[i].nclass = 4;
            if (allocate_class_metadata (&bmeta[i], 4) != SUCCESS)
            {
                sprintf (errmsg, "Cannot allocate memory for the DQF classes");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            bmeta[i].class_values[0].class = 0;
            bmeta[i].class_values[1].class = 1;
            bmeta[i].class_values[2].class = 2;
            bmeta[i].class_values[3].class = 3;
            strcpy (bmeta[i].class_values[0].description, "good pixel quality");
            strcpy (bmeta[i].class_values[1].description,
                "conditionally usable");
            strcpy (bmeta[i].class_values[2].description, "out of range");
            strcpy (bmeta[i].class_values[3].description, "no value");
        }
        else
        {
            strcpy (bmeta[i].category, "image");
            count = snprintf (bmeta[i].name, sizeof (bmeta[i].name), "%s",
                bandname);
            if (count < 0 || count >= sizeof (bmeta[i].name))
            {
                sprintf (errmsg, "Overflow of bmeta[].name string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            count = snprintf (bmeta[i].long_name,
                sizeof (bmeta[i].long_name), "band %c radiance", bandname[1]);
            if (count < 0 || count >= sizeof (bmeta[i].long_name))
            {
                sprintf (errmsg, "Overflow of bmeta[].long_name string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            count = snprintf (bmeta[i].data_units, sizeof (bmeta[i].data_units),
                "radiance");
            if (count < 0 || count >= sizeof (bmeta[i].data_units))
            {
                sprintf (errmsg, "Overflow of bmeta[].data_units string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* Set up the filename */
        count = snprintf (bmeta[i].file_name, sizeof (bmeta[i].file_name),
            "%s_%s.img", xml_basename, bmeta[i].name);
        if (count < 0 || count >= sizeof (bmeta[i].file_name))
        {
            sprintf (errmsg, "Overflow of bmeta[].file_name string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Convert the data type to the ESPA data type */
        switch (ncdf_attr[i].native_data_type)
        {
            case NC_BYTE:
                bmeta[i].data_type = ESPA_INT8; break;
            case NC_UBYTE:
                bmeta[i].data_type = ESPA_UINT8; break;
            case NC_SHORT:
                bmeta[i].data_type = ESPA_INT16; break;
            case NC_USHORT:
                bmeta[i].data_type = ESPA_UINT16; break;
            case NC_INT:
                bmeta[i].data_type = ESPA_INT32; break;
            case NC_UINT:
                bmeta[i].data_type = ESPA_UINT32; break;
            case NC_FLOAT:
                bmeta[i].data_type = ESPA_FLOAT32; break;
            default:
                sprintf (errmsg, "Unsupported GOES-R ABI netCDF data type");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
        }

        /* Compute the pixel size */
        bmeta[i].pixel_size[0] = ncdf_attr[GOES_CMI].pixel_size[0];
        bmeta[i].pixel_size[1] = ncdf_attr[GOES_CMI].pixel_size[1];
        strcpy (bmeta[i].pixel_units, "radians");

        /* Assign the scale, offset, min/max, and fill values */
        if (ncdf_attr[i].fill_defined)
            bmeta[i].fill_value = ncdf_attr[i].fill_value;
        if (ncdf_attr[i].scale_defined)
            bmeta[i].scale_factor = ncdf_attr[i].scale_fact;
        if (ncdf_attr[i].offset_defined)
            bmeta[i].add_offset = ncdf_attr[i].add_offset;
        if (ncdf_attr[i].valid_range_defined)
        {
            bmeta[i].valid_range[0] = ncdf_attr[i].valid_range[0];
            bmeta[i].valid_range[1] = ncdf_attr[i].valid_range[1];
        }

        /* Add production date/time and product version from the netCDF file */
        count = snprintf (bmeta[i].production_date,
            sizeof (bmeta[i].production_date), "%s", production_date);
        if (count < 0 || count >= sizeof (bmeta[i].production_date))
        {
            sprintf (errmsg, "Overflow of bmeta[].production_date string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        count = snprintf (bmeta[i].app_version, sizeof (bmeta[i].app_version),
            "Product version %s", product_version);
        if (count < 0 || count >= sizeof (bmeta[i].app_version))
        {
            sprintf (errmsg, "Overflow of bmeta[].app_version string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }  /* end for i (loop through the grids) */

    /* Successful read */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  setup_geo_xml

PURPOSE: Copies the GOES XML metadata and sets up the geographic projection
information for the reprojected data.  The UL/LR corners are set up to be the
same as the bounding coordinates.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error copying the GOES metadata
SUCCESS         Successfully set up the geographic metadata

NOTES:
******************************************************************************/
int setup_geo_xml
(
    Espa_internal_meta_t *goes_xml,  /* I: GOES metadata structure with input
                                           projection and corner information */
    Espa_internal_meta_t *xml_metadata  /* O: output metadata structure for
                                              reprojected data */
)
{
    char FUNC_NAME[] = "setup_geo_xml";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    Espa_global_meta_t *gmeta = &xml_metadata->global;  /* reprojected global
                                                           metadata */

    /* Copy the input GOES metadata structure */
    if (copy_metadata_struct (goes_xml, xml_metadata) != SUCCESS)
    {
        sprintf (errmsg, "Copying the GOES XML metadata");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Update the projection information to be geographic instead of the
       geostationary projection */
    gmeta->proj_info.proj_type = GCTP_GEO_PROJ;
    gmeta->proj_info.datum_type = ESPA_NAD83;
    strcpy (gmeta->proj_info.units, "degrees");

    /* Use the bounding lat/long coordinates as the UL and LR corners in our
       geographic product */
    gmeta->ul_corner[0] = gmeta->bounding_coords[ESPA_NORTH];  /* lat */
    gmeta->ul_corner[1] = gmeta->bounding_coords[ESPA_WEST];   /* lon */
    gmeta->lr_corner[0] = gmeta->bounding_coords[ESPA_SOUTH];  /* lat */
    gmeta->lr_corner[1] = gmeta->bounding_coords[ESPA_EAST];   /* lon */

    /* Use the bounding lat/long coordinates as the UL and LR corners in the
       projection information */
    gmeta->proj_info.ul_corner[0] = gmeta->bounding_coords[ESPA_WEST];  /* x */
    gmeta->proj_info.ul_corner[1] = gmeta->bounding_coords[ESPA_NORTH]; /* y */
    gmeta->proj_info.lr_corner[0] = gmeta->bounding_coords[ESPA_EAST];  /* x */
    gmeta->proj_info.lr_corner[1] = gmeta->bounding_coords[ESPA_SOUTH]; /* y */

    /* Clear out the geostationary projection parameters */
    gmeta->proj_info.semi_major_axis = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.semi_minor_axis = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.satellite_height = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.central_meridian = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.false_easting = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.false_northing = ESPA_FLOAT_META_FILL;

    /* Successful setup */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  convert_netcdf_to_img

PURPOSE: Convert the GOES-R ABI netCDF band to an ESPA raw binary (.img) file
and write the associated ENVI header for this band.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error converting the GOES-R ABI band
SUCCESS         Successfully converted GOES-R ABI band to raw binary

NOTES:
1. It is assumed the netCDF file is already opened.
******************************************************************************/
int convert_netcdf_to_img
(
    int ncid,            /* I: netCDF file ID */
    int primary_index,   /* I: index of the primary variable */
    nc_type native_data_type, /* I: native data type of the primary variable */
    int xml_band,        /* I: which band in the XML file is being processed */
    Espa_internal_meta_t *goes_xml,  /* I: GOES metadata structure with input
                                           projection and corner information */
    Espa_internal_meta_t *xml_metadata  /* I/O: output metadata structure for
                                                reprojected data */
)
{
    char FUNC_NAME[] = "convert_netcdf_to_img";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    char *cptr = NULL;        /* pointer to the file extension */
    char *img_file = NULL;    /* name of the output raw binary file */
    char envi_file[STR_SIZE]; /* name of the output ENVI header file */
    int nbytes;               /* number of bytes in the data type */
    int count;                /* number of chars copied in snprintf */
    void *file_buf = NULL;    /* pointer to input file buffer */
    void *reprojected_file_buf = NULL;  /* pointer to reprojected file buffer */
    FILE *fp_rb = NULL;       /* file pointer for the raw binary file */
    Envi_header_t envi_hdr;   /* output ENVI header information */
    Espa_band_meta_t *bmeta = NULL;  /* pointer to reprojected band metadata */
    Espa_band_meta_t *goes_bmeta = NULL;  /* pointer to GOES band metadata */
    Espa_global_meta_t *gmeta = &xml_metadata->global;  /* reprojected global
                                                           metadata */
    Espa_global_meta_t *goes_gmeta = &goes_xml->global; /* global GOES meta */

    /* Set up the band metadata pointer */
    bmeta = &xml_metadata->band[xml_band];
    goes_bmeta = &goes_xml->band[xml_band];

    /* Open the raw binary file for writing */
    img_file = bmeta->file_name;
    fp_rb = open_raw_binary (img_file, "wb");
    if (fp_rb == NULL)
    {
        sprintf (errmsg, "Opening the output raw binary file: %s", img_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the gridded data from the primary index in the netCDF file.
       Memory is allocated for the data buffer before the dataset is read. */
    if (ncdf_read_gridded_var (ncid, primary_index, goes_bmeta->nlines,
        goes_bmeta->nsamps, native_data_type, &file_buf) != SUCCESS)
    {
        sprintf (errmsg, "Reading the gridded netCDF data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Reproject the image from Geostationary to Geographic lat/long */
    if (reproject_goes (goes_gmeta, goes_bmeta, gmeta, bmeta, file_buf,
        &reprojected_file_buf) != SUCCESS)
    {
        sprintf (errmsg, "Reprojecting the gridded netCDF data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Determine the number of bytes per pixel for the output product */
    if (bmeta->data_type == ESPA_INT8)
        nbytes = sizeof (int8_t);
    else if (bmeta->data_type == ESPA_INT16)
        nbytes = sizeof (int16_t);
    else
    {
        sprintf (errmsg, "Unsupported GOES-R ABI data type %d.",
            bmeta->data_type);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write reprojected image to the raw binary file */
    if (write_raw_binary (fp_rb, bmeta->nlines, bmeta->nsamps, nbytes,
        reprojected_file_buf) != SUCCESS)
    {
        sprintf (errmsg, "Writing image to the raw binary file: %s", img_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the raw binary file */
    close_raw_binary (fp_rb);

    /* Free the memory */
    free (file_buf);
    free (reprojected_file_buf);

    /* Create the ENVI header file this band */
    if (create_envi_struct (bmeta, gmeta, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating the ENVI header structure for this file: "
            "%s", img_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the ENVI header */
    count = snprintf (envi_file, sizeof (envi_file), "%s", img_file);
    if (count < 0 || count >= sizeof (envi_file))
    {
        sprintf (errmsg, "Overflow of envi_file string");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    cptr = strrchr (envi_file, '.');
    strcpy (cptr, ".hdr");

    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing the ENVI header file: %s.", envi_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful conversion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  convert_goes_to_espa

PURPOSE: Converts the input GOES-R ABI netCDF file to the ESPA internal raw
binary file format (and associated XML file). The data are reprojected from
the input geostationary projection to Geographic lat/long.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error converting the GOES-R ABI file
SUCCESS         Successfully converted GOES-R ABI to ESPA format

NOTES:
  1. The ESPA raw binary band files will be generated from the ESPA XML
     filename.
******************************************************************************/
int convert_goes_to_espa
(
    char *goes_netcdf_file, /* I: input GOES-R ABI netCDF filename */
    char *espa_xml_file,    /* I: output ESPA XML metadata filename; if
                                  append_bands, then the XML metadata file
                                  will be appended to; otherwise the XML
                                  metadata file will be created */
    bool append_bands,      /* I: are the bands in this GOES-R file being
                                  appended?  If not, then the XML file will
                                  be created from scratch. */
    bool del_src            /* I: should the source .nc files be removed after
                                  conversion? */
)
{
    char FUNC_NAME[] = "convert_goes_to_espa";  /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int ncid;                /* netCDF file ID */
    int xml_band;            /* band number for the current file in the XML */
    Espa_internal_meta_t goes_xml; /* XML metadata structure for temporarily
                                      holding the GOES input metadata */
    Espa_internal_meta_t xml_metadata; /* XML metadata structure to be
                                          for the final Geographic product */
    Espa_ncdf_var_attr_t ncdf_attr[2]; /* file attributes for the netCDF primary
                                          variable ([0] - CMI, [1] - DQF) */

    /* Initialize the input and output metadata structures */
    init_metadata_struct (&goes_xml);
    init_metadata_struct (&xml_metadata);

    /* Open the GOES-R ABI netCDF file for reading */
    if (open_netcdf (goes_netcdf_file, &ncid) != SUCCESS)
    {
        sprintf (errmsg, "Opening the GOES-R ABI netCDF file: %s",
            goes_netcdf_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the attributes for CMI */
    printf ("Reading CMI variable attributes\n");
    if (ncdf_read_gridded_attr (ncid, "CMI", &ncdf_attr[GOES_CMI]) != SUCCESS)
    {
        sprintf (errmsg, "Reading attributes for CMI variable in the GOES-R "
            "ABI netCDF file: %s", goes_netcdf_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the attributes for DQF */
    printf ("Reading DQF variable attributes\n");
    if (ncdf_read_gridded_attr (ncid, "DQF", &ncdf_attr[GOES_DQF]) != SUCCESS)
    {
        sprintf (errmsg, "Reading attributes for DQF variable in the GOES-R "
            "ABI netCDF file: %s", goes_netcdf_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the attributes from the netCDF file into the GOES metadata struct */
    printf ("Reading global and related variable attributes\n");
    if (read_goes_netcdf (goes_netcdf_file, ncid, ncdf_attr, espa_xml_file,
        &goes_xml) != SUCCESS)
    {
        sprintf (errmsg, "Reading the GOES-R ABI netCDF file: %s",
            goes_netcdf_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Copy the input metadata to the output metadata, and set the geographic
       projection correctly.  The data are being reprojected from geostationary
       to geographic lat/long. */
    if (setup_geo_xml (&goes_xml, &xml_metadata) != SUCCESS)
    {
        /* Already printed error message */
        return (ERROR);
    }

    /* Convert each of the GOES-R ABI netCDF bands to raw binary.  If we are
       starting with the first GOES-R Channel (ch02/red) file, then the CMI
       and DQF will be the first two bands in the XML file.  If we are
       processing/appending the ch03/nir file, then the CMI and DQF will be
       the third and fourth bands in the output XML file. */
    /* CMI */
    printf ("Converting CMI variable to Geographic lat/long\n");
    xml_band = 0;
    if (convert_netcdf_to_img (ncid, ncdf_attr[GOES_CMI].primary_index,
        ncdf_attr[GOES_CMI].native_data_type, xml_band, &goes_xml,
        &xml_metadata) != SUCCESS)
    {
        sprintf (errmsg, "Converting %s CMI to raw binary", goes_netcdf_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* DQF */
    printf ("Converting DQF variable to Geographic lat/long\n");
    xml_band = 1;
    if (convert_netcdf_to_img (ncid, ncdf_attr[GOES_DQF].primary_index,
        ncdf_attr[GOES_DQF].native_data_type, xml_band, &goes_xml,
        &xml_metadata) != SUCCESS)
    {
        sprintf (errmsg, "Converting %s DQF to raw binary", goes_netcdf_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the output metadata from our internal metadata structure to the
       output XML filename */
    if (append_bands)
    {  /* append the additional two bands */
        if (append_metadata (2, xml_metadata.band, espa_xml_file) != SUCCESS)
        {  /* Error messages already written */
            return (ERROR);
        }
    }
    else
    {  /* create the XML metadata file from scratch */
        if (write_metadata (&xml_metadata, espa_xml_file) != SUCCESS)
        {  /* Error messages already written */
            return (ERROR);
        }
    }

    /* Validate the output metadata file */
    if (validate_xml_file (espa_xml_file) != SUCCESS)
    {  /* Error messages already written */
        return (ERROR);
    }

    /* Close the GOES-R ABI netCDF file */
    if (close_netcdf (ncid) != SUCCESS)
    {
        sprintf (errmsg, "Closing the GOES-R ABI netCDF file: %s",
            goes_netcdf_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Remove the source file if specified */
    if (del_src)
    {
        printf ("  Removing %s\n", goes_netcdf_file);
        if (unlink (goes_netcdf_file) != 0)
        {
            sprintf (errmsg, "Deleting source file: %s", goes_netcdf_file);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Free the metadata structures */
    free_metadata (&xml_metadata);
    free_metadata (&goes_xml);

    /* Successful conversion */
    return (SUCCESS);
}

