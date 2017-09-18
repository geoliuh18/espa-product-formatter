/*****************************************************************************
FILE: netcdf_io.h

PURPOSE: Contains defines and prototypes to read supported netCDF files.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/

#ifndef NETCDF_IO_H
#define NETCDF_IO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "netcdf.h"
#include "error_handler.h"

/* Defines */
typedef unsigned char ubyte;

/* Constant values for the generic gridded variables.  These routines will
   support two dimensions (x,y) that will be supported.  Buffer with extra
   just to support other datasets with possibly multiple numbers of primary
   variables.
   GOES-R ABI products have 5 dimensions and 36 variables */
#define MAX_GRIDDED_NDIMS 10
#define MAX_GRIDDED_NVARS 50

/* Structure for storing the gridded netCDF file attributes */
typedef struct
{
    int primary_index;     /* index of the primary variable */
    int nlines;            /* number of lines in y-dimension */
    int nsamps;            /* number of samples in x-dimension */
    bool pixel_size_defined; /* boolean to indicate the pixel size was
                              defined for the primary variable */
    double pixel_size[2];  /* size of the pixel in meters (assumes a projected
                              grid) for x (samples) and y (lines); NOTE:
                              pixel size is radians for Geostationary */
    nc_type native_data_type; /* native data type of the primary variable */
    bool fill_defined;     /* boolean to indicate if the fill value was
                              defined for the primary variable */
    float fill_value;      /* fill value of this variable if it exists */
    bool scale_defined;    /* boolean to indicate if the scale value was
                              defined for the primary variable */
    double scale_fact;     /* scale factor of this variable if it exists */
    bool offset_defined;   /* boolean to indicate if the offset value was
                              defined for the primary variable */
    double add_offset;     /* offset of this variable if it exists */
    bool valid_range_defined; /* boolean to indicate if the valid range was
                              defined for the primary variable */
    float valid_range[2];  /* valid range for this variable if it exists */
} Espa_ncdf_var_attr_t;

/* Prototypes */
int open_netcdf
(
    char *infile,      /* I: input filename for netCDF data */
    int *ncid          /* O: netCDF file ID */
);

int close_netcdf
(
    int ncid           /* I: netCDF file ID */
);

int ncdf_read_gridded_attr
(
    int ncid,          /* I: netCDF file ID */
    char *varname,     /* I: name of primary variable to be read */
    Espa_ncdf_var_attr_t *attr  /* O: file attributes for the netCDF primary
                                      variable */
);

int ncdf_read_gridded_var
(
    int ncid,            /* I: netCDF file ID */
    int primary_index,   /* I: index of the primary variable */
    int nlines,          /* I: number of lines in y-dimension */
    int nsamps,          /* I: number of samples in x-dimension */
    nc_type native_data_type, /* I: native data type of the primary variable */
    void **nc_data       /* O: address of 2D array of netcdf data values read */
);

#endif
