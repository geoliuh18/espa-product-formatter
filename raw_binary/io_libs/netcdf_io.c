/*****************************************************************************
FILE: netcdf_io.c

PURPOSE: Contains functions for reading GOES-R ABI and VIIRS netCDF products -
variables and associated attributes.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include "netcdf_io.h"

/******************************************************************************
MODULE:  open_netcdf

PURPOSE: Generic function to open a netCDF4 file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error opening the input netCDF file
SUCCESS         Successfully opened the netCDF file

NOTES:
1. The file can be closed by using close_ncdf.
******************************************************************************/
int open_netcdf
(
    char *infile,      /* I: input filename for netCDF data */
    int *ncid          /* O: netCDF file ID */
)
{
    char FUNC_NAME[] = "open_netcdf";  /* function name */
    char errmsg[STR_SIZE];             /* error message */
    int retval;                        /* function call return value */

    /* Open the netCDF file as read-only */
    if ((retval = nc_open (infile, NC_NOWRITE, ncid)))
    {
        nc_strerror (retval);
        sprintf (errmsg, "Error opening netCDF file: %s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  ncdf_read_gridded_attr

PURPOSE: Generic function to open a gridded netCDF4 file and return dimension
and variable information.  This function will support 2D datasets.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reading from the input netCDF files.
SUCCESS         Successfully read the provided netCDF file

NOTES:
1. This function supports [u]byte, [u]short, [u]int, and float data types.
   If other data types are needed, then support will need to be added.
2. The file can be closed by using close_ncdf.
3. The input variable data stored as 2D products are stored as [YDIM][XDIM].
4. The names of the x/y or lat/long dimensions are x, X, y, Y, lat, LAT,
   latitude, LATITUDE, lon, LON, long, LONG, longitude, LONGITUDE.
******************************************************************************/
int ncdf_read_gridded_attr
(
    int ncid,          /* I: netCDF file ID */
    char *varname,     /* I: name of primary variable to be read */
    Espa_ncdf_var_attr_t *attr  /* O: file attributes for the netCDF primary
                                      variable */
)
{
    char FUNC_NAME[] = "ncdf_read_gridded_attr";  /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char in_dimnames[MAX_GRIDDED_NDIMS][NC_MAX_NAME+1]; /* dim names as read */
    char in_varnames[MAX_GRIDDED_NDIMS][NC_MAX_NAME+1]; /* var names as read */
    char tmpstr[STR_SIZE];   /* temporary string for reading text */
    char ignore_str[STR_SIZE]; /* temporary string for reading text */
    int ndims;               /* number of input dimensions in netCDF file */
    int nvars;               /* number of input variables in netCDF file */
    int ngatts;              /* number of global attributes in netCDF file */
    int unlimdimid;          /* ID of the unlimited dimension */
    int primary_ndims=0;     /* number of dims for primary variable */
    int primary_dimids[2];   /* dimension IDs of the expected 2D dataset */
    int in_var_ndims[MAX_GRIDDED_NVARS];  /* num dims for each var as read */
    int in_var_dimids[MAX_GRIDDED_NVARS][MAX_GRIDDED_NDIMS];
                             /* array for the dimension IDs as read */
    int in_var_natts[MAX_GRIDDED_NDIMS];  /* number of var attributes as read */
    int i, d, retval;        /* loop indexes and error handling */
    int no_fill;             /* true if no_fill mode is set for this var */
    nc_type in_data_type[MAX_GRIDDED_NVARS];
                             /* data type for each variable as read */
    size_t in_dimsizes[MAX_GRIDDED_NDIMS];   /* dimension sizes as read
                                for the variable; [0]=z, [1]=y, [2]=x */

    signed char fill_char;   /* fill value for the primary variable */
    ubyte fill_ubyte;        /* fill value for the primary variable */
    short fill_short;        /* fill value for the primary variable */
    ushort fill_ushort;      /* fill value for the primary variable */
    int fill_int;            /* fill value for the primary variable */
    uint fill_uint;          /* fill value for the primary variable */
    float fill_float;        /* fill value for the primary variable */
    void *fill_value=NULL;   /* pointer to fill_{dtype} variable */

    signed char valid_range_char[2]; /* valid range for the primary variable */
    ubyte valid_range_ubyte[2];   /* valid range for the primary variable */
    short valid_range_short[2];   /* valid range for the primary variable */
    ushort valid_range_ushort[2]; /* valid range for the primary variable */
    int valid_range_int[2];       /* valid range for the primary variable */
    uint valid_range_uint[2];     /* valid range for the primary variable */
    float valid_range_float[2];   /* valid range for the primary variable */

    /* Determine how many netCDF variables, dimensions, and global attributes
       are in the file; also the dimension id of the unlimited dimension, if
       there is one. */
    if ((retval = nc_inq (ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
    {
        nc_strerror (retval);
        sprintf (errmsg, "Error inquiring about the variables, dimensions, "
            "global attributes, etc. for the netCDF file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* If the number of dimensions or variables is greater than what we've
       planned for, then flag it because the size of the inq arrays will be
       too small.  To fix, bump the size of the #defines in the netcdf_io.h
       file. */
    if (ndims > MAX_GRIDDED_NDIMS)
    {
        sprintf (errmsg, "Number of dimensions in the input netCDF4 file is "
            "greater than expected.  This software currently supports %d "
            "dimensions, but there are %d dimensions.", MAX_GRIDDED_NDIMS,
            ndims);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (nvars > MAX_GRIDDED_NVARS)
    {
        sprintf (errmsg, "Number of variables in the input netCDF4 file is "
            "greater than expected.  This software currently supports %d "
            "variables, but there are %d variables.", MAX_GRIDDED_NVARS, nvars);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get information about each of the dimensions */
    for (i = 0; i < ndims; i++)
    {
        if ((retval = nc_inq_dim (ncid, i, in_dimnames[i],
             &in_dimsizes[i])))
        {
            nc_strerror (retval);
            sprintf (errmsg, "Error inquiring about dimension %d", i);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Get information about the variables in the file */
    attr->primary_index = -1;     /* initialize index to invalid value */
    for (i = 0; i < nvars; i++)
    {
        if ((retval = nc_inq_var (ncid, i, in_varnames[i],
             &in_data_type[i], &in_var_ndims[i], in_var_dimids[i],
             &in_var_natts[i])))
        {
            nc_strerror (retval);
            sprintf (errmsg, "Error inquiring about variable %d", i);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* If this variable is our primary variable to be read then store
           the index and break out */
        if (!strcmp (in_varnames[i], varname))
        {
            attr->primary_index = i;
            primary_ndims = in_var_ndims[i];
            for (d = 0; d < primary_ndims; d++)
            {
                primary_dimids[d] = in_var_dimids[i][d];
            }
            break;
        }
    }

    /* Make sure the primary variable was found */
    if (attr->primary_index == -1)
    {
        sprintf (errmsg, "Primary variable %s was not found in the netCDF "
            "dataset.", varname);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Store the native data type of the primary variable */
    attr->native_data_type = in_data_type[attr->primary_index];

    /* Store the number of lines and samples in the desired primary variable */
    if (primary_ndims == 2)
    {  /* set up 2D dimensions */
        attr->nsamps = in_dimsizes[primary_dimids[1]];
        attr->nlines = in_dimsizes[primary_dimids[0]];
    }
    else
    {
        sprintf (errmsg, "Only 2D variables are supported.  %dD "
            "variable is not supported at this time.", primary_ndims);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get the pixel size for the specified primary variable */
    attr->pixel_size_defined = true;
    if ((retval = nc_get_att_text (ncid, attr->primary_index,
         "resolution", tmpstr)))
    {
        attr->pixel_size_defined = false;
    }
    else
    {
        /* Parse the text string to pull out the x and y pixel sizes
           Ex. - resolution = y: 0.000014 rad x: 0.000014 rad */
        sscanf (tmpstr, "%s %lf %s %s %lf %s", ignore_str, &attr->pixel_size[1],
        ignore_str, ignore_str, &attr->pixel_size[0], ignore_str);
    }

    /* Get the scale factor for the specified primary variable */
    attr->scale_defined = true;
    if ((retval = nc_get_att_double (ncid, attr->primary_index,
         "scale_factor", &attr->scale_fact)))
    {
        attr->scale_defined = false;
    }

    /* Get the offset for the specified primary variable */
    attr->offset_defined = true;
    if ((retval = nc_get_att_double (ncid, attr->primary_index,
         "add_offset", &attr->add_offset)))
    {
        attr->offset_defined = false;
    }


    /* Set up the void pointer for the fill value and read the valid range */
    attr->fill_defined = true;
    attr->fill_value = -99999.0;
    switch (in_data_type[attr->primary_index])
    {
        case NC_BYTE:
            fill_value = &fill_char;

            /* Get the valid range for the specified primary variable.  Valid
               range should have two values for this attribute. */
            if ((retval = nc_get_att_schar (ncid, attr->primary_index,
                 "valid_range", valid_range_char)))
            {
                attr->valid_range_defined = false;
            }
            else
            {
                attr->valid_range_defined = true;
                attr->valid_range[0] = valid_range_char[0];
                attr->valid_range[1] = valid_range_char[1];
            }
            break;

        case NC_UBYTE:
            fill_value = &fill_ubyte;

            /* Get the valid range for the specified primary variable.  Valid
               range should have two values for this attribute. */
            if ((retval = nc_get_att_ubyte (ncid, attr->primary_index,
                 "valid_range", valid_range_ubyte)))
            {
                attr->valid_range_defined = false;
            }
            else
            {
                attr->valid_range_defined = true;
                attr->valid_range[0] = valid_range_ubyte[0];
                attr->valid_range[1] = valid_range_ubyte[1];
            }
            break;

        case NC_SHORT:
            fill_value = &fill_short;

            /* Get the valid range for the specified primary variable.  Valid
               range should have two values for this attribute. */
            if ((retval = nc_get_att_short (ncid, attr->primary_index,
                 "valid_range", valid_range_short)))
            {
                attr->valid_range_defined = false;
            }
            else
            {
                attr->valid_range_defined = true;
                attr->valid_range[0] = valid_range_short[0];
                attr->valid_range[1] = valid_range_short[1];
            }
            break;

        case NC_USHORT:
            fill_value = &fill_ushort;

            /* Get the valid range for the specified primary variable.  Valid
               range should have two values for this attribute. */
            if ((retval = nc_get_att_ushort (ncid, attr->primary_index,
                 "valid_range", valid_range_ushort)))
            {
                attr->valid_range_defined = false;
            }
            else
            {
                attr->valid_range_defined = true;
                attr->valid_range[0] = valid_range_ushort[0];
                attr->valid_range[1] = valid_range_ushort[1];
            }
            break;

        case NC_INT:
            fill_value = &fill_int;

            /* Get the valid range for the specified primary variable.  Valid
               range should have two values for this attribute. */
            if ((retval = nc_get_att_int (ncid, attr->primary_index,
                 "valid_range", valid_range_int)))
            {
                attr->valid_range_defined = false;
            }
            else
            {
                attr->valid_range_defined = true;
                attr->valid_range[0] = valid_range_int[0];
                attr->valid_range[1] = valid_range_int[1];
            }
            break;

        case NC_UINT:
            fill_value = &fill_uint;

            /* Get the valid range for the specified primary variable.  Valid
               range should have two values for this attribute. */
            if ((retval = nc_get_att_uint (ncid, attr->primary_index,
                 "valid_range", valid_range_uint)))
            {
                attr->valid_range_defined = false;
            }
            else
            {
                attr->valid_range_defined = true;
                attr->valid_range[0] = valid_range_uint[0];
                attr->valid_range[1] = valid_range_uint[1];
            }
            break;

        case NC_FLOAT:
            fill_value = &fill_float;

            /* Get the valid range for the specified primary variable.  Valid
               range should have two values for this attribute. */
            if ((retval = nc_get_att_float (ncid, attr->primary_index,
                 "valid_range", valid_range_float)))
            {
                attr->valid_range_defined = false;
            }
            else
            {
                attr->valid_range_defined = true;
                attr->valid_range[0] = valid_range_float[0];
                attr->valid_range[1] = valid_range_float[1];
            }
            break;

        default:
            sprintf (errmsg, "Only byte, ubyte, ushort, short, uint, int, and "
                "float data types are supported at this time (%d).",
                in_data_type[attr->primary_index]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
    }

    /* Get the fill value, which is dependent upon the data type */
    if ((retval = nc_inq_var_fill (ncid, attr->primary_index, &no_fill,
        (void *) fill_value)) || (no_fill == 1))
    {
        sprintf (errmsg, "No _FillValue attribute value specified for the "
            "%s variable.", varname);
        error_handler (false, FUNC_NAME, errmsg);
        attr->fill_defined = false;
    }
    else
    {
        /* Assign the fill value */
        switch (in_data_type[attr->primary_index])
        {
            case NC_BYTE:
                attr->fill_value = fill_char;
                break;

            case NC_UBYTE:
                attr->fill_value = fill_ubyte;
                break;

            case NC_SHORT:
                attr->fill_value = fill_short;
                break;

            case NC_USHORT:
                attr->fill_value = fill_ushort;
                break;

            case NC_INT:
                attr->fill_value = fill_int;
                break;

            case NC_UINT:
                attr->fill_value = fill_uint;
                break;

            case NC_FLOAT:
                attr->fill_value = fill_float;
                break;

            default:
                sprintf (errmsg, "Only byte, ubyte, ushort, short, uint, int, "
                    "and float data types are supported at this time (%d).",
                    in_data_type[attr->primary_index]);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
        }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  close_netcdf

PURPOSE: Close the specified netCDF4 file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error closing the netCDF file
SUCCESS         Successfully closed the netCDF file

NOTES:
1. It is assumed the netCDF file (ncid) has been opened.
******************************************************************************/
int close_netcdf
(
    int ncid           /* I: netCDF file ID */
)
{
    int retval;                         /* function call return value */
    char errmsg[STR_SIZE];              /* error message */
    char FUNC_NAME[] = "close_netcdf";  /* function name */

    /* Close the netCDF file. This frees up any internal netCDF resources
       associated with the file and flushes any buffers. */
    retval = nc_close (ncid);
    if (retval)
    {
        nc_strerror (retval);
        sprintf (errmsg, "Error closing netCDF file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  ncdf_read_gridded_var

PURPOSE: Generic function to read the variable from a gridded netCDF4 file.
This function will support 2D datasets.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reading from the input netCDF files
SUCCESS         Successfully read the provided netCDF file

NOTES:
1. This function supports [u]byte, [u]short, [u]int, and float data types.
   If other data types are needed, then support will need to be added.
2. The input variable data stored as 2D products are stored as [YDIM][XDIM].
3. Memory is allocated for the nc_data array.  It will be up to the calling
   routine to free that memory.
******************************************************************************/
int ncdf_read_gridded_var
(
    int ncid,            /* I: netCDF file ID */
    int primary_index,   /* I: index of the primary variable */
    int nlines,          /* I: number of lines in y-dimension */
    int nsamps,          /* I: number of samples in x-dimension */
    nc_type native_data_type, /* I: native data type of the primary variable */
    void **nc_data      /* O: address of 2D array of netcdf data values read */
)
{
    char errmsg[STR_SIZE];                       /* error message */
    char FUNC_NAME[] = "ncdf_read_gridded_var";  /* function name */
    int retval;                       /* function call return value */
    size_t start[MAX_GRIDDED_NDIMS];  /* starting location for reading data */
    size_t count[MAX_GRIDDED_NDIMS];  /* how many values within each dimension
                                         will be read */

    /* Initialize start variable to start reading at 0,0 index */
    start[0] = 0;         /* y */
    start[1] = 0;         /* x */

    /* Read the entire variable */
    count[0] = nlines;    /* y */
    count[1] = nsamps;    /* x */

    /* Allocate space for reading the 2D data from the primary variable. Then
       read the data. */
    switch (native_data_type)
    {
        case NC_BYTE:
            /* Allocate signed bytes */
            *nc_data = (signed char *) calloc ((nlines * nsamps),
                sizeof (signed char));
            if (*nc_data == NULL)
            {
                sprintf (errmsg, "Error allocating memory for netCDF byte "
                    "array.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Read the data */
            if ((retval = nc_get_vara_schar (ncid, primary_index, start, count,
                 *nc_data)))
            {
                nc_strerror (retval);
                sprintf (errmsg, "Error reading data from primary variable");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            break;

        case NC_UBYTE:
            /* Allocate unsigned bytes */
            *nc_data = (ubyte *) calloc ((nlines * nsamps), sizeof (ubyte));
            if (*nc_data == NULL)
            {
                sprintf (errmsg, "Error allocating memory for netCDF ubyte "
                    "array.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Read the data */
            if ((retval = nc_get_vara_ubyte (ncid, primary_index, start, count,
                 *nc_data)))
            {
                nc_strerror (retval);
                sprintf (errmsg, "Error reading data from primary variable");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            break;

        case NC_USHORT:
            /* Allocate unsigned shorts */
            *nc_data = (ushort *) calloc ((nlines * nsamps), sizeof (ushort));
            if (*nc_data == NULL)
            {
                sprintf (errmsg, "Error allocating memory for netCDF ushort "
                    "array.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Read the data */
            if ((retval = nc_get_vara_ushort (ncid, primary_index, start, count,
                 *nc_data)))
            {
                nc_strerror (retval);
                sprintf (errmsg, "Error reading data from primary variable");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            break;

        case NC_SHORT:
            /* Allocate shorts */
            *nc_data = (short *) calloc ((nlines * nsamps), sizeof (short));
            if (*nc_data == NULL)
            {
                sprintf (errmsg, "Error allocating memory for netCDF short "
                    "array.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Read the data */
            if ((retval = nc_get_vara_short (ncid, primary_index, start, count,
                 *nc_data)))
            {
                nc_strerror (retval);
                sprintf (errmsg, "Error reading data from primary variable");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            break;

        case NC_UINT:
            /* Allocate unsigned ints */
            *nc_data = (uint *) calloc ((nlines * nsamps), sizeof (uint));
            if (*nc_data == NULL)
            {
                sprintf (errmsg, "Error allocating memory for netCDF uint "
                    "array.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Read the data */
            if ((retval = nc_get_vara_uint (ncid, primary_index, start, count,
                 *nc_data)))
            {
                nc_strerror (retval);
                sprintf (errmsg, "Error reading data from primary variable");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            break;

        case NC_INT:
            /* Allocate ints */
            *nc_data = (int *) calloc ((nlines * nsamps), sizeof (int));
            if (*nc_data == NULL)
            {
                sprintf (errmsg, "Error allocating memory for netCDF int "
                    "array.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Read the data */
            if ((retval = nc_get_vara_int (ncid, primary_index, start, count,
                 *nc_data)))
            {
                nc_strerror (retval);
                sprintf (errmsg, "Error reading data from primary variable");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            break;

        case NC_FLOAT:
            /* Allocate floats */
            *nc_data = (float *) calloc ((nlines * nsamps), sizeof (float));
            if (*nc_data == NULL)
            {
                sprintf (errmsg, "Error allocating memory for netCDF float "
                    "array.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Read the data */
            if ((retval = nc_get_vara_float (ncid, primary_index, start, count,
                 *nc_data)))
            {
                nc_strerror (retval);
                sprintf (errmsg, "Error reading data from primary variable");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            break;

        default:
            sprintf (errmsg, "Only byte, ubyte, ushort, short, uint, int, and "
                "float data types are supported at this time (%d).",
                native_data_type);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
    }

    return (SUCCESS);
}
