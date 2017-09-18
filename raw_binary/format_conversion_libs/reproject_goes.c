/*****************************************************************************
FILE: reproject_goes.c
  
PURPOSE: Contains functions for reprojecting GOES-R ABI geostationary products.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include <unistd.h>
#include <ctype.h>
#include "reproject_goes.h"


/******************************************************************************
MODULE:  goes_xy_to_latlon

PURPOSE: Converts the GOES-R ABI x/y coordinates to lat/long

RETURN VALUE:
Type = N/A

NOTES:
1. This conversion comes from the GOES-R ABI Product User Guide, L2, v5.
   Section 4.2.8.1 provides the equations for converting x, y coordinates to
   geodetic lat, long.
******************************************************************************/
void goes_xy_to_latlon
(
    double x,      /* I: GOES x coordinate - E/W scan angle (radians) */
    double y,      /* I: GOES y coordinate - N/S scan angle (radians) */
    Espa_proj_meta_t *proj_info, /* I: XML projection information */
    double *lat,    /* O: corresponding latitude value (degrees) */
    double *lon     /* O: corresponding longitude value (degrees) */
)
{
    double r_eq;     /* radius - use the GOES semi-major axis */
    double r_pol;    /* radius - use the GOES semi-minor axis */
    double r2_eq;    /* r_eq squared */
    double r2_pol;   /* r_pol squared */
    double lambda_0; /* use GOES longitude of projection origin (radians) */
    double H;        /* GOES perspective point height + GOES semi-major axis */
    double S_x, S_y, S_z;  /* coefficients to be used Sx, Sy, Sz */
    double S2_y;     /* S_y squared */
    double r_s;      /* distance from the satellite to point x,y */
    double a, b, c;  /* coefficients to be used */
    double sin_x, sin_y;     /* sine of x/y */
    double cos_x, cos_y;     /* cosine of x/y */
    double sin2_x, sin2_y;   /* sine of x/y, squared */
    double cos2_x, cos2_y;   /* cosine of x/y, squared */
    double H_minus_S_x;      /* H - Sx */
    double H2_minus_S_x;     /* (H - Sx) squared */
    double lat_rad, lon_rad; /* latitude/longitude in radians */

    /* Pull the variables from the XML metadata projection info */
    r_pol = proj_info->semi_minor_axis;
    r_eq = proj_info->semi_major_axis;
    lambda_0 = proj_info->central_meridian * RAD;  /* convert to radians */
    H = proj_info->satellite_height + proj_info->semi_major_axis;

    /* Compute coefficients a, b, and c */
    sin_x = sin (x);
    cos_x = cos (x);
    sin_y = sin (y);
    cos_y = cos (y);
    sin2_x = sin_x * sin_x;
    cos2_x = cos_x * cos_x;
    sin2_y = sin_y * sin_y;
    cos2_y = cos_y * cos_y;
    r2_eq = r_eq * r_eq;
    r2_pol = r_pol * r_pol;
    a = sin2_x + cos2_x * (cos2_y + (r2_eq / r2_pol) * sin2_y);
    b = -2.0 * H * cos_x * cos_y;
    c = H * H - r2_eq;

    /* Compute distance from the satellite to point x,y */
    r_s = (-b - (sqrt (b * b - 4.0 * a * c))) / (2.0 * a);

    /* Compute coefficients S_x, S_y, S_z */
    S_x = r_s * cos_x * cos_y;
    S_y = r_s * sin_x;
    S_z = r_s * cos_x * sin_y;
    S2_y = S_y * S_y;
    H_minus_S_x = H - S_x;
    H2_minus_S_x = H_minus_S_x * H_minus_S_x;

    /* Compute the latitude/longitude for point x,y in radians and convert to
       degrees */
    lat_rad = atan ((r2_eq / r2_pol) * (S_z / sqrt (H2_minus_S_x + S2_y)));
    *lat = lat_rad * DEG;
    lon_rad = lambda_0 + atan (S_y / H_minus_S_x);
    *lon = lon_rad * DEG;
}


/******************************************************************************
MODULE:  goes_latlon_to_xy

PURPOSE: Converts the GOES-R ABI lat/long coordinates to x/y

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Specified lat/long is not visible from the satellite
SUCCESS         Successfully converted lat/long to x/y

NOTES:
1. This conversion comes from the GOES-R ABI Product User Guide, L2, v5.
   Section 4.2.8.2 provides the equations for converting geodetic lat, long
   coordinates to x, y in radians.
******************************************************************************/
int goes_latlon_to_xy
(
    double lat,   /* I: latitude value (degrees) */
    double lon,   /* I: longitude value (degrees) */
    Espa_proj_meta_t *proj_info, /* I: XML projection information */
    double *x,    /* O: corresponding x coordinate - E/W scan angle (radians) */
    double *y     /* O: corresponding y coordinate - N/S scan angle (radians) */
)
{
    double r_eq;     /* radius - use the GOES semi-major axis */
    double r_pol;    /* radius - use the GOES semi-minor axis */
    double r2_eq;    /* r_eq squared */
    double r2_pol;   /* r_pol squared */
    double lambda_0; /* use GOES longitude of projection origin (radians) */
    double H;        /* GOES perspective point height + GOES semi-major axis */
    double S_x, S_y, S_z;    /* coefficients to be used Sx, Sy, Sz */
    double S2_x, S2_y, S2_z; /* S_x, S_y, S_z values squared */
    double r_c;      /* geocentric distance to the point on the ellipsoid */
    double lat_c;    /* geocentric latitude */
    double lat_rad, lon_rad; /* latitude/longitude in radians */
    double e = 0.0818191910435;  /* constant provided from the GOES UG */
    double e2;       /* e squared */
    double cos_lat_c;  /* cosine of the geocentric latitude */
    double cos2_lat_c; /* cosine squared of the geocentric latitude */

    /* Convert the lat/long from degrees to radians */
    lat_rad = lat * RAD;
    lon_rad = lon * RAD;

    /* Pull the variables from the XML metadata projection info */
    r_pol = proj_info->semi_minor_axis;
    r_eq = proj_info->semi_major_axis;
    r2_eq = r_eq * r_eq;
    r2_pol = r_pol * r_pol;
    lambda_0 = proj_info->central_meridian * RAD;  /* convert to radians */
    H = proj_info->satellite_height + proj_info->semi_major_axis;

    /* Compute the geocentric latitude */
    lat_c = atan ((r2_pol / r2_eq) * tan (lat_rad));

    /* Compute the geocentric distance to the point on the ellipsoid */
    e2 = e * e;
    cos_lat_c = cos (lat_c);
    cos2_lat_c = cos_lat_c * cos_lat_c;
    r_c = r_pol / sqrt (1.0 - (e2 * cos2_lat_c));

    /* Compute coefficients S_x, S_y, S_z */
    S_x = H - r_c * cos (lat_c) * cos (lon_rad - lambda_0);
    S_y = -r_c * cos (lat_c) * sin (lon_rad - lambda_0);
    S_z = r_c * sin (lat_c);
    S2_x = S_x * S_x;
    S2_y = S_y * S_y;
    S2_z = S_z * S_z;

    /* Determine if the given lat/long is visible from the satellite.  If
       not return an error and the calling routine will handle as desired. */
    if (H * (H - S_x) < S2_y + ((r2_eq / r2_pol) * S2_z))
        return (ERROR);

    /* Compute the x, y in radians from the input lat, lon */
    *x = asin (-S_y / sqrt (S2_x + S2_y + S2_z));
    *y = atan (S_z / S_x);

    return (SUCCESS);
}


/******************************************************************************
MODULE:  determine_pixsize_degs

PURPOSE: Determines the pixel size in degrees for the geographic projection

RETURN VALUE:
Type = N/A

NOTES:
******************************************************************************/
void determine_pixsize_degs
(
    Espa_global_meta_t *gmeta,  /* I: pointer to global metadata */
    Espa_band_meta_t *bmeta,    /* I: pointer to band metadata */
    double *pixsize_x,          /* O: pixel size in the x direction */
    double *pixsize_y           /* O: pixel size in the y direction */
)
{
    /* Compute the pixel size in degrees */
    *pixsize_y = (gmeta->ul_corner[0] - gmeta->lr_corner[0]) / bmeta->nlines;
    *pixsize_x = (gmeta->lr_corner[1] - gmeta->ul_corner[1]) / bmeta->nsamps;
}


/******************************************************************************
MODULE:  reproject_goes

PURPOSE: Convert the GOES-R ABI netCDF band to an ESPA raw binary (.img) file
and write the associated ENVI header for this band.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error converting the GOES-R ABI band
SUCCESS         Successfully converted GOES-R ABI band to raw binary

NOTES:
1. It is assumed the input goes_bmeta and goes_gmeta structures are populated
   from reading the input GOES data.
2. The output bmeta pixel size and units will be updated based on the
   reprojection.
3. The output gmeta LR corner is updated based on the reprojected data.
******************************************************************************/
int reproject_goes
(
    Espa_global_meta_t *goes_gmeta, /* I: GOES global metadata */
    Espa_band_meta_t *goes_bmeta,   /* I: GOES band metadata */
    Espa_global_meta_t *gmeta,  /* I/O: output/reprojected global metadata */
    Espa_band_meta_t *bmeta,    /* I/O: output/reprojected band metadata */
    void *in_file_buf,          /* I: input file buffer of GOES data,
                                      nlines x nsamps */
    void **out_file_buf         /* I: pointer to output file buffer for
                                      geographic data, nlines x nsamps */
)
{
    char FUNC_NAME[] = "reproject_goes";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int line, samp;           /* looping variable for line, sample of image */
    int pix;                  /* looping variable for the pixels in image */
    int x_loc, y_loc;         /* x,y location in the input buffer */
    double pixsize_x;         /* output pixel size for the longitude (degs) */
    double pixsize_y;         /* output pixel size for the latitude (degs) */
    double lat, lon;          /* lat,long in reprojected image */
    double x, y;              /* x,y coords for lat,long in GOES projection */
    int8_t *in_file_buf_int8 = NULL;  /* int8 pointer to the input buffer */
    int8_t *out_file_buf_int8 = NULL; /* int8 pointer to the output buffer */
    int16_t *in_file_buf_int16 = NULL;  /* int16 pointer to the input buffer */
    int16_t *out_file_buf_int16 = NULL; /* int16 pointer to the output buffer */

    /* Compute the output pixel size in degrees */
    determine_pixsize_degs (goes_gmeta, goes_bmeta, &pixsize_x, &pixsize_y);
    printf ("DEBUG: Output pixel size for geographic: %lf x %lf\n", pixsize_x, pixsize_y);

    /* Determine the number of bytes per pixel for the output product */
    if (bmeta->data_type == ESPA_INT8)
    {
        *out_file_buf = (int8_t *) calloc ((bmeta->nlines * bmeta->nsamps),
            sizeof (int8_t));
        if (*out_file_buf == NULL)
        {
            sprintf (errmsg, "Error allocating memory for reprojected int8 "
                "array.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        in_file_buf_int8 = in_file_buf;
        out_file_buf_int8 = *out_file_buf;
    }
    else if (bmeta->data_type == ESPA_INT16)
    {
        *out_file_buf = (int16_t *) calloc ((bmeta->nlines * bmeta->nsamps),
            sizeof (int16_t));
        if (*out_file_buf == NULL)
        {
            sprintf (errmsg, "Error allocating memory for reprojected short "
                "array.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        in_file_buf_int16 = in_file_buf;
        out_file_buf_int16 = *out_file_buf;
    }
    else
    {
        sprintf (errmsg, "Unsupported GOES-R ABI data type %d",
            bmeta->data_type);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Determine the new LR corner for the output product */
    gmeta->lr_corner[0] = gmeta->ul_corner[0] - bmeta->nlines * pixsize_y;
    gmeta->lr_corner[1] = gmeta->ul_corner[1] + bmeta->nsamps * pixsize_x;
    gmeta->proj_info.lr_corner[0] = gmeta->lr_corner[1]; /* x --> long */
    gmeta->proj_info.lr_corner[1] = gmeta->lr_corner[0]; /* y --> lat */

    /* Update the metadata for the current band in the output XML file.  Pixel
       size and pixel units need to be updated.  Resampling type is NN. */
    bmeta->pixel_size[0] = pixsize_x;
    bmeta->pixel_size[1] = pixsize_y;
    strcpy (bmeta->pixel_units, "degrees");
    bmeta->resample_method = ESPA_NN;

    /* Loop through the lines and samples, starting at the UL corner, and
       determine the closest value in the input geostationary image to
       represent the geographic lat/long image.  This will be a NN resample
       since we also need to support data quality flags. */
    pix = 0;
    lat = gmeta->ul_corner[0];
    for (line = 0; line < bmeta->nlines; line++, lat -= pixsize_y)
    {
        lon = gmeta->ul_corner[1];
        for (samp = 0; samp < bmeta->nsamps; samp++, lon += pixsize_x, pix++)
        {
            /* Determine the input x,y value for the current output lat,lon */
            if (goes_latlon_to_xy (lat, lon, &goes_gmeta->proj_info, &x, &y) !=
                SUCCESS)
            {
                /* This pixel is not visible from the satellite, so mark it
                   as fill */
                if (bmeta->data_type == ESPA_INT8)
                    out_file_buf_int8[pix] = bmeta->fill_value;
                else if (bmeta->data_type == ESPA_INT16)
                    out_file_buf_int16[pix] = bmeta->fill_value;
            }
            else
            {
                /* Determine the location of the x, y value in the input
                   buffer.  Just using a NN approach since this supports both
                   image as well as data quality bands. */
                x_loc = (int) roundf ((x - goes_gmeta->proj_info.ul_corner[0]) /
                    goes_bmeta->pixel_size[0]);
                y_loc = (int) roundf ((goes_gmeta->proj_info.ul_corner[1] - y) /
                    goes_bmeta->pixel_size[1]);

                /* If the x,y location doesn't fall within the valid input
                   image, then flag it as fill */
                if (x_loc < 0 || y_loc < 0 ||
                    x_loc >= goes_bmeta->nsamps || y_loc >= goes_bmeta->nlines)
                {
                    if (bmeta->data_type == ESPA_INT8)
                        out_file_buf_int8[pix] = bmeta->fill_value;
                    else if (bmeta->data_type == ESPA_INT16)
                        out_file_buf_int16[pix] = bmeta->fill_value;
                }
                else
                {
                    if (bmeta->data_type == ESPA_INT8)
                        out_file_buf_int8[pix] =
                            in_file_buf_int8[y_loc*bmeta->nsamps + x_loc];
                    else if (bmeta->data_type == ESPA_INT16)
                        out_file_buf_int16[pix] =
                            in_file_buf_int16[y_loc*bmeta->nsamps + x_loc];
                }
            }
        }
    }

    /* Successful conversion */
    return (SUCCESS);
}

