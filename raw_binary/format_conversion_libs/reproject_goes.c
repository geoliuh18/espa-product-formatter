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
Type = int
Value           Description
-----           -----------
ERROR           Specified x,y is not a valid lat,long
SUCCESS         Successfully converted x,y to lat,long

NOTES:
1. This conversion comes from the GOES-R ABI Product User Guide, L2, v5.
   Section 4.2.8.1 provides the equations for converting x, y coordinates to
   geodetic lat, long.
******************************************************************************/
int goes_xy_to_latlon
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
    double tmp_val;  /* used for storing intermediate results */
    double sin_x, sin_y;     /* sine of x/y */
    double cos_x, cos_y;     /* cosine of x/y */
    double sin2_x, sin2_y;   /* sine of x/y, squared */
    double cos2_x, cos2_y;   /* cosine of x/y, squared */
    double H_minus_S_x;      /* H - Sx */
    double H2_minus_S_x;     /* (H - Sx) squared */
    double lat_rad, lon_rad; /* latitude/longitude in radians */

printf ("DEBUG: x,y: %f, %f\n", x, y);
    /* Pull the variables from the XML metadata projection info */
    r_pol = proj_info->semi_minor_axis;
    r_eq = proj_info->semi_major_axis;
    lambda_0 = proj_info->central_meridian * RAD;  /* convert to radians */
    H = proj_info->satellite_height + proj_info->semi_major_axis;
printf ("DEBUG: r_pol = %f\n", r_pol);
printf ("DEBUG: r_eq = %f\n", r_eq);
printf ("DEBUG: lambda_0 = %f\n", lambda_0);
printf ("DEBUG: H = %f\n", H);

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

    /* Compute distance from the satellite to point x,y.  If the intermediate
       value of b^2 * 4ac is negative, then the lat/long cannot be computed
       and is not valid. */
    tmp_val = b * b - 4.0 * a * c;
    if (tmp_val < 0.0)
        return (ERROR);
    r_s = (-b - (sqrt (tmp_val))) / (2.0 * a);
printf ("DEBUG: result = %lf\n", tmp_val);
printf ("DEBUG: r_s = %f\n", r_s);

    /* Compute coefficients S_x, S_y, S_z */
    S_x = r_s * cos_x * cos_y;
    S_y = r_s * sin_x;
    S_z = r_s * cos_x * sin_y;
    S2_y = S_y * S_y;
    H_minus_S_x = H - S_x;
    H2_minus_S_x = H_minus_S_x * H_minus_S_x;
printf ("DEBUG: S_x = %f\n", S_x);
printf ("DEBUG: S_y = %f\n", S_y);
printf ("DEBUG: S_z = %f\n", S_z);

    /* Compute the latitude/longitude for point x,y in radians and convert to
       degrees */
    lat_rad = atan ((r2_eq / r2_pol) * (S_z / sqrt (H2_minus_S_x + S2_y)));
    *lat = lat_rad * DEG;
    lon_rad = lambda_0 + atan (S_y / H_minus_S_x);
    *lon = lon_rad * DEG;

    return (SUCCESS);
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
printf ("DEBUG: UL corner lat,long: %f, %f\n", gmeta->ul_corner[0], gmeta->ul_corner[1]);
printf ("DEBUG: LR corner lat,long: %f, %f\n", gmeta->lr_corner[0], gmeta->lr_corner[1]);
printf ("nlines, nsamps: %d, %d\n", bmeta->nlines, bmeta->nsamps);
    *pixsize_y = (gmeta->ul_corner[0] - gmeta->lr_corner[0]) / bmeta->nlines;
    *pixsize_x = (gmeta->lr_corner[1] - gmeta->ul_corner[1]) / bmeta->nsamps;
}


/******************************************************************************
MODULE:  reproject_goes

PURPOSE: Reproject the GOES data to the geographic lat/long projection, using
native resolutions in degrees.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error reprojecting the GOES-R ABI band
SUCCESS         Successfully reprojected GOES-R ABI band

NOTES:
1. It is assumed the input goes_bmeta and goes_gmeta structures are populated
   from reading the input GOES data. The projection and lat/long coordinates
   should reflect the overall bounding lat/long values.
2. The output bmeta pixel size and units will be updated to convert the GOES
   radian resolution into degrees.
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
    int tmp_percent;          /* current percentage for printing status */
    int curr_tmp_percent;     /* percentage for current line */
    int line, samp;           /* looping variable for line, sample of image */
    long pix;                 /* looping variable for the pixels in image */
    int x_loc, y_loc;         /* x,y location in the input buffer */
    double pixsize_x;         /* output pixel size for the longitude (degs) */
    double pixsize_y;         /* output pixel size for the latitude (degs) */
    double lat, lon;          /* lat,long in reprojected image */
    double x, y;              /* x,y coords for lat,long in GOES projection */
    int8_t *in_file_buf_int8 = NULL;  /* int8 pointer to the input buffer */
    int8_t *out_file_buf_int8 = NULL; /* int8 pointer to the output buffer */
    int16_t *in_file_buf_int16 = NULL;  /* int16 pointer to the input buffer */
    int16_t *out_file_buf_int16 = NULL; /* int16 pointer to the output buffer */

    /* Assign the output pixel size in degrees. If this is a 14 microradian
       band, then the output resolution will be 0.00449 degrees for both x&y
       directions. If this is a 28 microradian band, then the output resolution
       will be 0.00898 degrees for both x&y directions. These map to the
       documented 500m and 1000m resolutions for these bands, respectively.
       The conversion uses the fact that 1 deg = 111.325 km, on average.
       Use the x resolution for comparison, given that both the x and y
       resolutions are the same for the GOES geostationary bands. */
    if (fabs (goes_bmeta->pixel_size[0] - 0.000014) < 0.0000009)
        pixsize_x = pixsize_y = 0.00449;  /* 0.5 km */
    else if (fabs (goes_bmeta->pixel_size[0] - 0.000028) < 0.0000009)
        pixsize_x = pixsize_y = 0.00898;  /* 1 km */
    else
    {  /* not set up to support anything that isn't 14 or 28 microradian */
        sprintf (errmsg, "Unexpected resolution of the GOES band %f.  "
            "Currently only 14 and 28 microradian bands are supported.",
            goes_bmeta->pixel_size[0]);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    printf ("DEBUG: Output pixel size for geographic: %lf x %lf\n", pixsize_x, pixsize_y);

    /* Determine how many lines and samples are required at the current
       resolution to cover the lat/long boundaries. If there are any floating
       point values, round up to the next line/sample. */
    bmeta->nsamps = ceil ((gmeta->proj_info.lr_corner[0] -
        gmeta->proj_info.ul_corner[0]) / pixsize_x);
    bmeta->nlines = ceil ((gmeta->proj_info.ul_corner[1] -
        gmeta->proj_info.lr_corner[1]) / pixsize_y);
    printf ("DEBUG: Output nlines, nsamps for geographic: %d x %d\n", bmeta->nlines, bmeta->nsamps);

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
    gmeta->bounding_coords[ESPA_EAST] = gmeta->lr_corner[1];  /* long */
    gmeta->bounding_coords[ESPA_SOUTH] = gmeta->lr_corner[0]; /* lat */

    /* Adjust the projection corners for the center of the pixel */
    gmeta->proj_info.ul_corner[0] += pixsize_x * 0.5;
    gmeta->proj_info.ul_corner[1] -= pixsize_y * 0.5;
    gmeta->proj_info.lr_corner[0] -= pixsize_x * 0.5;
    gmeta->proj_info.lr_corner[1] += pixsize_y * 0.5;

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
    tmp_percent = 0;
    for (line = 0; line < bmeta->nlines; line++, lat -= pixsize_y)
    {
        /* update status */
        curr_tmp_percent = 100 * line / bmeta->nlines;
        if (curr_tmp_percent > tmp_percent)
        {
            tmp_percent = curr_tmp_percent;
            if (tmp_percent % 10 == 0)
            {
                printf ("%d%% ", tmp_percent);
                fflush (stdout);
            }
        }

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
                            in_file_buf_int8[y_loc*goes_bmeta->nsamps + x_loc];
                    else if (bmeta->data_type == ESPA_INT16)
                        out_file_buf_int16[pix] =
                            in_file_buf_int16[y_loc*goes_bmeta->nsamps + x_loc];
                }
            }
        }
    }

    /* Successful conversion */
    printf ("100%%\n");
    fflush (stdout);

    return (SUCCESS);
}

