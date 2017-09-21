/*****************************************************************************
FILE: reproject_goes.h
  
PURPOSE: Contains defines and prototypes to reproject GOES-R ABI data.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/

#ifndef REPROJECT_GOES_H
#define REPROJECT_GOES_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "error_handler.h"
#include "espa_metadata.h"
#include "espa_geoloc.h"


/* Prototypes */
int goes_xy_to_latlon
(
    double x,      /* I: GOES x coordinate - E/W scan angle (radians) */
    double y,      /* I: GOES y coordinate - N/S scan angle (radians) */
    Espa_proj_meta_t *proj_info, /* I: XML projection information */
    double *lat,    /* O: corresponding latitude value (degrees) */
    double *lon     /* O: corresponding longitude value (degrees) */
);

int goes_latlon_to_xy
(
    double lat,   /* I: latitude value (degrees) */
    double lon,   /* I: longitude value (degrees) */
    Espa_proj_meta_t *proj_info, /* I: XML projection information */
    double *x,    /* O: corresponding x coordinate - E/W scan angle (radians) */
    double *y     /* O: corresponding y coordinate - N/S scan angle (radians) */
);

void determine_pixsize_degs
(
    Espa_global_meta_t *gmeta,  /* I: pointer to global metadata */
    Espa_band_meta_t *bmeta,    /* I: pointer to band metadata */
    double *pixsize_x,          /* O: pixel size in the x direction */
    double *pixsize_y           /* O: pixel size in the y direction */
);

int reproject_goes
(
    Espa_global_meta_t *goes_gmeta, /* I: GOES global metadata */
    Espa_band_meta_t *goes_bmeta,   /* I: GOES band metadata */
    Espa_global_meta_t *gmeta,  /* I/O: output global metadata */
    Espa_band_meta_t *bmeta,    /* I/O: output band metadata */
    bool adjust_gmeta,          /* I: should the global metadata be adjusted */
    void *in_file_buf,          /* I: input file buffer of GOES data,
                                      nlines x nsamps */
    void **out_file_buf         /* I: pointer to output file buffer for
                                      geographic data, nlines x nsamps */
);

#endif
