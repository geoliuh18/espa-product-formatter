/*****************************************************************************
FILE: convert_goes_to_espa.h
  
PURPOSE: Contains defines and prototypes to read supported GOES-R ABI files,
create the XML metadata file, and convert from netCDF to raw binary file format.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/

#ifndef CONVERT_GOES_TO_ESPA_H
#define CONVERT_GOES_TO_ESPA_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "error_handler.h"
#include "espa_metadata.h"
#include "espa_geoloc.h"
#include "netcdf_io.h"
#include "raw_binary_io.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "reproject_goes.h"


/* Defines */
enum Espa_goes_band_types
{
    GOES_CMI=0, GOES_DQF, NGOES_BANDS
};

/* Prototypes */
int convert_goes_to_espa
(
    char *goes_netcdf_file, /* I: input GOES-R ABI netCDF filename */
    char *espa_xml_file,    /* I: output ESPA XML metadata filename */
    bool append_bands,      /* I: are the bands in this GOES-R file being
                                  appended?  If not, then the XML file will
                                  be created from scratch. */
    bool del_src            /* I: should the source .nc files be removed after
                                  conversion? */
);

#endif
