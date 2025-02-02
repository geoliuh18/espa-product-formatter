/*****************************************************************************
FILE: convert_sentinel_to_espa.c
  
PURPOSE: Contains functions for reading Sentinel-2 1C products and writing to
the ESPA raw binary format.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
  1. The XML metadata format written via this library follows the ESPA internal
     metadata format found in ESPA Raw Binary Format v1.0.doc.  The schema for
     the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/schema/espa_internal_metadata_v1_x.xsd.
*****************************************************************************/
#include <unistd.h>
#include <ctype.h>
#include "convert_sentinel_to_espa.h"

/* Band information for the Sentinel-2 L1C products */
char sentinel_bands[NUM_SENTINEL_BANDS][STR_SIZE] =
    {"B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A",
     "B09", "B10", "B11", "B12", "TCI"};
char sentinel_band_nums[NUM_SENTINEL_BANDS][STR_SIZE] =
    {"1", "2", "3", "4", "5", "6", "7", "8", "8A", "9", "10", "11", "12",
     "TCI"};  /* TCI is true colour image */

/******************************************************************************
MODULE:  rename_jp2

PURPOSE: Rename the Sentinel JPEG2000 files from the shortened JP2 filename to
the more informative granule name {product_id}_{band}.jp2.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error renaming the JP2 files
SUCCESS         Successfully renamed the JP2 files

NOTES:
******************************************************************************/
int rename_jp2
(
    Espa_internal_meta_t *xml_metadata /* I: valid ESPA metadata structure */
)
{
    char FUNC_NAME[] = "rename_jp2";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    char newfile[STR_SIZE];   /* name of the new Sentinel file */
    int i;                    /* looping variable for bands in XML file */
    int count;                /* number of chars copied in snprintf */
    Espa_band_meta_t *bmeta = NULL;  /* pointer to band metadata */
    Espa_global_meta_t *gmeta = &xml_metadata->global;  /* global metadata */

    /* Loop through the bands in the metadata file and convert each one to
       the ESPA format */
    for (i = 0; i < xml_metadata->nbands; i++)
    {
        /* Set up the band metadata pointer */
        bmeta = &xml_metadata->band[i];

        /* Rename the current JP2 filename to {product_id}_{bandname}.jp2 */
        sprintf (newfile, "%s_%s.jp2", gmeta->product_id, sentinel_bands[i]);
        if (rename (bmeta->file_name, newfile))
        {
            sprintf (errmsg, "Unable to rename the original Sentinel JP2 "
                "file (%s) to the new ESPA filename (%s)", bmeta->file_name,
                newfile);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Write the new filename to the ESPA metadata */
        count = snprintf (bmeta->file_name, sizeof (bmeta->file_name), "%s",
            newfile);
        if (count < 0 || count >= sizeof (bmeta->file_name))
        {
            sprintf (errmsg, "Overflow of bmeta->file_name string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Successful rename */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  convert_jp2_to_img

PURPOSE: Convert the Sentinel JP2 bands to an ESPA raw binary (.img) file,
and generate the ENVI header file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error converting the JP2 file
SUCCESS         Successfully converted JP2 file

NOTES:
******************************************************************************/
int convert_jp2_to_img
(
    Espa_internal_meta_t *xml_metadata /* I: valid ESPA metadata structure */
)
{
    char FUNC_NAME[] = "convert_jp2_to_img";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    char *cptr = NULL;        /* pointer to the file extension */
    char envi_file[STR_SIZE]; /* name of the output ENVI header file */
    char raw_file[STR_SIZE];  /* name of the output raw binary file (.raw) */
    char jp2_cmd[STR_SIZE];   /* command string for opj_decompress */
    int i;                    /* looping variable for bands in XML file */
    int count;                /* number of chars copied in snprintf */
    Envi_header_t envi_hdr;   /* output ENVI header information */
    Espa_band_meta_t *bmeta = NULL;  /* pointer to band metadata */
    Espa_global_meta_t *gmeta = &xml_metadata->global;  /* global metadata */

    /* Setup the opj_decompress command for converting all the bands in the
       current directory from JP2 to img.  This does not create an ENVI header
       file for the bands. */
    strcpy (jp2_cmd, "opj_decompress -ImgDir . -OutFor RAW -quiet");
    if (system (jp2_cmd) == -1)
    {
        sprintf (errmsg, "Decompressing JP2 files: %s. Make sure the current "
            "directory is writable and the openjpeg opj_decompress tool is in "
            "your system PATH", jp2_cmd);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through the bands in the metadata file and convert each one to
       the ESPA format */
    for (i = 0; i < xml_metadata->nbands; i++)
    {
        /* Set up the band metadata pointer */
        bmeta = &xml_metadata->band[i];

        /* Determine the name of the output raw binary file.  Replace the
           jp2 file extension with img in the Sentinel filenames. */
        cptr = strrchr (bmeta->file_name, '.');
        if (cptr == NULL)
        {
            sprintf (errmsg, "No file extension found in the Sentinel JP2 "
                "file: %s\n", bmeta->file_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        strcpy (cptr, ".img");

        /* Rename the .raw files from the opg_decompress to .img files */
        count = snprintf (raw_file, sizeof (raw_file), "%s", bmeta->file_name);
        if (count < 0 || count >= sizeof (raw_file))
        {
            sprintf (errmsg, "Overflow of raw_file string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        cptr = strrchr (raw_file, '.');
        strcpy (cptr, ".raw");

        if (rename (raw_file, bmeta->file_name))
        {
            sprintf (errmsg, "Unable to rename the decompressed Sentinel raw "
                "file (%s) to the new ESPA filename (%s)", raw_file,
                bmeta->file_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Create the ENVI header file for this band */
        if (create_envi_struct (bmeta, gmeta, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating the ENVI header structure for this "
                "file: %s", bmeta->file_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Write the ENVI header */
        count = snprintf (envi_file, sizeof (envi_file), "%s",
            bmeta->file_name);
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
    }  /* end for */

    /* Successful conversion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  convert_sentinel_to_espa

PURPOSE: Converts the input Sentinel-2 (A&B L1C) files to the ESPA internal raw
binary file format (and associated XML file).  The MTD_MSIL1C.xml and MTD_TL.xml
files are expected to be in the same directory as the Sentinel band data.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error converting the Sentinel-2 product
SUCCESS         Successfully converted Sentinel-2 product to ESPA format

NOTES:
  1. The Sentinel JP2 band files will be deciphered from the Sentinel XML file.
  2. The JP2 band files listed in the MTD_MSIL1C.xml file need to be available
     in the same directory as both the MTD_MSIL1C product and MTD_TL tile XML
     files.
******************************************************************************/
int convert_sentinel_to_espa
(
    bool del_src      /* I: should the source .jp2 files be removed after
                            conversion? */
)
{
    char FUNC_NAME[] = "convert_sentinel_to_espa";  /* function name */
    char errmsg[STR_SIZE];            /* error message */
    char espa_xml_file[STR_SIZE];     /* output ESPA XML metadata filename */
    char jp2_file[STR_SIZE];          /* jp2 image file to delete */
    char sentinel_xml_file[STR_SIZE]; /* current Sentinel XML filename */
    char *cptr = NULL;                /* pointer to the file extension */
    int i;                            /* looping variable */
    int count;                        /* number of chars copied in snprintf */
    Espa_internal_meta_t xml_metadata;  /* ESPA XML metadata structure to be
                                           populated by reading the Sentinel
                                           XML file */
    Espa_global_meta_t *gmeta;  /* global metadata structure */
    Espa_band_meta_t *bmeta;    /* band metadata pointer to all bands */

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);
    gmeta = &xml_metadata.global;

    /* Read the Sentinel MTD_MSIL1C product XML file and populate our internal
       ESPA metadata structure. The names of the image files/bands, the
       acquisition date/time, product generation date/time, lat/long coords,
       product type, and scale factor are all available in this XML file. */
    strcpy (sentinel_xml_file, "MTD_MSIL1C.xml");
    if (parse_sentinel_product_metadata (sentinel_xml_file, &xml_metadata) !=
        SUCCESS)
    {
        sprintf (errmsg, "Reading Sentinel product XML file: %s",
            sentinel_xml_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Make sure the number of bands that were found in the XML file matches
       the expected number of Sentinel L1C bands */
    if (xml_metadata.nbands != NUM_SENTINEL_BANDS)
    {
        sprintf (errmsg, "Number of bands read from %s (%d) does not match "
            "the expected number of bands %d", sentinel_xml_file,
            xml_metadata.nbands, NUM_SENTINEL_BANDS);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the Sentinel MTD_TL tile XML file and populate our internal ESPA
       metadata structure. The tile level datum, projection, and zone are
       available. The number of lines/samples for each resolution are
       available. The UL x/y position are also availble. */
    strcpy (sentinel_xml_file, "MTD_TL.xml");
    if (parse_sentinel_tile_metadata (sentinel_xml_file, &xml_metadata) !=
        SUCCESS)
    {
        sprintf (errmsg, "Reading Sentinel tile XML file: %s",
            sentinel_xml_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Add the data provider USGS/EROS (Sentinel products delivered via EE) */
    strcpy (gmeta->data_provider, "USGS/EROS");

    /* Add the instrument which is MSI - MultiSpectral Instrument */
    strcpy (gmeta->instrument, "MSI");

    /* Set the orientation angle to 0 */
    gmeta->orientation_angle = 0.0;

    /* The last band in the product is a 3-band browse. That band will be
       skipped from processing since it is not true science data. */
    xml_metadata.nbands--;

    /* Update remaining information for band metadata for each of the bands */
    for (i = 0; i < xml_metadata.nbands; i++)
    {
        bmeta = &xml_metadata.band[i];
        strcpy (bmeta->product, "MSIL1C");
        strcpy (bmeta->name, sentinel_bands[i]);
        strcpy (bmeta->category, "image");
        bmeta->data_type = ESPA_UINT16;
        bmeta->fill_value = 0;
        bmeta->saturate_value = 65535;
        bmeta->valid_range[0] = 0.0;
        bmeta->valid_range[1] = 65534.0;
        strcpy (bmeta->data_units, "reflectance");
        strcpy (bmeta->production_date, gmeta->level1_production_date);
        sprintf (bmeta->long_name, "band %s top-of-atmosphere reflectance",
            sentinel_band_nums[i]);

        /* Sentinel XML files don't indicate the application used to process
           the original image files, so set to "not available" */
        strcpy (bmeta->app_version, "not available");
    }

    /* Rename the current Sentinel JP2 bands to a new filename (using the
       product_id) to be used by ESPA */
    if (rename_jp2 (&xml_metadata) != SUCCESS)
    {
        sprintf (errmsg, "Renaming Sentinel JP2 image files");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Convert each of the Sentinel JP2 bands to raw binary, also create the
       ENVI header files using the XML metadata. Updates the filenames for
       each band to raw binary. */
    if (convert_jp2_to_img (&xml_metadata) != SUCCESS)
    {
        sprintf (errmsg, "Converting JP2 bands to raw binary");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the metadata from our internal metadata structure to the output
       XML filename */
    sprintf (espa_xml_file, "%s.xml", gmeta->product_id);
    if (write_metadata (&xml_metadata, espa_xml_file) != SUCCESS)
    {  /* Error messages already written */
        return (ERROR);
    }

    /* Validate the output metadata file */
    if (validate_xml_file (espa_xml_file) != SUCCESS)
    {  /* Error messages already written */
        return (ERROR);
    }

    /* Remove the source JP2 files if specified */
    if (del_src)
    {
        /* Remove the image band */
        for (i = 0; i < xml_metadata.nbands; i++)
        {
            bmeta = &xml_metadata.band[i];

            /* Remove the .jp2 files */
            count = snprintf (jp2_file, sizeof (jp2_file), "%s",
                bmeta->file_name);
            if (count < 0 || count >= sizeof (jp2_file))
            {
                sprintf (errmsg, "Overflow of jp2_file string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            cptr = strrchr (jp2_file, '.');
            strcpy (cptr, ".jp2");

            /* Remove the source file */
            printf ("  Removing %s\n", jp2_file);
            if (unlink (jp2_file) != 0)
            {
                sprintf (errmsg, "Deleting source file: %s", jp2_file);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* Remove the TCI jp2 band */
        bmeta = &xml_metadata.band[NUM_SENTINEL_BANDS-1];
        printf ("  Removing %s\n", bmeta->file_name);
        if (unlink (bmeta->file_name) != 0)
        {
            sprintf (errmsg, "Deleting source file: %s", bmeta->file_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Remove the TCI raw band */
        count = snprintf (jp2_file, sizeof (jp2_file), "%s", bmeta->file_name);
        if (count < 0 || count >= sizeof (jp2_file))
        {
            sprintf (errmsg, "Overflow of jp2_file string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        cptr = strrchr (jp2_file, '.');
        strcpy (cptr, ".raw");

        /* Remove the source file */
        printf ("  Removing %s\n", jp2_file);
        if (unlink (jp2_file) != 0)
        {
            sprintf (errmsg, "Deleting source file: %s", jp2_file);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Successful conversion */
    return (SUCCESS);
}

