/*****************************************************************************
FILE: espa_metadata.c
  
PURPOSE: Contains functions for reading/writing/appending the ESPA internal
metadata files along with inializing/freeing memory in the metadata structures.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
  1. The XML metadata format parsed or written via this library follows the
     ESPA internal metadata format found in ESPA Raw Binary Format v1.0.doc.
     The schema for the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/schema/espa_internal_metadata_v1_0.xsd.
  2. This code relies on the libxml2 library developed for the Gnome project.
*****************************************************************************/
#include <sys/stat.h>
#include "espa_metadata.h"

/******************************************************************************
MODULE:  validate_xml_file

PURPOSE:  Validates the specified XML file with the specified schema file/URL.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           XML does not validate against the specified schema
SUCCESS         XML validates

NOTES:
******************************************************************************/
int validate_xml_file
(
    char *meta_file           /* I: name of metadata file to be validated */
)
{
    char FUNC_NAME[] = "validate_xml_file";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    char *schema_file = NULL;     /* name of schema file or URL to be validated
                                     against */
    int status;                   /* return status */
    xmlDocPtr doc = NULL;         /* resulting document tree */
    xmlSchemaPtr schema = NULL;   /* pointer to the schema */
    xmlSchemaParserCtxtPtr ctxt = NULL;  /* parser context for the schema */
    xmlSchemaValidCtxtPtr valid_ctxt = NULL;  /* pointer to validate from the
                                                 schema */
    struct stat statbuf;          /* buffer for the file stat function */

    /* Get the ESPA schema environment variable which specifies the location
       of the XML schema to be used */
    schema_file = getenv ("ESPA_SCHEMA");
    if (schema_file == NULL)
    {  /* ESPA schema environment variable wasn't defined. Try the version in
          /usr/local... */
        schema_file = LOCAL_ESPA_SCHEMA;
        if (stat (schema_file, &statbuf) == -1)
        {  /* /usr/local ESPA schema file doesn't exist.  Try the version on
              the ESPA http site... */
            schema_file = ESPA_SCHEMA;
        }
    }
printf ("DEBUG: Using %s to validate schema\n", schema_file);

    /* Set up the schema parser and parse the schema file/URL */
    xmlLineNumbersDefault (1);
    ctxt = xmlSchemaNewParserCtxt (schema_file);
    xmlSchemaSetParserErrors (ctxt, (xmlSchemaValidityErrorFunc) fprintf,
        (xmlSchemaValidityWarningFunc) fprintf, stderr);
    schema = xmlSchemaParse (ctxt);

    /* Free the schema parser context */
    xmlSchemaFreeParserCtxt (ctxt);

    /* Load the XML file and parse it to the document tree */
    doc = xmlReadFile (meta_file, NULL, 0);
    if (doc == NULL)
    {
        sprintf (errmsg, "Could not parse %s", meta_file);
        error_handler (true, FUNC_NAME, errmsg);
        sprintf (errmsg, "Possible schema file not found.  ESPA_SCHEMA "
            "environment variable isn't defined.  The first default schema "
            "location of %s doesn't exist.  And the second default location of "
            "%s was used as the last default.", LOCAL_ESPA_SCHEMA, ESPA_SCHEMA);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Identify the schema file as the validation source */
    valid_ctxt = xmlSchemaNewValidCtxt (schema);
    xmlSchemaSetValidErrors (valid_ctxt, (xmlSchemaValidityErrorFunc) fprintf,
        (xmlSchemaValidityWarningFunc) fprintf, stderr);

    /* Validate the XML metadata against the schema */
    status = xmlSchemaValidateDoc (valid_ctxt, doc);
    if (status > 0)
    {
        sprintf (errmsg, "%s fails to validate", meta_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    else if (status != 0)
    {
        sprintf (errmsg, "%s validation generated an internal error",
            meta_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Free the resources and clean up the memory */
    xmlSchemaFreeValidCtxt (valid_ctxt);
    xmlFreeDoc (doc);
    if (schema != NULL)
        xmlSchemaFree (schema);
    xmlSchemaCleanupTypes();
    xmlCleanupParser();   /* cleanup the XML library */
    xmlMemoryDump();      /* for debugging */

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  init_metadata_struct

PURPOSE:  Initializes the ESPA internal metadata structure, particularly the
pointers within each sub-structure.  Assigns field values to fill to make it
easier to detect if the values were parsed from the input metadata file or
assigned by the user.

RETURN VALUE:
Type = None

NOTES:
******************************************************************************/
void init_metadata_struct
(
    Espa_internal_meta_t *internal_meta   /* I: pointer to internal metadata
                                                structure to be initialized */
)
{
    Espa_global_meta_t *gmeta = &internal_meta->global;
                                 /* pointer to the global metadata structure */

    /* Initialze the number of bands */
    internal_meta->nbands = 0;
    internal_meta->band = NULL;

    /* Initialize the global metadata values to fill for use by the write
       metadata routines */
    strcpy (gmeta->data_provider, ESPA_STRING_META_FILL);
    strcpy (gmeta->satellite, ESPA_STRING_META_FILL);
    strcpy (gmeta->instrument, ESPA_STRING_META_FILL);
    strcpy (gmeta->acquisition_date, ESPA_STRING_META_FILL);
    strcpy (gmeta->scene_center_time, ESPA_STRING_META_FILL);
    strcpy (gmeta->level1_production_date, ESPA_STRING_META_FILL);
    gmeta->solar_zenith = ESPA_FLOAT_META_FILL;
    gmeta->solar_azimuth = ESPA_FLOAT_META_FILL;
    strcpy (gmeta->solar_units, ESPA_STRING_META_FILL);
    gmeta->earth_sun_dist = ESPA_FLOAT_META_FILL;
    gmeta->wrs_system = ESPA_INT_META_FILL;
    gmeta->wrs_path = ESPA_INT_META_FILL;
    gmeta->wrs_row = ESPA_INT_META_FILL;
    gmeta->htile = ESPA_INT_META_FILL;
    gmeta->vtile = ESPA_INT_META_FILL;
    strcpy (gmeta->lpgs_metadata_file, ESPA_STRING_META_FILL);
    strcpy (gmeta->product_id, ESPA_STRING_META_FILL);
    gmeta->ul_corner[0] = gmeta->ul_corner[1] = ESPA_FLOAT_META_FILL;
    gmeta->lr_corner[0] = gmeta->lr_corner[1] = ESPA_FLOAT_META_FILL;
    gmeta->bounding_coords[0] = ESPA_FLOAT_META_FILL;
    gmeta->bounding_coords[1] = ESPA_FLOAT_META_FILL;
    gmeta->bounding_coords[2] = ESPA_FLOAT_META_FILL;
    gmeta->bounding_coords[3] = ESPA_FLOAT_META_FILL;
    gmeta->orientation_angle = ESPA_FLOAT_META_FILL;

    /* Initialize the projection information */
    gmeta->proj_info.proj_type = ESPA_INT_META_FILL;
    gmeta->proj_info.datum_type = ESPA_NODATUM;
    strcpy (gmeta->proj_info.units, ESPA_STRING_META_FILL);
    gmeta->proj_info.ul_corner[0] = gmeta->proj_info.ul_corner[1] =
        ESPA_FLOAT_META_FILL;
    gmeta->proj_info.lr_corner[0] = gmeta->proj_info.lr_corner[1] =
        ESPA_FLOAT_META_FILL;
    strcpy (gmeta->proj_info.grid_origin, ESPA_STRING_META_FILL);

    gmeta->proj_info.utm_zone = ESPA_INT_META_FILL;
    gmeta->proj_info.longitude_pole = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.latitude_true_scale = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.false_easting = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.false_northing = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.standard_parallel1 = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.standard_parallel2 = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.central_meridian = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.origin_latitude = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.sphere_radius = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.semi_major_axis = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.semi_minor_axis = ESPA_FLOAT_META_FILL;
    gmeta->proj_info.satellite_height = ESPA_FLOAT_META_FILL;
}


/******************************************************************************
MODULE:  allocate_band_metadata

PURPOSE:  Allocates memory in the ESPA internal metadata structure for nbands.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error allocating memory for the nbands
SUCCESS         Successfully allocated memory

NOTES:
  1. Initializes the bitmap_description and class_values for each band to NULL
     and sets the nbits, nclass, ncover to 0.
******************************************************************************/
int allocate_band_metadata
(
    Espa_internal_meta_t *internal_meta,  /* I: pointer to internal metadata
                                                structure */
    int nbands                            /* I: number of bands to allocate
                                                for the band field in the
                                                internal_meta */
)
{
    char FUNC_NAME[] = "allocate_band_metadata";   /* function name */
    char errmsg[STR_SIZE];          /* error message */
    Espa_band_meta_t *bmeta = NULL; /* pointer to array of bands metadata */
    int i;                          /* looping variable */

    /* Allocate the number of bands to nbands and the associated pointers */
    internal_meta->nbands = nbands;
    internal_meta->band = calloc (nbands, sizeof (Espa_band_meta_t));
    if (internal_meta->band == NULL)
    {
        sprintf (errmsg, "Allocating ESPA band metadata for %d bands", nbands);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    bmeta = internal_meta->band;

    /* Set the nbits, nclass, ncover fields in the band metadata to 0 for each
       band and initialize the pointers to NULL.  Initialize the other fields to
       fill to make it easy to distinguish if they were populated by reading
       an input metadata file or assigned directly. */
    for (i = 0; i < nbands; i++)
    {
        bmeta[i].nbits = 0;
        bmeta[i].bitmap_description = NULL;
        bmeta[i].nclass = 0;
        bmeta[i].class_values = NULL;
        bmeta[i].ncover = 0;
        bmeta[i].percent_cover = NULL;

        strcpy (bmeta[i].product, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].source, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].name, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].category, ESPA_STRING_META_FILL);
        bmeta[i].data_type = ESPA_UINT8;
        bmeta[i].nlines = ESPA_INT_META_FILL;
        bmeta[i].nsamps = ESPA_INT_META_FILL;
        bmeta[i].fill_value = ESPA_INT_META_FILL;
        bmeta[i].saturate_value = ESPA_INT_META_FILL;
        bmeta[i].scale_factor = ESPA_FLOAT_META_FILL;
        bmeta[i].add_offset = ESPA_FLOAT_META_FILL;
        bmeta[i].resample_method = ESPA_NONE;
        strcpy (bmeta[i].short_name, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].long_name, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].file_name, ESPA_STRING_META_FILL);
        bmeta[i].pixel_size[0] = bmeta[i].pixel_size[1] = ESPA_FLOAT_META_FILL;
        strcpy (bmeta[i].pixel_units, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].data_units, ESPA_STRING_META_FILL);
        bmeta[i].valid_range[0] = bmeta[i].valid_range[1] =
            ESPA_FLOAT_META_FILL;
        bmeta[i].rad_gain = ESPA_FLOAT_META_FILL;
        bmeta[i].rad_bias = ESPA_FLOAT_META_FILL;
        bmeta[i].refl_gain = ESPA_FLOAT_META_FILL;
        bmeta[i].refl_bias = ESPA_FLOAT_META_FILL;
        bmeta[i].k1_const = ESPA_FLOAT_META_FILL;
        bmeta[i].k2_const = ESPA_FLOAT_META_FILL;
        strcpy (bmeta[i].qa_desc, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].app_version, ESPA_STRING_META_FILL);
        strcpy (bmeta[i].production_date, ESPA_STRING_META_FILL);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  allocate_class_metadata

PURPOSE:  Allocates memory in the ESPA band metadata structure for nclasses.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error allocating memory for nclasses
SUCCESS         Successfully allocated memory

NOTES:
******************************************************************************/
int allocate_class_metadata
(
    Espa_band_meta_t *band_meta,  /* I: pointer to band metadata structure */
    int nclass                    /* I: number of classes to allocate for the
                                        band metadata */
)
{
    char FUNC_NAME[] = "allocate_class_metadata";   /* function name */
    char errmsg[STR_SIZE];        /* error message */

    /* Allocate the number of classes to nclass and the associated class_values
       pointer */
    band_meta->nclass = nclass;
    band_meta->class_values = calloc (nclass, sizeof (Espa_class_t));
    if (band_meta->class_values == NULL)
    {
        sprintf (errmsg, "Allocating ESPA band metadata for %d nclasses",
            nclass);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  allocate_percent_coverage_metadata

PURPOSE:  Allocates memory in the ESPA band metadata structure for ncover types.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error allocating memory for ncover types
SUCCESS         Successfully allocated memory

NOTES:
******************************************************************************/
int allocate_percent_coverage_metadata
(
    Espa_band_meta_t *band_meta,  /* I: pointer to band metadata structure */
    int ncover                    /* I: number of cover types to allocate for
                                        the band metadata */
)
{
    char FUNC_NAME[] = "allocate_percent_coverage_metadata"; /* function name */
    char errmsg[STR_SIZE];        /* error message */

    /* Allocate the number of cover types to ncover and the associated cover
       type descripts to the pointer */
    band_meta->ncover = ncover;
    band_meta->percent_cover = calloc (ncover, sizeof (Espa_percent_cover_t));
    if (band_meta->percent_cover == NULL)
    {
        sprintf (errmsg, "Allocating ESPA band metadata for %d cover types",
            ncover);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  allocate_bitmap_metadata

PURPOSE:  Allocates memory in the ESPA band metadata structure for nbits.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error allocating memory for nbits
SUCCESS         Successfully allocated memory

NOTES:
******************************************************************************/
int allocate_bitmap_metadata
(
    Espa_band_meta_t *band_meta,  /* I: pointer to band metadata structure */
    int nbits                     /* I: number of bits to allocate for the
                                        bitmap metadata */
)
{
    char FUNC_NAME[] = "allocate_bitmap_metadata";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    int i;                        /* looping variable */

    /* Allocate the number of bits to nbits and the associated bitmap pointer */
    band_meta->nbits = nbits;
    band_meta->bitmap_description = calloc (nbits, sizeof (char *));
    if (band_meta->bitmap_description == NULL)
    {
        sprintf (errmsg, "Allocating ESPA bitmap description");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    for (i = 0; i < nbits; i++)
    {
        band_meta->bitmap_description[i] = calloc (STR_SIZE, sizeof (char));
        if (band_meta->bitmap_description[i] == NULL)
        {
            sprintf (errmsg, "Allocating ESPA band metadata for %d nbits",
                nbits);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  free_metadata

PURPOSE:  Frees memory in the ESPA internal metadata structure.

RETURN VALUE: N/A

NOTES:
******************************************************************************/
void free_metadata
(
    Espa_internal_meta_t *internal_meta   /* I: pointer to internal metadata
                                                structure */
)
{
    int i, b;                      /* looping variables */

    /* Free the pointers in the band metadata */
    for (i = 0; i < internal_meta->nbands; i++)
    {
        if (internal_meta->band[i].nbits > 0)
        {
            for (b = 0; b < internal_meta->band[i].nbits; b++)
                free (internal_meta->band[i].bitmap_description[b]);
            free (internal_meta->band[i].bitmap_description);
        }

        free (internal_meta->band[i].class_values);
        free (internal_meta->band[i].percent_cover);
    }

    /* Free the band pointer itself */
    if (internal_meta->band)
        free (internal_meta->band);
}


/******************************************************************************
MODULE:  copy_metadata_struct

PURPOSE:  Copies the ESPA internal metadata structure.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error allocating memory for the new metadata structure
SUCCESS         Successfully copied the metadata structure

NOTES:
******************************************************************************/
int copy_metadata_struct
(
    Espa_internal_meta_t *in_meta,  /* I: input metadata struct to be copied */
    Espa_internal_meta_t *out_meta  /* I: metadata structure to be created */
)
{
    char FUNC_NAME[] = "copy_metadata_struct";   /* function name */
    char errmsg[STR_SIZE];          /* error message */
    int i, k;                       /* looping variables */
    Espa_global_meta_t *in_gmeta = &in_meta->global; /* input global metadata
                                                        structure */
    Espa_global_meta_t *out_gmeta = &out_meta->global; /* output global metadata
                                                          structure */
    Espa_band_meta_t *in_bmeta = NULL;  /* ptr to array of input bands meta */
    Espa_band_meta_t *out_bmeta = NULL; /* ptr to array of output bands meta */

    /* Copy the global metadata values */
    strcpy (out_gmeta->data_provider, in_gmeta->data_provider);
    strcpy (out_gmeta->satellite, in_gmeta->satellite);
    strcpy (out_gmeta->instrument, in_gmeta->instrument);
    strcpy (out_gmeta->acquisition_date, in_gmeta->acquisition_date);

    out_gmeta->ul_corner[0] = in_gmeta->ul_corner[0];
    out_gmeta->ul_corner[1] = in_gmeta->ul_corner[1];
    out_gmeta->lr_corner[0] = in_gmeta->lr_corner[0];
    out_gmeta->lr_corner[1] = in_gmeta->lr_corner[1];
    out_gmeta->bounding_coords[0] = in_gmeta->bounding_coords[0];
    out_gmeta->bounding_coords[1] = in_gmeta->bounding_coords[1];
    out_gmeta->bounding_coords[2] = in_gmeta->bounding_coords[2];
    out_gmeta->bounding_coords[3] = in_gmeta->bounding_coords[3];

    out_gmeta->wrs_system = in_gmeta->wrs_system;
    out_gmeta->wrs_path = in_gmeta->wrs_path;
    out_gmeta->wrs_row = in_gmeta->wrs_row;
    strcpy (out_gmeta->scene_center_time, in_gmeta->scene_center_time);
    strcpy (out_gmeta->product_id, in_gmeta->product_id);
    strcpy (out_gmeta->lpgs_metadata_file, in_gmeta->lpgs_metadata_file);
    out_gmeta->orientation_angle = in_gmeta->orientation_angle;
    out_gmeta->solar_zenith = in_gmeta->solar_zenith;
    out_gmeta->solar_azimuth = in_gmeta->solar_azimuth;
    strcpy (out_gmeta->solar_units, in_gmeta->solar_units);
    out_gmeta->earth_sun_dist = in_gmeta->earth_sun_dist;
    strcpy (out_gmeta->level1_production_date,
            in_gmeta->level1_production_date);
    out_gmeta->htile = in_gmeta->htile;
    out_gmeta->vtile = in_gmeta->vtile;

    /* Initialize the projection information */
    out_gmeta->proj_info.proj_type = in_gmeta->proj_info.proj_type;
    out_gmeta->proj_info.datum_type = in_gmeta->proj_info.datum_type;
    strcpy (out_gmeta->proj_info.units, in_gmeta->proj_info.units);
    out_gmeta->proj_info.ul_corner[0] = in_gmeta->proj_info.ul_corner[0];
    out_gmeta->proj_info.ul_corner[1] = in_gmeta->proj_info.ul_corner[1];
    out_gmeta->proj_info.lr_corner[0] = in_gmeta->proj_info.lr_corner[0];
    out_gmeta->proj_info.lr_corner[1] = in_gmeta->proj_info.lr_corner[1];
    strcpy (out_gmeta->proj_info.grid_origin, in_gmeta->proj_info.grid_origin);

    out_gmeta->proj_info.utm_zone = in_gmeta->proj_info.utm_zone;
    out_gmeta->proj_info.longitude_pole = in_gmeta->proj_info.longitude_pole;
    out_gmeta->proj_info.latitude_true_scale =
        in_gmeta->proj_info.latitude_true_scale;
    out_gmeta->proj_info.false_easting = in_gmeta->proj_info.false_easting;
    out_gmeta->proj_info.false_northing = in_gmeta->proj_info.false_northing;
    out_gmeta->proj_info.standard_parallel1 =
        in_gmeta->proj_info.standard_parallel1;
    out_gmeta->proj_info.standard_parallel2 = 
        in_gmeta->proj_info.standard_parallel2;
    out_gmeta->proj_info.central_meridian = 
        in_gmeta->proj_info.central_meridian;
    out_gmeta->proj_info.origin_latitude = in_gmeta->proj_info.origin_latitude;
    out_gmeta->proj_info.sphere_radius = in_gmeta->proj_info.sphere_radius;
    out_gmeta->proj_info.semi_major_axis = in_gmeta->proj_info.semi_major_axis;
    out_gmeta->proj_info.semi_minor_axis = in_gmeta->proj_info.semi_minor_axis;
    out_gmeta->proj_info.satellite_height = 
        in_gmeta->proj_info.satellite_height;

    /* Allocate the number of bands in the output metadata */
    out_meta->nbands = in_meta->nbands;
    out_meta->band = calloc (in_meta->nbands, sizeof (Espa_band_meta_t));
    if (out_meta->band == NULL)
    {
        sprintf (errmsg, "Allocating ESPA band metadata for %d bands",
            in_meta->nbands);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    in_bmeta = in_meta->band;
    out_bmeta = out_meta->band;

    /* Set the nbits, nclass, ncover fields in the band metadata to 0 for each
       band and initialize the pointers to NULL.  Initialize the other fields to
       fill to make it easy to distinguish if they were populated by reading
       an input metadata file or assigned directly. */
    for (i = 0; i < in_meta->nbands; i++)
    {
        strcpy (out_bmeta[i].product, in_bmeta[i].product);
        strcpy (out_bmeta[i].source, in_bmeta[i].source);
        strcpy (out_bmeta[i].name, in_bmeta[i].name);
        strcpy (out_bmeta[i].category, in_bmeta[i].category);
        out_bmeta[i].data_type = in_bmeta[i].data_type;
        out_bmeta[i].nlines = in_bmeta[i].nlines;
        out_bmeta[i].nsamps = in_bmeta[i].nsamps;
        out_bmeta[i].fill_value = in_bmeta[i].fill_value;
        out_bmeta[i].saturate_value = in_bmeta[i].saturate_value;
        out_bmeta[i].scale_factor = in_bmeta[i].scale_factor;
        out_bmeta[i].add_offset = in_bmeta[i].add_offset;
        out_bmeta[i].resample_method = in_bmeta[i].resample_method;

        strcpy (out_bmeta[i].short_name, in_bmeta[i].short_name);
        strcpy (out_bmeta[i].long_name, in_bmeta[i].long_name);
        strcpy (out_bmeta[i].file_name, in_bmeta[i].file_name);
        out_bmeta[i].pixel_size[0] = in_bmeta[i].pixel_size[0];
        out_bmeta[i].pixel_size[1] = in_bmeta[i].pixel_size[1];

        strcpy (out_bmeta[i].pixel_units, in_bmeta[i].pixel_units);
        strcpy (out_bmeta[i].data_units, in_bmeta[i].data_units);
        out_bmeta[i].valid_range[0] = in_bmeta[i].valid_range[0];
        out_bmeta[i].valid_range[1] = in_bmeta[i].valid_range[1];

        out_bmeta[i].rad_gain = in_bmeta[i].rad_gain;
        out_bmeta[i].rad_bias = in_bmeta[i].rad_bias;
        out_bmeta[i].refl_gain = in_bmeta[i].refl_gain;
        out_bmeta[i].refl_bias = in_bmeta[i].refl_bias;
        out_bmeta[i].k1_const = in_bmeta[i].k1_const;
        out_bmeta[i].k2_const = in_bmeta[i].k2_const;
        strcpy (out_bmeta[i].qa_desc, in_bmeta[i].qa_desc);
        strcpy (out_bmeta[i].app_version, in_bmeta[i].app_version);
        strcpy (out_bmeta[i].production_date, in_bmeta[i].production_date);

        /* bit description */
        out_bmeta[i].nbits = in_bmeta[i].nbits;
        if (in_bmeta[i].nbits != 0)
        {
            if (allocate_bitmap_metadata (&out_bmeta[i], out_bmeta[i].nbits) !=
                SUCCESS)
            {  /* Error messages already printed */
                return (ERROR);
            }

            for (k = 0; k < in_bmeta[i].nbits; k++)
            {
                strcpy (out_bmeta[i].bitmap_description[k],
                    in_bmeta[i].bitmap_description[k]);
            }
        }

        /* class description */
        out_bmeta[i].nclass = in_bmeta[i].nclass;
        if (in_bmeta[i].nclass != 0)
        {
            if (allocate_class_metadata (&out_bmeta[i], out_bmeta[i].nclass) !=
                SUCCESS)
            {  /* Error messages already printed */
                return (ERROR);
            }

            for (k = 0; k < in_bmeta[i].nclass; k++)
            {
                 out_bmeta[i].class_values[k].class =
                     in_bmeta[i].class_values[k].class;
                 strcpy (out_bmeta[i].class_values[k].description,
                     in_bmeta[i].class_values[k].description);
            }
        }

        /* percent coverage */
        out_bmeta[i].ncover = in_bmeta[i].ncover;
        if (in_bmeta[i].ncover != 0)
        {
            if (allocate_percent_coverage_metadata (&out_bmeta[i],
                out_bmeta[i].ncover) != SUCCESS)
            {  /* Error messages already printed */
                return (ERROR);
            }

            for (k = 0; k < in_bmeta[i].ncover; k++)
            {
                 out_bmeta[i].percent_cover[k].percent =
                     in_bmeta[i].percent_cover[k].percent;
                 strcpy (out_bmeta[i].percent_cover[k].description,
                     in_bmeta[i].percent_cover[k].description);
            }
        }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  print_element_names

PURPOSE:  Print the information for the elements in the document tree,
starting at the node provided.

RETURN VALUE:  N/A

NOTES:
  1. Prints to stdout.
******************************************************************************/
void print_element_names
(
    xmlNode *a_node   /* I: pointer to the current node in the tree to start
                            printing */
)
{
    xmlNode *cur_node = NULL;   /* pointer to the current node */

    /* Start at the input node and traverse the tree, visiting all the children
       and siblings */
    for (cur_node = a_node; cur_node;
         cur_node = xmlNextElementSibling (cur_node))
    {
        /* Only print the ELEMENT node types */
        if (cur_node->type == XML_ELEMENT_NODE) 
        {
            /* Print out the name of the element */
            xmlAttrPtr attr;     /* pointer to the element attributes */
            printf ("node type: Element, name: %s", cur_node->name);

            /* Print out the namespace info as well */
            xmlNsPtr ns = cur_node->nsDef;
            while (ns != 0)
            {
                printf (" with namespace: %s %p\n", ns->href, ns->prefix);
                ns = ns->next;
            }
            printf("\n");

            /* Print out the attribute properties for this element */
            for (attr = cur_node->properties; attr != NULL; attr = attr->next)
            {
                xmlChar *v = xmlGetProp (cur_node, attr->name);
                if (attr->ns != NULL)
                {
                    if (attr->ns->prefix != NULL)
                    {
                        printf (" with namespace: %s %p\n", attr->ns->href,
                            attr->ns->prefix);
                    }
                    else
                    {
                        printf (" with namespace: %s\n", attr->ns->href);
                    }
                }
                printf (" @%s=%s ", attr->name, v);
                xmlFree (v);
            }
            printf ("\n");
        }
        else if (cur_node->type == XML_TEXT_NODE) 
        {
            /* Print out the text for the element */
            printf ("   node type: Text, content: %s\n", cur_node->content);
        }

        print_element_names (cur_node->children);
    }
}

