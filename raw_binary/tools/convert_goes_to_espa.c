/*****************************************************************************
FILE: convert_goes_to_espa
  
PURPOSE: Contains functions for converting the GOES-R ABI products to the ESPA
internal raw binary file format.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
  1. The XML metadata format parsed or written via this library follows the
     ESPA internal metadata format found in ESPA Raw Binary Format v1.0.doc.
     The schema for the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/schema/espa_internal_metadata_v1_0.xsd.
*****************************************************************************/
#include <getopt.h>
#include "convert_goes_to_espa.h"

/******************************************************************************
MODULE: usage

PURPOSE: Prints the usage information for this application.

RETURN VALUE:
Type = None

NOTES:
******************************************************************************/
void usage ()
{
    printf ("convert_goes_to_espa converts the GOES-R ABI products to the ESPA "
            "internal format (XML metadata file and associated raw binary "
            "files).\n\n");
    printf ("usage: convert_goes_to_espa "
            "--red=input_red_filename "
            "--nir=input_nir_filename "
            "[--del_src_files]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -red: name of input GOES-R ABI netCDF band 2 (red) file\n");
    printf ("    -nir: name of input GOES-R ABI netCDF band 3 (NIR) file\n");
    printf ("    -del_src_files: if specified the source netCDF file will "
            "be removed.\n");
    printf ("\nExample: convert_goes_to_espa "
            "--red=OR_ABI-L2-CMIPC-M3C02_G16_s20171721702192_e20171721704565_c20171721705067.nc --nir=OR_ABI-L2-CMIPC-M3C03_G16_s20171721702192_e20171721704565_c20171721705037.nc\n");
}


/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

NOTES:
  1. Memory is allocated for the input and output files.  All of these should
     be character pointers set to NULL on input.  The caller is responsible
     for freeing the allocated memory upon successful return.
******************************************************************************/
short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **red_infile,    /* O: address of input GOES-R red filename */
    char **nir_infile,    /* O: address of input GOES-R NIR filename */
    char **xml_outfile,   /* O: address of output XML filename */
    bool *del_src         /* O: should source files be removed? */
)
{
    int c;                           /* current argument index */
    int option_index;                /* index for the command-line option */
    char *cptr = NULL;               /* pointer to .nc in netCDF filename */
    char errmsg[STR_SIZE];           /* error message */
    char tmpfile[STR_SIZE];          /* temporary filename */
    char FUNC_NAME[] = "get_args";   /* function name */
    static int del_flag = 0;         /* flag for removing the source files */
    static struct option long_options[] =
    {
        {"del_src_files", no_argument, &del_flag, 1},
        {"red", required_argument, 0, 'r'},
        {"nir", required_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Loop through all the cmd-line options */
    opterr = 0;   /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {   /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
     
            case 'h':  /* help */
                usage ();
                return (ERROR);
                break;

            case 'r':  /* GOES-R red infile */
                *red_infile = strdup (optarg);
                break;
     
            case 'n':  /* GOES-R NIR infile */
                *nir_infile = strdup (optarg);
                break;
     
            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind-1]);
                error_handler (true, FUNC_NAME, errmsg);
                usage ();
                return (ERROR);
                break;
        }
    }

    /* Make sure the input GOES-R files were specified */
    if (*red_infile == NULL)
    {
        sprintf (errmsg, "GOES-R red input file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    if (*nir_infile == NULL)
    {
        sprintf (errmsg, "GOES-R NIR input file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Generate the XML filename from the netCDF band 2 and band 3 filenames.
       Find the _c{blah}.nc and change that to .xml. Bands 2 and 3 should be
       the same except for (maybe) the _c{date} portion of the name which
       refers to the creation date. Strip off the creation date from band 2.
       Example: OR_ABI-L2-CMIPC-M3C02_G16_s20171721702192_e20171721704565_c20171721705067.nc */
    *xml_outfile = strdup (*red_infile);
    strcpy (tmpfile, *xml_outfile);
    cptr = strrchr (tmpfile, '_');
    *cptr = '\0';
    if (tmpfile == NULL)
    {
        sprintf (errmsg, "XML output file was not correctly generated");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Now remove the channel information from the XML filename, since it is
       different for bands 2 (C02) and 3 (C03). The XML filename needs to be
       generic. */
    strncpy (*xml_outfile, tmpfile, 18);
    strcpy (&(*xml_outfile)[18], &tmpfile[21]);
    sprintf (*xml_outfile, "%s.xml", *xml_outfile);
    if (*xml_outfile == NULL)
    {
        sprintf (errmsg, "XML output file was not correctly generated");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
printf ("XML filename: %s\n", *xml_outfile);

    /* Check the delete source files flag */
    if (del_flag)
        *del_src = true;

    return (SUCCESS);
}


/******************************************************************************
MODULE:  main

PURPOSE:  Converts the GOES-R ABI netCDF product to the ESPA internal format
(XML metadata file and associated raw binary files).

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error doing the conversion
SUCCESS         No errors encountered

NOTES:
******************************************************************************/
int main (int argc, char** argv)
{
    char *red_infile = NULL;    /* I: input GOES-R red filename */
    char *nir_infile = NULL;    /* I: input GOES-R NIR filename */
    char *xml_outfile = NULL;   /* output XML filename */
    bool append_bands;          /* should XML file be created or bands simply
                                   be appended to an existing XML file;  true
                                   will append the bands to an existing file */
    bool del_src = false;       /* should source files be removed? */

    /* Read the command-line arguments */
    if (get_args (argc, argv, &red_infile, &nir_infile, &xml_outfile,
        &del_src) != SUCCESS)
    {   /* get_args already printed the error message */
        exit (EXIT_FAILURE);
    }

    /* Convert the GOES netCDF red band to ESPA raw binary and XML.  The XML
       file will be created with this call to convert_goes_to_espa. */
    printf ("Converting red band from GOES to ESPA ... %s\n", red_infile);
    append_bands = false;
    if (convert_goes_to_espa (red_infile, xml_outfile, append_bands, del_src)
        != SUCCESS)
    {  /* Error messages already written */
        exit (EXIT_FAILURE);
    }

    /* Convert the GOES netCDF NIR band to ESPA raw binary and XML.  The XML
       file has already been created, so this time the bands will simply be
       appended. */
    printf ("\nConverting NIR band from GOES to ESPA ... %s\n", nir_infile);
    append_bands = true;
    if (convert_goes_to_espa (nir_infile, xml_outfile, append_bands, del_src)
        != SUCCESS)
    {  /* Error messages already written */
        exit (EXIT_FAILURE);
    }

    /* Free the pointers */
    free (red_infile);
    free (nir_infile);
    free (xml_outfile);

    /* Successful completion */
    exit (EXIT_SUCCESS);
}
