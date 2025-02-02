/*****************************************************************************
FILE: parse_sentinel_metadata.c
  
PURPOSE: Contains functions for parsing the Sentinel metadata files.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
  1. The XML metadata format parsed or written via this library follows the
     ESPA internal metadata format found in ESPA Raw Binary Format v1.0.doc.
     The schema for the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/schema/espa_internal_metadata_v1_0.xsd.
  2. This code relies on the libxml2 library developed for the Gnome project.
  3. The information on the Sentinel-2 L1C metadata files (MTD_MSIL1C.xml and
     MTD_TL.xml) can be found in the S2_MSI_Product_Specification.pdf file.
*****************************************************************************/
#include "espa_metadata.h"
#include "parse_sentinel_metadata.h"


/******************************************************************************
MODULE:  add_mean_solar_angles

PURPOSE: Add the current tile-based mean solar angles to the metadata
structure and process the children of this node.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the tile solar angles metadata
SUCCESS         Successful parse of the tile solar angles metadata

NOTES:
******************************************************************************/
int add_mean_solar_angles
(
    xmlNode *a_node,     /* I: pointer to the element node to process */
    Espa_global_meta_t *gmeta   /* I: global metadata structure */
)
{
    char FUNC_NAME[] = "add_mean_solar_angles";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    xmlNode *cur_node = NULL;     /* pointer to the current node */
    xmlNode *child_node = NULL;   /* pointer to the child node */

    /* Process the siblings in the Mean_Sun_Angle element */
    for (cur_node = a_node->children; cur_node;
         cur_node = xmlNextElementSibling (cur_node))
    {
        /* Process the zenith angle */
        if (xmlStrEqual (cur_node->name,
            (const xmlChar *) "ZENITH_ANGLE"))
        {
            /* Expect the child node to be a text node containing the value
               of this field */
            child_node = cur_node->children;
            if (child_node == NULL || child_node->type != XML_TEXT_NODE)
            {
                sprintf (errmsg, "Processing tile solar angle element: %s.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Copy the content of the child node into the value for this
               field */
            gmeta->solar_zenith = atof ((const char *) child_node->content);
        }

        /* Process the azimuth angle */
        if (xmlStrEqual (cur_node->name,
            (const xmlChar *) "AZIMUTH_ANGLE"))
        {
            /* Expect the child node to be a text node containing the value
               of this field */
            child_node = cur_node->children;
            if (child_node == NULL || child_node->type != XML_TEXT_NODE)
            {
                sprintf (errmsg, "Processing tile solar angle element: %s.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Copy the content of the child node into the value for this
               field */
            gmeta->solar_azimuth = atof ((const char *) child_node->content);
        }

    }  /* end for cur_node */

    /* Set the units to degrees for these angles */
    strcpy (gmeta->solar_units, "degrees");

    return (SUCCESS);
}


/******************************************************************************
MODULE:  add_mean_viewing_angles

PURPOSE: Add the current tile-based mean viewing angles to the metadata
structure and process the children of this node.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the tile viewing angles metadata
SUCCESS         Successful parse of the tile viewing angles metadata

NOTES:
******************************************************************************/
int add_mean_viewing_angles
(
    xmlNode *a_node,     /* I: pointer to the element node to process */
    Espa_global_meta_t *gmeta   /* I: global metadata structure */
)
{
    char FUNC_NAME[] = "add_mean_viewing_angles";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    xmlNode *cur_node = NULL;     /* pointer to the current node */
    xmlNode *child_node = NULL;   /* pointer to the child node */

    /* Process the siblings in the Mean_Viewing_Incidence_Angle element */
    for (cur_node = a_node->children; cur_node;
         cur_node = xmlNextElementSibling (cur_node))
    {
        /* Process the zenith angle */
        if (xmlStrEqual (cur_node->name,
            (const xmlChar *) "ZENITH_ANGLE"))
        {
            /* Expect the child node to be a text node containing the value
               of this field */
            child_node = cur_node->children;
            if (child_node == NULL || child_node->type != XML_TEXT_NODE)
            {
                sprintf (errmsg, "Processing tile viewing angle element: %s.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Copy the content of the child node into the value for this
               field */
            gmeta->view_zenith = atof ((const char *) child_node->content);
        }

        /* Process the azimuth angle */
        if (xmlStrEqual (cur_node->name,
            (const xmlChar *) "AZIMUTH_ANGLE"))
        {
            /* Expect the child node to be a text node containing the value
               of this field */
            child_node = cur_node->children;
            if (child_node == NULL || child_node->type != XML_TEXT_NODE)
            {
                sprintf (errmsg, "Processing tile view angle element: %s.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Copy the content of the child node into the value for this
               field */
            gmeta->view_azimuth = atof ((const char *) child_node->content);
        }

    }  /* end for cur_node */

    /* Set the units to degrees for these angles */
    strcpy (gmeta->view_units, "degrees");

    return (SUCCESS);
}


/******************************************************************************
MODULE:  add_tile_geocoding_metadata

PURPOSE: Add the current tile-based geocoding information to the metadata
structure and process the children of this node.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the tile geocoding metadata
SUCCESS         Successful parse of the tile geocoding metadata

NOTES:
******************************************************************************/
int add_tile_geocoding_metadata
(
    xmlNode *a_node,     /* I: pointer to the element node to process */
    Espa_global_meta_t *gmeta,  /* I: global metadata structure */
    int *nrows,          /* O: number of rows for 10m, 20m, 60m res */
    int *ncols           /* O: number of columns for 10m, 20m, 60m res */
)
{
    char FUNC_NAME[] = "add_tile_geocoding_metadata";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    char tmpstr[STR_SIZE];        /* temporary string */
    char tmp_zone[STR_SIZE];      /* temporary zone string */
    xmlNode *cur_node = NULL;     /* pointer to the current node */
    xmlNode *child_node = NULL;   /* pointer to the child node */
    xmlNode *sib_node = NULL;     /* pointer to the sibling node */
    xmlNode *sib_child_node = NULL; /* pointer to the sibling's child node */
    xmlAttrPtr attr = NULL;       /* pointer to the element attributes */
    xmlChar *attr_val = NULL;     /* attribute value */
    int count;                    /* number of chars copied in snprintf */
    int index;                    /* index for nrows/ncols arrays */
    int ulx[NUM_SENTINEL_RES];    /* ULx for sentinel resolutions */
    int uly[NUM_SENTINEL_RES];    /* ULy for sentinel resolutions */

    /* Process the siblings in the Geocoding element */
    for (cur_node = a_node->children; cur_node;
         cur_node = xmlNextElementSibling (cur_node))
    {
        /* Process the horizontal CS name for the datum, projection, and
           zone code */
        if (xmlStrEqual (cur_node->name,
            (const xmlChar *) "HORIZONTAL_CS_NAME"))
        {
            /* Expect the child node to be a text node containing the value
               of this field */
            child_node = cur_node->children;
            if (child_node == NULL || child_node->type != XML_TEXT_NODE)
            {
                sprintf (errmsg, "Processing tile metadata element: %s.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Copy the content of the child node into the value for this
               field */
            count = snprintf (tmpstr, sizeof (tmpstr), "%s",
                (const char *) child_node->content);
            if (count < 0 || count >= sizeof (tmpstr))
            {
                sprintf (errmsg, "Overflow of tmpstr string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Make sure this is UTM WGS84, as that's the promoted projection
               that we are going to use */
            if (strstr (tmpstr, "WGS84 / UTM") == NULL)
            {
                sprintf (errmsg, "Datum and projection should be "
                    "WGS84 / UTM but instead it is %s\n", tmpstr);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            gmeta->proj_info.proj_type = GCTP_UTM_PROJ;
            gmeta->proj_info.datum_type = ESPA_WGS84;
            strcpy (gmeta->proj_info.units, "meters");

            /* Grab the zone number. Negate the integer value if the zone is a
               southern zone. */
            strcpy (tmp_zone, &tmpstr[16]);
            gmeta->proj_info.utm_zone = atoi(tmp_zone);
            if (strstr (tmp_zone, "S"))
                gmeta->proj_info.utm_zone = -gmeta->proj_info.utm_zone;
        }

        /* Process the nrows/ncols */
        if (xmlStrEqual (cur_node->name, (const xmlChar *) "Size"))
        {
            /* Get the resolution attribute */
            index = -99;
            for (attr = cur_node->properties; attr != NULL; attr = attr->next)
            {
                attr_val = xmlGetProp (cur_node, attr->name);
                if (xmlStrEqual (attr->name, (const xmlChar *) "resolution"))
                {
                    if (xmlStrEqual (attr_val, (const xmlChar *) "10"))
                        index = 0;
                    else if (xmlStrEqual (attr_val, (const xmlChar *) "20"))
                        index = 1;
                    else if (xmlStrEqual (attr_val, (const xmlChar *) "60"))
                        index = 2;
                    else
                    {
                        sprintf (errmsg, "Unknown resolution for the Sentinel "
                            "tile geocoding specified (%s).", attr_val);
                        error_handler (false, FUNC_NAME, errmsg);
                    }
                }
                else
                {
                    sprintf (errmsg, "unknown attribute for element (%s): %s",
                        cur_node->name, attr->name);
                    error_handler (false, FUNC_NAME, errmsg);
                }
                xmlFree (attr_val);
            }

            /* Make sure we found the desired resolution */
            if (index == -99)
            {
                sprintf (errmsg, "resolution attribute not found for the "
                    "current Size element: %s", cur_node->name);
                error_handler (false, FUNC_NAME, errmsg);
            }

            /* Get the NROWS/NCOLS children for this resolution */
            for (sib_node = cur_node->children; sib_node;
                 sib_node = xmlNextElementSibling (sib_node))
            {
                /* Get the child node of the sibling */
                sib_child_node = sib_node->children;

                /* If this is the NROWS then store number of rows */
                if (xmlStrEqual (sib_node->name, (const xmlChar *) "NROWS"))
                {
                    nrows[index] =
                        atoi ((const char *) sib_child_node->content);
                }

                /* If this is the NCOLS then store number of columns */
                if (xmlStrEqual (sib_node->name, (const xmlChar *) "NCOLS"))
                {
                    ncols[index] =
                        atoi ((const char *) sib_child_node->content);
                }
            }
        }

        /* Process the UL x/y */
        if (xmlStrEqual (cur_node->name, (const xmlChar *) "Geoposition"))
        {
            /* Get the resolution attribute */
            index = -99;
            for (attr = cur_node->properties; attr != NULL; attr = attr->next)
            {
                attr_val = xmlGetProp (cur_node, attr->name);
                if (xmlStrEqual (attr->name, (const xmlChar *) "resolution"))
                {
                    if (xmlStrEqual (attr_val, (const xmlChar *) "10"))
                        index = 0;
                    else if (xmlStrEqual (attr_val, (const xmlChar *) "20"))
                        index = 1;
                    else if (xmlStrEqual (attr_val, (const xmlChar *) "60"))
                        index = 2;
                    else
                    {
                        sprintf (errmsg, "Unknown resolution for the Sentinel "
                            "tile geocoding specified (%s).", attr_val);
                        error_handler (false, FUNC_NAME, errmsg);
                    }
                }
                else
                {
                    sprintf (errmsg, "unknown attribute for element (%s): %s",
                        cur_node->name, attr->name);
                    error_handler (false, FUNC_NAME, errmsg);
                }
                xmlFree (attr_val);
            }

            /* Make sure we found the desired resolution */
            if (index == -99)
            {
                sprintf (errmsg, "resolution attribute not found for the "
                    "current Geoposition element: %s", cur_node->name);
                error_handler (false, FUNC_NAME, errmsg);
            }

            /* Get the ULX/ULY children for this resolution */
            for (sib_node = cur_node->children; sib_node;
                 sib_node = xmlNextElementSibling (sib_node))
            {
                /* Get the child node of the sibling */
                sib_child_node = sib_node->children;

                /* If this is the ULX then store the value */
                if (xmlStrEqual (sib_node->name, (const xmlChar *) "ULX"))
                    ulx[index] = atoi ((const char *) sib_child_node->content);

                /* If this is the ULY then store the value */
                if (xmlStrEqual (sib_node->name, (const xmlChar *) "ULY"))
                    uly[index] = atoi ((const char *) sib_child_node->content);
            }
        }
    }  /* end for cur_node */

    /* The global metadata only supports one set of corner coordinates. It
       is expected the center of the corners will be the same, so we will use
       the projection corners for the 10m resolution. */
    gmeta->proj_info.ul_corner[0] = ulx[0];
    gmeta->proj_info.ul_corner[1] = uly[0];
    gmeta->proj_info.lr_corner[0] = ulx[0] + ncols[0] * 10.0;
    gmeta->proj_info.lr_corner[1] = uly[0] - nrows[0] * 10.0;
    strcpy (gmeta->proj_info.grid_origin, "UL");

    return (SUCCESS);
}


/******************************************************************************
MODULE:  parse_sentinel_tile_xml_into_struct

PURPOSE: Parse the Sentinel L1C tile level XML document(MTD_TL.xml) data into
the ESPA internal metadata structure.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the metadata elements
SUCCESS         Successful parse of the metadata values

NOTES:
1. Uses a stack of character strings to keep track of the nodes that have
   been parsed.  The stack must be allocated before calling this routine.
2. The very first mean viewing angle is used in this case.  There is actually
   one mean view angle for each band, but the first one in the metadata is the
   one which is used.
******************************************************************************/
int parse_sentinel_tile_xml_into_struct
(
    xmlNode *a_node,     /* I: pointer to the current node */
    Espa_internal_meta_t *metadata,   /* I: ESPA internal metadata structure
                                            to be filled */
    int *top_of_stack,   /* I: pointer to top of the stack */
    char **stack,        /* I: stack to use for parsing */
    int *nrows,          /* O: number of rows for 10m, 20m, 60m res */
    int *ncols           /* O: number of columns for 10m, 20m, 60m res */
)
{
    char FUNC_NAME[] = "parse_sentinel_tile_xml_into_struct"; /* func name */
    char errmsg[STR_SIZE];       /* error message */
    char tmp_date[STR_SIZE];     /* temporary date string */
    int count;                   /* number of chars copied in snprintf */
    char *curr_stack_element = NULL;  /* element popped from the stack */
    xmlNode *cur_node = NULL;    /* pointer to the current node */
    xmlNode *child_node = NULL;  /* pointer to the child node */
    bool view_angle_found = false;  /* boolean to specify when the first view
                                    angle was processed */
    bool skip_child;             /* boolean to specify the children of this
                                    node should not be processed */
    Espa_global_meta_t *gmeta = &metadata->global;
                                 /* global metadata structure */

    /* Start at the input node and traverse the tree, visiting all the children
       and siblings */
    for (cur_node = a_node; cur_node;
         cur_node = xmlNextElementSibling (cur_node))
    {
        /* Process the children of this element unless otherwise specified */
        skip_child = false;

        /* Only print the ELEMENT node types */
        if (cur_node->type == XML_ELEMENT_NODE) 
        {
            /* Push the element to the stack and turn the boolean on if this
               is either Granule metadata */
            //printf ("***Pushed %s\n", cur_node->name); fflush (stdout);
            if (push (top_of_stack, stack, (const char *) cur_node->name))
            {
                sprintf (errmsg, "Pushing element '%s' to the stack.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Process the sensing time for the acquisition date/time.  The
               sensing time of the tile is defined as the time stamp of the
               first line of the Granule. */
            if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "SENSING_TIME"))
            {
                /* Expect the child node to be a text node containing the
                   value of this field */
                child_node = cur_node->children;
                if (child_node == NULL || child_node->type != XML_TEXT_NODE)
                {
                    sprintf (errmsg, "Processing tile metadata element: %s.",
                        cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Copy the content of the child node into the value for this
                   field */
                count = snprintf (tmp_date, sizeof (tmp_date), "%s",
                    (const char *) child_node->content);
                if (count < 0 || count >= sizeof (tmp_date))
                {
                    sprintf (errmsg, "Overflow of tmp_date string");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Copy the date into the acquisition date field */
                strncpy (gmeta->acquisition_date, tmp_date, 10);
                gmeta->acquisition_date[10] = '\0';
            }

            /* Process the geocoding metadata */
            else if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "Tile_Geocoding"))
            {
                /* Add the geocoding information to the band and global
                   metadata */
                if (add_tile_geocoding_metadata (cur_node, gmeta, nrows, ncols))
                {
                    sprintf (errmsg, "Consuming Tile Geocoding elements '%s'.",
                        cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Skip processing the children of this node, since they
                   will be handled by the global metadata parser */
                skip_child = true;
            }

            /* Process the mean solar angles */
            else if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "Mean_Sun_Angle"))
            {
                /* Add the solar angles to the global metadata */
                if (add_mean_solar_angles (cur_node, gmeta))
                {
                    sprintf (errmsg, "Consuming mean solar angle elements "
                        "'%s'.", cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Skip processing the children of this node, since they
                   will be handled by the global metadata parser */
                skip_child = true;
            }

            /* Process the mean viewing angles, only grab the first band in
               the list */
            else if (!view_angle_found && (xmlStrEqual (cur_node->name,
                (const xmlChar *) "Mean_Viewing_Incidence_Angle")))
            {
                /* Add the viewing angles to the global metadata */
                if (add_mean_viewing_angles (cur_node, gmeta))
                {
                    sprintf (errmsg, "Consuming mean viewing angle elements "
                        "'%s'.", cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Set to true once we have found a mean viewing incidence
                   angle in the metadata.  We will just read the first one. */
                view_angle_found = true;

                /* Skip processing the children of this node, since they
                   will be handled by the global metadata parser */
                skip_child = true;
            }

            /* Print out the name of the element */
            //xmlAttrPtr attr;     /* pointer to the element attributes */
            //printf ("node type: Element, name: %s\n", cur_node->name);

            ///* Print out the attribute properties for this element */
            //for (attr = cur_node->properties; attr != NULL; attr = attr->next)
            //{
            //    xmlChar *v = xmlGetProp (cur_node, attr->name);
            //    printf (" @%s=%s ", attr->name, v);
            //    xmlFree (v);
            //}
            //printf ("\n"); fflush (stdout);
        }
        else if (cur_node->type == XML_TEXT_NODE) 
        {
            /* Print out the text for the element */
            //printf ("   node type: Text, content: %s\n", cur_node->content);
        }

        /* Parse the children of this node if they haven't been consumed
           elsewhere */
        if (!skip_child)
        {
            if (parse_sentinel_tile_xml_into_struct (cur_node->children,
                metadata, top_of_stack, stack, nrows, ncols))
            {
                sprintf (errmsg, "Parsing the children of this element '%s'.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* Done with the element and its siblings so pop the element name off
           the stack */
        if (cur_node->type == XML_ELEMENT_NODE)
        {
            curr_stack_element = pop (top_of_stack, stack);
            if (curr_stack_element == NULL)
            {
                sprintf (errmsg, "Popping elements off the stack.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            //printf ("***Popped %s\n", curr_stack_element); fflush (stdout);
        }
    }  /* for cur_node */

    return (SUCCESS);
}


/******************************************************************************
MODULE:  parse_sentinel_tile_metadata

PURPOSE: Parse the Sentinel L1C tile metadata file (MTD_TL.xml) and populate
the associated ESPA internal metadata file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the Sentinel metadata elements
SUCCESS         Successful parse of the Sentinel metadata values

NOTES:
1. Uses a stack of character strings to keep track of the nodes that have been
   found in the metadata document.
2. For debugging purposes
   xmlDocDump (stderr, doc);
   can be used to dump/print the XML doc to the screen.
******************************************************************************/
int parse_sentinel_tile_metadata
(
    char *metafile,                 /* I: input Sentinel tile metadata file */
    Espa_internal_meta_t *metadata  /* I: input metadata structure which has
                                          been initialized via
                                          init_metadata_struct */
)
{
    char FUNC_NAME[] = "parse_sentinel_tile_metadata";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    xmlTextReaderPtr reader;  /* reader for the XML file */
    xmlDocPtr doc = NULL;     /* document tree pointer */
    xmlNodePtr current=NULL;  /* pointer to the current node */
    int i;                    /* looping variable */
    int status;               /* return status */
    int nodeType;             /* node type (element, text, attribute, etc.) */
    int top_of_stack;         /* top of the stack */
    int count;                /* number of chars copied in snprintf */
    int nrows[NUM_SENTINEL_RES]; /* num rows for each sentinel resolution */
    int ncols[NUM_SENTINEL_RES]; /* num columns for each sentinel resolution */
    char **stack = NULL;      /* stack to keep track of elements in the tree */
    Espa_band_meta_t *bmeta;  /* band metadata pointer to all bands */

    /* Establish the reader for this metadata file */
    reader = xmlNewTextReaderFilename (metafile);
    if (reader == NULL)
    {
        sprintf (errmsg, "Setting up reader for %s", metafile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Use the reader to parse the XML file, looking at each of the nodes,
       until the entire file has been parsed.  Start by reading the first
       node in the file. */
    status = xmlTextReaderRead (reader);
    while (status == 1)
    {
        /* Determine what kind of node the reader is at (element, end element,
           attribute, text/white space) and handle the information as desired */
        nodeType = xmlTextReaderNodeType (reader);
        if (nodeType == -1)
        {
            sprintf (errmsg, "Getting node type");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        switch (nodeType)
        {
            case XML_READER_TYPE_ELEMENT:
            {  /* Node is an element (ex. <global_metadata> */
                xmlNodePtr node=NULL;
                if (doc==NULL)
                {
                    doc=xmlNewDoc (BAD_CAST "1.0");
                }

                /* Get the URI defining the namespace for this node */
                if (xmlTextReaderConstNamespaceUri (reader) != NULL)
                {
                    /* Search the namespace for a document */
                    xmlNsPtr ns = xmlSearchNs (doc, current,
                        xmlTextReaderConstNamespaceUri(reader));

                    /* Create a tree node for this element in the XML file
                       using the element name */
                    node = xmlNewNode (ns, xmlTextReaderConstName(reader));

                    /* If the namespace is empty (i.e. root) then create a
                       new namespace pointer with this node */
                    if (ns == NULL)
                    {
                        ns = xmlNewNs (node,
                            xmlTextReaderConstNamespaceUri(reader),
                            xmlTextReaderConstPrefix(reader));
                    }
                }
                else
                {
                    /* Create a tree node for this element in the XML file
                       using the element name */
                    node = xmlNewNode (0, xmlTextReaderConstName(reader));
                }

                /* Set the element as the root if appropriate otherwise add
                   it as a child to the previous element */
                if (current == NULL)
                {
                    xmlDocSetRootElement (doc, node);
                }
                else
                {
                    xmlAddChild (current, node);
                }
                current = node;

                /* If the element has attributes, then handle them */
                if (xmlTextReaderHasAttributes (reader))
                {
                    /* Get the number of attributes and then process each one */
                    int i;
                    int n_att = xmlTextReaderAttributeCount (reader);
                    for (i = 0; i < n_att; i++)
                    {
                        /* Read each attribute, obtain the name and value,
                           then add it as a property for this node in the
                           tree */
                        const xmlChar *k = NULL;
                        xmlChar *v = NULL;
                        xmlTextReaderMoveToAttributeNo (reader, i);
                        k = xmlTextReaderConstName (reader);
                        v = xmlTextReaderValue (reader);
                        if (xmlTextReaderConstNamespaceUri (reader) != NULL)
                        {
                            if (!xmlStrEqual (
                                xmlTextReaderConstNamespaceUri(reader),
                                BAD_CAST "http://www.w3.org/2000/xmlns/"))
                            {
                                /* Search the namespace for the document */
                                xmlNsPtr ns = xmlSearchNs (doc, current,
                                    xmlTextReaderConstNamespaceUri(reader));
                                if (ns == NULL)
                                {
                                    ns = xmlNewNs (node,
                                        xmlTextReaderConstNamespaceUri(reader),
                                        xmlTextReaderConstPrefix(reader));
                                }

                                /* Create a new property tagged with this
                                   namespace and carried by this node */
                                xmlNewNsProp (current, ns,
                                    xmlTextReaderConstLocalName(reader), v);
                            }
                         }
                         else
                         {
                            /* Add the attribute as a property of the node
                               in the tree */
                            xmlNewProp (current, k, v);
                         }

                         /* Free the XML value pointer */
                         xmlFree (v);
                    }

                    /* We are done with the attributes so go to the current
                       attribute node */
                    xmlTextReaderMoveToElement (reader);
                }

                /* If this is an empty element, then return to the parent */
                if (xmlTextReaderIsEmptyElement(reader))
                    current = current->parent;
                break;
            }  /* End: Node is an element */

            case XML_READER_TYPE_END_ELEMENT:
            {  /* Node is an end element (ex. </global_metadata>, so return
                  to the parent */
                current = current->parent;
                break;
            }

            case XML_READER_TYPE_TEXT:
            {  /* Node is text or white space */
                /* Read the value of the text and add it as text for the
                   node, which is then added as a child to the tree */
                const xmlChar *v = xmlTextReaderConstValue (reader);
                xmlNodePtr node = xmlNewDocText (doc, v);
                xmlAddChild (current, node);
                break;
            }
        }  /* end switch */

        /* Read the next node */
        status = xmlTextReaderRead (reader);
    }  /* end while */
    if (status != 0)
    {
        sprintf (errmsg, "Failed to parse %s", metafile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* If the document tree is not NULL then send in the root node of the
       tree to be parsed and read into the ESPA metadata structure */
    if (doc != NULL)
    {
        /* Store the namespace for the overall metadata file */
        xmlNsPtr ns = xmlDocGetRootElement(doc)->nsDef;
        count = snprintf (metadata->meta_namespace,
            sizeof (metadata->meta_namespace), "%s", (const char *) ns->href);
        if (count < 0 || count >= sizeof (metadata->meta_namespace))
        {
            sprintf (errmsg, "Overflow of metadata->meta_namespace string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Initialize the stack to hold the elements */
        if (init_stack (&top_of_stack, &stack))
        {
            sprintf (errmsg, "Initializing the stack.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Parse the tile XML document into our ESPA internal metadata
           structure */
        if (parse_sentinel_tile_xml_into_struct (xmlDocGetRootElement(doc),
            metadata, &top_of_stack, stack, nrows, ncols))
        {
            sprintf (errmsg, "Parsing the tile metadata file into the internal "
                "metadata structure.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Clean up the XML document and the stack */
        xmlFreeDoc (doc);
        free_stack (&stack);
    }

    /* Free the reader and associated memory */
    xmlFreeTextReader (reader);
    xmlCleanupParser();
    xmlMemoryDump();

    /* The nrows/ncols need to be added to the band metadata for each of the
       bands */
    for (i = 0; i < metadata->nbands; i++)
    {
        bmeta = &metadata->band[i];
        switch (i)
        {
            /* 5 - 10m bands */
            case 1:  /* b2 */
            case 2:  /* b3 */
            case 3:  /* b4 */
            case 7:  /* b8 */
            case 13: /* TCI band */
                bmeta->nlines = nrows[0];
                bmeta->nsamps = ncols[0];
                bmeta->pixel_size[0] = 10.0;
                bmeta->pixel_size[1] = 10.0;
                strcpy (bmeta->pixel_units, "meters");
                break;

            /* 6 - 20m bands */
            case 4:  /* b5 */
            case 5:  /* b6 */
            case 6:  /* b7 */
            case 8:  /* b8a */
            case 11: /* b11 */
            case 12: /* b12 */
                bmeta->nlines = nrows[1];
                bmeta->nsamps = ncols[1];
                bmeta->pixel_size[0] = 20.0;
                bmeta->pixel_size[1] = 20.0;
                strcpy (bmeta->pixel_units, "meters");
                break;

            /* 3 - 60m bands */
            case 0:  /* b1 */
            case 9:  /* b9 */
            case 10: /* b10 */
                bmeta->nlines = nrows[2];
                bmeta->nsamps = ncols[2];
                bmeta->pixel_size[0] = 60.0;
                bmeta->pixel_size[1] = 60.0;
                strcpy (bmeta->pixel_units, "meters");
                break;
        }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  add_product_granule_metadata

PURPOSE: Add the current product granule-based information to the metadata
structure and process the children of this node.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the product granule-based metadata
SUCCESS         Successful parse of the product granule-based metadata

NOTES:
1. This MTD_MSIL1C.xml Sentinel product XML file does not provide the
   nlines, nsamps for the bands.
******************************************************************************/
int add_product_granule_metadata
(
    xmlNode *a_node,            /* I: pointer to the element node to process */
    Espa_global_meta_t *gmeta,  /* I: global metadata structure */
    int nbands,                 /* I: number of bands allocated in bmeta */
    Espa_band_meta_t *bmeta     /* I: band metadata pointer to all bands */
)
{
    char FUNC_NAME[] = "add_product_granule_metadata";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    char tmpfile[STR_SIZE];       /* temporary string for the filename */
    char *cptr = NULL;            /* pointer to image filename */
    int cur_band = 0;             /* current band being processed in the
                                     bands metadata section */
    xmlNode *cur_node = NULL;     /* pointer to the current node */
    xmlNode *child_node = NULL;   /* pointer to the child node */
    int count;                    /* number of chars copied in snprintf */

    /* Process the siblings in the Granule element */
    for (cur_node = a_node->children; cur_node;
         cur_node = xmlNextElementSibling (cur_node))
    {
        /* Look for the IMAGE_FILE elements and store them as the band name,
           which will then need to get cleaned up later */
        if (xmlStrEqual (cur_node->name, (const xmlChar *) "IMAGE_FILE"))
        {
            /* If we have more bands than were allocated, then we have an
               issue */
            if (cur_band >= nbands)
            {
                sprintf (errmsg, "Number of bands consumed already "
                    "reached the total number of bands allocated for this "
                    "XML file (%d).", nbands);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Expect the child node to be a text node containing the value
               of this field */
            child_node = cur_node->children;
            if (child_node == NULL || child_node->type != XML_TEXT_NODE)
            {
                sprintf (errmsg, "Processing Granule metadata element: %s.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Copy contents of the node into the file_name value for the
               current band */
            count = snprintf (tmpfile, sizeof (tmpfile), "%s",
                (const char *) child_node->content);
            if (count < 0 || count >= sizeof (tmpfile))
            {
                sprintf (errmsg, "Overflow of tmpfile string");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Strip the tile name of the directory/band and use that for the
               product ID; only need to do for the first band. */
            if (cur_band == 0)
            {
                strncpy (gmeta->product_id, &tmpfile[8], 34);
                gmeta->product_id[34] = '\0';
            }

            /* Grab the image filename for this band from the tmpfile string */
            cptr = strrchr (tmpfile, '/');
            if (cptr == NULL)
            {
                sprintf (errmsg, "Unsuspected format for the filename. "
                    "Expected GRANULE/L1C_{blah}/IMG_DATA/{image_filename} "
                    "however no directory separator ('/') was found in the "
                    "filename: %s.", tmpfile);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            cptr++;

            count = snprintf (bmeta[cur_band].file_name,
                sizeof (bmeta[cur_band].file_name), "%s.jp2",
                (const char *) cptr);
            if (count < 0 || count >= sizeof (bmeta[cur_band].file_name))
            {
                sprintf (errmsg, "Overflow of bmeta[%d].file_name string",
                    cur_band);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Increment the band count */
            cur_band++;
        }
    }  /* end for cur_node */

    return (SUCCESS);
}


/******************************************************************************
MODULE:  parse_sentinel_product_xml_into_struct

PURPOSE: Parse the Sentinel L1C product XML document(MTD_MSIL1C.xml) data
into the ESPA internal metadata structure.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the metadata elements
SUCCESS         Successful parse of the metadata values

NOTES:
1. Uses a stack of character strings to keep track of the nodes that have
   been parsed.  The stack must be allocated before calling this routine.
******************************************************************************/
int parse_sentinel_product_xml_into_struct
(
    xmlNode *a_node,                  /* I: pointer to the current node */
    Espa_internal_meta_t *metadata,   /* I: ESPA internal metadata structure
                                            to be filled */
    int *top_of_stack,                /* I: pointer to top of the stack */
    char **stack,                     /* I: stack to use for parsing */
    char *prodtype,                   /* O: product type string for all bands */
    float *scale_factor               /* O: scale factor for all bands */
)
{
    char FUNC_NAME[] = "parse_sentinel_product_xml_into_struct"; /* func name */
    char errmsg[STR_SIZE];       /* error message */
    char *curr_stack_element = NULL;  /* element popped from the stack */
    xmlNode *cur_node = NULL;    /* pointer to the current node */
    xmlNode *sib_node = NULL;    /* pointer to the sibling node */
    xmlNode *child_node = NULL;  /* pointer to the child node */
    static int nbands = 0;       /* number of bands in the XML structure */
    int count;                   /* number of chars copied in snprintf */
    float ul[2], ur[2];          /* UL and UR lat/long corner points */
    float ll[2], lr[2];          /* LL and LR lat/long corner points */
    bool skip_child;             /* boolean to specify the children of this
                                    node should not be processed */
    Espa_global_meta_t *gmeta = &metadata->global;
                                 /* global metadata structure */

    /* Start at the input node and traverse the tree, visiting all the children
       and siblings */
    for (cur_node = a_node; cur_node;
         cur_node = xmlNextElementSibling (cur_node))
    {
        /* Process the children of this element unless otherwise specified */
        skip_child = false;

        /* Only print the ELEMENT node types */
        if (cur_node->type == XML_ELEMENT_NODE) 
        {
            /* Push the element to the stack and turn the boolean on if this
               is either Granule metadata */
            //printf ("***Pushed %s\n", cur_node->name); fflush (stdout);
            if (push (top_of_stack, stack, (const char *) cur_node->name))
            {
                sprintf (errmsg, "Pushing element '%s' to the stack.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Process the product generation time */
            if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "GENERATION_TIME"))
            {
                /* Expect the child node to be a text node containing the
                   value of this field */
                child_node = cur_node->children;
                if (child_node == NULL || child_node->type != XML_TEXT_NODE)
                {
                    sprintf (errmsg, "Processing product metadata element: %s.",
                        cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Copy the content of the child node into the value for this
                   field */
                count = snprintf (gmeta->level1_production_date,
                    sizeof (gmeta->level1_production_date), "%s",
                    (const char *) child_node->content);
                if (count < 0 ||
                    count >= sizeof (gmeta->level1_production_date))
                {
                    sprintf (errmsg, "Overflow of gmeta->production_date "
                        "string");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
            }

            /* Process the spacecraft name */
            else if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "SPACECRAFT_NAME"))
            {
                /* Expect the child node to be a text node containing the
                   value of this field */
                child_node = cur_node->children;
                if (child_node == NULL || child_node->type != XML_TEXT_NODE)
                {
                    sprintf (errmsg, "Processing product metadata element: %s.",
                        cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Copy the content of the child node into the value for this
                   field */
                count = snprintf (gmeta->satellite, sizeof (gmeta->satellite),
                    "%s", (const char *) child_node->content);
                if (count < 0 ||
                    count >= sizeof (gmeta->level1_production_date))
                {
                    sprintf (errmsg, "Overflow of gmeta->production_date "
                        "string");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
            }

            /* Process the product type and save it for adding to the band
               metadata */
            else if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "PRODUCT_TYPE"))
            {
                /* Expect the child node to be a text node containing the
                   value of this field */
                child_node = cur_node->children;
                if (child_node == NULL || child_node->type != XML_TEXT_NODE)
                {
                    sprintf (errmsg, "Processing product metadata element: %s.",
                        cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Copy the content of the child node into the value for this
                   field */
                count = snprintf (prodtype, sizeof (prodtype), "%s",
                    (const char *) child_node->content);
                if (count < 0 || count >= sizeof (prodtype))
                {
                    sprintf (errmsg, "Overflow of prodtype string");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
            }

            /* Process the quantification value and store it as the scale
               factor for each band */
            else if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "QUANTIFICATION_VALUE"))
            {
                /* Expect the child node to be a text node containing the
                   value of this field */
                child_node = cur_node->children;
                if (child_node == NULL || child_node->type != XML_TEXT_NODE)
                {
                    sprintf (errmsg, "Error processing product metadata "
                        "element: %s", cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
    
                /* Copy the content of the child node into the value for this
                   field */
                *scale_factor = atof ((const char *) child_node->content);
            }

            /* Process the global footprint and store it as the corner UL/LR
               as well as the bounding coordinates */
            else if (xmlStrEqual (cur_node->name,
                (const xmlChar *) "EXT_POS_LIST"))
            {
                /* Expect the child node to be a text node containing the
                   value of this field */
                child_node = cur_node->children;
                if (child_node == NULL || child_node->type != XML_TEXT_NODE)
                {
                    sprintf (errmsg, "Error processing product metadata "
                        "element: %s", cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Copy the content of the child node into the value for this
                   field. According to the documentation, the footprint is a
                   closed polygon provided as a series of vertices (lat, lon)
                   counter-clockwise oriented. The last/fifth point can be
                   skipped in this case since it is a duplication of the UL
                   in the closed system. NOTE: After drawing out the corner
                   points in this XML metadata, they instead appear to be in
                   a clockwise direction, contradictory to the documentation. */
                sscanf ((const char *) child_node->content,
                    "%f %f %f %f %f %f %f %f", &ul[0], &ul[1], &ur[0], &ur[1],
                    &lr[0], &lr[1], &ll[0], &ll[1]);

                /* Assign the UL and LR lat/long values */
                gmeta->ul_corner[0] = ul[0];
                gmeta->ul_corner[1] = ul[1];
                gmeta->lr_corner[0] = lr[0];
                gmeta->lr_corner[1] = lr[1];

                /* Use the corners to determine the bounding coordinates */
                if (ul[1] < ll[1])
                    gmeta->bounding_coords[ESPA_WEST] = ul[1];
                else
                    gmeta->bounding_coords[ESPA_WEST] = ll[1];

                if (ur[1] > lr[1])
                    gmeta->bounding_coords[ESPA_EAST] = ur[1];
                else
                    gmeta->bounding_coords[ESPA_EAST] = lr[1];

                if (ul[0] > ur[0])
                    gmeta->bounding_coords[ESPA_NORTH] = ul[0];
                else
                    gmeta->bounding_coords[ESPA_NORTH] = ur[0];

                if (ll[0] < lr[0])
                    gmeta->bounding_coords[ESPA_SOUTH] = ll[0];
                else
                    gmeta->bounding_coords[ESPA_SOUTH] = lr[0];
            }

            /* Process the granule metadata. Count the number of band elements
               in this structure, then allocate memory for the nbands. */
            else if (xmlStrEqual (cur_node->name, (const xmlChar *) "Granule"))
            {
                /* Count the number of siblings which are band elements */
                nbands = 0;
                for (sib_node = cur_node->children; sib_node;
                     sib_node = xmlNextElementSibling (sib_node))
                {
                    /* If this is an IMAGE_FILE element then count it */
                    if (xmlStrEqual (sib_node->name,
                        (const xmlChar *) "IMAGE_FILE"))
                        nbands++;
                }

                /* Allocate band metadata */
                if (allocate_band_metadata (metadata, nbands) != SUCCESS)
                {   /* Error messages already printed */
                    return (ERROR);
                }

                /* Add the image filenames to the band metadata */
                if (add_product_granule_metadata (cur_node, gmeta, nbands,
                    metadata->band))
                {
                    sprintf (errmsg, "Consuming Granule image file elements "
                        "'%s'.", cur_node->name);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Skip processing the children of this node, since they
                   will be handled by the global metadata parser */
                skip_child = true;
            }

#ifdef DEBUG
            /* Print out the name of the element */
            xmlAttrPtr attr;     /* pointer to the element attributes */
            printf ("node type: Element, name: %s\n", cur_node->name);

            /* Print out the attribute properties for this element */
            for (attr = cur_node->properties; attr != NULL; attr = attr->next)
            {
                xmlChar *v = xmlGetProp (cur_node, attr->name);
                printf (" @%s=%s ", attr->name, v);
                xmlFree (v);
            }
            printf ("\n"); fflush (stdout);
#endif
        }
        else if (cur_node->type == XML_TEXT_NODE) 
        {
#ifdef DEBUG
            /* Print out the text for the element */
            printf ("   node type: Text, content: %s\n", cur_node->content);
#endif
        }

        /* Parse the children of this node if they haven't been consumed
           elsewhere */
        if (!skip_child)
        {
            if (parse_sentinel_product_xml_into_struct (cur_node->children,
                metadata, top_of_stack, stack, prodtype, scale_factor))
            {
                sprintf (errmsg, "Parsing the children of this element '%s'.",
                    cur_node->name);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* Done with the element and its siblings so pop the element name off
           the stack */
        if (cur_node->type == XML_ELEMENT_NODE)
        {
            curr_stack_element = pop (top_of_stack, stack);
            if (curr_stack_element == NULL)
            {
                sprintf (errmsg, "Popping elements off the stack.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            //printf ("***Popped %s\n", curr_stack_element); fflush (stdout);
        }
    }  /* for cur_node */

    return (SUCCESS);
}


/******************************************************************************
MODULE:  parse_sentinel_product_metadata

PURPOSE: Parse the Sentinel L1C product metadata file (MTD_MSIL1C.xml) and
populate the associated ESPA internal metadata file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error parsing the Sentinel metadata elements
SUCCESS         Successful parse of the Sentinel metadata values

NOTES:
1. Uses a stack of character strings to keep track of the nodes that have been
   found in the metadata document.
2. For debugging purposes
   xmlDocDump (stderr, doc);
   can be used to dump/print the XML doc to the screen.
******************************************************************************/
int parse_sentinel_product_metadata
(
    char *metafile,                 /* I: Sentinel product metadata file */
    Espa_internal_meta_t *metadata  /* I: input metadata structure which has
                                          been initialized via
                                          init_metadata_struct */
)
{
    char FUNC_NAME[] = "parse_sentinel_product_metadata";  /* function name */
    char errmsg[STR_SIZE];    /* error message */
    char prodtype[STR_SIZE];  /* product type string for all bands */
    xmlTextReaderPtr reader;  /* reader for the XML file */
    xmlDocPtr doc = NULL;     /* document tree pointer */
    xmlNodePtr current=NULL;  /* pointer to the current node */
    int i;                    /* looping variable */
    int status;               /* return status */
    int nodeType;             /* node type (element, text, attribute, etc.) */
    int top_of_stack;         /* top of the stack */
    int count;                /* number of chars copied in snprintf */
    float scale_factor;       /* scale factor for all bands */
    char **stack = NULL;      /* stack to keep track of elements in the tree */
    Espa_band_meta_t *bmeta;  /* band metadata pointer to all bands */

    /* Establish the reader for this metadata file */
    reader = xmlNewTextReaderFilename (metafile);
    if (reader == NULL)
    {
        sprintf (errmsg, "Setting up reader for %s", metafile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Use the reader to parse the XML file, looking at each of the nodes,
       until the entire file has been parsed.  Start by reading the first
       node in the file. */
    status = xmlTextReaderRead (reader);
    while (status == 1)
    {
        /* Determine what kind of node the reader is at (element, end element,
           attribute, text/white space) and handle the information as desired */
        nodeType = xmlTextReaderNodeType (reader);
        if (nodeType == -1)
        {
            sprintf (errmsg, "Getting node type");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        switch (nodeType)
        {
            case XML_READER_TYPE_ELEMENT:
            {  /* Node is an element (ex. <global_metadata>) */
                xmlNodePtr node=NULL;
                if (doc==NULL)
                {
                    doc=xmlNewDoc (BAD_CAST "1.0");
                }

                /* Get the URI defining the namespace for this node */
                if (xmlTextReaderConstNamespaceUri (reader) != NULL)
                {
                    /* Search the namespace for a document */
                    xmlNsPtr ns = xmlSearchNs (doc, current,
                        xmlTextReaderConstNamespaceUri(reader));

                    /* Create a tree node for this element in the XML file
                       using the element name */
                    node = xmlNewNode (ns, xmlTextReaderConstName(reader));

                    /* If the namespace is empty (i.e. root) then create a
                       new namespace pointer with this node */
                    if (ns == NULL)
                    {
                        ns = xmlNewNs (node,
                            xmlTextReaderConstNamespaceUri(reader),
                            xmlTextReaderConstPrefix(reader));
                    }
                }
                else
                {
                    /* Create a tree node for this element in the XML file
                       using the element name */
                    node = xmlNewNode (0, xmlTextReaderConstName(reader));
                }

                /* Set the element as the root if appropriate otherwise add
                   it as a child to the previous element */
                if (current == NULL)
                {
                    xmlDocSetRootElement (doc, node);
                }
                else
                {
                    xmlAddChild (current, node);
                }
                current = node;

                /* If the element has attributes, then handle them */
                if (xmlTextReaderHasAttributes (reader))
                {
                    /* Get the number of attributes and then process each one */
                    int i;
                    int n_att = xmlTextReaderAttributeCount (reader);
                    for (i = 0; i < n_att; i++)
                    {
                        /* Read each attribute, obtain the name and value,
                           then add it as a property for this node in the
                           tree */
                        const xmlChar *k = NULL;
                        xmlChar *v = NULL;
                        xmlTextReaderMoveToAttributeNo (reader, i);
                        k = xmlTextReaderConstName (reader);
                        v = xmlTextReaderValue (reader);
                        if (xmlTextReaderConstNamespaceUri (reader) != NULL)
                        {
                            if (!xmlStrEqual (
                                xmlTextReaderConstNamespaceUri(reader),
                                BAD_CAST "http://www.w3.org/2000/xmlns/"))
                            {
                                /* Search the namespace for the document */
                                xmlNsPtr ns = xmlSearchNs (doc, current,
                                    xmlTextReaderConstNamespaceUri(reader));
                                if (ns == NULL)
                                {
                                    ns = xmlNewNs (node,
                                        xmlTextReaderConstNamespaceUri(reader),
                                        xmlTextReaderConstPrefix(reader));
                                }

                                /* Create a new property tagged with this
                                   namespace and carried by this node */
                                xmlNewNsProp (current, ns,
                                    xmlTextReaderConstLocalName(reader), v);
                            }
                         }
                         else
                         {
                            /* Add the attribute as a property of the node
                               in the tree */
                            xmlNewProp (current, k, v);
                         }

                         /* Free the XML value pointer */
                         xmlFree (v);
                    }

                    /* We are done with the attributes so go to the current
                       attribute node */
                    xmlTextReaderMoveToElement (reader);
                }

                /* If this is an empty element, then return to the parent */
                if (xmlTextReaderIsEmptyElement(reader))
                    current = current->parent;
                break;
            }  /* End: Node is an element */

            case XML_READER_TYPE_END_ELEMENT:
            {  /* Node is an end element (ex. </global_metadata>, so return
                  to the parent */
                current = current->parent;
                break;
            }

            case XML_READER_TYPE_TEXT:
            {  /* Node is text or white space */
                /* Read the value of the text and add it as text for the
                   node, which is then added as a child to the tree */
                const xmlChar *v = xmlTextReaderConstValue (reader);
                xmlNodePtr node = xmlNewDocText (doc, v);
                xmlAddChild (current, node);
                break;
            }
        }  /* end switch */

        /* Read the next node */
        status = xmlTextReaderRead (reader);
    }  /* end while */
    if (status != 0)
    {
        sprintf (errmsg, "Failed to parse %s", metafile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* If the document tree is not NULL then send in the root node of the
       tree to be parsed and read into the ESPA metadata structure */
    if (doc != NULL)
    {
        /* Store the namespace for the overall metadata file */
        xmlNsPtr ns = xmlDocGetRootElement(doc)->nsDef;
        count = snprintf (metadata->meta_namespace,
            sizeof (metadata->meta_namespace), "%s", (const char *) ns->href);
        if (count < 0 || count >= sizeof (metadata->meta_namespace))
        {
            sprintf (errmsg, "Overflow of metadata->meta_namespace string");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Initialize the stack to hold the elements */
        if (init_stack (&top_of_stack, &stack))
        {
            sprintf (errmsg, "Initializing the stack.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        //print_element_names (xmlDocGetRootElement (doc));

        /* Parse the XML document into our ESPA internal metadata structure */
        if (parse_sentinel_product_xml_into_struct (xmlDocGetRootElement(doc),
            metadata, &top_of_stack, stack, prodtype, &scale_factor))
        {
            sprintf (errmsg, "Parsing the product metadata file into the "
                "internal metadata structure.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Clean up the XML document and the stack */
        xmlFreeDoc (doc);
        free_stack (&stack);
    }

    /* Free the reader and associated memory */
    xmlFreeTextReader (reader);
    xmlCleanupParser();
    xmlMemoryDump();

    /* The product type and scale factor need to be added to the band metadata
       for each of the bands */
    for (i = 0; i < metadata->nbands; i++)
    {
        bmeta = &metadata->band[i];
        strcpy (bmeta->short_name, prodtype);
        bmeta->scale_factor = 1.0 / scale_factor;
    }

    return (SUCCESS);
}

