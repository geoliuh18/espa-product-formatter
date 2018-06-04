# =======================================================
#   Binary Build Layer
# =======================================================
FROM usgseros/espa-processing:docker-devel-3.0rc1.dev1 as builder
LABEL maintainer="USGS EROS LSRD http://eros.usgs.gov" \
      description="ESPA Product Formatting Software"

ENV PREFIX=/usr/local \
    ESPAINC=/usr/local/include \
    ESPALIB=/usr/local/lib \
    ESPA_XML_SCHEMA=/usr/local/schema/espa_internal_metadata_v2_0.xsd \
    ESPA_LAND_MASS_POLYGON=/usr/local/auxiliaries/land_water_polygon/land_no_buf.ply

WORKDIR ${SRC_DIR}
COPY . ${SRC_DIR}

# `  Install the product-formatter applications
RUN cd ${SRC_DIR}/raw_binary \
    && make BUILD_STATIC=yes ENABLE_THREADING=yes \
    && make install \
    && cd ${SRC_DIR}/schema \
    && make install \
    && cd ${SRC_DIR} \
    && rm -rf *

