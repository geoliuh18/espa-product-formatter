# =======================================================
#   Binary Build Layer
# =======================================================
FROM jbrinkmann/lagoon-dragonfly:devel-1.0.0 as builder
LABEL maintainer="USGS EROS LSRD http://eros.usgs.gov" \
      description="ESPA Product Formatting Software"

ENV DIR_RAW_BINARY=raw_binary \
    DIR_SCHEMA=schema \
    BUILD_PREFIX=/usr/local

WORKDIR ${SRC_DIR}
COPY . ${SRC_DIR}

# `  Install the product-formatter applications
RUN cd ${SRC_DIR}/${DIR_RAW_BINARY} \
    && make BUILD_STATIC=yes ENABLE_THREADING=yes \
    && make install PREFIX=${BUILD_PREFIX} \
    && cd ${SRC_DIR}/${DIR_SCHEMA} \
    && make install PREFIX=${BUILD_PREFIX} \
    && cd ${SRC_DIR} \
    && rm -rf *
