# =======================================================
#   Binary Build Layer
# =======================================================
FROM jbrinkmann/lagoon-dragonfly:devel-1.0.0 as builder
LABEL maintainer="USGS EROS LSRD http://eros.usgs.gov" \
      description="ESPA Product Formatting Software"

ENV PREFIX=/usr/local \
    ESPAINC=/usr/local/include \
    ESPALIB=/usr/local/lib

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
