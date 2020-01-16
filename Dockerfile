# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=nrgljubljana_interface

COPY requirements.txt /src/$APPNAME/requirements.txt
RUN pip install -r /src/$APPNAME/requirements.txt

RUN apt-get install -y libgsl-dev || yum install -y gsl-devel

ENV BOOST_VERSION=1.72.0 BOOST_ROOT=/opt/boost
ADD https://github.com/boostorg/boost/archive/boost-$BOOST_VERSION.tar.gz /tmp/boost-boost.tar.gz
ADD https://github.com/boostorg/build/archive/boost-$BOOST_VERSION.tar.gz /tmp/boost-build.tar.gz
ADD https://github.com/boostorg/config/archive/boost-$BOOST_VERSION.tar.gz /tmp/boost-config.tar.gz
ADD https://github.com/boostorg/boost_install/archive/boost-$BOOST_VERSION.tar.gz /tmp/boost-boost_install.tar.gz
ADD https://github.com/boostorg/headers/archive/boost-$BOOST_VERSION.tar.gz /tmp/boost-headers.tar.gz
ADD https://github.com/boostorg/core/archive/boost-$BOOST_VERSION.tar.gz /tmp/boost-core.tar.gz
ADD https://github.com/boostorg/serialization/archive/boost-$BOOST_VERSION.tar.gz /tmp/boost-serialization.tar.gz
RUN chown build /tmp/boost-*.tar.gz
USER build
RUN mkdir -p /tmp/boost && \
    cd /tmp/boost && \
    for c in \
      boost:. \
      build:tools/build \
      config:libs/config \
      boost_install:tools/boost_install \
      headers:libs/headers \
      core:libs/core \
      serialization:libs/serialization \
    ; do \
      # https://github.com/boostorg/${c%:*}/archive/boost-$BOOST_VERSION.tar.gz
      tar -C ${c#*:} --strip-components=1 -xf /tmp/boost-${c%:*}.tar.gz && rm /tmp/boost-${c%:*}.tar.gz ; \
    done && \
    case $CC in (clang*) toolset=clang ;; (gcc*) toolset=gcc ;; esac ; \
    ./bootstrap.sh --prefix=$BOOST_ROOT --with-toolset=$toolset --with-libraries=serialization && \
    ./b2 cxxflags=$CXXFLAGS linkflags=$CXXFLAGS
USER root
RUN cd /tmp/boost && ./b2 install && \
    cd / && rm -rf /tmp/boost

COPY . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown -R build $SRC/$APPNAME .
USER build
ARG BUILD_DOC=0
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} -DBuild_Documentation=${BUILD_DOC} && make -j2 && make test CTEST_OUTPUT_ON_FAILURE=1
USER root
RUN make install
