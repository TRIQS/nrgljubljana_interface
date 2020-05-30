# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=nrgljubljana_interface

COPY requirements.txt /src/$APPNAME/requirements.txt
RUN pip3 install -r /src/$APPNAME/requirements.txt

RUN apt-get install -y libgsl-dev || yum install -y gsl-devel

ENV BOOST_ROOT=/opt/boost
ADD https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz /tmp/boost.tar.gz
RUN chown build /tmp/boost.tar.gz
USER build
RUN tar -C /tmp -xf /tmp/boost.tar.gz && \
    cd /tmp/boost_* && \
    case $CC in (clang*) toolset=clang ;; (gcc*) toolset=gcc ;; esac ; \
    ./bootstrap.sh --prefix=$BOOST_ROOT --with-toolset=$toolset --with-libraries=serialization && \
    ./b2 ${CXXFLAGS:+cxxflags=$CXXFLAGS linkflags=$CXXFLAGS}
USER root
RUN cd /tmp/boost_* && ./b2 install && \
    cd / && rm -rf /tmp/boost*

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_DOC=0
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} -DBuild_Documentation=${BUILD_DOC} && make -j2
USER root
RUN make install
