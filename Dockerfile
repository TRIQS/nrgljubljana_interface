# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=nrgljubljana_interface

RUN apt-get install -y libgsl-dev || yum install -y gsl-devel perl-Math-Complex

ENV BOOST_ROOT=/opt/boost
ENV LD_LIBRARY_PATH=/opt/boost/lib:$LD_LIBRARY_PATH
ADD https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz /tmp/boost.tar.gz
RUN chown build /tmp/boost.tar.gz
USER build
RUN tar -C /tmp -xf /tmp/boost.tar.gz && \
    cd /tmp/boost_* && \
    case $CC in (clang*) toolset=clang ;; (gcc*) toolset=gcc ;; (icx*) toolset=intel-linux ;; esac ; \
    ./bootstrap.sh --prefix=$BOOST_ROOT --with-toolset=$toolset --with-libraries=serialization && \
    ./b2 ${CXXFLAGS:+cxxflags=$CXXFLAGS linkflags=$CXXFLAGS} toolset=$toolset
USER root
RUN cd /tmp/boost_* && \
    case $CC in (clang*) toolset=clang ;; (gcc*) toolset=gcc ;; (icx*) toolset=intel-linux ;; esac ; ./b2 toolset=$toolset install && \
    cd / && rm -rf /tmp/boost*

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_ID
ARG CMAKE_ARGS
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} $CMAKE_ARGS && make -j4 || make -j1 VERBOSE=1
USER root
RUN make install
