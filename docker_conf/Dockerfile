FROM ubuntu:16.04

# Update base image ubuntu
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get upgrade -y
RUN apt-get install -y build-essential
RUN apt-get install -y git
RUN apt-get install -y wget
RUN apt-get install -y autoconf
RUN apt-get install -y gfortran
RUN apt-get install -y libblas-dev liblapack-dev  libpthread-stubs0-dev
RUN apt-get install -y pkg-config
RUN apt-get install -y libtool
RUN apt-get install -y libboost-all-dev

# Get dependencies
RUN wget https://cmake.org/files/v3.10/cmake-3.10.3-Linux-x86_64.sh; chmod +x cmake-3.10.3-Linux-x86_64.sh; ./cmake-3.10.3-Linux-x86_64.sh --skip-license; ln -s ./cmake-3.10.3-Linux-x86_64/bin/cmake /usr/bin/cmake;
RUN git clone https://github.com/BioPP/bpp-core; cd bpp-core; git checkout tags/v2.4.0 -b v240; mkdir build; cd build; cmake ..; make install; cd ../.. ; rm -r bpp-core;
RUN git clone https://github.com/BioPP/bpp-seq; cd bpp-seq; git checkout tags/v2.4.0 -b v240; mkdir build; cd build; cmake ..; make install; cd ../.. ; rm -r bpp-seq;
RUN git clone https://github.com/BioPP/bpp-phyl; cd bpp-phyl; git checkout tags/v2.4.0 -b v240; mkdir build; cd build; cmake ..; make install; cd ../.. ; rm -r bpp-phyl;
RUN git clone https://github.com/google/glog; cd glog; cmake -H. -Bbuild -G "Unix Makefiles"; cmake --build build --target install;
RUN git clone https://www.bitbucket.org/lorenzogatti89/tshlib.git; cd tshlib; git checkout master; git pull; cmake -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt; make install; cd .. ;
