[Back](./Index) | [Home](../Home)

---
# Compilation on a local environment
---

The user shoud have reading/writing rights on the system folders (i.e. /usr/local, /usr/local/include).


## Compiling and installing the dependencies


**bpp-core** http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-core
cd bpp-core
git checkout tags/v2.4.0 -b v240
mkdir build
cd build
cmake ..
make install
```

**bpp-seq** http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-seq
cd bpp-seq
git checkout tags/v2.4.0 -b v240
mkdir build
cd build
cmake ..
make install
```

**bpp-phyl**  http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-phyl
cd bpp-phyl
git checkout tags/v2.4.0 -b v240
mkdir build
cd build
cmake  ..
make install
```

**boost - C++ Libraries** http://www.boost.org/

```
#!bash
wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
tar xvf /path/to/boost_1_66_0.tar.gz
cd boost_1_66_0
./bootstrap.sh
./b2
./b2 install
```

**glog - Google Logging Library** https://github.com/google/glog

```
#!bash
git clone https://github.com/google/glog
cd glog
cmake -H. -Bbuild -G "Unix Makefiles"
cmake --build build --target install

```


**TSHLib - Tree Search Heuristics Library**

```
#!bash
git clone https://{username}@bitbucket.org/acg-team/tshlib.git
cd tshlib
git checkout develop
cmake -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt
make install
```


**Intel TBB - Intel(R) Threading Building Blocks 2018**

Under MacOS

```
#!bash

cd /opt
wget https://github.com/01org/tbb/releases/download/2018_U5/tbb2018_20180618oss_mac.tgz
tar -xvf tbb2018_20180618oss_mac.tgz
```

Under Linux

```
#!bash

cd /opt
wget https://github.com/01org/tbb/releases/download/2018_U5/tbb2018_20180618oss_lin.tgz
tar -xvf tbb2018_20180618oss_lin.tgz
```

**ete3  (optional for plots)** http://etetoolkit.org/

```
#!bash
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O Miniconda-latest-Linux-x86_64.sh
bash Miniconda-latest-Linux-x86_64.sh -b -p /var/www/html/efs/projects/bin/anaconda_ete
export PATH=/var/www/html/efs/projects/bin/anaconda_ete/bin:$PATH;
conda install anaconda-client
conda install -c etetoolkit ete3 ete_toolchain
xvfb-run ete3 build check
```

### Compiling Castor


*Dynamic linking*
```
#!bash
git clone https://bitbucket.org/lorenzogatti89/castor/
cd castor
cmake --target castor -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt
make
```

*Static linking* (executable is about 140mb)
```
#!bash
git clone https://bitbucket.org/lorenzogatti89/castor/
cd castor
cmake --target castor -- -DCMAKE_BUILD_TYPE=Release-static CMakeLists.txt
make
```


If *Intel TBB* is required, then cmake must be invoked as follows:

*Dynamic linking*
```
#!bash
git clone https://bitbucket.org/lorenzogatti89/castor/
cd castor
cmake --target castor -- -DCMAKE_BUILD_TYPE=Release -DPATH_INTELTBB=/opt/<path-to-tbb> CMakeLists.txt
make
```

*Static linking* (executable is about 140mb)
```
#!bash
git clone https://bitbucket.org/lorenzogatti89/castor/
cd castor
cmake --target castor -- -DCMAKE_BUILD_TYPE=Release-static -DPATH_INTELTBB=/opt/<path-to-tbb> CMakeLists.txt
make
```
