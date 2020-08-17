[Back](./Index.md) | [Home](../ProPIP/-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

---
# Compilation on a shared environment
---

### Preparing the environment


We define a unique location for installing the libraries (lib + include). This location should be accessible with read/write permissions by the user.


```
#!bash
export SharedLibraryPath=path/to/folder
export SharedIncludePath=${SharedLibraryPath}/include

```

For example:


```
SharedLibraryPath=$HOME/local
SharedIncludePath=$HOME/local/include
```




### Compiling and installing the dependencies


**bpp-core** http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-core
cd bpp-core
git checkout tags/v2.4.0 -b v240
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} ..
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
cmake -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} ..
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
cmake -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} ..
make install
```

**boost - C++ Libraries** http://www.boost.org/

```
#!bash
wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
tar xvf /path/to/boost_1_66_0.tar.gz
cd boost_1_66_0
./bootstrap.sh --libdir=${SharedLibraryPath}/lib --includedir=${SharedIncludePath}
./b2 --libdir=${SharedLibraryPath}/lib --includedir=${SharedIncludePath}
./b2 install --prefix=${SharedLibraryPath}
```

**glog - Google Logging Library** https://github.com/google/glog

```
#!bash
git clone https://github.com/google/glog
cd glog
cmake -H. -Bbuild -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath}
cmake --build build --target install

```


**TSHLib - Tree Search Heuristics Library**

```
#!bash
git clone https://{username}@bitbucket.org/acg-team/tshlib.git
cd tshlib
git checkout develop
cmake -- -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} -DCMAKE_PREFIX_PATH=${SharedLibraryPath} CMakeLists.txt
make install
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
git clone https://github.com/acg-team/ProPIP/
cd castor
cmake --target ProPIP -- -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${SharedIncludePath} -DCMAKE_PREFIX_PATH=${SharedLibraryPath} CMakeLists.txt
make
```
*Static linking*
```
#!bash
git clone https://github.com/acg-team/ProPIP/
cd castor
cmake --target ProPIP -- -DCMAKE_BUILD_TYPE=Release-static -DCMAKE_PREFIX_PATH=${SharedIncludePath} -DCMAKE_PREFIX_PATH=${SharedLibraryPath} CMakeLists.txt
make
```
