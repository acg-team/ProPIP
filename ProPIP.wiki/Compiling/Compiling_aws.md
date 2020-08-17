[Back](./Index.md) | [Home](https://github.com/acg-team/ProPIP/blob/master/ProPIP.wiki/ProPIP-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

---
# Compilation on AWS environment with EFS
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
SharedLibraryPath=/mnt/efs/projects/local
SharedLocalPath=/mnt/efs/projects/local/include
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
bash Miniconda-latest-Linux-x86_64.sh -b -p /mnt/efs/projects/bin/anaconda_ete
export PATH=/mnt/efs/projects/bin/anaconda_ete/bin:$PATH;
conda install anaconda-client
conda install -c etetoolkit ete3 ete_toolchain
xvfb-run ete3 build check
```



### Configuration AWS with load balancer with file setting


on AWS load balancer configure the slaves as follows:

```
#cloud-config
package_upgrade: true
packages:
- nfs-utils
- httpd
- php
- Xvfb
- xcb-util
- libXrender
- fontconfig
- urw-fonts
runcmd:
- echo "fs-f4b441ad.efs.eu-central-1.amazonaws.com:/    /var/www/html/efs   nfs4    defaults,acl" >> /etc/fstab
- mkdir /var/www/html/efs
- mount -a
- chmod go+rw /var/www/html/efs
- chmod 777 /var/www/html/efs
- touch /var/www/html/efs/test2.html
- service httpd start
- chkconfig httpd on
- export PATH=/mnt/efs/projects/bin/anaconda_ete/bin:/mnt/efs/projects/source/cmake-3.10.2-Linux-x86_64/bin:$PATH;
- export CPATH=/mnt/efs/projects/local/include;
- export LD_LIBRARY_PATH=/mnt/efs/projects/local/lib
```

### Configuration AWS with cnfcluster 

We use Ubuntu images for both computing nodes and login node (the cluster is dynamically scaled).

```

#!/bin/bash

sudo apt-get update
sudo apt-get -y upgrade
sudo apt-get install -y nfs-common
sudo apt-get install -y xvfb
sudo apt-get install -y xcb-util ***
sudo apt-get install -y apache2
sudo apt-get install -y php libapache2-mod-php php-mcrypt php-mysql
sudo apt-get install -y libxrender
sudo apt-get install -y fontconfig
sudo apt-get install -y fonts-texgyre
sudo apt-get install -y python3-pip
sudo apt-get install -y build-essential libssl-dev libffi-dev python3-dev
sudo apt-get install -y git
sudo apt-get install -y wget
sudo apt-get install -y autoconf
sudo apt-get install -y gfortran
sudo apt-get install -y libblas-dev liblapack-dev libpthread-stubs0-dev
sudo apt-get install -y pkg-config
sudo apt-get install -y libeigen3-dev
sudo apt-get install -y libtool
sudo apt-get install -y libboost-all-dev
# Mount permanent FileSystem
echo "fs-f4b441ad.efs.eu-central-1.amazonaws.com:/    /mnt/efs   nfs4    defaults,acl" >> /etc/fstab
mkdir /mnt/efs
mount -a
chmod go+rw /mnt/efs
chmod 777 /mnt/efs

# ---------------
#either this
export PATH=/mnt/efs/projects/source/cmake-3.10.2-Linux-x86_64/bin:$PATH;
# or this
cd /mnt/efs/projects/source/
wget https://cmake.org/files/v3.10/cmake-3.10.3-Linux-x86_64.sh; chmod +x cmake-3.10.3-Linux-x86_64.sh; ./cmake-3.10.3-Linux-x86_64.sh --skip-license; ln -s ./cmake-3.10.3-Linux-x86_64/bin/cmake /usr/bin/cmake;
# ---------------
# Bpp libraries
cd /mnt/efs/projects/source/bpp/
git clone https://github.com/BioPP/bpp-core; cd bpp-core; git checkout tags/v2.4.0 -b v240; mkdir build; cd build; cmake ..; sudo make install; cd ../..;
git clone https://github.com/BioPP/bpp-seq; cd bpp-seq; git checkout tags/v2.4.0 -b v240; mkdir build; cd build; cmake ..; sudo make install; cd ../..;
git clone https://github.com/BioPP/bpp-phyl; cd bpp-phyl; git checkout tags/v2.4.0 -b v240; mkdir build; cd build; cmake ..; sudo make install; cd ../..;
# ---------------
# GLOG library
cd /mnt/efs/projects/source/
git clone https://github.com/google/glog; cd glog; cmake -H. -Bbuild -G "Unix Makefiles"; sudo cmake --build build --target install;
# ---------------
# TSHLIB library
cd /mnt/efs/projects/source/
git clone https://${BITBUCKET_USERNAME}:${BITBUCKET_PASSWORD}@bitbucket.org/acg-team/tshlib.git; cd tshlib; git checkout master; cmake -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt; sudo make install;
# ---------------
# ProPIP executable
cd /mnt/efs/projects/source/
git clone https://${BITBUCKET_USERNAME}:${BITBUCKET_PASSWORD}@bitbucket.org/acg-team/minijati.git; cd minijati; git checkout master; cmake -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt; make;


```


### Update on AWS with shared environment policies


```
#!/bin/bash

export SharedLibraryPath=/mnt/efs/projects/local
export SharedIncludePath=${SharedLibraryPath}/include
export PATH=/mnt/efs/projects/source/cmake-3.10.2-Linux-x86_64/bin:$PATH;
cd /mnt/efs/projects/tshlib/
cmake -- -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} -DCMAKE_PREFIX_PATH=${SharedLibraryPath} CMakeLists.txt
sudo make install
cd /mnt/efs/projects/ProPIP/
cmake --target ProPIP -- -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${SharedLibraryPath} CMakeLists.txt
```
