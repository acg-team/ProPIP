[Back](./Index.md) | [Home](../ProPIP-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

---
# Updating ProPIP executable
---

The updating process is divided into two steps. When updating check first the branch to checkout on both TSHLib and ProPIP repositories and substitute the `<branch>` tag with an appropriate one (i.e. `master`).


### Update TSHLib

```
#!bash

cd /source/folder/tshlib
git checkout <branch>
git pull
cmake -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt; sudo make install;
```

### Update miniJATI

```
#!bash

cd /source/folder/castor
git checkout <branch>
git pull
cmake -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt; make;

```
