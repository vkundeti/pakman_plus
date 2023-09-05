# `pakman-mpi:` Improved PakMan Algorithm with Concurrent k-mer Counting 
## Summary
`pakman-mpi` is a PakMan executable with improved k-mer counting algorithm for better performance 

## Dependencies
Following are some dependencies to build `pakman-mpi`:
1. `oneTBB`: install [oneTBB](https://github.com/oneapi-src/oneTBB) and set `TBB_DIR` to the
directory inside the install directory pointing to `TBBConfig.cmake` (e.g. `oneTBB/build/install/lib/cmake/TBB`), set this to an environment variable `$TBB_DIR`.
Following is one example on how to build 

```
$git clone https://github.com/oneapi-src/oneTBB
$cd oneTBB
$mkdir build
$cd build
$cmake ../ -DCMAKE_INSTALL_PREFIX=./install/
$make -j10; make install
$export TBB_DIR=`pwd`/install/lib/cmake/TBB/
```

### Optional 
1. `icpx` Intel compiler: please use Intel compiler to avoid incompatibility issues with OMP (E.g. `source /opt/intel/oneapi/compiler/2021.1.1/env/vars.sh`)

2. `gtest`: Google unit test library [gtest](https://github.com/google/googletest).

## Build `pakman-mpi-exe` 

```
$ git clone git@github.com:intel-sandbox/pakman_intel_esc.git
$ cd pakman_intel_esc/
$ mkdir build
$ cd build
$ cmake ../src cmake ../src/ -DTBB_DIR=$TBB_DIR 
$ make -j10
```

### Build Unit Tests
If you want to build unit tests change the `cmake` command to the following: 
```
$cmake ../src/ -DTBB_DIR=$TBB_DIR -DBUILD_UNIT_TESTS=True
$make -j10
```
You should see `pakman-unit-tests-exe`


### Running `pakman-mpi`:
The executable `pakman-mpi` supports the same arguments as the original `pakman`
program. 
