#!/bin/bash

MPATH=.

export LD_LIBRARY_PATH=${MPATH}/visit/vtk/6.1.0/linux-x86_64_gcc-4.8/lib:${MPATH}/OSPRay/OSPRay-1.1.0/release:${MPATH}/visit/hdf5/1.8.14/linux-x86_64_gcc-4.8/lib:${MPATH}/visit/netcdf/4.1.1/linux-x86_64_gcc-4.8/lib:${MPATH}/visit/python/2.7.6/linux-x86_64_gcc-4.8/lib:${MPATH}/visit/qt/4.8.3/linux-x86_64_gcc-4.8/lib:${MPATH}/visit/silo/4.10.1/linux-x86_64_gcc-4.8/lib:${MPATH}/visit/szip/2.1/linux-x86_64_gcc-4.8/lib:${LD_LIBRARY_PATH}

export PATH=${MPATH}/visit/vtk/6.1.0/linux-x86_64_gcc-4.8/bin:${MPATH}/OSPRay/OSPRay-1.1.0/release:${MPATH}/visit/hdf5/1.8.14/linux-x86_64_gcc-4.8/bin:${MPATH}/visit/netcdf/4.1.1/linux-x86_64_gcc-4.8/bin:${MPATH}/visit/python/2.7.6/linux-x86_64_gcc-4.8/bin:${MPATH}/visit/qt/4.8.3/linux-x86_64_gcc-4.8/bin:${MPATH}/visit/silo/4.10.1/linux-x86_64_gcc-4.8/bin:${MPATH}/visit/szip/2.1/linux-x86_64_gcc-4.8/bin:${PATH}
