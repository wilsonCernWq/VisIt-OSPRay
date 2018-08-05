#/Users/qwu/Work/projects/visit/3rdparty/installs/trunk/cmake/3.8.1/i386-apple-darwin16_gcc/bin/cmake
##
## ./build_visit generated host.cmake
## created: Wed Aug  1 08:58:52 CDT 2018
## system: Darwin qwu.MacAir.local 17.7.0 Darwin Kernel Version 17.7.0: Thu Jun 21 22:53:14 PDT 2018; root:xnu-4570.71.2~1/RELEASE_X86_64 x86_64
## by: qwu

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /Users/qwu/Work/projects/visit/3rdparty/installs/trunk)
SET(VISITARCH i386-apple-darwin17_gcc)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_FORTRAN_COMPILER /usr/local/bin/gfortran TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS "-fno-common -fexceptions -O2 -g" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS "-fno-common -fexceptions -O2 -g -std=c++98" TYPE STRING)

##
## Parallel Build Setup.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER /usr/local/bin/mpicc TYPE FILEPATH)

##
## VisIt Thread Option
##
VISIT_OPTION_DEFAULT(VISIT_THREAD OFF TYPE BOOL)

##############################################################
##
## Database reader plugin support libraries
##
## The HDF4, HDF5 and NetCDF libraries must be first so that
## their libdeps are defined for any plugins that need them.
##
## For libraries with LIBDEP settings, order matters.
## Libraries with LIBDEP settings that depend on other
## Library's LIBDEP settings must come after them.
##############################################################
##

##
## Python
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.14/${VISITARCH})

##
## Qt
##
SETUP_APP_VERSION(QT 5.10.1)
VISIT_OPTION_DEFAULT(VISIT_QT_DIR ${VISITHOME}/qt/${QT_VERSION}/${VISITARCH})

##
## QWT
##
SETUP_APP_VERSION(QWT 6.1.2)
VISIT_OPTION_DEFAULT(VISIT_QWT_DIR ${VISITHOME}/qwt/${QWT_VERSION}/${VISITARCH})

##
## ISPC
##
VISIT_OPTION_DEFAULT(VISIT_ISPC_DIR /Users/qwu/Work/downloads/ispc-v1.9.2-osx)

##
## EMBREE
##
VISIT_OPTION_DEFAULT(VISIT_EMBREE_DIR /Users/qwu/Work/downloads/embree-3.2.0.x86_64.macosx)

##
## TBB
##
VISIT_OPTION_DEFAULT(TBB_ROOT /Users/qwu/Work/downloads/tbb2018_20180312oss)
VISIT_OPTION_DEFAULT(VISIT_TBB_DIR /Users/qwu/Work/downloads/tbb2018_20180312oss)

##
## OSPRay
##
VISIT_OPTION_DEFAULT(VISIT_OSPRAY ON TYPE BOOL)
SETUP_APP_VERSION(OSPRAY 1.7.0)
VISIT_OPTION_DEFAULT(VISIT_OSPRAY_DIR /Users/qwu/Work/projects/ospray/install/lib/cmake/ospray-1.7.0/../../../)

##
## VTK
##
SETUP_APP_VERSION(VTK 8.1.0)
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH})

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.14/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP  /usr/lib z TYPE STRING)

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/77c708f9090236b576669b74c53e9f105eedbd7e/${VISITARCH})
##

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.10.2/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} /usr/lib z TYPE STRING)

