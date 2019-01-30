#/home/qadwu/Work/projects/visit/3rdparty/installs/trunk/cmake/3.9.3/linux-x86_64_gcc-7.3/bin/cmake
##
## ./build_visit generated host.cmake
## created: Mon Jan 28 18:23:36 PST 2019
## system: Linux plenoptic 4.15.0-20-generic #21-Ubuntu SMP Tue Apr 24 06:16:15 UTC 2018 x86_64 x86_64 x86_64 GNU/Linux
## by: qadwu

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /home/qadwu/Work/projects/visit/3rdparty/installs/trunk)
SET(VISITARCH linux-x86_64_gcc-7.3)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_FORTRAN_COMPILER /usr/bin/gfortran TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS " -m64 -fPIC  -O2 -fvisibility=hidden" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS " -m64 -fPIC  -O2 -std=c++11 -fvisibility=hidden" TYPE STRING)

##
## Parallel Build Setup.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER /usr/bin/mpicc TYPE FILEPATH)

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
VISIT_OPTION_DEFAULT(VISIT_ISPC_DIR /home/qadwu/Work/softwares/ispc/ispc-1.9.2-install)

##
## EMBREE
##
VISIT_OPTION_DEFAULT(VISIT_EMBREE_DIR /home/qadwu/Work/softwares/embree/embree-3.2.0-install)

##
## TBB
##
VISIT_OPTION_DEFAULT(TBB_ROOT /home/qadwu/Work/softwares/tbb/tbb-2018_20180618-install)
VISIT_OPTION_DEFAULT(VISIT_TBB_DIR /home/qadwu/Work/softwares/tbb/tbb-2018_20180618-install)

##
## OSPRay
##
VISIT_OPTION_DEFAULT(VISIT_OSPRAY ON TYPE BOOL)
SETUP_APP_VERSION(OSPRAY 1.8.0)
VISIT_OPTION_DEFAULT(VISIT_OSPRAY_DIR /home/qadwu/Work/projects/ospray/builds/ospray-visit-install/lib/cmake/ospray-1.8.0/../../../)

##
## VTK
##
SETUP_APP_VERSION(VTK 8.1.0)
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH})

##
## SZIP
##
VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR ${VISITHOME}/szip/2.1/${VISITARCH})

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.14/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz /usr/lib/x86_64-linux-gnu z TYPE STRING)

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/77c708f9090236b576669b74c53e9f105eedbd7e/${VISITARCH})
##

##
## NetCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR ${VISITHOME}/netcdf/4.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_NETCDF_LIBDEP HDF5_LIBRARY_DIR hdf5_hl HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.10.2/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} /usr/lib/x86_64-linux-gnu z TYPE STRING)

