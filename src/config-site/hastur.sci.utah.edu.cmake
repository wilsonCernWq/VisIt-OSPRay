#/home/sci/qwu/software/Hastur/VisIt/visit/cmake/3.8.1/linux-x86_64_gcc-4.8/bin/cmake
##
## ./build_visit generated host.cmake
## created: Fri Jan 26 14:03:23 MST 2018
## system: Linux hastur.sci.utah.edu 3.10.0-327.36.1.el7.x86_64 #1 SMP Sun Sep 18 13:04:29 UTC 2016 x86_64 x86_64 x86_64 GNU/Linux
## by: qwu

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /home/sci/qwu/software/Hastur/VisIt/visit)
SET(VISITARCH linux-x86_64_gcc-4.8)

## Compiler flags.
##
IF (VISIT_USE_ICC)
  MESSAGE(STATUS "Using Intel ICC compiler")
  VISIT_OPTION_DEFAULT(VISIT_C_COMPILER icc TYPE FILEPATH)
  VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER icc TYPE FILEPATH)
  SET(VISIT_OPENMP_FLAG "-qopenmp")
ELSE (VISIT_USE_ICC)
  MESSAGE(STATUS "Using GCC compiler")
  VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
  VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
  SET(VISIT_OPENMP_FLAG "")
ENDIF (VISIT_USE_ICC)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS " -O3 -m64 -fPIC -fvisibility=hidden ${VISIT_OPENMP_FLAG}" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS " -O3 -m64 -fPIC -fvisibility=hidden ${VISIT_OPENMP_FLAG}" TYPE STRING)

##
## Parallel Build Setup.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER /opt/intel/impi/2017.2.174/bin64/mpicc TYPE FILEPATH)

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
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.11/${VISITARCH})

##
## Qt
##
SETUP_APP_VERSION(QT 4.8.3)
VISIT_OPTION_DEFAULT(VISIT_QT4 ON TYPE BOOL)
VISIT_OPTION_DEFAULT(VISIT_QT_DIR ${VISITHOME}/qt/${QT_VERSION}/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_QT_BIN ${VISIT_QT_DIR}/bin)

##
## QWT
##
SETUP_APP_VERSION(QWT 6.1.2)
VISIT_OPTION_DEFAULT(VISIT_QWT_DIR ${VISITHOME}/qwt/${QWT_VERSION}/${VISITARCH})

##
## VTK
##
SETUP_APP_VERSION(VTK 6.1.0)
VISIT_OPTION_DEFAULT(VISIT_VTK_DIR ${VISITHOME}/vtk/${VTK_VERSION}/${VISITARCH})

##
## BOOST
##
SETUP_APP_VERSION(BOOST 1_60_0)
VISIT_OPTION_DEFAULT(VISIT_BOOST_DIR ${VISITHOME}/boost/${BOOST_VERSION}/${VISITARCH})

##
## SZIP
##
VISIT_OPTION_DEFAULT(VISIT_SZIP_DIR ${VISITHOME}/szip/2.1/${VISITARCH})

##
## HDF5
##
VISIT_OPTION_DEFAULT(VISIT_HDF5_DIR ${VISITHOME}/hdf5/1.8.14/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_HDF5_LIBDEP ${VISITHOME}/szip/2.1/${VISITARCH}/lib sz /usr/lib z TYPE STRING)

##
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/1.0.0/${VISITARCH})
##

##
## NetCDF
##
VISIT_OPTION_DEFAULT(VISIT_NETCDF_DIR ${VISITHOME}/netcdf/4.1.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_NETCDF_LIBDEP HDF5_LIBRARY_DIR hdf5_hl HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} TYPE STRING)

##
## Silo
##
VISIT_OPTION_DEFAULT(VISIT_SILO_DIR ${VISITHOME}/silo/4.10.1/${VISITARCH})
VISIT_OPTION_DEFAULT(VISIT_SILO_LIBDEP HDF5_LIBRARY_DIR hdf5 ${VISIT_HDF5_LIBDEP} /usr/lib z TYPE STRING)

##
## Uintah
##
SETUP_APP_VERSION(UINTAH 2.5.0)
VISIT_OPTION_DEFAULT(VISIT_UINTAH_DIR ${VISITHOME}/uintah/${UINTAH_VERSION}/${VISITARCH})

##
## OSPRay
## -- VTK 6.1.0 forces to use python 2.7.6, need to manually change it
## -- recommend to remove tbb libraries inside embree and ospray binary folder
##                                                                             
SET(OSPRAY_USE_EXTERNAL_EMBREE ON)
SET(ospray_DIR /home/sci/qwu/OSPRay/Hastur/install-visit/lib64/cmake/ospray-1.4.3)
SET(embree_DIR /home/sci/qwu/software/embree-2.17.0.x86_64.linux)
SET(TBB_ROOT /home/sci/qwu/software/tbb2017_20160916oss)
SET(ISPC_EXECUTABLE /home/sci/qwu/software/ispc-v1.9.1-linux)
VISIT_OPTION_DEFAULT(VISIT_OSPRAY ON TYPE BOOL)

##
## PIDX
##
SET(PIDX_DIR /home/sci/qwu/software/Hastur/PIDX/install)
