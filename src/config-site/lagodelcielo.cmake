#/ssd/users/qwu/VisIt/cmake/3.8.1/linux-x86_64_gcc-7.2/bin/cmake
##
## ./build_visit generated host.cmake
## created: Thu Feb  1 19:12:30 MST 2018
## system: Linux lagodelcielo 4.4.104-39-default #1 SMP Thu Jan 4 08:11:03 UTC 2018 (7db1912) x86_64 x86_64 x86_64 GNU/Linux
## by: qwu

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /ssd/users/qwu/VisIt)
SET(VISITARCH linux-x86_64_gcc-7.2)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_FORTRAN_COMPILER /usr/bin/gfortran TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS " -m64 -fPIC  -O2 -g -fvisibility=hidden" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS " -m64 -fPIC  -O2 -g -std=c++98 -fvisibility=hidden" TYPE STRING)

##
## Parallel Build Setup.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER /opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpicc TYPE FILEPATH)

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
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR /usr)
VISIT_OPTION_DEFAULT(PYTHON_INCLUDE_PATH /usr/include/python2.7 )
VISIT_OPTION_DEFAULT(PYTHON_LIBRARY /usr/lib64/libpython2.7.so)
VISIT_OPTION_DEFAULT(PYTHON_LIBRARY_DIR /usr/lib64)
VISIT_OPTION_DEFAULT(PYTHON_VERSION 2.7)
SET(VISIT_PYTHON_SKIP_INSTALL ON)

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
VISIT_OPTION_DEFAULT(VISIT_BOOST_DIR /usr)

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
SET(ospray_DIR /ssd/users/qwu/ospray/install-visit/lib64/cmake/ospray-1.4.3)
SET(embree_DIR /home/sci/qwu/software/embree-2.17.0.x86_64.linux)
SET(TBB_ROOT /home/sci/qwu/software/tbb2017_20160916oss)
SET(ISPC_EXECUTABLE /home/sci/qwu/software/ispc-v1.9.1-linux)
VISIT_OPTION_DEFAULT(VISIT_OSPRAY ON TYPE BOOL)

##
## PIDX
##
SET(PIDX_DIR /ssd/users/qwu/PIDX/install)
