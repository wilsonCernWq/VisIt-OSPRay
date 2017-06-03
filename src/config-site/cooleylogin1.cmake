#/gpfs/mira-fs1/projects/SoPE_2/utah/users/allen/VisIt/2.12/visit/cmake/3.0.2/linux-rhel_6-x86_64_gcc-4.4/bin/cmake
##
## ./build_visit generated host.cmake
## created: Wed May  3 20:03:08 UTC 2017
## system: Linux cooleylogin1 2.6.32-642.11.1.el6.x86_64 #1 SMP Wed Oct 26 10:25:23 EDT 2016 x86_64 x86_64 x86_64 GNU/Linux
## by: allen

##
## Setup VISITHOME & VISITARCH variables.
##
SET(VISITHOME /gpfs/mira-fs1/projects/SoPE_2/utah/users/allen/VisIt/2.12/visit)
SET(VISITARCH linux-rhel_6-x86_64_gcc-4.4)

## Compiler flags.
##
VISIT_OPTION_DEFAULT(VISIT_C_COMPILER gcc TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_CXX_COMPILER g++ TYPE FILEPATH)
VISIT_OPTION_DEFAULT(VISIT_C_FLAGS " -m64 -fPIC -fvisibility=hidden" TYPE STRING)
VISIT_OPTION_DEFAULT(VISIT_CXX_FLAGS " -m64 -fPIC -fvisibility=hidden" TYPE STRING)

##
## Parallel Build Setup.
##
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)
## (configured w/ mpi compiler wrapper)
VISIT_OPTION_DEFAULT(VISIT_MPI_COMPILER /soft/libraries/mpi/mvapich2/gcc/bin/mpicc TYPE FILEPATH)

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
## Mesa
##
VISIT_OPTION_DEFAULT(VISIT_MESA_DIR ${VISITHOME}/mesa/7.10.2/${VISITARCH})

##
## Python
##
VISIT_OPTION_DEFAULT(VISIT_PYTHON_DIR ${VISITHOME}/python/2.7.11/${VISITARCH})

##
## Qt
##
SETUP_APP_VERSION(QT 4.8.3)
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
## Ice-T
##
VISIT_OPTION_DEFAULT(VISIT_ICET_DIR ${VISITHOME}/icet/1.0.0/${VISITARCH})
##

##
## Uintah
##
SETUP_APP_VERSION(UINTAH 1.6.0)
VISIT_OPTION_DEFAULT(VISIT_UINTAH_DIR ${VISITHOME}/uintah/${UINTAH_VERSION}/${VISITARCH})

##
## OSPRa
## -- VTK 6.1.0 forces to use python 2.7.6, need to manually change it
## -- recommend to remove tbb libraries inside embree and ospray binary folder
##                                                                             
SET(OSPRAY_USE_EXTERNAL_EMBREE ON)
SET(ospray_DIR /home/qiwu/software/ospray/install/lib64/cmake/ospray-1.3.0)
SET(embree_DIR /home/qiwu/software/embree-2.15.0.x86_64.linux/lib/cmake/embree-2.15.0)
SET(TBB_ROOT /home/qiwu/software/tbb2017_20170412oss)
SET(ISPC_EXECUTABLE /home/qiwu/software/ispc-v1.9.1-linux)
VISIT_OPTION_DEFAULT(VISIT_OSPRAY ON TYPE BOOL)

##
## PIDX
##
SET(PIDX_DIR /home/qiwu/software/PIDX/install)
