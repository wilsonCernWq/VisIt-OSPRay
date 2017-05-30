# VisIt-OSPRay

[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) is an Open Source, 
interactive, scalable, visualization, animation and analysis tool. 
From Unix, Windows or Mac workstations, users can interactively visualize and 
analyze data ranging in scale from small (<101 core) desktop-sized projects to 
large (>105 core) leadership-class computing facility simulation campaigns. 
Users can quickly generate visualizations, animate them through time, manipulate them 
with avariety of operators and mathematical expressions, and save the resulting images 
and animations for presentations. VisIt contains a rich set of visualization features 
to enable users to view a wide variety of data including scalar and vector fields 
defined on two- and three-dimensional (2D and 3D) structured, adaptive and unstructured 
meshes.

[OSPRay](http://www.ospray.org/index.html) is a ray-tracing based rendering engine for 
high-fidelity visualization, which features interactive CPU rendering capabilities geared 
towards Scientific Visualization applications. Advanced shading effects such as Ambient 
Occlusion, shadows, and transparency can be rendered interactively, enabling new insights
into data exploration.

This repository implements a ray-casting volume rendering backend for VisIt using OSPRay. 
The project is particularly sponsored by [CCMSC](http://ccmsc.utah.edu) & 
[PSAAP II](https://nnsa.energy.gov/mediaroom/pressreleases/psaap062713) Project for 
visualizing their exscale simulation data.
This is an on-going project, lease email me for more details if you want to use it.

## CONTENTS

System Information

Building OSPRay

Building VisIt

## System Information

1. Getting information about the platform
  * Which shell: `ps -p $$   or echo $SHELL`
  * How many cores: `cat /proc/cpuinfo` or `lscpu`
  * where is a library: e.g `which mpic++`

2. Platform Information: Location of MPI
* valhalla.sci.utah.edu _8_ _cores_ _bash_
```
export PAR_COMPILER=`which mpic++`
export PAR_INCLUDE="-I/usr/lib/openmpi/include"
export PAR_LIBS="-L/usr/lib/openmpi/lib"
```
* manasa.sci.utah.edu _8_ _cores_ _bash_
```
export PAR_COMPILER=`which mpic++`
export PAR_INCLUDE="-I/usr/lib/openmpi/include"
export PAR_LIBS="-L/usr/lib/openmpi/lib"
```
* maverick.tacc.utexas.edu _20_ _cores_ _bash_:
```
export PAR_COMPILER="/opt/apps/intel14/mvapich2/2.0b/bin/mpic++"
export PAR_INCLUDE="-I/opt/apps/intel14/mvapich2/2.0b/include"
export PAR_LIBS="-L/opt/apps/intel14/mvapich2/2.0b/lib"
```
* ash.chpc.utah.edu _16_ _cores_ _bash_
```
export PAR_COMPILER="/uufs/ash.peaks/sys/pkg/mvapich2/1.9i/bin/mpic++"
export PAR_INCLUDE="-I/uufs/ash.peaks/sys/pkg/mvapich2/1.9i/include"
export PAR_LIBS="-L/uufs/ash.peaks/sys/pkg/mvapich2/1.9i/lib"
```
* kingspeak.chpc.utah.edu
```
export PAR_COMPILER="/uufs/kingspeak.peaks/sys/pkg/mvapich2/1.9/bin/mpic++"
export PAR_INCLUDE="-I/uufs/kingspeak.peaks/sys/pkg/mvapich2/1.9/include"
```
3. Misc
* AVOID SUSE!!!

## Building OSPRay
## Building VisIt
## How to create patch
Add the files going to be removed into `removelist.txt`
Create a patch
```
./makeref.sh /path/to/visit/src/
diff -Naur src-ref src > patch.txt
cd /path/to/visit/src/
patch -p1  < ../working/patch.txt
```
Revert a patch
```
cd /path/to/visit/src/
patch -p1 -R  < ../working/patch.txt
```
Remember to revert the previous patch before applying a new one