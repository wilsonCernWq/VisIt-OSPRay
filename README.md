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

## Build

### Building OSPRay
Before the VisIt-OSPRay module is completely pushed into VisIt's official repository, you
need to build OSPRay yourself first. There is a bash [script](https://gist.github.com/wilsonCernWq/48a6aacd8c42ba45aae2b801a881d9e0#file-make_ospray-sh-L1) 
which should be able to build ospray with VisIt module automatically.
If it doesn't work, you can clone the [ospray](https://github.com/wilsonCernWq/ospray) from my github,
and build it as normal ospray. You might need to add an additional cmake flag `-DOSPRAY_MODULE_VISIT=ON`
in order to build the visit module. 

There is a different building script locates in the `ospray` repository called `premake.sh`. You can read
it or run `./premake.sh --help` to learn its usage.

### Building VisIt

```
Working Copy Root Path: /home/sci/qwu/software/Kepler/VisIt/visit-trunk/src
URL: http://visit.ilight.com/svn/visit/trunk/src
Repository Root: http://visit.ilight.com/svn/visit
Repository UUID: 18c085ea-50e0-402c-830e-de6fd14e8384
Revision: 32167
Node Kind: directory
Schedule: normal
Last Changed Author: bonnell
Last Changed Rev: 32166
Last Changed Date: 2018-01-05 18:12:50 -0700 (Fri, 05 Jan 2018)
```

## System Information

Getting information about the platform
  * Which shell: `ps -p $$   or echo $SHELL`
  * How many cores: `cat /proc/cpuinfo` or `lscpu`
  * where is a library: e.g `which mpic++`

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
