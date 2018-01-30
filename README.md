# VisIt-OSPRay

[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit) is an Open Source, interactive, scalable, visualization, animation and analysis tool. From Unix, Windows or Mac workstations, users can interactively visualize and analyze data ranging in scale from small (<101 core) desktop-sized projects to large (>105 core) leadership-class computing facility simulation campaigns. Users can quickly generate visualizations, animate them through time, manipulate them with avariety of operators and mathematical expressions, and save the resulting images and animations for presentations. VisIt contains a rich set of visualization features to enable users to view a wide variety of data including scalar and vector fields defined on two- and three-dimensional (2D and 3D) structured, adaptive and unstructured meshes.

[OSPRay](https://github.com/wilsonCernWq/ospray) is a ray-tracing based rendering engine for high-fidelity visualization, which features interactive CPU rendering capabilities geared towards Scientific Visualization applications. Advanced shading effects such as Ambient Occlusion, shadows, and transparency can be rendered interactively, enabling new insights into data exploration.

[IDX](https://github.com/sci-visus/PIDX) data format IDX is a binary format for storing structured data such as 3D Cartesian grids. A primary use of the format is to store data produced by large-scale scientific simulations. Often times different types of analysis require the same data to be organized in different data structures to be efficient. This quickly becomes a bottleneck when the data size increases to hundreds of gigabytes or terabytes. IDX facilitates progressively streaming and multi-resolution kinds of analysis and visualization, while at the same time performs decently at other kinds of task.

This repository integrates OSPRay as the scalable volume rendering back-end into VisIt. It also uses PIDX as the data reader. The project is particularly sponsored by CCMSC & PSAAP II Project for visualizing exscale simulation data. This is an on-going project, please email me for more details.

## Build

### Building OSPRay

Before the VisIt-OSPRay module is completely pushed into VisIt's official repository, you need to build OSPRay yourself first. There is a [fork](https://github.com/wilsonCernWq/ospray) of the official ospray repository which being actively
maintained for this project. You will also need to download a seperarte ospray module [`module_visit`](https://github.com/wilsonCernWq/module_visit) for extra functionalities being added into ospray.

Then should be able to build ospray following OSPRay's documentation. Please pass `-DOSPRAY_MODULE_VISIT=ON` flag to cmake
in order to have the visit module compiled. There is also a [script](https://github.com/wilsonCernWq/ospray/blob/master/premake.sh) that might be helpful for building ospray. You can read
it or run `./premake.sh --help` to learn its usage.

### Building VisIt

Currently this repository can be used with this visit revision. For details about building VisIt, please go
to the Wiki [page](https://github.com/wilsonCernWq/VisIt-OSPRay/wiki).

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

## How to use this repository

This repository only contains all the file I changed in VisIt source. Therefore the best way to 
apply those changes is to create a patch file. 

### What are the files will be removed
There are some files should be removed from VisIt source due to some renaming. You can find the file lists
in [`removelist.txt`](removelist.txt) file. However it will not cause any problem if those files are not removed.

### How to create a patch file

We can create a patch using `diff` on linux:

```
./makeref.sh /path/to/visit/src/
diff -Naur src-ref src > patch.txt
cd /path/to/visit/src/
patch -p1  < ../working/patch.txt
```

If you want to undo the change, you can revert a patch by this command:

```
cd /path/to/visit/src/
patch -p1 -R  < ../working/patch.txt
```
Remember to revert the previous patch before applying a new one
