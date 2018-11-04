# VisIt-OSPRay

This repository integrates [OSPRay](https://github.com/wilsonCernWq/ospray) as a CPU volume rendering backend in [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit). The details of this implementation has been described in our EGPGV paper: [VisIt-OSPRay: Toward an Exascale Volume Visualization System](http://www.sci.utah.edu/~wald/Publications/2018/visit-ospray-pgv18.pdf). In particular, this project was sponsored by Intel and the CCMSC/PSAAP II Project for visualizing exscale simulation data.

## Build

### Building OSPRay

You can use the `build_visit` script to build OSPRay and its dependencies. In order to enable ospray, you need to add the `--ospray` option to your `build_visit`:

`./build_visit --ospray ... `

Because OSPRay depends on TBB, ISPC and Embree, you might want to also enable `--ispc`, `--embree` and `--tbb`. If you have all these libraries downloaded alread, you can use alternative paths instead: 
```
--alt-tbb-dir <path-to-your-tbb>
--alt-ispc-dir <path-to-your-ispc>
--alt-embree-dir <path-to-your-embree>
```

### Building VisIt

Currently this repository is compatable with visit `trunk` revision `33734`. For details about building VisIt, please go
to the Wiki [page](https://github.com/wilsonCernWq/VisIt-OSPRay/wiki).

```
Path: .
URL: http://visit.ilight.com/svn/visit/trunk/src
Repository Root: http://visit.ilight.com/svn/visit
Repository UUID: 18c085ea-50e0-402c-830e-de6fd14e8384
Revision: 33734
Node Kind: directory
Schedule: normal
Last Changed Author: brugger
Last Changed Rev: 33734
Last Changed Date: 2018-10-30 03:06:51 -0700 (Tue, 30 Oct 2018)
```

## How to use this repository

This repository only contains all the file I changed in VisIt source. Therefore the best way to 
apply those changes is to create a patch file. 

### How to create a patch file

We can create a patch using `diff` on linux:

```
./exclude/tools/makeref.sh /path/to/visit/src/
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

### How to directly deploy or remove files

I have also provided a script for directly copying and removing files into the visit source tree. If you want to use this option, you can run the following script:

Copying all the `cmake` and `sh` files
```
./deploy.sh <visit-source> cmake
```


Copying all the source and `sh` files
```
./deploy.sh <visit-source> source
```

Remove unnecessary files
```
./deploy.sh <visit-source> clean
```

If there are some files should be removed from VisIt source due to some renaming. You can find the file lists
in [`exclude/tools/removelist.txt`](exclude/tools/removelist.txt) file. However it will not cause any problem if those files are not removed.

### Running the Development Version
Currently in this branch I have also implemented a version to test the new OSPRay distributed framebuffer. In order to test this functionality, you need to set the environmental variable `VTKOSPRAY_ARGS` to be `"--osp:mpi-distribute"`.
