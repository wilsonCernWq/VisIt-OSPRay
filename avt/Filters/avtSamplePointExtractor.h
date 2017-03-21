/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                           avtSamplePointExtractor.h                       //
// ************************************************************************* //

#ifndef AVT_SAMPLE_POINT_EXTRACTOR_H
#define AVT_SAMPLE_POINT_EXTRACTOR_H

#include "ospray/ospray.h"
#include "ospray/ospcommon/vec.h"

#include <filters_exports.h>

#include <avtDatasetToSamplePointsFilter.h>
#include <avtVolume.h>
#include <avtViewInfo.h>
#include <avtImgCommunicator.h>
#include <avtOpacityMap.h>
#include <imgMetaData.h>

#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include <vtkCamera.h>
#include <vtkMatrix4x4.h>

class  vtkDataArray;
class  vtkDataSet;
class  vtkHexahedron;
class  vtkQuadraticHexahedron;
class  vtkPixel;
class  vtkPyramid;
class  vtkQuad;
class  vtkTetra;
class  vtkTriangle;
class  vtkVoxel;
class  vtkWedge;
class  vtkPolygon;

class  avtHexahedronExtractor;
class  avtHexahedron20Extractor;
class  avtMassVoxelExtractor;
class  avtPointExtractor;
class  avtPyramidExtractor;
class  avtSamplePointArbitrator;
class  avtTetrahedronExtractor;
class  avtWedgeExtractor;

class  avtRayFunction;


// ****************************************************************************
//  Class: avtSamplePointExtractor
//
//  Purpose:
//      This is a component that will take an avtDataset as an input and find
//      all of the sample points from that dataset.
//
//  Programmer: Hank Childs
//  Creation:   December 5, 2000
//
//  Modifications:
//
//    Hank Childs, Sat Jan 27 15:09:34 PST 2001
//    Added support for sending cells when doing parallel volume rendering.
//
//    Kathleen Bonnell, Sat Apr 21, 13:09:27 PDT 2001 
//    Added recursive Execute method to walk down input data tree. 
//
//    Hank Childs, Tue Nov 13 15:51:15 PST 2001
//    Remove boolean argument to Extract<Cell> calls since it is no longer
//    necessary when all of the variables are being extracted.
//
//    Hank Childs, Sun Dec 14 11:07:56 PST 2003
//    Added mass voxel extractor.
//
//    Hank Childs, Fri Nov 19 13:41:56 PST 2004
//    Added view conversion option.
//
//    Hank Childs, Sat Jan 29 13:32:54 PST 2005
//    Added 2D extractors.
//
//    Hank Childs, Sun Dec  4 19:12:42 PST 2005
//    Added support for kernel-based sampling.
//
//    Hank Childs, Sun Jan  1 10:56:19 PST 2006
//    Added RasterBasedSample and KernelBasedSample.
//
//    Hank Childs, Tue Feb 28 08:25:33 PST 2006
//    Added PreExecute.
//
//    Jeremy Meredith, Thu Feb 15 11:44:28 EST 2007
//    Added support for rectilinear grids with an inherent transform.
//
//    Hank Childs, Fri Jun  1 11:47:56 PDT 2007
//    Add method GetLoadingInfoForArrays.
//
//    Hank Childs, Thu Sep 13 14:02:40 PDT 2007
//    Added support for hex-20s.
//
//    Hank Childs, Tue Jan 15 14:17:15 PST 2008
//    Have this class set up custom sample point arbitrators, since it has
//    the most knowledge.
//
//    Hank Childs, Fri Jan  9 14:09:57 PST 2009
//    Add support for jittering.
//
//    Kevin Griffin, Fri Apr 22 16:31:57 PDT 2016
//    Added support for polygons.
//
//    Qi Wu, to be determined
//    *) unified coding style
//
// ****************************************************************************

class AVTFILTERS_API avtSamplePointExtractor 
    : public avtDatasetToSamplePointsFilter
{
public:
    // ---
    // constructor & destructors
    //
                              avtSamplePointExtractor(int, int, int);
    virtual                  ~avtSamplePointExtractor();
    // ---
    // inherited functions 
    //
    virtual const char       *GetType(void)
    { return "avtSamplePointExtractor"; };
    virtual const char       *GetDescription(void)
    { return "Extracting sample points";};
    // ---
    // functions defined locally for this class
    //
    // coding style:
    //   _xx = local variables
    //    xx = class fields
    void                      RegisterRayFunction(avtRayFunction *_rf) { rayfoo = _rf; };
    void                      RestrictToTile(int, int, int, int);
    void                      StartTiling(void) { shouldDoTiling = true; }; // added by Qi 
    void                      StopTiling(void) { shouldDoTiling = false; };
    void                      SendCellsMode(bool);
    void                      SetRectilinearGridsAreInWorldSpace(bool, 
								 const avtViewInfo &,double);
    void                      Set3DMode(bool _m) { modeIs3D = _m; };
    void                      SetKernelBasedSampling(bool);
    void                      SetJittering(bool);
    void                      SetUpArbitrator(std::string &, bool);
    void                      SetTrilinear(bool _t) { trilinearInterpolation = _t; };
    void                      SetRayCastingSLIVR(bool _s) { rayCastingSLIVR = _s; };
    void                      SetRayCastingSLIVRParallel(bool _p)
    { rayCastingSLIVRParallel = _p; };
    void                      SetLighting(bool _l) { lighting = _l; };
    void                      SetLightPosition(double _lp[4])
    { for (int i=0;i<4;i++) { lightPosition[i] = _lp[i]; } }
    void                      SetLightDirection(double _ld[3])
    { for (int i=0;i<3;i++) { lightDirection[i] = _ld[i]; } }
    void                      SetMatProperties(double _matProp[4])
    { for (int i=0;i<4;i++) { materialProperties[i] = _matProp[i]; } }
    void                      SetTransferFn(avtOpacityMap *_transferFn1D)
    { transferFn1D = _transferFn1D; };
    void                      SetViewDirection(double *_vD) 
    { std::copy(_vD, _vD + 3, viewDirection); }
    void                      SetCameraPosition(double *_cp) 
    { std::copy(_cp, _cp + 3, cameraPosition); }
    void                      SetCameraUpVector(double *_cu) 
    { std::copy(_cu, _cu + 3, cameraUpVector); }
    void                      SetCameraAspect(double _a) { cameraAspect = _a; }
    void                      SetClipPlanes(double _camClip[2])
    { clipPlanes[0] = _camClip[0]; clipPlanes[1] = _camClip[1]; }
    void                      SetPanPercentages(double _pan[2])
    { panPercentage[0] = _pan[0]; panPercentage[1] = _pan[1]; }
    void                      SetImageZoom(double _zoom) { imageZoom = _zoom; }
    void                      SetDepthExtents(double _depthExtents[2])
    { depthExtents[0] = _depthExtents[0]; depthExtents[1] = _depthExtents[1]; }
    void                      SetMVPMatrix(vtkMatrix4x4 *_mvp)
    { modelViewProj->DeepCopy(_mvp); }
    void                      getSpatialExtents(double _spatialExtents[6])
    { for (int i=0; i<6; i++) _spatialExtents[i] = minMaxSpatialBounds[i]; }
    void                      getAvgPatchExtents(double _avgPatchExtents[6])
    { for (int i=0; i<3; i++) _avgPatchExtents[i] = avgPatchExtents[i]; }
    void                      getCellDimension(double _cellDimension[6])
    { for (int i=0; i<3; i++) _cellDimension[i] = cellDimension[i]; }
    void                      getProjectedExents(int _projectedExtents[4])
    { for (int i=0; i<4; i++) _projectedExtents[i]=projectedImageExtents[i]; }
    //
    // Getting image information
    //
    // gets the max number of patches it could have
    int                       getTotalAssignedPatches() { return totalAssignedPatches; } 
    // gets the number of patches
    int                       getImgPatchSize(){ return patchCount;};
    // gets the metadata
    imgMetaData               getImgMetaPatch(int patchId)
    { return imageMetaPatchVector.at(patchId);}
    // gets the image & erase its existence
    void                      getnDelImgData(int patchId, imgData &tempImgData);
    // deletes patches
    void                      delImgPatches();
    // Set background buffer
    void                      setDepthBuffer(float *_zBuffer, int _size)
    { depthBuffer = _zBuffer; }
    void                      setRGBBuffer(unsigned char  *_colorBuffer, 
					   int _width, int _height)
    { rgbColorBuffer = _colorBuffer; };
    void                      setBufferExtents(int _extents[4])
    { for (int i=0;i<4; i++) bufferExtents[i] = _extents[i]; }
    // Qi add
    void             SetOSPCamera(OSPCamera* _cam) { ospCamera = _cam; }
    void             SetOSPTransferFcn(OSPTransferFunction* _t) { ospTransferFcn = _t; }
    void             ActiveOSPData() { isDataDirty = true; }
    void             SetOSPVolumeList(std::vector<ospVolumeMeta>& l) { ospVolumeList = &l; }
    void             SetRendererSampleRate(double r) { rendererSampleRate = r; }
    void             SetOSPRayContext(OSPContext& o) { ospray = &o; }
    
public:
    typedef std::multimap<int, imgData>::iterator iter_t;

public:
    // Output data for RC SLIVR
    std::vector<imgMetaData>    imageMetaPatchVector;
    std::multimap<int, imgData> imgDataHashMap;

protected:
    int                       width,       height,       depth;
    int                       currentNode, totalNodes;
    int                       widthMin,    widthMax;
    int                       heightMin,   heightMax;
    bool                      shouldDoTiling;
    bool                      modeIs3D;
    bool                      kernelBasedSampling;
    double                    pointRadius;
    double                    minMaxSpatialBounds[6];
    double                    avgPatchExtents[3];
    double                    cellDimension[3];
    // background + other plots
    float                     *depthBuffer;             // depth buffer for the background and other plots
    unsigned char             *rgbColorBuffer;          // bounding box + pseudo color + ...
    int                       bufferExtents[4];         // extents of the buffer( minX, maxX, minY, maxY)
    // attributor
    bool                      shouldSetUpArbitrator;
    std::string               arbitratorVarName;
    bool                      arbitratorPrefersMinimum;
    avtSamplePointArbitrator *arbitrator;
    // different extractors
    avtHexahedronExtractor   *hexExtractor;
    avtHexahedron20Extractor *hex20Extractor;
    avtMassVoxelExtractor    *massVoxelExtractor;
    avtPointExtractor        *pointExtractor;
    avtPyramidExtractor      *pyramidExtractor;
    avtTetrahedronExtractor  *tetExtractor;
    avtWedgeExtractor        *wedgeExtractor;
    // miscellaneous
    bool                      sendCells;
    bool                      jitter;
    avtRayFunction           *rayfoo;
    bool                      rectilinearGridsAreInWorldSpace;
    avtViewInfo               viewInfo;
    double                    aspect;
    int                       projectedImageExtents[4];
    int                       patchCount;
    int                       totalAssignedPatches;
    // triliniear / raycastin SLIVR
    bool                      trilinearInterpolation;
    bool                      rayCastingSLIVR;
    bool                      rayCastingSLIVRParallel;
    // Camera stuff
    double                    viewDirection[3];  // this is camera direction also
    double                    cameraPosition[3]; // (Qi) camera location in world coordinate
    double                    cameraUpVector[3]; // (Qi) camera up vector direction
    double                    cameraAspect;
    double                    depthExtents[2];
    double                    clipPlanes[2];
    double                    panPercentage[2];
    double                    imageZoom;
    vtkMatrix4x4              *modelViewProj;
    // lighting & material
    bool                      lighting;
    double                    lightPosition[4];
    double                    lightDirection[3];
    double                    materialProperties[4];
    avtOpacityMap             *transferFn1D;
    virtual void              Execute(void);
    virtual void              PreExecute(void);
    virtual void              PostExecute(void);
    virtual void              ExecuteTree(avtDataTree_p);
    void                      SetUpExtractors(void);
    imgMetaData               initMetaPatch(int id);    // initialize a patch
    // Qi modification
    // ...
    // OSPRay stuffs
    //
    bool ospEmptyVolumeList;
    int  ospVolumeId;
    OSPCamera            *ospCamera;
    OSPTransferFunction  *ospTransferFcn;
    bool isDataDirty = false;
    std::vector<ospVolumeMeta> *ospVolumeList;
    OSPContext                 *ospray;
    double                      rendererSampleRate;
    
protected:
    typedef struct 
    {
      std::vector<int>                  cellDataIndex;
      std::vector<int>                  pointDataIndex;
      std::vector<int>                  cellDataSize;
      std::vector<int>                  pointDataSize;
      std::vector<vtkDataArray *>       cellArrays;
      std::vector<vtkDataArray *>       pointArrays;
      int                               nVars;
    } LoadingInfo;

    inline void               ExtractHex(vtkHexahedron*,vtkDataSet*, int, LoadingInfo &);
    inline void               ExtractHex20(vtkQuadraticHexahedron*,vtkDataSet*, int, LoadingInfo &);
    inline void               ExtractVoxel(vtkVoxel *, vtkDataSet *, int, LoadingInfo &);
    inline void               ExtractTet(vtkTetra *, vtkDataSet *, int, LoadingInfo &);
    inline void               ExtractPyramid(vtkPyramid *, vtkDataSet *, int, LoadingInfo &);
    inline void               ExtractWedge(vtkWedge *, vtkDataSet *, int, LoadingInfo &);
    inline void               ExtractTriangle(vtkTriangle *, vtkDataSet *, int, LoadingInfo &);
    inline void               ExtractQuad(vtkQuad *, vtkDataSet *, int, LoadingInfo &);
    inline void               ExtractPixel(vtkPixel *, vtkDataSet *, int, LoadingInfo &);
    inline void               ExtractPolygon(vtkPolygon *, vtkDataSet *, int, LoadingInfo &);
    void                      KernelBasedSample(vtkDataSet *);
    void                      RasterBasedSample(vtkDataSet *, int num = 0);
    virtual bool              FilterUnderstandsTransformedRectMesh();
    void                      GetLoadingInfoForArrays(vtkDataSet *, LoadingInfo &);
};

#endif


