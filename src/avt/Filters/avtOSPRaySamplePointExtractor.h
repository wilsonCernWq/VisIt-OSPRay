/*****************************************************************************
*
* Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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
//                      avtOSPRaySamplePointExtractor.h                       //
// ************************************************************************* //

#ifndef AVT_OSPRAY_SAMPLE_POINT_EXTRACTOR_H
#define AVT_OSPRAY_SAMPLE_POINT_EXTRACTOR_H

#include <filters_exports.h>

#include <avtSamplePointExtractorBase.h>
#include <avtOSPRayCommon.h> // this ensures VISIT_OSPRAY is defined
class     avtOSPRayVoxelExtractor;

#include <vtkMatrix4x4.h>

#include <vector>
#include <map>

// ****************************************************************************
//  Class: avtOSPRaySamplePointExtractor
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
// ****************************************************************************

class AVTFILTERS_API avtOSPRaySamplePointExtractor 
    : public avtSamplePointExtractorBase
{
  public:
                              avtOSPRaySamplePointExtractor(int, int, int);
    virtual                  ~avtOSPRaySamplePointExtractor();

    virtual const char       *GetType(void)
                                         { return "avtOSPRaySamplePointExtractor"; };
    virtual const char       *GetDescription(void)
                                         { return "Extracting sample points";};


    void                      SetLighting(bool l) {lighting = l; };
    void                      SetLightPosition(double lp[4])
                     { for (int i = 0; i < 4; ++i) lightPosition[i] = lp[i]; };
    void                      SetLightDirection(double ld[3]) 
                    { for (int i = 0; i < 3; ++i) lightDirection[i] = ld[i]; };
    void                      SetMatProperties(double _matProp[4]) 
                  { for (int i=0;i<4;i++) materialProperties[i]=_matProp[i]; };
    void                      SetTransferFn(avtOpacityMap *tfn1D) 
                                                     { transferFn1D = tfn1D; };
    void                      SetViewDirection(double *vD)
                         { for (int i=0; i<3; i++) viewDirection[i] = vD[i]; };
    void                      SetClipPlanes(double cp[2])
                             { clipPlanes[0] = cp[0]; clipPlanes[1] = cp[1]; };
    void                      SetPanPercentages(double p[2])
                         { panPercentage[0] = p[0]; panPercentage[1] = p[1]; };
    void                      SetDepthExtents(double de[2])
                         { depthExtents[0] = de[0]; depthExtents[1] = de[1]; };
    void                      SetMVPMatrix(vtkMatrix4x4 *mvp)
                                             { modelViewProj->DeepCopy(mvp); };

    void                      GetSpatialExtents(double se[6])
                   { for (int i=0; i<6; i++) se[i] = minMaxSpatialBounds[i]; };
    void                      GetAvgPatchExtents(double ae[6])
                       { for (int i=0; i<3; i++) ae[i] = avgPatchExtents[i]; };
    void                      GetCellDimension(double cd[6])
                         { for (int i=0; i<3; i++) cd[i] = cellDimension[i]; };
    //void                      GetProjectedExents(int pe[4])
    //               { for (int i=0; i<4; i++) pe[i]=projectedImageExtents[i]; };

    // Getting Image Information
    // Gets the max number of patches it could have
    //int                       GetTotalAssignedPatches() 
    //                                          { return totalAssignedPatches; };
    // Gets the number of patches
    int                       GetImgPatchSize() { return patchCount; };
    // Gets the metadata
    ospray::ImgMetaData       GetImgMetaPatch(int patchId)
                                  { return imageMetaPatchVector.at(patchId); };
    // Gets the image & erase its existence
    void                      GetAndDelImgData(int patchId,
					       ospray::ImgData &tempImgData);
    // Deletes patches
    void                      DelImgPatches();

    // Set background buffer
    void                      SetImageZoom(double z) { imageZoom = z; }
    void                      SetDepthBuffer(float *zBuffer, int size)
                                                     { depthBuffer= zBuffer; };
    void                      SetRGBBuffer(unsigned char *cb, int w, int h)
                                                      { rgbColorBuffer = cb; };
    void                      SetBufferExtents(int e[4])
                          { for (int i=0; i<4; i++) bufferExtents[i] = e[i]; };

    // Added by Qi (March 2018) for RayCasting:OSPRay  
    void SetOSPRay(OSPVisItContext* o) { ospray = o; }
    void SetRendererSampleRate(double r) { rendererSampleRate = r; }
    void SetFullImageExtents(int extents[4]) 
    {
	fullImageExtents[0] = extents[0];
	fullImageExtents[1] = extents[1];
	fullImageExtents[2] = extents[2];
	fullImageExtents[3] = extents[3];
    }

    // Output data for RayCasting OSPRay
    std::vector<ospray::ImgMetaData>    imageMetaPatchVector;
    std::multimap<int, ospray::ImgData> imgDataHashMap;
    typedef std::multimap<int, ospray::ImgData>::iterator iter_t;

  protected:
    virtual void              DoSampling(vtkDataSet *, int);
    virtual void              SetUpExtractors(void);
    virtual void              SendJittering(void);
    virtual bool              FilterUnderstandsTransformedRectMesh(void);

    double                    minMaxSpatialBounds[6];
    double                    avgPatchExtents[3];
    double                    cellDimension[3];

    // Background + other plots
    // depth buffer for the background and other plots
    float                    *depthBuffer; 
    // bounding box + pseudo color + ...
    unsigned char            *rgbColorBuffer; 
    // extents of the buffer (minX, maxX, minY, maxY)
    int                       bufferExtents[4];

    avtOSPRayVoxelExtractor  *osprayVoxelExtractor;
    int                       patchCount;

    // Camera stuff
    double                    viewDirection[3];
    double                    depthExtents[2];
    double                    clipPlanes[2];
    double                    panPercentage[2];
    double                    imageZoom;
    vtkMatrix4x4             *modelViewProj;

    // lighting & material
    bool                      lighting;
    double                    lightPosition[4];
    double                    lightDirection[3];
    double                    materialProperties[4];

    // OSPRay
    int                       fullImageExtents[4];
    OSPVisItContext          *ospray;
    double                    rendererSampleRate;

    // others
    ospray::ImgMetaData       InitMetaPatch(int id);    // initialize a patch
    void                      RasterBasedSample(vtkDataSet *, int num = 0);

};


#endif
