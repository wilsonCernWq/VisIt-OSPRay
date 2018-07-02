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
//                            avtOSPRayVoxelExtractor.h                      //
// ************************************************************************* //

#ifndef AVT_OSPRAY_VOXEL_EXTRACTOR_H
#define AVT_OSPRAY_VOXEL_EXTRACTOR_H

#include <filters_exports.h>

#include <avtVoxelExtractor.h>
#include <avtOSPRayCommon.h>
#include <avtOpacityMap.h>

#include <DebugStream.h>

#include <vtkMatrix3x3.h>
#include <vtkMatrix4x4.h>
#include <vtkCamera.h>

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

class     vtkRectilinearGrid;
class     vtkMatrix4x4;

// ****************************************************************************
//  Class: avtOSPRayVoxelExtractor
//
//  Purpose:
//      Extracts sample points from a collection of voxels.  It assumes that 
//      the voxels it has been given are in camera space and does not try to
//      populate points that are not in the cube [-1, 1], [-1, 1], [-1, 1].
//
//  Programmer: Hank Childs
//  Creation:   December 14, 2003
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files
//    because certain compilers have problems with them.
//
//    Hank Childs, Fri Nov 19 14:50:58 PST 2004
//    Added support for accepting grids that need to do a world space to
//    image space conversion as well.  Also changed API to AVTFILTERS_API.
//
//    Jeremy Meredith, Thu Feb 15 11:44:28 EST 2007
//    Added support for rectilinear grids with an inherent transform.
//
//    Hank Childs, Fri Jun  1 15:28:14 PDT 2007
//    Added support for non-scalars.
//
//    Hank Childs, Wed Aug 27 11:24:53 PDT 2008
//    Add support for non-floats.
//
//    Hank Childs, Wed Dec 24 11:24:47 PST 2008
//    Remove data member ProportionSpaceToZBufferSpace, as we now do our
//    sampling in even intervals (wbuffer).
//
//    Kathleen Biagas, Fri Jul 13 09:44:45 PDT 2012
//    Use double internally instead of float.
//
// ****************************************************************************
using namespace ospray;
class AVTFILTERS_API avtOSPRayVoxelExtractor : public avtVoxelExtractor
{
  public:
                     avtOSPRayVoxelExtractor(int, int, int, avtVolume *,
                                            avtCellList *);
    virtual         ~avtOSPRayVoxelExtractor();

    void             Extract(vtkRectilinearGrid *,
                             std::vector<std::string> &varnames,
                             std::vector<int> &varsize);

    // void             SetVariableInformation(std::vector<std::string> &names,
    //                                         std::vector<int> varsize);

    // RC OSPRay Specific
    //void             SetRayCastingSLIVR(bool s) { rayCastingSLIVR = s; };
    void             SetLighting(bool l) { lighting = l; };
    /* void             SetLightDirection(double ld[3]) */
    /*                 { for (int i=0;i <3; i++) { lightDirection[i] = ld[i]; } }; */
    /* void             SetLightPosition(double lp[4])  */
    /*                  { for (int i=0; i<4; i++) { lightPosition[i] = lp[i]; } }; */
    void             SetMatProperties(double matProp[4]) 
           { for (int i=0; i<4; i++) { materialProperties[i] = matProp[i]; } };
    void             SetScalarRange(double range[2])
                     { scalarRange[0] = range[0]; scalarRange[1] = range[1]; };
    void             SetTFVisibleRange(double tfRange[2])
           { tFVisibleRange[0] = tfRange[0]; tFVisibleRange[1] = tfRange[1]; };
    //void             SetTransferFn(avtOpacityMap *tf1D) 
    //                                                  { transferFn1D = tf1D; };
    void             SetViewDirection(double *vD)
                     { for (int i=0; i<3; i++) { viewDirection[i] = vD[i]; } };
    void             SetCameraPosition(double *cp) 
                                    { std::copy(cp, cp + 3, cameraPosition); };
    void             SetCameraUpVector(double *cu)
                                    { std::copy(cu, cu + 3, cameraUpVector); };
    void             SetCameraAspect(double a) { cameraAspect = a; };
    void             SetClipPlanes(double cc[2])
                             { clipPlanes[0] = cc[0]; clipPlanes[1] = cc[1]; };
    void             SetPanPercentages(double p[2])
                         { panPercentage[0] = p[0]; panPercentage[1] = p[1]; };
    void             SetImageZoom(double z) { imageZoom = z; };
    void             SetDepthExtents(double d[2])
       { fullVolumeDepthExtents[0] = d[0]; fullVolumeDepthExtents[1] = d[1]; };
    void             SetMVPMatrix(vtkMatrix4x4 *mvp)
    {
	model_to_screen_transform->DeepCopy(mvp); 
	vtkMatrix4x4::Invert(model_to_screen_transform, 
			     screen_to_model_transform); 
    }

    // Getting the image
    void             GetImageDimensions
    (int &, int dims[2], int screen_ll[2], int screen_ur[2], float &, float &);
    void             GetComputedImage(float *image);
    void             SetProcIdPatchID(int c, int p){ proc = c; patch = p; };

    // Set the background information
    //void             SetDepthBuffer(float *z, int size){ depthBuffer = z; };
    //void             SetRGBBuffer(unsigned char *cb, int width, int height)
    //                                                  { rgbColorBuffer = cb; };
    //void             SetBufferExtents(int e[4])
    //                       { for (int i=0;i<4; i++) bufferExtents[i] = e[i]; };
    void             SetRendererSampleRate(double r) 
                                                   { rendererSampleRate = r; };
    void             SetOSPRay(OSPVisItContext* o) { ospray = o; };
    void             SetRenderingExtents(int extents[4]) 
    {
	renderingExtents[0] = extents[0];
	renderingExtents[1] = extents[1];
	renderingExtents[2] = extents[2];	
	renderingExtents[3] = extents[3];
    }

  protected:
    bool            rayCastingOSPRay;

    vtkMatrix4x4    *model_to_screen_transform;
    vtkMatrix4x4    *screen_to_model_transform;

    double           clipPlanes[2];
    double           panPercentage[2];
    double           imageZoom;
    double           fullVolumeDepthExtents[2];
    double           viewDirection[3];
    double           cameraPosition[3]; // (Qi) camera location in world space
    double           cameraUpVector[3]; // (Qi) camera up vector direction
    double           cameraAspect;
    int              renderingExtents[4];

    // Color computation
    bool             lighting;
    double           materialProperties[4];
    float            gradient[3];
    double           scalarRange[2];
    double           tFVisibleRange[2];

    // Rendering
    double           renderingDepthsExtents[2];


    // Patch details for one image
    int              patchDrawn;       // whether the patch is drawn or not
    int              imgWidth;
    int              imgHeight;
    int              imgDims[2];       // size of the patch
    int              imgLowerLeft[2];  // coordinates in the whole image
    int              imgUpperRight[2]; // coordinates in the whole image
    float            eyeSpaceDepth;    // for blending patches
    float            clipSpaceDepth;   // clip space depth for blending with bg

    float            *imgArray;        // the final framebuffer
    int              proc;             // id of the processor
    int              patch;            // id of the patch

    int              fullImgWidth, fullImgHeight;
    int              xMin, xMax, yMin, yMax;

    // OSPRay stuffs
    OSPVisItContext *ospray;
    double           rendererSampleRate;

    // Added for RayCasting OSPRay
    void             ExtractWorldSpaceGridOSPRay(vtkRectilinearGrid *,  
                             std::vector<std::string> &varnames,
                             std::vector<int> &varsize);
};

#endif
