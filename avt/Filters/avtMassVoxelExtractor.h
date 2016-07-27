/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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
//                            avtMassVoxelExtractor.h                        //
// ************************************************************************* //

#ifndef AVT_MASS_VOXEL_EXTRACTOR_H
#define AVT_MASS_VOXEL_EXTRACTOR_H

#include <filters_exports.h>

#include <avtExtractor.h>
#include <avtViewInfo.h>

#include <avtOpacityMap.h>

class     vtkRectilinearGrid;
class     vtkMatrix4x4;

#include <vtkMatrix3x3.h>
#include <vtkMatrix4x4.h>
#include <vtkCamera.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm> 

// ****************************************************************************
//  Class: avtMassVoxelExtractor
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

class AVTFILTERS_API avtMassVoxelExtractor : public avtExtractor
{
  public:
                     avtMassVoxelExtractor(int, int, int, avtVolume *,
                                            avtCellList *);
    virtual         ~avtMassVoxelExtractor();

    void             Extract(vtkRectilinearGrid *,
                             std::vector<std::string> &varnames,
                             std::vector<int> &varsize);

    void             SetGridsAreInWorldSpace(bool, const avtViewInfo &,double,
                                             const double *);
    void             SetVariableInformation(std::vector<std::string> &names,
                                            std::vector<int> varsize);

    void             SetRayCastingSLIVR(bool s) {rayCastingSLIVR = s; };
    void             SetTrilinear(bool t) {trilinearInterpolation = t;   };

    void             SetLighting(bool l) {lighting = l; };
    void             SetLightDirection(double _lightDir[3]) { for (int i=0;i<3;i++) lightDirection[i]=_lightDir[i]; }
    void             SetLightPosition(double _lightPos[4]) { for (int i=0;i<4;i++) lightPosition[i]=_lightPos[i]; }
    void             SetMatProperties(double _matProp[4]) { for (int i=0;i<4;i++) materialProperties[i]=_matProp[i]; }

    void             SetTransferFn(avtOpacityMap *_transferFn1D) { transferFn1D = _transferFn1D; };

    void             SetViewDirection(double *vD){ for (int i=0; i<3; i++) viewDirection[i]=view_direction[i] = vD[i]; }
    void             SetClipPlanes(double _camClip[2]){ clipPlanes[0]=_camClip[0]; clipPlanes[1]=_camClip[1]; }
    void             SetMVPMatrix(vtkMatrix4x4 *_mvp){ modelViewProj->DeepCopy(_mvp); vtkMatrix4x4::Invert(modelViewProj, invModelViewProj); }



    // Getting the image
    void             getImageDimensions(int &inUse, int dims[2], int screen_ll[2], int screen_ur[2], float &eyeDepth, float &clipDepth);
    void             getComputedImage(float *image);
    void             getImageDepth(float *imageDepth);
    void             setProcIdPatchID(int _proc, int _patch){ proc = _proc; patch = _patch; }

    // TODO: Make that just a pointer instead of copy!!!
    void             setDepthBuffer(float *_zBuffer, int size){ depthBuffer=new float[size]; for (int i=0; i<size; i++) depthBuffer[i]=_zBuffer[i]; }
    void             setRGBBuffer(unsigned char  *_colorBuffer, int width, int height){ rgbColorBuffer=new unsigned char[width*height*3]; for (int i=0; i<width*height*3; i++) rgbColorBuffer[i]=_colorBuffer[i]; };
    void             setBufferExtents(int _extents[4]){ for (int i=0;i<4; i++) bufferExtents[i]=_extents[i]; }

 

  protected:
    bool             gridsAreInWorldSpace;
    bool             pretendGridsAreInWorldSpace;
    avtViewInfo      view;
    double           aspect;
    double           cur_clip_range[2];
    vtkMatrix4x4    *view_to_world_transform;
    vtkMatrix4x4    *world_to_view_transform;


    bool            trilinearInterpolation;
    bool            rayCastingSLIVR;

    vtkMatrix4x4    *modelViewProj;
    vtkMatrix4x4    *invModelViewProj;
    double           clipPlanes[2];
    double           viewDirection[3];
    double           view_direction[3];


    double           *X;
    double           *Y;
    double           *Z;
    unsigned char   *ghosts;
    int              dims[3];
    int              ncell_arrays;
    void            *cell_arrays[AVT_VARIABLE_LIMIT];
    int              cell_size[AVT_VARIABLE_LIMIT];
    int              cell_index[AVT_VARIABLE_LIMIT];
    int              cell_vartypes[AVT_VARIABLE_LIMIT];
    int              npt_arrays;
    void            *pt_arrays[AVT_VARIABLE_LIMIT];
    int              pt_size[AVT_VARIABLE_LIMIT];
    int              pt_index[AVT_VARIABLE_LIMIT];
    int              pt_vartypes[AVT_VARIABLE_LIMIT];


    double          *prop_buffer;
    int             *ind_buffer;
    bool            *valid_sample;

    // We repeatedly divide by the term (X[i+1]-X[i]).  In the interest of
    // performance, cache the term 1./(X[i+1]-X[i]) and use that for faster
    // multiplication.  This sped up total performance by about 5%.
    double           *divisors_X;
    double           *divisors_Y;
    double           *divisors_Z;

    bool             lighting;
    double           lightPosition[4];
    float            lightDirection[3];
    double           materialProperties[4];
    avtOpacityMap    *transferFn1D;
    float            gradient[3];

    
    bool             activeNow;
    int              tempCount;

    int debugOn;
    int countt;

    // Background + other plots
    float           *depthBuffer;           // depth buffer for the background and other plots
    unsigned char   *rgbColorBuffer;        // bounding box + pseudo color + ...      
    int              bufferExtents[4];      // extents of the buffer( minX, maxX, minY, maxY)    

    // Rendering
    int              renderingAreaExtents[4];
    float            renderingDepthsExtents[2];


    // Patch details for one image
    int              patchDrawn;            // whether the patch is drawn or not

    int              imgWidth, imgHeight;   
    int              imgDims[2];            // size of the patch

    int              imgLowerLeft[2];       // coordinates in the whole image
    int              imgUpperRight[2];      //


    float            eyeSpaceDepth;         // for blending patches
    float            clipSpaceDepth;        // clip space depth for blending with other visit stuff

    float            *imgArray;             // the image data
    float            *imgDepths;            // z depth

    int              proc;                  // id of the processor
    int              patch;                 // id of the patch

    
    int              fullImgWidth, fullImgHeight;
    int              xMin, xMax, yMin, yMax;

    void             ExtractImageSpaceGrid(vtkRectilinearGrid *,
                             std::vector<std::string> &varnames,
                             std::vector<int> &varsize);
    void             ExtractWorldSpaceGrid(vtkRectilinearGrid *,
                             std::vector<std::string> &varnames,
                             std::vector<int> &varsize);
    

    void             RegisterGrid(vtkRectilinearGrid*,
                                  std::vector<std::string>&,std::vector<int>&);
    void             SampleAlongSegment(const double *, const double*, int, int);
    void             SampleVariable(int, int, int, int);
    bool             FrustumIntersectsGrid(int, int, int, int) const;
    bool             FrustumIntersectsGridSLIVR(int, int, int, int) const;
    void             GetSegment(int, int, double *, double *) const;
    static void      FindPlaneNormal(const double *, const double *, 
                                     const double *, double *);
    bool             GridOnPlusSideOfPlane(const double *, const double *) const;
    bool             FindSegmentIntersections(const double *, const double *, 
                                              int &, int &);

    void             computePixelColor(double source_rgb[4], double dest_rgb[4]);
    double           trilinearInterpolate(double vals[8], float distRight, float distTop, float distBack);
    void             computeIndices(int dims[3], int indices[6], int returnIndices[8]);
    void             computeIndicesVert(int dims[3], int indices[6], int returnIndices[8]);
    void             getIndexandDistFromCenter(float dist, int index,    int &index_before, int &index_after,    float &dist_before, float &dist_after);

    void            normalize(float vec[3]);
    float           dot(float vecA[3], float vecB[3]);
    void            unProject(int _x, int _y, float _z, double _worldCoordinates[3], int _width, int _height);
    double          project(double _worldCoordinates[3], int pos2D[2], int _width, int _height);

    void            GetSegmentRCSLIVR(int x, int y, float depthsExtents[2], double *_origin, double *_terminus);
    void            SampleVariableRCSLIVR(int first, int last, int intersect, int x, int y);
    void            simpleExtractWorldSpaceGrid(vtkRectilinearGrid *,  // added for raycasting slivr
                             std::vector<std::string> &varnames,
                             std::vector<int> &varsize);
};

#endif
