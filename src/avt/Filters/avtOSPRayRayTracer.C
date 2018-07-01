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
//                       avtOSPRayRayTracer.C                                //
// ************************************************************************* //

#include <avtOSPRayRayTracer.h>

#include <visit-config.h>

#include <avtDataset.h>
#include <avtImage.h>
#include <avtParallel.h>
#include <avtRayCompositer.h>
#include <avtOSPRaySamplePointExtractor.h>
#include <avtWorldSpaceToImageSpaceTransform.h>

#include <DebugStream.h>
#include <ImproperUseException.h>
#include <TimingsManager.h>
#include <StackTimer.h>

#include <vtkImageData.h>
#include <vtkMatrix4x4.h>

#include <vector>

using     std::vector;

bool OSPRaySortImgMetaDataByDepth(ospray::ImgMetaData const& before, 
                                  ospray::ImgMetaData const& after)
{ return before.avg_z > after.avg_z; }
bool OSPRaySortImgMetaDataByEyeSpaceDepth(ospray::ImgMetaData const& before,
                                          ospray::ImgMetaData const& after)
{ return before.eye_z > after.eye_z; }

// ****************************************************************************
//  Method: avtOSPRayRayTracer constructor
//
//  Programmer: Hank Childs
//  Creation:   November 27, 2000
//
//  Modifications:
//    Brad Whitlock, Wed Dec 5 11:19:29 PDT 2001
//    Added gradient backgrounds.
//
//    Hank Childs, Mon Jan 16 11:11:47 PST 2006
//    Added kernel based sampling.
//
//    Pascal Grosset, Fri Sep 20 2013
//    Added ray casting slivr & trilinear interpolation
//
//    Qi WU, Wed Jun 20 2018
//    Added support for raycasting ospray
//
// ****************************************************************************

avtOSPRayRayTracer::avtOSPRayRayTracer() : avtRayTracerBase()
{
    panPercentage[0] = 0;
    panPercentage[1] = 0;
    lighting = false;

    materialProperties[0] = 0.4;
    materialProperties[1] = 0.75;
    materialProperties[2] = 0.0;
    materialProperties[3] = 15.0;

    ospray = NULL;
}


// ****************************************************************************
//  Method: avtOSPRayRayTracer destructor
//
//  Purpose:
//      Defines the destructor.  Note: this should not be inlined in the header
//      because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   February 5, 2004
//
// ****************************************************************************

avtOSPRayRayTracer::~avtOSPRayRayTracer()
{
    ;
}


// ****************************************************************************
//  Method: avtOSPRayRayTracer::Execute
//
//  Purpose:
//      Executes the ray tracer.  
//      This means:
//      - Put the input mesh through a transform so it is in camera space.
//      - Get the sample points.
//      - Communicate the sample points (parallel only).
//      - Composite the sample points along rays.
//      - Communicate the pixels from each ray (parallel only).
//      - Output the image.
//
//  Programmer: Hank Childs
//  Creation:   November 27, 2000
//
//  Modifications:
//
//    Hank Childs, Fri Jan 26 17:00:07 PST 2001
//    Added support for parallel ray tracing.
//
//    Hank Childs, Sat Jan 27 15:09:34 PST 2001
//    Added support for sending cells when doing parallel ray tracing.
//
//    Hank Childs, Sun Mar  4 21:50:36 PST 2001
//    Added image partition to aid parallel modules in minimizing
//    communication by dynamically assigning the partitions.
//
//    Eric Brugger, Mon Nov  5 13:47:19 PST 2001
//    Modified to always compile the timing code.
//
//    Brad Whitlock, Wed Dec 5 10:57:56 PDT 2001
//    Added code to set the gradient background colors, background mode in
//    the ray compositer.
//
//    Kathleen Bonnell, Tue Apr 30 12:34:54 PDT 2002  
//    Use new version of CopyTo. 
//
//    Hank Childs, Fri Sep 13 12:04:04 PDT 2002
//    Reverse arguments for CopyTo (we weren't sending in the input correctly).
//
//    Hank Childs, Fri Nov 19 13:47:20 PST 2004
//    Added option to have sampling of rectilinear grids done efficiently
//    by sample point extractor by not converting grid into image space.
//
//    Hank Childs, Thu Dec  2 09:26:28 PST 2004
//    No longer tighten clipping planes ['5699].
//
//    Hank Childs, Thu Dec  9 17:15:44 PST 2004
//    Cast rays in tiles for big images ['1948].
//
//    Hank Childs, Sun Mar 13 11:16:20 PST 2005
//    Fix memory leak.
//
//    Hank Childs, Tue Mar 29 16:19:19 PST 2005
//    If the image is large, force divisions of 512x512, even if we should
//    theoretically have enough memory to cover it.
//
//    Hank Childs, Mon Jan 16 11:11:47 PST 2006
//    Add support for kernel based sampling.
//
//    Hank Childs, Thu May 31 22:36:32 PDT 2007
//    Remove code for setting number of variables for processors that don't
//    have data, since that is now handled by the avtSamplePointExtractor's
//    base class.
//
//    Hank Childs, Wed Dec 24 09:46:56 PST 2008
//    Add code to convert the z-buffer of the background to the w-buffer.
//    This is being done so the samples lie more evenly.
//
//    Hank Childs, Wed Dec 24 14:16:04 PST 2008
//    Automatically tighten the clipping planes so we utilize our samples
//    more effectively.
//
//    Hank Childs, Fri Jan  9 14:17:10 PST 2009
//    Have each ray be randomly jittered forward or backwards along its
//    direction.
//
//    Hank Childs, Tue Jan 13 14:26:44 PST 2009
//    Fix oversight where parallel volume rendering was not being jittered.
//
//    Hank Childs, Sun Jan 24 15:35:50 PST 2010
//    Automatically use the kernel based resampling for point data.
//
//    Pascal Grosset & Manasa Prasad, Fri Aug 20 2016
//    Add the ray casting OSPRay code
//
//    Qi WU, Tue Aug 8 2017
//    Fix camera matrices multiplication order for ray casting OSPRay
//    Also fixed panning for ray casting OSPRay
//
//    Qi WU, Wed Jun 20 2018
//    Added support for raycasting ospray
//
// ****************************************************************************

void
avtOSPRayRayTracer::Execute()
{
    //======================================================================//
    // Initialization and Debug
    //======================================================================//
    // check memory in the beginning
    ospout << "[avrRayTracer] entering execute" << std::endl;
    ospray::CheckMemoryHere("[avtOSPRayRayTracer] Execute", "ospout");    

    // initialize current time
    StackTimer t0("Ray Tracing");

    //======================================================================//
    // Start of original pipeline
    //======================================================================//
    bool parallelOn = (imgComm.GetParSize() == 1) ? false : true;
    if (rayfoo == NULL)
    {
        debug1 << "Never set ray function for ray tracer." << endl;
        EXCEPTION0(ImproperUseException);
    }

    //
    // First we need to transform all of domains into camera space.
    //
    ospout << "[avrRayTracer] compute camera" << std::endl;
    double aspect = 1.;
    if (screen[1] > 0)
    {
        aspect = (double)screen[0] / (double)screen[1];
    }

    double scale[3] = {1,1,1};
    vtkMatrix4x4 *transform = vtkMatrix4x4::New();
    avtWorldSpaceToImageSpaceTransform::CalculateTransform(view, transform,
                                                           scale, aspect);
    double newNearPlane, newFarPlane, oldNearPlane, oldFarPlane;
    TightenClippingPlanes(view, transform, newNearPlane, newFarPlane);
    oldNearPlane = view.nearPlane;  oldFarPlane  = view.farPlane;
    view.nearPlane = newNearPlane;  view.farPlane  = newFarPlane;
    transform->Delete();

    avtWorldSpaceToImageSpaceTransform trans(view, aspect);
    trans.SetInput(GetInput());

    //
    // Extract all of the samples from the dataset.
    //
    avtOSPRaySamplePointExtractor extractor(screen[0], screen[1],
                                            samplesPerRay);
    extractor.SetJittering(true);
    extractor.SetTransferFn(transferFn1D);
    extractor.SetInput(trans.GetOutput());
	
    //
    // Before Rendering
    //
    double dbounds[6];  // Extents of the volume in world coordinates
    vtkMatrix4x4  *model_to_screen_transform = vtkMatrix4x4::New();
    vtkMatrix4x4  *screen_to_model_transform = vtkMatrix4x4::New();
    vtkMatrix4x4  *screen_to_camera_transform = vtkMatrix4x4::New();
    vtkImageData  *opaqueImageVTK = NULL;
    unsigned char *opaqueImageData = NULL;
    float         *opaqueImageZB = NULL;
    std::vector<float> opaqueImageDepth(screen[0] * screen[1], oldFarPlane);
    int            tileExtents[4];

    //
    // Camera Settings
    //
    vtkCamera *sceneCam = vtkCamera::New();
    sceneCam->SetPosition(view.camera[0],view.camera[1],view.camera[2]);
    sceneCam->SetFocalPoint(view.focus[0],view.focus[1],view.focus[2]);
    sceneCam->SetViewUp(view.viewUp[0],view.viewUp[1],view.viewUp[2]);
    sceneCam->SetViewAngle(view.viewAngle);
    sceneCam->SetClippingRange(oldNearPlane, oldFarPlane);
    if (view.orthographic) { sceneCam->ParallelProjectionOn(); }
    else { sceneCam->ParallelProjectionOff(); }
    sceneCam->SetParallelScale(view.parallelScale);	
    // Clip planes
    double oldclip[2] = {oldNearPlane, oldFarPlane};
    panPercentage[0] = view.imagePan[0];
    panPercentage[1] = view.imagePan[1];
    // Scaling
    vtkMatrix4x4 *matScale = vtkMatrix4x4::New();
    matScale->Identity(); 
    // Scale + Model + View Matrix
    vtkMatrix4x4 *matViewModelScale = vtkMatrix4x4::New();
    vtkMatrix4x4 *matViewModel = sceneCam->GetModelViewTransformMatrix();
    vtkMatrix4x4::Multiply4x4(matViewModel, matScale, matViewModelScale);
    // Zooming
    vtkMatrix4x4 *matZoomViewModelScale = vtkMatrix4x4::New();
    vtkMatrix4x4 *matZoom = vtkMatrix4x4::New();
    matZoom->Identity(); 
    matZoom->SetElement(0, 0, view.imageZoom); 
    matZoom->SetElement(1, 1, view.imageZoom);
    vtkMatrix4x4::Multiply4x4(matZoom, matViewModelScale, 
                              matZoomViewModelScale);
    // Projection:
    //
    // https://www.vtk.org/doc/release/6.1/html/classvtkCamera.html
    // HASH: #a4d9a509bf60f1555a70ecdee758c2753
    //
    // The Z buffer that is passed from visit is in clip scape with z 
    // limits of -1 and 1. However, using VTK 6.1.0, the z limits are 
    // wired. So, the projection matrix from VTK is hijacked here and
    // adjusted to be within -1 and 1 too
    //
    // Actually the correct way of using VTK GetProjectionTransformMatrix 
    // is to set near and far plane as -1 and 1
    //
    vtkMatrix4x4 *matProj = 
        sceneCam->GetProjectionTransformMatrix(aspect, -1, 1);
    double sceneSize[2];
    if (!view.orthographic) {
        sceneSize[0] = 2.0 * oldNearPlane / matProj->GetElement(0, 0);
        sceneSize[1] = 2.0 * oldNearPlane / matProj->GetElement(1, 1);
    }
    else {
        sceneSize[0] = 2.0 / matProj->GetElement(0, 0);
	    sceneSize[1] = 2.0 / matProj->GetElement(1, 1);
    }
    // Compute model_to_screen_transform matrix
    vtkMatrix4x4::Multiply4x4(matProj,matZoomViewModelScale,
                              model_to_screen_transform);
    vtkMatrix4x4::Invert(model_to_screen_transform,
                         screen_to_model_transform);
    vtkMatrix4x4::Invert(matProj,
                         screen_to_camera_transform);
    // Debug
    ospout << "[avrRayTracer] matZoom " << *matZoom << std::endl;
    ospout << "[avrRayTracer] matViewModel " << *matViewModel << std::endl;
    ospout << "[avrRayTracer] matScale " << *matScale << std::endl;
    ospout << "[avrRayTracer] matProj " << *matProj << std::endl;
    // Cleanup
    matScale->Delete();
    matViewModel->Delete();
    matViewModelScale->Delete();
    matZoom->Delete();
    matZoomViewModelScale->Delete();
    matProj->Delete();
    // Get the full image extents of the volume
    double depthExtents[2];
    GetSpatialExtents(dbounds);
    ospray::ProjectWorldToScreenCube(dbounds, screen[0], screen[1], 
                                     panPercentage, view.imageZoom,
                                     model_to_screen_transform,
                                     tileExtents, depthExtents);
    tileExtents[0] = std::max(tileExtents[0], 0);
    tileExtents[2] = std::max(tileExtents[2], 0);
    tileExtents[1] = std::min(1+tileExtents[1], screen[0]);
    tileExtents[3] = std::min(1+tileExtents[3], screen[1]);
    // Debug
    ospout << "[avrRayTracer] View settings: " << endl
           << "  inheriant view direction: "
           << viewDirection[0] << " "
           << viewDirection[1] << " "
           << viewDirection[2] << std::endl
           << "  camera: "       
           << view.camera[0] << ", " 
           << view.camera[1] << ", " 
           << view.camera[2] << std::endl
           << "  focus: "    
           << view.focus[0] << ", " 
           << view.focus[1] << ", " 
           << view.focus[2] << std::endl
           << "  viewUp: "    
           << view.viewUp[0] << ", " 
           << view.viewUp[1] << ", " 
           << view.viewUp[2] << std::endl
           << "  viewAngle: " << view.viewAngle << std::endl
           << "  eyeAngle:  " << view.eyeAngle  << std::endl
           << "  parallelScale: " << view.parallelScale  << std::endl
           << "  setScale: " << view.setScale << std::endl
           << "  scale:    " 
           << scale[0] << " " 
           << scale[1] << " " 
           << scale[2] << " " 
           << std::endl
           << "  nearPlane: " << view.nearPlane << std::endl
           << "  farPlane:  " << view.farPlane  << std::endl
           << "  imagePan[0]: " << view.imagePan[0] << std::endl 
           << "  imagePan[1]: " << view.imagePan[1] << std::endl
           << "  imageZoom:   " << view.imageZoom   << std::endl
           << "  orthographic: " << view.orthographic << std::endl
           << "  shear[0]: " << view.shear[0] << std::endl
           << "  shear[1]: " << view.shear[1] << std::endl
           << "  shear[2]: " << view.shear[2] << std::endl
           << "  oldNearPlane: " << oldNearPlane << std::endl
           << "  oldFarPlane:  " << oldFarPlane  << std::endl
           << "  aspect: " << aspect << std::endl
           << "[avrRayTracer] sceneSize: " 
           << sceneSize[0] << " " 
           << sceneSize[1] << std::endl
           << "[avrRayTracer] screen: " 
           << screen[0] << " " << screen[1] << std::endl
           << "[avrRayTracer] data bounds: " << std::endl
           << "\t" << dbounds[0] << " " << dbounds[1] << std::endl
           << "\t" << dbounds[2] << " " << dbounds[3] << std::endl
           << "\t" << dbounds[4] << " " << dbounds[5] << std::endl
           << "[avrRayTracer] full image extents: " << std::endl
           << "\t" << tileExtents[0] << " "
           << "\t" << tileExtents[1] << std::endl
           << "\t" << tileExtents[2] << " "
           << "\t" << tileExtents[3] << std::endl;
    ospout << "[avrRayTracer] model_to_screen_transform: " 
           << *model_to_screen_transform << std::endl;
    ospout << "[avrRayTracer] screen_to_model_transform: " 
           << *screen_to_model_transform << std::endl;
    ospout << "[avrRayTracer] screen_to_camera_transform: " 
           << *screen_to_camera_transform << std::endl;

    //===================================================================//
    // ospray stuffs
    //===================================================================//
    ospray::CheckMemoryHere("[avtOSPRayRayTracer] Execute before ospray", 
                            "ospout");
    // initialize ospray
    // -- multi-threading enabled
    ospray->InitOSP();
    // camera
    ospout << "[avrRayTracer] make ospray camera" << std::endl;
    ((ospray::visit::Camera)ospray->camera)
        .Set(view.orthographic, view.camera, view.focus, 
             view.viewUp,view.viewAngle, view.imagePan,
             view.imageZoom, sceneSize, screen, tileExtents);
    // transfer function
    ospout  << "[avrRayTracer] make ospray transfer function" 
            << std::endl;
    ((ospray::visit::TransferFunction)ospray->tfn)
        .Set(transferFn1D->GetTableFloat(), 
             transferFn1D->GetNumberOfTableEntries(),
             transferFn1D->GetMin(),
             transferFn1D->GetMax());
    // renderer
    ospout << "[avrRayTracer] make ospray renderer" << std::endl;
    ospray::visit::Renderer ren(ospray->renderer);
    ren.Init();
    ren.Set(0, 1, false, false, false);
    ren.ResetLights();
    ren.AddLight().Set(false, 1.5, 1.0, viewDirection); // sun light
    ren.AddLight().Set(true,  materialProperties[0], 1.0);  // ambient
    ren.AddLight().Set(false, materialProperties[1], 1.0, viewDirection); // diffuse
    for (int i = 0; i < 8; ++i) { // in VisIt there are only 8 lights
        const LightAttributes& la = lightList.GetLight(i);
        if (la.GetEnabledFlag()) {
            if (la.GetType() == LightAttributes::Ambient) {
                ren.AddLight().Set(true, la.GetBrightness(),
                                   (double)la.GetColor().Red() / 255.0,
                                   (double)la.GetColor().Green() / 255.0,
                                   (double)la.GetColor().Blue() / 255.0);
            } else {
                ren.AddLight().Set(false, la.GetBrightness(),
                                   (double)la.GetColor().Red() / 255.0,
                                   (double)la.GetColor().Green() / 255.0,
                                   (double)la.GetColor().Blue() / 255.0,
                                   la.GetDirection());        

            }          
        }
    }
    
    // others
    ospray->SetDataBounds(dbounds);
    ospray->SetScaling(scale);
    ospray->SetActiveVariable(activeVariable);
    // check memory
    ospray::CheckMemoryHere("[avtOSPRayRayTracer] Execute after ospray",
                            "ospout");    

    // 
    // Continuation of previous pipeline
    //
    extractor.SetJittering(false);
    extractor.SetLighting(lighting);
    //extractor.SetLightDirection(lightDirection);
    extractor.SetMatProperties(materialProperties);
    extractor.SetViewDirection(viewDirection);
    extractor.SetTransferFn(transferFn1D);
    extractor.SetClipPlanes(oldclip);
    extractor.SetPanPercentages(view.imagePan);
    extractor.SetImageZoom(view.imageZoom);
    extractor.SetRendererSampleRate(rendererSampleRate); 
    extractor.SetDepthExtents(depthExtents);
    extractor.SetMVPMatrix(model_to_screen_transform);
    extractor.SetFullImageExtents(tileExtents);
    extractor.SetOSPRay(ospray); // sending ospray

    //
    // Capture background
    //
    opaqueImageVTK  = opaqueImage->GetImage().GetImageVTK();
    opaqueImageData = 
        (unsigned char *)opaqueImageVTK->GetScalarPointer(0, 0, 0);
    opaqueImageZB   = opaqueImage->GetImage().GetZBuffer();
    int bufferScreenExtents[4] = {0,screen[0],0,screen[1]};
    extractor.SetDepthBuffer(opaqueImageZB,   screen[0]*screen[1]);
    extractor.SetRGBBuffer  (opaqueImageData, screen[0],screen[1]);
    extractor.SetBufferExtents(bufferScreenExtents);
    // Set the background to OSPRay
    for (int y = 0; y < screen[1]; ++y) {
        for (int x = 0; x < screen[0]; ++x) {
            int index = x + y * screen[0];
            int    screenCoord[2] = {x, y};
            double screenDepth = opaqueImageZB[index] * 2 - 1;
            double worldCoord[3];
            ospray::ProjectScreenToCamera
                (screenCoord, screenDepth, 
                 screen[0], screen[1],
                 screen_to_camera_transform, 
                 worldCoord);
            opaqueImageDepth[index] = -worldCoord[2];
        }
    }
    ospray->SetBgBuffer(opaqueImageData, opaqueImageDepth.data(), screen);
    
    // TODO We cannot delete camera here, why ?
    //sceneCam->Delete();
    
    //
    // For curvilinear and unstructured meshes, it makes sense to convert the
    // cells to image space.  But for rectilinear meshes, it is not the
    // most efficient strategy.  So set some flags here that allow the
    // extractor to do the extraction in world space.
    //
    trans.SetPassThruRectilinearGrids(true);
    extractor.SetRectilinearGridsAreInWorldSpace(true, view, aspect);

    // Qi debug
    ospray::CheckMemoryHere("[avtOSPRayRayTracer] Execute "
                            "raytracing setup done",
                            "ospout");

    // Execute raytracer
    avtDataObject_p samples = extractor.GetOutput();

    //
    // Ray casting: SLIVR ~ After Rendering
    //
    // Only required to force an update 
    // Need to find a way to get rid of that!!!!
    avtRayCompositer rc(rayfoo);
    rc.SetInput(samples);
    avtImage_p image  = rc.GetTypedOutput();
    
    // Execute rendering
    // This will call the execute function
    {
        StackTimer t1("AllPatchRendering");
        image->Update(GetGeneralContract()); 
    }

    //
    // Image Compositing
    //
    // Initialization
    int timingIdx;
    float *compositedData = NULL;
    int compositedW, compositedH;
    int compositedExtents[4];
    // Debug
    int numPatches = extractor.GetImgPatchSize();
    ospout << "[avtOSPRayRayTracer] Total num of patches " 
           << numPatches << std::endl;
    for (int i=0; i<numPatches; i++) {
        ospray::ImgMetaData currImgMeta = extractor.GetImgMetaPatch(i);
        ospout << "[avtOSPRayRayTracer] Rank " << PAR_Rank() << " "
               << "Idx " << i << " (" << currImgMeta.patchNumber << ") " 
               << " depth " << currImgMeta.eye_z << std::endl
               << "current patch size = " 
               << currImgMeta.dims[0] << ", " 
               << currImgMeta.dims[1] << std::endl
               << "current patch starting" 
               << " X = " << currImgMeta.screen_ll[0] 
               << " Y = " << currImgMeta.screen_ll[1] << std::endl
               << "current patch ending" 
               << " X = " << currImgMeta.screen_ur[0] 
               << " Y = " << currImgMeta.screen_ur[1] << std::endl;
    }
    //-------------------------------------------------------------------//
    // IceT: If each rank has only one patch, we use IceT to composite
    //-------------------------------------------------------------------//
    if (imgComm.IceTValid() && extractor.GetImgPatchSize() == 1) {
        //---------------------------------------------------------------//
        // Setup Local Tile
        ospray::ImgMetaData currMeta = extractor.GetImgMetaPatch(0);
        ospray::ImgData     currData;
        currData.imagePatch = NULL;
        extractor.GetAndDelImgData /* do shallow copy inside */
            (currMeta.patchNumber, currData);
        //---------------------------------------------------------------//
        //---------------------------------------------------------------//
        // First Composition
        if (PAR_Size() > 1) { 
            compositedW = tileExtents[1] - tileExtents[0];
            compositedH = tileExtents[3] - tileExtents[2];
            compositedExtents[0] = tileExtents[0];
            compositedExtents[1] = tileExtents[1];
            compositedExtents[2] = tileExtents[2];
            compositedExtents[3] = tileExtents[3];
            if (PAR_Rank() == 0) {
                compositedData = 
                    new float[4 * compositedW * compositedH]();
            }
            int currExtents[4] = 
                {std::max(currMeta.screen_ll[0]-tileExtents[0], 0), 
                 std::min(currMeta.screen_ur[0]-tileExtents[0], 
                          compositedW), 
                 std::max(currMeta.screen_ll[1]-tileExtents[2], 0),
                 std::min(currMeta.screen_ur[1]-tileExtents[2],
                          compositedH)};
            imgComm.IceTInit(compositedW, compositedH);
            imgComm.IceTSetTile(currData.imagePatch, 
                                currExtents,
                                currMeta.eye_z);
            imgComm.IceTComposite(compositedData);
            if (currData.imagePatch != NULL) {
                delete[] currData.imagePatch;
                currData.imagePatch = NULL;
            }
        } else {
            compositedW = currMeta.dims[0];
            compositedH = currMeta.dims[1];
            compositedExtents[0] = tileExtents[0];
            compositedExtents[1] = tileExtents[0] + compositedW;
            compositedExtents[2] = tileExtents[2];
            compositedExtents[3] = tileExtents[2] + compositedH;
            compositedData = currData.imagePatch;
            currData.imagePatch = NULL;
        }
        //---------------------------------------------------------------//
        //---------------------------------------------------------------//
        // Memory
        ospray::CheckMemoryHere("[avtOSPRayRayTracer] Execute "
                                "IceT Compositing Done", 
                                "ospout");
        //---------------------------------------------------------------//
    }
    //-------------------------------------------------------------------//
    // SERIAL: Image Composition
    //-------------------------------------------------------------------//
    else if (parallelOn == false) {
        //---------------------------------------------------------------//
        // Get the Metadata for All Patches
        ospray::CheckSectionStart("avtOSPRayRayTracer", "Execute", timingIdx,
                                  "Serial-Composite: Get the Metadata for "
                                  "All Patches");
        // contains the metadata to composite the image
        std::vector<ospray::ImgMetaData> allPatchMeta;
        std::vector<ospray::ImgData>     allPatchData;
        // get the number of patches
        int numPatches = extractor.GetImgPatchSize();
        for (int i=0; i<numPatches; i++)
        {
            allPatchMeta.push_back(extractor.GetImgMetaPatch(i));
        }
        ospray::CheckSectionStop("avtOSPRayRayTracer", "Execute", timingIdx,
                                 "Serial-Composite: Get the Metadata for "
                                 "All Patches");
        //---------------------------------------------------------------//
        //---------------------------------------------------------------//
        // Sort with the Largest z First
        ospray::CheckSectionStart("avtOSPRayRayTracer", "Execute", timingIdx,
                                  "Serial-Composite: Sort with the Largest "
                                  "z First");
        std::sort(allPatchMeta.begin(), allPatchMeta.end(), 
                  &OSPRaySortImgMetaDataByEyeSpaceDepth);
        ospray::CheckSectionStop("avtOSPRayRayTracer", "Execute", timingIdx,
                                 "Serial-Composite: Sort with the Largest "
                                 "z First");
        //---------------------------------------------------------------//
        //---------------------------------------------------------------//
        // Blend Images
        ospray::CheckSectionStart("avtOSPRayRayTracer", "Execute", timingIdx,
                                  "Serial-Composite: Blend Images");
        compositedW = tileExtents[1] - tileExtents[0];
        compositedH = tileExtents[3] - tileExtents[2];
        compositedExtents[0] = tileExtents[0];
        compositedExtents[1] = tileExtents[0] + compositedW;
        compositedExtents[2] = tileExtents[2];
        compositedExtents[3] = tileExtents[2] + compositedH;	    
        if (PAR_Rank() == 0) {
            compositedData = new float[compositedW * compositedH * 4]();
        }
        for (int i=0; i<numPatches; i++)
        {
            ospray::ImgMetaData currImgMeta = allPatchMeta[i];
            ospray::ImgData     currImgData;
            currImgData.imagePatch = NULL;
            extractor.GetAndDelImgData /* do shallow copy inside */
                (currImgMeta.patchNumber, currImgData);
            const float* currData = currImgData.imagePatch;
            const int currExtents[4] = 
                {currImgMeta.screen_ll[0], currImgMeta.screen_ur[0], 
                 currImgMeta.screen_ll[1], currImgMeta.screen_ur[1]};
            // ospray::WriteArrayToPPM
            // 	("/home/qwu/work/visit/build/img-"+std::to_string(i)+".ppm",
            // 	 currData,
            // 	 currImgMeta.screen_ur[0] - currImgMeta.screen_ll[0],
            // 	 currImgMeta.screen_ur[1] - currImgMeta.screen_ll[1]);
            avtOSPRayImageCompositor::BlendBackToFront(currData,
                                                       currExtents,
                                                       compositedData, 
                                                       compositedExtents);
            // Clean up data
            if (currImgData.imagePatch != NULL) {
                delete[] currImgData.imagePatch;
            }
            currImgData.imagePatch = NULL;
        }
        allPatchMeta.clear();
        allPatchData.clear();
        ospray::CheckSectionStop("avtOSPRayRayTracer", "Execute", timingIdx,
                                 "Serial-Composite: Blend Images");
        //---------------------------------------------------------------//
        //---------------------------------------------------------------//
        // Memory
        ospray::CheckMemoryHere("[avtOSPRayRayTracer] Execute "
                                "Sequential Compositing Done", 
                                "ospout");
        //---------------------------------------------------------------//
    } 
    //
    // PARALLEL: Image Composition
    //
    else { 
        //---------------------------------------------------------------//
        // Parallel Direct Send
        ospray::CheckSectionStart("avtOSPRayRayTracer", "Execute", timingIdx,
                                  "Parallel-Composite: "
                                  "Parallel Direct Send");
        int tags[2] = {1081, 1681};
        int tagGather = 2681;
        int *regions = NULL;
        imgComm.RegionAllocation(regions);
        int myRegionHeight =
            imgComm.ParallelDirectSendManyPatches
            (extractor.imgDataHashMap, extractor.imageMetaPatchVector,
             numPatches, regions, imgComm.GetParSize(), tags, 
             tileExtents);
        imgComm.gatherImages(regions, imgComm.GetParSize(), 
                             imgComm.intermediateImage, 
                             imgComm.intermediateImageExtents, 
                             imgComm.intermediateImageExtents, 
                             tagGather, tileExtents, myRegionHeight);

        ospray::CheckSectionStop("avtOSPRayRayTracer", "Execute", timingIdx,
                                 "Parallel-Composite: "
                                 "Parallel Direct Send");
        //---------------------------------------------------------------//
        //---------------------------------------------------------------//
        // Some Cleanup
        ospray::CheckSectionStart("avtOSPRayRayTracer", "Execute", timingIdx,
                                  "Parallel-Composite: Some Cleanup");
        if (regions != NULL)
            delete [] regions;
        regions = NULL;
        if (imgComm.intermediateImage != NULL)
            delete [] imgComm.intermediateImage;
        imgComm.intermediateImage = NULL;		
        imgComm.Barrier();
        ospray::CheckSectionStop("avtOSPRayRayTracer", "Execute", timingIdx,
                                 "Parallel-Composite: Some Cleanup");
        //---------------------------------------------------------------//
        //---------------------------------------------------------------//
        // Setup for Final Composition
        compositedW = 
            imgComm.finalImageExtents[1] -
            imgComm.finalImageExtents[0];
        compositedH = 
            imgComm.finalImageExtents[3] -
            imgComm.finalImageExtents[2];
        compositedExtents[0] = imgComm.finalImageExtents[0];
        compositedExtents[1] = imgComm.finalImageExtents[1];
        compositedExtents[2] = imgComm.finalImageExtents[2];
        compositedExtents[3] = imgComm.finalImageExtents[3];
        if (PAR_Rank() == 0) {
            compositedData = imgComm.GetFinalImageBuffer();
        }
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        // Memory
        ospray::CheckMemoryHere("[avtOSPRayRayTracer] Execute "
                                "Parallel Compositing Done", 
                                "ospout");
        //--------------------------------------------------------------//
    }	

    ///////////////////////////////////////////////////////////////////
    //
    // Final Composition for Displaying
    //
    ///////////////////////////////////////////////////////////////////
    if (PAR_Rank() == 0) {
        avtImage_p finalImage = new avtImage(this);
        vtkImageData *finalVTKImage = 
            avtImageRepresentation::NewImage(screen[0], screen[1]);
        finalImage->GetImage() = finalVTKImage;
        unsigned char *finalImageBuffer = 
            finalImage->GetImage().GetRGBBuffer();
        ospray::CompositeBackground(screen,
                                    compositedExtents,
                                    compositedW,
                                    compositedH,
                                    compositedData,
                                    opaqueImageData,
                                    opaqueImageZB,
                                    finalImageBuffer);
        finalVTKImage->Delete();
        SetOutput(finalImage);
    }
    if (compositedData != NULL) { 
        delete [] compositedData;
    }
    compositedData = NULL; 
    ospout << "[avtOSPRayRayTracer] Raycasting OSPRay is Done !" << std::endl;

    //
    // Clean up
    //
    screen_to_model_transform->Delete();
    model_to_screen_transform->Delete();
    screen_to_camera_transform->Delete();
    ospray->Finalize();
}
