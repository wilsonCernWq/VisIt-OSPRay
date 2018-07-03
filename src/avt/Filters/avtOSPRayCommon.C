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

#include <avtOSPRayCommon.h>

#include <avtParallel.h>
#include <avtMemory.h>

#include <DebugStream.h>
#include <StackTimer.h>
#include <TimingsManager.h>
#include <ImproperUseException.h>

#include <vtkCamera.h>
#include <vtkMatrix4x4.h>

#include <ospray/ospray.h>
#include <ospray/visit/VisItWrapper.h>
#include <ospray/visit/VisItModuleCommon.h>
#include <ospray/visit/VisItImageComposite.h>

#ifdef __unix__
# include <unistd.h>
#endif

inline bool CheckThreadedBlend_MetaData() {
    bool use = true;
    const char* env_use = std::getenv("OSPRAY_SERIAL_BLEND");
    if (env_use) { use = atoi(env_use) <= 0; }
    return use;
}
static bool UseThreadedBlend_MetaData = CheckThreadedBlend_MetaData();

std::ostream *ospray::osp_out = (ospray::visit::CheckVerbose()) ?
    &std::cout : &DebugStream::Stream5();
std::ostream *ospray::osp_err = (ospray::visit::CheckVerbose()) ?
    &std::cerr : &DebugStream::Stream1();

// ****************************************************************************
//
// OSPRay
//
// ****************************************************************************

void OSPContext_ErrorFunc(OSPError, const char* msg)
{ 
    osperr << "#osp: (rank " << PAR_Rank() << ")" << msg; 
}
void OSPContext_StatusFunc(const char* msg)
{ 
    osperr << "#osp: (rank " << PAR_Rank() << ")" << msg; 
}
static bool ospray_initialized = false;
void ospray::Finalize()
{
}
void ospray::InitOSP(int numThreads) 
{   
    if (!ospray_initialized) 
    {
#ifdef __unix__
        // check hostname
        char hname[200];
        gethostname(hname, 200);
        ospout << "[ospray] on host >> " << hname << "<<" << std::endl;;
#endif
        // load ospray device
        ospout << "[ospray] Initialize OSPRay" << std::endl;	
        OSPDevice device = ospGetCurrentDevice();
	// check if ospray has been initialized already
        if (!device) {
	    ospout << "[ospray] device not found, creating one" << std::endl;
            device = ospNewDevice("default"); 
            if (DebugStream::Level5()) { 
                ospout << "[ospray] --> debug mode" << std::endl;
                ospDeviceSet1i(device, "debug", 0);
            }	
            if (numThreads > 0) {
                ospout << "[ospray] --> numThreads: " << numThreads 
                       << std::endl;
                ospDeviceSet1i(device, "numThreads", numThreads);
            }
            ospDeviceSetErrorFunc(device, OSPContext_ErrorFunc);
            ospDeviceSetStatusFunc(device, OSPContext_StatusFunc);
            ospDeviceCommit(device);
            ospSetCurrentDevice(device);
        }
        // load ospray module
        OSPError err = ospLoadModule("visit");
        if (err != OSP_NO_ERROR) {
	    osperr << "[Error] can't load visit module" << std::endl;
        }
        ospray_initialized = true;
    }
}

// ****************************************************************************
//
// OSPRay
//
// ****************************************************************************











// We use this function to minimize interface
void OSPVisItContext::Render(float xMin, float xMax, float yMin, float yMax,
                             int imgWidth, int imgHeight,
                             float*& dest, OSPVisItVolume* volume) 
{

    ospray::visit::Camera   cam(camera);
    ospray::visit::Renderer ren(renderer);
    
    cam.SetScreen(xMin, xMax, yMin, yMax);
    
    ren.Set(*(volume->patch.model));
    ren.Set(*camera);
        
    ospray::visit::FrameBuffer fb(volume->patch.fb);
    fb.Render(imgWidth, imgHeight,
	      cam.GetWindowExts(0),
	      cam.GetWindowExts(2),
	      ren->bgSize[0],
	      ren->bgDepthBuffer,
	      *ren,
	      dest);
}

void OSPVisItContext::InitPatch(int id) 
{
    if (volumes.find(id) == volumes.end()) {
        OSPVisItVolume v;
        v.parent = this;
        volumes[id] = v;
    }
}

// ****************************************************************************
//
// OSPVolume
//
// ****************************************************************************

void OSPVisItVolume::Set(int type, void *ptr, double *X, double *Y, double *Z, 
                         int nX, int nY, int nZ,
                         double volumePBox[6], 
                         double volumeBBox[6], 
                         double mtl[4], float sr,
                         bool shading)
{
    specularKs    = (float)mtl[2];
    specularNs    = (float)mtl[3];
    enableShading = shading;
    samplingRate  = sr;

    std::string str_type;
    OSPDataType osp_type;
    ospray::CheckVolumeFormat(type, str_type, osp_type);
    
    ospray::visit::Volume volume(patch.volume);
    volume.Init("visit_shared_structured_volume", osp_type, str_type,
		(size_t)nX * (size_t)nY * (size_t)nZ, ptr);
    volume.Set(false, false, false, false, shading, samplingRate,
	       specularKs, specularNs, X, Y, Z, nX, nY, nZ,
	       volumePBox, volumeBBox,
	       parent->bbox.upper,
	       parent->bbox.lower,
	       parent->regionScaling,
	       *(parent->tfn));
    
    ospray::visit::Model model(patch.model);
    model.Reset();
    model.Init();
    model.Set(*volume);
    
    finished = true;
}

/*
// ospFrameBuffer component     
void OSPVisItVolume::InitFB(unsigned int width, unsigned int height)
{
    // preparation
    imageSize.x = width;
    imageSize.y = height;

    std::vector<float> maxDepth(width * height);
    //
    // The reason I use round(r * (N-1)) instead of floor(r * N) is that 
    // during the composition phase, there will be a wired offset between
    // rendered image and the background, which is about one pixel in size.
    // Using round(r * (N - 1)) can remove the problem
    //
    // It seems this is the correct way of doing it
    //
    // It seems we need to also fix pan and zoom also
    //
    const int Xs = ((ospray::visit::Camera)parent->camera).GetWindowExts(0);
    const int Ys = ((ospray::visit::Camera)parent->camera).GetWindowExts(2);
    for (int i = 0; i < width; ++i) {
    	for (int j = 0; j < height; ++j) {
    	    maxDepth[i + j * width] = parent->renderer.bgDepthBuffer
                [Xs + i + (Ys + j) * parent->renderer.bgSize[0]];
    	}
    }
    framebufferBg = ospNewTexture2D(imageSize, OSP_TEXTURE_R32F, 
                                    maxDepth.data(),
                                    OSP_TEXTURE_FILTER_NEAREST);
    ospCommit(framebufferBg);
    ospSetObject(*(parent->renderer), "maxDepthTexture", framebufferBg);
    ospCommit(*(parent->renderer));
    ospRelease(framebufferBg);
    framebufferBg = NULL;

    CleanFB();
    framebuffer = ospNewFrameBuffer(imageSize, 
                                    OSP_FB_RGBA32F,
                                    OSP_FB_COLOR);
}
void OSPVisItVolume::RenderFB() {
    static int i = 0;
    ospRenderFrame(framebuffer, *(parent->renderer), OSP_FB_COLOR);
    framebufferData = (float*) ospMapFrameBuffer(framebuffer, OSP_FB_COLOR);
}
float* OSPVisItVolume::GetFBData() {
    return framebufferData;
}
*/
// ****************************************************************************
//
//
//
//  Extra Functions Defined here
//
//
//
// ****************************************************************************

void
ospray::CheckVolumeFormat(const int dt,
			  std::string& str_type,
			  OSPDataType& osp_type)
{
    if (dt == VTK_UNSIGNED_CHAR) {
	str_type = "uchar";
	osp_type = OSP_UCHAR;
    } else if (dt == VTK_SHORT) {
	str_type = "short";
	osp_type = OSP_SHORT;
    } else if (dt == VTK_UNSIGNED_SHORT) {
	str_type = "ushort";
	osp_type = OSP_USHORT;
    } else if (dt == VTK_FLOAT) {
	str_type = "float";
	osp_type = OSP_FLOAT;
    } else if (dt == VTK_DOUBLE) {
	str_type = "double";
	osp_type = OSP_DOUBLE;
    } else {
	ospray::Exception("ERROR: Unsupported ospray volume type");
    }
    ospout << "[ospray] data type " << str_type << std::endl;
}

void
ospray::ComputeProjections(const avtViewInfo &view, 
			   const double &aspect,
			   const int     screen[2],
			   const double  scale[3],
			   const double &oldNearPlane,
			   const double &oldFarPlane,
			   vtkMatrix4x4 *model_to_screen_transform, 
			   vtkMatrix4x4 *screen_to_model_transform, 
			   vtkMatrix4x4 *screen_to_camera_transform,
			   int           renderingExtents[4],
			   double        sceneSize[2],
			   double        dbounds[6])       
{
    vtkCamera *sceneCam = vtkCamera::New();
    sceneCam->SetPosition(view.camera[0],view.camera[1],view.camera[2]);
    sceneCam->SetFocalPoint(view.focus[0],view.focus[1],view.focus[2]);
    sceneCam->SetViewUp(view.viewUp[0],view.viewUp[1],view.viewUp[2]);
    sceneCam->SetViewAngle(view.viewAngle);
    sceneCam->SetClippingRange(oldNearPlane, oldFarPlane);
    if (view.orthographic) { sceneCam->ParallelProjectionOn(); }
    else { sceneCam->ParallelProjectionOff(); }
    sceneCam->SetParallelScale(view.parallelScale);	
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
    //sceneCam->Delete();
    
    // Get the full image extents of the volume
    double depthExtents[2];
    ospray::ProjectWorldToScreenCube(dbounds, screen[0], screen[1], 
				     view.imagePan, view.imageZoom,
				     model_to_screen_transform,
				     renderingExtents, depthExtents);
    ospray::ProjectWorldToScreenCube(dbounds, screen[0], screen[1], 
				     view.imagePan, view.imageZoom,
				     model_to_screen_transform,
				     renderingExtents, depthExtents);
    renderingExtents[0] = std::max(renderingExtents[0], 0);
    renderingExtents[2] = std::max(renderingExtents[2], 0);
    renderingExtents[1] = std::min(1+renderingExtents[1], screen[0]);
    renderingExtents[3] = std::min(1+renderingExtents[3], screen[1]);
    // Debug
    ospout << "[avrRayTracer] View settings: " << endl
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
	   << "[avrRayTracer] data bounds: "
	   << dbounds[0] << " " << dbounds[1] << std::endl
	   << "               data bounds  "
	   << dbounds[2] << " " << dbounds[3] << std::endl
	   << "               data bounds  "
	   << dbounds[4] << " " << dbounds[5] << std::endl
	   << "[avrRayTracer] full image extents: " 
	   << renderingExtents[0] << " " << renderingExtents[1] << std::endl
	   << "               full image extents: "
	   << renderingExtents[2] << " " << renderingExtents[3] << std::endl;
    
    ospout << "[avrRayTracer] model_to_screen_transform: " 
	   << *model_to_screen_transform << std::endl;
    ospout << "[avrRayTracer] screen_to_model_transform: " 
	   << *screen_to_model_transform << std::endl;
    ospout << "[avrRayTracer] screen_to_camera_transform: " 
	   << *screen_to_camera_transform << std::endl;

}

void
ospray::CheckMemoryHere(const std::string& message, std::string debugN)
{
    if (debugN.compare("ospout") == 0) {	
        ospray::CheckMemoryHere(message, *osp_out);
    }
    else if (debugN.compare("debug5") == 0) {
        if (DebugStream::Level5()) {
            ospray::CheckMemoryHere(message, DebugStream::Stream5());
        }       
    }
    else if (debugN.compare("debug4") == 0) {
        if (DebugStream::Level4()) {
            ospray::CheckMemoryHere(message, DebugStream::Stream4());
        }       
    }
    else if (debugN.compare("debug3") == 0) {
        if (DebugStream::Level3()) {
            ospray::CheckMemoryHere(message, DebugStream::Stream3());
        }       
    }
    else if (debugN.compare("debug2") == 0) {
        if (DebugStream::Level2()) {
            ospray::CheckMemoryHere(message, DebugStream::Stream2());
        }       
    }
    else if (debugN.compare("debug1") == 0) {
        if (DebugStream::Level1()) {
            ospray::CheckMemoryHere(message, DebugStream::Stream1());
        }       
    }
}

void
ospray::CheckMemoryHere(const std::string& message, std::ostream& out)
{
    unsigned long m_size, m_rss;
    avtMemory::GetMemorySize(m_size, m_rss);
    out << message << std::endl << "\t"
        << " Rank " << PAR_Rank()
        << " Memory use begin " << m_size 
        << " rss " << m_rss/(1024*1024) << " (MB)"
        << std::endl;
}

double
ospray::ProjectWorldToScreen(const double worldCoord[3], 
			     const int screenWidth, 
			     const int screenHeight,
			     const double panPercentage[2], 
			     const double imageZoom,
			     vtkMatrix4x4 *mvp, int screenCoord[2])
{
    // world space coordinate in homogeneous coordinate
    double worldHCoord[4] = {worldCoord[0],worldCoord[1],worldCoord[2],1.0};
    // world to clip space (-1 ~ 1)
    double clipHCoord[4];
    mvp->MultiplyPoint(worldHCoord, clipHCoord);
    // check error
    if (clipHCoord[3] == 0.0)
    {
        std::cerr << "world coordinates: (" 
                  << worldHCoord[0] << ", " 
                  << worldHCoord[1] << ", " 
                  << worldHCoord[2] << ", " 
                  << worldHCoord[3] << ")" << std::endl
                  << "clip space coordinate: ("
                  << clipHCoord[0] << ", " 
                  << clipHCoord[1] << ", " 
                  << clipHCoord[2] << ", "
                  << clipHCoord[3] << std::endl
		  << "Matrix: " << *mvp << std::endl;
	ospray::Exception("Zero Division During Projection");
    }
    // screen coordinates (int integer)
    screenCoord[0] =
	round((clipHCoord[0] / clipHCoord[3] + 1) * screenWidth  * 0.5) +
	round(screenWidth  * panPercentage[0]);
    screenCoord[1] =
	round((clipHCoord[1] / clipHCoord[3] + 1) * screenHeight * 0.5) +
	round(screenHeight * panPercentage[1]); 
    // return point depth
    return clipHCoord[2]/clipHCoord[3];
}

void
ospray::ProjectScreenToWorld(const int screenCoord[2], const double z,
                             const int screenWidth, const int screenHeight, 
                             const double panPercentage[2], 
                             const double imageZoom,
                             vtkMatrix4x4 *imvp, double worldCoord[3])
{
    // remove panning
    const int x = screenCoord[0] - round(screenWidth*panPercentage[0]);
    const int y = screenCoord[1] - round(screenHeight*panPercentage[1]);   
    // do projection
    double worldHCoord[4] = {0,0,0,1};
    double clipHCoord[4] = {
        (x - screenWidth  / 2.0) / (screenWidth  / 2.0),
        (y - screenHeight / 2.0) / (screenHeight / 2.0), z, 1.0};
    imvp->MultiplyPoint(clipHCoord, worldHCoord);
    // check error
    if (worldHCoord[3] == 0) {
        std::cerr << "world coordinates: (" 
                  << worldHCoord[0] << ", " 
                  << worldHCoord[1] << ", " 
                  << worldHCoord[2] << ", " 
                  << worldHCoord[3] << ")" << std::endl
                  << "clip space coordinate: ("
                  << clipHCoord[0] << ", " 
                  << clipHCoord[1] << ", " 
                  << clipHCoord[2] << ", "
                  << clipHCoord[3] << std::endl
		  << "Matrix: " << *imvp << std::endl;
	ospray::Exception("Zero Division During Projection");
    }    
    // normalize world space coordinate	
    worldCoord[0] = worldHCoord[0]/worldHCoord[3];
    worldCoord[1] = worldHCoord[1]/worldHCoord[3];
    worldCoord[2] = worldHCoord[2]/worldHCoord[3];
}

void
ospray::ProjectScreenToCamera(const int screenCoord[2], const double z,
                              const int screenWidth, const int screenHeight, 
                              vtkMatrix4x4 *imvp, double cameraCoord[3])
{
    // remove panning
    const int x = screenCoord[0];
    const int y = screenCoord[1];
    // do projection
    double cameraHCoord[4] = {0,0,0,1};
    double clipHCoord[4] = {
        (x - screenWidth /2.0)/(screenWidth /2.0),
        (y - screenHeight/2.0)/(screenHeight/2.0), z, 1.0};
    imvp->MultiplyPoint(clipHCoord, cameraHCoord);
    // check error
    if (cameraHCoord[3] == 0) {
        std::cerr << "world coordinates: (" 
                  << cameraHCoord[0] << ", " 
                  << cameraHCoord[1] << ", " 
                  << cameraHCoord[2] << ", " 
                  << cameraHCoord[3] << ")" << std::endl
                  << "clip space coordinate: ("
                  << clipHCoord[0] << ", " 
                  << clipHCoord[1] << ", " 
                  << clipHCoord[2] << ", "
                  << clipHCoord[3] << std::endl
		  << "Matrix: " << *imvp << std::endl;
	ospray::Exception("Zero Division During Projection");
    }
    // normalize world space coordinate	
    cameraCoord[0] = cameraHCoord[0]/cameraHCoord[3];
    cameraCoord[1] = cameraHCoord[1]/cameraHCoord[3];
    cameraCoord[2] = cameraHCoord[2]/cameraHCoord[3];
}

void
ospray::ProjectWorldToScreenCube(const double cube[6],
                                 const int screenWidth, 
                                 const int screenHeight, 
                                 const double panPercentage[2], 
                                 const double imageZoom,
                                 vtkMatrix4x4 *mvp, 
                                 int screenExtents[4], 
                                 double depthExtents[2])
{
    int xMin = std::numeric_limits<int>::max();
    int xMax = std::numeric_limits<int>::min();
    int yMin = std::numeric_limits<int>::max();
    int yMax = std::numeric_limits<int>::min();
    double zMin = std::numeric_limits<double>::max();
    double zMax = std::numeric_limits<double>::min();

    float coordinates[8][3];
    coordinates[0][0] = cube[0];   
    coordinates[0][1] = cube[2];   
    coordinates[0][2] = cube[4];	

    coordinates[1][0] = cube[1];   
    coordinates[1][1] = cube[2];   
    coordinates[1][2] = cube[4];	

    coordinates[2][0] = cube[1];  
    coordinates[2][1] = cube[3];
    coordinates[2][2] = cube[4];	

    coordinates[3][0] = cube[0]; 
    coordinates[3][1] = cube[3]; 
    coordinates[3][2] = cube[4];

    coordinates[4][0] = cube[0];
    coordinates[4][1] = cube[2];
    coordinates[4][2] = cube[5];

    coordinates[5][0] = cube[1]; 
    coordinates[5][1] = cube[2]; 
    coordinates[5][2] = cube[5];	

    coordinates[6][0] = cube[1]; 
    coordinates[6][1] = cube[3];
    coordinates[6][2] = cube[5];

    coordinates[7][0] = cube[0]; 
    coordinates[7][1] = cube[3]; 
    coordinates[7][2] = cube[5];

    double worldCoord[3];
    int screenCoord[2]; double depth;
    for (int i=0; i<8; i++)
    {
        worldCoord[0] = coordinates[i][0];
        worldCoord[1] = coordinates[i][1];
        worldCoord[2] = coordinates[i][2];
        depth = ospray::ProjectWorldToScreen
            (worldCoord, screenWidth, screenHeight, 
             panPercentage, imageZoom, mvp, screenCoord);
        // clamp values
        screenCoord[0] = CLAMP(screenCoord[0], 0, screenWidth);
        screenCoord[1] = CLAMP(screenCoord[1], 0, screenHeight);
        screenExtents[0] = xMin = std::min(xMin, screenCoord[0]);
        screenExtents[1] = xMax = std::max(xMax, screenCoord[0]);
        screenExtents[2] = yMin = std::min(yMin, screenCoord[1]);
        screenExtents[3] = yMax = std::max(yMax, screenCoord[1]);
        depthExtents[0] = zMin = std::min(zMin, depth);
        depthExtents[1] = zMax = std::max(zMax, depth);
    }
}

void
ospray::CompositeBackground(int screen[2],
                            int compositedImageExtents[4],
                            int compositedImageWidth,
                            int compositedImageHeight,
                            float *compositedImageBuffer,
                            unsigned char *opaqueImageColor,
                            float         *opaqueImageDepth,
                            unsigned char *&imgFinal)
{
    if (UseThreadedBlend_MetaData) {
        ospray::visit::CompositeBackground(screen,
                                           compositedImageExtents,
                                           compositedImageWidth,
                                           compositedImageHeight,
                                           compositedImageBuffer,
                                           opaqueImageColor,
                                           opaqueImageDepth,
                                           imgFinal);
	return;
    } 
    for (int y = 0; y < screen[1]; y++)
    {
	for (int x = 0; x < screen[0]; x++)
	{
	    int indexScreen     = y * screen[0] + x;
	    int indexComposited =
		(y - compositedImageExtents[2]) * compositedImageWidth +
		(x - compositedImageExtents[0]);

	    bool insideComposited = 
		((x >= compositedImageExtents[0] && 
		  x < compositedImageExtents[1]) &&
		 (y >= compositedImageExtents[2] && 
		  y < compositedImageExtents[3]));

	    if (insideComposited)
	    {
		if (compositedImageBuffer[indexComposited*4 + 3] == 0)
		{
		    // No data from rendering here! - Good
		    imgFinal[indexScreen * 3 + 0] = 
			opaqueImageColor[indexScreen * 3 + 0];
		    imgFinal[indexScreen * 3 + 1] = 
			opaqueImageColor[indexScreen * 3 + 1];
		    imgFinal[indexScreen * 3 + 2] = 
			opaqueImageColor[indexScreen * 3 + 2];
		}
		else
		{
		    // Volume in front
		    float alpha = 
			(1.0 - compositedImageBuffer[indexComposited * 4 + 3]);
		    imgFinal[indexScreen * 3 + 0] = 
			CLAMP(opaqueImageColor[indexScreen * 3 + 0] * alpha +
			      compositedImageBuffer[indexComposited * 4 + 0] *
			      255.f,
			      0.f, 255.f);
		    imgFinal[indexScreen * 3 + 1] = 
			CLAMP(opaqueImageColor[indexScreen * 3 + 1] * alpha +
			      compositedImageBuffer[indexComposited * 4 + 1] *
			      255.f,
			      0.f, 255.f);
		    imgFinal[indexScreen * 3 + 2] =
			CLAMP(opaqueImageColor[indexScreen * 3 + 2] * alpha +
			      compositedImageBuffer[indexComposited * 4 + 2] *
			      255.f,
			      0.f, 255.f);
		}
	    }
	    else
	    {
		// Outside bounding box: Use the background : Good
		imgFinal[indexScreen * 3 + 0] = 
		    opaqueImageColor[indexScreen * 3 + 0];
		imgFinal[indexScreen * 3 + 1] =
		    opaqueImageColor[indexScreen * 3 + 1];
		imgFinal[indexScreen * 3 + 2] =
		    opaqueImageColor[indexScreen * 3 + 2];
	    }
	}
    }
}

void
ospray::WriteArrayToPPM(std::string filename, const float * image,
                        int dimX, int dimY)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
                             std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n"; 
    for (int y=dimY-1; y>=0; --y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x)*4;
            char color[3];
            float alpha = image[index + 3];
            color[0] = CLAMP(image[index + 0]*alpha, 0.0f, 1.0f) * 255;
            color[1] = CLAMP(image[index + 1]*alpha, 0.0f, 1.0f) * 255;
            color[2] = CLAMP(image[index + 2]*alpha, 0.0f, 1.0f) * 255;
            outputFile.write(color,3);
        }
    } 
    outputFile.close();
}

void
ospray::WriteArrayToPPM(std::string filename,  const unsigned char *image, 
                        int dimX, int dimY)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
                             std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n"; 
    for (int y=dimY-1; y>=0; --y)
    {
        outputFile.write(reinterpret_cast<const char*>(&image[y * dimX * 3]), 
                         dimX * 3);
    } 
    outputFile.close();
}

void
ospray::WriteArrayGrayToPPM(std::string filename, const float* image, 
                            int dimX, int dimY)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
                             std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n"; 
    for (int y=dimY-1; y>=0; --y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x);
            char var = CLAMP(image[index], 0.f, 1.f) * 255;
            char color[3];
            color[0] = var;
            color[1] = var;
            color[2] = var;
            outputFile.write(color,3);
        }
    } 
    outputFile.close();
}
