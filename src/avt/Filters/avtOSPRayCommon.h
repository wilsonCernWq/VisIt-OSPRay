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

#ifndef AVT_OSPRAY_COMMON_H
#define AVT_OSPRAY_COMMON_H

#ifndef VISIT_OSPRAY /* make sure VISIT_OSPRAY is defined */
# error "VISIT_OSPRAY is not defined but ospray is used"
#endif

#include <avtParallel.h>

#include <vtkType.h>
#include <vtkMatrix4x4.h>

#include <ospray/ospray.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>

// some constants
#define OSP_PERSPECTIVE              1
#define OSP_ORTHOGRAPHIC             2
#define OSP_BLOCK_BRICKED_VOLUME     3
#define OSP_SHARED_STRUCTURED_VOLUME 4
#define OSP_INVALID                  5
#define OSP_VALID                    6

// ostreams customized for ospray
#ifdef ospout
#undef ospout
#endif
#define ospout \
    if (!ospray::CheckVerbose() && !DebugStream::Level5()) ; \
    else (*ospray::osp_out)
#ifdef osperr
#undef osperr
#endif
#define osperr \
    if (!ospray::CheckVerbose() && !DebugStream::Level1()) ; \
    else (*ospray::osp_err)
namespace ospray {
    extern std::ostream *osp_out;
    extern std::ostream *osp_err;
    //////////////////////////////////////////////////////
    //                                                  //
    // Those function has to be inline, otherwise we    //
    // need to link this library to other components    //
    // manually                                         //
    //                                                  //
    //////////////////////////////////////////////////////
    //                                                  //
    // OSPRay defines following environmental variables //
    //                                                  //
    // OSPRAY_DEBUG                                     //
    // OSPRAY_THREADS                                   //
    // OSPRAY_LOG_LEVEL                                 //
    // OSPRAY_LOG_OUTPUT                                //
    // OSPRAY_SET_AFFINITY                              //
    //                                                  //
    // We define one more environmental variable here   //
    //                                                  //
    // OSPRAY_VERBOSE                                   //
    //                                                  //
    //////////////////////////////////////////////////////
    inline bool InitVerbose() {
        const char* env_verbose   = std::getenv("OSPRAY_VERBOSE");
        const char* env_debug     = std::getenv("OSPRAY_DEBUG");
        const char* env_log_level = std::getenv("OSPRAY_LOG_LEVEL");	
        bool verbose = false;
        if (env_verbose) { if (atoi(env_verbose) > 0) { verbose = true; } }
        if (env_debug) { if (atoi(env_debug) > 0) { verbose = true; } }
        if (env_log_level) { if (atoi(env_log_level) > 0) { verbose = true; } }
        if (verbose) {
            ospray::osp_out = &std::cout;
            ospray::osp_err = &std::cerr;
            return true;
	} else { return false; }
    }
    inline int InitOSPRaySpp() {
        int spp = 1;
        const char* env_spp = std::getenv("OSPRAY_SPP");
        if (env_spp) { if (atoi(env_spp) > 0) { spp = atoi(env_spp); } }
        return spp;
    }
    inline bool CheckVerbose() { // initialize OSPRAY_VERBOSE    
        static bool OSPRAY_VERBOSE = ospray::InitVerbose();
        return OSPRAY_VERBOSE;
    }
    inline int CheckOSPRaySpp() {
        static int spp = InitOSPRaySpp();
        return spp;
    }
};

// ****************************************************************************
//  Struct:  OSPVisItVolume
//
//  Purpose:
//    
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

class OSPVisItContext;
class OSPVisItVolume 
{
 private:
    friend class OSPVisItContext;

    OSPVisItContext *parent;

    // objects owned by the struct
    // -- ospray model ---
    OSPModel            world;
    unsigned char       worldType;
    // --- ospray framebuffer ---
    OSPFrameBuffer      framebuffer;
    float              *framebufferData;
    OSPTexture2D        framebufferBg;
    osp::vec2i imageSize;
    // --- ospray volume ---
    OSPVolume           volume;
    unsigned char       volumeType;
    // --- ospray data ---
    OSPDataType         voxelDataType;
    OSPData             voxelData;
    size_t              voxelSize;
    void*               dataPtr;
    std::string         dataType;

    // metadata for volume
    int                 patchId;       // volume patch id
    bool                finished;      // check if this volume is initialized
    bool                enableShading;
    bool                enableDVR;     // Distributed Volume Renderer
    float               specularKs;
    float               specularNs;
    float               samplingRate;

    // geometric parameters for volume
    osp::vec3i          regionSize;
    osp::vec3f          regionStart;
    osp::vec3f          regionStop;
    osp::vec3f          regionSpacing;
    osp::vec3f          regionUpperClip;
    osp::vec3f          regionLowerClip;
    osp::vec3f          regionScaling;
    
 public:
    // constructor
    OSPVisItVolume() {
	// objects owned by the struct
	world           = NULL;
	worldType       = OSP_INVALID;
	framebuffer     = NULL;
	framebufferData = NULL;
	framebufferBg   = NULL;
	volume          = NULL;
	volumeType      = OSP_INVALID;
	voxelDataType   = OSP_VOID_PTR;
	voxelData       = NULL;
	voxelSize       = 0;
	dataPtr         = NULL;
	dataType        = "";
	// metadata for volume
	//patchId = id;    
	finished      = false; 
	enableShading = false;
	enableDVR     = false;
	specularKs    = 1.0f;
	specularNs    = 15.0f;
	samplingRate  = 3.0f;
	// geometric parameters for volume
	regionSize.x  = regionSize.y  = regionSize.z  = 0;
	regionStart.x = regionStart.y = regionStart.z = 0.0f;
	regionStop.x  = regionStop.y  = regionStop.z  = 0.0f;
	regionSpacing.x   = regionSpacing.y   = regionSpacing.z   = 0.0f;
	regionUpperClip.x = regionUpperClip.y = regionUpperClip.z = 0.0f;
	regionLowerClip.x = regionLowerClip.y = regionLowerClip.z = 0.0f;
	regionScaling.x   = regionScaling.y   = regionScaling.z   = 1.0f;
    }

    // destructor
    ~OSPVisItVolume() { Clean(); }    
    void Clean() {
        CleanFB();
        CleanVolume();	
        CleanWorld();
    }
    
    // other function
    void Set(int type, void *ptr, 
             double *X, double *Y, double *Z, 
             int nX, int nY, int nZ, 
             double volumePBox[6], double volumeBBox[6],
             double mtl[4], float sr, bool shading);
    bool GetDVRFlag() { return enableDVR; }
    void SetDVRFlag(bool mode) { enableDVR = mode; }
    bool GetFinishedFlag() { return finished; }
    void SetFinishedFlag(bool f) { finished = f; } 

    // ospModel component
    OSPModel GetWorld() { return world; }
    void InitWorld();
    void SetWorld();
    void CleanWorld() {
        if (world != NULL) {	    
            ospRelease(world);
            world = NULL;
        }
        worldType = OSP_INVALID;
    }
	
    // ospVolume component
    void InitVolume(unsigned char type = OSP_SHARED_STRUCTURED_VOLUME); 
    OSPVolume GetVolume() { return volume; }
    void SetVolume(int type, void *ptr,
                   double *X, double *Y, double *Z, 
                   int nX, int nY, int nZ,
                   double volumePBox[6], 
                   double volumeBBox[6]);
    void CleanVolume() {	
        if (volume != NULL) { ospRelease(volume); volume = NULL; }
        if (voxelData != NULL) { 
            ospRelease(voxelData);
            voxelData = NULL; 
        }
        volumeType = OSP_INVALID;
    }

    // framebuffer component     
    void InitFB(unsigned int width, unsigned int height);
    void RenderFB();
    float* GetFBData();
    void CleanFB() {
        if (framebufferData != NULL) { 
            ospUnmapFrameBuffer(framebufferData, framebuffer); 
            framebufferData = NULL;
        }
        if (framebuffer != NULL) { 
            ospRelease(framebuffer); 	    
            framebuffer = NULL;
        }
        if (framebufferBg != NULL) {
            ospRelease(framebufferBg); 	    
            framebufferBg = NULL;
        }
    }
};

// ****************************************************************************
//  Struct:  OSPVisItLight
//
//  Purpose:
//
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

struct OSPVisItLight
{
    OSPLight aLight;
    OSPLight dLight;
    OSPLight sLight; // constant sun light
    OSPData  lightdata;
    OSPVisItLight() {
        aLight = NULL;
        dLight = NULL;
        sLight = NULL;
        lightdata = NULL;
    }
    ~OSPVisItLight() { Clean(); }
    void Clean() {/* TODO should we delete them? */}
    void Init(const OSPRenderer& renderer);
    void Set(double materialProperties[4], double viewDirection[3]);
};

// ****************************************************************************
//  Struct:  OSPVisItRenderer
//
//  Purpose:
//
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

struct OSPVisItRenderer
{
public:
    enum State {
        INVALID, /* TODO do we need this actually ? */
        SCIVIS,
    } rendererType;
    OSPRenderer renderer;
    OSPVisItLight lights;
    // properties
    int aoSamples;
    int spp; //!< samples per pixel
    bool flagOneSidedLighting;
    bool flagShadowsEnabled;
    bool flagAoTransparencyEnabled;
    float       *maxDepthBuffer;  // depth buffer (shared, never delete)
    osp::vec2i   maxDepthSize;    // buffer extents (minX, maxX, minY, max)  
public:
    OSPVisItRenderer() {
        renderer = NULL;
        rendererType = INVALID;
        aoSamples = 0;
        spp = ospray::CheckOSPRaySpp();
        flagOneSidedLighting = false;
        flagShadowsEnabled = false;
        flagAoTransparencyEnabled = false;
    }
    ~OSPVisItRenderer() { Clean(); }
    void Clean() {
        if (renderer != NULL) {
            lights.Clean();
            ospRelease(renderer);
            renderer = NULL;
            rendererType = INVALID;
        }
    }
    void Init();
    void Set(double materialProperties[4], double viewDirection[3], bool);
    void SetCamera(const OSPCamera& camera);
    void SetModel(const OSPModel& world);
};

// ****************************************************************************
//  Struct:  OSPVisItCamera
//
//  Purpose:
//
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

struct OSPVisItCamera
{
public:
    enum State {
        INVALID,
        PERSPECTIVE,
        ORTHOGRAPHIC,
    } cameraType;
    OSPCamera camera;
    float panx; // this is a ratio [0, 1]
    float pany; // this is a ratio [0, 1]
    float r_xl; 
    float r_yl;
    float r_xu;
    float r_yu;
    float zoom; 
    int   size[2];
    osp::vec2f imgS, imgE;
public:
    OSPVisItCamera() {
        camera = NULL;
        cameraType = INVALID;
        panx = 0.0f;
        pany = 0.0f;
        zoom = 1.0f;
        size[0] = size[1] = 0.0f;
        imgS.x = 0.f;
        imgS.y = 0.f;
        imgE.x = 0.f;
        imgE.y = 0.f;
    }
    ~OSPVisItCamera() { Clean(); }
    void Clean() {
        if (camera != NULL) {
            ospRelease(camera);
            camera = NULL;
            cameraType = INVALID;
        }
    }
    void Init(State type);
    void Set(const double camp[3], 
             const double camf[3], 
             const double camu[3], 
             const double camd[3],
             const double sceneSize[2],
             const double aspect, 
             const double fovy, 
             const double zoom_ratio, 
             const double pan_ratio[2],
             const int bufferExtents[4],
             const int screenExtents[2]);
    void SetScreen(float xMin, float xMax, float yMin, float yMax);
};

// ****************************************************************************
//  Struct:  OSPVisItColor
//
//  Purpose:
//
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

struct OSPVisItColor { float R,G,B, A; };

// ****************************************************************************
//  Struct:  OSPVisItTransferFunction
//
//  Purpose:
//
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

struct OSPVisItTransferFunction
{
public:
    enum State { INVALID, PIECEWISE_LINEAR, } transferfcnType;
    OSPTransferFunction  transferfcn;
public:
    OSPVisItTransferFunction() {
        transferfcn = NULL;
        transferfcnType = INVALID;
    }
    void Clean() {
        if (transferfcn != NULL) {
            ospRelease(transferfcn);
            transferfcn = NULL;
            transferfcnType = INVALID;
        }
    }
    void Init();
    void Set(const OSPVisItColor* table,
             const unsigned int size, 
             const float datamin,
             const float datamax);
};

// ****************************************************************************
//  Struct:  OSPVisItContext
//
//  Purpose:
//
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

class OSPVisItContext
{
public:
    // ************************************************************************
    // We expose this in header because iy will be called in other components
    // where we dont have direct library linkage
    // ************************************************************************
    OSPVisItContext() 
    {
        regionScaling.x = regionScaling.y = regionScaling.z = 1.0f;
    }
    ~OSPVisItContext() {	
        volumes.clear();
        renderer.Clean();
        camera.Clean();
        transferfcn.Clean();
    }

    // helper
    void Render(float xMin, float xMax, float yMin, float yMax,
                int imgWidth, int imgHeight, 
                float*& dest, OSPVisItVolume* volume);
    void InitOSP(int numThreads = 0);
    void Finalize();
    void InitPatch(int id);
    OSPVisItVolume* GetPatch(int id) { return &volumes[id]; }

    // parameters
    void SetDataBounds(double dbounds[6]) {
        for (int i = 0; i < 6; ++i) { bounds[i] = dbounds[i]; }
    }
    void SetBgBuffer(float* depth, int extents[4]) {
        renderer.maxDepthBuffer = depth;
        renderer.maxDepthSize.x = extents[1] - extents[0];
        renderer.maxDepthSize.y = extents[3] - extents[2];
    }
    void SetScaling(double s[3]) { 
        regionScaling.x = (float)s[0];
        regionScaling.y = (float)s[1];
        regionScaling.z = (float)s[2]; 
    }

    OSPVisItRenderer renderer;
    OSPVisItCamera   camera;
    OSPVisItTransferFunction transferfcn;
    std::map<int, OSPVisItVolume> volumes;

private:
    
    friend class OSPVisItVolume;
    osp::vec3f     regionScaling;
    double bounds[6];
    static bool initialized;
    
};

#endif//AVT_OSPRAY_COMMON_H

// ****************************************************************************
//
//
//
//  Extra Functions Defined here
//
//
//
// ****************************************************************************

#ifndef VISIT_OSPRAY_CONTEXT_ONLY

#ifndef AVT_OSPRAY_COMMON_EXTRA_H
#define AVT_OSPRAY_COMMON_EXTRA_H

#include <DebugStream.h>
#include <StackTimer.h>
#include <TimingsManager.h>
#include <ImproperUseException.h>
#include <ospray/visit/VisItImageComposite.h>

#ifndef CLAMP
# define CLAMP(x, l, h) (x > l ? x < h ? x : h : l)
#endif
#ifndef M_MIN
# define M_MIN(x, r) (x < r ? x : r)
#endif
#ifndef M_MAX
# define M_MAX(x, r) (x > r ? x : r)
#endif

// ****************************************************************************
//  Namespace:  ospray
//
//  Purpose:
//    
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************

namespace ospray
{
    // ************************************************************************
    //  Struct:  ImgMetaData
    //
    //  Purpose:
    //    Holds information about patches but not the image 
    //
    //  Programmer:  
    //  Creation:   
    //
    // ************************************************************************

    struct ImgMetaData
    {
        int procId;       // processor that produced the patch
        int patchNumber;  // id of the patch on that processor
        int destProcId;   // destination proc where this patch gets composited
        int inUse;        // whether the patch is composed locally or not
        int dims[2];      // height, width
        int screen_ll[2]; // (lower left)  position in the final image
        int screen_ur[2]; // (upper right)
        float avg_z;      // camera space depth of the patch (average)
        float eye_z;      // camera space z
        float clip_z;     // clip space z
    };

    // ************************************************************************
    //  Struct:  ImgData
    //
    //  Purpose:
    //    Holds the image data generated
    //
    //  Programmer:  
    //  Creation:    
    //
    // ************************************************************************
    
    struct ImgData
    {
        // acts as a key
        int procId;        // processor that produced the patch
        int patchNumber;   // id of the patch on that processor
        float *imagePatch; // the image data - RGBA
        ImgData() { imagePatch = NULL; }
        bool operator==(const ImgData &a) {
            return (patchNumber == a.patchNumber);
        }
    };

    // ************************************************************************
    //
    //  Helper Functions
    //
    // ************************************************************************
        
    void CheckMemoryHere(const std::string& message, 
                         std::string debugN = "debug5");
    void CheckMemoryHere(const std::string& message, 
                         std::ostream& out);

    typedef int timestamp;
    inline void CheckSectionStart(const std::string& c,
                                  const std::string& f,
                                  timestamp& timingDetail,
                                  const std::string& str) 
    {
        debug5 << c << "::" << f << " " << str << " Start" << std::endl;
        timingDetail = visitTimer->StartTimer();	    
    }
    
    inline void CheckSectionStop(const std::string& c,
                                 const std::string& f, 
                                 timestamp& timingDetail,
                                 const std::string& str) 
    {
        visitTimer->StopTimer(timingDetail, 
                              (c + "::" + f + " " + str).c_str());
        ospray::CheckMemoryHere(("[" + c + "]" + " " + f + " " + str).c_str(),
                                "debug5");
        debug5 << c << "::" << f << " " << str << " Done" << std::endl;
    }

    inline void Exception(const std::string str)
    {
        std::cerr << str << std::endl;
        debug1    << str << std::endl;
        EXCEPTION1(VisItException, str.c_str()); 
    }

    double ProjectWorldToScreen
        (const double worldCoord[3], 
         const int screenWidth, const int screenHeight,	 
         const double panPercentage[2], const double imageZoom,
         vtkMatrix4x4 *mvp, int screenCoord[2]);
    
    void ProjectScreenToWorld
        (const int screenCoord[2], const double z,
         const int screenWidth, const int screenHeight, 
         const double panPercentage[2], const double imageZoom,
         vtkMatrix4x4 *imvp, double worldCoord[3]);

    void ProjectScreenToCamera
        (const int screenCoord[2], const double z,
         const int screenWidth, const int screenHeight,
         vtkMatrix4x4 *imvp, double cameraCoord[3]);

    inline void ProjectScreenToWorld
        (const int x, const int y, const double z,
         const int screenWidth, const int screenHeight, 
         const double panPercentage[2], const double imageZoom,
         vtkMatrix4x4 *imvp, double worldCoord[3]) 
    {
        int screen_coord[2] = {x, y};
        ProjectScreenToWorld(screen_coord, z, screenWidth, screenHeight, 
                             panPercentage, imageZoom, imvp, worldCoord);
    }

    void ProjectWorldToScreenCube
        (const double cube[6], const int screenWidth, const int screenHeight, 
         const double panPercentage[2], const double imageZoom, 
         vtkMatrix4x4 *mvp,int screenExtents[4], double depthExtents[2]);

    void CompositeBackground(int screen[2],
                             int compositedImageExtents[4],
                             int compositedImageWidth,
                             int compositedImageHeight,
                             float *compositedImageBuffer,
                             unsigned char *opaqueImageColor,
                             float         *opaqueImageDepth,
                             unsigned char *&imgFinal);
    
    void WriteArrayToPPM
        (std::string filename, const float *image, int dimX, int dimY);

    void WriteArrayToPPM
        (std::string filename, const unsigned char *image, int dimX, int dimY);

    void WriteArrayGrayToPPM
        (std::string filename, const float * image, int dimX, int dimY);
};

#endif//AVT_OSPRAY_COMMON_EXTRA_H

#endif//VISIT_OSPRAY_CONTEXT_ONLY

