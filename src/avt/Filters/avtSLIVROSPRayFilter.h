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

#ifndef AVT_OSPRAY_FILTER_H
#define AVT_OSPRAY_FILTER_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include <vtkType.h>
#include <DebugStream.h>

#ifdef VISIT_OSPRAY
# include "ospray/ospray.h"
# define OSP_PERSPECTIVE              1
# define OSP_ORTHOGRAPHIC             2
# define OSP_BLOCK_BRICKED_VOLUME     3
# define OSP_SHARED_STRUCTURED_VOLUME 4
# define OSP_INVALID                  5
# define OSP_VALID                    6
#endif

// ****************************************************************************
// Debug ostreams customized for ospray
// ****************************************************************************
namespace slivr 
{
    extern std::ostream *osp_out;
    extern std::ostream *osp_err;
};

// ****************************************************************************
// Those function has to be inline, otherwise we need to link this library 
// to other components manually
// ****************************************************************************
namespace slivr 
{
    // ************************************************************************
    //
    // Detect verbose from environmental variable
    //
    // ************************************************************************
    inline bool InitVerbose() 
    {
#ifndef VISIT_OSPRAY
	return false;
#else
        // ********************************************************************
	//
	// OSPRay defines following environmental variables
	//
	// OSPRAY_THREADS
	// OSPRAY_SET_AFFINITY
	// OSPRAY_DEBUG
	// OSPRAY_LOG_LEVEL
	// OSPRAY_LOG_OUTPUT
	//
	// We define one more environmental variable here
	//
	// OSPRAY_VERBOSE
	//
        // ********************************************************************
	const char* env_verbose   = std::getenv("OSPRAY_VERBOSE");
	const char* env_debug     = std::getenv("OSPRAY_DEBUG");
	const char* env_log_level = std::getenv("OSPRAY_LOG_LEVEL");	
	bool verbose = false;
	if (env_verbose) {
	    if (atoi(env_verbose) > 0) { verbose = true; }
	}
	if (env_debug) {
	    if (atoi(env_debug) > 0) { verbose = true; }
	}
	if (env_log_level) {
	    if (atoi(env_log_level) > 0) { verbose = true; }
	}
	if (verbose) {
	    slivr::osp_out = &std::cout;
	    slivr::osp_err = &std::cerr;
	    return true;
	} else {
	    return false;
	}
#endif
    }

    // ************************************************************************
    //
    // Detect sample per pixel from environmental variable
    //
    // ************************************************************************
    inline int InitOSPRaySpp() {
#ifndef VISIT_OSPRAY
	return 1;
#else
	int spp = 1;
	const char* env_spp = std::getenv("OSPRAY_SPP");
	if (env_spp) {
	    if (atoi(env_spp) > 0) { 
		spp = atoi(env_spp); 
	    }
	}	
	return spp;
#endif
    }
    inline bool CheckVerbose() // initialize OSPRAY_VERBOSE
    {
	static bool OSPRAY_VERBOSE = slivr::InitVerbose();
	return OSPRAY_VERBOSE;
    }
    inline int CheckOSPRaySpp()
    {
	static int spp = InitOSPRaySpp();
	return spp;
    }
};

// ****************************************************************************
//
// Over-write ostream marcos
//
// ****************************************************************************
#define ospout \
    if (!slivr::CheckVerbose() && !DebugStream::Level5()) ; \
    else (*slivr::osp_out)
#define osperr \
    if (!slivr::CheckVerbose() && !DebugStream::Level1()) ; \
    else (*slivr::osp_err)

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
#ifdef VISIT_OSPRAY
class OSPVisItContext;
class OSPVisItVolume 
{
 private:
    friend class OSPVisItContext;
 private:
    OSPVisItContext *parent;
    /* // object references */
    /* // (those objects are created and handled by other parts of the program) */
    /* // (shouldnt be deleted in this struct) */
    /* OSPTransferFunction transferfcn; */
    /* OSPRenderer         renderer; */
    /* float         *bgDepthBuffer; // depth buffer for the background and other plots */
    /* unsigned char *bgColorBuffer; // bounding box + pseudo color + ... */
    /* // int            bgExtents[4];  // extents of the buffer(minX, maxX, minY, max) */

    // objects owned by the struct
    // -- ospray model ---
    OSPModel            world;
    unsigned char       worldType;
    // --- ospray framebuffer ---
    OSPFrameBuffer      framebuffer;
    float              *framebufferData;
    OSPTexture2D        framebufferBg;
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
    bool                finished;    // check if this volume is initialized
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
    OSPVisItVolume(int id) {
	// object references 
	//transferfcn = NULL;
	//renderer    = NULL;
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
	patchId = id;    
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
    void Set(int type, void *ptr, double *X, double *Y, double *Z, 
	     int nX, int nY, int nZ, double volumePBox[6], double volumeBBox[6],
	     double mtl[4], float sr, bool shading);
    bool GetDVRFlag() { return enableDVR; }
    void SetDVRFlag(bool mode) { enableDVR = mode; }
    bool GetFinishedFlag() { return finished; }
    void SetFinishedFlag(bool f) { finished = f; } 
    //void SetScaling(const osp::vec3f& s) { regionScaling = s; }
    //void SetTransferFunction(const OSPTransferFunction& t) { transferfcn = t; }
    //void SetRenderer(const OSPRenderer& r) { renderer = r; }
    //void SetBgBuffer(unsigned char* color, float* depth, int extents[4]) {
	//bgColorBuffer = color;
	//bgDepthBuffer = depth;
	/* bgExtents[0] = extents[0]; */
	/* bgExtents[1] = extents[1]; */
	/* bgExtents[2] = extents[2]; */
	/* bgExtents[3] = extents[3]; */
    //}

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
    /* void CleanFBData() { */
    /* 	if (framebufferData != NULL) {  */
    /* 	    ospUnmapFrameBuffer(framebufferData, framebuffer);  */
    /* 	    framebufferData = NULL; */
    /* 	}	 */
    /* } */
    void CleanFB() {
	if (framebufferData != NULL) { 
	    ospUnmapFrameBuffer(framebufferData, framebuffer); 
	    framebufferData = NULL;
	}
	if (framebuffer != NULL) { 
	    ospFreeFrameBuffer(framebuffer); 	    
	    framebuffer = NULL;
	}
	if (framebufferBg != NULL) {
	    ospRelease(framebufferBg); 	    
	    framebufferBg = NULL;
	}
    }
};
#endif//VISIT_OSPRAY


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

#ifdef VISIT_OSPRAY
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
#endif//VISIT_OSPRAY


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

#ifdef VISIT_OSPRAY
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
    OSPTexture2D maxDepthTexture;
public:
    OSPVisItRenderer() {
	renderer = NULL;
	rendererType = INVALID;
	aoSamples = 0;
	spp = slivr::CheckOSPRaySpp();
	flagOneSidedLighting = false;
	flagShadowsEnabled = false;
	flagAoTransparencyEnabled = false;
	maxDepthTexture = NULL;
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
#endif//VISIT_OSPRAY


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

#ifdef VISIT_OSPRAY
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
    float zoom; 
    int   size[2];

public:
    OSPVisItCamera() {
	camera = NULL;
	cameraType = INVALID;
	panx = 0.0f;
	pany = 0.0f;
	zoom = 1.0f;
	size[0] = size[1] = 0.0f;
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
#endif//VISIT_OSPRAY


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

#ifdef VISIT_OSPRAY
struct OSPVisItColor { float R,G,B, A; };
#endif//VISIT_OSPRAY

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

#ifdef VISIT_OSPRAY
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
#endif//VISIT_OSPRAY


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
#ifdef VISIT_OSPRAY
	regionScaling.x = regionScaling.y = regionScaling.z = 1.0f;
	bgDepthBuffer = NULL;
	bgColorBuffer = NULL;
	initialized = false;
#endif//VISIT_OSPRAY
    }
    ~OSPVisItContext() {	
#ifdef VISIT_OSPRAY
	// clean stuffs
	volumes.clear();
	renderer.Clean();
	camera.Clean();
	transferfcn.Clean();
#endif//VISIT_OSPRAY
    }

#ifdef VISIT_OSPRAY
private:
    friend class OSPVisItVolume;
 public:
    OSPVisItRenderer renderer;
    OSPVisItCamera   camera;
    OSPVisItTransferFunction transferfcn;
    std::vector<OSPVisItVolume> volumes;
private:
    // class parameters
    osp::vec3f regionScaling;
    float         *bgDepthBuffer; // depth buffer for the background and other plots
    unsigned char *bgColorBuffer; // bounding box + pseudo color + ...
    int            bgExtents[4];  // extents of the buffer(minX, maxX, minY, max)
    // ospray mode
    bool initialized;
 public:
    // helper
    void Render(float xMin, float xMax, float yMin, float yMax,
		int imgWidth, int imgHeight, 
		float*& dest, OSPVisItVolume* volume);
    // parameters
    void SetBgBuffer(unsigned char* color, float* depth, int extents[4]) {
	bgColorBuffer = color;
	bgDepthBuffer = depth;
	/* bgExtents[0] = extents[0]; */
	/* bgExtents[1] = extents[1]; */
	/* bgExtents[2] = extents[2]; */
	/* bgExtents[3] = extents[3]; */
    }
    void SetScaling(double s[3]) { 
	regionScaling.x = (float)s[0];
	regionScaling.y = (float)s[1];
	regionScaling.z = (float)s[2]; 
    }

    // patch 
    void InitOSP(int numThreads = 0);
    void InitPatch(int id);
    OSPVisItVolume* GetPatch(int id) { return &volumes[id]; }
#endif//VISIT_OSPRAY
};


// ****************************************************************************
//  Namespace:  slivr
//
//  Purpose:
//    
//
//  Programmer: Qi WU
//  Creation:   
//
// ****************************************************************************
namespace slivr
{
    double deg2rad (double degrees);
    double rad2deg (double radins);
    void CheckMemoryHere(const std::string& message, 
			 std::string debugN = "debug5");
    void CheckMemoryHere(const std::string& message, 
			 std::ostream& out);
};

#endif//AVT_OSPRAY_FILTER_H

