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

namespace slivr {
    // output stream
    extern std::ostream *osp_out;
    extern std::ostream *osp_err;
    //! this function has to be inline, otherwise we need to 
    //! modify library linkages
    // detect environmental variables
    inline bool InitVerbose() 
    {
#ifndef VISIT_OSPRAY
	return false;
#else
	/* OSPRay defined environmental variables */
	// OSPRAY_THREADS, OSPRAY_SET_AFFINITY
	// OSPRAY_DEBUG, OSPRAY_LOG_LEVEL
	// OSPRAY_LOG_OUTPUT
	/* visit shortcuts*/
	// OSPRAY_VERBOSE
	const char* env_debug = std::getenv("OSPRAY_DEBUG");
	const char* env_verbose = std::getenv("OSPRAY_VERBOSE");
	const char* env_log_level = std::getenv("OSPRAY_LOG_LEVEL");	
	bool verbose = false;
	if (env_debug) {
	    if (atoi(env_debug) > 0) { verbose = true; }
	}
	if (env_verbose) {
	    if (atoi(env_verbose) > 0) { verbose = true; }
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
    //! this function has to be inline, otherwise we need to 
    //! modify library linkages
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
#define ospout \
    if (!slivr::CheckVerbose() && !DebugStream::Level5()) ; \
    else (*slivr::osp_out)
#define osperr \
    if (!slivr::CheckVerbose() && !DebugStream::Level1()) ; \
    else (*slivr::osp_err)

// ****************************************************************************
//  Struct:  OSPVolumePatch
//
//  Purpose:
//    Holds information about patches but not the image 
//
//  Programmer:  
//  Creation:   
//
// ****************************************************************************
// this struct will contain informations related to one volume
#ifdef VISIT_OSPRAY
class OSPVolumePatch 
{    
 private:
    // object references 
    // (those objects are created and handled by other parts of the program)
    // (shouldnt be deleted in this struct)
    OSPTransferFunction transferfcn;
    OSPRenderer         renderer;
    
    // objects owned by the struct
    // -- ospray model ---
    OSPModel            world;
    unsigned char       worldType;
    // --- ospray framebuffer ---
    OSPFrameBuffer      framebuffer;
    float              *framebufferData;
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
    OSPVolumePatch(int id) {
	// object references 
	transferfcn = NULL;
	renderer    = NULL;
	// objects owned by the struct
	world           = NULL;
	worldType       = OSP_INVALID;
	framebuffer     = NULL;
	framebufferData = NULL;    
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
    ~OSPVolumePatch() { Clean(); }    
    void Clean() {
	CleanFBData();
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
    void SetScaling(const osp::vec3f& s) { regionScaling = s; }
    void SetTransferFunction(const OSPTransferFunction& t) { transferfcn = t; }
    void SetRenderer(const OSPRenderer& r) { renderer = r; }

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
    void CleanFBData() {
	if (framebufferData != NULL) { 
	    ospUnmapFrameBuffer(framebufferData, framebuffer); 
	    framebufferData = NULL;
	}
    }
    void CleanFB() {	   
	if (framebuffer != NULL) { 
	    ospFreeFrameBuffer(framebuffer); 	    
	    framebuffer = NULL;
	}	
    }

};

// ****************************************************************************
//  Struct:  OSPContext
//
//  Purpose:
//    Holds information about patches but not the image 
//
//  Programmer:  
//  Creation:   
//
// ****************************************************************************
class OSPContext
{
 public:
    struct OSPColor
    {
	float R;
	float G;
	float B;
	float A;
    };
 private:
    // ospray objects
    std::vector<OSPVolumePatch> volumePatch;
    OSPLight aLight;
    OSPLight dLight;
    OSPLight sLight; // constant sun light
    OSPData  lightdata;
    OSPRenderer             renderer;
    unsigned char           rendererType;
    OSPCamera               camera;
    unsigned char           cameraType;
    OSPTransferFunction     transferfcn;
    unsigned char           transferfcnType;
    // -- DVR mode variable --
    OSPModel                modelDVR;
    OSPFrameBuffer          fbDVR;
    float                  *fbDataDVR;
    std::vector<osp::box3f> volumeRegionsDVR;
    // class parameters
    osp::vec3f regionScaling;
    float r_panx;
    float r_pany;
    float zoom;
    int   screenSize[2];
    // ospray mode
    bool  refreshData;
    bool  enableDVR; // (not used yet) Distributed Volume Renderer
 public:
    OSPContext() {
	// ospray objects
	aLight = NULL;
	dLight = NULL;
	sLight = NULL;
	lightdata = NULL;
	renderer        = NULL;
	rendererType    = OSP_INVALID;
	camera          = NULL;
	cameraType      = OSP_INVALID;
	transferfcn     = NULL;
	transferfcnType = OSP_INVALID;
	//--- DVR ---
	modelDVR        = NULL;
	fbDVR           = NULL;
	fbDataDVR       = NULL;
	// class parameters
	regionScaling.x = regionScaling.y = regionScaling.z = 1.0f;
	r_panx = 0.0f;
	r_pany = 0.0f;
	zoom   = 1.0f;
	screenSize[0] = screenSize[1] = 0.0f;
	// ospray mode
	refreshData = false;
	enableDVR   = false; // Distributed Volume Renderer
    }
    // expose this in header
    // because this will be called in other libraries
    ~OSPContext() {	
	// clean stuffs
	volumePatch.clear();
	if (renderer    != NULL) {
	    ospRelease(lightdata);
	    ospRelease(renderer); 
	}
	if (camera      != NULL) { ospRelease(camera); }
	if (transferfcn != NULL) { ospRelease(transferfcn); }
	rendererType    = OSP_INVALID;
	cameraType      = OSP_INVALID;
	transferfcnType = OSP_INVALID;
    }

    // helper
    void Render(float xMin, float xMax, float yMin, float yMax,
		int imgWidth, int imgHeight, 
		float*& dest, OSPVolumePatch* volume);
    // flags
    bool IsDVRMode() { return enableDVR; }
    void SetDVRMode(bool mode) { enableDVR = mode; } 

    // parameters
    void SetScaling(double s[3]) { 
	regionScaling.x = (float)s[0];
	regionScaling.y = (float)s[1];
	regionScaling.z = (float)s[2]; 
    }

    // patch 
    void InitOSP(bool flag, int numThreads = 0);
    void InitPatch(int id);
    OSPVolumePatch* GetPatch(int id) { return &volumePatch[id]; }

    // ospRenderer component
    void InitRenderer();
    void SetRenderer(bool shading, double mtl[4], double dir[3]);
    void SetModel(OSPModel world);

    // ospCamera component
    void InitCamera(unsigned char type);
    void SetCamera
    (const double campos[3], const double camfocus[3], const double camup [3], 
     const double camdir[3], const double sceneSize[2], const double aspect,
     const double viewAngle, const double zoomratio, const double imagepan[2], 
     const int imageExtents[4], const int screenExtents[2]);
    void SetSubCamera(float xMin, float xMax, float yMin, float yMax);

    // ospTransferFunction component  
    void InitTransferFunction();
    void SetTransferFunction
    (const OSPColor* table, const unsigned int size, 
     const float datamin, const float datamax);
};
#endif//VISIT_OSPRAY

// ****************************************************************************
//  Namespace:  slivr
//
//  Purpose:
//    
//
//  Programmer:  
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

