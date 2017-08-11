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
	const char* env_verbose = std::getenv("OSPRAY_VERBOSE");
	if (env_verbose) {
	    if (atoi(env_verbose) > 0) {
		slivr::osp_out = &std::cout;
		slivr::osp_err = &std::cerr;
		return true;
	    }
	    else {
		return false;
	    }
	}
	else {
	    return false;
	}
#endif
    }
    //! this function has to be inline, otherwise we need to 
    //! modify library linkages
    inline bool CheckVerbose() // initialize OSPRAY_VERBOSE
    {
	static bool OSPRAY_VERBOSE = slivr::InitVerbose();
	return OSPRAY_VERBOSE;
    }
};
#define ospout \
    if (!slivr::CheckVerbose() && !DebugStream::Level5()) ; \
    else (*slivr::osp_out)
#define osperr \
    if (!slivr::CheckVerbose() && !DebugStream::Level1()) ; \
    else (*slivr::osp_err)

// ****************************************************************************
//  Struct:  VolumeInfo
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
class VolumeInfo 
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
    bool                isComplete; // check if this volume is initialized
    bool                lightingFlag;
    float               specularColor;
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
    VolumeInfo(int id) {
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
	isComplete = false; 
	lightingFlag = false;
	specularColor = 0.0f;
	specularNs    = 0.0f;
	samplingRate = -1.0f;
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
    ~VolumeInfo() { Clean(); }    
    void Clean() {
	CleanFBData();
	CleanFB();
	CleanVolume();	
	CleanWorld();
    }
    
    // other function
    void Set(int type, void *ptr, unsigned char* ghost,
	     double *X, double *Y, double *Z, int nX, int nY, int nZ,
	     double volumePBox[6], double volumeBBox[6], double mtl[4],
	     float sr, bool lighting, bool cellDataFormat);
    bool GetCompleteFlag() { return isComplete; }
    void SetCompleteFlag(bool f) { isComplete = f; } 
    void SetTransferFunction(const OSPTransferFunction& tf) { transferfcn = tf; }
    void SetRenderer(const OSPRenderer& r) { renderer = r; }
    void SetScaling(const osp::vec3f& s) { regionScaling = s; }

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
    void SetVolume(int type, void *ptr, unsigned char* ghost,
		   double *X, double *Y, double *Z, 
		   int nX, int nY, int nZ,
		   double volumePBox[6], 
		   double volumeBBox[6],
		   bool cellDataFormat);
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
    std::vector<VolumeInfo> volumePatch;
    OSPLight aLight;
    OSPLight dLight;
    OSPData  lightdata;
    OSPRenderer             renderer;
    unsigned char           rendererType;
    OSPCamera               camera;
    unsigned char           cameraType;
    OSPTransferFunction     transferfcn;
    unsigned char           transferfcnType;
    // class parameters
    bool refreshData;
    osp::vec3f regionScaling;
    float r_panx;
    float r_pany;
    float zoom;
    int   screenSize[2];
    bool  enabledOSPRay;
    bool  enabledDVR; // (not used yet) Distributed Volume Renderer
 public:
    OSPContext() {
	// ospray objects
	aLight = NULL;
	dLight = NULL;
	lightdata = NULL;
	renderer        = NULL;
	rendererType    = OSP_INVALID;
	camera          = NULL;
	cameraType      = OSP_INVALID;
	transferfcn     = NULL;
	transferfcnType = OSP_INVALID;
	// class parameters
	refreshData = false;
	regionScaling.x = regionScaling.y = regionScaling.z = 1.0f;
	r_panx = 0.0f;
	r_pany = 0.0f;
	zoom = 1.0f;
	screenSize[0] = screenSize[1] = 0.0f;
	enabledOSPRay = false;
	enabledDVR = false; // Distributed Volume Renderer
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

    // flags
    bool IsEnabled() { return enabledOSPRay; }
    bool IsDVRMode() { return enabledDVR; }
    void SetDVRMode(bool mode) { enabledDVR = mode; } 
    void SetScaling(double s[3]) { 
	regionScaling.x = (float)s[0];
	regionScaling.y = (float)s[1];
	regionScaling.z = (float)s[2]; 
    }

    // patch 
    void InitOSP(bool flag, int numThreads = 0);
    void InitPatch(int id);
    VolumeInfo* GetPatch(int id) { return &volumePatch[id]; }

    // ospRenderer component
    void InitRenderer();
    void SetRenderer(bool lighting, double mtl[4], double dir[3]);
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

