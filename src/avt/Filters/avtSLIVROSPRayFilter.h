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
#include <string>
#include <vector>
#include <cmath>

#include <vtkType.h>
#include <DebugStream.h>

#include "ospray/ospray.h"
#include "ospray/ospcommon/vec.h"
#include "ospray/ospcommon/math.h"
#include "ospray/ospcommon/common.h"

#define OSP_PERSPECTIVE              1
#define OSP_ORTHOGRAPHIC             2
#define OSP_BLOCK_BRICKED_VOLUME     3
#define OSP_SHARED_STRUCTURED_VOLUME 4
#define OSP_INVALID                  5
#define OSP_VALID                    6

typedef ospcommon::vec2f vec2f;
typedef ospcommon::vec2i vec2i;
typedef ospcommon::vec3f vec3f;
typedef ospcommon::vec3i vec3i;
typedef ospcommon::vec4f vec4f;
typedef ospcommon::vec4i vec4i;


namespace slivr {
    // output stream
    extern std::ostream *osp_out;
    extern std::ostream *osp_err;
    // detect environmental variables
    inline bool InitVerbose() 
    {
	auto OSPRAY_VERBOSE_PAIR = 
	    ospcommon::getEnvVar<int>("OSPRAY_VERBOSE");
	if (OSPRAY_VERBOSE_PAIR.first) {
	    slivr::osp_out = &std::cout;
	    slivr::osp_err = &std::cerr;
	    return OSPRAY_VERBOSE_PAIR.second > 0;
	} else {
	    return false;
	}
    }
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
struct VolumeInfo 
{
    // meta info
    int                 patchId = -1;
    bool                isComplete = false;

    // object references (shouldnt be deleted in this struct)
    OSPTransferFunction transferfcn = nullptr;
    OSPRenderer         renderer    = nullptr;

    // objects owned by the struct
    OSPModel            world           = nullptr;
    unsigned char       worldType       = OSP_INVALID;
    OSPFrameBuffer      framebuffer     = nullptr;
    float              *framebufferData = nullptr;
    OSPVolume           volume          = nullptr;
    unsigned char       volumeType      = OSP_INVALID;
    void*               dataPtr         = nullptr;
    std::string         dataType        = "";
    OSPDataType         voxelDataType   = OSP_VOID_PTR;
    OSPData             voxelData       = nullptr;
    size_t              voxelSize       = 0;
    OSPData             ghostData       = nullptr;
    size_t              ghostSize       = 0;

    float               samplingRate = -1.0f;
    vec3f               regionStart;
    vec3f               regionStop;
    vec3f               regionSpacing;
    vec3i               regionSize;
    vec3f               regionUpperClip;
    vec3f               regionLowerClip;

    bool                lightingFlag = false;
    float               specularColor = 0.0f;

    // constructor
    VolumeInfo(int id) : patchId(id) {}
    // destructor
    ~VolumeInfo() { Clean(); }

    void Clean() {
	CleanFBData();
	CleanFB();
	CleanVolume();	
	CleanWorld();
    }
    
    // other function
    void Set
    (void *ptr, int type, unsigned char* ghost,
     double *X, double *Y, double *Z,
     int nX, int nY, int nZ,
     bool cellDataFormat, float sr, 
     double volumePBox[6], double volumeBBox[6],		  
     bool lighting, double mtl[4]);

    bool GetCompleteFlag() { return isComplete; }
    void SetCompleteFlag(bool f) { isComplete = f; } 
    void SetTransferFunction(OSPTransferFunction tf) { transferfcn = tf; }
    void SetRenderer(OSPRenderer r) { renderer = r; }

    // ospModel component
    OSPModel GetWorld() { return world; }
    void InitWorld();
    void SetWorld();
    void CleanWorld() {
	if (world != nullptr) {	    
	    ospRelease(world);
	    world = nullptr;
	}
	worldType = OSP_INVALID;
    }
	
    // ospVolume component
    void InitVolume(unsigned char type = OSP_SHARED_STRUCTURED_VOLUME); 
    OSPVolume GetVolume() { return volume; }
    void SetVolume(void *ptr, int type, unsigned char* ghost,
		   double *X, double *Y, double *Z, 
		   int nX, int nY, int nZ,
		   bool cellDataFormat,
		   double volumePBox[6], 
		   double volumeBBox[6]);
    void SetSamplingRate(float r);
    void SetLighting(bool lighting, float Ks);
    void CleanVolume() {	
	if (voxelData != nullptr) { 
	    ospRelease(voxelData);
	    voxelData = nullptr; 
	}
	if (ghostData != nullptr) {
	    ospRelease(ghostData);
	    ghostData = nullptr; 
	}
	if (volume != nullptr) { ospRelease(volume); volume = nullptr; }
	volumeType = OSP_INVALID;
    }

    // framebuffer component     
    void InitFB(unsigned int width, unsigned int height);
    void RenderFB();
    float* GetFBData();
    void CleanFBData() {
	if (framebufferData != nullptr) { 
	    ospUnmapFrameBuffer(framebufferData, framebuffer); 
	    framebufferData = nullptr;
	}
    }
    void CleanFB() {	   
	if (framebuffer != nullptr) { 
	    ospFreeFrameBuffer(framebuffer); 	    
	    framebuffer = nullptr;
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
struct OSPContext
{
    struct OSPColor
    {
	float R;
	float G;
	float B;
	float A;
    };

    bool refreshData = false;

    std::vector<VolumeInfo> volumePatch;
    OSPRenderer             renderer        = nullptr;
    unsigned char           rendererType    = OSP_INVALID;
    OSPCamera               camera          = nullptr;
    unsigned char           cameraType      = OSP_INVALID;
    OSPTransferFunction     transferfcn     = nullptr;
    unsigned char           transferfcnType = OSP_INVALID;

    float r_panx;
    float r_pany;
    int   screenSize[2];
    bool  enabledOSPRay = false;
    bool  enabledDVR = false; // Distributed Volume Renderer

    // expose this in header
    // because this will be called in other libraries
    ~OSPContext() {	
	// clean stuffs
	volumePatch.clear();
	if (camera      != nullptr) { ospRelease(camera); }
	if (renderer    != nullptr) { ospRelease(renderer); }
	if (transferfcn != nullptr) { ospRelease(transferfcn); }
    }

    // flags
    bool IsEnabled() { return enabledOSPRay; }
    bool IsDVRMode() { return enabledDVR; }
    void setDVRMode(bool mode) { enabledDVR = mode; } 

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
     const double fovy, const double zoomratio, const double imagepan[2], 
     const int imageExtents[4], const int screenExtents[2]);
    void SetSubCamera(float xMin, float xMax, float yMin, float yMax);

    // ospTransferFunction component  
    void InitTransferFunction();
    void SetTransferFunction
    (const OSPColor* table, const unsigned int size, 
     const float datamin, const float datamax);

};

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
    void CheckMemoryHere(const std::string& message, 
			 std::string debugN = "debug5");
    void CheckMemoryHere(const std::string& message, 
			 std::ostream& out);
};

#endif//AVT_OSPRAY_FILTER_H

