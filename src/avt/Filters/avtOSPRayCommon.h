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

#ifndef VISIT_OSPRAY /* make sure VISIT_OSPRAY is defined */
# error "VISIT_OSPRAY is not defined but ospray is used"
#endif

#ifndef AVT_OSPRAY_COMMON_H
#define AVT_OSPRAY_COMMON_H

#include <ospray/visit/VisItWrapperCore.h>
#include <ospray/visit/VisItWrapper.h>

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
public:
    ospray::visit::PatchCore patch;
private:
    friend class OSPVisItContext;

    OSPVisItContext *parent;
    
    // objects owned by the struct
    // -- ospray model ---
    //OSPModel            world;
    //unsigned char       worldType;
    // --- ospray framebuffer ---
    OSPFrameBuffer      framebuffer;
    float              *framebufferData;
    OSPTexture2D        framebufferBg;
    osp::vec2i imageSize;
    /* // --- ospray volume --- */
    /* OSPVolume           volume; */
    /* unsigned char       volumeType; */
    /* // --- ospray data --- */
    /* OSPDataType         voxelDataType; */
    /* OSPData             voxelData; */
    /* size_t              voxelSize; */
    /* void*               dataPtr; */
    /* std::string         dataType; */

    // metadata for volume
    int                 patchId;       // volume patch id
    bool                finished;      // check if this volume is initialized
    bool                enableShading;
    bool                enableDVR;     // Distributed Volume Renderer
    float               specularKs;
    float               specularNs;
    float               samplingRate;

    /* // geometric parameters for volume */
    /* osp::vec3i          regionSize; */
    /* osp::vec3f          regionStart; */
    /* osp::vec3f          regionStop; */
    /* osp::vec3f          regionSpacing; */
    /* osp::vec3f          regionUpperClip; */
    /* osp::vec3f          regionLowerClip; */
    /* osp::vec3f          regionScaling; */
    

public:
    // constructor
    OSPVisItVolume() {
        // objects owned by the struct
        //world           = NULL;
        //worldType       = OSP_INVALID;
        framebuffer     = NULL;
        framebufferData = NULL;
        framebufferBg   = NULL;
        /* volume          = NULL; */
        /* volumeType      = OSP_INVALID; */
        /* voxelDataType   = OSP_VOID_PTR; */
        /* voxelData       = NULL; */
        /* voxelSize       = 0; */
        /* dataPtr         = NULL; */
        /* dataType        = ""; */
        // metadata for volume
        finished      = false; 
        enableShading = false;
        enableDVR     = false;
        specularKs    = 1.0f;
        specularNs    = 15.0f;
        samplingRate  = 3.0f;
        /* // geometric parameters for volume */
        /* regionSize.x  = regionSize.y  = regionSize.z  = 0; */
        /* regionStart.x = regionStart.y = regionStart.z = 0.0f; */
        /* regionStop.x  = regionStop.y  = regionStop.z  = 0.0f; */
        /* regionSpacing.x   = regionSpacing.y   = regionSpacing.z   = 0.0f; */
        /* regionUpperClip.x = regionUpperClip.y = regionUpperClip.z = 0.0f; */
        /* regionLowerClip.x = regionLowerClip.y = regionLowerClip.z = 0.0f; */
        /* regionScaling.x   = regionScaling.y   = regionScaling.z   = 1.0f; */
    }

    // destructor
    ~OSPVisItVolume() { Clean(); }    
    void Clean() {
        CleanFB();
        //CleanVolume();	
        //CleanWorld();
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
    /* OSPModel GetWorld() { return world; } */
    /* void InitWorld(); */
    /* void SetWorld(); */
    /* void CleanWorld() { */
    /*     if (world != NULL) {	     */
    /*         ospRelease(world); */
    /*         world = NULL; */
    /*     } */
    /*     worldType = OSP_INVALID; */
    /* } */
	
    /* // ospVolume component */
    /* void InitVolume(int type, void *ptr, */
    /* 		    int nX, int nY, int nZ, */
    /* 		    unsigned char volumeType);  */
    /* OSPVolume GetVolume() { return volume; } */
    /* void SetVolume(double *X, double *Y, double *Z,  */
    /*                int nX, int nY, int nZ, */
    /*                double volumePBox[6],  */
    /*                double volumeBBox[6]); */
    /* void CleanVolume() {	 */
    /*     if (volume != NULL) { ospRelease(volume); volume = NULL; } */
    /*     if (voxelData != NULL) {  */
    /*         ospRelease(voxelData); */
    /*         voxelData = NULL;  */
    /*     } */
    /*     volumeType = OSP_INVALID; */
    /* } */

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
    void SetScaleAndDataBounds(double s[3], double d[6]) {
	
	regionScaling.x = (float)s[0];
        regionScaling.y = (float)s[1];
        regionScaling.z = (float)s[2];
	
	bbox.lower.x = d[0] * regionScaling.x;
	bbox.upper.x = d[1] * regionScaling.x;
	bbox.lower.y = d[2] * regionScaling.y;
	bbox.upper.y = d[3] * regionScaling.y;
	bbox.lower.z = d[4] * regionScaling.z;
	bbox.upper.z = d[5] * regionScaling.z;

        for (int i = 0; i < 6; ++i) { bounds[i] = d[i]; }
    }
    

    ospray::visit::CameraCore           camera;
    ospray::visit::RendererCore         renderer;
    ospray::visit::TransferFunctionCore tfn;


    void SetBgBuffer(unsigned char* color, float* depth, int size[2]) 
    {
        ((ospray::visit::Renderer)renderer)
            .SetBackgroundBuffer(color, depth, size);
    }


    std::map<int, OSPVisItVolume> volumes;

    void SetActiveVariable(const char* str) { varname = str; }
    const std::string& GetActiveVariable() const { return varname; }
    
private:
    
    friend class OSPVisItVolume;

    osp::vec3f     regionScaling;
    
    osp::box3f     bbox;
    
    double bounds[6];
    static bool initialized;
    std::string varname;
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

#include <avtParallel.h>
#include <avtViewInfo.h>

#include <DebugStream.h>
#include <StackTimer.h>
#include <TimingsManager.h>
#include <ImproperUseException.h>

#include <vtkType.h>
#include <vtkMatrix4x4.h>

#include <ospray/ospray.h>
#include <ospray/visit/VisItModuleCommon.h>
#include <ospray/visit/VisItWrapper.h>
#include <ospray/visit/VisItImageComposite.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdio.h>

#ifndef CLAMP
# define CLAMP(x, l, h) (x > l ? x < h ? x : h : l)
#endif
#ifndef M_MIN
# define M_MIN(x, r) (x < r ? x : r)
#endif
#ifndef M_MAX
# define M_MAX(x, r) (x > r ? x : r)
#endif

// ostreams customized for ospray
#ifdef ospout
#undef ospout
#endif
#define ospout                                                      \
    if (!ospray::visit::CheckVerbose() && !DebugStream::Level5()) ; \
    else (*ospray::osp_out)
#ifdef osperr
#undef osperr
#endif
#define osperr                                                      \
    if (!ospray::visit::CheckVerbose() && !DebugStream::Level1()) ; \
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
};

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
    
    void CheckVolumeFormat(const int dt,
			   std::string& str_type,
			   OSPDataType& osp_type);

    void ComputeProjections(const avtViewInfo &view, 
			    const double &aspect,
			    const int screen[2],
			    const double scale[3],
			    const double &oldNearPlane,
			    const double &oldFarPlane,
			    vtkMatrix4x4  *model_to_screen_transform, 
			    vtkMatrix4x4  *screen_to_model_transform, 
			    vtkMatrix4x4  *screen_to_camera_transform,
			    int            renderingExtents[4],
			    double         sceneSize[2],
			    double         dbounds[6]);
    
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

