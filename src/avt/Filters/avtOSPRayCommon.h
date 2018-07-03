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

#include <string>
#include <vector>
#include <map>

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

typedef ospray::visit::ContextCore OSPVisItContext;
namespace ospray {
    void InitOSP(int numThreads = 0);
    void Finalize();
    struct Context : public ospray::visit::ContextCore {
    public:
	void SetBackgroundBuffer(const unsigned char* color,
				 const float* depth, const int size[2]);
	void SetSpecular(const double& k, const double& n) { Ks = k; Ns = n; }
	void SetScaleAndDataBounds(const double v[3], const double d[6])
	{
	    scale[0] = v[0]; scale[1] = v[1]; scale[2] = v[2];
	    gbbox[0] = d[0] * scale[0];
	    gbbox[3] = d[1] * scale[0];
	    gbbox[1] = d[2] * scale[1];
	    gbbox[4] = d[3] * scale[1];
	    gbbox[2] = d[4] * scale[2];
	    gbbox[5] = d[5] * scale[2];
	}	
	void SetSamplingRate(const double& v) { samplingRate = v; }
	void SetAoSamples(const int v) { aoSamples = v; } 
	void SetSpp(const int v) { spp = v; }
	void SetOneSidedLighting(bool v) { oneSidedLighting = v; }
	void SetShadowsEnabled(bool v) { shadowsEnabled = v; }
	void SetAoTransparencyEnabled(bool v) { aoTransparencyEnabled = v; }
	void SetUseGridAccelerator(bool v) { useGridAccelerator = v; }
	void SetAdaptiveSampling(bool v) { adaptiveSampling = v; }
	void SetPreIntegration(bool v) { preIntegration = v; }
	void SetSingleShade(bool v) { singleShade = v; }
	void SetGradientShadingEnabled(bool v) { gradientShadingEnabled = v; }
	void SetVariableName(const std::string& str) { varname = str; }
	const std::string& GetVariableName() const { return varname; }
	void InitPatch(const int patchID);
	void SetupPatch(const int patchID, const int vtk_type,
			const size_t data_size, const void* data_ptr,
			const double *X, const double *Y, const double *Z, 
			const int nX, const int nY, const int nZ,
			const double dbox[6], const double cbox[6]);
	void RenderPatch(const int patchID,
			 const float xMin, const float xMax, 
			 const float yMin, const float yMax,
			 const int tile_w, const int tile_h,
			 float*& dest); 
    };
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
#include <ospray/visit/VisItWrapper.h>
#include <ospray/visit/VisItModuleCommon.h>
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
    typedef ospray::visit::TransferFunction TransferFunction;
    typedef ospray::visit::Camera Camera;
    typedef ospray::visit::Renderer Renderer;
    typedef ospray::visit::Volume Volume;
    typedef ospray::visit::Model Model;
    typedef ospray::visit::FrameBuffer FrameBuffer;
    typedef ospray::visit::Patch Patch;
    
    void CheckVolumeFormat(const int dt,
			   std::string& str_type,
			   OSPDataType& osp_type);

    void ComputeProjections(const avtViewInfo &view, 
			    const double      &aspect,
			    const double      &old_near_plane,
			    const double      &old_far_plane,
			    const double       data_scale[3],
			    const double       data_bound[6],
			    const int          screen_size[2],
			    vtkMatrix4x4 *model_to_screen_transform, 
			    vtkMatrix4x4 *screen_to_model_transform, 
			    vtkMatrix4x4 *screen_to_camera_transform,
			    double        canvas_size[2],
			    int           rendering_extents[4]);
    
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

