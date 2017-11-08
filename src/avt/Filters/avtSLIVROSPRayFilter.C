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

#include "avtSLIVROSPRayFilter.h"

#include <avtMemory.h>
#include <avtParallel.h>
#include <ImproperUseException.h>
#include <TimingsManager.h>

#ifdef __unix__
# include <unistd.h>
#endif

// helper
namespace slivr {
    // output stream
    std::ostream *osp_out = &DebugStream::Stream5();
    std::ostream *osp_err = &DebugStream::Stream1();
};

double slivr::deg2rad (double degrees) {
    return degrees * 4.0 * atan (1.0) / 180.0;
}
double slivr::rad2deg (double radins) {
    return radins / 4.0 / atan (1.0) * 180.0;
}

// other function
void 
OSPVolumePatch::Set
(int type, void *ptr, double *X, double *Y, double *Z, 
 int nX, int nY, int nZ, double volumePBox[6], double volumeBBox[6], 
 double mtl[4], float sr, bool shading)
{
    /* OSPRay Volume */
    specularKs    = (float)mtl[2];
    specularNs    = (float)mtl[3];
    enableShading = shading;
    samplingRate  = sr;
    // TODO: It seems if a volume is recovered from a session
    // ospray will crash during zooming ...
    // So we refresh volume everytime to fix the bug
    // which means we need to disable grid accelerator
    // to speed things up. Until I found the reason of crashing
    if (ptr != dataPtr) {
       ospout << "[ospray] update data" << std::endl;
    };
    if (true /*!finished*/) {
	// Because we initialized the volume each frame
	// we need to removed the old volume from model first
	// if (finished) {
	//     ospRemoveVolume(world, volume);
	//     ospCommit(world);	    
	// }
	volumeType = OSP_INVALID;
	InitVolume();
	SetVolume(type, ptr, X, Y, Z, nX, nY, nZ,
		  volumePBox, volumeBBox);
	// now we add the volume back
	// if (finished) { SetWorld(); }
    }
    /* OSPRay Model */
    if (true/*!finished*/) {
	worldType = OSP_INVALID; 
	InitWorld();
	SetWorld();
    }
    /* update volume */
    finished = true;
}

// ospModel component
void OSPVolumePatch::InitWorld() {
    if (worldType == OSP_INVALID) {
	CleanWorld();
	worldType = OSP_VALID;
	world = ospNewModel();
    }
}
void OSPVolumePatch::SetWorld() {
    if (world != NULL) { 
	ospAddVolume(world, volume);
	ospCommit(world);
    }
}

// ospVolume component
void OSPVolumePatch::InitVolume(unsigned char type) {
    if (volumeType != type) { // only initialize once
	CleanVolume();
	volumeType = type;
	switch (type) {
	case (OSP_BLOCK_BRICKED_VOLUME):
	    volume = ospNewVolume("block_bricked_volume"); 
	    break;
	case (OSP_SHARED_STRUCTURED_VOLUME):
	    volume = ospNewVolume("visit_shared_structured_volume"); 
	    break;
	default:
	    debug1 << "ERROR: ospray volume not initialized"
		   << std::endl;
	    volumeType = OSP_INVALID;
	    EXCEPTION1(VisItException, 
		       "ERROR: ospray volume not initialized");
	}
    }
}
void 
OSPVolumePatch::SetVolume
(int type, void *ptr, double *X, double *Y, double *Z, 
 int nX, int nY, int nZ, double volumePBox[6], double volumeBBox[6]) 
{
    // calculate volume data type
    if (type == VTK_UNSIGNED_CHAR) {
	dataType = "uchar";
	voxelDataType = OSP_UCHAR;
    } else if (type == VTK_SHORT) {
	dataType = "short";
	voxelDataType = OSP_SHORT;
    } else if (type == VTK_UNSIGNED_SHORT) {
	dataType = "ushort";
	voxelDataType = OSP_USHORT;
    } else if (type == VTK_FLOAT) {
	dataType = "float";
	voxelDataType = OSP_FLOAT;
    } else if (type == VTK_DOUBLE) {
	dataType = "double";
	voxelDataType = OSP_DOUBLE;
    } else {
	debug1 << "ERROR: Unsupported ospray volume type" << std::endl;
	EXCEPTION1(VisItException, "ERROR: Unsupported ospray volume type");
    }
    ospout << "[ospray] data type " << dataType << std::endl;
    // assign data pointer
    dataPtr = ptr;
    // assign structure
    regionStart.x   = volumePBox[0];
    regionStart.y   = volumePBox[1];
    regionStart.z   = volumePBox[2];
    regionStop.x    = volumePBox[3];
    regionStop.y    = volumePBox[4];
    regionStop.z    = volumePBox[5];
    regionSize.x    = nX;
    regionSize.y    = nY;
    regionSize.z    = nZ;
    regionSpacing.x = (regionStop.x-regionStart.x)/((float)regionSize.x-1.0f);
    regionSpacing.y = (regionStop.y-regionStart.y)/((float)regionSize.y-1.0f);
    regionSpacing.z = (regionStop.z-regionStart.z)/((float)regionSize.z-1.0f);
    regionLowerClip.x = volumeBBox[0];
    regionLowerClip.y = volumeBBox[1];
    regionLowerClip.z = volumeBBox[2];
    regionUpperClip.x = volumeBBox[3];
    regionUpperClip.y = volumeBBox[4];
    regionUpperClip.z = volumeBBox[5];

    // other objects
    ospSetString(volume, "voxelType", dataType.c_str());
    ospSetObject(volume, "transferFunction", transferfcn);

    // commit voxel data
    if (voxelData != NULL) { 
	debug1 << "ERROR: Found VoxelData to be non-empty "
	       << "while creating new volume" << std::endl;
	EXCEPTION1(VisItException, 
		   "ERROR: Found VoxelData to be non-empty "
		   "while creating new volume");
    }
    voxelSize = nX * nY * nZ;
    voxelData = ospNewData(voxelSize, voxelDataType,
			   dataPtr, OSP_DATA_SHARED_BUFFER);
    ospSetData(volume, "voxelData", voxelData);

    // commit volume
    // -- no lighting by default
    ospout << "[ospray] setting specular value to " << specularKs << std::endl;
    osp::vec3f Ks; Ks.x = Ks.y = Ks.z = specularKs;
    ospSetVec3f(volume, "specular", Ks);
    ospSet1f(volume, "Ns", specularNs);
    ospSet1i(volume, "gradientShadingEnabled", (int)enableShading);
    // -- other properties
    osp::vec3f scaledBBoxLower;
    osp::vec3f scaledBBoxUpper;
    osp::vec3f scaledSpacing;
    osp::vec3f scaledOrigin;

    scaledBBoxLower.x = regionLowerClip.x * regionScaling.x;
    scaledBBoxUpper.x = regionUpperClip.x * regionScaling.x;
    scaledSpacing.x = regionSpacing.x * regionScaling.x;
    scaledOrigin.x  = regionStart.x * regionScaling.x;

    scaledBBoxLower.y = regionLowerClip.y * regionScaling.y;
    scaledBBoxUpper.y = regionUpperClip.y * regionScaling.y;
    scaledSpacing.y = regionSpacing.y * regionScaling.y;
    scaledOrigin.y  = regionStart.y * regionScaling.y;

    scaledBBoxLower.z = regionLowerClip.z * regionScaling.z;
    scaledBBoxUpper.z = regionUpperClip.z * regionScaling.z;
    scaledSpacing.z = regionSpacing.z * regionScaling.z;
    scaledOrigin.z  = regionStart.z * regionScaling.z;

    ospSet1i(volume, "useGridAccelerator", 0);
    ospSetVec3f(volume, "volumeClippingBoxLower", scaledBBoxLower);
    ospSetVec3f(volume, "volumeClippingBoxUpper", scaledBBoxUpper);
    ospSetVec3f(volume, "gridSpacing", scaledSpacing);
    ospSetVec3f(volume, "gridOrigin",  scaledOrigin);
    ospSetVec3i(volume, "dimensions",  regionSize);
    ospSet1f(volume, "samplingRate", samplingRate);
    ospSet1i(volume, "adaptiveSampling", 0);
    ospSet1i(volume, "preIntegration", 0);
    ospSet1i(volume, "singleShade", 0);
    ospCommit(volume);
}

// framebuffer component     
void OSPVolumePatch::InitFB(unsigned int width, unsigned int height) {
    osp::vec2i imageSize;
    imageSize.x = width;
    imageSize.y = height;
    CleanFBData(); CleanFB();	    
    framebuffer = ospNewFrameBuffer(imageSize, 
				    OSP_FB_RGBA32F,
				    OSP_FB_COLOR);	    
}
void OSPVolumePatch::RenderFB() {
    // ospFrameBufferClear(framebuffer, OSP_FB_COLOR);
    ospRenderFrame(framebuffer, renderer, OSP_FB_COLOR);
    framebufferData = (float*) ospMapFrameBuffer(framebuffer, OSP_FB_COLOR);
}
float* OSPVolumePatch::GetFBData() {
    return framebufferData;
}

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
void OSPContext_ErrorFunc(OSPError, const char* msg) { osperr << msg; }
void OSPContext_StatusFunc(const char* msg) { ospout << msg; }
void OSPContext::InitOSP(bool flag, int numThreads) 
{ 
    OSPDevice device = ospGetCurrentDevice();
    if (device == NULL) 
    {
	// check hostname
#ifdef __unix__
	char hname[200];
	gethostname(hname, 200);
        ospout << "[ospray] on host >> " << hname << "<<" << std::endl;;
#endif
	// load ospray module
	if (enableDVR) {
	    OSPError err = ospLoadModule("mpi");
	    if (err != OSP_NO_ERROR) {
		osperr << "[Error] can't load ospray MPI module" << std::endl;
	    }
	}
	OSPError err = ospLoadModule("visit");
	if (err != OSP_NO_ERROR) {
	    osperr << "[Error] can't load visit module" << std::endl;
	}
	// initialize ospray
        ospout << "[ospray] Initialize OSPRay";
	if (enableDVR) {
	    device = ospNewDevice("mpi_distributed");
	    ospDeviceSet1i(device, "masterRank", 0);
	} else {
	    device = ospNewDevice();
	}
	// setup debug 
	if (DebugStream::Level5()) {
	    ospout << " debug mode";
	    ospDeviceSet1i(device, "debug", 0);
	}
	// setup number of threads (this can only be hard-coded)
	if (numThreads > 0) {
	    ospout << " numThreads: " << numThreads;
	    ospDeviceSet1i(device, "numThreads", numThreads);
	}
	ospout << std::endl;
	ospDeviceSetErrorFunc(device, OSPContext_ErrorFunc);
	ospDeviceSetStatusFunc(device, OSPContext_StatusFunc);
	ospDeviceCommit(device);
	ospSetCurrentDevice(device);
    }
    refreshData = flag;
    ospout << "[ospray] Initialize OSPRay (new data " << flag << ")" 
	   << std::endl;    
}

// We use this function to minimize interface
void OSPContext::Render
(float xMin, float xMax, float yMin, float yMax,
 int imgWidth, int imgHeight, float*& dest, OSPVolumePatch* volume) 
{
    // render frame
    int timing_setsubcamera = visitTimer->StartTimer();
    SetSubCamera(xMin, xMax, yMin, yMax);
    visitTimer->StopTimer(timing_setsubcamera, "[OSPRay] Calling OSPContext::SetSubCamera");

    int timing_setmodel = visitTimer->StartTimer();
    SetModel(volume->GetWorld());
    visitTimer->StopTimer(timing_setmodel, "[OSPRay] Calling OSPContext::SetModel");

    int timing_initfb = visitTimer->StartTimer();
    volume->InitFB(imgWidth, imgHeight);
    visitTimer->StopTimer(timing_initfb, "[OSPRay] Calling OSPContext::InitFB");

    int timing_renderfb = visitTimer->StartTimer();
    volume->RenderFB();
    visitTimer->StopTimer(timing_renderfb, "[OSPRay] Calling OSPContext::RenderFB");

    // copy data
    int timing_stdcopy = visitTimer->StartTimer();
    std::copy(volume->GetFBData(), 
	      volume->GetFBData() + (imgWidth * imgHeight) * 4, 
	      dest);
    visitTimer->StopTimer(timing_stdcopy, "[OSPRay] Calling OSPContext::std::copy");

    int timing_cleanfbdata = visitTimer->StartTimer();
    volume->CleanFBData();
    visitTimer->StopTimer(timing_cleanfbdata, "[OSPRay] Calling OSPContext::CleanFBData");
}

void OSPContext::InitPatch(int id) 
{
    if (volumePatch.size() < id) {
	debug1 << "ERROR: wrong patch index " << id << std::endl;
	EXCEPTION1(VisItException, "ERROR: wrong patch index"); 
	return;
    }
    if (volumePatch.size() == id) { 
	volumePatch.push_back(id); 
    }
    volumePatch[id].SetScaling(regionScaling);
    volumePatch[id].SetTransferFunction(transferfcn);
    volumePatch[id].SetRenderer(renderer);
    volumePatch[id].SetDVRFlag(enableDVR);
    volumePatch[id].SetFinishedFlag(!refreshData); // reset volume for new data
    // if the data is refreshed -> not complete
}

// ospRenderer component
void OSPContext::InitRenderer() 
{
    if (rendererType == OSP_INVALID) {
	renderer = ospNewRenderer("scivis");
	aLight = ospNewLight(renderer, "ambient");
	dLight = ospNewLight(renderer, "distant");
	sLight = ospNewLight(renderer, "distant");
	OSPLight lights[3] = { aLight, dLight, sLight };
	lightdata = ospNewData(3, OSP_OBJECT, lights);
	rendererType = OSP_VALID;
    }
}

void OSPContext::SetRenderer(bool shading, double mtl[4], double dir[3]) 
{
    ospSetObject(renderer, "camera", camera);
    ospSet1i(renderer, "backgroundEnabled", 0);
    ospSet1i(renderer, "oneSidedLighting", 0);
    ospSet1i(renderer, "aoSamples", 0);
    ospSet1i(renderer, "spp", slivr::CheckOSPRaySpp());
    if (shading)
    {
	ospout << "[ospray] use lighting " 
	       << "use mtl " 
	       << mtl[0] << " "
	       << mtl[1] << " "
	       << mtl[2] << " "
	       << mtl[3] << std::endl;
	ospSet1i(renderer, "shadowsEnabled", 0);
	// light direction
	osp::vec3f lightDir;
	lightDir.x = (float)dir[0];
	lightDir.y = (float)dir[1];
	lightDir.z = (float)dir[2];
	// ambient light
	ospSet1f(aLight, "intensity", (float)mtl[0]);
	ospSet1i(aLight, "isVisible", 0);
	ospCommit(aLight);
	// directional light
	ospSet1f(dLight, "intensity", (float)mtl[1]);
	ospSet1f(dLight, "angularDiameter", 0.53f);
	ospSet1i(dLight, "isVisible", 0);
	ospSetVec3f(dLight, "direction", lightDir);
	ospCommit(dLight);
	// directional light
	ospSet1f(sLight, "intensity", 1.5f);
	ospSet1f(sLight, "angularDiameter", 0.53f);
	ospSet1i(sLight, "isVisible", 0);
	ospSetVec3f(sLight, "direction", lightDir);
	ospCommit(sLight);
	ospSetData(renderer, "lights", lightdata);
    }
    ospCommit(renderer);
}

void OSPContext::SetModel(OSPModel world)
{
    ospSetObject(renderer, "camera", camera);
    ospSetObject(renderer, "model",  world);
    ospCommit(renderer);
}

// ospCamera component
void OSPContext::InitCamera(unsigned char type) 
{
    if (cameraType != type) {
	cameraType = type;
	if (camera != NULL) { ospRelease(camera); }
	switch (type) {
	case (OSP_PERSPECTIVE):
	    camera = ospNewCamera("perspective");
	    break;
	case (OSP_ORTHOGRAPHIC):
	    camera = ospNewCamera("orthographic");
	    break;
	default:
	    debug1 << "ERROR: wrong ospray camera type" << std::endl;
	    cameraType = OSP_INVALID;
	    EXCEPTION1(VisItException, "ERROR: wrong ospray camera type"); 
	}
    }
}

void OSPContext::SetCamera(const double campos[3], 
			   const double camfocus[3], 
			   const double camup [3], 
			   const double camdir[3],
			   const double sceneSize[2],
			   const double aspect, 
			   const double viewAngle, 
			   const double zoomratio, 
			   const double imagepan[2],
			   const int imageExtents[4],
			   const int screenExtents[2]) 
{
    osp::vec3f camPos, camDir, camUp;
    camPos.x = campos[0]; camPos.y = campos[1]; camPos.z = campos[2];    
    camDir.x = camdir[0]; camDir.y = camdir[1]; camDir.z = camdir[2];
    camUp.x = camup[0];   camUp.y = camup[1];   camUp.z = camup[2];
    ospSetVec3f(camera, "pos", camPos);
    ospSetVec3f(camera, "dir", camDir);
    ospSetVec3f(camera, "up",  camUp);
    if (cameraType == OSP_PERSPECTIVE) {
	ospSet1f(camera, "aspect", aspect);
	ospSet1f(camera, "fovy", viewAngle);
    }
    else if (cameraType == OSP_ORTHOGRAPHIC) {
	ospSet1f(camera, "aspect", aspect);
	ospSet1f(camera, "height", sceneSize[1]);
    }
    r_panx = imagepan[0] * zoomratio;
    r_pany = imagepan[1] * zoomratio;
    this->SetSubCamera(imageExtents[0], imageExtents[1],
		       imageExtents[2], imageExtents[3]);
    screenSize[0] = screenExtents[0];
    screenSize[1] = screenExtents[1];
    zoom = zoomratio;
}

void OSPContext::SetSubCamera(float xMin, float xMax, float yMin, float yMax) 
{
    osp::vec2f imgS, imgE;
    float r_xl = xMin/screenSize[0] - r_panx; 
    float r_yl = yMin/screenSize[1] - r_pany; 
    float r_xu = xMax/screenSize[0] - r_panx;
    float r_yu = yMax/screenSize[1] - r_pany;	
    imgS.x = (r_xl - 0.5f) / zoom + 0.5f;
    imgS.y = (r_yl - 0.5f) / zoom + 0.5f;
    imgE.x = (r_xu - 0.5f) / zoom + 0.5f;
    imgE.y = (r_yu - 0.5f) / zoom + 0.5f;
    ospSetVec2f(camera, "imageStart", imgS);
    ospSetVec2f(camera, "imageEnd",   imgE);
    ospCommit(camera);
}

// ospTransferFunction component  
void OSPContext::InitTransferFunction() 
{
    if (transferfcnType == OSP_INVALID) {
	if (transferfcn != NULL) { ospRelease(transferfcn); }
	transferfcn = ospNewTransferFunction("piecewise_linear");
	transferfcnType = OSP_VALID;
    }
}

void OSPContext::SetTransferFunction(const OSPColor *table,
				     const unsigned int size, 
				     const float datamin, 
				     const float datamax) 
{
    std::vector<osp::vec3f> colors;
    std::vector<float>      opacities;
    for (int i = 0; i < size; ++i) {
	osp::vec3f color;
	color.x = table[i].R;
	color.y = table[i].G;
	color.z = table[i].B;
	colors.push_back(color);
	opacities.push_back(table[i].A);
    }

    // printf("\n\n");
    // for (int i = 0; i < size; ++i) {
    // 	printf("vec3f(%f,%f,%f),\n",colors[i].x,colors[i].y,colors[i].z);
    // }
    // printf("\n{");
    // for (int i = 0; i < size; ++i) {
    // 	printf("%f,",opacities[i]);
    // }
    // printf("}\n");

    OSPData colorData   = 
	ospNewData(colors.size(), OSP_FLOAT3, colors.data());
    OSPData opacityData = 
	ospNewData(opacities.size(), OSP_FLOAT, opacities.data());
    osp::vec2f range;
    range.x = datamin;
    range.y = datamax;
    ospSetData(transferfcn, "colors",      colorData);
    ospSetData(transferfcn, "opacities",   opacityData);
    ospSetVec2f(transferfcn, "valueRange", range);
    ospCommit(transferfcn);
    ospRelease(colorData);
    ospRelease(opacityData);
}

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

void slivr::CheckMemoryHere(const std::string& message, std::string debugN)
{
    if (debugN.compare("debug5") == 0) {
	if (DebugStream::Level5()) {
	    slivr::CheckMemoryHere(message, DebugStream::Stream5());
	}       
    }
    else if (debugN.compare("debug4") == 0) {
	if (DebugStream::Level4()) {
	    slivr::CheckMemoryHere(message, DebugStream::Stream4());
	}       
    }
    else if (debugN.compare("debug3") == 0) {
	if (DebugStream::Level3()) {
	    slivr::CheckMemoryHere(message, DebugStream::Stream3());
	}       
    }
    else if (debugN.compare("debug2") == 0) {
	if (DebugStream::Level2()) {
	    slivr::CheckMemoryHere(message, DebugStream::Stream2());
	}       
    }
    else if (debugN.compare("debug1") == 0) {
	if (DebugStream::Level1()) {
	    slivr::CheckMemoryHere(message, DebugStream::Stream1());
	}       
    }
}

void slivr::CheckMemoryHere(const std::string& message, std::ostream& out)
{
    unsigned long m_size, m_rss;
    avtMemory::GetMemorySize(m_size, m_rss);
    out << message << std::endl << "\t"
	<< " Rank " << PAR_Rank()
	<< " Memory use begin " << m_size 
	<< " rss " << m_rss/(1024*1024) << " (MB)"
	<< std::endl;
}
