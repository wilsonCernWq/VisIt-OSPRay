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
#include <DebugStream.h>
#include <ImproperUseException.h>
#include <TimingsManager.h>

#include <cmath>

// helper
double slivr::deg2rad (double degrees) {
    return degrees * 4.0 * atan (1.0) / 180.0;
}

// other function
void 
VolumeInfo::Set
(void *ptr, int type, double *X, double *Y, double *Z, 
 int nX, int nY, int nZ, float sr, double volumePBox[6], double volumeBBox[6])
{
    if (!isComplete) {
	worldType = OSP_INVALID;
	volumeType = OSP_INVALID;
    }
    InitWorld();
    InitVolume();
    if (!isComplete) { 
	SetVolume(ptr, type, X, Y, Z, nX, nY, nZ, volumePBox, volumeBBox); 
    }
    if (samplingRate != sr) {
	samplingRate = sr;
	SetSamplingRate(samplingRate);
    }
    if (!isComplete) { 
	SetWorld();
    }
    isComplete = true;	
}

// ospModel component
void VolumeInfo::InitWorld() {
    if (worldType == OSP_INVALID) {
	CleanWorld();
	worldType = OSP_VALID;
	world = ospNewModel();
    }
}
void VolumeInfo::SetWorld() {
    if (world != nullptr) { 
	ospAddVolume(world, volume);
	ospCommit(world);
    }
}
void VolumeInfo::CleanWorld() {
    if (world != nullptr) {	    
	ospRelease(world);
	world = nullptr;
    }
    worldType = OSP_INVALID;
}

// ospVolume component
void VolumeInfo::InitVolume(unsigned char type) {
    if (volumeType != type) { // only initialize once
	CleanVolume();
	volumeType = type;
	switch (type) {
	case (OSP_BLOCK_BRICKED_VOLUME):
	    volume = ospNewVolume("block_bricked_volume"); 
	    break;
	case (OSP_SHARED_STRUCTURED_VOLUME):
	    volume = ospNewVolume("shared_structured_volume"); 
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
void VolumeInfo::SetVolume(void *ptr, int type, 
			   double *X, double *Y, double *Z, 
			   int nX, int nY, int nZ,
			   double volumePBox[6], 
			   double volumeBBox[6]) {
    // refresh existing data
    if (voxelData != nullptr) { 
	ospRelease(voxelData); 
	voxelData = nullptr; 
    }
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
    // assign data pointer
    dataPtr = ptr;
    // assign structure
    regionStart   = vec3f(volumePBox[0], volumePBox[1], volumePBox[2]);
    regionStop    = vec3f(volumePBox[3], volumePBox[4], volumePBox[5]);
    regionSize    = vec3i(nX, nY, nZ);
    // regionSpacing = vec3f(X[1]-X[0], Y[1]-Y[0], Z[1]-Z[0]);
    regionSpacing = (regionStop-regionStart)/
	((ospcommon::vec3f)regionSize-1.0f);

    regionLowerClip.x = volumeBBox[0];
    regionLowerClip.y = volumeBBox[1];
    regionLowerClip.z = volumeBBox[2];

    regionUpperClip.x = volumeBBox[3];
    regionUpperClip.y = volumeBBox[4];
    regionUpperClip.z = volumeBBox[5];

    // commit data
    voxelSize = nX * nY * nZ;
    voxelData = ospNewData(voxelSize, voxelDataType,
			   dataPtr, OSP_DATA_SHARED_BUFFER);
    ospSetData(volume, "voxelData", voxelData);
    ospSetString(volume, "voxelType", dataType.c_str());
    ospSetObject(volume, "transferFunction", transferfcn);

    // commit volume
    ospSetVec3f(volume, "specular", osp::vec3f{0.3f,0.3f,0.3f});
    ospSetVec3f(volume,
		"volumeClippingBoxLower",
		(const osp::vec3f&)regionLowerClip);
    ospSetVec3f(volume,
		"volumeClippingBoxUpper",
		(const osp::vec3f&)regionUpperClip);
    ospSetVec3f(volume, "gridSpacing", (const osp::vec3f&)regionSpacing);
    ospSetVec3f(volume, "gridOrigin",  (const osp::vec3f&)regionStart);
    ospSetVec3i(volume, "dimensions",  (const osp::vec3i&)regionSize);
    ospSet1i(volume, "gradientShadingEnabled", 1);
    ospSet1i(volume, "adaptiveSampling", 0);
    ospSet1i(volume, "preIntegration", 0);
    ospSet1i(volume, "singleShade", 1);
    int volumeInitIndex = visitTimer->StartTimer();
    ospCommit(volume);
    visitTimer->StopTimer(volumeInitIndex, "Commit OSPRay patch");
    visitTimer->DumpTimings();
}
void VolumeInfo::SetSamplingRate(float r) {
    ospSet1f(volume, "samplingRate", r);
    ospCommit(volume);
}
void VolumeInfo::CleanVolume() {	
    if (voxelData != nullptr) { 
	ospRelease(voxelData);
	voxelData = nullptr; 
    }
    if (volume != nullptr) { ospRelease(volume); volume = nullptr; }
    volumeType = OSP_INVALID;
}

// framebuffer component     
void VolumeInfo::InitFB(unsigned int width, unsigned int height) {
    vec2i imageSize(width, height);
    CleanFBData(); CleanFB();	    
    framebuffer = ospNewFrameBuffer((osp::vec2i&)imageSize, 
				    OSP_FB_RGBA32F,
				    OSP_FB_COLOR | OSP_FB_ACCUM);	    
}
void VolumeInfo::RenderFB() {
    ospFrameBufferClear(framebuffer, OSP_FB_COLOR | OSP_FB_ACCUM);
    ospRenderFrame(framebuffer, renderer, OSP_FB_COLOR | OSP_FB_ACCUM);
    framebufferData = (float*) ospMapFrameBuffer(framebuffer, OSP_FB_COLOR);
}
float* VolumeInfo::GetFBData() {
    return framebufferData;
}
void VolumeInfo::CleanFBData() {
    if (framebufferData != nullptr) { 
	ospUnmapFrameBuffer(framebufferData, framebuffer); 
	framebufferData = nullptr;
    }
}
void VolumeInfo::CleanFB() {	   
    if (framebuffer != nullptr) { 
	ospFreeFrameBuffer(framebuffer); 	    
	framebuffer = nullptr;
    }	
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
OSPContext::~OSPContext() 
{
    std::cout << "deleting ospray" << std::endl;
    // check memory
    slivr::CheckMemoryHere("OSPContext::~OSPContext before deleting OSPRay");
    // clean stuffs
    volumePatch.clear();
    if (camera      != nullptr) { ospRelease(camera); }
    if (renderer    != nullptr) { ospRelease(renderer); }
    if (transferfcn != nullptr) { ospRelease(transferfcn); }
    // check memory again
    slivr::CheckMemoryHere("OSPContext::~OSPContext after deleting OSPRay");
}

void OSPContext::InitOSP(bool flag, bool debug, int numThreads) 
{ 
    std::cout << "Initialize OSPRay (new data " << flag << ")" << std::endl;
    enabledOSPRay = true;
    OSPDevice device = ospGetCurrentDevice();
    if (device == nullptr) {
	debug5 << "Initializing OSPRay" 
	       << " debug: " << debug  
	       << " numThreads: " << numThreads
	       << std::endl;
	device = ospCreateDevice();
	ospDeviceSet1i(device, "debug", debug ? 1 : 0);
	if (numThreads != -1) {
	    ospDeviceSet1i(device, "numThreads", numThreads);
	}
	ospDeviceCommit(device);
	ospSetCurrentDevice(device);
	ospDeviceSetStatusFunc
	    (device, [](const char *msg) { debug5 << msg; });
	ospDeviceCommit(device);
    }
    refreshData = flag;
}

void OSPContext::InitPatch(int id) 
{
    if (volumePatch.size() < id) {
	debug1 << "ERROR: wrong patch index " << id << std::endl;
	EXCEPTION1(VisItException, "ERROR: wrong patch index"); 
	return;
    }
    if (volumePatch.size() == id) { 
	volumePatch.emplace_back(id); 
    }
    volumePatch[id].SetTransferFunction(transferfcn);
    volumePatch[id].SetRenderer(renderer);
    volumePatch[id].SetCompleteFlag(!refreshData); 
    // if the data is refreshed -> not complete
}

// ospRenderer component
void OSPContext::InitRenderer() 
{
    if (rendererType == OSP_INVALID) {
	renderer = ospNewRenderer("scivis");
	rendererType = OSP_VALID;
    }
}

void OSPContext::SetRenderer(bool lighting, double material[4], double dir[3]) 
{
    ospSetObject(renderer, "camera", camera);
    ospSet1i(renderer, "backgroundEnabled", 0);
    ospSet1i(renderer, "oneSidedLighting", 0);
    ospSet1i(renderer, "aoSamples", 16);
    if (lighting)
    {
	std::cout << "use lighting" << std::endl;
	std::cout << "use material " 
		  << material[0] << " "
		  << material[1] << " "
		  << material[2] << " "
		  << material[3] << std::endl;
	ospSet1i(renderer, "shadowsEnabled", 1);
	OSPLight aLight = ospNewLight(renderer, "AmbientLight");
	ospSet1f(aLight, "intensity", material[0]);
	ospCommit(aLight);
	OSPLight dLight = ospNewLight(renderer, "DirectionalLight");
	ospSet1f(dLight, "intensity", material[1]);
	ospSetVec3f(dLight, "direction", 
		    osp::vec3f{(float)-dir[0],(float)-dir[1],(float)-dir[2]});
	ospCommit(dLight);
	OSPLight lights[2] = { aLight, dLight };
	ospSetData(renderer,"lights", ospNewData(2, OSP_OBJECT, lights));
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
	if (camera != nullptr) { ospRelease(camera); }
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
			   const double fovy, 
			   const double zoomratio, 
			   const double imagepan[2],
			   const int imageExtents[4],
			   const int screenExtents[2]) 
{
    float current[3];
    for (int i = 0; i < 3; ++i) {
	current[i] = (campos[i] - camfocus[i]) / zoomratio + camfocus[i];
    }
    const ospcommon::vec3f camPos(current[0], current[1], current[2]);
    const ospcommon::vec3f camDir(camdir[0], camdir[1], camdir[2]);
    const ospcommon::vec3f camUp (camup[0], camup[1], camup[2]);
    ospSetVec3f(camera, "pos", (osp::vec3f&)camPos);
    ospSetVec3f(camera, "dir", (osp::vec3f&)camDir);
    ospSetVec3f(camera, "up",  (osp::vec3f&)camUp);
    if (cameraType == OSP_PERSPECTIVE) {
	ospSet1f(camera, "aspect", aspect);
	ospSet1f(camera, "fovy", fovy);
    }
    else if (cameraType == OSP_ORTHOGRAPHIC) {
	ospSet1f(camera, "aspect", aspect);
	ospSet1f(camera, "height", sceneSize[1]);
    }
    r_panx = imagepan[0] * zoomratio;
    r_pany = imagepan[1] * zoomratio;
    this->SetSubCamera(imageExtents[0], imageExtents[1], imageExtents[2], imageExtents[3]);
    screenSize[0] = screenExtents[0];
    screenSize[1] = screenExtents[1];
}

void OSPContext::SetSubCamera(float xMin, float xMax, float yMin, float yMax) 
{
    float r_xl = xMin/screenSize[0] - r_panx; 
    float r_yl = yMin/screenSize[1] - r_pany; 
    float r_xu = xMax/screenSize[0] - r_panx;
    float r_yu = yMax/screenSize[1] - r_pany;	
    ospSetVec2f(camera, "imageStart", osp::vec2f{r_xl, r_yl});
    ospSetVec2f(camera, "imageEnd",   osp::vec2f{r_xu, r_yu});
    ospCommit(camera);
}

// ospTransferFunction component  
void OSPContext::InitTransferFunction() 
{
    if (transferfcnType == OSP_INVALID) {
	if (transferfcn != nullptr) { ospRelease(transferfcn); }
	transferfcn = ospNewTransferFunction("piecewise_linear");
	transferfcnType = OSP_VALID;
    }
}

void OSPContext::SetTransferFunction(const OSPColor *table,
				     const unsigned int size, 
				     const float datamin, 
				     const float datamax) 
{
    std::vector<ospcommon::vec3f> colors;
    std::vector<float>            opacities;
    for (int i = 0; i < size; ++i) {
	colors.emplace_back(table[i].R, table[i].G, table[i].B);
	opacities.emplace_back(table[i].A);
    }
    OSPData colorData   = 
	ospNewData(colors.size(), OSP_FLOAT3, colors.data());
    OSPData opacityData = 
	ospNewData(opacities.size(), OSP_FLOAT, opacities.data());
    const ospcommon::vec2f range(datamin, datamax);
    ospSetData(transferfcn, "colors",    colorData);
    ospSetData(transferfcn, "opacities", opacityData);
    ospSetVec2f(transferfcn, "valueRange", (osp::vec2f&)range);
    ospCommit(transferfcn);
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
