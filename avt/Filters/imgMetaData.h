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
#ifndef IMG_METADATA_H
#define IMG_METADATA_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <limits>

#include "ospray/ospray.h"
#include "ospray/ospcommon/vec.h"
#include "ospray/ospcommon/math.h"

#include <vtkType.h>
#include <avtOpacityMap.h>

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

//! this struct will contain informations related to one volume
struct VolumeInfo 
{

    // meta info
    int           patchId = -1;
    bool          isComplete = false;

    // object references (shouldnt be deleted in this struct)
    OSPTransferFunction transferfcn = nullptr;
    OSPRenderer         renderer    = nullptr;

    // objects owned by the struct
    OSPModel                world = nullptr;
    unsigned char           worldType = OSP_INVALID;
    OSPFrameBuffer          framebuffer = nullptr;
    unsigned char           framebufferType = OSP_INVALID;
    float                  *framebufferData = nullptr;
    OSPVolume               volume = nullptr;
    unsigned char           volumeType = OSP_INVALID;
    void*                   dataPtr;
    std::string             dataType;
    size_t                  voxelSize;
    OSPDataType             voxelDataType;
    OSPData                 voxelData = nullptr;

    vec3f regionStart;
    vec3f regionStop;
    vec3f regionSpacing;
    vec3i regionSize;
    vec3f regionUpperClip;
    vec3f regionLowerClip;

    //! constructor
    VolumeInfo(int id) : patchId(id) {}
    ~VolumeInfo() {
	if (framebuffer != nullptr) { 
	    ospRelease(framebuffer); 
	}
	if (world != nullptr) { 
	    ospRelease(world); 
	}
	if (volume != nullptr) { 
	    ospRelease(volume); 
	}
    }

    //! ospTransfer function
    void SetTransferFunction(OSPTransferFunction tf) { transferfcn = tf; }
    void SetRenderer(OSPRenderer r) { renderer = r; }
    void SetCompleteFlag(bool f) { isComplete = f; }

    //! ospModel component
    void InitWorld() {
	if (worldType == OSP_INVALID) {
	    std::cout << "-- initializing world " << patchId << std::endl;
	    world = ospNewModel();
	    worldType = OSP_VALID;
	}
    }
    void SetWorld() {
	ospAddVolume(world, volume);
	ospCommit(world);
    }
    OSPModel GetWorld() { return world; }
	
    //! ospVolume component
    void InitVolume(unsigned char type = OSP_SHARED_STRUCTURED_VOLUME) {
	if (volumeType != type) { // only initialize once
	    std::cout << "-- initializing volume " << patchId << std::endl;
	    volumeType = type;	    
	    ospRelease(volume);
	    switch (type) {
	    case (OSP_BLOCK_BRICKED_VOLUME):
		volume = ospNewVolume("block_bricked_volume"); 
		break;
	    case (OSP_SHARED_STRUCTURED_VOLUME):
		volume = ospNewVolume("shared_structured_volume"); 
		break;
	    default:
		cout << "ERROR: wrong ospray volume type" << std::endl;
		volumeType = OSP_INVALID;
	    }
	}
    }
    void SetVolume
    (void* ptr, int type, double *X, double *Y, double *Z, int nX, int nY, int nZ) {
	if (isComplete) { return; }
	//! calculate volume data type
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
	    EXCEPTION1(VisItException, "ERROR: Unsupported ospray volume type");
	}
	//! assign data pointer
	dataPtr = ptr;
	//! assign structure
	regionStart   = vec3f(X[0],    Y[0],    Z[0]);
	regionStop    = vec3f(X[nX-1], Y[nY-1], Z[nZ-1]);
	regionSize    = vec3i(nX, nY, nZ);
	regionSpacing = (regionStop - regionStart)/((ospcommon::vec3f)regionSize - 1.0f);
	regionUpperClip = vec3f(X[0],Y[0],Z[0]);
	regionLowerClip = vec3f(X[nX-2], Y[nY-2], Z[nZ-2]);
	//! commit data
	voxelSize = nX * nY * nZ;
	voxelData = ospNewData(voxelSize, voxelDataType, dataPtr, OSP_DATA_SHARED_BUFFER);
	ospSetData(volume, "voxelData", voxelData);
	ospSetString(volume, "voxelType", dataType.c_str());
	//! commit volume
	ospSetVec3f(volume, "volumeClippingBoxLower", (const osp::vec3f&)regionLowerClip);
	ospSetVec3f(volume, "volumeClippingBoxUpper", (const osp::vec3f&)regionUpperClip);
	ospSetVec3f(volume, "gridOrigin",  (const osp::vec3f&)regionStart);
	ospSetVec3f(volume, "gridSpacing", (const osp::vec3f&)regionSpacing);
	ospSetVec3i(volume, "dimensions", (const osp::vec3i&)regionSize);
	ospSetObject(volume, "transferFunction", transferfcn);
	ospSet1i(volume, "gradientShadingEnabled", 0);
	ospSet1i(volume, "preIntegration", 0);
	ospSet1i(volume, "singleShade", 1);
	ospSet1i(volume, "adaptiveSampling", 0);
	ospSetVec3f(volume, "specular", osp::vec3f{1.0f,1.0f,1.0f});
	ospCommit(volume);
	isComplete = true;
    }
    void SetSamplingRate(float r) {
	ospSet1f(volume, "samplingRate", r);
	ospCommit(volume);
    }
    void FreeVolume() {	
	ospRelease(voxelData);
	ospRelease(volume);
    }

    //! framebuffer component     
    void InitFB(unsigned int width, unsigned int height) {
	std::cout << "-- initializing fb " << patchId << std::endl;
	vec2i imageSize(width, height);
	if (framebuffer != nullptr) {
	    // ospFreeFrameBuffer(framebuffer); 
	    // ospRelease(framebuffer); 
	}
	framebuffer = ospNewFrameBuffer((osp::vec2i&)imageSize, 
					OSP_FB_RGBA32F, OSP_FB_COLOR | OSP_FB_ACCUM);	    
    }
    void RenderFB() {
	ospFrameBufferClear(framebuffer, OSP_FB_COLOR | OSP_FB_ACCUM);
	ospRenderFrame(framebuffer, renderer, OSP_FB_COLOR | OSP_FB_ACCUM);
	framebufferData = (float*) ospMapFrameBuffer(framebuffer, OSP_FB_COLOR);
    }
    float* GetFB() {
	return framebufferData;
    }
    void CleanFBData() {
	ospUnmapFrameBuffer(framebufferData, framebuffer); 
    }
    void CleanFB() {	   
	ospFreeFrameBuffer(framebuffer); 
    }
};

struct OSPContext 
{
    bool refreshData = false;

    std::vector<VolumeInfo> volumePatch;
    OSPRenderer             renderer = nullptr;
    unsigned char           rendererType = OSP_INVALID;
    OSPCamera               camera = nullptr;
    unsigned char           cameraType = OSP_INVALID;
    OSPTransferFunction     transferfcn = nullptr;
    unsigned char           transferfcnType = OSP_INVALID;

    float r_panx;
    float r_pany;
    int   screenSize[2];
    bool  enabledOSPRay = false;

    ~OSPContext() {
	std::cout << "deleting ospray" << std::endl;    
	if (camera != nullptr) { ospRelease(camera); }
	if (transferfcn != nullptr) { ospRelease(transferfcn); }
	if (renderer != nullptr) { ospRelease(renderer); }
    }
    
    bool IsEnabled() { return enabledOSPRay; }

    void InitOSP(bool flag, bool debug = false, int numThreads = -1) { 
	enabledOSPRay = true;
	OSPDevice device = ospGetCurrentDevice();
	if (device == nullptr) {
	    std::cout << "Initializing OSPRay" << std::endl;
	    device = ospCreateDevice();
	    ospDeviceSet1i(device, "debug", debug ? 1 : 0);
	    if (numThreads != -1) {
		ospDeviceSet1i(device, "numThreads", numThreads);
	    }
	    ospDeviceCommit(device);
	    ospSetCurrentDevice(device);
	    ospDeviceSetErrorMsgFunc
		(device, [](const char *msg) { std::cout << msg; });
	}
	refreshData = flag;   
    }

    void InitPatch(int id) {
        if (volumePatch.size() < id) {
	    std::cerr << "ERROR: wrong patch index " << id << std::endl;
	    EXCEPTION1(VisItException, "ERROR: wrong patch index"); 
	    return;
	}
	if (volumePatch.size() == id) { 
	    volumePatch.emplace_back(id); 
	}
	volumePatch[id].SetTransferFunction(transferfcn);
	volumePatch[id].SetRenderer(renderer);
	volumePatch[id].SetCompleteFlag(!refreshData); // if the data is refreshed -> not complete
    }
    VolumeInfo* GetPatch(int id) {
	return &volumePatch[id];
    }

    //! ospRenderer component
    void InitRenderer() {
	if (rendererType == OSP_INVALID) {
	    renderer = ospNewRenderer("scivis");
	    rendererType = OSP_VALID;
	}
    }
    void SetRenderer(bool lighting, double material[4], double dir[3]) {
	ospSetObject(renderer, "camera", camera);
	ospSet1i(renderer, "backgroundEnabled", 0);
	ospSet1i(renderer, "oneSidedLighting", 0);
	ospSet1i(renderer, "aoSamples", 0);
	if (lighting == true)
	{
	    ospSet1i(renderer, "shadowsEnabled", 0);
	    OSPLight aLight = ospNewLight(renderer, "AmbientLight");
	    ospSet1f(aLight, "intensity", material[0]);
	    ospCommit(aLight);
	    OSPLight dLight = ospNewLight(renderer, "DirectionalLight");
	    ospSet1f(dLight, "intensity", material[2]);
	    ospSetVec3f(dLight, "direction", osp::vec3f{(float)dir[0],(float)dir[1],(float)dir[2]});
	    ospCommit(dLight);
	    OSPLight lights[2] = { aLight, dLight };
	    ospSetData(renderer,"lights",
		       ospNewData(2, OSP_OBJECT, lights));
	}
	ospCommit(renderer);
    }
    void SetModel(OSPModel world) {
	ospSetObject(renderer, "model",  world);
	ospCommit(renderer);
    }

    //! ospCamera component
    void InitCamera(unsigned char type) {
	if (cameraType != type) {
	    cameraType = type;
	    ospRelease(camera);
	    switch (type) {
	    case (OSP_PERSPECTIVE):
                camera = ospNewCamera("perspective");
		break;
	    case (OSP_ORTHOGRAPHIC):
                camera = ospNewCamera("orthographic");
		break;
	    default:
		cout << "ERROR: wrong ospray camera type" << std::endl;
	        cameraType = OSP_INVALID;
	    }
	}
    }
    void SetCamera(const double campos[3], 
		   const double camfocus[3], 
		   const double camup [3], 
		   const double camdir[3],
		   const double aspect, 
		   const double fovy, 
		   const double zoomratio, 
		   const double imagepan[2],
		   const int imageExtents[4],
		   const int screenExtents[2]) {
	float current[3];
	for (int i = 0; i < 3; ++i) {
	    current[i] = (campos[i] - camfocus[i]) / zoomratio + camfocus[i];
	}
	const ospcommon::vec3f camPos(current[0], current[1], current[2]);
	const ospcommon::vec3f camDir(camdir[0], camdir[1], camdir[2]);
	const ospcommon::vec3f camUp (camup[0], camup[1], camup[2]);
	ospSetf(camera, "aspect", aspect);
	ospSetVec3f(camera, "pos", (osp::vec3f&)camPos);
	ospSetVec3f(camera, "dir", (osp::vec3f&)camDir);
	ospSetVec3f(camera, "up",  (osp::vec3f&)camUp);
	ospSet1f(camera, "fovy", fovy);
	r_panx = imagepan[0] * zoomratio;
	r_pany = imagepan[1] * zoomratio;
	float r_xl = (float)imageExtents[0]/(float)screenExtents[0] - r_panx; 
	float r_yl = (float)imageExtents[2]/(float)screenExtents[1] - r_pany; 
	float r_xu = (float)imageExtents[1]/(float)screenExtents[0] - r_panx;
	float r_yu = (float)imageExtents[3]/(float)screenExtents[1] - r_pany;
	ospSetVec2f(camera, "imageStart", osp::vec2f{r_xl, r_yl});
	ospSetVec2f(camera, "imageEnd",   osp::vec2f{r_xu, r_yu});
	ospCommit(camera);
	screenSize[0] = screenExtents[0];
	screenSize[1] = screenExtents[1];
    }
    void SetSubCamera(float xMin, float xMax, float yMin, float yMax) {
	float r_xl = xMin/screenSize[0] - r_panx; 
	float r_yl = yMin/screenSize[1] - r_pany; 
	float r_xu = xMax/screenSize[0] - r_panx;
	float r_yu = yMax/screenSize[1] - r_pany;	
	ospSetVec2f(camera, "imageStart", osp::vec2f{r_xl, r_yl});
	ospSetVec2f(camera, "imageEnd",   osp::vec2f{r_xu, r_yu});
	ospCommit(camera);
    }

    //! ospTransferFunction component  
    void InitTransferFunction() {
	if (transferfcnType == OSP_INVALID) {
	    transferfcn = ospNewTransferFunction("piecewise_linear");
	    transferfcnType = OSP_VALID;
	}
    }
    void SetTransferFunction(const RGBAF* table, 
			     const unsigned int size, 
			     const float datamin, 
			     const float datamax) {
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
};

// ****************************************************************************
//  Struct:  imgMetaData
//
//  Purpose:
//    Holds information about patches but not the image 
//
//  Programmer:  
//  Creation:   
//
// ****************************************************************************
struct imgMetaData
{
    int procId;       // processor that produced the patch
    int patchNumber;  // id of the patch on that processor - with procId, acts as a key
    int destProcId;   // destination proc where this patch gets composited
    int inUse;        // whether the patch is composed locally or not
    int dims[2];      // height, width
    int screen_ll[2]; // position in the final image
    int screen_ur[2];
    float avg_z;      // camera space z = depth of the patch - used for compositing
    float eye_z;      // camera space z
    float clip_z;     // clip space z
};


// ****************************************************************************
//  Struct:  imgData
//
//  Purpose:
//    Holds the image data generated
//
//  Programmer:  
//  Creation:    
//
// ****************************************************************************
struct imgData
{
    int procId;         // processor that produced the patch
    int patchNumber;    // id of the patch on that processor  - with procId, acts as a key

    float *imagePatch;  // the image data - RGBA

    bool operator==(const imgData &a){
        return (patchNumber == a.patchNumber);
    }
};


// ****************************************************************************
//  Struct:  convexHull
//
//  Purpose:
//    Holds the image data generated
//
//  Programmer:  
//  Creation:    
//
// ****************************************************************************
struct convexHull
{
    int numPatches;
    int arrangement[3];     // [0] rows along x axis, [1] rows along y axis, [2] rows along z axis

    float extents[6];       // minX, maxX   minY, maxY   minZ, maxZ
    float cellDims[3];      // x, y, z
    float tolerance;        // amount of overlap that is considered ok - typically 2 cells for cell centered data

    // 0: no overlap    1: overlpa in Z    2: overlap in Y    3: overlap in Z
    int overlap(convexHull _hull)
    {

        if ( (_hull.extents[1] < extents[0]) || (_hull.extents[0] > extents[1]) )   // No overlap in X
        {
            if ( (_hull.extents[3] < extents[2]) || (_hull.extents[2] > extents[3]) )   // No overlap in Y
            {
                if ( (_hull.extents[5] < extents[4]) || (_hull.extents[4] > extents[5]) )   // No overlap in Z
                {
                    return 0;
                }
                else
                    return 3;
            }
            else
                return 2;
        }
        else
            return 1;
    }
};


template <class T> 
inline std::string toStr(T x){
    std::ostringstream ss;
    ss << x;
    return ss.str();
}


inline void createColorPPM(std::string filename, unsigned char *data, int width, int height){
    std::ofstream outputFile( (filename+ ".ppm").c_str(), std::ios::out | std::ios::binary);
    outputFile << "P6\n" << width << "\n" << height << "\n" << 255 << "\n";
    
    for (int y=0; y<height; ++y){
        for (int x=0; x<width; ++x){
            int index = ((y * width) + x)*3;
            
            char color[3];
            color[0] = data[index+0];  // red
            color[1] = data[index+1];  // green   
            color[2] = data[index+2];  // blue
            outputFile.write(color,3);
        }
    }
    
    outputFile.close();
}


inline void writeOutputToFile( std::string filename, float * data, int dimX, int dimY )
{
    std::ofstream outputFile( (filename+ ".txt").c_str(), std::ios::out);
    outputFile << "Dims: " << dimX << ", " << dimY << "\n";
 
    for (int y=0; y<dimY; ++y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x);
            outputFile << data[index] << ", ";
        }
        outputFile  << "\n";
    }
 
    outputFile.close();
}


inline void writeOutputToFileByLine(std::string filename, float * data, int dimX, int dimY )
{
    std::ofstream outputFile( (filename+ ".txt").c_str(), std::ios::out);
    outputFile << "Dims: " << dimX << ", " << dimY << "\n";
 
    for (int y=0; y<dimY; ++y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x);
            outputFile << x  << ", " << y << " : " << data[index] << "\n";
        }
        outputFile  << "\n";
    }
 
    outputFile.close();
}

inline void writeDepthBufferToPPM( std::string filename , float * data, int dimX, int dimY )
{
    std::ofstream outputFile( (filename+ ".ppm").c_str(), std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n";
 
    for (int y=0; y<dimY; ++y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x);
 
            char color[3];
            color[0] = data[index] * 255;  // red
            color[1] = data[index] * 255;  // green
            color[2] = data[index] * 255;  // blue
            outputFile.write(color,3);
        }
    }
 
    outputFile.close();
}


inline void writeArrayToPPM( std::string filename , float * image, int dimX, int dimY )
{
    std::ofstream outputFile( (filename+ ".ppm").c_str(), std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n";
 
    for (int y=0; y<dimY; ++y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x)*4;
 
            char color[3];
            float alpha = image[index + 3];
            color[0] = std::max( std::min(image[index + 0]*alpha , 1.0f), 0.0f) * 255;  // red
            color[1] = std::max( std::min(image[index + 1]*alpha , 1.0f), 0.0f) * 255;  // green
            color[2] = std::max( std::min(image[index + 2]*alpha , 1.0f), 0.0f) * 255;  // blue
            outputFile.write(color,3);
        }
    }
 
    outputFile.close();
}

#endif
