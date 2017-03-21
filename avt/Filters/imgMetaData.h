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

#include "ospray/ospray.h"
#include "ospray/ospcommon/vec.h"

#include <vtkType.h>
#include <avtOpacityMap.h>

#define OSP_PERSPECTIVE              1
#define OSP_ORTHOGRAPHIC             2
#define OSP_BLOCK_BRICKED_VOLUME     3
#define OSP_SHARED_STRUCTURED_VOLUME 4
#define OSP_INVALID                  5
#define OSP_VALID                    6

struct OSPContext {

    ~OSPContext() {
	ospUnmapFrameBuffer(framebufferData, framebuffer);
	ospRelease(camera);
	ospRelease(transferfcn);
	ospRelease(world);
	ospRelease(renderer);
	ospRelease(framebuffer);
    }

    //! framebuffer component
    OSPFrameBuffer framebuffer;
    float         *framebufferData = NULL;
    void InitFB(unsigned int width, unsigned int height) {
	ospcommon::vec2i imageSize(width, height);	
	framebuffer = ospNewFrameBuffer((osp::vec2i&)imageSize, 
					OSP_FB_RGBA32F, 
					OSP_FB_COLOR | OSP_FB_ACCUM);
	ospFrameBufferClear(framebuffer, OSP_FB_COLOR | OSP_FB_ACCUM);
    }
    void RenderFB() {
	if (framebufferData != NULL) { ospUnmapFrameBuffer(framebufferData, framebuffer); }
	ospRenderFrame(framebuffer, renderer, OSP_FB_COLOR | OSP_FB_ACCUM);
        framebufferData = (float*) ospMapFrameBuffer(framebuffer, OSP_FB_COLOR);
    }
    float* GetData() {
	return framebufferData;
    }

    //! ospRenderer component
    OSPRenderer         renderer;
    unsigned char       rendererType = OSP_INVALID;
    void InitRenderer() {
	if (rendererType == OSP_INVALID) {
	    renderer = ospNewRenderer("scivis");
	    rendererType = OSP_VALID;
	}
    }
    void SetRenderer(bool lighting, float material[4], float dir[3]) {
	ospSetObject(renderer, "camera", camera);
	ospSetObject(renderer, "model",  world);
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

    //! ospModel component
    OSPModel            world;
    unsigned char       worldType = OSP_INVALID;
    void InitWorld() {
	if (worldType == OSP_INVALID) {
	    world = ospNewModel();
	    worldType = OSP_VALID;
	}
    }
    void SetWorld() {
	ospAddVolume(world, volume);
	ospCommit(world);
    }

    //! ospVolume component
    OSPVolume           volume;
    unsigned char       volumeType = OSP_INVALID;
    void InitVolume(unsigned char type) {
	if (volumeType != type) {
	    volumeType = type;	    
            ospRelease(volume);
	    switch (type) {
	    case (OSP_BLOCK_BRICKED_VOLUME):
		volume = ospNewVolume("block_bricked_volume"); 
		break;
	    case (OSP_SHARED_STRUCTURED_VOLUME):
		volume = ospNewVolume("block_bricked_volume"); 
		break;
	    default:
		cout << "ERROR: wrong ospray volume type" << std::endl;
		volumeType = OSP_INVALID;
	    }
	}
    }
    void SetVolume(void* dataPtr, int dataType) {
/* 	//! calculate volume data type */
/* 	std::string voxelType; */
/* 	OSPDataType voxelDataType; */
/* 	if (dataType == VTK_UNSIGNED_CHAR) { */
/* 	    voxelType = "uchar"; */
/* 	    voxelDataType = OSP_UCHAR; */
/* 	} else if (dataType == VTK_SHORT) { */
/* 	    voxelType = "short"; */
/* 	    voxelDataType = OSP_SHORT; */
/* 	} else if (dataType == VTK_UNSIGNED_SHORT) { */
/* 	    voxelType = "ushort"; */
/* 	    voxelDataType = OSP_USHORT; */
/* 	} else if (dataType == VTK_FLOAT) { */
/* 	    voxelType = "float"; */
/* 	    voxelDataType = OSP_FLOAT; */
/* 	} else if (dataType == VTK_DOUBLE) {	 */
/* 	    voxelType = "double"; */
/* 	    voxelDataType = OSP_DOUBLE; */
/* 	} else { */
/* 	    EXCEPTION1(VisItException, "ERROR: Unsupported ospray volume type"); */
/* 	} */
    }

    //! ospCamera component
    OSPCamera           camera;
    unsigned char       cameraType = OSP_INVALID;
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
    void SetCamera(const float campos[3], 
		   const float camfocus[3], 
		   const float camup [3], 
		   const float camdir[3],
		   const float aspect, 
		   const float fovy, 
		   const float zoomratio, 
		   const float imagepan[2],
		   const int imageExtents[2],
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
	float r_panx = imagepan[0] * zoomratio;
	float r_pany = imagepan[1] * zoomratio;
	float r_xl = (float)imageExtents[0]/(float)screenExtents[0] - r_panx; 
	float r_yl = (float)imageExtents[2]/(float)screenExtents[1] - r_pany; 
	float r_xu = (float)imageExtents[1]/(float)screenExtents[0] - r_panx;
	float r_yu = (float)imageExtents[3]/(float)screenExtents[1] - r_pany;
	ospSetVec2f(camera, "imageStart", osp::vec2f{r_xl, r_yl});
	ospSetVec2f(camera, "imageEnd",   osp::vec2f{r_xu, r_yu});
	ospCommit(camera);
    }

    //! ospTransferFunction component
    OSPTransferFunction transferfcn;
    unsigned char       transferfcnType = OSP_INVALID;  
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

struct ospVolumeMeta {

    std::string ospVoxelType;
    OSPVolume volume;
    void* ospVolumePointer;
    OSPDataType ospVoxelDataType;
    size_t ospVolumeSize;
    OSPData ospVoxelData;

    ospcommon::vec3i volumeDims;
    ospcommon::vec3f volumeLbox;
    ospcommon::vec3f volumeMbox;
    ospcommon::vec3f volumeSpac;

    void init() { 
	//volume = ospNewVolume("shared_structured_volume"); 
	volume = ospNewVolume("block_bricked_volume"); 
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
