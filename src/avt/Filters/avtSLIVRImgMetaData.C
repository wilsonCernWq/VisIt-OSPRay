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

#include "avtSLIVRImgMetaData.h"

#include <avtMemory.h>
#include <avtParallel.h>
#include <ImproperUseException.h>
#include <DebugStream.h>
#include <TimingsManager.h>

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
double slivr::ProjectWorldToScreen
(const double worldCoord[3], const int screenWidth, const int screenHeight,
 const double panPercentage[2], const double imageZoom,
 vtkMatrix4x4 *mvp, int screenCoord[2])
{
    // world space coordinate in homogeneous coordinate
    double worldHCoord[4] = {
	worldCoord[0],
	worldCoord[1],
	worldCoord[2],
	1.0
    };

    // world to clip space (-1 ~ 1)
    double clipHCoord[4];
    mvp->MultiplyPoint(worldHCoord, clipHCoord);
    if (clipHCoord[3] == 0.0)
    {
	std::cerr << "slivr::ProjectWorldToScreen "
		  << "Zero Division During Projection" 
		  << std::endl;
	std::cerr << "world coordinates: (" 
		  << worldHCoord[0] << ", " 
		  << worldHCoord[1] << ", " 
		  << worldHCoord[2] << ", " 
		  << worldHCoord[3] << ")" << std::endl
		  << "clip space coordinate: ("
		  << clipHCoord[0] << ", " 
		  << clipHCoord[1] << ", " 
		  << clipHCoord[2] << ", "
		  << clipHCoord[3] << std::endl;
	std::cerr << "Matrix: " << *mvp << std::endl;
	EXCEPTION1(VisItException, "Zero Division During Projection");
    }

    // normalize clip space coordinate
    double clipCoord[3] = {
	clipHCoord[0]/clipHCoord[3],
	clipHCoord[1]/clipHCoord[3],
	clipHCoord[2]/clipHCoord[3]
    };

    // screen coordinates (int integer)
    screenCoord[0] = round(clipCoord[0]*(screenWidth /2.0)+(screenWidth /2.0));
    screenCoord[1] = round(clipCoord[1]*(screenHeight/2.0)+(screenHeight/2.0));

    // add panning
    screenCoord[0] += round(screenWidth  * panPercentage[0] * imageZoom);
    screenCoord[1] += round(screenHeight * panPercentage[1] * imageZoom); 

    // return point depth
    return clipCoord[2];
}

void
slivr::ProjectScreenToWorld
(const int screenCoord[2], const double z,
 const int screenWidth, const int screenHeight, 
 const double panPercentage[2], const double imageZoom,
 vtkMatrix4x4 *imvp, double worldCoord[3])
{
    // remove panning
    const int x = 
	screenCoord[0] - round(screenWidth*panPercentage[0]*imageZoom);
    const int y = 
	screenCoord[1] - round(screenHeight*panPercentage[1]*imageZoom);
    
    // do projection
    double worldHCoord[4] = {0,0,0,1};
    double clipHCoord[4] = {
	(x - screenWidth/2.0)/(screenWidth/2.0),
	(y - screenHeight/2.0)/(screenHeight/2.0),
	z, 1.0};
    imvp->MultiplyPoint(clipHCoord, worldHCoord);
    if (worldHCoord[3] == 0) {
	debug5 << "slivr::ProjectScreenToWorld "
	       << "Zero Division During Projection" 
	       << std::endl;
	std::cerr << "world coordinates: (" 
		  << worldHCoord[0] << ", " 
		  << worldHCoord[1] << ", " 
		  << worldHCoord[2] << ", " 
		  << worldHCoord[3] << ")" << std::endl
		  << "clip space coordinate: ("
		  << clipHCoord[0] << ", " 
		  << clipHCoord[1] << ", " 
		  << clipHCoord[2] << ", "
		  << clipHCoord[3] << std::endl;
	std::cerr << "Matrix: " << *imvp << std::endl;
	EXCEPTION1(VisItException, "Zero Division During Projection");
    }
    
    // normalize world space coordinate	
    worldCoord[0] = worldHCoord[0]/worldHCoord[3];
    worldCoord[1] = worldHCoord[1]/worldHCoord[3];
    worldCoord[2] = worldHCoord[2]/worldHCoord[3];
}

void
slivr::ProjectWorldToScreenCube
(const double cube[6], const int screenWidth, const int screenHeight, 
 const double panPercentage[2], const double imageZoom, vtkMatrix4x4 *mvp, 
 int screenExtents[4], double depthExtents[2])
{
    int xMin = std::numeric_limits<int>::max();
    int xMax = std::numeric_limits<int>::min();
    int yMin = std::numeric_limits<int>::max();
    int yMax = std::numeric_limits<int>::min();
    double zMin = std::numeric_limits<double>::max();
    double zMax = std::numeric_limits<double>::min();

    float coordinates[8][3];
    coordinates[0][0] = cube[0];   
    coordinates[0][1] = cube[2];   
    coordinates[0][2] = cube[4];	

    coordinates[1][0] = cube[1];   
    coordinates[1][1] = cube[2];   
    coordinates[1][2] = cube[4];	

    coordinates[2][0] = cube[1];  
    coordinates[2][1] = cube[3];
    coordinates[2][2] = cube[4];	

    coordinates[3][0] = cube[0]; 
    coordinates[3][1] = cube[3]; 
    coordinates[3][2] = cube[4];

    coordinates[4][0] = cube[0];
    coordinates[4][1] = cube[2];
    coordinates[4][2] = cube[5];

    coordinates[5][0] = cube[1]; 
    coordinates[5][1] = cube[2]; 
    coordinates[5][2] = cube[5];	

    coordinates[6][0] = cube[1]; 
    coordinates[6][1] = cube[3];
    coordinates[6][2] = cube[5];

    coordinates[7][0] = cube[0]; 
    coordinates[7][1] = cube[3]; 
    coordinates[7][2] = cube[5];

    double worldCoord[3];
    int screenCoord[2]; double depth;
    for (int i=0; i<8; i++)
    {
	worldCoord[0] = coordinates[i][0];
	worldCoord[1] = coordinates[i][1];
	worldCoord[2] = coordinates[i][2];
	depth = slivr::ProjectWorldToScreen
	    (worldCoord, screenWidth, screenHeight, 
	     panPercentage, imageZoom, mvp, screenCoord);
	// clamp values
	screenCoord[0] = CLAMP(screenCoord[0], 0, screenWidth);
	screenCoord[1] = CLAMP(screenCoord[1], 0, screenHeight);
	screenExtents[0] = xMin = std::min(xMin, screenCoord[0]);
	screenExtents[1] = xMax = std::max(xMax, screenCoord[0]);
	screenExtents[2] = yMin = std::min(yMin, screenCoord[1]);
	screenExtents[3] = yMax = std::max(yMax, screenCoord[1]);
	depthExtents[0] = zMin = std::min(zMin, depth);
	depthExtents[1] = zMax = std::max(zMax, depth);
    }
}


// ****************************************************************************
//  Function:  
//
//  Purpose:
//    
//
//  Programmer:  
//  Creation:    
//
// ****************************************************************************
void 
CreateColorPPM
(std::string filename, unsigned char *data, int width, int height)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
			     std::ios::out | std::ios::binary);
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

void 
WriteOutputToFile( std::string filename, float * data, int dimX, int dimY )
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

void 
WriteOutputToFileByLine(std::string filename, float * data, int dimX, int dimY)
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

void WriteDepthBufferToPPM
(std::string filename , float * data, int dimX, int dimY)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
			     std::ios::out | std::ios::binary);
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

void WriteArrayToPPM(std::string filename , float * image, int dimX, int dimY)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
			     std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n"; 
    for (int y=0; y<dimY; ++y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x)*4;
            char color[3];
            float alpha = image[index + 3];
            color[0] = CLAMP(image[index + 0]*alpha, 0.0f, 1.0f) * 255;
            color[1] = CLAMP(image[index + 1]*alpha, 0.0f, 1.0f) * 255;
            color[2] = CLAMP(image[index + 2]*alpha, 0.0f, 1.0f) * 255;
            outputFile.write(color,3);
        }
    } 
    outputFile.close();
}

void WriteArrayToPPM(std::string filename , unsigned char * image, 
		     int dimX, int dimY)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
			     std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n"; 
    // for (int y=0; y<dimY; ++y)
    // {
    //     for (int x=0; x<dimX; ++x)
    //     {
    //         int index = (y * dimX + x)*3;
    //         char color[3];
    //         color[0] = CLAMP(image[index + 0], 0, 255);
    //         color[1] = CLAMP(image[index + 1], 0, 255);
    //         color[2] = CLAMP(image[index + 2], 0, 255);
    //         outputFile.write(color,3);
    //     }
    // } 
    outputFile.write((const char*)image, dimX * dimY * 3);
    outputFile.close();
}

void WriteArrayGrayToPPM(std::string filename , float * image, 
			 int dimX, int dimY)
{
    std::ofstream outputFile((filename+ ".ppm").c_str(), 
			     std::ios::out | std::ios::binary);
    outputFile <<  "P6\n" << dimX << "\n" << dimY << "\n" << 255 << "\n"; 
    for (int y=0; y<dimY; ++y)
    {
        for (int x=0; x<dimX; ++x)
        {
            int index = (y * dimX + x);
            char color[3];
            color[0] = CLAMP(image[index], 0.f, 1.f) * 255;
            color[1] = CLAMP(image[index], 0.f, 1.f) * 255;
            color[2] = CLAMP(image[index], 0.f, 1.f) * 255;
            outputFile.write(color,3);
        }
    } 
    outputFile.close();
}

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
int slivr::ConvexHull::Overlap(ConvexHull _hull) 
{
    if ( (_hull.extents[1] < extents[0]) || 
	 (_hull.extents[0] > extents[1]) )   // No overlap in X
    {
	if ( (_hull.extents[3] < extents[2]) || 
	     (_hull.extents[2] > extents[3]) )   // No overlap in Y
	{
	    if ( (_hull.extents[5] < extents[4]) ||
		 (_hull.extents[4] > extents[5]) )   // No overlap in Z
	    {
		return 0;
	    }
	    else
	    {
		return 3;
	    }
	}
	else
	{
	    return 2;
	}
    }
    else
    {
	return 1;
    }
}
