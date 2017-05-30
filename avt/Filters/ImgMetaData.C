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

#include "ImgMetaData.h"

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
(const double world_pos[3], const int screenWidth, const int screenHeight,
 const double panPercentage[2], const double imageZoom,
 vtkMatrix4x4 *mvp, int screen_pos[2])
{
    // world space coordinate in homogeneous coordinate
    double worldHCoord[4] = {
	world_pos[0],
	world_pos[1],
	world_pos[2],
	1.0
    };

    // world to clip space (-1 ~ 1)
    double clipHCoord[4];
    mvp->MultiplyPoint(worldHCoord, clipHCoord);
    if (clipHCoord[3] == 0.0)
    {
	std::cerr << "avtMassVoxelExtractor::Zero Division During Projection" 
		  << endl;
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
	std::cerr << "Matrix: " << *mvp << endl;
	EXCEPTION1(VisItException, "Zero Division During Projection");
    }

    // normalize clip space coordinate
    double clip_pos[3] = {
	clipHCoord[0]/clipHCoord[3],
	clipHCoord[1]/clipHCoord[3],
	clipHCoord[2]/clipHCoord[3]
    };

    // screen coordinates (int integer)
    screen_pos[0] = round(clip_pos[0]*(screenWidth /2.0)+(screenWidth /2.0));
    screen_pos[1] = round(clip_pos[1]*(screenHeight/2.0)+(screenHeight/2.0));

    // add panning
    screen_pos[0] += round(screenWidth  * panPercentage[0] * imageZoom);
    screen_pos[1] += round(screenHeight * panPercentage[1] * imageZoom); 

    // return point depth
    return clip_pos[2];
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
createColorPPM
(std::string filename, unsigned char *data, int width, int height){
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
writeOutputToFile( std::string filename, float * data, int dimX, int dimY )
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
writeOutputToFileByLine(std::string filename, float * data, int dimX, int dimY)
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

void writeDepthBufferToPPM(std::string filename , float * data, int dimX, int dimY)
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

void writeArrayToPPM(std::string filename , float * image, int dimX, int dimY)
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
            color[0] = std::max(std::min(image[index + 0]*alpha, 1.0f), 0.0f) * 255;  // red
            color[1] = std::max(std::min(image[index + 1]*alpha, 1.0f), 0.0f) * 255;  // green
            color[2] = std::max(std::min(image[index + 2]*alpha, 1.0f), 0.0f) * 255;  // blue
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
int convexHull::overlap(convexHull _hull) 
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
