/*****************************************************************************
*
* Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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

// ************************************************************************* //
//                            avtSLIVRVoxelExtractor.C                        //
// ************************************************************************* //

#include <avtSLIVRVoxelExtractor.h>

#include <float.h>

#include <avtAccessor.h>
#include <avtCellList.h>
#include <avtVolume.h>
#include <avtSLIVRImgMetaData.h>
#include <avtMemory.h>
#include <avtCallback.h>

#include <vtkDataArray.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkTemplateAliasMacro.h>
#include <vtkUnsignedCharArray.h>

#include <DebugStream.h>
#include <TimingsManager.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>
#include <math.h>

#if defined (_MSC_VER) && (_MSC_VER < 1800) && !defined(round)
inline double round(double x) {return (x-floor(x)) > 0.5 ? ceil(x) : floor(x);}
#endif

// ****************************************************************************
//  Method: Dot
//
//  Purpose:
//      Dot product for vectors
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//  Modifications:
//
// ****************************************************************************

float Dot (float vecA[3], float vecB[3])
{
    return ((vecA[0] * vecB[0]) + (vecA[1] * vecB[1]) + (vecA[2] * vecB[2]));
}


// ****************************************************************************
//  Method: Normalize
//
//  Purpose:
//      Normalize vector
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//  Modifications:
//
// ****************************************************************************

void Normalize(float vec[3])
{
    float inverse_sqrt_sum_squared = 
	sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (inverse_sqrt_sum_squared != 0)
	inverse_sqrt_sum_squared = 1.0/inverse_sqrt_sum_squared;
    for (int i=0;i<3; i++)
	vec[i] = vec[i]*inverse_sqrt_sum_squared;
}


// ****************************************************************************
//  Method: Reflect
//
//  Purpose:
//      Reflect vector
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//  Modifications:
//
// ****************************************************************************

void Reflect(float vec[3], float normal[3], float refl[3])
{
    // Vnew = -2*(V dot N)*N + V
    // out = incidentVec - 2.f * Dot(incidentVec, normal) * normal;
    float vDotN = Dot(vec,normal);
    for (int i=0; i<3; i++)
    {
        refl[i] = -2.0 * vDotN * normal[i] + vec[i];
    }
    Normalize(refl);
}


// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor constructor
//
//  Arguments:
//     w     The number of sample points in the x direction (width).
//     h     The number of sample points in the y direction (height).
//     d     The number of sample points in the z direction (depth).
//     vol   The volume to put samples into.
//     cl    The cell list to put cells whose sampling was deferred.
//
//  Programmer: Hank Childs
//  Creation:   December 14, 2003
//
//  Modifications:
//
//    Hank Childs, Fri Nov 19 14:50:58 PST 2004
//    Initialize gridsAreInWorldSpace.
//
//    Jeremy Meredith, Thu Feb 15 13:11:34 EST 2007
//    Added an ability to extract voxels using the world-space version
//    even when they're really in image space.
//
//    Hank Childs, Wed Aug 27 11:11:28 PDT 2008
//    Initialize spatial coordinates array.
//
//    Hank Childs, Wed Dec 24 11:22:43 PST 2008
//    Remove reference to ProportionSpaceToZBufferSpace data member.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

avtSLIVRVoxelExtractor::avtSLIVRVoxelExtractor(int w, int h, int d,
                                             avtVolume *vol, avtCellList *cl)
    : avtVoxelExtractor(w, h, d, vol, cl)
{
    fullImgWidth = w;
    fullImgHeight = h;

    debug5 << "fullImgWidth: "  << fullImgWidth << "  "
	   << "fullImgHeight: " << fullImgHeight << std::endl;

    rayCastingSLIVR = true;

    model_to_screen_transform = vtkMatrix4x4::New();
    screen_to_model_transform = vtkMatrix4x4::New();

    prop_buffer   = new double[3*depth];
    ind_buffer    = new int[3*depth];
    valid_sample  = new bool[depth];

    lighting = false;
    lightPosition[0] = lightPosition[1] = lightPosition[2] = 0.0;   
    lightPosition[3] = 1.0;
    materialProperties[0] = 0.4; materialProperties[1] = 0.75;
    materialProperties[3] = 0.0; materialProperties[3] = 15.0;
    gradient[0] = gradient[1] = gradient[2] = 0;

    proc = patch = 0;
    patchDrawn = 0;
    imgDims[0] = imgDims[1] = 0; // size of the patch
    // coordinates in the whole image
    imgLowerLeft[0]  = imgLowerLeft[1]  = 0;   
    imgUpperRight[0] = imgUpperRight[1] = 0;
    eyeSpaceDepth = -1;
    clipSpaceDepth = -1;
    imgArray = NULL;                         

    depthBuffer = NULL;
    rgbColorBuffer = NULL;

    // RayCasting:OSPRay renderer
    ospray = NULL;
}


// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor destructor
//
//  Purpose:
//      Defines the destructor.  Note: this should not be inlined in the header
//      because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   February 5, 2004
//
//  Modifications:
//
//    Hank Childs, Sun Nov 21 10:35:40 PST 2004
//    Delete the view to world transform.
//
//    Hank Childs, Wed Aug 27 11:10:51 PDT 2008
//    Delete the spatial coordinate arrays.
//
//    Hank Childs, Wed Dec 24 11:22:43 PST 2008
//    Remove reference to ProportionSpaceToZBufferSpace data member.
//
// ****************************************************************************

avtSLIVRVoxelExtractor::~avtSLIVRVoxelExtractor()
{
    model_to_screen_transform->Delete();
    screen_to_model_transform->Delete();

    if (prop_buffer != NULL)
        delete [] prop_buffer;
    if (ind_buffer != NULL)
        delete [] ind_buffer;
    if (valid_sample != NULL)
        delete [] valid_sample;
    if (imgArray != NULL)
        delete []imgArray;

    imgArray = NULL;
}

// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::Extract
//
//  Purpose:
//      Extracts the grid into the sample points.
//
//  Programmer: Hank Childs
//  Creation:   November 19, 2004
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 13:11:34 EST 2007
//    Added an ability to extract voxels using the world-space version
//    even when they're really in image space.
//
//    Hank Childs, Fri Jun  1 16:40:10 PDT 2007
//    Added support for non-scalars.
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::Extract(vtkRectilinearGrid *rgrid,
                std::vector<std::string> &varnames, std::vector<int> &varsizes)
{
    if (gridsAreInWorldSpace || pretendGridsAreInWorldSpace)
        ExtractWorldSpaceGridRCSLIVR(rgrid, varnames, varsizes);
    else
        ExtractImageSpaceGrid(rgrid, varnames, varsizes);
}

// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::SampleVariable
//
//  Purpose:
//      Actually samples the variable into our temporaray structure.
//
//  Programmer: Hank Childs
//  Creation:   November 20, 2004
//
//  Modifications:
//
//    Hank Childs, Fri Jun  1 15:45:58 PDT 2007
//    Add support for non-scalars.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************
void
avtSLIVRVoxelExtractor::SampleVariable(int first, int last, int w, int h)
{
    bool inrun = false;
    int  count = 0;

    avtRay *ray = volume->GetRay(w, h);
    bool calc_cell_index = ((ncell_arrays > 0) || (ghosts != NULL));

    double dest_rgb[4] = {0.0,0.0,0.0, 0.0};     // to store the computed color
    for (int i = first ; i < last ; i++)
    {
        const int *ind = ind_buffer + 3*i;
        const double *prop = prop_buffer + 3*i;

        int index = 0;
        if (calc_cell_index)
            index = ind[2]*((dims[0]-1)*(dims[1]-1)) +
		    ind[1]*(dims[0]-1) +
		    ind[0];

        if (ghosts != NULL)
        {
            if (ghosts[index] != 0)
               valid_sample[i] = false;
        }

        if (!valid_sample[i] && inrun)
        {
            ray->SetSamples(i-count, i-1, tmpSampleList);
            inrun = false;
            count = 0;
        }


        //
        // Trilinear SETUP
        int indices[6];
        float dist_from_left, 
	      dist_from_right,
	      dist_from_top,
	      dist_from_bottom,
	      dist_from_front, dist_from_back;
        float x_left, y_bottom, z_front;

        if (trilinearInterpolation)
        {
            int index_left,  index_right,
		index_top,   index_bottom, 
		index_front, index_back;

            int newInd[3];
            newInd[0] = ind[0];
            newInd[1] = ind[1];
            newInd[2] = ind[2];

            float x_right = prop[0];        x_left   = 1. - x_right;
            float y_top   = prop[1];        y_bottom = 1. - y_top;
            float z_back  = prop[2];        z_front  = 1. - z_back;

            // get the index and distance from the center of the neighbouring
	    // cells
            GetIndexandDistFromCenter(x_right, newInd[0], 
				      index_left, index_right,
				      dist_from_left, dist_from_right);
            GetIndexandDistFromCenter(y_top,   newInd[1], 
				      index_bottom,index_top,
				      dist_from_bottom,dist_from_top);
            GetIndexandDistFromCenter(z_back,  newInd[2],
				      index_front, index_back, 
				      dist_from_front, dist_from_back);


            indices[4] = index_front;       indices[5] = index_back;
            indices[2] = index_bottom;      indices[3] = index_top;
            indices[0] = index_left;        indices[1] = index_right;


            if (indices[0] < 0 || indices[0]>dims[0]-2)
                valid_sample[i] = false;

            if (indices[1] < 0 || indices[1]>dims[0]-2)
                valid_sample[i] = false;


            if (indices[2] < 0 || indices[2]>dims[1]-2)
                valid_sample[i] = false;

            if (indices[3] < 0 || indices[3]>dims[1]-2)
                valid_sample[i] = false;


            if (indices[4] < 0 || indices[4]>dims[2]-2)
                valid_sample[i] = false;

            if (indices[5] < 0 || indices[5]>dims[2]-2)
                valid_sample[i] = false;
        }

        if (!valid_sample[i])
            continue;

        //
        // Trilinear RUN
        if (trilinearInterpolation)
        {
            //
            // Cell centered data
            //
            if (ncell_arrays > 0){
                int indexT[8];
                ComputeIndices(dims, indices, indexT);

                for (int l = 0 ; l < ncell_arrays ; l++)            // ncell_arrays: usually 1
                {

                    void  *cell_array = cell_arrays[l];
                    double values[8];
                    for (int m = 0 ; m < cell_size[l] ; m++){       // cell_size[l] usually 1
                        AssignEight(cell_vartypes[l], values, indexT, cell_size[l], m, cell_array);
                        double scalarValue = TrilinearInterpolate(values, dist_from_left, dist_from_bottom, dist_from_front);

                        tmpSampleList[count][cell_index[l]+m] = scalarValue;
                    }
                }
            }

            //
            // Node centered data
            //
            if (npt_arrays > 0)
            {
                int indexT[8];
                ComputeIndicesVert(dims, indices, indexT);

                for (int l = 0 ; l < npt_arrays ; l++)
                {
                    void  *pt_array = pt_arrays[l];
                    double values[8];
                    for (int m = 0 ; m < pt_size[l] ; m++)
                    {
                        AssignEight(pt_vartypes[l], values, indexT, pt_size[l], m, pt_array);
                        double scalarValue = TrilinearInterpolate(values, x_left, y_bottom, z_front);

                        tmpSampleList[count][pt_index[l]+m] = scalarValue;
                    }
                }
            }
        }
        else
        {
            if (ncell_arrays > 0)
            {
                for (int l = 0 ; l < ncell_arrays ; l++)
                {
                    for (int m = 0 ; m < cell_size[l] ; m++)
                        tmpSampleList[count][cell_index[l]+m] =
                                         ConvertToDouble(cell_vartypes[l], index,
                                                      cell_size[l], m, cell_arrays[l]);
                }
            }

            if (npt_arrays > 0)
            {
                int index[8];
                index[0] = (ind[2])  *dims[0]*dims[1] + (ind[1])  *dims[0] + (ind[0]);
                index[1] = (ind[2])  *dims[0]*dims[1] + (ind[1])  *dims[0] + (ind[0]+1);
                index[2] = (ind[2])  *dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]);
                index[3] = (ind[2])  *dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]+1);
                index[4] = (ind[2]+1)*dims[0]*dims[1] + (ind[1])  *dims[0] + (ind[0]);
                index[5] = (ind[2]+1)*dims[0]*dims[1] + (ind[1])  *dims[0] + (ind[0]+1);
                index[6] = (ind[2]+1)*dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]);
                index[7] = (ind[2]+1)*dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]+1);

                double x_right = prop[0];
                double x_left = 1. - prop[0];
                double y_top = prop[1];
                double y_bottom = 1. - prop[1];
                double z_back = prop[2];
                double z_front = 1. - prop[2];
                for (int l = 0 ; l < npt_arrays ; l++)
                {
                    void  *pt_array = pt_arrays[l];
                    int    s = pt_size[l];
                    for (int m = 0 ; m < s ; m++)
                    {
                        double vals[8];
                        AssignEight(pt_vartypes[l], vals, index, s, m, pt_array);
                        double val =
                          x_left*y_bottom*z_front*vals[0] +
                          x_right*y_bottom*z_front*vals[1] +
                          x_left*y_top*z_front*vals[2] +
                          x_right*y_top*z_front*vals[3] +
                          x_left*y_bottom*z_back*vals[4] +
                          x_right*y_bottom*z_back*vals[5] +
                          x_left*y_top*z_back*vals[6] +
                          x_right*y_top*z_back*vals[7];

                        tmpSampleList[count][pt_index[l]+m] = val;
                    }
                }
            }
        }

        inrun = true;
        count++;
    }

    //
    // Make sure we get runs at the end.
    //
    if (inrun)
        ray->SetSamples(last-count, last-1, tmpSampleList);
}

// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::SampleAlongSegment
//
//  Purpose:
//      Samples the grid along a line segment.
//
//  Programmer: Hank Childs
//  Creation:   November 20, 2004
//
//  Modifications:
//
//    Hank Childs, Tue Jan  3 17:26:11 PST 2006
//    Fix bug that ultimately led to UMR where sampling occurred along
//    invalid values.
//
//    Hank Childs, Wed Dec 24 11:22:43 PST 2008
//    No longer use the ProportionSpaceToZBufferSpace data member, as we now
//    do our sampling in even intervals (wbuffer).
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::SampleAlongSegment(const double *origin,
                                          const double *terminus, int w, int h)
{
    int first = 0;
    int last = 0;
    bool hasIntersections = FindSegmentIntersections(origin, terminus,
                                                     first, last);

    if (!hasIntersections) { return; }

    //
    // Determine if there is intersection with buffer
    //
    double worldOpaqueCoordinates[3];
    bool   intersecting = false;
    int    intersect = -1;
    float  propAlongVector = 0;

    bool foundHit = false;
    int curX = -1;
    int curY = -1;
    int curZ = -1;
    bool xGoingUp = (terminus[0] > origin[0]);
    bool yGoingUp = (terminus[1] > origin[1]);
    bool zGoingUp = (terminus[2] > origin[2]);

    double x_dist = (terminus[0]-origin[0]);
    double y_dist = (terminus[1]-origin[1]);
    double z_dist = (terminus[2]-origin[2]);

    double pt[3];
    bool hasSamples = false;

    if (rayCastingSLIVR)
    {
	int screenX = bufferExtents[1] - bufferExtents[0];
	int screenY = bufferExtents[3] - bufferExtents[2];
	int screenIndex = h * screenX + w;
	if (screenIndex < 0 || screenIndex >= screenX * screenY) { return; }

        if (depthBuffer[screenIndex] != 1)  
        {
	    // There is some other things to blend with at this location ...
            // -- switching from (0 - 1) to (-1 - 1)
	    float normalizedDepth = depthBuffer[screenIndex]*2 - 1;  

	    debug5 << "normalizedDepth: " << normalizedDepth << " "
	    	   << "renderingDepthsExtents[0]: " 
		   << renderingDepthsExtents[0] << " " 
	    	   << "renderingDepthsExtents[1]: " 
		   << renderingDepthsExtents[1] << std::endl;

	    // ... and it's within this patch
	    if ((normalizedDepth >= renderingDepthsExtents[0]) && 
		(normalizedDepth <= renderingDepthsExtents[1]))  
	    {

		slivr::ProjectScreenToWorld(w, h, depthBuffer[screenIndex],
					    fullImgWidth, fullImgHeight, 
					    panPercentage, imageZoom, 
					    screen_to_model_transform,
					    worldOpaqueCoordinates);

		// debug5 << "Location: " << w << ", " << h << std::endl;
		// double start[3];
		// start[0] = terminus[0];
		// start[1] = terminus[1];
		// start[2] = terminus[2];
		// if (xGoingUp)
		//     start[0] = origin[0];
		// if (yGoingUp)
		//     start[1] = origin[1];
		// if (zGoingUp)
		//     start[1] = origin[2];

		double start[3] = {
		    xGoingUp ? origin[0] : terminus[0],
		    yGoingUp ? origin[1] : terminus[1],
		    zGoingUp ? origin[2] : terminus[2]
		};

		float distOriginTerminus_Squared = 
		    (origin[0]-terminus[0])*(origin[0]-terminus[0]) +
		    (origin[1]-terminus[1])*(origin[1]-terminus[1]) + 
		    (origin[2]-terminus[2])*(origin[2]-terminus[2]);
		float distCoordStart_Squared = 
		    (worldOpaqueCoordinates[0]-start[0])*(worldOpaqueCoordinates[0]-start[0]) +
		    (worldOpaqueCoordinates[1]-start[1])*(worldOpaqueCoordinates[1]-start[1]) + 
		    (worldOpaqueCoordinates[2]-start[2])*(worldOpaqueCoordinates[2]-start[2]);

		// lies along the vector
		if (distCoordStart_Squared < distOriginTerminus_Squared)    
		{
		    intersecting = true;
		    propAlongVector = 
			(float)sqrt(distCoordStart_Squared)/sqrt(distOriginTerminus_Squared);
		}
	    }
	}
    }

    for (int i = first ; i < last ; i++)
    {
        int *ind = ind_buffer + 3*i;
        double *dProp = prop_buffer + 3*i;
        valid_sample[i] = false;

        double proportion = ((double)i)/((double)depth);
        pt[0] = origin[0] + proportion*x_dist;
        pt[1] = origin[1] + proportion*y_dist;
        pt[2] = origin[2] + proportion*z_dist;

        ind[0] = -1;
        ind[1] = -1;
        ind[2] = -1;

        if (!foundHit)
        {
            //
            // We haven't found any hits previously.  Exhaustively search
            // through arrays and try to find a hit.
            //
            ind[0] = FindMatch(X, pt[0], dims[0]);
            if (ind[0] >= 0)
                dProp[0] = (pt[0] - X[ind[0]]) * divisors_X[ind[0]];
            ind[1] = FindMatch(Y, pt[1], dims[1]);
            if (ind[1] >= 0)
                dProp[1] = (pt[1] - Y[ind[1]]) * divisors_Y[ind[1]];
            ind[2] = FindMatch(Z, pt[2], dims[2]);
            if (ind[2] >= 0)
                dProp[2] = (pt[2] - Z[ind[2]]) * divisors_Z[ind[2]];
        }
        else
        {
            //
            // We have found a hit before.  Try to locate the next sample
            // based on what we already found.
            //
            if (xGoingUp)
            {
                for ( ; curX < dims[0]-1 ; curX++)
                {
                    if (pt[0] >= X[curX] && pt[0] <= X[curX+1])
                    {
                        dProp[0] = (pt[0] - X[curX]) * divisors_X[curX];
                        ind[0] = curX;
                        break;
                    }
                }
            }
            else
            {
                for ( ; curX >= 0 ; curX--)
                {
                    if (pt[0] >= X[curX] && pt[0] <= X[curX+1])
                    {
                        dProp[0] = (pt[0] - X[curX]) * divisors_X[curX];
                        ind[0] = curX;
                        break;
                    }
                }
            }
            if (yGoingUp)
            {
                for ( ; curY < dims[1]-1 ; curY++)
                {
                    if (pt[1] >= Y[curY] && pt[1] <= Y[curY+1])
                    {
                        dProp[1] = (pt[1] - Y[curY]) * divisors_Y[curY];
                        ind[1] = curY;
                        break;
                    }
                }
            }
            else
            {
                for ( ; curY >= 0 ; curY--)
                {
                    if (pt[1] >= Y[curY] && pt[1] <= Y[curY+1])
                    {
                        dProp[1] = (pt[1] - Y[curY]) * divisors_Y[curY];
                        ind[1] = curY;
                        break;
                    }
                }
            }
            if (zGoingUp)
            {
                for ( ; curZ < dims[2]-1 ; curZ++)
                {
                    if (pt[2] >= Z[curZ] && pt[2] <= Z[curZ+1])
                    {
                        dProp[2] = (pt[2] - Z[curZ]) * divisors_Z[curZ];
                        ind[2] = curZ;
                        break;
                    }
                }
            }
            else
            {
                for ( ; curZ >= 0 ; curZ--)
                {
                    if (pt[2] >= Z[curZ] && pt[2] <= Z[curZ+1])
                    {
                        dProp[2] = (pt[2] - Z[curZ]) * divisors_Z[curZ];
                        ind[2] = curZ;
                        break;
                    }
                }
            }
        }

        bool intersectedDataset = !(ind[0] < 0 || ind[1] < 0 || ind[2] < 0);
        if (!intersectedDataset)
        {
            if (!foundHit)
            {
                // We still haven't found the start.  Keep looking.
                continue;
            }
            else
            {
                // This is the true terminus.
                last = i;
                break;
            }
        }
        else  // Did intersect data set.
        {
            if (!foundHit)
            {
                // This is the first true sample.  "The true start"
                first = i;
            }
        }

        valid_sample[i] = true;
        foundHit = true;
        hasSamples = true;

        curX = ind[0];
        curY = ind[1];
        curZ = ind[2];
    }

    debug5 << "First: " << first << "  last: " << last << std::endl;
    if (intersecting){
	intersect = floor(propAlongVector * (last-first) + first);
	debug5 << "intersect: " << intersect
	       << " first: " << first << " last: " << last << std::endl;
    }

    if (hasSamples) {
	if (rayCastingSLIVR) {
	    SampleVariableRCSLIVR(first, last, intersect, w, h);
	}
	else {
	    SampleVariable(first, last, w, h);
	}
    }
}


// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::ExtractWorldSpaceGridRCSLIVR
//
//  Purpose:
//      Compute region that patch covers
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//  Modifications:
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::ExtractWorldSpaceGridRCSLIVR(vtkRectilinearGrid *rgrid,
                 std::vector<std::string> &varnames, std::vector<int> &varsize)
{
    int timing_ExtractWorldSpaceGridRCSLIVR = visitTimer->StartTimer();
    //=======================================================================//
    // Initialization
    //=======================================================================//
    // Flag to indicate if the patch is drawn
    patchDrawn = 0;
    
    //=======================================================================//
    // Register data and early skipping
    //=======================================================================//
    int timing_register_data = visitTimer->StartTimer();
    // Some of our sampling routines need a chance to pre-process the data.
    // Register the grid here so we can do that.
    // Stores the values in a structure so that it can be used
    RegisterGrid(rgrid, varnames, varsize);
    // Determine what range we are dealing with on this iteration.
    int w_min = restrictedMinWidth;
    int w_max = restrictedMaxWidth + 1;
    int h_min = restrictedMinHeight;
    int h_max = restrictedMaxHeight + 1;
    imgWidth = imgHeight = 0;
    // Let's find out if this range can even intersect the dataset.
    // If not, just skip it.
    if (!FrustumIntersectsGrid(w_min, w_max, h_min, h_max)) { return; }
    // Timing
    visitTimer->StopTimer(timing_register_data, 
			  "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			  "Register Data (VisIt preparation)");

    //=======================================================================//
    // obtain data pointers & ghost region information
    //=======================================================================//
    int timing_get_metadata = visitTimer->StartTimer();
    // Calculate patch dimensions for point array and cell array
    //   This is to check if the patch is a cell data or a point data
    //   I have to assume cell dataset has a higher priority
    void* volumePointer = NULL;
    int   volumeDataType;
    int nX = 0, nY = 0, nZ = 0;
    if (ncell_arrays > 0) {
	ospout << "[avtMassVoxelExtractor] Cell Dataset " << std::endl;
	if (ncell_arrays != 1 || cell_size[0] != 1) {
	    EXCEPTION1(VisItException, 
		       "Trying to plot more than one field, "
		       "which is not supported by OSPRay SLIVR. "
		       "Use other render type instead");
	}
	nX = dims[0] - 1;
	nY = dims[1] - 1;
	nZ = dims[2] - 1;
	volumePointer = cell_arrays[0];
	volumeDataType = cell_vartypes[0];
    }
    else if (npt_arrays > 0) {
	ospout << "[avtMassVoxelExtractor] Point Dataset " << std::endl;
	if (npt_arrays != 1 || pt_size[0] != 1) {
	    EXCEPTION1(VisItException, 
		       "Trying to plot more than one field, "
		       "which is not supported by OSPRay SLIVR. "
		       "Use other render type instead");
	}
	nX = dims[0];
	nY = dims[1];
	nZ = dims[2];
	volumePointer = pt_arrays[0];
	volumeDataType = pt_vartypes[0];	
    } else {
	std::cerr << "WARNING: Empty dataset " << std::endl;
    }
    ospout << "[avtMassVoxelExtractor] patch dimension "
	   << nX << " " << nY << " " << nZ << std::endl;
    // Calculate ghost region boundaries
    //   ghost_boundaries is an array to indicate if the patch contains
    //   any ghost regions in six different directions
    // Here I assume the patch is larger than 3-cube
    // If not then you might want to dig into this code and see if
    // there will be any special boundary cases
    //
    // debug5 << "VAR: ghost value " << (int)ghosts[0] << std::endl;
    //
    bool ghost_bound[6] = {false};
    if (ghosts != NULL)
    {
	int gnX = 0, gnY = 0, gnZ = 0;
	gnX = dims[0] - 1;
	gnY = dims[1] - 1;
        gnZ = dims[2] - 1;

	// debug the meaning of the ghost zoom
	// for (int z = 0; z < gnZ; ++z) {
	//     for (int y = 0; y < gnY; ++y) {
	// 	for (int x = 0; x < gnX; ++x) {
	// 	    std::cout << (int)ghosts[z*gnY*gnX+y*gnX+x] << " ";
	// 	}
	// 	std::cout << std::endl;
	//     }
	//     std::cout << std::endl << std::endl;
	// }
	
	for (int y = 1; y < (gnY-1); ++y) {
	    for (int z = 1; z < (gnZ-1); ++z) {
		if (!ghost_bound[0]) {
		    if (ghosts[z*gnY*gnX+y*gnX        ] != 0)
		    { ghost_bound[0] = true; }
		}
		if (!ghost_bound[3]) {
		    if (ghosts[z*gnY*gnX+y*gnX+(gnX-1)] != 0)
		    { ghost_bound[3] = true; }
		}
		if (ghost_bound[0] && ghost_bound[3]) { break; }
	    }
	}
	for (int x = 1; x < (gnX-1); ++x) {
	    for (int z = 1; z < (gnZ-1); ++z) {
		if (!ghost_bound[1]) {
		    if (ghosts[z*gnY*gnX            +x] != 0)
		    { ghost_bound[1] = true; }
		}
		if (!ghost_bound[4]) {
		    if (ghosts[z*gnY*gnX+(gnY-1)*gnX+x] != 0)
		    { ghost_bound[4] = true; }
		}
		if (ghost_bound[1] && ghost_bound[4]) { break; }
	    }
	}
	for (int x = 1; x < (gnX-1); ++x) {
	    for (int y = 1; y < (gnY-1); ++y) {
		if (!ghost_bound[2]) {
		    if (ghosts[                y*gnX+x] != 0) 
		    { ghost_bound[2] = true; }
		}
		if (!ghost_bound[5]) {
		    if (ghosts[(gnZ-1)*gnY*gnX+y*gnX+x] != 0)
		    { ghost_bound[5] = true; }
		}
		if (ghost_bound[2] && ghost_bound[5]) { break; }
	    }
	}
    }
    // Data bounding box
    double volumeCube[6] = {
	X[0], X[nX-1],
	Y[0], Y[nY-1],
	Z[0], Z[nZ-1]
    };
    // Timing
    visitTimer->StopTimer(timing_get_metadata , 
			  "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			  "Compute metadata & ghost boundary (Pre-OSPRay preparation)");

    //=======================================================================//
    // Determine the screen size of the patch being processed
    //=======================================================================//
    int timing_get_screen_projection = visitTimer->StartTimer();
    int patchScreenExtents[4];
    slivr::ProjectWorldToScreenCube(volumeCube, w_max, h_max, 
				    panPercentage, imageZoom,
				    model_to_screen_transform, 
				    patchScreenExtents, 
				    renderingDepthsExtents);
    xMin = patchScreenExtents[0];
    xMax = patchScreenExtents[1];
    yMin = patchScreenExtents[2];
    yMax = patchScreenExtents[3];

    ospout << "[avtMassVoxelExtractor] patch ghost bounds:"
	   << "   " << ghost_bound[0] << " " << ghost_bound[3] 
	   << " | " << ghost_bound[1] << " " << ghost_bound[4] 
	   << " | " << ghost_bound[2] << " " << ghost_bound[5]
	   << std::endl;   
    
    // calculate patch depth
    double patch_center[3];
    patch_center[0] = (volumeCube[0] + volumeCube[1])/2.0;
    patch_center[1] = (volumeCube[2] + volumeCube[3])/2.0;
    patch_center[2] = (volumeCube[4] + volumeCube[5])/2.0;        
    double patch_depth = // use the norm of patch center as patch depth
	std::sqrt((patch_center[0]-view.camera[0])*
		  (patch_center[0]-view.camera[0])+
		  (patch_center[1]-view.camera[1])*
		  (patch_center[1]-view.camera[1])+
		  (patch_center[2]-view.camera[2])*
		  (patch_center[2]-view.camera[2]));
    eyeSpaceDepth = patch_depth;
    clipSpaceDepth = renderingDepthsExtents[0];
    // Timing
    visitTimer->StopTimer(timing_get_screen_projection, 
			  "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			  "Get screen size of the patch (Pre-OSPRay preparation)");

    //=======================================================================//
    // create framebuffer
    //=======================================================================//
    int timing_create_imgarray = visitTimer->StartTimer();
    // assign data to the class
    //xMax+=1; yMax+=1;
    ospout << "[avtMassVoxelExtractor] patch extents " 
	   << xMin << " " << xMax << " "
	   << yMin << " " << yMax << std::endl;
    if (xMin < fullImageExtents[0]) { xMin = fullImageExtents[0]; }
    if (yMin < fullImageExtents[2]) { yMin = fullImageExtents[2]; }    
    if (xMax > fullImageExtents[1]) { xMax = fullImageExtents[1]; }
    if (yMax > fullImageExtents[3]) { yMax = fullImageExtents[3]; }
    imgWidth  = xMax-xMin;
    imgHeight = yMax-yMin;

    // Initialize memory (framebuffer)
    if (avtCallback::UseOSPRay()) {
	// framebuffer
	imgArray = new float[((imgWidth)*4) * imgHeight];   
    } else {
        // framebuffer initialized
	imgArray = new float[((imgWidth)*4) * imgHeight](); 
    };
    // Timing
    visitTimer->StopTimer(timing_create_imgarray, 
			  "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			  "Create ImgArray (Pre-OSPRay preparation)");

    //=======================================================================//
    // Render using OSPRay
    //=======================================================================//
    int timing_using_ospray = visitTimer->StartTimer();
    if (avtCallback::UseOSPRay()) {
	int timing_setup_ospray = visitTimer->StartTimer();
	if (!((npt_arrays == 1)^(ncell_arrays == 1)))
	{
	    std::cerr << "WARNING: Multiple data found within one patch, " 
		      << " We don't know what to do !! " 
		      << std::endl
		      << "         One of the dataset might be missing "
		      << std::endl;
	}
	// shift grid and make it cel centered for cell data
	double volumePBox[6] = {
	    // for cell centered data, we put the voxel on its left boundary
	    X[0], Y[0], Z[0], 
	    X[nX-1], Y[nY-1], Z[nZ-1]
	};
	// compute boundingbox and clipping plane for ospray
	double volumeBBox[6];
	if (ncell_arrays > 0) {
	    volumeBBox[0] = ghost_bound[0] ? (X[0]+X[1])/2. : volumePBox[0];
	    volumeBBox[1] = ghost_bound[1] ? (Y[0]+Y[1])/2. : volumePBox[1];
	    volumeBBox[2] = ghost_bound[2] ? (Z[0]+Z[1])/2. : volumePBox[2];
	    volumeBBox[3] = 
		ghost_bound[3] ? (X[nX-1]+X[nX-2])/2. : volumePBox[3];
	    volumeBBox[4] = 
		ghost_bound[4] ? (Y[nY-1]+Y[nY-2])/2. : volumePBox[4];
	    volumeBBox[5] = 
		ghost_bound[5] ? (Z[nZ-1]+Z[nZ-2])/2. : volumePBox[5];
	}
	else {
	    volumeBBox[0] = ghost_bound[0] ? X[1] : volumePBox[0];
	    volumeBBox[1] = ghost_bound[1] ? Y[1] : volumePBox[1];
	    volumeBBox[2] = ghost_bound[2] ? Z[1] : volumePBox[2];
	    volumeBBox[3] = ghost_bound[3] ? X[nX-2] : volumePBox[3];
	    volumeBBox[4] = ghost_bound[4] ? Y[nY-2] : volumePBox[4];
	    volumeBBox[5] = ghost_bound[5] ? Z[nZ-2] : volumePBox[5];
	}
	ospout << "[avtMassVoxelExtractor] patch data position:" 
	       << " " << volumePBox[0]
	       << " " << volumePBox[1]
	       << " " << volumePBox[2]
	       << " |"
	       << " " << volumePBox[3]
	       << " " << volumePBox[4]
	       << " " << volumePBox[5]
	       << std::endl;  
	ospout << "[avtMassVoxelExtractor] patch data bbox:" 
	       << " " << volumeBBox[0]
	       << " " << volumeBBox[1]
	       << " " << volumeBBox[2]
	       << " |"
	       << " " << volumeBBox[3]
	       << " " << volumeBBox[4]
	       << " " << volumeBBox[5]
	       << std::endl; 
	// Timing
	visitTimer->StopTimer(timing_setup_ospray, 
			      "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			      "OSPRay bbox and clip (OSPRay preparation)");
	
	// Create volume and model
	int timing_create_volume = visitTimer->StartTimer();
	ospray->GetPatch(patch)->Set(volumeDataType, volumePointer,
				     X, Y, Z, nX, nY, nZ, volumePBox, volumeBBox, 
				     materialProperties, (float)rendererSampleRate,
				     lighting);
	visitTimer->StopTimer(timing_create_volume, 
			      "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			      "OSPRay Create Volume");

	// Render Volume
	int timing_render_volume = visitTimer->StartTimer();
	if ((scalarRange[1] >= tFVisibleRange[0]) &&
	    (scalarRange[0] <= tFVisibleRange[1]))
	{
	    ospray->Render(xMin, xMax, yMin, yMax,
			   imgWidth, imgHeight, imgArray,
			   ospray->GetPatch(patch));
	    patchDrawn = 1;

	}
	visitTimer->StopTimer(timing_render_volume, 
			      "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			      "OSPRay Render Volume");	
    }
    // Timing
    visitTimer->StopTimer(timing_using_ospray, 
			  "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			  "Using OSPRay");

    //=======================================================================//
    // Send rays
    //=======================================================================//
    int timing_using_pascal = visitTimer->StartTimer();
    imgDims[0] = imgWidth;
    imgDims[1] = imgHeight;
    imgLowerLeft[0] = xMin;
    imgLowerLeft[1] = yMin;
    imgUpperRight[0] = xMax; 
    imgUpperRight[1] = yMax;

    if (!avtCallback::UseOSPRay()) {
	ospout << "[avtMassiveVoxelExtractor] "
	       << "Using CPU version raytracer" << std::endl;
    	for (int patchX = xMin; patchX < xMax ; patchX++) {
	    for (int patchY = yMin; patchY < yMax ; patchY++) {

		const int pIndex = (patchY-yMin)*imgWidth + (patchX-xMin);
		const int fIndex = ((patchY-bufferExtents[2])*
				    (bufferExtents[1]-bufferExtents[0])+
				    (patchX-bufferExtents[0]));
		// outside visible range
		if ((scalarRange[1] < tFVisibleRange[0]) ||
		    (scalarRange[0] > tFVisibleRange[1]))	
		{
		    if (depthBuffer[fIndex] != 1)
		    {
			const double clipDepth = depthBuffer[fIndex]*2 - 1;
			if (clipDepth >= renderingDepthsExtents[0] && 
			    clipDepth <= renderingDepthsExtents[1])
			{
			    patchDrawn = 1;
			    const int imgID = pIndex * 4;
			    imgArray[imgID + 0] = 
				rgbColorBuffer[fIndex*3 + 0]/255.0;
			    imgArray[imgID + 1] = 
				rgbColorBuffer[fIndex*3 + 1]/255.0;
			    imgArray[imgID + 2] = 
				rgbColorBuffer[fIndex*3 + 2]/255.0;
			    imgArray[imgID + 3] = 1.0;
			}
		    }
		}
		else
		{
		    // std::cout << "draw" << std::endl;
		    patchDrawn = 1;		    
		    // starting point where we start sampling
		    double origin[4]   = {0,0,0,1};
		    // ending point where we stop sampling 
		    double terminus[4] = {0,0,0,1};
		    // find the starting point & ending point of the ray
		    GetSegmentRCSLIVR(patchX, patchY, fullVolumeDepthExtents,
				      origin, terminus); 
		    // Go get the segments along this ray and store them in
		    SampleAlongSegment(origin, terminus, patchX, patchY);
		}
	    }
	}
    }
    // Timing
    visitTimer->StopTimer(timing_using_pascal, 
			  "avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR "
			  "Using Pascal (Post OSPRay)");

    //=======================================================================//
    // Deallocate memory if not used
    //=======================================================================//
    if (patchDrawn == 0)
    { 
	if (imgArray != NULL) 
	{ 
	    delete []imgArray; imgArray = NULL; 
	} 
    } 
    // else {
    // 	WriteArrayToPPM("patch-after-render"+
    // 			std::to_string(proc), imgArray, 
    // 			imgWidth, imgHeight);
    // }
    visitTimer->StopTimer(timing_ExtractWorldSpaceGridRCSLIVR, "Calling avtMassVoxelExtractor::ExtractWorldSpaceGridRCSLIVR");
}


// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::GetSegmentRCSLIVR
//
//  Purpose:
//      Determine ray start and end
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//  Modifications:
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::GetSegmentRCSLIVR(int x, int y, double depthsExtents[2],
					  double *origin, double *terminus)
{
   slivr::ProjectScreenToWorld(x, y, depthsExtents[0], 
				fullImgWidth, fullImgHeight, 
				panPercentage, imageZoom, 
				screen_to_model_transform,
				origin);
    slivr::ProjectScreenToWorld(x, y, depthsExtents[1], 
				fullImgWidth, fullImgHeight, 
				panPercentage, imageZoom, 
				screen_to_model_transform,
				terminus);
}




// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::SampleVariableRCSLIVR
//
//  Purpose:
//      Sample each ray
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::SampleVariableRCSLIVR(int first, int last, 
					      int intersect, int x, int y)
{
    int  count = 0;
    bool calc_cell_index = ((ncell_arrays > 0) || (ghosts != NULL));

    double dest_rgb[4] = {0.0,0.0,0.0, 0.0};     // to store the computed color
    for (int i = first ; i < last ; i++)
    {
        // If we intersect a value in the z buffer
        if (i == intersect)
        {
            int fullIndex = 
		(y * (bufferExtents[1]-bufferExtents[0]) + x) * 3.0;

            float bufferColor[4];
            bufferColor[0] = rgbColorBuffer[fullIndex + 0] / 255.0;
            bufferColor[1] = rgbColorBuffer[fullIndex + 1] / 255.0;
            bufferColor[2] = rgbColorBuffer[fullIndex + 2] / 255.0;
            bufferColor[3] = 1.0;

            for (int j=0; j<4; j++) {
                dest_rgb[j] = 
		    bufferColor[j] * (1.0 - dest_rgb[3]) + dest_rgb[j];
	    }

 	    // debug5 << x << ", " << y 
	    // 	   << " ~ First: " << first 
	    // 	   << " i:  " << i 
	    // 	   << " intersect: " << intersect 
	    // 	   << " bufferColor: " 
	    // 	   << bufferColor[0] << ", " 
	    // 	   << bufferColor[1] << ", " 
	    // 	   << bufferColor[2] 
	    // 	   << " dest_rgb: " 
	    // 	   << dest_rgb[0] << ", " 
	    // 	   << dest_rgb[1] << ", " 
	    // 	   << dest_rgb[2] << ", " 
	    // 	   << dest_rgb[3] << std::endl;
	    break;
       }

        const int *ind = ind_buffer + 3*i;
        const double *prop = prop_buffer + 3*i;

        int index = 0;
        if (calc_cell_index)
            index = ind[2]*((dims[0]-1)*(dims[1]-1)) +
		    ind[1]*(dims[0]-1) + 
		    ind[0];

        if (ghosts != NULL)
        {
            if (ghosts[index] != 0)
               valid_sample[i] = false;
        }

        int index_left,  index_right,
            index_top,   index_bottom,
	    index_front, index_back;
        float dist_from_left, dist_from_right,
	    dist_from_top,    dist_from_bottom,
	    dist_from_front,  dist_from_back;

        int newInd[3];
        newInd[0] = ind[0];
        newInd[1] = ind[1];
        newInd[2] = ind[2];

        float x_right = prop[0];        float x_left   = 1. - x_right;
        float y_top   = prop[1];        float y_bottom = 1. - y_top;
        float z_back  = prop[2];        float z_front  = 1. - z_back;       

        // get the index and distance from the center of the neighbouring cells
        GetIndexandDistFromCenter(x_right, newInd[0], index_left, index_right,
				  dist_from_left, dist_from_right);
        GetIndexandDistFromCenter(y_top,   newInd[1], index_bottom,index_top, 
				  dist_from_bottom,dist_from_top);
        GetIndexandDistFromCenter(z_back,  newInd[2], index_front, index_back,
				  dist_from_front, dist_from_back);

        int indices[6];
        indices[4] = index_front;       indices[5] = index_back;
        indices[2] = index_bottom;      indices[3] = index_top;
        indices[0] = index_left;        indices[1] = index_right;


        if (indices[0] < 0 || indices[0]>dims[0]-2)
            valid_sample[i] = false;

        if (indices[1] < 0 || indices[1]>dims[0]-2)
            valid_sample[i] = false;


        if (indices[2] < 0 || indices[2]>dims[1]-2)
            valid_sample[i] = false;

        if (indices[3] < 0 || indices[3]>dims[1]-2)
            valid_sample[i] = false;


        if (indices[4] < 0 || indices[4]>dims[2]-2)
            valid_sample[i] = false;

        if (indices[5] < 0 || indices[5]>dims[2]-2)
            valid_sample[i] = false;


        if (!valid_sample[i])
            continue;


        //
        // Cell centered data
        //
        if (ncell_arrays > 0)
        {
            int indexT[8];
            ComputeIndices(dims, indices, indexT);

            // ncell_arrays: usually 1
            for (int l = 0 ; l < ncell_arrays ; l++)
            {
                void  *cell_array = cell_arrays[l];
                double values[8];
		// cell_size[l] usually 1
                for (int m = 0 ; m < cell_size[l] ; m++)
                {
                    AssignEight(cell_vartypes[l], values, indexT,
				cell_size[l], m, cell_array);
                    double scalarValue = 
			TrilinearInterpolate(values, 
					     dist_from_left, 
					     dist_from_bottom,
					     dist_from_front);
                    double source_rgb[4];
                    int retVal = transferFn1D->QueryTF(scalarValue,source_rgb);

                    if (((retVal == 0) || (source_rgb[3]==0)) || 
			(source_rgb[0]==0 && 
			 source_rgb[1]==0 && 
			 source_rgb[2]==0))
		    {
                        // no need to do anything more if there will be no
			// color
                    }
                    else
                    {
                        //
                        // Compute Lighting (if needed)
                        //
                        if (lighting == true)
                        {
                            double vals[6];

                            // h = offset = 1/2 the distance between grids
                            // grad = 1/2*h * ( f(x+h,y,z)-f(x-h,y,z) 
			    // f(x,y+h,z)-f(x,y-h,z)  
			    // f(x,y,z-h)-f(x,y,z-h)  )
			    float distFromRight, distFromLeft, 
				  distFromTop, distFromBottom, 
				  distFromFront, distFromBack;
			    int indexLeft, indexRight, indexTop, 
				indexBottom, indexFront, indexBack;
			    float gradientOffset = 0.25;

			    double gradVals[8];
			    int indexGrad[8], gradInd[3], gradIndices[6];
			    float xRight, yTop, zBack;

			    void  *cell_array = cell_arrays[0];

                            //
                            // X
                            //
                            for (int i=0; i<6; i++)
                                gradIndices[i] = indices[i];

                            //
                            // find x-h
                            //
                            if (x_right - gradientOffset < 0.0){
                                xRight = (x_right - gradientOffset)+1.0;
                                gradInd[0] = ind[0]-1;
                            }
                            else{
                                xRight = x_right - gradientOffset;
                                gradInd[0] = ind[0];
                            }

                            GetIndexandDistFromCenter(xRight, gradInd[0], 
						      indexLeft, indexRight,
						  distFromLeft, distFromRight);
                            gradIndices[0] = indexLeft;
			    gradIndices[1] = indexRight;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(cell_vartypes[0], gradVals,
					indexGrad, 1, 0, cell_array);
                            vals[0] = TrilinearInterpolate(gradVals, 
							   distFromLeft,
							   dist_from_bottom,
							   dist_from_front);

                            //
                            // find x+h
                            //
                            if (x_right + gradientOffset > 1.0){
                                xRight = (x_right + gradientOffset)-1.0;
                                gradInd[0] = ind[0]+1;
                            }else{
                                xRight = x_right + gradientOffset;
                                gradInd[0] = ind[0];
                            }

                            GetIndexandDistFromCenter(xRight, gradInd[0],  indexLeft, indexRight,  distFromLeft, distFromRight);
                            gradIndices[0] = indexLeft;    gradIndices[1] = indexRight;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                            vals[1] = TrilinearInterpolate(gradVals, distFromLeft, dist_from_bottom, dist_from_front);



                            //
                            // Y
                            //
                            for (int i=0; i<6; i++)
                                gradIndices[i] = indices[i];

                            //
                            // find y-h
                            //
                            if (y_top - gradientOffset < 0.0){
                                yTop = (y_top - gradientOffset)+1.0;
                                gradInd[1] = ind[1]-1;
                            }
                            else{
                                yTop = y_top - gradientOffset;
                                gradInd[1] = ind[1];
                            }

                            GetIndexandDistFromCenter(yTop, gradInd[1],  indexBottom, indexTop,  distFromBottom, distFromTop);
                            gradIndices[2] = indexBottom ;    gradIndices[3] = indexTop;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                            vals[2] = TrilinearInterpolate(gradVals, dist_from_left, distFromBottom, dist_from_front);

                            //
                            // find y+h
                            //
                            yTop = y_top;
                            if (y_top + gradientOffset > 1.0){
                                yTop = (y_top + gradientOffset)-1.0;
                                gradInd[1] = ind[1]+1;
                            }else{
                                yTop = y_top + gradientOffset;
                                gradInd[1] = ind[1];
                            }

                            GetIndexandDistFromCenter(yTop, gradInd[1],  indexBottom, indexTop,  distFromBottom, distFromTop);
                            gradIndices[2] = indexBottom;    gradIndices[3] = indexTop;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                            vals[3] = TrilinearInterpolate(gradVals, dist_from_left, distFromBottom, dist_from_front);


                            //
                            // Z
                            //
                            for (int i=0; i<6; i++)
                                gradIndices[i] = indices[i];

                            //
                            // z-h
                            //
                            if (z_back - gradientOffset < 0.0){
                                zBack = (z_back - gradientOffset)+1.0;
                                gradInd[2] = ind[2]-1;
                            }
                            else{
                                zBack = z_back - gradientOffset;
                                gradInd[2] = ind[2];
                            }

                            GetIndexandDistFromCenter(zBack, gradInd[2],  indexFront, indexBack,  distFromFront, distFromBack);
                            gradIndices[4] = indexFront;    gradIndices[5] = indexBack;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                            vals[4] = TrilinearInterpolate(gradVals, dist_from_left, dist_from_bottom, distFromFront);

                            //
                            // z+h
                            //
                            if (z_back + gradientOffset > 1.0){
                                zBack = (z_back + gradientOffset)-1.0;
                                gradInd[2] = ind[2]+1;
                            }else{
                                zBack = z_back + gradientOffset;
                                gradInd[2] = ind[2];
                            }

                            GetIndexandDistFromCenter(zBack, gradInd[2],  indexFront, indexBack,  distFromFront, distFromBack);
                            gradIndices[4] = indexFront;    gradIndices[5] = indexBack;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                            vals[5] = TrilinearInterpolate(gradVals, dist_from_left, dist_from_bottom, distFromFront);


                            gradient[0] = (1.0/(2.0*gradientOffset)) * (vals[1] - vals[0]);
                            gradient[1] = (1.0/(2.0*gradientOffset)) * (vals[3] - vals[2]);
                            gradient[2] = (1.0/(2.0*gradientOffset)) * (vals[5] - vals[4]);

                            Normalize(gradient);
                        }

                        //
                        // Compute the color
                        //
                        ComputePixelColor(source_rgb, dest_rgb, gradient);
                    }

                }
            }
        }

        //
        // Node centered data
        //
        if (npt_arrays > 0)
        {
            int indexT[8];
            ComputeIndicesVert(dims, indices, indexT);

            for (int l = 0 ; l < npt_arrays ; l++)
            {
                void  *pt_array = pt_arrays[l];
                double values[8];
                for (int m = 0 ; m < pt_size[l] ; m++)
                {
                    AssignEight(pt_vartypes[l], values, indexT, pt_size[l], m, pt_array);
                    double scalarValue = TrilinearInterpolate(values, x_left, y_bottom, z_front);
                    double source_rgb[4];
                    int retVal = transferFn1D->QueryTF(scalarValue,source_rgb);
                    if ( ((retVal == 0)||(source_rgb[3]==0)) || (source_rgb[0]==0 && source_rgb[1]==0 && source_rgb[2]==0) )
                    {
                        // no need to do anything more if there will be no color
                    }
                    else
                    {
                        //
                        // Compute Lighting (if needed)
                        //
                        if (lighting == true)
                        {
                            double vals[6];

                            // h = offset = 1/2 the distance between grids
                            // grad = 1/2*h * ( f(x+h,y,z)-f(x-h,y,z)    f(x,y+h,z)-f(x,y-h,z)   f(x,y,z-h)-f(x,y,z-h)  )

                            float distFromRight, distFromLeft, distFromTop, distFromBottom, distFromFront, distFromBack;
                            int indexLeft, indexRight, indexTop, indexBottom, indexFront, indexBack;
                            float gradientOffset = 0.5;

                            double gradVals[8];
                            int indexGrad[8], gradInd[3], gradIndices[6];
                            float xRight, yTop, zBack = 0;



                            //
                            // X
                            //
                            for (int i=0; i<6; i++)
                                gradIndices[i] = indices[i];

                            //
                            // find x-h
                            //
                            if (x_right - gradientOffset < 0.0){
                                xRight = (x_right - gradientOffset)+1.0;
                                gradInd[0] = ind[0]-1;
                            }
                            else{
                                xRight = x_right - gradientOffset;
                                gradInd[0] = ind[0];
                            }

                            GetIndexandDistFromCenter(xRight, gradInd[0],  indexLeft, indexRight,  distFromLeft, distFromRight);
                            gradIndices[0] = indexLeft;    gradIndices[1] = indexRight;
                            ComputeIndicesVert(dims, gradIndices, indexGrad);
                            AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                            vals[0] = TrilinearInterpolate(gradVals, x_left-gradientOffset, y_bottom, z_front);

                            //
                            // find x+h
                            //
                            if (x_right + gradientOffset > 1.0){
                                xRight = (x_right + gradientOffset)-1.0;
                                gradInd[0] = ind[0]+1;
                            }else{
                                xRight = x_right + gradientOffset;
                                gradInd[0] = ind[0];
                            }

                            GetIndexandDistFromCenter(xRight, gradInd[0],  indexLeft, indexRight,  distFromLeft, distFromRight);
                            gradIndices[0] = indexLeft;    gradIndices[1] = indexRight;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                            vals[1] = TrilinearInterpolate(gradVals, x_left+gradientOffset, y_bottom, z_front);



                            //
                            // Y
                            //
                            for (int i=0; i<6; i++)
                                gradIndices[i] = indices[i];

                            //
                            // find y-h
                            //
                            if (y_top - gradientOffset < 0.0){
                                yTop = (y_top - gradientOffset)+1.0;
                                gradInd[1] = ind[1]-1;
                            }
                            else{
                                yTop = y_top - gradientOffset;
                                gradInd[1] = ind[1];
                            }

                            GetIndexandDistFromCenter(yTop, gradInd[1],  indexBottom, indexTop,  distFromBottom, distFromTop);
                            gradIndices[2] = indexBottom ;    gradIndices[3] = indexTop;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                            vals[2] = TrilinearInterpolate(gradVals, x_left, y_bottom-gradientOffset, z_front);

                            //
                            // find y+h
                            //
                            yTop = y_top;
                            if (y_top + gradientOffset > 1.0){
                                yTop = (y_top + gradientOffset)-1.0;
                                gradInd[1] = ind[1]+1;
                            }else{
                                yTop = y_top + gradientOffset;
                                gradInd[1] = ind[1];
                            }

                            GetIndexandDistFromCenter(yTop, gradInd[1],  indexBottom, indexTop,  distFromBottom, distFromTop);
                            gradIndices[2] = indexBottom;    gradIndices[3] = indexTop;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                            vals[3] = TrilinearInterpolate(gradVals, x_left, y_bottom+gradientOffset, z_front);


                            //
                            // Z
                            //
                            for (int i=0; i<6; i++)
                                gradIndices[i] = indices[i];

                            //
                            // z-h

                            if (z_back - gradientOffset < 0.0){
                                zBack = (z_back - gradientOffset)+1.0;
                                gradInd[2] = ind[2]-1;
                            }
                            else
                            {
                                zBack = z_back - gradientOffset;
                                gradInd[2] = ind[2];
                            }

                            GetIndexandDistFromCenter(zBack, gradInd[2],  indexFront, indexBack,  distFromFront, distFromBack);
                            gradIndices[4] = indexFront;    gradIndices[5] = indexBack;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                            vals[4] = TrilinearInterpolate(gradVals, x_left, y_bottom, z_front-gradientOffset);

                            //
                            // z+h
                            //
                            if (z_back + gradientOffset > 1.0){
                                zBack = (z_back + gradientOffset)-1.0;
                                gradInd[2] = ind[2]+1;
                            }else{
                                zBack = z_back + gradientOffset;
                                gradInd[2] = ind[2];
                            }

                            GetIndexandDistFromCenter(zBack, gradInd[2],  indexFront, indexBack,  distFromFront, distFromBack);
                            gradIndices[4] = indexFront;    gradIndices[5] = indexBack;
                            ComputeIndices(dims, gradIndices, indexGrad);
                            AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                            vals[5] = TrilinearInterpolate(gradVals, x_left, y_bottom, z_front+gradientOffset);


                            gradient[0] = (1.0/(2.0*gradientOffset)) * (vals[1] - vals[0]);
                            gradient[1] = (1.0/(2.0*gradientOffset)) * (vals[3] - vals[2]);
                            gradient[2] = (1.0/(2.0*gradientOffset)) * (vals[5] - vals[4]);

                            Normalize(gradient);
                        }

                        //
                        // Compute the color
                        //
                        ComputePixelColor(source_rgb, dest_rgb, gradient);
                    }
                }
            }
        }
        count++;
    }

    //
    // Set color
    imgArray[(y-yMin)*(imgWidth*4) + (x-xMin)*4 + 0] = std::min(std::max(dest_rgb[0],0.0),1.0);
    imgArray[(y-yMin)*(imgWidth*4) + (x-xMin)*4 + 1] = std::min(std::max(dest_rgb[1],0.0),1.0);
    imgArray[(y-yMin)*(imgWidth*4) + (x-xMin)*4 + 2] = std::min(std::max(dest_rgb[2],0.0),1.0);
    imgArray[(y-yMin)*(imgWidth*4) + (x-xMin)*4 + 3] = std::min(std::max(dest_rgb[3],0.0),1.0);
}


// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::getImageDimensions
//
//  Purpose:
//      Transfers the metadata of the patch
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//  Modifications:
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::GetImageDimensions(int &inUse, int dims[2], int screen_ll[2], int screen_ur[2], float &eyeDepth, float &clipDepth)
{
    inUse = patchDrawn;

    dims[0] = imgDims[0];    dims[1] = imgDims[1];

    screen_ll[0] = imgLowerLeft[0];     screen_ll[1] = imgLowerLeft[1];
    screen_ur[0] = imgUpperRight[0];    screen_ur[1] = imgUpperRight[1];

    eyeDepth = eyeSpaceDepth;
    clipDepth = clipSpaceDepth;
}


// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::getComputedImage
//
//  Purpose:
//      Allocates space to the pointer address and copy the image generated to it
//
//  Programmer: Pascal Grosset
//  Creation:   August 14, 2016
//
//  Modifications:
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::GetComputedImage(float *image)
{
    memcpy(image, imgArray, imgDims[0]*4*imgDims[1]*sizeof(float));

    if (imgArray != NULL)
        delete []imgArray;
    imgArray = NULL;
}


// ****************************************************************************
//  Method: avtSLIVRVoxelExtractor::ComputePixelColor
//
//  Purpose:
//      Computes color
//      By replicating avtPhong::AddLightingHeadlight
//
//  Programmer: Pascal Grosset
//  Creation:   June 10, 2013
//
//  Modifications:
//      Need to take into account lighting
//      Need to take into accoujnt multiple light sources
//
// ****************************************************************************

void
avtSLIVRVoxelExtractor::ComputePixelColor(double source_rgb[4], 
					  double dest_rgb[4],
					  float gradient[3])
{
    if (dest_rgb[3] >= 0.99)
    {
        patchDrawn = 1;
        return;
    }

    // Phong Shading
    if (lighting == true)
    {
        float dir[3];           
        dir[0] = -viewDirection[0];
        dir[1] = -viewDirection[1];
        dir[2] = -viewDirection[2];

        Normalize(dir);

        double temp_rgb[4];
        temp_rgb[0] = source_rgb[0];
        temp_rgb[1] = source_rgb[1];
        temp_rgb[2] = source_rgb[2];
        temp_rgb[3] = source_rgb[3];

        // cos(angle) = a.b;  angle between normal and light
        float normal_dot_light = Dot(gradient,dir);   // angle between light and normal;
        if (normal_dot_light < 0)
            normal_dot_light = -normal_dot_light;

        debug5 << "normal_dot_light: " << normal_dot_light << "   gradient: " << gradient[0] << ", " << gradient[1] << ", " << gradient[2] << std::endl;
        // Calculate color using phong shading
        // I = (I  * ka) + [ (I_i  * kd * (L.N)) + (Ia_i * ks * (R.V)^ns) ]_for each light source i
        for (int i=0; i<3; i++)
        {
            source_rgb[i] =  ( (materialProperties[0] + materialProperties[1] * normal_dot_light)           * source_rgb[i] ) +     // I  * ( ka + kd*abs(cos(angle)) )
                               (materialProperties[2] * pow((double)normal_dot_light,materialProperties[3]) * source_rgb[3] )  ;    // I  * kd*abs(cos(angle))

        }
    }
    for (int i=0; i<4; i++)
    {
        // front to back compositing
        dest_rgb[i] = source_rgb[i] * (1.0 - dest_rgb[3]) + dest_rgb[i];

        // // back to front
        //  dest_rgb[i] = std::min( dest_rgb[i] * (1.0 - intermediate_rgb[3]) + intermediate_rgb[i], 1.0);
    }


    patchDrawn = 1;
}


