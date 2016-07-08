/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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
//                         avtImgCommunicator.C                              //
// ************************************************************************* //
#include <cmath>
#include <avtParallel.h>
#include <ImproperUseException.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include <avtImgCommunicator.h>
#include <fstream>
#include <DebugStream.h>
#include <limits>


#ifdef PARALLEL

MPI_Datatype createImgDataType(){ 
  MPI_Datatype _img_mpi;
  const int numItems = 8;
  int blockLengths[numItems] = {1, 1,   1, 1,   2, 2, 2, 1};
  MPI_Datatype type[numItems] = { MPI_INT, MPI_INT,    MPI_INT, MPI_INT,   MPI_INT, MPI_INT, MPI_INT,   MPI_FLOAT};
  MPI_Aint offsets[numItems] = {0, sizeof(int), sizeof(int)*2, sizeof(int)*3, sizeof(int)*4, sizeof(int)*6, sizeof(int)*8, sizeof(int)*10};
  MPI_Type_struct(numItems, blockLengths,  offsets, type, &_img_mpi);
  
  return _img_mpi;
}

#endif


bool value_comparer(const std::pair<int,int> &before, const std::pair<int,int> &after){ return before.second < after.second; }
bool sortByVecSize(const std::vector<iotaMeta> &before, const std::vector<iotaMeta> &after){return before.size() > after.size();}
bool sortImgByCoordinatesIota(iotaMeta const& before, iotaMeta const& after){
  if(before.screen_ll[0] != after.screen_ll[0]) 
    return before.screen_ll[0] < after.screen_ll[0];
  else 
    return before.screen_ll[1] < after.screen_ll[1];
}
bool sortImgByDepthIota(iotaMeta const& before, iotaMeta const& after){ 
  if(before.avg_z != after.avg_z) 
    return (before.avg_z < after.avg_z);
  else 
    return (before.procId < after.procId); 
}

enum blendDirection {FRONT_TO_BACK = 0, BACK_TO_FRONT = 1};

// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
avtImgCommunicator::avtImgCommunicator()
{
  #ifdef PARALLEL
    MPI_Comm_size(VISIT_MPI_COMM, &num_procs);
    MPI_Comm_rank(VISIT_MPI_COMM, &my_id);

    MPI_Type_commit(&_img_mpi);
  #else
    num_procs = 1;
    my_id = 0;
  #endif
    
    totalPatches = 0;
    numPatchesToCompose = 0;

    processorPatchesCount = NULL;
    imgBuffer = NULL;

    allRecvIotaMeta = NULL;
    patchesToSendArray = NULL;
    patchesToRecvArray = NULL;
    numPatchesToSendArray = NULL;
    numPatchesToRecvArray = NULL;
    recvDisplacementForProcs = NULL;
    sendDisplacementForProcs = NULL;
    numPatchesToSendRecvArray = NULL;
    boundsPerBlockArray = NULL;
    blockDisplacementForProcs = NULL;
    numBlocksPerProc = NULL;
    patchesToCompositeLocallyArray = NULL;
    numPatchesToCompositeLocally = NULL;
    compositeDisplacementForProcs = NULL;

    compressedSizePerDiv = NULL;

    all_avgZ_proc0.clear();
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
avtImgCommunicator::~avtImgCommunicator(){
  if (my_id == 0){
    if (processorPatchesCount != NULL)
      delete []processorPatchesCount;

    if (allRecvIotaMeta != NULL)
      delete []allRecvIotaMeta;

    if (imgBuffer != NULL)
      delete []imgBuffer;

    if (compressedSizePerDiv != NULL)
      delete []compressedSizePerDiv;
    compressedSizePerDiv = NULL;


    if (patchesToSendArray!=NULL) delete[] patchesToSendArray;
        if (patchesToRecvArray!=NULL) delete[] patchesToRecvArray;
        if (numPatchesToSendArray!=NULL) delete[] numPatchesToSendArray;
        if (numPatchesToRecvArray!=NULL) delete[] numPatchesToRecvArray;
        if (recvDisplacementForProcs!=NULL) delete[] recvDisplacementForProcs;
        if (sendDisplacementForProcs!=NULL) delete[] sendDisplacementForProcs;
        if (blockDisplacementForProcs!=NULL) delete[] blockDisplacementForProcs;
        if (numPatchesToSendRecvArray!=NULL) delete[] numPatchesToSendRecvArray;
        if (boundsPerBlockArray!=NULL) delete[] boundsPerBlockArray;
        if (numBlocksPerProc!=NULL) delete[] numBlocksPerProc;
        if (patchesToCompositeLocallyArray!=NULL) delete[] patchesToCompositeLocallyArray;
        if (numPatchesToCompositeLocally!=NULL) delete[] numPatchesToCompositeLocally;
        if (compositeDisplacementForProcs!=NULL) delete[] compositeDisplacementForProcs;

        patchesToSendArray = NULL;
    patchesToRecvArray = NULL;
    numPatchesToSendArray = NULL;
    numPatchesToRecvArray = NULL;
    recvDisplacementForProcs = NULL;
    sendDisplacementForProcs = NULL;
    numPatchesToSendRecvArray = NULL;
    boundsPerBlockArray = NULL;
    blockDisplacementForProcs = NULL;
    numBlocksPerProc = NULL;
    patchesToCompositeLocallyArray = NULL;
    numPatchesToCompositeLocally = NULL;
    compositeDisplacementForProcs = NULL;
  }
  
  all_avgZ_proc0.clear();
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
int avtImgCommunicator::getDataPatchID(int procID, int patchID){
    int sumPatches = 0;
    for (int i=0; i<procID; i++)
        sumPatches += processorPatchesCount[i];
  
    return (sumPatches+patchID);
}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//    initialize 
//
//  Programmer: Pascal Grosset
//  Creation: July 2013  
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::init(){
}


// ****************************************************************************
//  Method: avtImgCommunicator::waitToSync
//
//  Purpose:
//      Wait for all processors to hearch here before continuing
//
//  Programmer: Pascal Grosset
//  Creation: July 2013  
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::syncAllProcs()
{
  #ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
  #endif
}



// ****************************************************************************
//  Method: avtImgCommunicator::gatherNumPatches
//
//  Purpose:
//    Get the number of patches each processor has
//
//  Arguments:
//    numPatches: number of patches a processor has
//
//  Programmer: Pascal Grosset
//  Creation: July 2013
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::gatherNumPatches(int numPatches)
{
  #ifdef PARALLEL
    int patchesProc[2];
    patchesProc[0] = my_id; patchesProc[1] = numPatches;
    int *tempRecvBuffer = NULL;

    if (my_id == 0){  
        processorPatchesCount = new int[num_procs];   // creates a buffer to hold the number of patches per processor
        tempRecvBuffer = new int[num_procs*2];
    }

    MPI_Gather(patchesProc, 2, MPI_INT,   tempRecvBuffer, 2,MPI_INT,         0, MPI_COMM_WORLD);    // all send to proc 0

    if (my_id == 0){    
        for (int i=0; i<num_procs; i++){
            processorPatchesCount[tempRecvBuffer[i*2]] = tempRecvBuffer[i*2 + 1]; // enter the number of patches for each processor
            totalPatches += processorPatchesCount[tempRecvBuffer[i*2]];       // count the number of patches
        }

        if (tempRecvBuffer != NULL)
            delete []tempRecvBuffer;
        tempRecvBuffer = NULL;
    }

  #endif
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose: 
//    Send the metadata needed by the root node to make decisions
//
//  Arguments:
//    arraySize   : the number of elements being sent
//    allIotaMetadata : the metadata bieng sent
//
//  Programmer: Pascal Grosset
//  Creation: July 2013
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::gatherIotaMetaData(int arraySize, float *allIotaMetadata)
{
  #ifdef PARALLEL
    int *recvSizePerProc = NULL;
    float *tempRecvBuffer = NULL;
    int *offsetBuffer = NULL;

    if (my_id == 0)
    {
        tempRecvBuffer = new float[totalPatches*8]; // x8: procId, patchNumber, dims[0], dims[1], screen_ll[0], screen_ll[1], eye_z, clip_z
        recvSizePerProc = new int[num_procs]; 
        offsetBuffer = new int[num_procs];

        for (int i=0; i<num_procs; i++)
        {
            if (i == 0)
                offsetBuffer[i] = 0;
            else
                offsetBuffer[i] = offsetBuffer[i-1] + recvSizePerProc[i-1];

            recvSizePerProc[i] = processorPatchesCount[i]*8;
        }
    }


    MPI_Gatherv(allIotaMetadata, arraySize, MPI_FLOAT,   tempRecvBuffer, recvSizePerProc, offsetBuffer,MPI_FLOAT,    0, MPI_COMM_WORLD);// all send to proc 0


    if (my_id == 0)
    {
        allRecvIotaMeta = new iotaMeta[totalPatches]; // allocate space to receive the many patches

        iotaMeta tempPatch;
        for (int i=0; i<totalPatches; i++)
        {
            tempPatch.procId        = (int) tempRecvBuffer[i*8 + 0];
            tempPatch.patchNumber   = (int) tempRecvBuffer[i*8 + 1];
            tempPatch.dims[0]       = (int) tempRecvBuffer[i*8 + 2];
            tempPatch.dims[1]       = (int) tempRecvBuffer[i*8 + 3];
            tempPatch.screen_ll[0]  = (int) tempRecvBuffer[i*8 + 4];
            tempPatch.screen_ll[1]  = (int) tempRecvBuffer[i*8 + 5];
            tempPatch.avg_z = tempPatch.eye_z = tempRecvBuffer[i*8 + 6];
            tempPatch.clip_z        =       tempRecvBuffer[i*8 + 7];

            int patchIndex = getDataPatchID(tempPatch.procId, tempPatch.patchNumber);
            allRecvIotaMeta[patchIndex] = setIota(tempPatch.procId, tempPatch.patchNumber, tempPatch.dims[0], tempPatch.dims[1], tempPatch.screen_ll[0], tempPatch.screen_ll[1], tempPatch.avg_z, tempPatch.clip_z);
            all_avgZ_proc0.insert(tempPatch.avg_z); //insert avg_zs into the set to keep a count of the total number of avg_zs
        }


        if (recvSizePerProc != NULL)
            delete []recvSizePerProc;
        recvSizePerProc = NULL;

        if (offsetBuffer != NULL)
            delete []offsetBuffer;
        offsetBuffer = NULL;

        if (tempRecvBuffer != NULL)
            delete []tempRecvBuffer;
        tempRecvBuffer = NULL;
    }
  #endif  
}





float avtImgCommunicator::clamp(float x){
    if (x > 1.0)
        x = 1.0;

    if (x < 0.0)
        x = 0.0;

    return x;
}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Pascal Grosset
//  Creation: July 2013  
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::getcompositedImage(int imgBufferWidth, int imgBufferHeight, unsigned char *wholeImage)
{
    for (int i=0; i< imgBufferHeight; i++)
        for (int j=0; j<imgBufferWidth; j++){
            int bufferIndex = (imgBufferWidth*4*i) + (j*4);
            int wholeImgIndex = (imgBufferWidth*3*i) + (j*3);

            wholeImage[wholeImgIndex+0] = (imgBuffer[bufferIndex+0] ) * 255;
            wholeImage[wholeImgIndex+1] = (imgBuffer[bufferIndex+1] ) * 255;
            wholeImage[wholeImgIndex+2] = (imgBuffer[bufferIndex+2] ) * 255;
        }


    if (imgBuffer != NULL)
        delete []imgBuffer;
    imgBuffer = NULL;
}








// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Pascal Grosset
//  Creation: July 2013  
//
//  Modifications:
//
// ****************************************************************************
imgMetaData avtImgCommunicator::setImg(int _inUse, int _procId, int _patchNumber, float dim_x, float dim_y, float screen_ll_x, float screen_ll_y, float screen_ur_x, float screen_ur_y, float _avg_z, float _clip_z){
  imgMetaData temp;
  temp.inUse = _inUse;
  temp.procId = _procId;
  temp.destProcId = _procId;
  temp.patchNumber = _patchNumber;
  temp.dims[0] = dim_x;         temp.dims[1] = dim_y;
  temp.screen_ll[0] = screen_ll_x;    temp.screen_ll[1] = screen_ll_y;
  temp.screen_ur[0] = screen_ur_x;    temp.screen_ur[1] = screen_ur_y;
  temp.avg_z = _avg_z;
  temp.eye_z = _avg_z;
  temp.clip_z = _clip_z;
  
  return temp;
}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Pascal Grosset
//  Creation: July 2013  
//
//  Modifications:
//
// ****************************************************************************
iotaMeta avtImgCommunicator::setIota(int _procId, int _patchNumber, int dim_x, int dim_y, int screen_ll_x, int screen_ll_y, float _avg_z, float _clip_z){
  iotaMeta temp;
  temp.procId = _procId;
  temp.patchNumber = _patchNumber;
  temp.dims[0] = dim_x;               temp.dims[1] = dim_y;
  temp.screen_ll[0] = screen_ll_x;    temp.screen_ll[1] = screen_ll_y;
  temp.avg_z  = _avg_z;
  temp.eye_z  = _avg_z;
  temp.clip_z = _clip_z;
  
  return temp;
}





// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************


void 
avtImgCommunicator::colorImage(float *& srcImage, int widthSrc, int heightSrc, float _color[4])
{
    for (int _y=0; _y<heightSrc; _y++)
        for (int _x=0; _x<widthSrc; _x++)
        {
            int srcIndex = widthSrc*_y*4 + _x*4;

            srcImage[srcIndex+0] = _color[0];
            srcImage[srcIndex+1] = _color[1];
            srcImage[srcIndex+2] = _color[2];
            srcImage[srcIndex+3] = _color[3];
        }
}


void 
avtImgCommunicator::placeInImage(float * srcImage, int srcExtents[4], float *& dstImage, int dstExtents[4])
{
    int widthSrc, heightSrc, widthDst;
    widthSrc  = srcExtents[1] - srcExtents[0];
    heightSrc = srcExtents[3] - srcExtents[2];

    widthDst  = dstExtents[1] - dstExtents[0];

    for (int _y=0; _y<heightSrc; _y++)
        for (int _x=0; _x<widthSrc; _x++)
        {
            int startingX = srcExtents[0];
            int startingY = srcExtents[2]; 

            int srcIndex = widthSrc*_y*4 + _x*4;                                                                  // index in the subimage 
            int dstIndex = ( (startingY+_y - dstExtents[2])*widthDst*4  + (startingX+_x - dstExtents[0])*4 );     // index in the big buffer

            dstImage[dstIndex+0] = srcImage[srcIndex+0];
            dstImage[dstIndex+1] = srcImage[srcIndex+1];
            dstImage[dstIndex+2] = srcImage[srcIndex+2];
            dstImage[dstIndex+3] = srcImage[srcIndex+3];
        }
}



void 
avtImgCommunicator::blendWithBackground(float *_image, int extents[4], float backgroundColor[4])
{
    int numPixels = (extents[3]-extents[2]) * (extents[1]-extents[0]);

    for (int index=0; index<numPixels; index++)      // estimated potential speedup: 2.240
    {
        int indexSrc = index*4;
        float alpha = (1.0 - _image[indexSrc+3]);

        _image[indexSrc+0] = backgroundColor[0] * alpha +  _image[indexSrc+0];
        _image[indexSrc+1] = backgroundColor[1] * alpha +  _image[indexSrc+1];
        _image[indexSrc+2] = backgroundColor[2] * alpha +  _image[indexSrc+2];
        _image[indexSrc+3] = backgroundColor[3] * alpha +  _image[indexSrc+3];
    }
}


void 
avtImgCommunicator::blendFrontToBack(float * srcImage, int srcExtents[4], int blendExtents[4], float *& dstImage, int dstExtents[4])
{
    int widthSrc, heightSrc, widthDst;
    widthSrc  = srcExtents[1] - srcExtents[0];
    heightSrc = srcExtents[3] - srcExtents[2];

    widthDst  = dstExtents[1] - dstExtents[0];

    for (int _y=blendExtents[2]; _y<blendExtents[3]; _y++)
        for (int _x=blendExtents[0]; _x<blendExtents[1]; _x++)
        {

            int srcIndex = (_y-srcExtents[2]) * widthSrc * 4 + (_x-srcExtents[0]) * 4;    
            int dstIndex = (_y-dstExtents[2]) * widthDst * 4 + (_x-dstExtents[0]) * 4;

            // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
            float alpha = 1.0 - dstImage[dstIndex+3];
            dstImage[dstIndex+0] = clamp( (srcImage[srcIndex+0] * alpha) + dstImage[dstIndex+0] );
            dstImage[dstIndex+1] = clamp( (srcImage[srcIndex+1] * alpha) + dstImage[dstIndex+1] );
            dstImage[dstIndex+2] = clamp( (srcImage[srcIndex+2] * alpha) + dstImage[dstIndex+2] );
            dstImage[dstIndex+3] = clamp( (srcImage[srcIndex+3] * alpha) + dstImage[dstIndex+3] );
        }
}

void 
avtImgCommunicator::blendBackToFront(float * srcImage, int srcExtents[4], int blendExtents[4], float *& dstImage, int dstExtents[4])
{
    int widthSrc, heightSrc, widthDst;
    widthSrc  = srcExtents[1] - srcExtents[0];
    heightSrc = srcExtents[3] - srcExtents[2];

    widthDst  = dstExtents[1] - dstExtents[0];

    for (int _y=blendExtents[2]; _y<blendExtents[3]; _y++)
        for (int _x=blendExtents[0]; _x<blendExtents[1]; _x++)
        {

            int srcIndex = (_y-srcExtents[2]) * widthSrc * 4 + (_x-srcExtents[0]) * 4;    
            int dstIndex = (_y-dstExtents[2]) * widthDst * 4 + (_x-dstExtents[0]) * 4;

            // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
            float alpha = 1.0 - srcImage[srcIndex+3];
            dstImage[dstIndex+0] = clamp( (dstImage[dstIndex+0] * alpha) + srcImage[srcIndex+0] );
            dstImage[dstIndex+1] = clamp( (dstImage[dstIndex+1] * alpha) + srcImage[srcIndex+1] );
            dstImage[dstIndex+2] = clamp( (dstImage[dstIndex+2] * alpha) + srcImage[srcIndex+2] );
            dstImage[dstIndex+3] = clamp( (dstImage[dstIndex+3] * alpha) + srcImage[srcIndex+3] );
        }
}


void 
avtImgCommunicator::blendFrontToBack(float * srcImage, int srcExtents[4], float *& dstImage, int dstExtents[4])
{
    int widthSrc, heightSrc, widthDst;
    widthSrc  = srcExtents[1] - srcExtents[0];
    heightSrc = srcExtents[3] - srcExtents[2];

    widthDst  = dstExtents[1] - dstExtents[0];

    for (int _y=0; _y<heightSrc; _y++)
        for (int _x=0; _x<widthSrc; _x++)
        {
            int startingX = srcExtents[0];
            int startingY = srcExtents[2]; 

            if ((startingX + _x) > dstExtents[1])
                continue;

            if ((startingY + _y) > dstExtents[3])
                continue;
            
            int srcIndex = widthSrc*_y*4 + _x*4;                                                                  // index in the subimage 
            int dstIndex = ( (startingY+_y - dstExtents[2])*widthDst*4  + (startingX+_x - dstExtents[0])*4 );     // index in the big buffer

            // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
            float alpha = 1.0 - dstImage[dstIndex+3];
            dstImage[dstIndex+0] = clamp( (srcImage[srcIndex+0] * alpha) + dstImage[dstIndex+0] );
            dstImage[dstIndex+1] = clamp( (srcImage[srcIndex+1] * alpha) + dstImage[dstIndex+1] );
            dstImage[dstIndex+2] = clamp( (srcImage[srcIndex+2] * alpha) + dstImage[dstIndex+2] );
            dstImage[dstIndex+3] = clamp( (srcImage[srcIndex+3] * alpha) + dstImage[dstIndex+3] );
        }
}

void 
avtImgCommunicator::blendBackToFront(float * srcImage, int srcExtents[4], float *& dstImage, int dstExtents[4])
{
    int widthSrc, heightSrc, widthDst;
    widthSrc  = srcExtents[1] - srcExtents[0];
    heightSrc = srcExtents[3] - srcExtents[2];

    widthDst  = dstExtents[1] - dstExtents[0];

    for (int _y=0; _y<heightSrc; _y++)
        for (int _x=0; _x<widthSrc; _x++)
        {
            int startingX = srcExtents[0];
            int startingY = srcExtents[2]; 

            if ((startingX + _x) > dstExtents[1])
                continue;

            if ((startingY + _y) > dstExtents[3])
                continue;
            
            int srcIndex = widthSrc*_y*4 + _x*4;                                                                  // index in the subimage 
            int dstIndex = ( (startingY+_y - dstExtents[2])*widthDst*4  + (startingX+_x - dstExtents[0])*4 );     // index in the big buffer

            // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
            float alpha = 1.0 - srcImage[srcIndex+3];
            dstImage[dstIndex+0] = clamp( (dstImage[dstIndex+0] * alpha) + srcImage[srcIndex+0] );
            dstImage[dstIndex+1] = clamp( (dstImage[dstIndex+1] * alpha) + srcImage[srcIndex+1] );
            dstImage[dstIndex+2] = clamp( (dstImage[dstIndex+2] * alpha) + srcImage[srcIndex+2] );
            dstImage[dstIndex+3] = clamp( (dstImage[dstIndex+3] * alpha) + srcImage[srcIndex+3] );
        }
}


void 
avtImgCommunicator::regionAllocation(int numMPIRanks, int *& regions)
{
    regions = new int[numMPIRanks];

    // Initial allocation: partition for section rank
    for (int i=0; i<numMPIRanks; i++)
        regions[i] = i;
}


void 
avtImgCommunicator::updateBoundingBox(int currentBoundingBox[4], int imageExtents[4])
{
    if ( (currentBoundingBox[0] == 0 && currentBoundingBox[1] == 0) && (currentBoundingBox[2] == 0 && currentBoundingBox[3] == 0))
    {
        currentBoundingBox[0]=imageExtents[0];
        currentBoundingBox[1]=imageExtents[1];
        currentBoundingBox[2]=imageExtents[2];
        currentBoundingBox[3]=imageExtents[3];

        return;
    }

    if (imageExtents[0] < currentBoundingBox[0])
        currentBoundingBox[0] = imageExtents[0];

    if (imageExtents[2] < currentBoundingBox[2])
        currentBoundingBox[2] = imageExtents[2];


    if (imageExtents[1] > currentBoundingBox[1])
        currentBoundingBox[1] = imageExtents[1];

    if (imageExtents[3] > currentBoundingBox[3])
        currentBoundingBox[3] = imageExtents[3];
}



void 
avtImgCommunicator::gatherDepthAtRoot(int numlocalPatches, float *localPatchesDepth, int &totalPatches, int *& patchCountPerRank, float *& allPatchesDepth)
{
  #ifdef PARALLEL
    //
    // Get how many patches are coming from each MPI rank
    totalPatches = 0;
    int *patchesOffset = NULL;


    if (my_id == 0) // root!
        patchCountPerRank = new int[num_procs]();

    MPI_Gather(&numlocalPatches, 1, MPI_INT,   patchCountPerRank, 1, MPI_INT,    0, MPI_COMM_WORLD);


    //
    // Gather number of patch group
    if (my_id == 0)
    {
        patchesOffset = new int[num_procs]();
        patchesOffset[0] = 0;

        for (int i=0; i<num_procs; i++)
        {
            totalPatches += patchCountPerRank[i];

            if (i == 0)
                patchesOffset[i] = 0;
            else
                patchesOffset[i] = patchesOffset[i-1] + patchCountPerRank[i-1]; 
        }

        allPatchesDepth = new float[totalPatches];
    }
    
    MPI_Gatherv(localPatchesDepth, numlocalPatches, MPI_FLOAT,   allPatchesDepth, patchCountPerRank, patchesOffset,    MPI_FLOAT, 0, MPI_COMM_WORLD);

    //
    // Cleanup
    if (my_id == 0)
        if (patchesOffset != NULL)
            delete []patchesOffset;
        
    patchesOffset = NULL;
  #endif
} 



//
// Serial Direct Send
//
void 
avtImgCommunicator::serialDirectSend(int numPatches, float *localPatchesDepth, int *extents, float *imgData, float backgroundColor[4], int width, int height)
{
  #ifdef PARALLEL
    //debug5 << "serialDirectSend" << std::endl;

    float *recvImage = NULL; 

    int tags[2] = {5781, 5782};

    int totalPatches;
    int *patchCountPerRank = NULL;
    float *patchesDepth = NULL;
    gatherDepthAtRoot(numPatches, localPatchesDepth, totalPatches, patchCountPerRank, patchesDepth);

    // debug5 << "Gather done!  totalPatches: " << totalPatches << std::endl;
    // if (my_id == 0)
    //     for (int i=0; i<totalPatches; i++)
    //         debug5 << "patchesDepth[" << i << "]: " << patchesDepth[i] << std::endl; 


    if (my_id == 0)
    {
        //
        // Root
        int srcSize[2], srcPos[2], dstSize[2], dstPos[2];
        srcSize[0] = width;  srcSize[1] = height;  
        srcPos[0] = 0;       srcPos[1] = 0;  

        //
        // Sort patches we will receive
        std::multimap<float,int> depthRankPatches;

        int index = 0;
        for (int i=0; i<num_procs; i++)
            for (int j=0; j<patchCountPerRank[i]; j++)
            {
                depthRankPatches.insert( std::pair<float,int>(patchesDepth[index],i) );
                index++;
            }
        

        //
        // Create space for buffers
        int recvParams[4];                          // minX, maxX, minY, maxY
        int imgExtents[4];
        imgExtents[0] = 0;  imgExtents[1] = width;
        imgExtents[2] = 0;  imgExtents[3] = height;

        recvImage = new float[width*height*4]();
        imgBuffer = new float[width*height*4]();

        int localIndex = 0;

        //
        // Compositing
        int __index = 0;
        for (std::multimap<float,int>::iterator it=depthRankPatches.begin(); it!=depthRankPatches.end(); ++it)
        {
            int rank = (*it).second;

            //debug5 << "\nRecv and blend from " << rank << " depth: " << (*it).first << std::endl;

            if (rank != my_id)
            {
                MPI_Recv(recvParams,             4, MPI_INT,   rank, tags[0],  MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // recv image info
                MPI_Recv(recvImage, width*height*4, MPI_FLOAT, rank, tags[1],  MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // recv image

                dstPos[0]  = dstPos[0];                      dstPos[1]  = dstPos[1];
                dstSize[0] = recvParams[2]-recvParams[0];    dstSize[1] = recvParams[3]-recvParams[1];

                //writeArrayToPPM("/home/pascal/Desktop/debugImages/recv_from_" + toStr(rank), recvImage, recvParams[1]-recvParams[0], recvParams[3]-recvParams[2]);
            }
            else
            {
                // It's local
                recvParams[0] = extents[ localIndex*4 + 0];
                recvParams[1] = extents[ localIndex*4 + 1];
                recvParams[2] = extents[ localIndex*4 + 2];
                recvParams[3] = extents[ localIndex*4 + 3];

                recvImage = &imgData[ localIndex*(width*height*4) ];
                localIndex++;
            }

            //debug5 << "Original extents: " << imgExtents[0] << ", " << imgExtents[1] << ", " << imgExtents[2] << ", " << imgExtents[3] << std::endl;
            //debug5 << "recvParams extents: " << recvParams[0] << ", " << recvParams[1] << ", " << recvParams[2] << ", " << recvParams[3] << std::endl;

            blendFrontToBack(recvImage, recvParams, imgBuffer, imgExtents);

            //writeArrayToPPM("/home/pascal/Desktop/debugImages/blended_with_" + toStr(__index), imgBuffer, width, height);
            //__index++;
        }


        //writeArrayToPPM("/home/pascal/Desktop/debugImages/full_" + toStr(__index), imgBuffer, imgExtents[1]-imgExtents[0], imgExtents[3]-imgExtents[2]);

        blendWithBackground(imgBuffer, imgExtents, backgroundColor);

        writeArrayToPPM("/home/pascal/Desktop/debugImages/full_with back_" + toStr(__index), imgBuffer, imgExtents[1]-imgExtents[0], imgExtents[3]-imgExtents[2]);

    }
    else
    {
        //
        // Sender
        for (int i=0; i<numPatches; i++)
        {
            int imgSize = (extents[i*4 + 1] - extents[i*4 + 0]) * (extents[i*4 + 3] - extents[i*4 + 2]) * 4;

            if (imgSize > 0)
            {
                //debug5 << "Sending: Extents " <<  extents[i*4 + 0] << ", " << extents[i*4 + 1] << ", " << extents[i*4 + 2] << ", " << extents[i*4 + 3] << std::endl;
                //writeArrayToPPM("/home/pascal/Desktop/debugImages/sending_to_root_from_" + toStr(my_id), &imgData[i*(width*height*4)], (extents[i*4 + 1] - extents[i*4 + 0]), (extents[i*4 + 3] - extents[i*4 + 2]) );

                MPI_Send( &extents[i*4],                       4, MPI_INT,   0, tags[0], MPI_COMM_WORLD);
                MPI_Send( &imgData[i*(width*height*4)],  imgSize, MPI_FLOAT, 0, tags[1], MPI_COMM_WORLD);
            }
        }
    }

    //debug5 << "Free memory" << std::endl; 

    //
    // Cleanup
    if (patchesDepth != NULL)
      delete []patchesDepth;
        
    if (patchCountPerRank != NULL)
        delete []patchCountPerRank; 

    if (recvImage != NULL)
        delete []recvImage;

    recvImage = NULL;
    patchCountPerRank = NULL;
    patchesDepth = NULL;


  #endif
}



// numInRegion / numRegions
void 
avtImgCommunicator::parallelDirectSend(float *imgData, int imgExtents[4], int region[], int numRegions, int tags[3], float backgroundColor[4], int width, int height)
{
  #ifdef PARALLEL
    // 
    // Determine position in region (myPositionInRegion)
    compositingDone = false;
    int myPositionInRegion = -1;
    bool inRegion = true;
    std::vector<int> regionVector(region, region+numRegions);
    std::vector<int>::iterator it = std::find(regionVector.begin(), regionVector.end(), my_id);

    if (it == regionVector.end())
    {
        inRegion = false;
        debug5 << my_id << " ~ SHOULD NOT HAPPEN: Not found " << my_id <<  " !!!" << std::endl;
    }
    else{
        myPositionInRegion = it - regionVector.begin();
    }


    //
    // Region boundaries
    int regionHeight = height/numRegions;
    int lastRegionHeight = height - regionHeight*(numRegions-1);


    // Extents of my region
    int myStartingHeight = myPositionInRegion * regionHeight;   
    int myEndingHeight = myStartingHeight + regionHeight;       
    if (myPositionInRegion == numRegions-1) 
        myEndingHeight = height;

    int myRegionHeight = myEndingHeight-myStartingHeight;


    // Size of one buffer
    int sizeOneBuffer = std::max(regionHeight,lastRegionHeight) * width * 4;




    //
    // MPI Async

    // Recv
    MPI_Request *recvMetaRq = new MPI_Request[ numRegions-1 ];
    MPI_Request *recvImageRq = new MPI_Request[ numRegions-1 ];

    MPI_Status *recvMetaSt = new MPI_Status[ numRegions-1 ];
    MPI_Status *recvImageSt = new MPI_Status[ numRegions-1 ];

    // Send
    MPI_Request *sendMetaRq = new MPI_Request[ numRegions-1 ];
    MPI_Request *sendImageRq = new MPI_Request[ numRegions-1 ];

    MPI_Status *sendMetaSt = new MPI_Status[ numRegions-1 ];
    MPI_Status *sendImageSt = new MPI_Status[ numRegions-1 ];



    //
    // Create Buffers

    // Create buffer for receiving images
    float *recvDataBuffer;
    recvDataBuffer = new float[ sizeOneBuffer * numRegions];


    // Create buffer for receiving messages
    std::vector<int> msgBuffer;
    msgBuffer.clear();
    msgBuffer.resize(5 * numRegions); 

    // Create buffer for sending messages
    int *sendExtents = new int[numRegions*5];

    //
    // Async Recv
    if (inRegion)
    {
        int recvCount=0;
        for (int i=0; i<numRegions; i++)
        {
            if ( regionVector[i] == my_id )
                continue;

            int src = regionVector[i];
            MPI_Irecv(&msgBuffer[i*5],                              5, MPI_INT,   src, tags[0], MPI_COMM_WORLD,  &recvMetaRq[recvCount] );
            MPI_Irecv(&recvDataBuffer[i*sizeOneBuffer], sizeOneBuffer, MPI_FLOAT, src, tags[1], MPI_COMM_WORLD,  &recvImageRq[recvCount] );
            recvCount++;
        }
    }


    
    //
    // Async Send
    int sendCount = 0;
    int sendingOffset;
    for (int i=0; i<numRegions; i++)
    {
        int regionStart, regionEnd, imgSize, dest;
        dest = regionVector[i];

        if ( dest == my_id )
            continue;

        regionStart = i*regionHeight;
        regionEnd = regionStart + regionHeight;
        if (i == numRegions-1) // the last one in region
            regionEnd = height;


        if (regionStart < imgExtents[2])
            regionStart = imgExtents[2];

        if (regionEnd > imgExtents[3])
            regionEnd = imgExtents[3];


        bool hasData = true;
        if (regionEnd - regionStart <= 0 || imgExtents[1]-imgExtents[0] <= 0)
        {
            hasData = false;

            sendingOffset = 0;
            imgSize = sendExtents[i*5 + 0] = sendExtents[i*5 + 1] = sendExtents[i*5 + 2] = sendExtents[i*5 + 3] =  sendExtents[i*5 + 4] = 0;
        }
        else
        {
            imgSize = (regionEnd-regionStart) * (imgExtents[1]-imgExtents[0]) * 4;
            sendingOffset = (regionStart-imgExtents[2]) * (imgExtents[1]-imgExtents[0]) * 4;

            sendExtents[i*5 + 0] = imgExtents[0];
            sendExtents[i*5 + 1] = imgExtents[1];
            sendExtents[i*5 + 2] = regionStart;
            sendExtents[i*5 + 3] = regionEnd;
            sendExtents[i*5 + 4] = 0;
        }

        //std::cout << my_id << " ~ i: " << i << "   regionVector[index]: " << regionVector[index] << "  extents: " <<  sendExtents[index*5 + 0] << ", " << sendExtents[index*5 + 1]  << ", " << sendExtents[index*5 + 2] << ", " << sendExtents[index*5 + 3] << "  sending ... " << std::endl;
        MPI_Isend(&sendExtents[i*5],             5,   MPI_INT, dest, tags[0], MPI_COMM_WORLD, &sendMetaRq[sendCount]);
        MPI_Isend(&imgData[sendingOffset], imgSize, MPI_FLOAT, dest, tags[1], MPI_COMM_WORLD, &sendImageRq[sendCount]);
        
        sendCount++;
    }



    //
    // Create buffer for region
    int myCompositedRegionExtents[4];
    myCompositedRegionExtents[0] = 0;                 myCompositedRegionExtents[1] = width;
    myCompositedRegionExtents[2] = myStartingHeight;  myCompositedRegionExtents[3] = myEndingHeight;

    float *myCompositedRegionImg = new float[width * (myEndingHeight-myStartingHeight)]();


    int recvImageExtents[4];
    float *recvImageData;
    

    //
    // Blend
    int numBlends = 0;
    int countBlend = 0;
    int boundingBox[4];
    boundingBox[0] = boundingBox[2] = std::numeric_limits<int>::max();
    boundingBox[1] = boundingBox[3] = std::numeric_limits<int>::min();

    if (inRegion)
    {
        for (int i=0; i<numRegions; i++)
        {
            int index = i;

            if (regionVector[index] == my_id)
            {
                int regionStart = myStartingHeight;
                int regionEnd = myEndingHeight;

                if (regionStart < imgExtents[2])
                    regionStart = imgExtents[2];

                if (regionEnd > imgExtents[3])
                    regionEnd = imgExtents[3];


                bool hasData = true;
                if (regionEnd - regionStart <= 0)
                {
                    hasData = false;
                    regionEnd = regionStart = 0;
                }

                if (hasData == true)
                {
                    int extentsSectionRecv[4];
                    extentsSectionRecv[0] = imgExtents[0];
                    extentsSectionRecv[1] = imgExtents[1];
                    extentsSectionRecv[2] = regionStart;
                    extentsSectionRecv[3] = regionEnd;

                    blendFrontToBack(imgData, imgExtents, extentsSectionRecv, myCompositedRegionImg, myCompositedRegionExtents);
                   
                    updateBoundingBox(boundingBox, extentsSectionRecv);
                    numBlends++;
                }
            }
            else
            {
                MPI_Wait(&recvMetaRq[countBlend], &recvMetaSt[countBlend]);

                for (int j=0; j<4; j++)
                    recvImageExtents[j] = msgBuffer[index*5 + j];


                bool hasData =  false;
                if (recvImageExtents[1]-recvImageExtents[0] > 0 && recvImageExtents[3]-recvImageExtents[2] > 0)
                {
                    hasData = true;
                    MPI_Wait(&recvImageRq[countBlend], &recvImageSt[countBlend]);
                    recvImageData = &recvDataBuffer[index*sizeOneBuffer];
                }

                
                if (hasData)
                {
                    blendFrontToBack(recvImageData, recvImageExtents, myCompositedRegionImg, myCompositedRegionExtents);

                    updateBoundingBox(boundingBox, recvImageExtents);
                    numBlends++;
                }

                countBlend++;
            }
        }
    }
    else
        compositingDone = true;

    
    // Gather
    gatherImages(region, num_procs, myCompositedRegionImg, myCompositedRegionExtents, tags[2], width, height);


    msgBuffer.clear();

    delete []recvDataBuffer;
    recvDataBuffer = NULL;


    if (numBlends == 0)
        boundingBox[0]=boundingBox[1]=boundingBox[2]=boundingBox[3] = 0;

    delete []recvMetaRq;
    delete []recvImageRq;
    delete []recvMetaSt;
    delete []recvImageSt;

    delete []sendMetaRq;
    delete []sendImageRq;
    delete []sendMetaSt;
    delete []sendImageSt;

    recvMetaRq = NULL;
    recvImageRq = NULL;
    recvMetaSt = NULL;
    recvImageSt = NULL;

    sendMetaRq = NULL;
    sendImageRq = NULL;
    sendMetaSt = NULL;
    sendImageSt = NULL;
  #endif
}


//imgBuffer
void avtImgCommunicator::gatherImages(int regionGather[], int numToRecv, float * inputImg, int imgExtents[4], int tag, int width, int height)
{
  #ifdef PARALLEL
    if (my_id == 0)
    {
        //
        // Receive at root/display node!
        imgBuffer = new float[width*height*4];
        int imgBufferExtents[4] = {0,0,0,0};
        imgBufferExtents[1] = width;
        imgBufferExtents[3] = height;


        int regionHeight      = height / numToRecv;
        int lastRegionHeight  = height - regionHeight*(numToRecv-1);
        int lastBufferSize    = lastRegionHeight * width * 4; 
        int regularBufferSize = regionHeight * width * 4;  

        for (int i=0; i<numToRecv; i++)
            if (regionGather[i] == my_id){
                numToRecv--;
                break;
            }
        
        //
        // Create buffers for async reciving
        MPI_Request *recvImageRq = new MPI_Request[ numToRecv ];
        MPI_Status  *recvImageSt = new MPI_Status[ numToRecv ];

        // Async Recv
        int recvCount=0;
        for (int i=0; i<numToRecv; i++)
        {
            int src = regionGather[i];
            if (src == my_id)
                continue;

            if (i == numToRecv-1)
                MPI_Irecv(&imgBuffer[i*regularBufferSize], lastBufferSize,     MPI_FLOAT, src, tag, MPI_COMM_WORLD,  &recvImageRq[recvCount] ); 
            else
                MPI_Irecv(&imgBuffer[i*regularBufferSize], regularBufferSize,  MPI_FLOAT, src, tag, MPI_COMM_WORLD,  &recvImageRq[recvCount] );
            recvCount++;
        }


        if (compositingDone == false)   // If root has data for the final image
            placeInImage(inputImg, imgBufferExtents, imgBuffer, imgExtents);

        MPI_Waitall(numToRecv, recvImageRq, recvImageSt);

        delete []recvImageRq;
        recvImageRq = NULL;
        delete []recvImageSt;
        recvImageSt = NULL;
    }
    else
    {
        if (compositingDone == false)   // If root has data for the final image
        {
            int imgSize = (imgExtents[1]-imgExtents[0]) * (imgExtents[3]-imgExtents[2]) * 4;
            MPI_Send(inputImg, imgSize, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);    
          
            compositingDone = true;
        }
    }
  #endif
}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************


void avtImgCommunicator::barrier(){
  #ifdef PARALLEL
    MPI_Barrier( MPI_COMM_WORLD );
  #endif
}
