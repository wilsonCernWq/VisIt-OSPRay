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
avtImgCommunicator::avtImgCommunicator(){

#ifdef PARALLEL
  MPI_Comm_size(VISIT_MPI_COMM, &num_procs);
  MPI_Comm_rank(VISIT_MPI_COMM, &my_id);

  _img_mpi = createMetaDataType();
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

#ifdef PARALLEL
  MPI_Type_free(&_img_mpi);
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
void avtImgCommunicator::syncAllProcs(){
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
void avtImgCommunicator::gatherNumPatches(int numPatches){

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


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Manasa Prasad
//  Creation: August 2013
//
//  Modifications:
//
// ****************************************************************************
void determinePatchesToCompositeLocally(const std::vector<iotaMeta>& all_patches_sorted_avgZ_proc0, std::vector<std::vector<iotaMeta> >& patchesToCompositeLocallyVector, const int& procToSend, std::vector<std::vector<int> >& divisions){

  if(all_patches_sorted_avgZ_proc0.size() == 0) return;
  iotaMeta prev_data = *all_patches_sorted_avgZ_proc0.begin();
  bool already_in = false;

  iotaMeta delimiter;
  delimiter.patchNumber = -1;

  for (size_t i=1; i < all_patches_sorted_avgZ_proc0.size(); ++i){
    if(prev_data.procId == all_patches_sorted_avgZ_proc0[i].procId && all_patches_sorted_avgZ_proc0[i].procId != procToSend){
      if(!already_in) {
        patchesToCompositeLocallyVector[prev_data.procId].push_back(delimiter);
        patchesToCompositeLocallyVector[prev_data.procId].push_back(prev_data);
        divisions[prev_data.procId].push_back(patchesToCompositeLocallyVector[prev_data.procId].size() - 1);
        already_in = true;
      }
      patchesToCompositeLocallyVector[prev_data.procId].push_back(all_patches_sorted_avgZ_proc0[i]);
    }else{
      prev_data = all_patches_sorted_avgZ_proc0[i];
      already_in = false;
    }
  }

  for (size_t i=0; i < patchesToCompositeLocallyVector.size(); ++i){
    divisions[i].push_back(patchesToCompositeLocallyVector[i].size() + 1);
  }
}

// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Manasa Prasad
//  Creation: August 2013
//
//  Modifications:
//
// ****************************************************************************

bool adjacencyTest(const iotaMeta& patch_1, const iotaMeta& patch_2){
  if( (patch_2.screen_ll[0] == patch_1.screen_ll[0] /*+ patch_1.dims[0]*/) && 
    (patch_2.screen_ll[1] <= patch_1.screen_ll[1] + patch_1.dims[1])) 
    return true;
  return false;
}

void determinePatchAdjacency(std::vector<iotaMeta>& allPatchesSorted, std::vector<std::vector<iotaMeta> >& patchesToComposite, std::vector<std::vector<int> >& divisions ){

  if(allPatchesSorted.size() == 0) return;

  iotaMeta delimiter; 
  delimiter.patchNumber = -1;

  std::vector<iotaMeta>::iterator it;

  iotaMeta prevPatch, currentPatch, nextPatch;
  int lower, upper;

  for (size_t i=0; i < patchesToComposite.size(); ++i){
    for(size_t k=0; k<divisions[i].size()-1; k++){

      lower = divisions[i][k];
      upper = divisions[i][k+1]-1;

      std::sort(  patchesToComposite[i].begin()+lower, 
            patchesToComposite[i].begin()+upper,  &sortImgByCoordinatesIota);

      for(int j = lower; j < upper; j++){

        currentPatch = patchesToComposite[i][j];
        if(j < upper-1)   nextPatch = patchesToComposite[i][j+1];
        if(j > lower)   prevPatch = patchesToComposite[i][j-1];
      
        //
        // if patch is adjacent to the preceding patch, remove it from the allPatchesSorted list
        //
        if (j > lower && adjacencyTest(prevPatch, currentPatch)) {
          it = std::find(allPatchesSorted.begin(), allPatchesSorted.end(), currentPatch); 
          allPatchesSorted.erase(it);
        }

        //
        // if patch is adjacent to the succeeding patch, check if:
        //          patch lies somewhere in the middle => insert a delimiter before patch
        //
        else if (j < upper-1 && adjacencyTest(currentPatch, nextPatch)){
          if (j != lower){ //insert a -1      
            it = std::find(patchesToComposite[i].begin(), patchesToComposite[i].end(), currentPatch);
            patchesToComposite[i].insert(it, delimiter);

            for(size_t m = k; m < divisions[i].size() - 1; m++)
              divisions[i][m+1]++;
            j++; upper++;

            int position = it - patchesToComposite[i].begin();
            lower = position + 1;
          }
          // else do nothing. Let the first patch remain
        }

        //
        // if patch does not lie beside any adjoining patch, remove it from patchesToComposite
        //
        else{
          it = std::find(patchesToComposite[i].begin(), patchesToComposite[i].end(), currentPatch);
          patchesToComposite[i].erase(it);

          // if this is the only remaining element, also remove the -1
          if(upper == lower+1){ 
            patchesToComposite[i].erase(--it);
            for(size_t m = k; m < divisions[i].size() - 1; m++)
              divisions[i][m+1] -= 2;
          }else{
            for(size_t m = k; m < divisions[i].size() - 1; m++)
              divisions[i][m+1]--;
            j--; upper--;

          }
        }
        // end of j loop
      }
      // end of k loop
    }
    // end of i loop
  }
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Manasa Prasad
//  Creation: August 2013
//
//  Modifications:
//
// ****************************************************************************
int calculatePatchDivision (const std::vector<iotaMeta>& imgVector, std::vector<int>& procsAlreadyInList){
std::map<int,int> numPatchesPerProc;
  std::pair<std::map<int,int>::iterator, bool> isPresent;

  if (imgVector.size() == 0) return 0;

  for (size_t i = 0; i < imgVector.size(); ++i){
    const bool is_in = ((std::find(procsAlreadyInList.begin(), procsAlreadyInList.end(), imgVector[i].procId)) != (procsAlreadyInList.end()));
    if(!is_in){
      isPresent = numPatchesPerProc.insert (  std::pair<int,int>(imgVector[i].procId, imgVector[i].dims[0] * imgVector[i].dims[1]) );
      if(isPresent.second == false)     
        numPatchesPerProc[imgVector[i].procId] += imgVector[i].dims[0] * imgVector[i].dims[1];
    }
  }

  if(numPatchesPerProc.size() == 0){
    for (size_t i = 0; i < imgVector.size(); ++i){
      isPresent = numPatchesPerProc.insert (  std::pair<int,int>(imgVector[i].procId, imgVector[i].dims[0] * imgVector[i].dims[1]) );
      if(isPresent.second == false)     
        numPatchesPerProc[imgVector[i].procId] += imgVector[i].dims[0] * imgVector[i].dims[1];
    }
  }

  return std::max_element(numPatchesPerProc.begin(), numPatchesPerProc.end(), value_comparer)->first;
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Manasa Prasad
//  Creation: August 2013
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::patchAllocationLogic(){

  // 1. do not send the numpatchespreproc vector!
  // 2. do away with numcompositedPatches test and replace with totalRecvPatches

  std::vector<std::vector<iotaMeta> > all_patches_sorted_avgZ_proc0(num_procs); 
  std::vector<int> procToSend;
  std::vector<int> numAvgZEachBlock (num_procs);

  numBlocksPerProc = new int[num_procs]();
  boundsPerBlockVec.resize(num_procs);

  int num_avgZ_proc0 = all_avgZ_proc0.size();

  // calculate the range of avg_zs
  int num_avgZ_perBlock = (num_avgZ_proc0 / num_procs);
  int rem_avgZ = num_avgZ_proc0 % num_procs;

  //printf("num_avg_z: %d num_procs: %d num_avgZ_perBlock: %d rem_avgZ: %d \n\n", num_avgZ_proc0, num_procs, num_avgZ_perBlock, rem_avgZ);
  
  //procToSend.resize           (num_procs);
  //all_patches_sorted_avgZ_proc0.resize(num_procs);

  for(int i = 0; i < num_procs; i++){
    numAvgZEachBlock[i] = num_avgZ_perBlock + (rem_avgZ-- > 0 ? 1 : 0);
  }

  // Sorting the patches
  std::sort(allRecvIotaMeta, allRecvIotaMeta + totalPatches, &sortImgByDepthIota);

  // Populate the all_patches_sorted_avgZ_proc0 vector with data from allPatchesMeta
  int current_avgZ_Block = 0;
  float prevAvgZ = totalPatches > 0 ? allRecvIotaMeta[0].avg_z : 0;
  int num_avgZ_thisBlock = 1;

  for (int i = 0; i < totalPatches; ++i){
    float currentAvgZ = allRecvIotaMeta[i].avg_z;

    if (currentAvgZ != prevAvgZ && num_avgZ_thisBlock == numAvgZEachBlock[current_avgZ_Block]){
      num_avgZ_thisBlock = 1;
      current_avgZ_Block++;
      prevAvgZ = currentAvgZ;
    } else if (currentAvgZ != prevAvgZ){
      prevAvgZ = currentAvgZ;
      num_avgZ_thisBlock++;
    }

    all_patches_sorted_avgZ_proc0[current_avgZ_Block].push_back(allRecvIotaMeta[i]);
  }

  //sort the avg_z vector by size
  std::sort (all_patches_sorted_avgZ_proc0.begin(), all_patches_sorted_avgZ_proc0.end(), sortByVecSize);
  std::vector<std::vector<iotaMeta> > patchesToCompositeLocallyVector(num_procs);

  // Calculate which block of avg_z to send to which processor
  // ith block is sent to the processor in procToSend[i]

  //printf("num_avg_z: %d num_procs: %d\n\n", num_avgZ_proc0, num_procs);
  for (int i = 0; i < num_procs; ++i){
    
    procToSend.push_back(calculatePatchDivision(all_patches_sorted_avgZ_proc0[i], procToSend));
    std::vector<std::vector<int> > divisions(num_procs);
    determinePatchesToCompositeLocally(all_patches_sorted_avgZ_proc0[i], patchesToCompositeLocallyVector, procToSend[i], divisions); 
    determinePatchAdjacency(all_patches_sorted_avgZ_proc0[i], patchesToCompositeLocallyVector, divisions);

    int size = all_patches_sorted_avgZ_proc0[i].size();
    if(size > 0){ 
      boundsPerBlockVec[procToSend[i]].push_back(all_patches_sorted_avgZ_proc0[i][0].avg_z);
      boundsPerBlockVec[procToSend[i]].push_back(all_patches_sorted_avgZ_proc0[i][size-1].avg_z);
    }
    //printf("division: %d, procToSend: %d \n", i, procToSend[i]);
  }

  // Printing for error check
  for(int i = 0; i < num_procs; ++i){
    debug5 <<  "Block: " << i << " \t size:" << all_patches_sorted_avgZ_proc0[i].size() << endl;
  //  printf("Block %d\n", i );
    for(std::vector<iotaMeta>::iterator it = all_patches_sorted_avgZ_proc0[i].begin(); it != all_patches_sorted_avgZ_proc0[i].end(); ++it){
      debug5 <<  it->avg_z << "\t " << it->procId << endl;
  //    printf("\t %.5f \t %d\n", it->avg_z, it->procId);
    }
  }


  std::vector< std::vector<int> >  patchesToSendVec (num_procs);
  std::vector< std::map<int,int> > patchesToRecvMap (num_procs); 
  std::pair< std::map<int,int>::iterator, bool > isPresent;

  for (int procId = 0; procId < num_procs; procId++)  patchesToRecvMap[procId].clear();
  
  // Populate recvMap and sendVec
  for (size_t i = 0; i < all_patches_sorted_avgZ_proc0.size(); i++){

    int destProcId = procToSend[i];
    for (size_t j = 0; j < all_patches_sorted_avgZ_proc0[i].size(); j++){

      int originProcId = all_patches_sorted_avgZ_proc0[i][j].procId;
      if(originProcId != destProcId){

        // Array of the form [0 1 2 3 ...]
        //           even numbers(0,2..): patchNumber 
        //           odd numbers(1,3...): destProcId
          patchesToSendVec[originProcId].push_back( all_patches_sorted_avgZ_proc0[i][j].patchNumber );
        patchesToSendVec[originProcId].push_back( destProcId );

        // Array of the form [0 1 2 3 ...]
        //           even numbers(0,2..): procId 
        //           odd numbers(1,3...): numPatches
          isPresent = patchesToRecvMap[destProcId].insert( std::pair<int,int>(originProcId, 1) );
          if(isPresent.second == false)   ++patchesToRecvMap[destProcId][originProcId];
          
          //printf("Inserted to dest: %d the origin: %d avg_z: %.2f\n", destProcId, originProcId, allPatchesMeta[patchIndex].avg_z );
      }
    }

    debug5 <<  "------division: " << i << " procToSend:" << procToSend[i] << endl;
  }


  recvDisplacementForProcs  = new int[num_procs]();
  sendDisplacementForProcs  = new int[num_procs]();
  numPatchesToSendArray     = new int[num_procs]();
  numPatchesToRecvArray     = new int[num_procs]();
  blockDisplacementForProcs   = new int[num_procs]();
  numPatchesToSendRecvArray   = new int[4*num_procs]();

  compositeDisplacementForProcs = new int[num_procs]();
  numPatchesToCompositeLocally = new int[num_procs]();

  int totalRecvSize = 0;
  int totalSendSize = 0;
  int totalBlockSize = 0;
  int totalCompositeSize = 0;

  for (int currentProcId = 0; currentProcId < num_procs; currentProcId++){

    numPatchesToSendArray[currentProcId] = patchesToSendVec[currentProcId].size();
    numPatchesToRecvArray[currentProcId] = 2*patchesToRecvMap[currentProcId].size();
    numBlocksPerProc[currentProcId]    = boundsPerBlockVec[currentProcId].size();
    numPatchesToCompositeLocally[currentProcId] = patchesToCompositeLocallyVector[currentProcId].size();

    numPatchesToSendRecvArray[currentProcId*4]    = numPatchesToSendArray[currentProcId];
    numPatchesToSendRecvArray[currentProcId*4 + 1]  = numPatchesToRecvArray[currentProcId];
    numPatchesToSendRecvArray[currentProcId*4 + 2]  = numBlocksPerProc[currentProcId];
    numPatchesToSendRecvArray[currentProcId*4 + 3]  = numPatchesToCompositeLocally[currentProcId];

    totalRecvSize += numPatchesToRecvArray[currentProcId];
    totalSendSize += numPatchesToSendArray[currentProcId];
    totalBlockSize+= numBlocksPerProc[currentProcId];
    totalCompositeSize += numPatchesToCompositeLocally[currentProcId];

    if(currentProcId!=0) {
      recvDisplacementForProcs[currentProcId] = recvDisplacementForProcs[currentProcId-1] + numPatchesToRecvArray[currentProcId-1];
      sendDisplacementForProcs[currentProcId] = sendDisplacementForProcs[currentProcId-1] + numPatchesToSendArray[currentProcId-1];
      blockDisplacementForProcs[currentProcId] = blockDisplacementForProcs[currentProcId-1] + numBlocksPerProc[currentProcId-1];
      compositeDisplacementForProcs[currentProcId] = compositeDisplacementForProcs[currentProcId-1] + numPatchesToCompositeLocally[currentProcId-1];
    }

  }

  patchesToRecvArray  = new int[totalRecvSize];
  patchesToSendArray  = new int[totalSendSize];
  boundsPerBlockArray = new float[totalBlockSize];
  patchesToCompositeLocallyArray = new int[totalCompositeSize];

  for (int currentProcId = 0; currentProcId < num_procs; currentProcId++){

    //printf("\n\n### dest: %d size: %ld\n ", currentProcId, patchesToRecvMap[currentProcId].size());
    int count = 0;
    for (std::map<int,int>::iterator it = patchesToRecvMap[currentProcId].begin(); it!=patchesToRecvMap[currentProcId].end(); ++it){
        patchesToRecvArray[recvDisplacementForProcs[currentProcId] + (count++)] = it->first;
        patchesToRecvArray[recvDisplacementForProcs[currentProcId] + (count++)] = it->second;
        //printf("\t numPatches %d origin: %d\n", it->second, it->first);
    }

    count = 0;
    //printf("\n\n### origin: %d size: %ld\n ", currentProcId, patchesToSendVec[currentProcId].size());
    for (std::vector<int>::iterator it = patchesToSendVec[currentProcId].begin(); it!=patchesToSendVec[currentProcId].end(); ++it){
      patchesToSendArray[sendDisplacementForProcs[currentProcId] + (count++)] = *it;

      //if(count%2!=0) printf("\t patchNumber %d", *it);
      //else       printf(" dest: %d\n", *it);

    }

    count = 0;
    for (std::vector<float>::iterator it = boundsPerBlockVec[currentProcId].begin(); it!=boundsPerBlockVec[currentProcId].end(); ++it)
      boundsPerBlockArray[blockDisplacementForProcs[currentProcId] + (count++)] = *it;

    count = 0;
    for (std::vector<iotaMeta>::iterator it = patchesToCompositeLocallyVector[currentProcId].begin(); it!=patchesToCompositeLocallyVector[currentProcId].end(); ++it)
      patchesToCompositeLocallyArray[compositeDisplacementForProcs[currentProcId] + (count++)] = it->patchNumber;

  }

  //printf("\n ..................................................\n");

}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Manasa Prasad
//  Creation: July 2013
//
//  Modifications:
//
// ****************************************************************************

void
avtImgCommunicator::scatterNumDataToCompose(int &totalSendData, int &totalRecvData,
    int &numDivisions, int &totalPatchesToCompositeLocally)
{
  int totalSendRecvData[4] = {0,0,0,0};
  #ifdef PARALLEL
    MPI_Scatter(numPatchesToSendRecvArray, 4, MPI_INT, totalSendRecvData, 4, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

    totalSendData = totalSendRecvData[0];
    totalRecvData = totalSendRecvData[1];
    numDivisions = totalSendRecvData[2]; // twice the  number of divisions
    totalPatchesToCompositeLocally = totalSendRecvData[3];
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: Manasa Prasad
//  Creation: July 2013
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::scatterDataToCompose(  int& totalSendData, int* informationToSendArray, 
                        int& totalRecvData, int* informationToRecvArray, 
                        int& numDivisions, float* blocksPerProc,
                        int& totalPatchesToCompositeLocally, int* patchesToCompositeLocally){
  #ifdef PARALLEL
        MPI_Scatterv(patchesToSendArray, numPatchesToSendArray, sendDisplacementForProcs, MPI_INT, informationToSendArray, totalSendData, MPI_INT, 0,MPI_COMM_WORLD);
        MPI_Scatterv(patchesToRecvArray, numPatchesToRecvArray, recvDisplacementForProcs, MPI_INT, informationToRecvArray, totalRecvData, MPI_INT, 0,MPI_COMM_WORLD);
        MPI_Scatterv(boundsPerBlockArray, numBlocksPerProc, blockDisplacementForProcs, MPI_FLOAT, blocksPerProc, numDivisions, MPI_FLOAT, 0,MPI_COMM_WORLD);
        MPI_Scatterv(patchesToCompositeLocallyArray, numPatchesToCompositeLocally, compositeDisplacementForProcs, MPI_INT, patchesToCompositeLocally, totalPatchesToCompositeLocally, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

    if(my_id == 0)
    {
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
void avtImgCommunicator::sendPointToPoint(imgMetaData toSendMetaData, imgData toSendImgData, int tag){
  #ifdef PARALLEL
    // Commit the datatype
    MPI_Datatype _TestImg_mpi;
    _TestImg_mpi = createImgDataType();
    MPI_Type_commit(&_TestImg_mpi);

      MPI_Send(&toSendMetaData, 1, _TestImg_mpi, toSendMetaData.destProcId, tag, MPI_COMM_WORLD); // 2
    MPI_Send(toSendImgData.imagePatch, toSendMetaData.dims[0]*toSendMetaData.dims[1]*4, MPI_FLOAT, toSendMetaData.destProcId, tag-1, MPI_COMM_WORLD); //send the image data // 1
   
      MPI_Type_free(&_TestImg_mpi);
    #endif
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
void avtImgCommunicator::recvPointToPoint(imgMetaData &recvMetaData, imgData &recvImgData){
  #ifdef PARALLEL
    // Commit the datatype
    MPI_Datatype _TestImg_mpi;
    _TestImg_mpi = createImgDataType();
    MPI_Type_commit(&_TestImg_mpi);

        MPI_Recv (&recvMetaData, 1, _TestImg_mpi, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);

        recvImgData.procId = recvMetaData.procId;
        recvImgData.patchNumber = recvMetaData.patchNumber;
        recvImgData.imagePatch = new float[recvMetaData.dims[0]*recvMetaData.dims[1] * 4];

        MPI_Recv(recvImgData.imagePatch, recvMetaData.dims[0]*recvMetaData.dims[1]*4, MPI_FLOAT, status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);

        MPI_Type_free(&_TestImg_mpi);
    #endif
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
void avtImgCommunicator::recvPointToPointMetaData(imgMetaData &recvMetaData, int tag){
  #ifdef PARALLEL
    // Commit the datatype
    MPI_Datatype _TestImg_mpi;
    _TestImg_mpi = createImgDataType();
    MPI_Type_commit(&_TestImg_mpi);

        MPI_Recv (&recvMetaData, 1, _TestImg_mpi, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status); // 2

        MPI_Type_free(&_TestImg_mpi);
    #endif
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
void avtImgCommunicator::recvPointToPointImgData(imgMetaData recvMetaData, imgData &recvImgData, int tag){
  #ifdef PARALLEL
        MPI_Recv(recvImgData.imagePatch, recvMetaData.dims[0]*recvMetaData.dims[1]*4, MPI_FLOAT, status.MPI_SOURCE, tag-1, MPI_COMM_WORLD, &status); // 1
    #endif
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
void avtImgCommunicator::gatherEncodingSizes(int *sizeEncoding, int numDivisions){

  #ifdef PARALLEL
  int *offsetBuffer = NULL;
  int *recvSizePerProc = NULL;
  int totalDivisions = 0;

    // Only proc 0 receives data
    if (my_id == 0){

      recvSizePerProc = new int[num_procs];
      offsetBuffer = new int[num_procs];
      
      for (int i=0; i<num_procs; i++){
        int numBoundsPerBlock = boundsPerBlockVec[i].size()/2;
        totalDivisions += numBoundsPerBlock;
        recvSizePerProc[i] = numBoundsPerBlock;

        if (i == 0)
          offsetBuffer[i] = 0;
        else
          offsetBuffer[i] = offsetBuffer[i-1] + recvSizePerProc[i-1];

      }
      compressedSizePerDiv = new int[totalDivisions];
    }

        //  send   recv  others
    MPI_Gatherv(sizeEncoding, numDivisions, MPI_INT,    compressedSizePerDiv, recvSizePerProc, offsetBuffer,MPI_INT,      0, MPI_COMM_WORLD); // all send to proc 0



    if (my_id == 0){
      for (int i=0; i<totalDivisions; i++)
        debug5 <<  "  0 div: " << i << " : " << compressedSizePerDiv[i] << endl;

      if (offsetBuffer != NULL) 
        delete []offsetBuffer;
      offsetBuffer = NULL;

      if (recvSizePerProc != NULL)  
        delete []recvSizePerProc;
      recvSizePerProc = NULL;
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
void avtImgCommunicator::gatherAndAssembleEncodedImages(int sizex, int sizey, int sizeSending, float *images, int numDivisions){
  #ifdef PARALLEL
    float *tempRecvBuffer = NULL;
    int *recvSizePerProc = NULL;
    int *offsetBuffer = NULL;
    int totalDivisions = 0;
    std::map<float,int> depthPartitions;  //float: z-value, int: division index

    // Only proc 0 receives data
    if (my_id == 0){
      int divIndex = 0;
      int totalSize = 0;
      recvSizePerProc = new int[num_procs]; 
      offsetBuffer = new int[num_procs];  

      for (int i=0; i<num_procs; i++){
        int numBoundsPerBlock = boundsPerBlockVec[i].size()/2;
        totalDivisions += numBoundsPerBlock;
        
        int sizeEncoded = 0;
        for (size_t j=0; j<boundsPerBlockVec[i].size(); j+=2){
          depthPartitions.insert( std::pair< float,int > ( (boundsPerBlockVec[i][j] + boundsPerBlockVec[i][j+1]) * 0.5, divIndex ) );   
          sizeEncoded += compressedSizePerDiv[divIndex]*5;
          divIndex++;
        }
        totalSize += sizeEncoded;
        recvSizePerProc[i] = sizeEncoded;

        if (i == 0)
          offsetBuffer[i] = 0;
        else
          offsetBuffer[i] = offsetBuffer[i-1] + recvSizePerProc[i-1];
      }

      tempRecvBuffer = new float[ totalSize ];
    }

        //  send   recv  others
    MPI_Gatherv(images, sizeSending, MPI_FLOAT,    tempRecvBuffer, recvSizePerProc, offsetBuffer,MPI_FLOAT,        0, MPI_COMM_WORLD);    // all send to proc 0

    if (my_id == 0){
      // Create a buffer to store the composited image
      imgBuffer = new float[sizex * sizey * 4];
      
      // Front to back
      for (int i=0; i<(sizex * sizey * 4); i+=4){
        imgBuffer[i+0] = background[0]/255.0; 
        imgBuffer[i+1] = background[1]/255.0; 
        imgBuffer[i+2] = background[2]/255.0; 
        imgBuffer[i+3] = 1.0;
      }

      int count = 0;
      int offset = 0;
      int index = 0;

      if (depthPartitions.size() > 0){
        std::map< float,int >::iterator it;
        it=depthPartitions.end();
        it=depthPartitions.end();
        do{
          --it;
          debug5 << it->first << " => " << it->second << "    compressedSizePerDiv[count]: " << compressedSizePerDiv[count] << std::endl;

          imageBuffer temp;
          temp.image = new float[sizex*sizey*4];
          index = it->second;

          if (index == 0)
            offset = 0;
          else{
            offset = 0;
            for (int k=0; k<index; k++)
              offset += compressedSizePerDiv[k];
          }

          rleDecode(compressedSizePerDiv[index], tempRecvBuffer, offset*5, temp.image);

              for (int j=0; j<sizey; j++){
            for (int k=0; k<sizex; k++){
              int imgIndex = sizex*4*j + k*4;                   // index in the image 

              // Back to front compositing
              imgBuffer[imgIndex+0] = clamp((imgBuffer[imgIndex+0] * (1.0 - temp.image[imgIndex+3])) + temp.image[imgIndex+0]);
              imgBuffer[imgIndex+1] = clamp((imgBuffer[imgIndex+1] * (1.0 - temp.image[imgIndex+3])) + temp.image[imgIndex+1]);
              imgBuffer[imgIndex+2] = clamp((imgBuffer[imgIndex+2] * (1.0 - temp.image[imgIndex+3])) + temp.image[imgIndex+2]);
            }
          }

              delete []temp.image;
              temp.image = NULL;

              count++;
        }while( it!=depthPartitions.begin() );
      }

      if (tempRecvBuffer != NULL)
        delete []tempRecvBuffer;
      tempRecvBuffer = NULL;

      if (offsetBuffer != NULL)
        delete []offsetBuffer;
      offsetBuffer = NULL;

      if (recvSizePerProc != NULL)
        delete []recvSizePerProc;
      recvSizePerProc = NULL;
    }

    syncAllProcs();

  #endif
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

bool compareColor(code x, float r, float g, float b, float a){
  if ((x.color[0] == r && x.color[1]==g) &&  (x.color[2] == b && x.color[3]==a))
    return true;
  else
    return false;
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

code initCode(int count, float r, float g, float b, float a){
  code temp;
  temp.count = count;
  temp.color[0] = r;  temp.color[1]=g;  temp.color[2]=b;  temp.color[3]=a;
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

code incrCode(code x){
  x.count += 1;
  return x;
}



// ****************************************************************************
//  Method: avtImgCommunicator::rleEncodeAll
//
//  Purpose:
//    Encodes an image using RLE: Runlength, R, G, B, A
//
//  Arguments:
//    dimsX & dimsY : width & height of the image
//    imgArray    : array containing all the division composited by this processor
//    numDivs     : number of Z divisions
//
//    encoding      : the image compressed using RLE; it's a float so that it can be sent directly using MPI
//    sizeOfEncoding  : the compressed size of each image
//
//
//  Programmer: Pascal Grosset
//  Creation: August 23, 2013
//
//  Modifications:
//
// ****************************************************************************
int avtImgCommunicator::rleEncodeAll(int dimsX, int dimsY, int numDivs, float *imgArray,  float *& encoding, int *& sizeOfEncoding){
  std::vector<code> encodingVec;
  encodingVec.clear();
  code tempCode;
  sizeOfEncoding = new int[numDivs];
  int prev = 0;

  // Compress the data
  for (int j=0; j<numDivs; j++){
    int offset = dimsX*dimsY*4  *  j; // to get the next image in the array of images

    int i=0;
    tempCode = initCode(1, imgArray[offset + (i*4+0)],imgArray[offset + (i*4+1)],imgArray[offset + (i*4+2)],  imgArray[offset + (i*4+3)]);
    for (i=1; i<dimsX*dimsY; i++){
      if ( compareColor(tempCode, imgArray[offset + (i*4+0)],imgArray[offset + (i*4+1)],imgArray[offset + (i*4+2)],imgArray[offset + (i*4+3)]) )
        tempCode = incrCode(tempCode);
      else{
        encodingVec.push_back(tempCode);
        tempCode = initCode(1, imgArray[offset + (i*4+0)],imgArray[offset + (i*4+1)],imgArray[offset + (i*4+2)],imgArray[offset + (i*4+3)]);
      }
    }
    encodingVec.push_back(tempCode);

    if (j == 0)
      prev = sizeOfEncoding[j] = encodingVec.size();
    else{
      sizeOfEncoding[j] = encodingVec.size() - prev;
      prev = encodingVec.size();
    }

    debug5 <<  my_id << "  ~  encoding.size(): " << sizeOfEncoding[j] << "   offset: " << offset << "    original size: " << (dimsX * dimsY * 4) << "   size: " << encodingVec.size() << std::endl;
  }


  // Transfer the data to the encoding array
  int encSize = encodingVec.size();
  encoding = new float[encSize*5];
  
  int index = 0;
  for (int j=0; j<encSize; j++){
    encoding[index] = encodingVec[j].count; index++;
    encoding[index] = encodingVec[j].color[0];  index++;
    encoding[index] = encodingVec[j].color[1];  index++;
    encoding[index] = encodingVec[j].color[2];  index++;
    encoding[index] = encodingVec[j].color[3];  index++;  
  }
  encodingVec.clear();

  return encSize;   // size of the array
}



// ****************************************************************************
//  Method: avtImgCommunicator::rleDecode
//
//  Purpose:
//    Decodes rle encoded image to an RGBA image the size of the screen
//
//  Arguments:
//    encSize : size of compressed image
//    encoding: the image to be decoded
//    offset  : offset in the encoding array (all images from all procs are stored in 1 array)
//    img     : an array containing the RGBA decompressed image
//
//  Programmer: Pascal Grosset
//  Creation: August 23, 2013
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::rleDecode(int encSize, float *encoding, int offset, float *& img){
  int imgIndex = 0; // index into the image to be decompressed

  //debug5 << " 0 ~ offset: " << offset  << " encSize: " << encSize << std::endl;

  for (int i=0; i<encSize; i++)
    for (int j=0; j<(int)encoding[offset + i*5 + 0]; j++) // *5 to offset into image; +0 run count
    {
      img[imgIndex*4 + 0] = encoding[offset + i*5 + 1]; // R
      img[imgIndex*4 + 1] = encoding[offset + i*5 + 2]; // G
      img[imgIndex*4 + 2] = encoding[offset + i*5 + 3]; // B
      img[imgIndex*4 + 3] = encoding[offset + i*5 + 4]; // A

      imgIndex++;
    }
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
void avtImgCommunicator::initImage(int sizeX, int sizeY, float _color[4])
{
    imgBuffer = new float[sizeX * sizeY * 4];

    for (int _y=0; _y<sizeY; _y++)
        for (int _x=0; _x<sizeX; _x++)
        {
            int index = (_y*4*sizeX) + (_x*4);

            imgBuffer[index+0] = _color[0];
            imgBuffer[index+1] = _color[1];
            imgBuffer[index+2] = _color[2];
            imgBuffer[index+3] = _color[3];
        }
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
void avtImgCommunicator::ParallelDirectSend()
{

}



1447
1448
1449
1450
1451
1452
1453
1454
1455
1456
1457
1458
1459
1460
1461
1462
1463
1464
1465
1466
1467
1468
1469
1470
1471
1472
1473
1474
1475
1476
1477
1478
1479
1480
1481
1482
1483
1484
1485
1486
1487
1488
1489
1490
1491
1492
1493
1494
1495
1496
1497
1498
1499
1500
1501
1502
1503
1504
1505
1506
1507
1508
1509
1510
1511
1512
1513
1514
1515
1516
1517
1518
1519
1520
1521
1522
1523
1524
1525
1526
1527
1528
1529
1530
1531
1532
1533
1534
1535
1536
1537
1538
1539
1540
1541
1542
1543
1544
1545
1546
1547
1548
1549
1550
1551
1552
1553
1554
1555
1556
1557
1558
1559
1560
1561
1562
1563
1564
1565
1566
1567
1568
1569
1570
1571
1572
1573
1574
1575
1576
1577
1578
1579
1580
1581
1582
1583
1584
1585
1586
1587
1588
1589
1590
1591
1592
1593
1594
1595
1596
1597
1598
1599
1600
1601
1602
1603
1604
1605
1606
1607
1608
1609
1610
1611
1612
1613
1614
1615
1616
1617
1618
1619
1620
1621
1622
1623
1624
1625
1626
1627
1628
1629
1630
1631
1632
1633
1634
1635
1636
1637
1638
1639
1640
1641
1642
1643
1644
1645
1646
1647
1648
1649
1650
1651
1652
1653
1654
1655
1656
1657
1658
1659
1660
1661
1662
1663
1664
1665
1666
1667
1668
1669
1670
1671
1672
1673
1674
1675
1676
1677
1678
1679
1680
1681
1682
1683
1684
1685
1686
1687
1688
1689
1690
1691
1692
1693
1694
1695
1696
1697
1698
1699
1700
1701
1702
1703
1704
1705
1706
1707
1708
1709
1710
1711
1712
1713
1714
1715
1716
1717
1718
1719
1720
1721
1722
1723
1724
1725
1726
1727
1728
1729
1730
1731
1732
1733
1734
1735
1736
1737
1738
1739
1740
1741
1742
1743
1744
1745
1746
1747
1748
1749
1750
1751
1752
1753
1754
1755
1756
1757
1758
1759
1760
1761
1762
1763
1764
1765
1766
1767
1768
1769
1770
1771
1772
1773
1774
1775
1776
1777
1778
1779
1780
1781
1782
1783
1784
1785
1786
1787
1788
1789
1790
1791
1792
1793
1794
1795
1796
1797
1798
1799
1800
1801
1802
1803
1804
1805
1806
1807
1808
1809
1810
1811
1812
1813
1814
1815
1816
1817
1818
1819
1820
1821
1822
1823
1824
1825
1826
1827
1828
1829
1830
1831
1832
1833
1834
1835
1836
1837
1838
1839
1840
1841
1842
1843
1844
1845
1846
1847
1848
1849
1850
1851
1852
1853
1854
1855
1856
1857
1858
1859
1860
1861
1862
1863
1864
1865
1866
1867
1868
1869
1870
1871
1872
1873
1874
1875
1876
1877
1878
1879
1880
1881
1882
1883
1884
1885
1886
1887
1888
1889
1890
1891
1892
1893
1894
1895
1896
1897
1898
1899
1900
1901
1902
1903
1904
1905
1906
1907
1908
1909
1910
1911
1912
1913
1914
1915
1916
1917
1918
1919
1920
1921
1922
1923
1924
1925
1926
1927
1928
1929
1930
1931
1932
1933
1934
1935
1936
1937
1938
1939
1940
1941
1942
1943
1944
1945
1946
1947
1948
1949
1950
1951
1952
1953
1954
1955
1956
1957
1958
1959
1960
1961
1962
1963
1964
1965
1966
1967
1968
1969
1970
1971
1972
1973
1974
1975
1976
1977
1978
1979
1980
1981
1982
1983
1984
1985
1986
1987
1988
1989
1990
1991
1992
1993
1994
1995
1996
1997
1998
1999
2000
2001
2002
2003
2004
2005
2006
2007
2008
2009
2010
2011
2012
2013
2014
2015
2016
2017
2018
2019
2020
2021
2022
2023
2024
2025
2026
2027
2028
2029
2030
2031
2032
2033
2034
2035
2036
2037
2038
2039
2040
2041
2042
2043
2044
2045
2046
2047
2048
2049
2050
2051
2052
2053
2054
2055
2056
2057
2058
2059
2060
2061
2062
2063
2064
2065
2066
2067
2068
2069
2070
2071
2072
2073
2074
2075
2076
2077
2078
2079
2080
2081
2082
2083
2084
2085
2086
2087
2088
2089
2090
2091
2092
2093
2094
2095
2096
2097
2098
2099
2100
2101
2102
2103
2104
2105
2106
2107
2108
2109
2110
2111
2112
2113
2114
2115
2116
2117
2118
2119
2120
2121
2122
2123
2124
2125
2126
2127
2128
2129
2130
2131
2132
2133
2134
2135
2136
2137
2138
2139
2140
2141
2142
2143
2144
2145
2146
2147
2148
2149
2150
2151
2152
2153
2154
2155
2156
2157
2158
2159
2160
2161
2162
2163
2164
2165
2166
2167
2168
2169
2170
2171
2172
2173
2174
2175
2176
2177
2178
2179
2180
2181
2182
2183
2184
2185
2186
2187
2188
2189
2190
2191
2192
2193
2194
2195
2196
2197
2198
2199
2200
2201
2202
2203
2204
2205
2206
2207
2208
2209
2210
2211
2212
2213
2214
2215
2216
2217
2218
2219
2220
2221
2222
2223
2224
2225
2226
2227
2228
2229
2230
2231
2232
2233
2234
2235
2236
2237
2238
2239
2240
2241
2242
2243
2244
2245
2246
2247
2248
2249
2250
2251
2252
2253
2254
2255
2256
2257
2258
2259
2260
2261
2262
2263
2264
2265
2266
2267
2268
2269
2270
2271
2272
2273
2274
2275
2276
2277
2278
2279
2280
2281
2282
2283
2284
2285
2286
2287
2288
2289
2290
2291
2292
2293
2294
2295
2296
2297
2298
2299
2300
2301
2302
2303
2304
2305
2306
2307
2308
2309
2310
2311
2312
2313
2314
2315
2316
2317
2318
2319
2320
2321
2322
2323
2324
2325
2326
2327
2328
2329
2330
2331
2332
2333
2334
2335
2336
2337
2338
2339
2340
2341
2342
2343
2344
2345
2346
2347
2348
2349
2350
2351
2352
2353
2354
2355
2356
2357
2358
2359
2360
2361
2362
2363
2364
2365
2366
2367
2368
2369
2370
2371
2372
2373
2374
2375
2376
2377
2378
2379
2380
2381
2382
2383
2384
2385
2386
2387
2388
2389
2390
2391
2392
2393
2394
2395
2396
2397
2398
2399
2400
2401
2402
2403
2404
2405
2406
2407
2408
2409
2410
2411
2412
2413
2414
2415
2416
2417
2418
2419
2420
2421
2422
2423
2424
2425
2426
2427
2428
2429
2430
2431
2432
2433
2434
2435
2436
2437
2438
2439
2440
2441
2442
2443
2444
2445
2446
2447
2448
2449
2450
2451
2452
2453
2454
2455
2456
2457
2458
2459
2460
2461
2462
2463
2464
2465
2466
2467
2468
2469
2470
2471
2472
2473
2474
2475
2476
2477
2478
2479
2480
2481
2482
2483
2484
2485
2486
2487
2488
2489
2490
2491
2492
2493
2494
2495
2496
2497
2498
2499
2500
2501
2502
2503
2504
2505
2506
2507
2508
2509
2510
2511
2512
2513
2514
2515
2516
2517
2518
2519
2520
2521
2522
2523
2524
2525
2526
2527
2528
2529
2530
2531
2532
2533
2534
2535
2536
2537
2538
2539
2540
2541
2542
2543
2544
2545
2546
2547
2548
2549
2550
2551
2552
2553
2554
2555
2556
2557
2558
2559
2560
2561
2562
2563
2564
2565
2566
2567
2568
2569
2570
2571
2572
2573
2574
2575
2576
2577
2578
2579
2580
2581
2582
2583
2584
2585
2586
2587
2588
2589
2590
2591
2592
2593
2594
2595
2596
2597
2598
2599
2600
2601
2602
2603
2604
2605
2606
2607
2608
2609
2610
2611
2612
2613
2614
2615
2616
2617
2618
2619
2620
2621
2622
2623
2624
2625
2626
2627
2628
2629
2630
2631
2632
2633
2634
2635
2636
2637
2638
2639
2640
2641
2642
2643
2644
2645
2646
2647
2648
2649
2650
2651
2652
2653
2654
2655
2656
2657
2658
2659
2660
2661
2662
2663
2664
2665
2666
2667
2668
2669
2670
2671
2672
2673
2674
2675
2676
2677
2678
2679
2680
2681
2682
2683
2684
2685
2686
2687
2688
2689
2690
2691
2692
2693
2694
2695
2696
2697
2698
2699
2700
2701
2702
2703
2704
2705
2706
2707
2708
2709
2710
2711
2712
2713
2714
2715
2716
2717
2718
2719
2720
2721
2722
2723
2724
2725
2726
2727
2728
2729
2730
2731
2732
2733
2734
2735
2736
2737
2738
2739
2740
2741
2742
2743
2744
2745
2746
2747
2748
2749
2750
2751
2752
2753
2754
2755
2756
2757
2758
2759
2760
2761
2762
2763
2764
2765
2766
2767
2768
2769
2770
2771
2772
2773
2774
2775
2776
2777
2778
2779
2780
2781
2782
2783
2784
2785
2786
2787
2788
2789
2790
2791
2792
2793
2794
2795
2796
2797
2798
2799
2800
2801
2802
2803
2804
2805
2806
2807
2808
2809
2810
2811
2812
2813
2814
2815
2816
2817
2818
2819
2820
2821
2822
2823
2824
2825
2826
2827
2828
2829
2830
2831
2832
2833
2834
2835
2836
2837
2838
2839
2840
2841
2842
2843
2844
2845
2846
2847
2848
2849
2850
2851
2852
2853
2854
2855
2856
2857
2858
2859
2860
2861
2862
2863
2864
2865
2866
2867
2868
2869
2870
2871
2872
2873
2874
2875
2876
2877
2878
2879
2880
2881
2882
2883
2884
2885
2886
2887
2888
2889
2890
2891
2892
2893
2894
2895
2896
2897
2898
2899
2900
2901
2902
2903
2904
2905
2906
2907
2908
2909
2910
2911
2912
2913
2914
2915
2916
2917
2918
2919
2920
2921
2922
2923
2924
2925
2926
2927
2928
2929
2930
2931
2932
2933
2934
2935
2936
2937
2938
2939
2940
2941
2942
2943
2944
2945
2946
2947
2948
2949
2950
2951
2952
2953
2954
2955
2956
2957
2958
2959
2960
2961
2962
2963
2964
2965
2966
2967
2968
2969
2970
2971
2972
2973
2974
2975
2976
2977
2978
2979
2980
2981
2982
2983
2984
2985
2986
2987
2988
2989
2990
2991
2992
2993
2994
2995
2996
2997
2998
2999
3000
3001
3002
3003
3004
3005
3006
3007
3008
3009
3010
3011
3012
3013
3014
3015
3016
3017
3018
3019
3020
3021
3022
3023
3024
3025
3026
3027
3028
3029
3030
3031
3032
3033
3034
3035
3036
3037
3038
3039
3040
3041
3042
3043
3044
3045
3046
3047
3048
3049
3050
3051
3052
3053
3054
3055
3056
3057
3058
3059
3060
3061
3062
3063
3064
3065
3066
3067
3068
3069
3070
3071
3072
3073
3074
3075
3076
3077
3078
3079
3080
3081
3082
3083
3084
3085
3086
3087
3088
3089
3090
3091
3092
3093
3094
3095
3096
3097
3098
3099
3100
3101
3102
3103
3104
3105
3106
3107
3108
3109
3110
3111
3112
3113
3114
3115
3116
3117
3118
3119
3120
3121
3122
3123
3124
3125
3126
3127
3128
3129
3130
3131
3132
3133
3134
3135
3136
3137
3138
3139
3140
3141
3142
3143
3144
3145
3146
3147
3148
3149
3150
3151
3152
3153
3154
3155
3156
3157
3158
3159
3160
3161
3162
3163
3164
3165
3166
3167
3168
3169
3170
3171
3172
3173
3174
3175
3176
3177
3178
3179
3180
3181
3182
3183
3184
3185
3186
3187
3188
3189
3190
3191
3192
3193
3194
3195
3196
3197
3198
3199
3200
3201
3202
3203
3204
3205
3206
3207
3208
3209
3210
3211
3212
3213
3214
3215
3216
3217
3218
3219
3220
3221
3222
3223
3224
3225
3226
3227
3228
3229
3230
3231
3232
3233
3234
3235
3236
3237
3238
3239
3240
3241
3242
3243
3244
3245
3246
3247
3248
3249
3250
3251
3252
3253
3254
3255
3256
3257
3258
3259
3260
3261
3262
3263
3264
3265
3266
3267
3268
3269
3270
3271
3272
3273
3274
3275
3276
3277
3278
3279
3280
3281
3282
3283
3284
3285
3286
3287
3288
3289
3290
3291
3292
3293
3294
3295
3296
3297
3298
3299
3300
3301
3302
3303
3304
3305
3306
3307
3308
3309
3310
3311
3312
3313
3314
3315
3316
3317
3318
3319
3320
3321
3322
3323
3324
3325
3326
3327
3328
3329
3330
3331
3332
3333
3334
3335
3336
3337
3338
3339
3340
3341
3342
3343
3344
3345
3346
3347
3348
3349
3350
3351
3352
3353
3354
3355
3356
3357
3358
3359
3360
3361
3362
3363
3364
3365
3366
3367
3368
3369
3370
3371
3372
3373
3374
3375
3376
3377
3378
3379
3380
3381
3382
3383
3384
3385
3386
3387
3388
3389
3390
3391
3392
3393
3394
3395
3396
3397
3398
3399
3400
3401
3402
3403
3404
3405
3406
3407
3408
3409
3410
3411
3412
3413
3414
3415
3416
3417
3418
3419
3420
3421
3422
3423
3424
3425
3426
3427
3428
3429
3430
3431
3432
3433
3434
3435
3436
3437
3438
3439
3440
3441
3442
3443
3444
3445
3446
3447
3448
3449
3450
3451
3452
3453
3454
3455
3456
3457
3458
3459
3460
3461
3462
3463
3464
3465
3466
3467
3468
3469
3470
3471
3472
3473
3474
3475
3476
3477
3478
3479
3480
3481
3482
3483
3484
3485
3486
3487
3488
3489
3490
3491
3492
3493
3494
3495
3496
3497
3498
3499
3500
3501
3502
3503
3504
3505
3506
3507
3508
3509
3510
3511
3512
3513
3514
3515
3516
3517
3518
3519
3520
3521
3522
3523
3524
3525
3526
3527
3528
3529
3530
3531
3532
3533
3534
3535
3536
3537
3538
3539
3540
3541
3542
3543
3544
3545
3546
3547
3548
3549
3550
3551
3552
3553
3554
3555
3556
3557
3558
3559
3560
3561
3562
3563
3564
3565
3566
3567
3568
3569
3570
3571
3572
3573
3574
3575
3576
3577
3578
3579
3580
3581
3582
3583
3584
3585
3586
3587
3588
3589
3590
3591
3592
3593
3594
3595
3596
3597
3598
3599
3600
3601
3602
3603
3604
#ifndef _CPU_COMPOSITING_H_
#define _CPU_COMPOSITING_H_

#include <iostream>
#include <limits>       // std::numeric_limits

#include <vector>
#include <algorithm>
#include <map>

#include <mpi.h>
#include <omp.h>

#include <math.h>


#include "cpuBlending.h"
#include "../utils/timer.h"
#include "../utils/utils.h"
#include "../scheduler/scheduler.h"

#include "../scheduler/regionScheduler.h"


struct imgBuff{
    int partition;
    int extents[4];
    bool dataInBuffer;
    float *dataBuffer;
};

class CPUCompositing
{
    int myId;                   // mpi Id; 0 => root
    int numProcesses;           // mpi: numProcs

    bool hasData;
    bool compositingDone;
    bool renderingDone;
    float depth;                // z distance from the camera

    bool compositorOnSeperateRank;


    // Timing logs
    std::string outputTiming;
    std::string outputProfilingTiming;


    // Buffer
    std::vector<int> msgBuffer;         // (0/1, min-pos_x, min-pos_y, max-pos_x, max-pos_y) 0-not received, 1-received      
    float *dataBuffer;                  // instances of the data
    float *dataBufferStage1;            // instances of the data
    float *dataBufferRecv;
    float *alphaArray;

    

    Image bufferImg, recvImg, interimImg, interimImg1, interimImg2, currentImage;
    Image recvImage, sendImage, backupImg, newImage, compositedImg;   // temporary storages for other images


    bool compressionOn;
    int compositingStrategy;
    int numIterations;


    // the graph based load balancer scheduler
    Scheduler schCompositor;
    RegionScheduler schRegionCompositor;

  public:
    Image img;                               // img     - the data from rendering
    Image fullImg;                           // fullImg - the composited image
    std::vector<int> compositingOrder;       // sorted in z, from the camera to world's end (closest to farthest)


    CPUCompositing();
    ~CPUCompositing();

    void init(int _myRank, int _worldSize);


    //
    // Order needed for blending
    int exchangeDepthForPartition();
    void exchangeDepth();
    void gatherDepth();
    void updateBoundingBox(int currentBoundingBox[4], int imageExtents[4]);


    //
    // Getters and Setters
    void setId(int _id){ myId = _id; }
    void setNumProcesses(int _numProcesses){ numProcesses = _numProcesses; }
    void setDepth(float _depth){ depth = _depth; }
    void setHasData(bool _hasData){ hasData = _hasData; }
    void setCompositingDone(bool _compositingDone){ compositingDone = _compositingDone; }
    void setNumIterations(int _n){ numIterations = _n;}
    void setSeperateCompositor(bool _status){ compositorOnSeperateRank = _status; }

    int getId(){ return myId; }
    float getDepth(){ return depth; }
    bool getHasData(){ return hasData; }
    bool getCompositingDone(){ return compositingDone; }


    //
    // Log messages
    void clearTimingMsg(){ outputTiming = ""; }
    void addTimingMsg(std::string _msg){ outputTiming += _msg; }
    void outputTimingMessage(){ std::cout << outputTiming.c_str();}
    std::string getOutputTiming(){ return outputTiming; }

    //
    // Timing messages
    void clearProfilingTimingMsg(){ outputProfilingTiming = ""; }
    void addProfilingTimingMsg(std::string _msg){ outputProfilingTiming += _msg; }
    void outputProfilingTimingMessage(){ std::cout << myId << "," << outputProfilingTiming.c_str() << "-1\n";}
    std::string getOutputProfilingTiming(){ return outputProfilingTiming; }


    //
    // Sending
    void sendImgNoPackOpaqueFinal(int tags[2], int dest);
    void sendSImgNoPack(int tags[2], int dest, int offset, int extents[4], bool hasData);

    void gatherOpaqueImages(int regionGather[], int numToRecv, Image inputImg, Image outputImg, int tag, int width, int height);
    void gatherImages(int regionGather[], int numToRecv, Image inputImg, Image outputImg, int tag[2], int width, int height, float backgroundColor[4]);


    //
    // Compositing 
    //

    //
    // TOD-Tree
    void todTreePDS(int width, int height, int tag, Image img,                                                              // image input
                    int indexInLocalRegion, int r, int myRegionSize, int myRegionBB[4], std::vector<int> pdsProcesses,      // region information
                    MPI_Request *sendMetaRq, MPI_Request *sendImageRq, int *sendExtents,                                    // async buffers
                    Image &interimImg1);

    void todTree(Image _img, int r, int k, int tag[3], int width, int height, float backgroundColor[4]);


    //
    // Direct Send
    void serialDS(Image _img, int tags[2], float backgroundColor[4], int width, int height);    // Serial
    void parallelDS(Image _img, int region[], int numInRegion, int tags[3], float backgroundColor[4], int width, int height);  // Parallel


    //
    // Scheduler
    void schedulerInitGather(int worldSize, int width, int height, float imageInfo[5]);
    void launchScheduler(int tag[2], int width, int height, float backgroundColor[4]);
    void runCompositor(Image _img, float depth, int tags[3], float backgroundColor[4], int width, int height);

    //imgBuff *storage;
    void schedulerInitRegionGather(int numRegions, int worldSize, int width, int height, float imageInfo[5]);
    double launchRegionScheduler(int numRegions, int tag[2], int width, int height);
    double runRegionCompositor(Image _img, int tags[3], float backgroundColor[4], int numRegions, int width, int height);
    double runRegionCompositorASYNC(Image _img, int tags[3], float backgroundColor[4], int numRegions, int width, int height);

    bool validateExtents(int extents[4], int width, int height, std::string err);
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline CPUCompositing::CPUCompositing()
{
    compositorOnSeperateRank = false;

    dataBuffer = NULL;
    dataBufferStage1 = NULL;
    dataBufferRecv = NULL;
    alphaArray = NULL;
}


inline CPUCompositing::~CPUCompositing()
{
    #ifdef __INTEL_COMPILER
        if (dataBuffer != NULL)
            _mm_free(dataBuffer);

        if (dataBufferStage1 != NULL)
            _mm_free(dataBufferStage1);

        if (dataBufferRecv != NULL)
            _mm_free(dataBufferRecv);

        if (alphaArray != NULL)
            _mm_free(alphaArray);
    #else
        if (dataBuffer != NULL)
            delete []dataBuffer;

        if (dataBufferRecv != NULL)
            delete []dataBufferRecv;

        if (dataBufferStage1 != NULL)
            delete []dataBufferStage1;

        if (alphaArray != NULL)
            delete []alphaArray;
    #endif

    dataBuffer = NULL;
    dataBufferRecv = NULL;
    dataBufferStage1 = NULL;
    alphaArray = NULL;

    // //
    // // Local Cleanup
    // if (storage != NULL)
    //     delete []storage;
    // storage=NULL;
}





inline void CPUCompositing::init(int _myRank, int _worldSize)
{
    myId = _myRank;
    numProcesses = _worldSize;

    hasData = true;
    renderingDone = true;
    compositingDone = false;
    depth = -1;
    numIterations = 1;

    img.data = NULL;
    
    recvImage.data = NULL;
    sendImage.data = NULL;
    newImage.data = NULL;
    compositedImg.data = NULL;
    
    outputTiming = "";
    outputProfilingTiming = "";
}




void CPUCompositing::exchangeDepth()
{
    // Gather all the depths
    float *allDepths = new float[numProcesses];
    MPI_Allgather(&depth, 1, MPI_FLOAT, allDepths, 1, MPI_FLOAT, MPI_COMM_WORLD);


    // Sort the rank in compositingOrder
    std::multimap<float,int> depthsMap;
    for (int i=0; i<numProcesses; i++)
        //depthsMap.insert( std::pair<float,int>(1.0-allDepths[i],i) );
        depthsMap.insert( std::pair<float,int>(allDepths[i],i) );


    // Store the rank
    for (auto it=depthsMap.begin(); it!=depthsMap.end(); ++it)
        compositingOrder.push_back( (*it).second );


    delete []allDepths;
    allDepths = NULL;

    depthsMap.clear();
}


int CPUCompositing::exchangeDepthForPartition()
{
    // Gather all the depths
    float *allDepths = new float[numProcesses];
    MPI_Allgather(&depth, 1, MPI_FLOAT, allDepths, 1, MPI_FLOAT, MPI_COMM_WORLD);


    // Sort the rank in compositingOrder
    std::multimap<float,int> depthsMap;
    for (int i=0; i<numProcesses; i++)
        //depthsMap.insert( std::pair<float,int>(1.0-allDepths[i],i) );
        depthsMap.insert( std::pair<float,int>(allDepths[i],i) );


    // Store the rank
    for (auto it=depthsMap.begin(); it!=depthsMap.end(); ++it)
        compositingOrder.push_back( (*it).second );


    delete []allDepths;
    allDepths = NULL;


    int rank = (std::find(compositingOrder.begin(), compositingOrder.end(), myId)) - compositingOrder.begin();

    depthsMap.clear();
    compositingOrder.clear();


    return rank;
}


void CPUCompositing::gatherDepth()
{
    // Gather all the depths
    float *allDepths = new float[numProcesses];
    MPI_Gather(&depth, 1, MPI_FLOAT, allDepths, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);


    // Sort the rank in compositingOrder
    std::multimap<float,int> depthsMap;
    for (int i=1; i<numProcesses; i++)
        //depthsMap.insert( std::pair<float,int>(1.0-allDepths[i],i) );
        depthsMap.insert( std::pair<float,int>(allDepths[i],i) );


    // Store the rank
    for (auto it=depthsMap.begin(); it!=depthsMap.end(); ++it)
        compositingOrder.push_back( (*it).second );


    delete []allDepths;
    allDepths = NULL;

    depthsMap.clear();
}


void CPUCompositing::updateBoundingBox(int currentBoundingBox[4], int imageExtents[4])
{
    if ( (currentBoundingBox[0] == 0 && currentBoundingBox[1] == 0) && (currentBoundingBox[2] == 0 && currentBoundingBox[3] == 0)){
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



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void CPUCompositing::sendSImgNoPack(int tags[2], int dest, int offset, int extents[4], bool hasData)
{
    Timer send;
    MPI_Request myRequest;


  send.start();
    MPI_Send(extents, 4, MPI_INT, dest, tags[0], MPI_COMM_WORLD);

    if (hasData)
        MPI_Send(&sendImage.data[offset], (extents[3]-extents[2])*(extents[1]-extents[0])*4, MPI_FLOAT, dest, tags[1], MPI_COMM_WORLD);    // Send the actual image
  send.stop();

    addTimingMsg( toString(myId) + " ~  Elapsed time sending to " + toString(dest) + " : " + toString(extents[1]-extents[0]) + " x " + toString(extents[3]-extents[2]) + " : " + toString(send.getDuration()) + "\n" );
    addProfilingTimingMsg("e|to " + toString(dest) + " ~ "  + toString(extents[1]-extents[0]) + " x " + toString(extents[3]-extents[2]) + ":" + toString(send.getDuration()) + ",");
}


inline void CPUCompositing::sendImgNoPackOpaqueFinal(int tags[2], int dest)
{
    Timer send;
    MPI_Request myRequest;

  send.start();
        int imgSize = (sendImage.extents[1]-sendImage.extents[0])*(sendImage.extents[3]-sendImage.extents[2]);
        MPI_Send(sendImage.data, imgSize*3, MPI_FLOAT, dest, tags[1], MPI_COMM_WORLD);    // Send the actual image
  send.stop();

    addTimingMsg( toString(myId) + " ~  Elapsed time sending to " + toString(dest) + " : " + toString(sendImage.extents[1]-sendImage.extents[0]) + " x " + toString(sendImage.extents[3]-sendImage.extents[2]) + " : " + toString(send.getDuration()) + "\n" );
    addProfilingTimingMsg("v|to " + toString(dest) + " ~ "  + toString(sendImage.extents[1]-sendImage.extents[0]) + " x " + toString(sendImage.extents[3]-sendImage.extents[2]) + ":" + toString(send.getDuration()) + ",");
}



// Gathering
inline void CPUCompositing::gatherOpaqueImages(int regionGather[], int numToRecv, Image inputImg, Image outputImg, int tag, int width, int height)
{
    if (myId == 0)
    {
        //
        // Receive at root/display node!
        Timer recvTimer, blendTimer, setupTimer;
      setupTimer.start();
      
        fullImg.createOpaqueImage(0,width, 0,height);

        int regionHeight     = height / numToRecv;
        int lastRegionHeight = height - regionHeight*(numToRecv-1);
        int lastOpaqueBufferSize = lastRegionHeight * width * 3; 
        int opaqueRegularBuffer  = regionHeight * width * 3;  

        for (int i=0; i<numToRecv; i++)
            if (regionGather[i] == myId){
                numToRecv--;
                break;
            }
        
        //
        // Create buffers for async reciving
        MPI_Request *recvImageRq = new MPI_Request[ numToRecv ];
        MPI_Status  *recvImageSt = new MPI_Status[ numToRecv ];
      setupTimer.stop();

      recvTimer.start();

        // Async Recv
        int recvCount=0;
        for (int i=0; i<numToRecv; i++){
            int src = regionGather[i];
            if (src == myId)
                continue;

            if (i == numToRecv-1)
                MPI_Irecv(&fullImg.data[i*opaqueRegularBuffer], lastOpaqueBufferSize, MPI_FLOAT, src, tag, MPI_COMM_WORLD,  &recvImageRq[recvCount] ); 
            else
                MPI_Irecv(&fullImg.data[i*opaqueRegularBuffer], opaqueRegularBuffer, MPI_FLOAT, src, tag, MPI_COMM_WORLD,  &recvImageRq[recvCount] );
            recvCount++;
        }

      
        if (compositingDone == false)   // If root has data for the final image
        {
          blendTimer.start();
            fullImg.placeInOpaqueImage(inputImg.data, inputImg.extents, inputImg.extents[0], inputImg.extents[1], inputImg.extents[2], inputImg.extents[3]);
          blendTimer.stop();
          
            addProfilingTimingMsg("z| from " + toString(myId) + " ~ " + toString(inputImg.extents[1]-inputImg.extents[0]) + "x"  + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
            addTimingMsg(toString(myId) + " ~ Final Blending " + toString(inputImg.extents[1]-inputImg.extents[0]) + "x"  + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(blendTimer.getDuration()) + "\n");            

            inputImg.deleteImage();
        }

        MPI_Waitall(numToRecv, recvImageRq, recvImageSt);
      recvTimer.stop();

        addProfilingTimingMsg("s:" + toString(setupTimer.getDuration()) + ",");
        addProfilingTimingMsg("u:" + toString(recvTimer.getDuration()) + ",");

        addTimingMsg(toString(myId) + " ~ Final setup " + toString(setupTimer.getDuration()) + "\n");
        addTimingMsg(toString(myId) + " ~ Final recv " + toString(recvTimer.getDuration()) + "\n");


        delete []recvImageRq;
        recvImageRq = NULL;
        delete []recvImageSt;
        recvImageSt = NULL;
    }
    else
    {
        if (compositingDone == false)   // If root has data for the final image
        {
            Timer sendTimer;

          sendTimer.start();
            MPI_Send(inputImg.data, inputImg.getNumPixels()*inputImg.getNumComponents(), MPI_FLOAT, 0, tag, MPI_COMM_WORLD);    // Send the actual image
          sendTimer.stop();

            addTimingMsg( toString(myId) + " ~  Elapsed time sending to 0 : " + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + " : " + toString(sendTimer.getDuration()) + "\n" );
            addProfilingTimingMsg("v|to 0 ~ "  + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(sendTimer.getDuration()) + ",");

            inputImg.deleteImage();
            compositingDone = true;
        }
    }
}


inline void CPUCompositing::gatherImages(int regionGather[], int numToRecv, Image inputImg, Image outputImg, int tag[2], int width, int height, float backgroundColor[4])
{
    if (myId == 0)  // Display process
    {
        //
        // Receive at root/display node!
        Timer recvTimer, blendTimer, bkgTimer, setupTimer;

        int *msgBuffer = new int[numToRecv*5];
        int regionHeight = height/numToRecv;
        int lastRegionHeight = height - regionHeight*(numToRecv-1);
        if ( lastRegionHeight > regionHeight)
            regionHeight = lastRegionHeight;
        int bufferSize = width * regionHeight * 4;

        Image tempImg;
        float *gatherBuffer = NULL;

        #ifdef __INTEL_COMPILER
            gatherBuffer = (float *)_mm_malloc( bufferSize * numToRecv * sizeof(float), _ALIGN_);
        #else
            gatherBuffer = new float[ bufferSize * numToRecv];
        #endif


        // Adjust numtoRecv
        for (int i=0; i<numToRecv; i++)
            if (regionGather[i] == myId){
                numToRecv--;
                break;
            }
        
        //
        // Create buffers for async reciving
        MPI_Request *recvMetaRq = new MPI_Request[ numToRecv ];
        MPI_Request *recvImageRq = new MPI_Request[ numToRecv ];

        MPI_Status *recvMetaSt = new MPI_Status[ numToRecv ];
        MPI_Status *recvImageSt = new MPI_Status[ numToRecv ];


        // Async Recv
        int recvCount=0;
        for (int i=0; i<numToRecv; i++)
        {
            int src = regionGather[i];
            if (src == myId)
                continue;

            MPI_Irecv(&msgBuffer[i*5],                         5*sizeof(int),   MPI_INT, src, tag[0], MPI_COMM_WORLD,  &recvMetaRq[recvCount] );
            MPI_Irecv(&gatherBuffer[i*bufferSize],  bufferSize*sizeof(float), MPI_FLOAT, src, tag[1], MPI_COMM_WORLD,  &recvImageRq[recvCount] );
            recvCount++;
        }



        // create image with background
      bkgTimer.start();
        outputImg.colorImage(backgroundColor[0], backgroundColor[1], backgroundColor[2], backgroundColor[3]);
      bkgTimer.stop();
        addProfilingTimingMsg("b| with bkg ~ "  + toString(width) + " x " + toString(height) + ":" + toString(bkgTimer.getDuration()) + ",");


        // If root has data for the final image
        if (compositingDone == false)
        {
          blendTimer.start();
            blendInto(FRONT_TO_BACK, outputImg, inputImg);
          blendTimer.stop();

            //outputImg.outputPPM("output/blendLocal_" + toString(myId));
            addProfilingTimingMsg("z| local ~ "  + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
        }

       
        int recvBlend = 0;
        for (int i=0; i<numToRecv; i++)
        {
            MPI_Wait(&recvMetaRq[recvBlend], &recvMetaSt[recvBlend]);

            if ((msgBuffer[recvBlend*5+1]-msgBuffer[recvBlend*5+0]) * (msgBuffer[recvBlend*5+3]-msgBuffer[recvBlend*5+2]))
            {
                // Has Data!
              recvTimer.start();
                MPI_Wait(&recvImageRq[recvBlend], &recvImageSt[recvBlend]);
              recvTimer.stop();
                
                // Blend
              blendTimer.start();
                tempImg.extents[0] = msgBuffer[recvBlend*5 + 0];
                tempImg.extents[1] = msgBuffer[recvBlend*5 + 1];
                tempImg.extents[2] = msgBuffer[recvBlend*5 + 2];
                tempImg.extents[3] = msgBuffer[recvBlend*5 + 3];
                tempImg.data = &gatherBuffer[recvBlend*bufferSize];

                blendInto(FRONT_TO_BACK, outputImg, tempImg);

              blendTimer.stop();

                addProfilingTimingMsg("f| " + toString(i) + " :" + toString(recvTimer.getDuration()) + ",");
                addProfilingTimingMsg("z| " + toString(i) + " ~ " + toString(tempImg.extents[1]-tempImg.extents[0]) + " x " + toString(tempImg.extents[3]-tempImg.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
            }

            recvBlend++;
        }

 
         #ifdef __INTEL_COMPILER
            _mm_free(gatherBuffer);
            gatherBuffer = NULL;
        #else
            delete []gatherBuffer;
            gatherBuffer = NULL;
        #endif

        if (recvMetaRq != NULL){
            delete []recvMetaRq;
            recvMetaRq = NULL;
        }

        if (recvImageRq != NULL){
           delete []recvImageRq;
           recvImageRq = NULL;
        }

        if (recvMetaSt != NULL){
            delete []recvMetaSt;
            recvMetaSt = NULL;
        }

        if (recvImageSt != NULL){
           delete []recvImageSt;
           recvImageSt = NULL;
        }

        delete []msgBuffer;
        msgBuffer = NULL;
    }
    else
    {
        if (compositingDone == false)   // If root has data for the final image
        {
            Timer sendImageTimer;

          sendImageTimer.start();
            int sendParams[5];
            for (int i=0; i<4; i++)
                sendParams[i] = inputImg.extents[i];
            sendParams[4]=0;

            MPI_Send(sendParams,                                                      5, MPI_INT,   0, tag[0], MPI_COMM_WORLD);    // Send meta data
            MPI_Send(inputImg.data, inputImg.getNumPixels()*inputImg.getNumComponents(), MPI_FLOAT, 0, tag[1], MPI_COMM_WORLD);    // Send the actual image
          sendImageTimer.stop();

            addTimingMsg(toString(myId) +" ~ Gather Send ~ " + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(sendImageTimer.getDuration()) + "\n");
            addProfilingTimingMsg("v| " + toString(myId) + " ~ " + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(sendImageTimer.getDuration()) + ",");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//
// TOD-Tree
//


void CPUCompositing::todTreePDS(int width, int height, int tag, Image img,                                                              // image input
                                int indexInLocalRegion, int r, int myRegionSize, int myRegionBB[4], std::vector<int> pdsProcesses,      // region information
                                MPI_Request *sendMetaRq, MPI_Request *sendImageRq, int *sendExtents,                                    // async buffers
                                Image &interimImg1)                                                                                      // output image
{
    Timer blendTimer, prepareTimer, waitTimer, recvTimer, otherTimer;

  prepareTimer.start();

    int localExTags[2];
    localExTags[0] = tag;
    localExTags[1] = tag+1;

    MPI_Request *recvMetaRq = new MPI_Request[ myRegionSize-1 ];
    MPI_Request *recvImageRq = new MPI_Request[ myRegionSize-1 ];

    MPI_Status *recvMetaSt = new MPI_Status[ myRegionSize-1 ];
    MPI_Status *recvImageSt = new MPI_Status[ myRegionSize-1 ];

    // std::cout << myId << " ~ indexInLocalRegion: " << indexInLocalRegion << "  r: " << r <<  "  myRegionSize: " << myRegionSize << "   myRegionBB: " << myRegionBB[0] << ", " << myRegionBB[1] << "   " << myRegionBB[2] << ", " << myRegionBB[3] << std::endl;
    // std::cout << "Rank: ";
    // for (int i=0; i<pdsProcesses.size(); i++)
    //     std::cout << pdsProcesses[i] << " ";
    // std::cout << "\n" << std::endl;
    // MPI_Barrier(MPI_COMM_WORLD);

    // 
    // Determine if the process will be receiving data and recv order
    bool inRegion = true;
    if (indexInLocalRegion >= r)
        inRegion = false;

    bool firstHalf = true;
    if (indexInLocalRegion >= r/2)
        firstHalf = false;

    // size of buffer needed
    int sizeOneBuffer = (myRegionBB[1]-myRegionBB[0])*(myRegionBB[3]-myRegionBB[2])*4;
    int regionHeight = height/r;


    //
    // Create buffer for recv
    #ifdef __INTEL_COMPILER
        dataBufferStage1 = (float *)_mm_malloc( sizeOneBuffer * myRegionSize * sizeof(float), _ALIGN_);
    #else
        dataBufferStage1 = new float[ sizeOneBuffer * myRegionSize];
    #endif
    msgBuffer.clear();
    msgBuffer.resize(5 * myRegionSize); 

    //
    // Async Recv
    if (inRegion)
    {
        int recvCount=0;
        for (int i=0; i<myRegionSize; i++)
        {
            int index = i;
            if (!firstHalf)
                index = (myRegionSize-1)-i;

            if ( pdsProcesses[index] == myId )
                continue;

            int src = pdsProcesses[index];
            MPI_Irecv(&msgBuffer[index*5],                                5,   MPI_INT, src, localExTags[0], MPI_COMM_WORLD,  &recvMetaRq[recvCount] );
            MPI_Irecv(&dataBufferStage1[index*sizeOneBuffer], sizeOneBuffer, MPI_FLOAT, src, localExTags[1], MPI_COMM_WORLD,  &recvImageRq[recvCount] );
            recvCount++;
        }
    }


    //
    // Async Send
    int sendCount = 0;
    int sendingOffset;
    for (int i=0; i<r; i++)
    {
        int regionStart, regionEnd, imgSize, index, dest;
        index = i;

        if (!firstHalf)
            index = (r-1)-i;
        dest = pdsProcesses[index];

        if ( dest == myId )
            continue;

        regionStart = index*regionHeight;
        regionEnd = regionStart + regionHeight;
        if (index == r-1) // the last one in region
            regionEnd = height;


        if (regionStart < img.extents[2])
            regionStart = img.extents[2];

        if (regionEnd > img.extents[3])
            regionEnd = img.extents[3];


        bool hasData = true;
        if (regionEnd - regionStart <= 0 || img.extents[1]-img.extents[0] <= 0){
            hasData = false;

            sendingOffset = 0;
            imgSize = sendExtents[index*5 + 0] = sendExtents[index*5 + 1] = sendExtents[index*5 + 2] = sendExtents[index*5 + 3] = sendExtents[index*5 + 4] =  0;
            //std::cout << myId << " ~  sending to  " << pdsProcesses[index] << " ... NO DATA!!!!" << std::endl;
        }
        else
        {
            imgSize = (regionEnd-regionStart) * (img.extents[1]-img.extents[0]) * 4;
            sendingOffset = (regionStart-img.extents[2]) * (img.extents[1]-img.extents[0]) * 4;

            sendExtents[index*5 + 0] = img.extents[0];
            sendExtents[index*5 + 1] = img.extents[1];
            sendExtents[index*5 + 2] = regionStart;
            sendExtents[index*5 + 3] = regionEnd;
            sendExtents[index*5 + 4] = 0;

            //std::cout << myId << " ~  sending to  " << pdsProcesses[index] << " ... HAS DATA!!!!" << std::endl;
        }

        // std::cout << myId << " ~ index: " << index << "   pdsProcesses[index]: " << pdsProcesses[index] << "  extents: " <<  sendExtents[index*5 + 0] << ", " << sendExtents[index*5 + 1]  << ", " << sendExtents[index*5 + 2] << ", " << sendExtents[index*5 + 3] << ", " << sendExtents[index*5 + 4] << "  sending ... " << std::endl;
        MPI_Isend(&sendExtents[index*5],          5,   MPI_INT, dest, localExTags[0], MPI_COMM_WORLD, &sendMetaRq[sendCount]);
        MPI_Isend(&img.data[sendingOffset], imgSize, MPI_FLOAT, dest, localExTags[1], MPI_COMM_WORLD, &sendImageRq[sendCount]);
        
        sendCount++;
    }


    //
    // Create buffer for region
    interimImg1.initializeZero();
    interimImg1.setBB(2*width,-2*width, 2*height,-2*height);

  prepareTimer.stop();
    addProfilingTimingMsg("p| from " + toString(myId) + ":" + toString(prepareTimer.getDuration()) + ",");
    addTimingMsg( toString(myId) + " ~ Prepare " + toString(myId) +  " : " + toString(prepareTimer.getDuration()) + "\n" );



    //
    // Blend
    int numBlends = 0;
    int countBlend = 0;
    if (inRegion)
    {
        for (int i=0; i<r; i++)
        {
            int index;
            if (firstHalf)
                index = i;
            else
                index = (r-1)-i;

            if (pdsProcesses[index] == myId)
            {
              blendTimer.start();

                int regionStart = myRegionBB[2];
                int regionEnd = myRegionBB[3];

                if (regionStart < img.extents[2])
                    regionStart = img.extents[2];

                if (regionEnd > img.extents[3])
                    regionEnd = img.extents[3];


                bool hasData = true;
                if (regionEnd - regionStart <= 0){
                    hasData = false;
                    regionEnd = regionStart = 0;
                }

                if (hasData == true)
                {
                    int extentsSectionRecv[4];
                    extentsSectionRecv[0] = img.extents[0];
                    extentsSectionRecv[1] = img.extents[1];
                    extentsSectionRecv[2] = regionStart;
                    extentsSectionRecv[3] = regionEnd;

                    if (firstHalf)
                        blendInto(FRONT_TO_BACK, interimImg1, img, extentsSectionRecv);
                    else
                        blendInto(BACK_TO_FRONT, interimImg1, img, extentsSectionRecv);
                    updateBoundingBox(interimImg1.boundingBox, extentsSectionRecv);

                    //interimImg1.outputPPM("output/localBlend_" + toString(myId));
                    numBlends++;
                }
              blendTimer.stop();

                addProfilingTimingMsg("b| from " + toString(myId) + " ~ " + toString(img.extents[1]-img.extents[0]) + "x"  + toString(img.extents[3]-img.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
                addTimingMsg( toString(myId) + " ~ Blending " + toString(myId) + " size: " + toString(img.extents[1]-img.extents[0]) + "x"  + toString(img.extents[3]-img.extents[2]) + " : " + toString(blendTimer.getDuration()) + "\n" ); 
            }
            else
            {
              recvTimer.start();
                MPI_Wait(&recvMetaRq[countBlend], &recvMetaSt[countBlend]);

                for (int j=0; j<4; j++)
                    recvImage.extents[j] = msgBuffer[index*5 + j];


                bool hasData =  false;
                if (recvImage.extents[1]-recvImage.extents[0] > 0 && recvImage.extents[3]-recvImage.extents[2] > 0){
                    hasData = true;
                    MPI_Wait(&recvImageRq[countBlend], &recvImageSt[countBlend]);
                    recvImage.data = &dataBufferStage1[index*sizeOneBuffer];
                }
                // else
                // {
                //     std::cout << myId << " ~ index: " << index << "  recv has no data!" << std::endl;
                // }
                
              recvTimer.stop();
                
              blendTimer.start();
                if (hasData)
                {
                    //recvImage.outputPPM("output/recv_" + toString(myId));
                    if (firstHalf)
                        blendInto(FRONT_TO_BACK, interimImg1, recvImage);
                    else
                        blendInto(BACK_TO_FRONT, interimImg1, recvImage);
                    updateBoundingBox(interimImg1.boundingBox, recvImage.extents);
                    //interimImg1.outputPPM("output/localRecv_" + toString(myId));
                    
                    numBlends++;
                }
              blendTimer.stop();
                countBlend++;

                addProfilingTimingMsg("f| from " + toString(pdsProcesses[index]) + " ~ " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + ":" + toString(recvTimer.getDuration()) + ",");
                addTimingMsg( toString(myId) + " ~ Receiving " + toString(pdsProcesses[index]) + " size: " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + " : " + toString(recvTimer.getDuration()) + "\n" );

                addProfilingTimingMsg("b| from " + toString(pdsProcesses[index]) + " ~ " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
                addTimingMsg( toString(myId) + " ~ Blending " + toString(pdsProcesses[index]) + " size: " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + " : " + toString(blendTimer.getDuration()) + "\n" );
            }
        }
    }
    else
        compositingDone = true;

    // std::cout <<  myId << " ~ Blend " << std::endl;
    // MPI_Barrier(MPI_COMM_WORLD);

  waitTimer.start();
  
    msgBuffer.clear();
    #ifdef __INTEL_COMPILER
        _mm_free(dataBufferStage1);
        dataBufferStage1 = NULL;
    #else
        delete []dataBufferStage1;
        dataBufferStage1 = NULL;
    #endif

    if (numBlends == 0)
        interimImg1.boundingBox[0]=interimImg1.boundingBox[1]=interimImg1.boundingBox[2]=interimImg1.boundingBox[3] = 0;

    delete []recvMetaRq;
    recvMetaRq = NULL;
    delete []recvImageRq;
    recvImageRq = NULL;
    delete []recvMetaSt;
    recvMetaSt = NULL;
    delete []recvImageSt;
    recvImageSt = NULL;


  waitTimer.stop();
    addProfilingTimingMsg("i| from " + toString(myId) + ":" + toString(waitTimer.getDuration()) + ",");
    addTimingMsg( toString(myId) + " ~ Other wait " + toString(myId) +  " : " + toString(waitTimer.getDuration()) + "\n" );
}


// numProcessesInRegion = r
// groupSize = k
void CPUCompositing::todTree(Image _img, int r, int k, int tag[3], int width, int height, float backgroundColor[4])
{
    Timer overallTimer;
    Timer send, recv, wait, blend, other, prepare;
    Timer setupTimer, sendTimer,  recvTimer, waitTimer, blendTimer, otherTimer, prepareTimer;
    //fullImg.createImage(0,width, 0,height);

    for (int iter=0; iter<numIterations; iter++)
    {
      overallTimer.start();
      setupTimer.start();
        
        compositingDone = false;

        // Default for even regions
        int numRegions = numProcesses / r;              // Get the number of regions
        int lastRegionID = numRegions - 1;              // Determine the ID of last region
        int myRegionSize = r;                           // Typical size of each region
        bool extraneousNodes = false;                   //   


        // Find position among all ranks
        std::vector<int>::iterator it = std::find(compositingOrder.begin(), compositingOrder.end(), myId);
        int myIndex = it - compositingOrder.begin();    // Index among all processes      
        
        // Find position in region
        int myRegionIndex = myIndex / r;                // Index of the region i'm in
        int indexInLocalRegion  = myIndex % r;          // Index in region for Parallel Direct Send
        int _indexInLocalRegion = myIndex % r;          //  Index in region for Parallel Direct Send (temp to avoid if)

        if (numProcesses%r != 0)                        // Adjust for last region of uneven number of processes
        {
            if ( myRegionIndex >= lastRegionID )        // Adjust region size if in last n before to last
                myRegionSize = numProcesses - lastRegionID * r;

            if ( myRegionIndex > lastRegionID )         // Adjust region index + position if in last region
            {
                myRegionIndex = lastRegionID;
                indexInLocalRegion = indexInLocalRegion + r;
                extraneousNodes = true;
            }      
        }
      
        

        //
        // Compute local regions
        std::vector<int> pdsProcesses;                  // Processes for parallel Direct Send
        std::vector<int> karyProcesses;                 // Processes for k-ary compositing
        std::vector<int> gatherProcesses;               // Processes for gather

        // Processes for parallel Direct Send
        for (int i=0; i<myRegionSize; i++)
            pdsProcesses.push_back( compositingOrder[ myRegionIndex*r + i] );
        
        // Processes for k-ary compositing
        for (int i=0; i<numRegions; i++)
            karyProcesses.push_back( compositingOrder[_indexInLocalRegion + i*r] );

        // Processes for gather
        for (int i=0; i<r; i++)
            gatherProcesses.push_back( compositingOrder[i] );



        // 
        // Compute local region dimensions
        int regionHeight = height / r;                           // Generic ( 501/4 = 125.25 = 125; we take the min )
        int lastRegionHeight = height - regionHeight*(r-1); 

        int myStartingHeight = _indexInLocalRegion * regionHeight;   
        int myEndingHeight   = myStartingHeight + regionHeight;       

        if (_indexInLocalRegion == r-1) // the last one in region
            myEndingHeight = height;

        int myRegionHeight = myEndingHeight-myStartingHeight;


        //
        // Buffer Sizes
        int myBufferSize = myRegionHeight * width * 4;

        int myRegionBB[4];
        myRegionBB[0] = 0;                myRegionBB[1] = width;
        myRegionBB[2] = myStartingHeight;   myRegionBB[3] = myEndingHeight;

        

        // Gather opaque buffers
        int lastOpaqueBufferSize = lastRegionHeight * width * 3; 
        int opaqueRegularBuffer  = regionHeight * width * 3;  

    
        // Buffer Size
        interimImg1.boundingBox[0] =  width*2;
        interimImg1.boundingBox[1] = -width*2;
        interimImg1.boundingBox[2] =  height*2;
        interimImg1.boundingBox[3] = -height*2;


        // Async Send Buffers
        int *sendExtents;;
        MPI_Request *sendMetaRq, *sendImageRq;
        MPI_Status *sendMetaSt, *sendImageSt;


        if (extraneousNodes)
        {
            // Send and do not receive
            sendMetaRq = new MPI_Request[ r ];
            sendImageRq = new MPI_Request[ r ];

            sendMetaSt = new MPI_Status[ r ];
            sendImageSt = new MPI_Status[ r ];

            sendExtents = new int[r * 5];
        }
        else
        {
            sendMetaRq = new MPI_Request[ r-1 ];
            sendImageRq = new MPI_Request[ r-1 ];

            sendMetaSt = new MPI_Status[ r-1 ];
            sendImageSt = new MPI_Status[ r-1 ];

            sendExtents = new int[r * 5];
        }

      setupTimer.stop();
        addProfilingTimingMsg("o| from " + toString(myId) + ":" + toString(setupTimer.getDuration()) + ",");
        addTimingMsg( toString(myId) + " ~ initialization other Prepare " + toString(myId) +  " : " + toString(setupTimer.getDuration()) + "\n" );




        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // Stage 1 - exchange to make one region authoritative
        //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        //
        // Create output image
        interimImg1.createImage(0, width, myRegionBB[2], myRegionBB[3]);

        todTreePDS(width, height, tag[0], _img,                                         // image input
                   indexInLocalRegion, r, myRegionSize, myRegionBB, pdsProcesses,       // region information
                   sendMetaRq, sendImageRq, sendExtents,                                // async buffers
                   interimImg1);                                                        // output image

        // // Debug
        // MPI_Barrier(MPI_COMM_WORLD);    
        // std::cout << myId << " ~ done PDS: Extents: " << interimImg1.extents[0]     <<", " << interimImg1.extents[1] << "   " << interimImg1.extents[2] << ", " << interimImg1.extents[3] <<
        //                                      "  BB: " << interimImg1.boundingBox[0] <<", " << interimImg1.boundingBox[1] << "   " << interimImg1.boundingBox[2] << ", " << interimImg1.boundingBox[3] << std::endl;
        // interimImg1.outputPPM("output/afterPDS_" + toString(myId));

        // std::cout << myId << " ~  : " << karyProcesses[0] << ", " << karyProcesses[1] <<  ", " << karyProcesses[2] <<  ", " << karyProcesses[3] << std::endl;
        // MPI_Barrier(MPI_COMM_WORLD);    


        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // Stage 2 - tree like compositing
        //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        

      otherTimer.start();
        int rounds = -1;
        int localExTags[2];
        float *alphaArray = NULL;
        bool lastRound = false;

        
        //
        // Create buffer for recv only for those who will need it
        std::vector<int>::iterator _it = std::find(karyProcesses.begin(), karyProcesses.end(), myId);
        int _myKIndex = _it - karyProcesses.begin();    // myIndex = my position in the array
        
        if (_myKIndex % k == 0)       
        {
            #ifdef __INTEL_COMPILER
                dataBuffer = (float *)_mm_malloc( myBufferSize * k * sizeof(float), _ALIGN_);
            #else
                dataBuffer = new float[ myBufferSize * k];
            #endif
            msgBuffer.resize(5 * k);
        }
       
      otherTimer.stop();

        addTimingMsg( toString(myId) + " ~ Elapsed time preparing: " + toString(otherTimer.getDuration()) + "\n" );
        addProfilingTimingMsg("o:" + toString(otherTimer.getDuration()) + ",");


        while (karyProcesses.size() > 1)
        {
            if (compositingDone == false)
            {
              prepareTimer.start();
                rounds++;
                localExTags[0] = tag[1] + 2*rounds;         localExTags[1] = tag[1]+ 2*rounds + 1;

                std::vector<int>::iterator it = std::find(karyProcesses.begin(), karyProcesses.end(), myId);
                int myKIndex = it - karyProcesses.begin();    // myKIndex = my position in the array

                
                //
                // Sending
                //
                if (myKIndex % k != 0 )
                {
                    int destProc = karyProcesses[ ( myKIndex/k )  *  k ];

                    //std::cout << myId << " ~ interimImg1.boundingBox: " << interimImg1.boundingBox[0] << ", " <<interimImg1.boundingBox[1] << "   " << interimImg1.boundingBox[2] << ", " << interimImg1.boundingBox[3] << "   my extents height: " << myStartingHeight << ", " << endingHeight << std::endl;

                    int sendingOffset = (interimImg1.boundingBox[2] - interimImg1.extents[2]) * (interimImg1.getWidth()) * 4;
                    bool hasData = true;
                    int sendParams[5];

                    if (interimImg1.getBBwidth() <= 0 || interimImg1.getBBHeight() <= 0){
                        hasData = false;
                        sendingOffset = 0;
                        sendParams[0] = sendParams[1] = sendParams[2] = sendParams[3] = sendParams[4] = 0;
                    }
                    else
                    {
                        sendParams[0] = interimImg1.extents[0];
                        sendParams[1] = interimImg1.extents[1];
                        sendParams[2] = interimImg1.boundingBox[2];
                        sendParams[3] = interimImg1.boundingBox[3];
                        sendParams[4] = 0;
                    }
                  prepareTimer.stop();

                  sendTimer.start();
                    MPI_Send(sendParams, 5, MPI_INT, destProc, localExTags[0], MPI_COMM_WORLD);

                    if (hasData)
                        MPI_Send(&interimImg1.data[sendingOffset], (sendParams[3]-sendParams[2])*(sendParams[1]-sendParams[0])*4, MPI_FLOAT, destProc, localExTags[1], MPI_COMM_WORLD);    // Send the actual image
                  sendTimer.stop();


                    addProfilingTimingMsg("p| from " + toString(myId) + ":" + toString(prepareTimer.getDuration()) + ",");
                    addTimingMsg( toString(myId) + " ~ Prepare " + toString(myId) +  " : " + toString(prepareTimer.getDuration()) + "\n" );

                    addTimingMsg( toString(myId) + " ~  Elapsed time sending to " + toString(destProc) + " : " + toString(sendParams[1]-sendParams[0]) + " x " + toString(sendParams[3]-sendParams[2]) + " : " + toString(sendTimer.getDuration()) + "\n" );
                    addProfilingTimingMsg("e|to " + toString(destProc) + " ~ "  + toString(sendParams[1]-sendParams[0]) + " x " + toString(sendParams[3]-sendParams[2]) + ":" + toString(sendTimer.getDuration()) + ",");

                    interimImg1.deleteImage();
                    compositingDone = true;
                }
                else
                {
                    //
                    // Receive
                    //

                    //
                    // Processes to receive from
                    std::vector<int> localExchange;
                    for (int i=0; i<k; i++)
                        if ( (myKIndex/k) * k  + i <= karyProcesses.size()-1)
                            localExchange.push_back( karyProcesses[ (myKIndex/k) * k + i] );

                    //std::cout << myId << " ~ myKIndex: " << myKIndex << " : " << localExchange[0] << ", " << localExchange[1] << "   round: " << rounds << "   karyProcesses.size(): " << karyProcesses.size() <<  "   localExchange.size(): " << localExchange.size() << std::endl;


                    //
                    // Create buffers for async reciving
                    MPI_Request *recvMetaRq = new MPI_Request[ localExchange.size()-1 ];
                    MPI_Request *recvImageRq = new MPI_Request[ localExchange.size()-1 ];

                    MPI_Status *recvMetaSt = new MPI_Status[ localExchange.size()-1 ];
                    MPI_Status *recvImageSt = new MPI_Status[ localExchange.size()-1 ];

                    
                    int recvCount=0;
                    for (int i=0; i<localExchange.size(); i++){
                        if ( localExchange[i] == myId )
                            continue;

                        int src = localExchange[i];

                        MPI_Irecv(&msgBuffer[i*5],                        5,   MPI_INT, src, localExTags[0], MPI_COMM_WORLD,  &recvMetaRq[recvCount] );
                        MPI_Irecv(&dataBuffer[i*myBufferSize], myBufferSize, MPI_FLOAT, src, localExTags[1], MPI_COMM_WORLD,  &recvImageRq[recvCount] );
                        recvCount++;
                    }


                    //
                    // Last round, blend with background
                    if (karyProcesses.size() <= k)
                    {
                        lastRound = true;
                        interimImg2.createOpaqueImage(0, width, myStartingHeight, myEndingHeight);

                        #ifdef __INTEL_COMPILER
                            alphaArray = (float *)_mm_malloc( width*(myEndingHeight-myStartingHeight) * sizeof(float), _ALIGN_);
                        #else
                            alphaArray = new float[ width*(myEndingHeight-myStartingHeight) ];
                        #endif

                        blendWithBackground(backgroundColor, interimImg2, interimImg1, alphaArray);
                        //finalImage.outputOpaquePPM("debug/final_with_bkg" + toString(myId) + "_.ppm");
                        interimImg1.deleteImage();
                    }

                  prepare.stop();
                    addProfilingTimingMsg("p| from " + toString(myId) + ":" + toString(prepare.getDuration()) + ",");
                    addTimingMsg( toString(myId) + " ~ Prepare " + toString(myId) +  " : " + toString(prepare.getDuration()) + "\n" );

                    //
                    // Receive
                    for (int receiveFrom = 1; receiveFrom < k; receiveFrom++)
                    {
                      recvTimer.start();
                        if (receiveFrom > localExchange.size()-1)
                            break;

                        if (myKIndex == localExchange.size()-1)   // the last one in the list; no neighbors
                            continue;

                        int sourceProc = karyProcesses[myKIndex + receiveFrom];


                        MPI_Wait(&recvMetaRq[receiveFrom-1], &recvMetaSt[receiveFrom-1]);
                        for (int j=0; j<4; j++)
                            recvImage.extents[j] = msgBuffer[receiveFrom*5 + j];

                        bool hasData =  false;
                        if (recvImage.extents[1]-recvImage.extents[0] > 0 && recvImage.extents[3]-recvImage.extents[2] > 0){
                            hasData = true;
                            MPI_Wait(&recvImageRq[receiveFrom-1], &recvImageSt[receiveFrom-1]);
                            recvImage.data = &dataBuffer[receiveFrom*myBufferSize];
                        }

                      recvTimer.stop();
                        addProfilingTimingMsg("f| from " + toString(sourceProc) + " ~ " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + ":" + toString(recvTimer.getDuration()) + ",");
                        addTimingMsg( toString(myId) + " ~ Receiving " + toString(sourceProc) + " size: " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + " : " + toString(recvTimer.getDuration()) + "\n" );
                

                        if (hasData == true)
                        {
                            //recvImage.outputPPM("output/kary_recv" + toString(myId) + "_" + toString(localExchange[receiveFrom]) );
                          blendTimer.start();
                            if (lastRound)
                            {
                                blendInto(FRONT_TO_BACK, interimImg2, recvImage, alphaArray);    // incorporate new image into old
                            }
                            else
                            {
                                blendInto(FRONT_TO_BACK, interimImg1, recvImage);    // incorporate new image into old
                                updateBoundingBox(interimImg1.boundingBox, recvImage.extents);
                            }
                          blendTimer.stop();
                            
                            addProfilingTimingMsg("b| from " + toString(sourceProc) + " ~ " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
                            addTimingMsg( toString(myId) + " ~ Blending " + toString(sourceProc) + " size: " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + " : " + toString(blendTimer.getDuration()) + "\n" );
                        }
                    }

                    delete []recvMetaRq;
                    recvMetaRq = NULL;
                    delete []recvImageRq;
                    recvImageRq = NULL;
                    delete []recvMetaSt;
                    recvMetaSt = NULL;
                    delete []recvImageSt;
                    recvImageSt = NULL;

                    #ifdef __INTEL_COMPILER
                        if ( lastRound )
                            _mm_free(alphaArray);
                        alphaArray = NULL;
                    #else
                        if ( lastRound )
                           delete []alphaArray;

                        alphaArray = NULL;
                    #endif
                }
            }

            //
            if (compositingDone || lastRound)
                break;
            
            //
            // Generate new karyProcesses list
          otherTimer.start();
            std::vector<int> _karyProcessesTemp;
            _karyProcessesTemp = karyProcesses;

            karyProcesses.clear();
            for (int i=0; i<_karyProcessesTemp.size(); i++)
                if (i%k == 0)
                    karyProcesses.push_back(_karyProcessesTemp[i]);

            _karyProcessesTemp.clear();
          otherTimer.stop();

            addTimingMsg( toString(myId) + " ~ Create new list: " + toString(otherTimer.getDuration()) + "\n" );
            addProfilingTimingMsg("o:" + toString(otherTimer.getDuration()) + ",");
        }


        #ifdef __INTEL_COMPILER
            _mm_free(dataBuffer);
            dataBuffer = NULL;
        #else
            delete []dataBuffer;
            dataBuffer = NULL;
        #endif


        // MPI_Barrier(MPI_COMM_WORLD);
        // std::cout << myId << " ~ ^^^^^^^^^ tree like compositing!!!" << width << ", " << height << std::endl;
        // if (compositingDone == false)
        //     interimImg2.outputOpaquePPM("output/karyDone_" + toString(myId));
        // MPI_Barrier(MPI_COMM_WORLD);


        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // Stage 3 - gather on root
        //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        gatherOpaqueImages(gatherProcesses.data(), gatherProcesses.size(), interimImg2, fullImg, tag[2], width, height);
      overallTimer.stop();

        // // Output
        // if (myId == 0)
        // {
        //     outputTimingMessage();
        //     clearTimingMsg();
        //     std::cout << "Total TOD compositing time: " << overallTimer.getDuration() << " s" << std::endl;
        //     compositingDone = false;
        // }

        //
        // Cleanup
        //

        // Async Send cleanup
        MPI_Waitall(r-1, sendMetaRq, sendMetaSt);
        MPI_Waitall(r-1, sendImageRq, sendImageSt);

        
        // Delete buffers
        delete []sendExtents;
        sendExtents = NULL;

        delete []sendMetaRq;
        sendMetaRq = NULL;
        delete []sendImageRq;
        sendImageRq = NULL;
        delete []sendMetaSt;
        sendMetaSt = NULL;
        delete []sendImageSt;
        sendImageSt = NULL;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Display the image
    // if (myId == 0)
    // {
    //     fullImg.outputOpaquePPM("output/tod_" +  toString(width) + "__" + toString(numProcesses) + "__" + toString(r)+ "_" + toString(k) + "_");
    // }

    // Point img to the new buffer
    //_img.deleteImage();

    // if (myId == 0){
    //     // Point img to the new buffer
    //     for (int i=0; i<4; i++)
    //         img.extents[i] = fullImg.extents[i];
    //     img.data = fullImg.data;
    // }
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
avtImgCommunicator::blendBackToFront(float *src, int dimsSrc[2], int posSrc[2],  float *dst, int dimsDst[2], int posDst[2])
{
    for (int _y=0; _y<dimsSrc[1]; _y++)
        for (int _x=0; _x<dimsSrc[0]; _x++)
        {
            int startingX = posSrc[0];
            int startingY = posSrc[1]; 

            if ((startingX + _x) > (posDst[0]+dimsDst[0]))
                continue;

            if ((startingY + _y) > (posDst[1]+dimsDst[1]))
                continue;
            
            int srcIndex = dimsSrc[0]*_y*4 + _x*4;                                     // index in the subimage 
            int dstIndex = ( (startingY+_y - posDst[1])*dimsDst[0]*4  + (startingX+_x - posDst[0])*4 );    // index in the big buffer

            // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
            float alpha = src[srcIndex+3];
            dst[dstIndex+0] = clamp( (dst[dstIndex+0] * (1.0 - alpha)) + src[srcIndex+0] );
            dst[dstIndex+1] = clamp( (dst[dstIndex+1] * (1.0 - alpha)) + src[srcIndex+1] );
            dst[dstIndex+2] = clamp( (dst[dstIndex+2] * (1.0 - alpha)) + src[srcIndex+2] );
            dst[dstIndex+3] = clamp( (dst[dstIndex+3] * (1.0 - alpha)) + src[srcIndex+3] );
        }
}





void 
avtImgCommunicator::gatherDepthAtRoot(int _numPatchGroup, int &totalPatches, int *patchCountPerRank, int *patchesDepth)
{
    //
    // Get how many patches are coming from each MPI rank
    totalPatches = 0;
    int *patchesOffset = NULL;

    if (my_id == 0) // root!
        patchCountPerRank = new int[num_procs];

    MPI_Gather(&_numPatchGroup, _numPatchGroup, MPI_INT,  &patchCountPerRank, num_procs, MPI_INT,    0, MPI_COMM_WORLD);



    //
    // Gather number of patch group
    if (my_id == 0)
    {
        patchesOffset = new int[num_procs];
        patchesOffset[0] = 0;

        for (int i=0; i<num_procs; i++)
        {
            totalPatches += patchCountPerRank[i];

            if (i == 0)
                patchesOffset[i] = 0;
            else
                patchesOffset[i] = patchesOffset[i-1] + patchCountPerRank[i-1]; 
        }

        patchesDepth = new float[totalPatches];
    }

    MPI_Gatherv(&_numPatchGroup, _numPatchGroup, MPI_FLOAT,  patchesDepth, patchCountPerRank, patchesOffset,  MPI_FLOAT, 0, MPI_COMM_WORLD);


    //
    // Cleanup
    if (i == 0)
        if (patchesOffset != NULL)
            delete []patchesOffset;
        
    patchesOffset = NULL;
} 



//
// Serial Direct Send
//
inline void avtImgCommunicator::serialDirectSend(int numPatches, float backgroundColor[4], int width, int height)
{
    float *recvImage = NULL; 
    float *fullImage = NULL;

    int tags[2] = {5781, 5782};

    int totalPatches;
    int *patchCountPerRank;
    int *patchesDepth;
    gatherDepthAtRoot(numPatches, totalPatches, patchCountPerRank, patchesDepth);

    //int myId = my_id;
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
                depthRankPatches.insert( std::pair<float,int>(patchCountPerRank[index],i) );
                index++;
            }
        

        //
        // Create space for buffers
        int recvParams[4];                          // minX, minY, maxX, maxY
        recvImage = new float[width*height*4]();
        fullImage = new float[width*height*4]();


        //
        // Compositing
        for (std::multimap<float,int>::iterator it=depthRankPatches.begin(); it!=depthRankPatches.end(); ++it)
        {
            int rank = (*it).second;

            if (rank != my_id)
            {
                MPI_Recv(recvParams,             4, MPI_INT,   rank, tags[0],  MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // recv image info
                MPI_Recv(recvImage, width*height*4, MPI_FLOAT, rank, tags[1],  MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // recv image

                dstPos[0]  = dstPos[0];                      dstPos[1]  = dstPos[1];
                dstSize[0] = recvParams[2]-recvParams[0];    dstSize[1] = recvParams[3]-recvParams[1];
            }
            else
            {
                // It's local

            }

            blendBackToFront(fullImage, srcSize, srcPos,  recvImage, dstSize, dstPos);
            //blendWithBackground();
        }
    }
    else
    {
        //
        // Sender
        for (int i=0; i<numPatches; i++)
        {
          //MPI_Send(sendParams, 5, MPI_INT, 0, tags[0], MPI_COMM_WORLD);             // send params

            MPI_Send( , 4, MPI_INT, 0, tags[0], MPI_COMM_WORLD);
            MPI_Send( ,  , MPI_INT, 0, tags[1], MPI_COMM_WORLD);
        }
    }
}

/*
inline void avtImgCommunicator::parallelDirectSend(Image _img, int region[], int numInRegion, int tags[3], float backgroundColor[4], int width, int height)
{
    Timer overallTimer;

    // Create final image for Display node
    // if (myId == 0)
    //     fullImg.createImage(0,width, 0,height);

    //for (int iter=0; iter<numIterations; iter++)
    //{
        //MPI_Barrier(MPI_COMM_WORLD);


      overallTimer.start();
        Timer blendTimer, prepareTimer, waitTimer, recvTimer, otherTimer;


      prepareTimer.start();

        // 
        // Determine position in region (positionInRegion)
        int positionInRegion = -1;
        bool inRegion = true;
        std::vector<int> regionVector (region, region+numInRegion);
        std::vector<int>::iterator it = std::find(regionVector.begin(), regionVector.end(), myId);

        if (it == regionVector.end()){
            inRegion = false;
            std::cout << myId << " ~ SHOULD NOT HAPPEN: Not found " << myId <<  " !!!" << std::endl;
        }
        else{
            positionInRegion = it - regionVector.begin();
        }


        //
        // Region boundaries
        int regionHeight = height/numInRegion;
        int lastRegionHeight = height - regionHeight*(numInRegion-1);

        int startingHeight = positionInRegion * regionHeight;   //  1*128 = 128
        int endingHeight = startingHeight + regionHeight;       // 128+128 = 256


        // Extents of my region
        if (positionInRegion == numInRegion-1) // the last one in region
            endingHeight = height;

        int myRegionHeight = endingHeight-startingHeight;


        // Size of one buffer
        int sizeOneBuffer = regionHeight*width * 4;


        //
        // MPI Async

        // Recv
        MPI_Request *recvMetaRq = new MPI_Request[ numInRegion-1 ];
        MPI_Request *recvImageRq = new MPI_Request[ numInRegion-1 ];

        MPI_Status *recvMetaSt = new MPI_Status[ numInRegion-1 ];
        MPI_Status *recvImageSt = new MPI_Status[ numInRegion-1 ];

        // Send
        MPI_Request *sendMetaRq = new MPI_Request[ numInRegion-1 ];
        MPI_Request *sendImageRq = new MPI_Request[ numInRegion-1 ];

        MPI_Status *sendMetaSt = new MPI_Status[ numInRegion-1 ];
        MPI_Status *sendImageSt = new MPI_Status[ numInRegion-1 ];


        //
        // Create Buffers

        // Create buffer for receiving images
        float *dataBuffer;
        #ifdef __INTEL_COMPILER
            dataBuffer = (float *)_mm_malloc( sizeOneBuffer * numInRegion * sizeof(float), _ALIGN_);
        #else
            dataBuffer = new float[ sizeOneBuffer * numInRegion];
        #endif

        // Create buffer for receiving messages
        std::vector<int> msgBuffer;
        msgBuffer.clear();
        msgBuffer.resize(5 * numInRegion); 

        // Create buffer for sending messages
        int *sendExtents = new int[numInRegion*5];


        //
        //  Position in region
        bool firstHalf = true;
        if (positionInRegion >= numInRegion/2)
            firstHalf = false;

        //
        // Async Recv
        if (inRegion)
        // MPI_Status *recvImageSt = new MPI_Status[ numInRegion-1 ];
        {
            int recvCount=0;
            for (int i=0; i<numInRegion; i++)
            {
                int index = i;
                if (!firstHalf)
                    index = (numInRegion-1)-i;

                if ( regionVector[index] == myId )
                    continue;

                int src = regionVector[index];
                MPI_Irecv(&msgBuffer[index*5],                          5, MPI_INT,   src, tags[0], MPI_COMM_WORLD,  &recvMetaRq[recvCount] );
                MPI_Irecv(&dataBuffer[index*sizeOneBuffer], sizeOneBuffer, MPI_FLOAT, src, tags[1], MPI_COMM_WORLD,  &recvImageRq[recvCount] );
                recvCount++;
            }
        }

        
        //
        // Async Send
        int sendCount = 0;
        int sendingOffset;
        for (int i=0; i<numInRegion; i++)
        {
            int regionStart, regionEnd, imgSize, index, dest;
            index = i;
            if (!firstHalf)
                index = (numInRegion-1)-i;
            dest = regionVector[index];

            if ( dest == myId )
                continue;

            regionStart = index*regionHeight;
            regionEnd = regionStart + regionHeight;
            if (index == numInRegion-1) // the last one in region
                regionEnd = height;


            if (regionStart < _img.extents[2])
                regionStart = _img.extents[2];

            if (regionEnd > _img.extents[3])
                regionEnd = _img.extents[3];


            bool hasData = true;
            if (regionEnd - regionStart <= 0 || _img.extents[1]-_img.extents[0] <= 0){
                hasData = false;

                sendingOffset = 0;
                imgSize = sendExtents[index*5 + 0] = sendExtents[index*5 + 1] = sendExtents[index*5 + 2] = sendExtents[index*5 + 3] =  sendExtents[index*5 + 4] = 0;
            }
            else
            {
                imgSize = (regionEnd-regionStart) * (_img.extents[1]-_img.extents[0]) * 4;
                sendingOffset = (regionStart-_img.extents[2]) * (_img.extents[1]-_img.extents[0]) * 4;

                sendExtents[index*5 + 0] = _img.extents[0];
                sendExtents[index*5 + 1] = _img.extents[1];
                sendExtents[index*5 + 2] = regionStart;
                sendExtents[index*5 + 3] = regionEnd;
                sendExtents[index*5 + 4] = 0;

                //std::cout << myId << " ~  sending to  " << locality[index] << " ... HAS DATA!!!!" << std::endl;
            }

            //std::cout << myId << " ~ index: " << index << "   regionVector[index]: " << regionVector[index] << "  extents: " <<  sendExtents[index*5 + 0] << ", " << sendExtents[index*5 + 1]  << ", " << sendExtents[index*5 + 2] << ", " << sendExtents[index*5 + 3] << "  sending ... " << std::endl;
            MPI_Isend(&sendExtents[index*5],           5,   MPI_INT, dest, tags[0], MPI_COMM_WORLD, &sendMetaRq[sendCount]);
            MPI_Isend(&_img.data[sendingOffset], imgSize, MPI_FLOAT, dest, tags[1], MPI_COMM_WORLD, &sendImageRq[sendCount]);
            
            sendCount++;
        }


        //
        // Create buffer for region
        interimImg.createImage(0, width, startingHeight, endingHeight);
        interimImg.initializeZero();


      prepareTimer.stop();
        addProfilingTimingMsg("p| from " + toString(myId) + ":" + toString(prepareTimer.getDuration()) + ",");
        addTimingMsg( toString(myId) + " ~ Prepare " + toString(myId) +  " : " + toString(prepareTimer.getDuration()) + "\n" );



        //
        // Blend
        int numBlends = 0;
        int countBlend = 0;
        int boundingBox[4];
        boundingBox[0] = boundingBox[2] = std::numeric_limits<int>::max();
        boundingBox[1] = boundingBox[3] = std::numeric_limits<int>::min();

        if (inRegion)
        {
            for (int i=0; i<numInRegion; i++)
            {
                int index;
                if (firstHalf)
                    index = i;
                else
                    index = (numInRegion-1)-i;

                if (regionVector[index] == myId)
                {
                  blendTimer.start();

                    int regionStart = startingHeight;
                    int regionEnd = endingHeight;

                    if (regionStart < _img.extents[2])
                        regionStart = _img.extents[2];

                    if (regionEnd > _img.extents[3])
                        regionEnd = _img.extents[3];


                    bool hasData = true;
                    if (regionEnd - regionStart <= 0){
                        hasData = false;
                        regionEnd = regionStart = 0;
                    }

                    if (hasData == true)
                    {
                        int extentsSectionRecv[4];
                        extentsSectionRecv[0] = _img.extents[0];
                        extentsSectionRecv[1] = _img.extents[1];
                        extentsSectionRecv[2] = regionStart;
                        extentsSectionRecv[3] = regionEnd;

                        if (firstHalf)
                            blendInto(FRONT_TO_BACK, interimImg, _img, extentsSectionRecv);
                        else
                            blendInto(BACK_TO_FRONT, interimImg, _img, extentsSectionRecv);
                        updateBoundingBox(boundingBox, extentsSectionRecv);
                        numBlends++;
                    }
                  blendTimer.stop();

                    addProfilingTimingMsg("b| from " + toString(myId) + " ~ " + toString(_img.extents[1]-_img.extents[0]) + "x"  + toString(_img.extents[3]-_img.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
                    addTimingMsg( toString(myId) + " ~ Blending " + toString(myId) + " size: " + toString(_img.extents[1]-_img.extents[0]) + "x"  + toString(_img.extents[3]-_img.extents[2]) + " : " + toString(blendTimer.getDuration()) + "\n" ); 
                }
                else
                {
                  recvTimer.start();
                    MPI_Wait(&recvMetaRq[countBlend], &recvMetaSt[countBlend]);

                    for (int j=0; j<4; j++)
                        recvImage.extents[j] = msgBuffer[index*5 + j];


                    bool hasData =  false;
                    if (recvImage.extents[1]-recvImage.extents[0] > 0 && recvImage.extents[3]-recvImage.extents[2] > 0){
                        hasData = true;
                        MPI_Wait(&recvImageRq[countBlend], &recvImageSt[countBlend]);
                        recvImage.data = &dataBuffer[index*sizeOneBuffer];
                    }
                  recvTimer.stop();
                    
                  blendTimer.start();
                    if (hasData)
                    {
                        if (firstHalf)
                            blendInto(FRONT_TO_BACK, interimImg, recvImage);
                        else
                            blendInto(BACK_TO_FRONT, interimImg, recvImage);
                        updateBoundingBox(boundingBox, recvImage.extents);
                        numBlends++;
                    }
                  blendTimer.stop();
                    countBlend++;

                    addProfilingTimingMsg("f| from " + toString(regionVector[index]) + " ~ " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + ":" + toString(recvTimer.getDuration()) + ",");
                    addTimingMsg( toString(myId) + " ~ Receiving " + toString(regionVector[index]) + " size: " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + " : " + toString(recvTimer.getDuration()) + "\n" );

                    addProfilingTimingMsg("b| from " + toString(regionVector[index]) + " ~ " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
                    addTimingMsg( toString(myId) + " ~ Blending " + toString(regionVector[index]) + " size: " + toString(recvImage.extents[1]-recvImage.extents[0]) + "x"  + toString(recvImage.extents[3]-recvImage.extents[2]) + " : " + toString(blendTimer.getDuration()) + "\n" );
                }
            }


            for (int i=0; i<4; i++)
                interimImg.boundingBox[i] = boundingBox[i];
        }
        else
            compositingDone = true;

        
        // Gather
        int gatherTags[2];
        gatherTags[0] = tags[2];     gatherTags[1] = tags[2]+1;
        gatherImages(region, numInRegion, interimImg, fullImg, gatherTags, width, height, backgroundColor);

      overallTimer.stop();

        // if (myId == 0)
        // {
        //     outputTimingMessage();
        //     clearTimingMsg();
        //     std::cout << "Total Parallel Direct Send compositing time: " << overallTimer.getDuration() << " s" << std::endl;
        //     compositingDone = false;
        // }


        msgBuffer.clear();
        #ifdef __INTEL_COMPILER
            _mm_free(dataBuffer);
            dataBuffer = NULL;
        #else
            delete []dataBuffer;
            dataBuffer = NULL;
        #endif

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

        interimImg.deleteImage();
    //}

    //_img.deleteImage();

  
    // Display the image
    // if (myId == 0)
    // {
    //     fullImg.outputPPM("output/pds_" +  toString(width) + "__" + toString(numProcesses) + "_");
    //     compositingDone = true;
    // }
}


inline void avtImgCommunicator::gatherImages(int regionGather[], int numToRecv, Image inputImg, Image outputImg, int tag[2], int width, int height, float backgroundColor[4])
{
    if (myId == 0)  // Display process
    {
        //
        // Receive at root/display node!
        Timer recvTimer, blendTimer, bkgTimer, setupTimer;

        int *msgBuffer = new int[numToRecv*5];
        int regionHeight = height/numToRecv;
        int lastRegionHeight = height - regionHeight*(numToRecv-1);
        if ( lastRegionHeight > regionHeight)
            regionHeight = lastRegionHeight;
        int bufferSize = width * regionHeight * 4;

        Image tempImg;
        float *gatherBuffer = NULL;

        #ifdef __INTEL_COMPILER
            gatherBuffer = (float *)_mm_malloc( bufferSize * numToRecv * sizeof(float), _ALIGN_);
        #else
            gatherBuffer = new float[ bufferSize * numToRecv];
        #endif


        // Adjust numtoRecv
        for (int i=0; i<numToRecv; i++)
            if (regionGather[i] == myId){
                numToRecv--;
                break;
            }
        
        //
        // Create buffers for async reciving
        MPI_Request *recvMetaRq = new MPI_Request[ numToRecv ];
        MPI_Request *recvImageRq = new MPI_Request[ numToRecv ];

        MPI_Status *recvMetaSt = new MPI_Status[ numToRecv ];
        MPI_Status *recvImageSt = new MPI_Status[ numToRecv ];


        // Async Recv
        int recvCount=0;
        for (int i=0; i<numToRecv; i++)
        {
            int src = regionGather[i];
            if (src == myId)
                continue;

            MPI_Irecv(&msgBuffer[i*5],                         5*sizeof(int),   MPI_INT, src, tag[0], MPI_COMM_WORLD,  &recvMetaRq[recvCount] );
            MPI_Irecv(&gatherBuffer[i*bufferSize],  bufferSize*sizeof(float), MPI_FLOAT, src, tag[1], MPI_COMM_WORLD,  &recvImageRq[recvCount] );
            recvCount++;
        }



        // create image with background
      bkgTimer.start();
        outputImg.colorImage(backgroundColor[0], backgroundColor[1], backgroundColor[2], backgroundColor[3]);
      bkgTimer.stop();
        addProfilingTimingMsg("b| with bkg ~ "  + toString(width) + " x " + toString(height) + ":" + toString(bkgTimer.getDuration()) + ",");


        // If root has data for the final image
        if (compositingDone == false)
        {
          blendTimer.start();
            blendInto(FRONT_TO_BACK, outputImg, inputImg);
          blendTimer.stop();

            //outputImg.outputPPM("output/blendLocal_" + toString(myId));
            addProfilingTimingMsg("z| local ~ "  + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
        }

       
        int recvBlend = 0;
        for (int i=0; i<numToRecv; i++)
        {
            MPI_Wait(&recvMetaRq[recvBlend], &recvMetaSt[recvBlend]);

            if ((msgBuffer[recvBlend*5+1]-msgBuffer[recvBlend*5+0]) * (msgBuffer[recvBlend*5+3]-msgBuffer[recvBlend*5+2]))
            {
                // Has Data!
              recvTimer.start();
                MPI_Wait(&recvImageRq[recvBlend], &recvImageSt[recvBlend]);
              recvTimer.stop();
                
                // Blend
              blendTimer.start();
                tempImg.extents[0] = msgBuffer[recvBlend*5 + 0];
                tempImg.extents[1] = msgBuffer[recvBlend*5 + 1];
                tempImg.extents[2] = msgBuffer[recvBlend*5 + 2];
                tempImg.extents[3] = msgBuffer[recvBlend*5 + 3];
                tempImg.data = &gatherBuffer[recvBlend*bufferSize];

                blendInto(FRONT_TO_BACK, outputImg, tempImg);

              blendTimer.stop();

                addProfilingTimingMsg("f| " + toString(i) + " :" + toString(recvTimer.getDuration()) + ",");
                addProfilingTimingMsg("z| " + toString(i) + " ~ " + toString(tempImg.extents[1]-tempImg.extents[0]) + " x " + toString(tempImg.extents[3]-tempImg.extents[2]) + ":" + toString(blendTimer.getDuration()) + ",");
            }

            recvBlend++;
        }

 
         #ifdef __INTEL_COMPILER
            _mm_free(gatherBuffer);
            gatherBuffer = NULL;
        #else
            delete []gatherBuffer;
            gatherBuffer = NULL;
        #endif

        if (recvMetaRq != NULL){
            delete []recvMetaRq;
            recvMetaRq = NULL;
        }

        if (recvImageRq != NULL){
           delete []recvImageRq;
           recvImageRq = NULL;
        }

        if (recvMetaSt != NULL){
            delete []recvMetaSt;
            recvMetaSt = NULL;
        }

        if (recvImageSt != NULL){
           delete []recvImageSt;
           recvImageSt = NULL;
        }

        delete []msgBuffer;
        msgBuffer = NULL;
    }
    else
    {
        if (compositingDone == false)   // If root has data for the final image
        {
            Timer sendImageTimer;

          sendImageTimer.start();
            int sendParams[5];
            for (int i=0; i<4; i++)
                sendParams[i] = inputImg.extents[i];
            sendParams[4]=0;

            MPI_Send(sendParams,                                                      5, MPI_INT,   0, tag[0], MPI_COMM_WORLD);    // Send meta data
            MPI_Send(inputImg.data, inputImg.getNumPixels()*inputImg.getNumComponents(), MPI_FLOAT, 0, tag[1], MPI_COMM_WORLD);    // Send the actual image
          sendImageTimer.stop();

            addTimingMsg(toString(myId) +" ~ Gather Send ~ " + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(sendImageTimer.getDuration()) + "\n");
            addProfilingTimingMsg("v| " + toString(myId) + " ~ " + toString(inputImg.extents[1]-inputImg.extents[0]) + " x " + toString(inputImg.extents[3]-inputImg.extents[2]) + ":" + toString(sendImageTimer.getDuration()) + ",");
        }
    }
}
*/

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
#ifdef PARALLEL
MPI_Datatype avtImgCommunicator::createMetaDataType(){
  MPI_Datatype _imgMeta_mpi;
  const int numItems = 8;
  int blockLengths[numItems] = {1, 1, 1, 2, 2, 2, 1, 1};  // 11
  MPI_Datatype type[numItems] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT};
  MPI_Aint offsets[numItems] = {0, sizeof(int), sizeof(int)*2, sizeof(int)*3, sizeof(int)*5, sizeof(int)*7, sizeof(int)*9, sizeof(int)*10 };
  MPI_Type_struct(numItems, blockLengths,  offsets, type, &_imgMeta_mpi);
  
  return _imgMeta_mpi;
}
#endif

void avtImgCommunicator::barrier(){
  #ifdef PARALLEL
    MPI_Barrier( MPI_COMM_WORLD );
  #endif
}
