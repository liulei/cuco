
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>

#include	<cutil_inline.h>

#include	<cutil_math.h>

#include	"allvars.h"
#include	"proto.h"

#include	"radixsort.h"
#include	"proto.cuh"
#include	"gravlist_kernel.cu"

extern "C"{

RadixSort	*sorter;

static uint	maxThreads	=	16384;
static uint	numThreads	=	256;

void cudaInit(){
	
	printf("Init cuda device...\n");

	cudaSetDevice( cutGetMaxGflopsDeviceId());

}

void allocateArray(void **devPtr, size_t size){

	cutilSafeCall(cudaMalloc(devPtr, size));

}

void freeArray(void *devPtr){

	cutilSafeCall(cudaFree(devPtr));

}

void copyArrayFromDevice(void *host, void *device, int size){

	cutilSafeCall(cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost));

}

void copyArrayToDevice(void *device, void *host, int offset, int size){
	cutilSafeCall(cudaMemcpy((char *) device + offset, host, size, cudaMemcpyHostToDevice));

}

void copyPosToDevice(){
	
	float4	*hPos4	=	(float4 *)hPos;

	uint	i;
	for(i = 0; i < NumPart; ++i){

		hPos4[i].x	=	P[i].Pos[0];
		hPos4[i].y	=	P[i].Pos[1];
		hPos4[i].z	=	P[i].Pos[2];
	}
	copyArrayToDevice(dOldPos, hPos, 0, numParticles * sizeof(float4));
}

void copyAccelFromDevice(){
	
	copyArrayFromDevice(hGravAccel, dGravAccel, numParticles * sizeof(float4));

	float4	*hGravAccel4	=	(float4 *)hGravAccel;

	int	i;
	for(i = 0; i < NumPart; ++i){

		P[i].GravAccel[0]	=	hGravAccel4[i].x;
		P[i].GravAccel[1]	=	hGravAccel4[i].y;
		P[i].GravAccel[2]	=	hGravAccel4[i].z;
	}
}

void calcHash(uint * gridParticleHash,
			  uint * gridParticleIndex,
			  float *pos,
			  uint	numParticles){

	dim3	dimBlock(numThreads, 1);
	dim3	dimGrid;
	if(numParticles > maxThreads){
		dimGrid.y	=	numParticles / maxThreads;
		dimGrid.x	=	maxThreads / numThreads;
	}else{
		dimGrid.y	=	1;
		dimGrid.x	=	numParticles / numThreads;
	}

	calcHashD<<<dimGrid, dimBlock>>>(gridParticleHash,
										 gridParticleIndex,
										 (float4 *) pos,
										 numParticles);
}



void reorderDataAndFindCellStart(uint 	*cellStart,
								 uint 	*cellEnd,
								 float	*sortedPos,
								 uint	*gridParticleHash,
								 uint	*gridParticleIndex,
								 float	*oldPos,
								 uint	numParticles,
								 uint	numCells){
	
	dim3	dimBlock(numThreads, 1);
	dim3	dimGrid;
	if(numParticles > maxThreads){
		dimGrid.y	=	numParticles / maxThreads;
		dimGrid.x	=	maxThreads / numThreads;
	}else{
		dimGrid.y	=	1;
		dimGrid.x	=	numParticles / numThreads;
	}
	cutilSafeCall(cudaMemset(cellStart, 0xffffffff, numCells * sizeof(uint)));

	cutilSafeCall(cudaBindTexture(0, oldPosTex, oldPos, numParticles * sizeof(float4)));

	uint	smemSize	=	sizeof(uint) * (numThreads + 1);

	reorderDataAndFindCellStartD<<< dimGrid, dimBlock, smemSize>>>(
		cellStart,
		cellEnd,
		(float4 *) sortedPos,
		gridParticleHash,
		gridParticleIndex,
		(float4 *) oldPos,
		numParticles);

	cutilSafeCall(cudaUnbindTexture(oldPosTex));
}

void cudaForceEvaluateShortrange(float	*gravAccel,
								 float	*sortedPos,
								 uint	*gridParticleIndex,
								 uint	*cellStart,
								 uint	*cellEnd,
								 uint	numParticles,
								 uint	numCells){
	cutilSafeCall(cudaBindTexture(0, oldPosTex, sortedPos, numParticles * sizeof(float4)));
	cutilSafeCall(cudaBindTexture(0, cellStartTex, cellStart, numCells * sizeof(uint)));
	cutilSafeCall(cudaBindTexture(0, cellEndTex, cellEnd, numCells * sizeof(uint)));

	dim3	dimBlock(numThreads, 1);
	dim3	dimGrid;
	if(numParticles > maxThreads){
		dimGrid.y	=	numParticles / maxThreads;
		dimGrid.x	=	maxThreads / numThreads;
	}else{
		dimGrid.y	=	1;
		dimGrid.x	=	numParticles / numThreads;
	}

	cudaForceEvaluateShortrangeD<<<dimGrid, dimBlock>>>((float4 *)gravAccel,
															(float4 *)sortedPos,
															gridParticleIndex,
															cellStart,
															cellEnd,
															numParticles);
}

void setSimParam(SIMPARAM *pHostSimParam){

	cutilSafeCall(cudaMemcpyToSymbol(dSimParam, pHostSimParam, sizeof(SIMPARAM)));

}

void setSoftParam(SOFTPARAM *pHostSoftParam){
	cutilSafeCall(cudaMemcpyToSymbol(dSoftParam, pHostSoftParam, sizeof(SOFTPARAM)));
}

}
