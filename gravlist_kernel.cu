#ifndef	_GRAVLIST_KERNEL_H_
#define	_GRAVLIST_KERNEL_H_

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#define	FETCH(t, i) tex1Dfetch(t##Tex, i)

texture<float4, 1, cudaReadModeElementType> oldPosTex;

texture<uint, 1, cudaReadModeElementType> gridParticleHashTex;
texture<uint, 1, cudaReadModeElementType> cellStartTex;
texture<uint, 1, cudaReadModeElementType> cellEndTex;

__constant__ SIMPARAM	dSimParam;

__constant__ SOFTPARAM	dSoftParam;

__device__ float3 calcForce(float3 pos, float3 pos2){
	
	float3	acc	=	make_float3(0.0f);

	float	dx	=	pos2.x - pos.x;
	float	dy	=	pos2.y - pos.y;
	float	dz	=	pos2.z - pos.z;
	
	float	boxsize	=	dSimParam.boxsize;
	float	boxhalf	=	dSimParam.boxhalf;
	
	if(dx > boxhalf)
		dx	-=	boxsize;
	if(dx < -boxhalf)
		dx	+=	boxsize;

	if(dy > boxhalf)
		dy	-=	boxsize;
	if(dy < -boxhalf)
		dy	+=	boxsize;

	if(dz > boxhalf)
		dz	-=	boxsize;
	if(dz < -boxhalf)
		dz	+=	boxsize;

	float	r2	=	dx * dx + dy * dy + dz * dz;

	float	rcut	=	dSimParam.rcut;
	float	rcut2	=	dSimParam.rcut2;
	float	asmth	=	dSimParam.asmth;
	float	asmthfac	=	dSimParam.asmthfac;
	float	mass	=	dSimParam.mass;

	if(r2 > rcut2)
		return acc;
	
	float	h	=	dSoftParam.h;
	float	h_inv	=	dSoftParam.h_inv;
	float	h3_inv	=	dSoftParam.h3_inv;

	float	r	=	sqrtf(r2);
	float	fac, u;

	if(r >= h){
		fac	=	mass / (r2 * r);
	}else{
		u	=	r * h_inv;
		if(u < 0.5)
			fac	=	mass * h3_inv * (10.66667 + u * u * (32.0 * u - 38.4));
		else
			fac	=	mass * h3_inv * (21.33333 - 48.0 * u 
			+ 38.4 * u * u - 10.66667 * u * u * u - 0.06667 / (u * u * u));
	}

	int	tabindex	=	(int) (asmthfac * r);

	if(tabindex < NTAB){

		fac	*=	dSimParam.shortrange_table[tabindex];

		acc.x	+=	dx * fac;
		acc.y	+=	dy * fac;
		acc.z	+=	dz * fac;
	}

	acc	*=	dSimParam.G;
	return acc;
}



__global__
void calcHashD(uint *gridParticleHash,
			   uint *gridParticleIndex,
			   float4	*pos,
			   uint	numParticles){

	uint	index	=	blockIdx.y * gridDim.x * blockDim.x 
						+ blockIdx.x * blockDim.x + threadIdx.x;
	if(index >= numParticles)
		return;
	
	float	gridFac	=	dSimParam.to_grid_fac;

	float4	p	=	pos[index];
	
	int3	gridPos;

	gridPos.x	=	gridFac * p.x;
	gridPos.y	=	gridFac * p.y;
	gridPos.z	=	gridFac * p.z;

	uint	hash;
	hash	=	(gridPos.x * PMGRID + gridPos.y) * PMGRID + gridPos.z;

	gridParticleHash[index]	=	hash;
	gridParticleIndex[index]	=	index;
}


__global__
void reorderDataAndFindCellStartD(uint	*cellStart,
								  uint	*cellEnd,
								  float4	*sortedPos,
								  uint	*gridParticleHash,
								  uint	*gridParticleIndex,
								  float4	*oldPos,
								  uint	numParticles){
	extern __shared__ uint	sharedHash[];
	
	uint	index	=	blockIdx.y * gridDim.x * blockDim.x 
						+ blockIdx.x * blockDim.x + threadIdx.x;
	if(index >= numParticles)
		return;
	
	uint hash;
	hash	=	gridParticleHash[index];

	sharedHash[threadIdx.x + 1]	=	hash;

	if(index > 0 && threadIdx.x == 0){
		sharedHash[0]	=	gridParticleHash[index - 1];
	}

	__syncthreads();

	if(index == 0 || hash != sharedHash[threadIdx.x]){
		cellStart[hash]	=	index;
		if(index > 0)
			cellEnd[sharedHash[threadIdx.x]] =	index;
	}

	if(index == numParticles - 1){
		cellEnd[hash]	=	index + 1;
	}

	uint sortedIndex	=	gridParticleIndex[index];
	float4	pos	=	FETCH(oldPos, sortedIndex);
	sortedPos[index]	=	pos;
}

__global__
void cudaForceEvaluateShortrangeD(float4 *gravAccel,
								  float4 *oldPos,
								  uint	*gridParticleIndex,
								  uint	*cellStart,
								  uint	*cellEnd,
								  uint	numParticles){
	uint	index	=	blockIdx.y * gridDim.x * blockDim.x 
						+ blockIdx.x * blockDim.x + threadIdx.x;
	if(index >= numParticles)
		return;

	float3	pos	=	make_float3(FETCH(oldPos, index));
	float3	acc	=	make_float3(0.0f);

	float	gridFac	=	dSimParam.to_grid_fac;
	float	rcut	=	dSimParam.rcut;

	int	xl, xr, yl, yr, zl, zr;

	xl	=	gridFac * (pos.x - rcut) + PMGRID;
	xl	-=	PMGRID;
	xr	=	gridFac * (pos.x + rcut);

	yl	=	gridFac * (pos.y - rcut) + PMGRID;
	yl	-=	PMGRID;
	yr	=	gridFac * (pos.y + rcut);

	zl	=	gridFac * (pos.z - rcut) + PMGRID;
	zl	-=	PMGRID;
	zr	=	gridFac * (pos.z + rcut);

	int	ix, iy, iz, iix, iiy, iiz;
	for(ix = xl; ix <= xr; ++ix){
		for(iy = yl; iy <= yr; ++iy){
			for(iz = zl; iz <= zr; ++iz){
				
				iix	=	ix;
				iiy	=	iy;
				iiz	=	iz;

				if(iix < 0)
					iix	+=	PMGRID;
				if(iiy < 0)
					iiy	+=	PMGRID;
				if(iiz < 0)
					iiz	+=	PMGRID;

				if(iix >= PMGRID)
					iix	-=	PMGRID;
				if(iiy >= PMGRID)
					iiy	-=	PMGRID;
				if(iiz >= PMGRID)
					iiz	-=	PMGRID;

				uint gridHash	=	(iix * PMGRID + iiy) * PMGRID + iiz;

				uint startIndex	=	FETCH(cellStart, gridHash);
				
				if(startIndex != 0xffffffff){
					uint endIndex	=	FETCH(cellEnd, gridHash);
					for(uint j = startIndex; j < endIndex; ++j){
						float3	pos2	=	make_float3(FETCH(oldPos, j));
						acc	+=	calcForce(pos, pos2);
					}
				}
			}
		}
	}

	uint originalIndex	=	gridParticleIndex[index];
	gravAccel[originalIndex]	=	make_float4(acc, 0.0f);
}

#endif
