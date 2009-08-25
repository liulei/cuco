#ifdef __cplusplus
extern "C"{
#endif

#include	"radixsort.h"

typedef struct tagSimParam{
	float	shortrange_table[NTAB];
	float	boxsize;
	float	boxhalf;
	float	to_grid_fac;
	float	rcut;
	float	rcut2;
	float	asmth;
	float	asmthfac;
	float	mass;
	float	G;
}SIMPARAM;

typedef struct tagSoftParam{
	float	h;
	float	h_inv;
	float	h3_inv;
}SOFTPARAM;

void	cudaInit(void);
void	allocateArray(void **devPtr, size_t size);
void	freeArray(void *devPtr);
void	copyArrayFromDevice(void *host, void *device, int size);
void	copyArrayToDevice(void *device, void *host, int offset, int size);
void	setSimParam(SIMPARAM *pHostSimParam);
void	setSoftParam(SOFTPARAM *pHostSoftParam);
void	copyPosToDevice(void);
void	copyAccelFromDevice(void);

void	calcHash(uint	*gridParticleHash,
				 uint	*gridParticleIndex,
				 float	*pos,
				 uint	numParticles);

void	reorderDataAndFindCellStart(
				uint	*cellStart,
				uint	*cellEnd,
				float	*sortedPos,
				uint	*gridParticleHash,
				uint	*gridParticleIndex,
				float	*oldPos,
				uint	numParticles,
				uint	numCells);

void	cudaForceEvaluateShortrange(
				float	*gravAccel,
				float	*sortedPos,
				uint	*gridParticleIndex,
				uint	*cellStart,
				uint	*cellEnd,
				uint	numParticles,
				uint	numCells);

extern	RadixSort	*sorter;

#ifdef __cplusplus
}
#endif
