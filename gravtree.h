#ifndef	_GRAVTREE_KERNEL_H_
#define	_GRAVTREE_KERNEL_H_

__constant__ SIMPARAM	dSimParam;

__constant__ SOFTPARAM	dSoftParam;

__global__ void force_treeevaluate_shortrange_device(int numParticles){

	uint	index	=	blockIdx.y * gridDim.x * blockDim.x
						+ blockIdx.x * blockDim.x + threadIdx.x;
	if(index >= numParticles){
		return;
	}

}

#endif
