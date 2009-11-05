#ifndef	_GRAVTREE_KERNEL_H_
#define	_GRAVTREE_KERNEL_H_

#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))

#define	BUSY	(-250)

texture<float4, 1, cudaReadModeElementType> dPosTex;
texture<float4, 1, cudaReadModeElementType> dNodesTex;
texture<int, 1, cudaReadModeElementType> dNextnodeTex;

__constant__ SIMPARAM	dSimParam;

__constant__ SOFTPARAM	dSoftParam;

__device__ int	nfree, dNumNodes;

__global__ void force_treebuild_device(
					float4	*dPos,
					NODE	*dNodes,
					SUNS	*dSuns,
					int		numParticles){
	
	int	index	=	blockIdx.y * gridDim.x * blockDim.x
						+ blockIdx.x * blockDim.x + threadIdx.x;
	if(index >= numParticles){
		return;
	}

	__shared__ int	count;

	if(threadIdx.x == 0){
		count	=	0;
		__threadfence_block();
	}

	NODE	node, pnode;
	SUNS	suns;
	int		subnode, nn, thisNfree, thSubnode;
	float	lenhalf;
//	float	mass	=	dSimParam.mass;

	float4	pos		=	dPos[index];
	float4	thPos;
	
	int	th	=	numParticles;
	int	parent	=	-1;

//	int	run_flag	=	1;

	while(1){
//	while(count < blockDim.x){
//	while(threadIdx.x == count){
//		printf("%d: new loop, entry th: %d...\n", index, th);

//		__syncthreads();
		__threadfence_block();
		
//		if(threadIdx.x == count)
//		if(run_flag	== 1)
		if(th >= numParticles){
			
			node	=	dNodes[th];

			subnode	=	0;
			if(pos.x > node.center[0])
				subnode	+=	1;
			if(pos.y > node.center[1])
				subnode	+=	2;
			if(pos.z > node.center[2])
				subnode	+=	4;
			
			nn	=	atomicExch(&dSuns[th].suns[subnode], BUSY);
//			__threadfence();
			if(nn != BUSY){
				if(nn >= 0){
//					if(nn >= numParticles){
//						dNodes[nn].u.d.mass	+=	mass;

//						dNodes[nn].u.d.s[0]	+=	mass * dPos[nn].x;
//						dNodes[nn].u.d.s[1]	+=	mass * dPos[nn].y;
//						dNodes[nn].u.d.s[2]	+=	mass * dPos[nn].z;
//					}
					atomicExch(&dSuns[th].suns[subnode], nn);
//					__threadfence();
//					printf("%d: subnode %d of node %d: reset to %d, continue...\n", index, subnode, th, nn);
					parent	=	th;
					th	=	nn;
				}
				if(nn < 0){
					atomicExch(&dSuns[th].suns[subnode], index);
//					__threadfence();
//					printf("%d: subnode %d of node %d: set to %d, inserting done!\n", index, subnode, th, index);
					count++;
//					run_flag	=	0;
//					__threadfence_block();
//					__syncthreads();
					break;
				}
//			}else{
//				printf("%d: waiting for subnode %d of node %d...\n", index, subnode, th);
//				continue;
			}
		}
			
//		if(threadIdx.x == count)
//		if(run_flag	== 1)
		if(th < numParticles){
			
			nn	=	atomicExch(&dSuns[parent].suns[subnode], BUSY);
//			__threadfence();
			if(nn != BUSY){
				
				thisNfree	=	atomicAdd(&nfree, 1);
//				__threadfence();
				thisNfree++;

//				printf("%d: subnode %d of node %d: generating new node %d for particle %d...\n", index, subnode, parent, thisNfree, nn);

				pnode		=	dNodes[parent];
				node.len	=	0.5 * pnode.len;
				lenhalf		=	0.25 * pnode.len;
				
				if(subnode & 1)
					node.center[0]	=	pnode.center[0] + lenhalf;
				if(!(subnode & 1))
					node.center[0]	=	pnode.center[0] - lenhalf;

				if(subnode & 2)
					node.center[1]	=	pnode.center[1] + lenhalf;
				if(!(subnode & 2))
					node.center[1]	=	pnode.center[1] - lenhalf;

				if(subnode & 4)
					node.center[2]	=	pnode.center[2] + lenhalf;
				if(!(subnode & 4))
					node.center[2]	=	pnode.center[2] - lenhalf;

				node.u.d.s[0]	=	0.0;
				node.u.d.s[1]	=	0.0;
				node.u.d.s[2]	=	0.0;

				node.u.d.mass	=	0.0;

				suns.suns[0]	=	-1;
				suns.suns[1]	=	-1;
				suns.suns[2]	=	-1;
				suns.suns[3]	=	-1;
				suns.suns[4]	=	-1;
				suns.suns[5]	=	-1;
				suns.suns[6]	=	-1;
				suns.suns[7]	=	-1;

				thSubnode	=	0;
				thPos	=	dPos[th];
				if(thPos.x > node.center[0])
					thSubnode	+=	1;
				if(thPos.y > node.center[1])
					thSubnode	+=	2;
				if(thPos.z > node.center[2])
					thSubnode	+=	4;

				suns.suns[thSubnode]	=	th;
/*
				node.u.d.s[0]	=	mass * thPos.x;
				node.u.d.s[1]	=	mass * thPos.y;
				node.u.d.s[2]	=	mass * thPos.z;

				node.u.d.mass	=	mass;
*/

				node.u.d.bitflags	=	0;

				dNodes[thisNfree]	=	node;
				dSuns[thisNfree]	=	suns;

				th	=	thisNfree;

				atomicAdd(&dNumNodes, 1);
//				__threadfence();
				atomicExch(&dSuns[parent].suns[subnode], thisNfree);
//				__threadfence();
//				printf("%d: subnode %d of node %d set to new node %d !\n", index, subnode, parent, thisNfree);
//			}else{
//				printf("%d: waiting for particle in subnode %d of node %d...\n", index, subnode, parent);
//				continue;
			}
		}
	}

/*	Don't forget update dNodes[numParticles]! Also every node need be 
 *	devided by mass!
 */
//	count++;
//	printf("%d: return, node number: %d\n", index, dNumNodes);
}

__global__ void read_tree_information(int *pNumNodes){
	*pNumNodes	=	dNumNodes;
}

__global__ void update_tree_device(
							float4	*dPos, 
							NODE	*dNodes, 
							SUNS	*dSuns, 
							int		numParticles){
	__shared__ int	flag[256];

	flag[threadIdx.x]	=	-1;
	__syncthreads();

	int	numBlocks	=	gridDim.x * gridDim.y;
	int	blockIndex	=	blockIdx.y * gridDim.x + blockIdx.x;
	int	index	=	(numBlocks - 1 - blockIndex) * blockDim.x
						+ threadIdx.x;

//	int	index	=	blockIdx.y * gridDim.x * blockDim.x
//						+ blockIdx.x * blockDim.x + threadIdx.x;
	
	if(index >= dNumNodes){
//		return;
	}

	if(blockIdx.x > 0){
//		printf("index: %d: (%d, %d)\n", index, blockIdx.x, blockIdx.y);
	}

	if(index < dNumNodes)
		flag[threadIdx.x]	=	0;
	__syncthreads();

	index	+=	numParticles;

	int	i, ss;
	int	ready	=	0;
	int	activeThreads	=	0;
	int	count	=	0;
	int	done	=	0;
	float	particleMass	=	dSimParam.mass;
	float	mass	=	0.0;
	float	s[3];
	float4	pos;
	NODE	node;


	s[0] = s[1] = s[2] = 0.0;
	
	for(i = 0; i < blockDim.x; ++i){
		if(flag[i] == 0){
			activeThreads++;
		}
	}

//	printf("%d: count: %d\tactive threads: %d\n", index, count, activeThreads);

	__syncthreads();
	SUNS suns	=	dSuns[index];
	while(count	< activeThreads){
//	while(blockIdx.x == 1 && count < activeThreads){

		__syncthreads();
		__threadfence_block();

//		printf("%d: done: %d\tready: %d\n", index, done, ready);

		if(done == 0){
			ready	=	1;
			for(i = 0; i < 8; ++i){
				if(suns.suns[i] >= numParticles){
					ss	=	suns.suns[i];
					if(dNodes[ss].u.d.bitflags == 0){
						ready	=	0;
					}
				}
			}


			if(ready == 1){
				for(i = 0; i < 8; ++i){
					ss	=	suns.suns[i];
					if(ss >= 0){
						if(ss < numParticles){
							pos	=	dPos[ss];
							mass	+=	particleMass;
							s[0]	+=	particleMass * pos.x;
							s[1]	+=	particleMass * pos.y;
							s[2]	+=	particleMass * pos.z;
						}
						if(ss >= numParticles){
							node	=	dNodes[ss];
							mass	+=	node.u.d.mass;
							s[0]	+=	node.u.d.s[0] * node.u.d.mass;
							s[1]	+=	node.u.d.s[1] * node.u.d.mass;
							s[2]	+=	node.u.d.s[2] * node.u.d.mass;
						}
					}
				}

				s[0]	/=	mass;
				s[1]	/=	mass;
				s[2]	/=	mass;

				dNodes[index].u.d.mass	=	mass;
				
				dNodes[index].u.d.s[0]	=	s[0];
				dNodes[index].u.d.s[1]	=	s[1];
				dNodes[index].u.d.s[2]	=	s[2];

				dNodes[index].u.d.bitflags	=	1;

				done	=	1;
					
				flag[threadIdx.x]	=	1;
			}
		}
		
		__syncthreads();
		count	=	0;
		for(i = 0; i < activeThreads; ++i){
			if(index == numParticles){
				printf("%d", flag[i]);
			}
			if(flag[i] == 1){
				count++;
			}
		}
		if(index == numParticles){
			printf("\n");
		}
//		printf("%d: count(after while): %d\n", index, count);
	}
}

__global__ void update_treenext_device(
							float4	*dPos, 
							NODE	*dNodes, 
							SUNS	*dSuns, 
							int		*dNextnode,
							int		numParticles){
	__shared__ int	flag[256];

	flag[threadIdx.x]	=	-1;
	__syncthreads();

	int	index	=	blockIdx.y * gridDim.x * blockDim.x
						+ blockIdx.x * blockDim.x 
						+ blockDim.x - 1 - threadIdx.x;
	
	if(index < dNumNodes){
		flag[blockDim.x - 1 - threadIdx.x]	=	0;
	}
	__syncthreads();

	index	+=	numParticles;

	int	i, j, jj, p, pp, nextsib, father, ready;
	int	activeThreads	=	0;
	int	count	=	0;
	int	done	=	0;
	SUNS	suns;

	if(index == numParticles){
		dNodes[index].u.d.bitflags	=	0;
		dNodes[index].u.d.sibling	=	-1;
		dNodes[index].u.d.father	=	numParticles;
		suns	=	dSuns[index];
		for(j = 0; j < 8; ++j){
			if(suns.suns[j] >= 0){
				break;
			}
		}
		dNodes[index].u.d.nextnode	=	suns.suns[j];
	}
	__syncthreads();

	if(index > numParticles && index < numParticles + dNumNodes){
		dNodes[index].u.d.bitflags	=	0;
	}

	__syncthreads();

	for(i = 0; i < blockDim.x; ++i){
		if(flag[i] == 0){
			activeThreads++;
		}
	}

//	printf("%d: count: %d\tactive threads: %d\n", index, count, activeThreads);

	__syncthreads();
	if(index < numParticles + dNumNodes){
		suns	=	dSuns[index];
		for(i = 0; i < 8; ++i){
			p	=	suns.suns[i];
			if(p >= numParticles){
				dNodes[p].u.d.father	=	index;
			}
		}
	}
	__syncthreads();
	while(count	< activeThreads){
//	while(blockIdx.x == 1 && count < activeThreads){
		__syncthreads();
		__threadfence();
		
		if((index < numParticles + dNumNodes) && (done == 0)){
			father	=	dNodes[index].u.d.father;
			ready	=	0;
			if(dNodes[father].u.d.bitflags == 1){
				ready	=	1;
			}
			if(index == numParticles){
				ready	=	1;
			}

			if(ready == 1){

//				printf("%d: working node: bitflag: %d\n", index, dNodes[index].u.d.bitflags);

				for(j = 0; j < 8; ++j){
					if(suns.suns[j] >= 0){
						break;
					}
				}
				dNodes[index].u.d.nextnode	=	suns.suns[j];

				for(j = 0; j < 8; ++j){
					p	=	suns.suns[j];
					if(p >= 0){
						for(jj = j + 1; jj < 8; ++jj){
							pp	=	suns.suns[jj];
							if(pp >= 0)
								break;
						}

						if(jj < 8)
							nextsib	=	pp;
						else
							nextsib	=	dNodes[index].u.d.sibling;

						if(p < numParticles){
							dNextnode[p]	=	nextsib;
						}

						if(p >= numParticles){
							dNodes[p].u.d.sibling	=	nextsib;
//							printf("%d: set bitflag of node %d\n", index, p);
//							dNodes[p].u.d.bitflags	=	p;
						}
					}
				}
				dNodes[index].u.d.bitflags	=	1;
				done	=	1;
				flag[blockDim.x - 1- threadIdx.x]	=	1;
			}
		}
		
		__syncthreads();
		count	=	0;
		for(i = 0; i < activeThreads; ++i){
			if(index == numParticles){
				printf("%d", flag[i]);
			}
			if(flag[i] == 1){
				count++;
			}
		}
		if(index == numParticles){
			printf("\n");
		}
	}
}

__global__ void set_numnodes(SUNS *dSuns, int numParticles){
	dSuns[2 * numParticles - 1].suns[0]	=	dNumNodes;
	dSuns[2 * numParticles - 1].suns[1]	=	nfree;
}


__global__ void set_variable_device(NODE *dNodes, SUNS *dSuns, int numParticles){
	
	NODE	node;
	SUNS	suns;
	int	i;

	float	boxsize	=	dSimParam.boxsize;
	float	boxhalf	=	boxsize / 2.0;

	node.len	=	boxsize;
	for(i = 0; i < 3; ++i){
		node.center[i]	=	boxhalf;
//		node.u.d.s[i]	=	0.0;
	}
//	node.u.d.mass	=	0.0;

	for(i = 0; i < 8; ++i)
		suns.suns[i]	=	-1;
	
	nfree	=	numParticles;
	dNodes[nfree]	=	node;
	dSuns[nfree]	=	suns;

	for(i = 0; i < 8; ++i){
//		printf("subnode %d of node %d: %d\n", i, nfree, dSuns[nfree].suns[i]);
	}

	dNumNodes	=	1;
	__threadfence();
}

__global__ void force_treeevaluate_shortrange_device(
						float4	*dPos,
						float4	*dGravAccel,
						NODE	*dNodes,
						int		*dNextnode,
						int 	numParticles){

	uint	index	=	blockIdx.y * gridDim.x * blockDim.x
						+ blockIdx.x * blockDim.x + threadIdx.x;
	if(index >= numParticles){
		return;
	}

//	NODE	node;
	NODE_1	node_1;
	NODE_2	node_2;
	NODE_3	node_3;
	int		no, tabindex;
	float	r2, dx, dy, dz, r, fac, u;
	float	eff_dist, dist;
	float4	pos, pos_no, acc;

	float4	temp;
	
	float	rcut	=	dSimParam.rcut;
	float	rcut2	=	dSimParam.rcut2;
	float	asmthfac	=	dSimParam.asmthfac;
	float	mass	=	dSimParam.mass;
	float	ErrTolTheta	=	dSimParam.ErrTolTheta;
	float	G		=	dSimParam.G;
	float	boxsize	=	dSimParam.boxsize;
	float	boxhalf	=	boxsize / 2.0;

/*
	printf("device parameter: \n");
	printf("rcut: %g\n", rcut);
	printf("rcut2: %g\n", rcut2);
	printf("asmthfac: %g\n", asmthfac);
	printf("mass: %g\n", mass);
	printf("ErrTolTheta: %g\n", ErrTolTheta);
	printf("G: %g\n", G);
	printf("boxsize: %g\n", boxsize);
	printf("boxhalf: %g\n", boxhalf);
*/
	float	h	=	dSoftParam.h;
	float	h_inv	=	dSoftParam.h_inv;
	float	h3_inv	=	dSoftParam.h3_inv;
/*
	printf("h: %g\n", h);
	printf("h_inv: %g\n", h_inv);
	printf("h3_inv: %g\n", h3_inv);
*/
	acc.x	=	0.0;
	acc.y	=	0.0;
	acc.z	=	0.0;

	pos	=	dPos[index];
	
	no	=	numParticles;

	while(no >= 0){

		if(no < numParticles){

//			pos_no	=	dPos[no];
			pos_no	=	tex1Dfetch(dPosTex, no);
			
			dx	=	pos_no.x - pos.x;
			dy	=	pos_no.y - pos.y;
			dz	=	pos_no.z - pos.z;

			dx	=	NEAREST(dx);
			dy	=	NEAREST(dy);
			dz	=	NEAREST(dz);

			r2	=	dx * dx + dy * dy + dz * dz;

			mass	=	dSimParam.mass;

//			no	=	dNextnode[no];
			no	=	tex1Dfetch(dNextnodeTex, no);
		
		}else{

//			node	=	dNodes[no];
//			node	=	tex1Dfetch(dNodesTex, no);

			temp	=	tex1Dfetch(dNodesTex, 3 * (no - numParticles));
			node_1	=	*(NODE_1 *) &temp;

			temp	=	tex1Dfetch(dNodesTex, 3 * (no - numParticles) + 1);
			node_2	=	*(NODE_2 *) &temp;

			temp	=	tex1Dfetch(dNodesTex, 3 * (no - numParticles) + 2);
			node_3	=	*(NODE_3 *) &temp;

			mass	=	node_2.mass;
			
			dx	=	node_2.s[0] - pos.x;
			dy	=	node_2.s[1] - pos.y;
			dz	=	node_2.s[2] - pos.z;

			dx	=	NEAREST(dx);
			dy	=	NEAREST(dy);
			dz	=	NEAREST(dz);

			r2	=	dx * dx + dy * dy + dz * dz;

			if(r2 > rcut2){

				eff_dist	=	rcut + 0.5 * node_1.len;

				dist	=	NEAREST(node_1.center[0] - pos.x);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	node_3.sibling;
					continue;
				}

				dist	=	NEAREST(node_1.center[1] - pos.y);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	node_3.sibling;
					continue;
				}

				dist	=	NEAREST(node_1.center[2] - pos.z);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	node_3.sibling;
					continue;
				}
			}
			
			if(node_1.len * node_1.len > r2 * ErrTolTheta * ErrTolTheta){
				no	=	node_3.nextnode;
				continue;
			}

			no	=	node_3.sibling;
		}

		r	=	__fsqrt_rn(r2);

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

		tabindex	=	(int) (asmthfac * r);

		if(tabindex < NTAB){

			fac	*=	dSimParam.shortrange_table[tabindex];

			acc.x	+=	dx * fac;
			acc.y	+=	dy * fac;
			acc.z	+=	dz * fac;
		}
	}
	
	acc.x	*=	G;
	acc.y	*=	G;
	acc.z	*=	G;

	dGravAccel[index]	=	acc;

}

#endif
