#ifndef	_GRAVTREE_KERNEL_H_
#define	_GRAVTREE_KERNEL_H_

#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))

texture<float4, 1, cudaReadModeElementType> dPosTex;
texture<float4, 1, cudaReadModeElementType> dNodesTex;
texture<int, 1, cudaReadModeElementType> dNextnodeTex;

__constant__ SIMPARAM	dSimParam;

__constant__ SOFTPARAM	dSoftParam;

__device__	float4	*Pos;

//__global__ void force_treebuild_device(
//					float4)

__global__ void force_treeevaluate_shortrange_device(
						float4	*dPos,
						float4	*dGravAccel,
						NODE	*dNodes,
						EXTNODE	*dExtnodes,
						int		*dNextnode,
						int		*dFather,
						int 	numParticles){

	uint	index	=	blockIdx.y * gridDim.x * blockDim.x
						+ blockIdx.x * blockDim.x + threadIdx.x;
	if(index >= numParticles){
		return;
	}
	
	Pos	=	dPos;

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

	pos	=	Pos[index];
	
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
