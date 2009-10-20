#ifndef	_GRAVTREE_KERNEL_H_
#define	_GRAVTREE_KERNEL_H_

#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))

__constant__ SIMPARAM	dSimParam;

__constant__ SOFTPARAM	dSoftParam;

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

	NODE	node;
	int		no, tabindex;
	float	r2, dx, dy, dz, r, fac, u;
	float	eff_dist, dist;
	float4	pos, pos_no, acc;
	
	float	rcut	=	dSimParam.rcut;
	float	rcut2	=	dSimParam.rcut2;
	float	asmthfac	=	dSimParam.asmthfac;
	float	mass	=	dSimParam.mass;
	float	ErrTolTheta	=	dSimParam.ErrTolTheta;
	float	G		=	dSimParam.G;
	float	boxsize	=	dSimParam.boxsize;
	float	boxhalf	=	boxsize / 2.0;
	
	float	h	=	dSoftParam.h;
	float	h_inv	=	dSoftParam.h_inv;
	float	h3_inv	=	dSoftParam.h3_inv;

	acc.x	=	0.0;
	acc.y	=	0.0;
	acc.z	=	0.0;

	pos	=	dPos[index];
	
	no	=	numParticles;

/*
	int	i;
	i	=	no + 250;
	printf("no: %g|%g|%g\t%g\t%d\t%d\t%d\t%d\n", 
			dNodes[i].u.d.s[0], 
			dNodes[i].u.d.s[1], 
			dNodes[i].u.d.s[2], 
			dNodes[i].u.d.mass, 
			dNodes[i].u.d.bitflags, 
			dNodes[i].u.d.sibling,
			dNodes[i].u.d.nextnode,
			dNodes[i].u.d.father);
*/

	while(no >= 0){
		
		if(no < numParticles){

			pos_no	=	dPos[no];

			dx	=	pos_no.x - pos.x;
			dy	=	pos_no.y - pos.y;
			dz	=	pos_no.z - pos.z;

			dx	=	NEAREST(dx);
			dy	=	NEAREST(dy);
			dz	=	NEAREST(dz);

			r2	=	dx * dx + dy * dy + dz * dz;


			no	=	dNextnode[no];
		
		}else{

			node	=	dNodes[no];

			mass	=	node.u.d.mass;
			
			dx	=	node.u.d.s[0] - pos.x;
			dy	=	node.u.d.s[1] - pos.y;
			dz	=	node.u.d.s[2] - pos.z;

			dx	=	NEAREST(dx);
			dy	=	NEAREST(dy);
			dz	=	NEAREST(dz);

			r2	=	dx * dx + dy * dy + dz * dz;

			if(r2 > rcut2){

				eff_dist	=	rcut + 0.5 * node.len;

				dist	=	NEAREST(node.center[0] - pos.x);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	node.u.d.sibling;
					continue;
				}

				dist	=	NEAREST(node.center[1] - pos.y);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	node.u.d.sibling;
					continue;
				}

				dist	=	NEAREST(node.center[2] - pos.z);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	node.u.d.sibling;
					continue;
				}
			}
			
			if(node.len * node.len > r2 * ErrTolTheta * ErrTolTheta){
				no	=	node.u.d.nextnode;
				continue;
			}

			no	=	node.u.d.sibling;
		}

		r	=	sqrtf(r2);

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
