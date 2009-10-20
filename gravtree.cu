#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	<cuda.h>
#include	<cuda_runtime.h>

#include	"allvars.h"
#include	"proto.h"

#include	"gravtree_kernel.cu"


extern	"C"{

static	int	maxThreads	=	32768;
static	int	numThreads	=	256;

static	SIMPARAM	hSimParam;

static	SOFTPARAM	hSoftParam;

static	int		last;

static	int		numnodes;

static float	shortrange_table[NTAB];

static float	boxsize, boxhalf;

static double	to_grid_fac;

static double	rcut, rcut2, asmth, asmthfac;

static double	h, h_inv, h3_inv;

void cudaInit(){
	cudaSetDevice(0);
}

void set_softenings(void){

	if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys){
		All.SofteningTable	=	All.SofteningHaloMaxPhys / All.Time;
	}else{
		All.SofteningTable	=	All.SofteningHalo;
	}

	All.ForceSoftening	=	2.8 * All.SofteningTable;
	
	h	=	All.ForceSoftening;
	h_inv	=	1.0 / h;
	h3_inv	=	h_inv * h_inv * h_inv;

	hSoftParam.h	=	h;
	hSoftParam.h_inv	=	h_inv;
	hSoftParam.h3_inv	=	h3_inv;
}


void gravity_tree(void){

	set_softenings();

	cudaMemcpyToSymbol(dSoftParam, &hSoftParam, sizeof(SOFTPARAM));

	printf("Tree building\n");

	if(TreeReconstructFlag){

		force_treebuild(NumPart);

//		TreeReconstructFlag	=	0;
	}

	printf("Copy data to device\n");
	copyTreeToDevice();
	copyPosToDevice();

	dim3	dimBlock(numThreads, 1);
	dim3	dimGrid;

	if(NumPart > maxThreads){
		dimGrid.y	=	NumPart / maxThreads;
		dimGrid.x	=	maxThreads / numThreads;
	}else{
		dimGrid.y	=	1;
		dimGrid.x	=	NumPart / numThreads;
	}

	printf("calc force\n");
	force_treeevaluate_shortrange_device<<<dimGrid, dimBlock>>>((float4 *)dPos, (float4 *)dGravAccel, dNodes, dExtnodes, dNextnode, dFather, NumPart);

	printf("copy back\n");
	copyAccelFromDevice();

	printf("done!\n");

	cudaError_t cudaError;
	cudaError = cudaGetLastError();
	if( cudaError != cudaSuccess )
	{
		fprintf(stderr, "CUDA Runtime API Error reported : %s\n", cudaGetErrorString(cudaError));
		exit(EXIT_FAILURE);
	}


}

void force_treeevaluate_shortrange(int target){

	NODE	* nop;
	int		no, tabindex;
	double	r2, dx, dy, dz, mass, r, fac, u;
	double	acc_x, acc_y, acc_z, pos_x, pos_y, pos_z;
	double	eff_dist, dist;

	int	count	=	0;

	acc_x	=	0;
	acc_y	=	0;
	acc_z	=	0;

	pos_x	=	P[target].Pos[0];
	pos_y	=	P[target].Pos[1];
	pos_z	=	P[target].Pos[2];

	no	=	All.NumPart;	/* root node */
	while(no >= 0){

		if(no < All.NumPart){

			dx	=	P[no].Pos[0] - pos_x;
			dy	=	P[no].Pos[1] - pos_y;
			dz	=	P[no].Pos[2] - pos_z;

			dx	=	NEAREST(dx);
			dy	=	NEAREST(dy);
			dz	=	NEAREST(dz);

			r2	=	dx * dx + dy * dy + dz * dz;

			mass	=	P[no].Mass;

			no	=	Nextnode[no];

		}else{
			
			nop	=	&Nodes[no];

			mass	=	nop->u.d.mass;

			dx	=	nop->u.d.s[0] - pos_x;
			dy	=	nop->u.d.s[1] - pos_y;
			dz	=	nop->u.d.s[2] - pos_z;

			dx	=	NEAREST(dx);
			dy	=	NEAREST(dy);
			dz	=	NEAREST(dz);

			r2	=	dx * dx + dy * dy + dz * dz;

			if(r2 > rcut2){

				eff_dist	=	rcut + 0.5 * nop->len;

				dist	=	NEAREST(nop->center[0] - pos_x);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	nop->u.d.sibling;
					continue;
				}

				dist	=	NEAREST(nop->center[1] - pos_y);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	nop->u.d.sibling;
					continue;
				}
	
				dist	=	NEAREST(nop->center[2] - pos_z);
				if(dist < -eff_dist || dist > eff_dist){
					no	=	nop->u.d.sibling;
					continue;
				}
			}
			
			if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta){
				no	=	nop->u.d.nextnode;
				continue;
			}

			no	=	nop->u.d.sibling;
		}

					
		r	=	sqrt(r2);
					
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

			fac	*=	shortrange_table[tabindex];

			acc_x	+=	dx * fac;
			acc_y	+=	dy * fac;
			acc_z	+=	dz * fac;
		}

		count++;
	}
	
	P[target].GravAccel[0]  =   acc_x * All.G;
	P[target].GravAccel[1]  =   acc_y * All.G;
	P[target].GravAccel[2]  =   acc_z * All.G;
	
}

void force_treeallocate(int maxnodes, int maxpart){

	MaxNodes	=	maxnodes;
	
	Nodes_base	=	(NODE_BASE *) malloc((MaxNodes + 1) * sizeof(NODE_BASE));

	Extnodes_base	=	(EXTNODE_BASE *) malloc((MaxNodes + 1) * sizeof(EXTNODE_BASE));

	Nodes		=	Nodes_base - NumPart;

	Extnodes	=	Extnodes_base - NumPart;

	Nextnode	=	(int *) malloc(NumPart * sizeof(int));

	Father	=	(int *)malloc(NumPart * sizeof(int));

	int		i;
	double	u;

	for(i = 0; i < NTAB; ++i){
		u	=	3.0 / NTAB * (i + 0.5);
		shortrange_table[i]	=	erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
	}

	printf("boxsize in cu: %g\n", All.BoxSize);

	to_grid_fac	=	PMGRID / All.BoxSize;

	rcut	=	All.Rcut;
	rcut2	=	rcut * rcut;

	asmth		=	All.Asmth;
	asmthfac	=	0.5 / asmth * (NTAB / 3.0);

	boxsize	=	All.BoxSize;
	boxhalf	=	boxsize / 2.0;

	printf("to_grid_fac: %g\tasmth: %g\trcut: %g\n", 
			to_grid_fac, All.Asmth, All.Rcut);

	hPos	=	(float *)malloc(sizeof(float) * 4 * NumPart);
	if(hPos == NULL){
		printf("failed to allocate memory for hPos!\n");
		exit(1);
	}

	hGravAccel	=	(float *)malloc(sizeof(float) * 4 * NumPart);

	int	memSize	=	sizeof(float) * 4 * NumPart;

	cudaMalloc((void **) &dPos, memSize);
	cudaMalloc((void **) &dGravAccel, memSize);

	cudaMalloc((void **) &dNodes, sizeof(NODE) * MaxNodes);

	dNodes	-=	NumPart;

	cudaMalloc((void **) &dExtnodes, sizeof(EXTNODE) * MaxNodes);
	dExtnodes	-=	NumPart;

	cudaMalloc((void **) &dNextnode, sizeof(int) * NumPart);
	cudaMalloc((void **) &dFather, sizeof(int) * NumPart);

	for(i = 0; i < NTAB; ++i){
		hSimParam.shortrange_table[i]	=	shortrange_table[i];
	}

	hSimParam.boxsize	=	boxsize;
	hSimParam.rcut		=	rcut;
	hSimParam.rcut2		=	rcut2;
	hSimParam.asmthfac	=	asmthfac;
	hSimParam.mass		=	header.Mass;
	hSimParam.ErrTolTheta	=	All.ErrTolTheta;
	hSimParam.G			=	All.G;

	printf("setting param...\n");

	printf("ErrTolTheta: %g\n", hSimParam.ErrTolTheta);

	cudaMemcpyToSymbol(dSimParam, &hSimParam, sizeof(SIMPARAM));

	printf("param done\n");

}

void force_treebuild(int npart){
	
	int		i, subnode, parent;
	int		nfree, th, nn;
	NODE	* nfreep;
	double	lenhalf;
	
	subnode	=	0;

/*	create empty root node, size is just boxsize
 */
	nfree	=	NumPart;	/* index of first node */
	nfreep	=	&Nodes[nfree];	/* pointer of first node */

	nfreep->len	=	boxsize;	
	for(i = 0; i < 3; ++i)
		nfreep->center[i]	=	boxhalf;
	for(i = 0; i < 8; ++i)
		nfreep->u.suns[i]	=	-1;
	
	numnodes	=	1;
	nfreep++;
	nfree++;

	nfreep	=	&Nodes[nfree];
	parent	=	-1;

/*	next insert all particles one by one
 */
	for(i = 0; i < NumPart; ++i){
		
		th	=	NumPart;

		while(1){

			if(th >= NumPart){
			/* This is a node, continue.
			 */
				subnode	=	0;
				if(P[i].Pos[0] > Nodes[th].center[0])
					subnode	+=	1;
				if(P[i].Pos[1] > Nodes[th].center[1])
					subnode	+=	2;
				if(P[i].Pos[2] > Nodes[th].center[2])
					subnode	+=	4;

				nn	=	Nodes[th].u.suns[subnode];

				if(nn >= 0){
					parent	=	th;
					th	=	nn;
				}else{
					/* OK! We have find an empty slot, insert particle i here.
					 */
					Nodes[th].u.suns[subnode]	=	i;
					break;
				}
			}else{
			/* This is a particle, we have to generate a new node and insert 
			 * the old one
			 */
				Nodes[parent].u.suns[subnode]	=	nfree;

				nfreep->len	=	0.5 * Nodes[parent].len;
				lenhalf	=	0.25 * Nodes[parent].len;

				if(subnode & 1)
					nfreep->center[0]	=	Nodes[parent].center[0] + lenhalf;
				else
					nfreep->center[0]	=	Nodes[parent].center[0] - lenhalf;

				if(subnode & 2)
					nfreep->center[1]	=	Nodes[parent].center[1] + lenhalf;
				else
					nfreep->center[1]	=	Nodes[parent].center[1] - lenhalf;

				if(subnode & 4)
					nfreep->center[2]	=	Nodes[parent].center[2] + lenhalf;
				else
					nfreep->center[2]	=	Nodes[parent].center[2] - lenhalf;

				/* I think here does not use loop is for the sake of speed
				 */
				nfreep->u.suns[0]	=	-1;
				nfreep->u.suns[1]	=	-1;
				nfreep->u.suns[2]	=	-1;
				nfreep->u.suns[3]	=	-1;
				nfreep->u.suns[4]	=	-1;
				nfreep->u.suns[5]	=	-1;
				nfreep->u.suns[6]	=	-1;
				nfreep->u.suns[7]	=	-1;

				subnode	=	0;
				if(P[th].Pos[0] > nfreep->center[0])
					subnode	+=	1;
				if(P[th].Pos[1] > nfreep->center[1])
					subnode	+=	2;
				if(P[th].Pos[2] > nfreep->center[2])
					subnode	+=	4;

				nfreep->u.suns[subnode]	=	th;

				th	=	nfree;

				numnodes++;
				nfree++;
				nfreep++;

				if(numnodes >= MaxNodes){
					printf("maximun number of tree nodes reached\n");
					printf("we'd better stop!\n");
					printf("numnodes: %d\nMaxNodes: %d\n", numnodes, MaxNodes);
					exit(1);
				}
			}
		}
	}

	/* next computer multipole moments recursively.
	 */
	last	=	-1;
	force_update_node_recursive(All.NumPart, -1, -1);

	printf("last: %d\n", last);
	printf("numnodes: %d\n", numnodes);

	if(last >= All.NumPart)
		Nodes[last].u.d.nextnode	=	-1;
	else
		Nextnode[last]	=	-1;

	i	=	All.NumPart + 250;
	printf("build: %g|%g|%g\t%g\t%d\t%d\t%d\t%d\n", 
			Nodes[i].u.d.s[0], 
			Nodes[i].u.d.s[1], 
			Nodes[i].u.d.s[2], 
			Nodes[i].u.d.mass, 
			Nodes[i].u.d.bitflags, 
			Nodes[i].u.d.sibling,
			Nodes[i].u.d.nextnode,
			Nodes[i].u.d.father);
}

void force_update_node_recursive(int no, int sib, int father){

	int	j, jj, p, pp, nextsib, suns[8];
	PARTICLE	* pa;
	double	s[3], vs[3], mass;

	if(no >= All.NumPart){
	/* This is a node. */
		for(j = 0; j < 8; ++j)
			suns[j]	=	Nodes[no].u.suns[j];
		
		if(last >= 0){

			if(last >= All.NumPart)
				Nodes[last].u.d.nextnode	=	no;
			else
				Nextnode[last]	=	no;
		}

		last	=	no;

		mass	=	0;
		s[0]	=	0;
		s[1]	=	0;
		s[2]	=	0;
		vs[0]	=	0;
		vs[1]	=	0;
		vs[2]	=	0;

		for(j = 0; j < 8; ++j){

			p	=	suns[j];
			if(p >= 0){
					
				for(jj = j + 1; jj < 8; ++jj){
					pp	=	suns[jj];
					if(pp >= 0)
						break;
				}

				if(jj < 8)
					nextsib	=	pp;
				else
					nextsib	=	sib;

				force_update_node_recursive(p, nextsib, no);

				if(p >= All.NumPart){
						
					mass	+=	Nodes[p].u.d.mass;
					s[0]	+=	Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
					s[1]	+=	Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
					s[2]	+=	Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
					vs[0]	+=	Nodes[p].u.d.mass * Extnodes[p].vs[0];
					vs[0]	+=	Nodes[p].u.d.mass * Extnodes[p].vs[1];
					vs[0]	+=	Nodes[p].u.d.mass * Extnodes[p].vs[2];
					
				}else{

					pa	=	&P[p];

					mass	+=	pa->Mass;
					s[0]	+=	pa->Mass * pa->Pos[0];
					s[1]	+=	pa->Mass * pa->Pos[1];
					s[2]	+=	pa->Mass * pa->Pos[2];
					vs[0]	+=	pa->Mass * pa->Vel[0];
					vs[1]	+=	pa->Mass * pa->Vel[1];
					vs[2]	+=	pa->Mass * pa->Vel[2];
				}
			}
		}

		if(mass){
			
			s[0]	/=	mass;
			s[1]	/=	mass;
			s[2]	/=	mass;
			vs[0]	/=	mass;
			vs[1]	/=	mass;
			vs[2]	/=	mass;
		
		}else{

			s[0]	=	Nodes[no].center[0];
			s[1]	=	Nodes[no].center[1];
			s[2]	=	Nodes[no].center[2];

		}

		Nodes[no].u.d.s[0]	=	s[0];
		Nodes[no].u.d.s[1]	=	s[1];
		Nodes[no].u.d.s[2]	=	s[2];
		Nodes[no].u.d.mass	=	mass;

		Extnodes[no].vs[0]	=	vs[0];
		Extnodes[no].vs[1]	=	vs[1];
		Extnodes[no].vs[2]	=	vs[2];

		Nodes[no].u.d.sibling	=	sib;
		Nodes[no].u.d.father	=	father;
	
	}else{
		
		if(last >= 0){
			
			if(last >= All.NumPart)
				Nodes[last].u.d.nextnode	=	no;
			else
				Nextnode[last]	=	no;
		}

		last	=	no;

		if(no < All.NumPart)
			Father[no]	=	father;
	}
}

void copyTreeToDevice(){

	printf("copy numnodes: %d\n", numnodes);

	cudaMemcpy((char *) &dNodes[NumPart], (void *) &Nodes[NumPart], numnodes * sizeof(NODE), cudaMemcpyHostToDevice);
	
	cudaMemcpy((char *) &dExtnodes[NumPart], (void *) &Extnodes[NumPart], numnodes * sizeof(EXTNODE), cudaMemcpyHostToDevice);

	cudaMemcpy((char *) dNextnode, (void *)Nextnode, NumPart * sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy((char *) dFather, (void *)Father, NumPart * sizeof(int), cudaMemcpyHostToDevice);

}

void copyPosToDevice(){
	
	float4	*hPos4	=	(float4 *)hPos;
	
	int	i;
	for(i = 0; i < NumPart; ++i){
		hPos4[i].x	=	P[i].Pos[0];
		hPos4[i].y	=	P[i].Pos[1];
		hPos4[i].z	=	P[i].Pos[2];
	}

	cudaMemcpy((char *)dPos, (void *)hPos, NumPart * sizeof(float4), cudaMemcpyHostToDevice);
}

void copyAccelFromDevice(){

	cudaMemcpy((void *)hGravAccel, (char *)dGravAccel, NumPart * sizeof(float4), cudaMemcpyDeviceToHost);

	float4	*hGravAccel4	=	(float4 *)hGravAccel;

	int	i;
	for(i = 0; i < NumPart; ++i){
		P[i].GravAccel[0]	=	hGravAccel4[i].x;
		P[i].GravAccel[1]	=	hGravAccel4[i].y;
		P[i].GravAccel[2]	=	hGravAccel4[i].z;
	}
}

}
