#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	<cuda.h>
#include	<cuda_runtime.h>
#include	<cufft.h>
//#include	<cutil_inline.h>

#include	"allvars.h"
#include	"proto.h"

#define	PMGRID2	(2*(PMGRID/2+1))

static	cufftHandle	cufft_forward_plan, cufft_inverse_plan;

static	int	fftsize, gridsize;

static	cufftReal		*rhogrid, *forcegrid;
static	cufftComplex	*fft_of_rhogrid;
static	cufftComplex	*data;

static	FLOAT	to_slab_fac;

void pm_init_periodic(void){

	All.Asmth	=	ASMTH * All.BoxSize / PMGRID;
	All.Rcut	=	RCUT * All.Asmth;

	fftsize		=	PMGRID * PMGRID * PMGRID2;

	cudaMalloc((void **)&data, sizeof(cufftComplex) * (fftsize / 2));

	cufftPlan3d(&cufft_forward_plan, PMGRID, PMGRID, PMGRID, CUFFT_R2C);
	cufftPlan3d(&cufft_inverse_plan, PMGRID, PMGRID, PMGRID, CUFFT_C2R);

/* Since we will make a in-place fft, the z dim should be PMGRID2
 */

	gridsize	=	PMGRID * PMGRID * PMGRID;

	to_slab_fac	=	PMGRID / All.BoxSize;
}

void pm_init_periodic_allocate(void){

		if(!(rhogrid = (cufftReal *) malloc(fftsize * sizeof(cufftReal)))){
			printf("failed to allocate memory for FFT-rhogrid!\n");
			exit(1);
		}

		if(!(forcegrid = (cufftReal *) malloc(gridsize * sizeof(cufftReal)))){
			printf("failed to allocate memory for FFT-forcegrid!\n");
			exit(1);
		}

		fft_of_rhogrid	=	(cufftComplex *) & rhogrid[0];
}

void pm_init_periodic_free(void){

	free(forcegrid);
	free(rhogrid);

}

void pmforce_periodic(void){
	
	double	k2, kx, ky, kz, smth;
	double	dx, dy, dz;
	double	fx, fy, fz, ff;
	double	asmth2, fac, acc_dim;

	int	i, ip;
	int	x, y, z, xl, yl, zl, xr, yr, zr, xll, yll, zll, xrr, yrr, zrr, dim;
	int	slab_x, slab_y, slab_z;
	int	slab_xx, slab_yy, slab_zz;
	int	dimx, dimy, dimz;

	printf("Starting long range force calculation...");

	asmth2	=	(2 * M_PI) * All.Asmth / All.BoxSize;
	asmth2	*=	asmth2;

	fac	=	All.G / (M_PI * All.BoxSize);
	fac	*=	1 / (2 * All.BoxSize / PMGRID);
	
	dimx	=	PMGRID;
	dimy	=	PMGRID;
	dimz	=	PMGRID2;

	pm_init_periodic_allocate();

	for(i = 0; i < dimx * dimy * dimz; ++i)
		rhogrid[i]	=	0;
	
	for(i = 0; i < NumPart; ++i){

		slab_x	=	to_slab_fac * P[i].Pos[0];
		if(slab_x >= PMGRID)
			slab_x	=	PMGRID -1;
		dx	=	to_slab_fac * P[i].Pos[0] - slab_x;
		slab_xx	=	slab_x + 1;
		if(slab_xx >= PMGRID)
			slab_xx	-=	PMGRID;
			
		slab_y	=	to_slab_fac * P[i].Pos[1];
		if(slab_y >= PMGRID)
			slab_y	=	PMGRID -1;
		dy	=	to_slab_fac * P[i].Pos[1] - slab_y;
		slab_yy	=	slab_y + 1;
		if(slab_yy >= PMGRID)
			slab_yy	-=	PMGRID;
		
		slab_z	=	to_slab_fac * P[i].Pos[2];
		if(slab_z >= PMGRID)
			slab_z	=	PMGRID -1;
		dz	=	to_slab_fac * P[i].Pos[2] - slab_z;
		slab_zz	=	slab_z + 1;
		if(slab_zz >= PMGRID)
			slab_zz	-=	PMGRID;

		rhogrid[(slab_x * dimy + slab_y) * dimz + slab_z]	+=	P[i].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
		rhogrid[(slab_x * dimy + slab_yy) * dimz + slab_z]	+=	P[i].Mass * (1.0 - dx) * dy * (1.0 - dz);
		rhogrid[(slab_x * dimy + slab_y) * dimz + slab_zz]	+=	P[i].Mass * (1.0 - dx) * (1.0 - dy) * dz;
		rhogrid[(slab_x * dimy + slab_yy) * dimz + slab_zz]	+=	P[i].Mass * (1.0 - dx) * dy * dz;

		rhogrid[(slab_xx * dimy + slab_y) * dimz + slab_z]	+=	P[i].Mass * dx * (1.0 - dy) * (1.0 - dz);
		rhogrid[(slab_xx * dimy + slab_yy) * dimz + slab_z]	+=	P[i].Mass * dx * dy * (1.0 - dz);
		rhogrid[(slab_xx * dimy + slab_y) * dimz + slab_zz]	+=	P[i].Mass * dx * (1.0 - dy) * dz;
		rhogrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz]+=	P[i].Mass * dx * dy * dz;
	}

/* Do the FFT of the density field, since we do not use mpi, the workspace is not necessary, here we make a 
 * real to complex in place fft, the results are stored in rhogrid.
 */

	cudaMemcpy(data, rhogrid, fftsize * sizeof(cufftReal), cudaMemcpyHostToDevice);

	cufftExecR2C(cufft_forward_plan, (cufftReal *) data, data);

	cudaMemcpy(rhogrid, data, fftsize * sizeof(cufftReal), cudaMemcpyDeviceToHost);

/* Multiply with Green's function
 */
	for(x = 0; x < PMGRID; ++x){
		for(y = 0; y < PMGRID; ++y){
			for(z = 0; z < PMGRID / 2 + 1; ++z){
				if(x > PMGRID / 2)
					kx	=	x - PMGRID;
				else
					kx	=	x;
				if(y > PMGRID / 2)
					ky	=	y - PMGRID;
				else
					ky	=	y;
				if(z > PMGRID / 2)
					kz	=	z - PMGRID;
				else
					kz	=	z;

				k2	=	kx * kx + ky * ky + kz * kz;

				if(k2 > 0){

					smth	=	-exp(-k2 * asmth2) / k2;

					fx = fy = fz = 1;
					
					if(kx != 0){
						fx	=	(M_PI * kx) / PMGRID;
						fx	=	sin(fx) / fx;
					}

					if(ky != 0){
						fy	=	(M_PI * ky) / PMGRID;
						fy	=	sin(fy) / fy;
					}

					if(kz != 0){
						fz	=	(M_PI * kz) / PMGRID;
						fz	=	sin(fz) / fz;
					}

					ff	=	1 / (fx * fy * fz);
					smth	*=	ff * ff * ff * ff;

					ip	=	PMGRID * (PMGRID / 2 + 1) * x + (PMGRID / 2 + 1) * y + z;

					fft_of_rhogrid[ip].x	*=	smth;
					fft_of_rhogrid[ip].y	*=	smth;
				}
			}
		}
	}
	
	fft_of_rhogrid[0].x	=	fft_of_rhogrid[0].y	=	0.0;

	cudaMemcpy(data, rhogrid, fftsize * sizeof(cufftReal), cudaMemcpyHostToDevice);

	cufftExecC2R(cufft_inverse_plan, data, (cufftReal *) data);

	cudaMemcpy(rhogrid, data, fftsize * sizeof(cufftReal), cudaMemcpyDeviceToHost);

	for(dim = 0; dim < 3; dim++){
		for(x = 0; x < PMGRID; ++x){
			for(y = 0; y < PMGRID; ++y){
				for(z = 0; z < PMGRID; ++z){
					xrr = xll = xr = xl = x;
					yrr = yll = yr = yl = y;
					zrr = zll = zr = zl = z;

					switch(dim){
						case 0:
							xr	=	x + 1;
							if(xr >= PMGRID)
								xr	-=	PMGRID;
							xrr	=	x + 2;
							if(xrr >= PMGRID)
								xrr	-=	PMGRID;
							xl	=	x - 1;
							if(xl < 0)
								xl	+=	PMGRID;
							xll	=	x - 2;
							if(xll < 0)
								xll +=	PMGRID;
							break;

						case 1:
							yr	=	y + 1;
							if(yr >= PMGRID)
								yr	-=	PMGRID;
							yrr	=	y + 2;
							if(yrr >= PMGRID)
								yrr	-=	PMGRID;
							yl	=	y - 1;
							if(yl < 0)
								yl	+=	PMGRID;
							yll	=	y - 2;
							if(yll < 0)
								yll +=	PMGRID;
							break;
							
						case 2:
							zr	=	z + 1;
							if(zr >= PMGRID)
								zr	-=	PMGRID;
							zrr	=	z + 2;
							if(zrr >= PMGRID)
								zrr	-=	PMGRID;
							zl	=	z - 1;
							if(zl < 0)
								zl	+=	PMGRID;
							zll	=	z - 2;
							if(zll < 0)
								zll +=	PMGRID;
							break;
					}

					forcegrid[(x * PMGRID + y) * PMGRID + z]
						=	
						fac * ((4.0 / 3) *
							   (rhogrid[(xl * dimy + yl) * dimz + zl] 
								- rhogrid[(xr * dimy + yr) * dimz + zr])
							   -(1.0 / 6) *
							   (rhogrid[(xll * dimy + yll) * dimz + zll] 
							    - rhogrid[(xrr * dimy + yrr) * dimz + zrr]));
				}
			}
		}

		for(i = 0; i < NumPart; ++i){
			slab_x	=	to_slab_fac * P[i].Pos[0];
			if(slab_x >= PMGRID)
				slab_x	=	PMGRID -1;
			dx	=	to_slab_fac * P[i].Pos[0] - slab_x;
			slab_xx	=	slab_x + 1;
			if(slab_xx >= PMGRID)
				slab_xx	-=	PMGRID;

			slab_y	=	to_slab_fac * P[i].Pos[1];
			if(slab_y >= PMGRID)
				slab_y	=	PMGRID -1;
			dy	=	to_slab_fac * P[i].Pos[1] - slab_y;
			slab_yy	=	slab_y + 1;
			if(slab_yy >= PMGRID)
				slab_yy	-=	PMGRID;

			slab_z	=	to_slab_fac * P[i].Pos[2];
			if(slab_z >= PMGRID)
				slab_z	=	PMGRID -1;
			dz	=	to_slab_fac * P[i].Pos[2] - slab_z;
			slab_zz	=	slab_z + 1;
			if(slab_zz >= PMGRID)
				slab_zz	-=	PMGRID;

			acc_dim	=	
				forcegrid[(slab_x * PMGRID + slab_y) * PMGRID + slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
			acc_dim	+=	forcegrid[(slab_x * PMGRID + slab_yy) * PMGRID + slab_z] * (1.0 - dx) * dy * (1.0 - dz);
			acc_dim	+=	forcegrid[(slab_x * PMGRID + slab_y) * PMGRID + slab_zz] * (1.0 - dx) * (1.0 - dy) * dz;
			acc_dim	+=	forcegrid[(slab_x * PMGRID + slab_yy) * PMGRID + slab_zz] * (1.0 - dx) * dy * dz;

			acc_dim	+=	forcegrid[(slab_xx * PMGRID + slab_y) * PMGRID + slab_z] * (dx) * (1.0 - dy) * (1.0 - dz);
			acc_dim	+=	forcegrid[(slab_xx * PMGRID + slab_yy) * PMGRID + slab_z] * (dx) * dy * (1.0 - dz);
			acc_dim	+=	forcegrid[(slab_xx * PMGRID + slab_y) * PMGRID + slab_zz] * (dx) * (1.0 - dy) * dz;
			acc_dim	+=	forcegrid[(slab_xx * PMGRID + slab_yy) * PMGRID + slab_zz] * (dx) * dy * dz;

			P[i].GravPM[dim]	=	acc_dim;
		}
	}
	
	pm_init_periodic_free();
	printf("done!\n");
}					


