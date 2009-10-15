#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	<rfftw.h>

#include	"allvars.h"
#include	"proto.h"

#define	PMGRID2	(2*(PMGRID/2+1))

static	rfftwnd_plan	fft_forward_plan, fft_inverse_plan;

static	int	fftsize, gridsize;

static	fftw_real	*rhogrid, *forcegrid;
static	fftw_complex	*fft_of_rhogrid;

static	FLOAT	to_slab_fac;

void pm_init_periodic(void){

	All.Asmth	=	ASMTH * All.BoxSize / PMGRID;
	All.Rcut	=	RCUT * All.Asmth;

	fft_forward_plan = rfftw3d_create_plan(PMGRID, PMGRID, PMGRID,
						FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
	fft_inverse_plan = rfftw3d_create_plan(PMGRID, PMGRID, PMGRID,
						FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

/* Since we will make a in-place fft, the z dim should be PMGRID2
 */
	fftsize		=	PMGRID * PMGRID * PMGRID2;

	gridsize	=	PMGRID * PMGRID * PMGRID;

	to_slab_fac	=	PMGRID / All.BoxSize;
}

void pm_init_periodic_allocate(void){

		if(!(rhogrid = (fftw_real *) malloc(fftsize * sizeof(fftw_real)))){
			printf("failed to allocate memory for FFT-rhogrid!\n");
			exit(1);
		}

		if(!(forcegrid = (fftw_real *) malloc(gridsize * sizeof(fftw_real)))){
			printf("failed to allocate memory for FFT-forcegrid!\n");
			exit(1);
		}

		fft_of_rhogrid	=	(fftw_complex *) & rhogrid[0];
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

	rfftwnd_one_real_to_complex(fft_forward_plan, rhogrid, NULL);

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

					fft_of_rhogrid[ip].re	*=	smth;
					fft_of_rhogrid[ip].im	*=	smth;
				}
			}
		}
	}
	
	fft_of_rhogrid[0].re	=	fft_of_rhogrid[0].im	=	0.0;

	rfftwnd_one_complex_to_real(fft_inverse_plan, fft_of_rhogrid, NULL);

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


