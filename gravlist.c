#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

#define	NTAB	1000

static float	shortrange_table[NTAB];

static float	boxsize, boxhalf;

static double	to_grid_fac;

static double	rcut, rcut2, asmth, asmthfac;

static double	h, h_inv, h3_inv;

static int *	ll;
static int *	ihoc;

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
}

void linklist_init(void){
	
	int		i;
	double	u;
	
	ll		=	(int *)malloc(sizeof(int) * NumPart);
	ihoc	=	(int *)malloc(sizeof(int) * PMGRID * PMGRID * PMGRID);

	for(i = 0; i < NTAB; ++i){
		u	=	3.0 / NTAB * (i + 0.5);
		shortrange_table[i]	=	erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
	}

	to_grid_fac	=	PMGRID / All.BoxSize;

	rcut	=	All.Rcut;
	rcut2	=	rcut * rcut;

	asmth		=	All.Asmth;
	asmthfac	=	0.5 / asmth * (NTAB / 3.0);

	boxsize	=	All.BoxSize;
	boxhalf	=	boxsize / 2.0;

	printf("Initiate linklist...\n");
	printf("to_grid_fac: %g\tasmth: %g\trcut: %g\n", 
			to_grid_fac, asmth, rcut);
}

void build_linklist(void){
	
	printf("Building linklist...\n");

	int	i, ix, iy, iz, ibox;

	for(i = 0; i < NumPart; ++i)
		ll[i]	=	-1;
	
	for(i = 0; i < PMGRID * PMGRID * PMGRID; ++i){
		ihoc[i]	=	-1;
	}
	
	for(i = 0; i < NumPart; ++i){

		ix	=	to_grid_fac * P[i].Pos[0];
		iy	=	to_grid_fac * P[i].Pos[1];
		iz	=	to_grid_fac * P[i].Pos[2];

		ibox	=	ix * PMGRID * PMGRID + iy * PMGRID + iz;

		ll[i]		=	ihoc[ibox];
		ihoc[ibox]	=	i;
	}
	
	ibox	=	PMGRID * PMGRID * PMGRID / 2;

	printf("done!\n");
}

void gravity_linklist(){

	set_softenings();

	build_linklist();

	int	i;

	printf("Calculating short range force...\n");

	for(i = 0; i < NumPart; ++i){

		if(P[i].Ti_endstep == All.Ti_Current){
			force_evaluate_shortrange(i);
		}
	}

	printf("done!\n");
}

void force_evaluate_shortrange(int target){
	
	double	acc_x, acc_y, acc_z, pos_x, pos_y, pos_z;
	double	r, r2, dx, dy, dz, fac, u, mass;
	int	xl, xr, yl, yr, zl, zr, ix, iy, iz, ibox, tabindex;
	int	no, iix, iiy, iiz;

	acc_x	=	0;
	acc_y	=	0;
	acc_z	=	0;

	pos_x	=	P[target].Pos[0];
	pos_y	=	P[target].Pos[1];
	pos_z	=	P[target].Pos[2];

	mass	=	P[target].Mass;

	xl	=	to_grid_fac * (pos_x - rcut) + PMGRID;
	xl	-=	PMGRID;
	xr	=	to_grid_fac * (pos_x + rcut);

	yl	=	to_grid_fac * (pos_y - rcut) + PMGRID;
	yl	-=	PMGRID;
	yr	=	to_grid_fac * (pos_y + rcut);

	zl	=	to_grid_fac * (pos_z - rcut) + PMGRID;
	zl	-=	PMGRID;
	zr	=	to_grid_fac * (pos_z + rcut);

	for(ix = xl; ix <= xr; ++ix){
		for(iy = yl; iy <= yr; ++iy){
			for(iz = zl; iz <= zr; ++iz){

				iix	=	ix;
				iiy	=	iy;
				iiz	=	iz;

				while(iix < 0)
					iix	+=	PMGRID;
				while(iix >= PMGRID)
					iix	-=	PMGRID;

				while(iiy < 0)
					iiy	+=	PMGRID;
				while(iiy >= PMGRID)
					iiy	-=	PMGRID;
					
				while(iiz < 0)
					iiz	+=	PMGRID;
				while(iiz >= PMGRID)
					iiz	-=	PMGRID;

				ibox	=	iix * PMGRID * PMGRID + iiy * PMGRID + iiz;
			
				no	=	ihoc[ibox];
				
				while(no != -1){

					if(no == target){
						no	=	ll[no];
						continue;
					}

					dx	=	P[no].Pos[0] - pos_x;
					dy	=	P[no].Pos[1] - pos_y;
					dz	=	P[no].Pos[2] - pos_z;
					
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
					
					r2	=	dx * dx + dy * dy + dz * dz;

					if(r2 > rcut2){
						no	=	ll[no];
						continue;
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

					no	=	ll[no];
				}
			}
		}
	}
	
	P[target].GravAccel[0]	=	acc_x * All.G;
	P[target].GravAccel[1]	=	acc_y * All.G;
	P[target].GravAccel[2]	=	acc_z * All.G;
}
