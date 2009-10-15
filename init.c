#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

void init(void){
	
	double	a3;
	int		i,j;

	read_ic(All.InitCondFile);

	All.Time	=	All.TimeBegin;
	All.Ti_Current	=	0;
	All.Timebase_interval	=	(log(All.TimeMax)-log(All.TimeBegin)) / TIMEBASE;

	a3	=	All.Time * All.Time * All.Time;

	set_softenings();

	All.NumCurrentTiStep	=	0;
	All.SnapshotFileCount	=	0;

	check_omega();

	for(i = 0; i < NumPart; ++i){
		for(j = 0; j < 3; ++j){
			P[i].Vel[j]	*=	sqrt(All.Time)*All.Time;
			P[i].GravAccel[j]	=	0;
			P[i].GravPM[j]	=	0;
		}
		P[i].Ti_endstep	=	0;
		P[i].Ti_begstep	=	0;
	}

	All.PM_Ti_endstep	=	All.PM_Ti_begstep	=	0;

#ifdef	TREE
	force_treeallocate(0.8 * NumPart, NumPart);
#endif
	TreeReconstructFlag	=	1;

}

void check_omega(void){
	
	double	mass, omega;
	int		i;
	
	mass	=	0.0;
	
	for(i = 0; i < NumPart; ++i){
		mass	+=	P[i].Mass;
	}

	omega	=	mass / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

	if(fabs(omega - All.Omega0) > 0.001){
		printf("The mass in ic file accounts for Omega = %g, but you specified Omega = %g in parameter file!\n", omega, All.Omega0);
		exit(1);
	}
}
