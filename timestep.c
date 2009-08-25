#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

static double	fac1, hubble_a, atime;
static double	dt_displacement;

void advance_and_find_timesteps(void){

	int		i, j, ti_step, ti_min, tend, tstart;
	double	dt_gravkick;
	float	dv[3];

	fac1	=	1 / (All.Time * All.Time);
	atime	=	All.Time;
	hubble_a=	All.Omega0 / (atime * atime * atime) 
				+ (1 - All.Omega0 - All.OmegaLambda) 
				/ (atime * atime) + All.OmegaLambda;
	hubble_a	=	All.Hubble * sqrt(hubble_a);

	dt_displacement	=	All.MaxSizeTimestep;

	for(i = 0; i < NumPart; ++i){

		if(P[i].Ti_endstep == All.Ti_Current){

			ti_step =	get_timestep(i);
			
			ti_min	=	TIMEBASE;
			while(ti_min > ti_step)
				ti_min	>>= 1;
			ti_step	=	ti_min;

			if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep)){
				if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
					ti_step	=	P[i].Ti_endstep - P[i].Ti_begstep;
			}
			
			if(All.Ti_Current == TIMEBASE)
				ti_step	=	0;

			if((TIMEBASE - All.Ti_Current) < ti_step)
				ti_step	=	TIMEBASE - All.Ti_Current;

			tstart	=	(P[i].Ti_begstep + P[i].Ti_endstep) / 2;
			tend	=	P[i].Ti_endstep + ti_step / 2;

			dt_gravkick	=	get_gravkick_factor(tstart, tend);

			P[i].Ti_begstep	=	P[i].Ti_endstep;
			P[i].Ti_endstep	=	P[i].Ti_begstep + ti_step;

			for(j = 0; j < 3; ++j){
				dv[j]	=	P[i].GravAccel[j] * dt_gravkick;
				P[i].Vel[j]	+=	dv[j];
			}
		}
	}

	if(All.PM_Ti_endstep == All.Ti_Current){
	
/* the dt_displacement parameter should be set properly! */
		ti_step	=	TIMEBASE;
		while(ti_step > (dt_displacement / All.Timebase_interval))
			ti_step	>>=	1;

		if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep)){
			if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
				ti_step	=	All.PM_Ti_endstep - All.PM_Ti_begstep;
		}

		if(All.Ti_Current == TIMEBASE)
			ti_step	=	0;
		
		tstart	=	(All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
		tend	=	All.PM_Ti_endstep + ti_step / 2;

		dt_gravkick	=	get_gravkick_factor(tstart, tend);

		All.PM_Ti_begstep	=	All.PM_Ti_endstep;
		All.PM_Ti_endstep	=	All.PM_Ti_begstep + ti_step;

		for(i = 0; i < NumPart; ++i){
			for(j = 0; j < 3; ++j){
				P[i].Vel[j]	+=	P[i].GravPM[j] * dt_gravkick;
//				printf("%g|%g\n", P[i].GravPM[j], P[i].GravAccel[j]);
			}
//			printf("\n");
		}
	}
}

int get_timestep(int p){

	double	ax, ay, az, ac;
	double	dt;
	int		ti_step;

	ax	=	fac1 * P[p].GravAccel[0];
	ay	=	fac1 * P[p].GravAccel[1];
	az	=	fac1 * P[p].GravAccel[2];

	ax	+=	fac1 * P[p].GravPM[0];
	ay	+=	fac1 * P[p].GravPM[1];
	az	+=	fac1 * P[p].GravPM[2];
	
	ac	=	sqrt(ax * ax + ay * ay + az * az);
	if(ac == 0)
		ac	=	1.0e-30;

	dt	=	sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable / ac);
	
	dt	*=	hubble_a;

	if(dt >= dt_displacement)
		dt	=	dt_displacement;

	if(dt < All.MinSizeTimestep){
		printf("warning: Timestep will be below the limit of MinSizeTimestep!\n");
		printf("Part ID: %d	dt: %g	ac: %g\n", p, dt, ac);
		dt	=	All.MinSizeTimestep;
	}
	
	ti_step	=	dt / All.Timebase_interval;

	if(!(ti_step > 0 && ti_step < TIMEBASE)){
		printf("Error: a timestep of size zero was assigned on the integer timeline!\n");
		printf("ti_step: %d\n", ti_step);
		printf("timebase interval: %g\n", All.Timebase_interval);
		printf("ax|ay|az: %g|%g|%g\n", ax, ay, az);
		exit(1);
	}

	return	ti_step;
}
