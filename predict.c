#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	<gsl/gsl_math.h>

#include	"allvars.h"
#include	"proto.h"

void move_particles(int time0, int time1){
	
	int		i, j;
	double	dt_drift;

	dt_drift	=	get_drift_factor(time0, time1);

	for(i = 0; i < NumPart; ++i){
		for(j = 0; j < 3; ++j){
			P[i].Pos[j]	+=	P[i].Vel[j] * dt_drift;
		}
	}
}

void do_box_wrapping(void){

	int		i, j;
	double	boxsize;

	boxsize	=	All.BoxSize;
	
	for(i = 0; i < NumPart; ++i){
		for(j = 0; j < 3; ++j){
			
			while(P[i].Pos[j] < 0)
				P[i].Pos[j]	+=	boxsize;
			
			while(P[i].Pos[j] >= boxsize)
				P[i].Pos[j]	-=	boxsize;
		}
	}
}
