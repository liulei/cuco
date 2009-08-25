#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

void run(void){
	
	do{
		find_next_sync_point_and_drift();

		every_timestep_stuff();
		
		compute_accelerations();

		advance_and_find_timesteps();

		All.NumCurrentTiStep++;

	}while(All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

	savepositions(All.SnapshotFileCount++);

}

void find_next_sync_point_and_drift(void){

	int		i;
	double	timeold;
	double	min, min_glob;

	timeold	=	All.Time;

	min		=	P[0].Ti_endstep;
	for(i = 1; i < NumPart; ++i){
		if(min > P[i].Ti_endstep){
			min	=	P[i].Ti_endstep;
		}
	}

	min_glob	=	min;

	if(min_glob >= All.PM_Ti_endstep){
		min_glob	=	All.PM_Ti_endstep;
	}

	while(min_glob >= All.Ti_nextoutput && All.Ti_nextoutput >= 0){
		
		move_particles(All.Ti_Current, All.Ti_nextoutput);

		All.Ti_Current	=	All.Ti_nextoutput;

		All.Time	=	All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
		savepositions(All.SnapshotFileCount++);

		All.Ti_nextoutput	=	find_next_outputtime(All.Ti_nextoutput + 1);
	}

	move_particles(All.Ti_Current, min_glob);

	All.Ti_Current	=	min_glob;

	All.Time	=	All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);

	All.TimeStep	=	All.Time - timeold;

}

/*	This function ensures that next output time in the proper time range and
 *	must be the one after current output time. Also it allowed a output time 
 *	list with shuffled order.
 */
int find_next_outputtime(int ti_curr){

	int	i, ti, ti_next;
	double	time;
	
	ti_next	=	-1;

	for(i = 0; i < All.OutputListLength; ++i){
		
		time	=	All.OutputListTimes[i];

		if(time >= All.TimeBegin && time <= All.TimeMax){

			ti	=	log(time / All.TimeBegin) / All.Timebase_interval;
			
			if(ti > ti_curr){
				if(ti_next == -1)
					ti_next	=	ti;
				if(ti_next > ti)
					ti_next	=	ti;
			}
		}
	}

	if(ti_next == -1){
		ti_next	=	2 * TIMEBASE;
	}

	return(ti_next);
}

void every_timestep_stuff(void){
	
	double	z;
	
	z	=	1.0 / (All.Time) - 1;

	printf("\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time)- log(All.Time - All.TimeStep));
}
