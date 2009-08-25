#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

void long_range_init(void){

	pm_init_periodic();

}

void long_range_force(void){
	
	int	i;
	for(i = 0; i < NumPart; ++i){
		P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;
	}

	pmforce_periodic();
}
