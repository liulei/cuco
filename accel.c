#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

void compute_accelerations(void){

	do_box_wrapping();

/* If long range kick needed, computer the long range force
 */
	if(All.PM_Ti_endstep == All.Ti_Current){
		long_range_force();
	}

/* Calculate the short range force by direct summation with linklist 
 * method to enhance calculation.
 */
#ifdef	TREE
	gravity_tree();
#endif

#ifdef	LINKLIST
	gravity_linklist();
#endif
}
