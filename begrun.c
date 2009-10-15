#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<gsl/gsl_rng.h>

#include	"allvars.h"
#include	"proto.h"

void begrun(void){

	read_parameter_file();

	set_units();

	long_range_init();

	init();

	init_drift_table();

#ifdef	LINKLIST
	linklist_init();
#endif

	All.Ti_nextoutput	=	find_next_outputtime(All.Ti_Current);
}

/*	At present I just set parameters in this function instead of from a file.
 */
void read_parameter_file(void){
	
	All.UnitLength_in_m	=	3.085678e19;	/* 1.0 kpc in m */
	All.UnitMass_in_kg	=	1.98892e40;	/* 1.0e10 solar mass in kg */
	All.UnitVelocity_in_m_per_s	=	1e3; /* 1 km/s in m/s */

	All.Omega0		=	0.3;
	All.OmegaLambda	=	0.7;
	All.OmegaBaryon	=	0.04;
	All.HubbleParam	=	0.7;

	All.BoxSize		=	10000.0;

	All.ErrTolTheta	=	0.5;

	All.ErrTolIntAccuracy	=	0.025;

	All.MaxSizeTimestep	=	0.03;
	All.MinSizeTimestep	=	0.0;

	All.SofteningHalo	=	18.0;
	All.SofteningHaloMaxPhys	=	1.8;

/* The NumPart, Mass will be set according to the initial condition file, 
 */
	All.TimeBegin	=	0.0909; /* z=10 */
	All.TimeMax		=	1.0;

	strcpy(All.InitCondFile, "/home/liulei/program/N-body/ic/ic_32.cuco");
	strcpy(All.OutputDir, "output_32");
	strcpy(All.SnapshotFileBase, "cuco_32");
	strcpy(All.OutputListFilename, "output_32/list_32.txt");

	read_outputlist(All.OutputListFilename);
}

void set_units(void){
	
	All.UnitTime_in_s	=	All.UnitLength_in_m / All.UnitVelocity_in_m_per_s;
	All.UnitTime_in_Megayears	=	All.UnitTime_in_s / SEC_PER_MEGAYEAR;
	All.G	=	GRAVITY / pow(All.UnitLength_in_m, 3) * All.UnitMass_in_kg * pow(All.UnitTime_in_s, 2);
	All.Hubble	=	HUBBLE * All.UnitTime_in_s;

	printf("Hubble (internal units) = %g\n", All.Hubble);
	printf("G (internal units) = %g\n", All.G);
	printf("UnitTime_in_kg = %g\n", All.UnitMass_in_kg);
	printf("UnitTime_in_s = %g\n", All.UnitTime_in_s);
	printf("UnitVelocity_in_m_per_s = %g\n",All.UnitVelocity_in_m_per_s);

}

int read_outputlist(char * fname){
	
	FILE *	fp;
	if(!(fp = fopen(fname, "r"))){
		printf("cannot read output list file %s!\n", fname);
		return(1);
	}

	All.OutputListLength	=	0;

	do{
		if(fscanf(fp, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
			All.OutputListLength++;
		else
			break;
	}while(All.OutputListLength < MAXLEN_OUTPUTLIST);

	fclose(fp);

	printf("find %d output times in output list file %s\n", 
			All.OutputListLength, fname);
	return(0);
}

