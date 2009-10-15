#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

void read_ic(char * fname){
	
	FILE *	fp;
	if(!(fp=fopen(fname, "r"))){
		printf("Cannot open initial condition file %s!\n",fname);
		exit(1);
	}

	fread(&header,sizeof(HEADER),1,fp);
	
	NumPart	=	header.NumPart;
	P	=	(PARTICLE *)malloc(sizeof(PARTICLE)*NumPart);
	fread(P,sizeof(PARTICLE),NumPart,fp);

	fclose(fp);

	All.NumPart	=	NumPart;
	All.BoxSize	=	header.BoxSize;
	All.Mass	=	header.Mass;

	printf("Particle number: %d\n", NumPart);
	printf("Boxsize(kpc/h): %f\n", header.BoxSize);
	printf("Mass(10^10/h solar mass): %f\n", header.Mass);
	printf("Redshift: %f\n", header.Redshift);
}

	
