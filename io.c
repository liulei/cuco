#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>

#include	"allvars.h"
#include	"proto.h"

void savepositions(int num){
	
	FILE *	fp;
	
	char	buffer[250];
	sprintf(buffer, "%s/%s_%d.cuco", 
			All.OutputDir, All.SnapshotFileBase, num);

	if(!(fp = fopen(buffer, "w"))){
		printf("Cannot open snapshot file %d!\n", num);
		exit(1);
	}

	HEADER	header1;

	header1.NumPart	=	header.NumPart;
	header1.Mass	=	header.Mass;
	header1.Time	=	All.Time;
	header1.Redshift=	1 / All.Time -1;
	header1.BoxSize	=	header.BoxSize;
	header1.Omega0	=	header.Omega0;
	header1.OmegaLambda	=	header.OmegaLambda;
	header1.HubbleParam	=	header.HubbleParam;

	fwrite(&header1, sizeof(HEADER), 1, fp);
	fwrite(P, sizeof(PARTICLE), NumPart, fp);

	fclose(fp);
	
	printf("write output file %d: a = %g z = %g\n", 
			num, All.Time, header1.Redshift);
}
