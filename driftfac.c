#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#include	<gsl/gsl_math.h>
#include	<gsl/gsl_integration.h>

#include	"allvars.h"
#include	"proto.h"

static	double	logTimeBegin;
static	double	logTimeMax;

void init_drift_table(void){

#define	WORKSIZE	100000

	int		i;
	double	result, abserr;
	gsl_function	F;
	gsl_integration_workspace *	workspace;

	logTimeBegin	=	log(All.TimeBegin);
	logTimeMax		=	log(All.TimeMax);

	workspace	=	gsl_integration_workspace_alloc(WORKSIZE);

	for(i = 0; i < DRIFT_TABLE_LENGTH; ++i){

		F.function	=	&drift_integ;
		gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,
			1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
		DriftTable[i]	=	result;

		F.function	=	&gravkick_integ;
		gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,
			1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
		GravKickTable[i]	=	result;
	}

	gsl_integration_workspace_free(workspace);
}

double drift_integ(double a, void * param){
	
	double	h;

	h	=	All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a)
			+All.OmegaLambda;
	h	=	All.Hubble * sqrt(h);

	return(1 / (h * a * a * a));
}

double gravkick_integ(double a, void * param){
	
	double	h;

	h	=	All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a)
			+ All.OmegaLambda;
	h	=	All.Hubble * sqrt(h);

	return(1 / (h * a * a));
}

double get_drift_factor(int time0, int time1){
	
	double	a1, a2, df1, df2, u1, u2;
	int		i1, i2;

	a1	=	logTimeBegin + time0 * All.Timebase_interval;
	a2	=	logTimeBegin + time1 * All.Timebase_interval;

	u1	=	(a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
	i1	=	(int) u1;
	if(i1 >= DRIFT_TABLE_LENGTH)
		i1	=	DRIFT_TABLE_LENGTH - 1;
	
	if(i1 < 1)
		df1	=	u1 * DriftTable[0];
	else
		df1	=	DriftTable[i1 - 1] + (DriftTable[i1] - DriftTable[i1 - 1]) * (u1 - i1);

	u2	=	(a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
	i2	=	(int) u2;
	if(i2 >= DRIFT_TABLE_LENGTH)
		i2	=	DRIFT_TABLE_LENGTH - 1;
	
	if(i2 < 1)
		df2	=	u2 * DriftTable[0];
	else
		df2	=	DriftTable[i2 - 1] + (DriftTable[i2] - DriftTable[i2 - 1]) * (u2 - i2);

	return df2 - df1;
}

double get_gravkick_factor(int time0, int time1){
	
	double	a1, a2, df1, df2, u1, u2;
	int		i1, i2;

	a1	=	logTimeBegin + time0 * All.Timebase_interval;
	a2	=	logTimeBegin + time1 * All.Timebase_interval;

	u1	=	(a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
	i1	=	(int) u1;
	if(i1 >= DRIFT_TABLE_LENGTH)
		i1	=	DRIFT_TABLE_LENGTH - 1;
	
	if(i1 < 1)
		df1	=	u1 * GravKickTable[0];
	else
		df1	=	GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);

	u2	=	(a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
	i2	=	(int) u2;
	if(i2 >= DRIFT_TABLE_LENGTH)
		i2	=	DRIFT_TABLE_LENGTH - 1;
	
	if(i2 < 1)
		df2	=	u2 * GravKickTable[0];
	else
		df2	=	GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

	return df2 - df1;
}
