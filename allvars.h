/*	This file declares all global variables. 
 */

#ifndef	ALLVARS_H
#define	ALLVARS_H

#include	<stdio.h>

#define	PMGRID			64

#define	TIMEBASE		(1<<28)

/*	Some physical constants in mks units
 */
#define	GRAVITY			6.673e-11	/* in mks units */
#define	SOLAR_MASS		1.98892e30	/* in kg */
#define	C				2.99792458e8	/* in m/s */
#define	HUBBLE			3.2407789e-18	/* in h/sec */

#define	METER_PER_MPC		3.08568025e22	
#define	SEC_PER_MEGAYEAR	3.155e13
#define	SEC_PER_YEAR		3.155e7

#define	ASMTH	1.25
#define	RCUT	4.5

#define	NTAB	1000

#define	DRIFT_TABLE_LENGTH	1000
#define	MAXLEN_OUTPUTLIST	500
#define	MAXLEN_FILENAME		250

typedef	float	FLOAT;

typedef	struct	particle_data{
	
	FLOAT	Pos[3];
	FLOAT	Mass;
	FLOAT	Vel[3];
	FLOAT	GravAccel[3];
	FLOAT	GravPM[3];

	int		Ti_endstep;
	int		Ti_begstep;
}PARTICLE;

typedef struct	io_header{

	int		NumPart;
	float	Mass;
	float	Time;
	float	Redshift;
	float	BoxSize;
	float	Omega0;
	float	OmegaLambda;
	float	HubbleParam;
}HEADER;

typedef	struct	global_data{
	int		NumPart;
	float	BoxSize;
	float	Mass;

/*	Some parameter in system units
 */
	double	G;
	double	UnitTime_in_s;
	double	UnitMass_in_kg;
	double	UnitVelocity_in_m_per_s;
	double	UnitLength_in_m;
	double	UnitTime_in_Megayears;

	float	Hubble;
	float	Omega0;
	float	OmegaLambda;
	float	OmegaBaryon;
	float	HubbleParam;

/*	Time parameter
 */
	float	Time;
	float	TimeBegin;
	float	TimeStep;
	float	TimeMax;

	double	Timebase_interval;
	int		Ti_Current;
	int		Ti_nextoutput;

	int		PM_Ti_endstep;
	int		PM_Ti_begstep;

	float	Asmth;	/* Long-range/short-range split */
	float	Rcut;	/* Maximum radius for which short range force is evaluated */

	float	ErrTolIntAccuracy;

	float	MaxSizeTimestep;
	float	MinSizeTimestep;

	float	SofteningHalo;
	float	SofteningHaloMaxPhys;
	float	SofteningTable;
	float	ForceSoftening;

	int		NumCurrentTiStep;

	int		SnapshotFileCount;

	float	OutputListTimes[MAXLEN_OUTPUTLIST];
	int		OutputListLength;

	char	InitCondFile[MAXLEN_FILENAME];
	char	OutputDir[MAXLEN_FILENAME];
	char	SnapshotFileBase[MAXLEN_FILENAME];
	char	OutputListFilename[MAXLEN_FILENAME];

}ALL;

extern	unsigned int			NumPart;
extern	PARTICLE *	P;
extern	HEADER		header;
extern	ALL			All;
extern	float		DriftTable[DRIFT_TABLE_LENGTH];
extern	float		GravKickTable[DRIFT_TABLE_LENGTH];

extern	unsigned int	numParticles;
extern	unsigned int	numGridCells;

/*	CPU data	*/
extern	float	*hPos;
extern	float	*hGravAccel;

/*	GPU data	*/
extern	float	*dOldPos;
extern	float	*dSortedPos;
extern	float	*dGravAccel;

extern	unsigned int	*dGridParticleHash;
extern	unsigned int	*dGridParticleIndex;
extern	unsigned int	*dCellStart;
extern	unsigned int	*dCellEnd;

//extern	RadixSort		*sorter;

#endif
