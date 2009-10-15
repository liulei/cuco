/*	This file declares all global variables. 
 */

#ifndef	ALLVARS_H
#define	ALLVARS_H

#include	<stdio.h>

#define	PMGRID		128

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
	double	BoxSize;
	double	Mass;

/*	Some parameter in system units
 */
	double	G;
	double	UnitTime_in_s;
	double	UnitMass_in_kg;
	double	UnitVelocity_in_m_per_s;
	double	UnitLength_in_m;
	double	UnitTime_in_Megayears;

	double	Hubble;
	double	Omega0;
	double	OmegaLambda;
	double	OmegaBaryon;
	double	HubbleParam;

/*	Time parameter
 */
	double	Time;
	double	TimeBegin;
	double	TimeStep;
	double	TimeMax;

	double	Timebase_interval;
	int		Ti_Current;
	int		Ti_nextoutput;

	int		PM_Ti_endstep;
	int		PM_Ti_begstep;

	double	Asmth;	/* Long-range/short-range split */
	double	Rcut;	/* Maximum radius for which short range force is evaluated */
	
	double	ErrTolTheta;

	double	ErrTolIntAccuracy;

	double	MaxSizeTimestep;
	double	MinSizeTimestep;

	double	SofteningHalo;
	double	SofteningHaloMaxPhys;
	double	SofteningTable;
	double	ForceSoftening;

	int		NumCurrentTiStep;

	int		SnapshotFileCount;

	double	OutputListTimes[MAXLEN_OUTPUTLIST];
	int		OutputListLength;

	char	InitCondFile[MAXLEN_FILENAME];
	char	OutputDir[MAXLEN_FILENAME];
	char	SnapshotFileBase[MAXLEN_FILENAME];
	char	OutputListFilename[MAXLEN_FILENAME];

}ALL;

typedef struct tagNODE{
	FLOAT	len;
	FLOAT	center[3];

	union{
		int	suns[8];
		struct{
			FLOAT	s[3];
			FLOAT	mass;
			int		bitflags;
			int		sibling;
			int		nextnode;
			int		father;
		}d;
	}u;
} NODE, NODE_BASE;

typedef struct tagEXTNODE{
	FLOAT	vs[3];
}EXTNODE, EXTNODE_BASE;

extern	int			NumPart;
extern	PARTICLE *	P;
extern	HEADER		header;
extern	ALL			All;
extern	double		DriftTable[DRIFT_TABLE_LENGTH];
extern	double		GravKickTable[DRIFT_TABLE_LENGTH];
extern	int			TreeReconstructFlag;
extern	NODE		* Nodes;
extern	NODE_BASE	* Nodes_base;
extern	int			MaxNodes;
extern	int			* Nextnode;
extern	int			* Father;
extern	EXTNODE		* Extnodes;
extern	EXTNODE_BASE	* Extnodes_base;
extern	int			* Father;

#endif
