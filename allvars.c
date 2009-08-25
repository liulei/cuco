/*	Provides instances of all global variables;
 */

#include	<stdio.h>
#include	"allvars.h"

unsigned int			NumPart;
PARTICLE *	P;
HEADER		header;
ALL			All;
float		DriftTable[DRIFT_TABLE_LENGTH];
float		GravKickTable[DRIFT_TABLE_LENGTH];

unsigned int	numParticles;
unsigned int	numGridCells;

float	*hPos;
float	*hGravAccel;

float	*dOldPos;
float	*dSortedPos;
float	*dGravAccel;

unsigned int	*dGridParticleHash;
unsigned int	*dGridParticleIndex;
unsigned int	*dCellStart;
unsigned int	*dCellEnd;

//RadixSort	*sorter;
