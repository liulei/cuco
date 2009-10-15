/*	Provides instances of all global variables;
 */

#include	<stdio.h>
#include	"allvars.h"

int			NumPart;
PARTICLE *	P;
HEADER		header;
ALL			All;
double		DriftTable[DRIFT_TABLE_LENGTH];
double		GravKickTable[DRIFT_TABLE_LENGTH];
int			TreeReconstructFlag;
NODE		* Nodes;
NODE_BASE	* Nodes_base;
int			MaxNodes;
int			* Nextnode;
int			* Father;
EXTNODE		* Extnodes;
EXTNODE_BASE	* Extnodes_base;
int			* Father;
