/******************************************************************************
 *                                                                            *
 * This is the header file for traces() version 2.0, which is included into   *
 *   nauty() version 2.5.                                                     *
 *  nauty is Copyright (1984-2013) Brendan McKay.  All rights reserved.      *
 *  Subject to the waivers and disclaimers in nauty.h.                       *
 *  Traces is Copyright Adolfo Piperno, 2008-2013.  All rights reserved.     *
 *                                                                            *
 *   CHANGE HISTORY                                                           *
 *       28-Dec-12 : final changes for version 2.0                            *
 *       20-Jan-13 - add code for ^C catching in Traces                       *
 *****************************************************************************/

#include "gtools.h"
#include "schreier.h" 

#define NAUABORTED   4      /* Traces/nauty is terminated early under program control */
#define NAUKILLED    5      /* Traces/nauty is terminated early by caught signal */

typedef struct TracesOptions {
	boolean getcanon;
	boolean writeautoms;
	boolean cartesian;
	boolean digraph;
	boolean defaultptn;
	int linelength;
	FILE* outfile;
	int strategy;
	int verbosity;
	permnode **generators;
    void (*userautomproc)(int,int*,int);
    int  (*usercanonproc)(graph*,int*,graph*,int,int,int,int);
} TracesOptions;

#define DEFAULTOPTIONS_TRACES(opts) TracesOptions opts \
= { FALSE, FALSE, FALSE, FALSE, TRUE, 0, NULL, 0, 0, NULL, NULL, NULL }

typedef struct TracesStats {
	double grpsize1;
	int grpsize2;
	int numgenerators;
	int numorbits;
	int treedepth;
	int canupdates;
	int errstatus;
	unsigned long numnodes;
	unsigned long interrupted;
	unsigned long peaknodes;
} TracesStats;

#ifdef __cplusplus
extern "C" {
#endif

extern void Traces(sparsegraph*,int*,int*,int*,TracesOptions*,
				   TracesStats*,sparsegraph*);									
extern void refine_tr(sparsegraph*,int*,int*,int*,int*,TracesOptions*);		
extern void traces_freedyn(void);

#ifdef __cplusplus
}
#endif
