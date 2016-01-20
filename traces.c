/******************************************************************************
 *                                                                            *
 * This is the main file for traces() version 2.0, which is included into     *
 *   nauty() version 2.5.                                                     *
 *                                                                            *
 *   nauty is Copyright (1984-2013) Brendan McKay.  All rights reserved.      *
 *   Subject to the waivers and disclaimers in nauty.h.                       *
 *   Traces is Copyright Adolfo Piperno, 2008-2013.  All rights reserved.     *
 *                                                                            *
 *   CHANGE HISTORY                                                           *
 *       28-Dec-12 : final changes for version 2.0                            *
 *       20-Jan-13 : add code for ^C catching in Traces                       *
 *       29-Mar-13 : bug correction in automorphism mode                      *
 *       02-Apr-13 : add preprocessing                                        *
 *****************************************************************************/

#include "traces.h"

#define SORT_OF_SORT 2
#define SORT_NAME sort2ints
#define SORT_TYPE1 int
#define SORT_TYPE2 int
#include "sorttemplates.c"

#ifdef NAUTY_IN_MAGMA
#include "cleanup.e"
#endif

#define NAUTY_ABORTED (-11)
#define NAUTY_KILLED (-12)

typedef struct Candidate {
    boolean sortedlab;
    int *invlab;
    int *lab;
    int code;
    int do_it;
    int indnum;
    int name;
    int vertex;
    struct Candidate *next;
    struct searchtrie *stnode;
    unsigned int firstsingcode;
    unsigned int pathsingcode;
    unsigned int singcode;
} Candidate;

typedef struct Partition {
    int *cls;
    int *inv;
    int active;
    int cells;
    int code;
} Partition;

typedef struct trielist {
    struct searchtrie *triearray;
    struct trielist *prev;
    struct trielist *next;
} trielist;

typedef struct TracesVars {
    char digstring[25];
	double autchk;
	double expaths;
	double schreier1;
	double schreier2;
	double schreier3;
    int build_autom;
    int *currorbit;
    int *orbits;
    int answ;
    int brkstpcount;
    int compstage;
    int cand_level;
    int canlist;
    int digits;
    int expathlength;
    int firstpathlength;
    int fromlevel;
    int indivend;
    int indivstart;
    int indiv_vtx;
    int lastcell;
    int lastlev;
    int levelfromCS0;
    int linelgth;
    int mark;
    int treemark;
    int autmark;
    int markcell1;
    int markcell2;
    int maxdeg;
    int maxtreelevel;
    int maxspineorblevel;
    int name;
    struct searchtrie *gotonode;
    struct searchtrie *newgotonode;
    struct searchtrie *newst_stage1;
    int newindex;
    int nextlevel;
    int nfix;
    int finalnumcells;
    int permInd;
    int preprocessed;
    int samepref;
    int specialgens;
    int stackmark;
    int steps;
    int strategy;
    trielist *strielist;
    int strienext;
    int tcell_sz;
    int tcell;
    int tcellevel;
    int tcellexpath_sz;
    int tcellexpath;
    int tolevel_tl;
    int tolevel;
    int treedepth;
    int trienext;
    int triepos;
    TracesOptions *options;
    TracesStats *stats;
    unsigned int singlongcode;
    sparsegraph *graph;
    sparsegraph *cangraph;
    sparsegraph *input_graph;
    int conta0;
    int conta1;
    int conta2;
    int conta3;
    int conta4;
    int conta5;
    int conta6;
    int conta7;
    int contatc;
    //    double auxtime1;
    //    double auxtime2;
    //    double auxtime3;
    //    double auxtime4;
    //    double auxtime5;
} TracesVars;

typedef struct TracesInfo {
    boolean autofound;
    boolean deg_one;
    boolean exitfromref;
    boolean identitygroup;
    boolean minimalinorbits;
    boolean newtc;
    boolean thegraphisparse;
    boolean thegrouphaschanged;
    boolean thereisnextlevel;
    boolean useTempOrbits1;
    boolean useTempOrbits2;
} TracesInfo;

typedef struct TracesSpine {
    boolean thetracexists;
    boolean thexpathexists;
    Candidate *listend;
    Candidate *liststart;
    int ccend;
    int ccstart;
    int idx;
    int listcounter;
    int stpend;
    int stpstart;
    int tgtcell;
    int tgtend;
    int tgtfrom;
    int tgtpos;
    int tgtsize;
    int trcend;
    int trcstart;
    int singend;
    int singstart;
    int updates;
    unsigned long keptcounter;
    unsigned long levelcounter;
    Partition *part;
    unsigned int singcode;
} TracesSpine;

typedef struct trie {
	int value;
	struct trie *first_child;
	struct trie *next_sibling;
} trie;

typedef struct searchtrie {
	int index;
	int name;
	int vtx;
    int level;
	struct searchtrie *father;
	struct searchtrie *first_child;
	struct searchtrie *last_child;
	struct searchtrie *next_sibling;
	struct searchtrie *goes_to;
} searchtrie;

typedef struct pair {
	int arg;
	int val;
} pair;

typedef struct grph_strct {
	int* e;
	int d;
	boolean one;
} grph_strct;

static boolean traces_degree_refine(sparsegraph*, Candidate*, int, Partition*, 
						  struct TracesVars*, struct TracesInfo*, int);
static int  traces_refine(Candidate*, int, Partition*, 
						  struct TracesVars*, struct TracesInfo*, int, boolean);
static void traces_refine_notrace(Candidate*, int, Partition*, 
								  struct TracesVars*, struct TracesInfo*);
static void traces_refine_maketrie(Candidate*, int, Partition*, 
								   struct TracesVars*, struct TracesInfo*);
static int  traces_refine_comptrie(Candidate*, int, Partition*, 
								   struct TracesVars*, struct TracesInfo*);
static int  traces_refine_sametrace(Candidate*, int, Partition*, 
									struct TracesVars*, struct TracesInfo*);
static int  traces_refine_refine(sparsegraph*, Candidate*, int, Partition*, 
								 struct TracesVars*, struct TracesInfo*);
static void quickSort(int*, int);
static struct Candidate* NewCandidate(int, Candidate**, int);
static int FreeList(Candidate*, int);
static int FixBase(int*, struct TracesVars*, Candidate*, int, int);
static void factorial(double*, int*, int);
static void factorial2(double*, int*, int);
static int CheckForAutomorphisms(Candidate*, Candidate*, struct TracesVars*, struct TracesInfo*, int, int, Partition*);
static int CheckForSingAutomorphisms(Candidate*, Partition*, Candidate*, struct TracesVars*, struct TracesInfo*, int, int);
static int CheckForMatching(Candidate*, Candidate*, Partition*, struct TracesVars*, struct TracesInfo*, int, int);
static void Individualize(Partition*, Candidate*, int, int, int, int);
static boolean TreeFyTwo(int, Candidate*, Candidate*, Partition*, int, struct TracesVars*, struct TracesInfo*);
static void ExperimentalStep(Partition*, Candidate*, struct TracesVars*, struct TracesInfo*, int, int);
static boolean TargetCell(Candidate*, Partition*, int, struct TracesVars*, int);
static boolean TargetCellFirstPath(Candidate*, Partition*, struct TracesVars*);
static int TargetCellExpPath(Candidate*, Partition*, struct TracesVars*);
static boolean SelectNextLevel(int, struct TracesVars*);
static void CopyCand(Candidate*, Candidate*, int, int*, int*);
static struct trie* trie_new(int, struct TracesVars*);
static struct trie* trie_make(trie*, int, int, struct TracesVars*);
static struct trie* trie_comp(trie*, int);
static void RemoveFromLevel(int, int, int, boolean);
static int CompStage0(Partition*, Partition*, Candidate*, Candidate*, int, int, struct TracesVars*, struct TracesInfo*);
static int CompStage1(Partition*, Partition*, Candidate*, Candidate*, int, int, struct TracesVars*, struct TracesInfo*);
static int CompStage2(Partition*, Partition*, Candidate*, Candidate*, int, int, struct TracesVars*, struct TracesInfo*);
static void grouporderplus(sparsegraph*, Candidate*, Partition*, permnode**, double*, int*, int, struct TracesVars*, struct TracesInfo*);
static boolean Prefix(Candidate*, Candidate*, int);
static boolean findperm(permnode*, int*, int);
static int spinelementorbsize(int*, int*, int, int);
static trielist* searchtrie_new(int, struct TracesVars*);
static searchtrie* searchtrie_make(Candidate*, Candidate*, int, struct TracesVars*);
static boolean lookup(searchtrie*);
static int* findcurrorbits(schreier*, int);
static int Preprocess(sparsegraph*, permnode**, Candidate*, int, Partition*, struct TracesVars*);
static void MakeTree(int, int, sparsegraph*, int, struct TracesVars*, boolean);
static void MakeCanTree(int, sparsegraph*, int, Candidate*, Partition*, struct TracesVars*);
static int max(int, int);
static int min(int, int);
static void orbjoin_sp_perm(int*, int*, int*, int, int*);
static void orbjoin_sp_pair(int*, int*, int, int, int, int*);
static boolean isautom_sg_pair(graph*, int*, boolean, int, int, struct TracesVars*);
static void SetAutom(int, int, struct TracesVars*);
static void ResetAutom(int, int, struct TracesVars*);
static void PrintVect(int*, int, int, int);
static void putgraphplus_sg(FILE*, sparsegraph*, sparsegraph*, int);
static boolean VerifyId(int *p, int n);
static void PrintPartition(int*, int*, int, int, int);
static void Place(int, Candidate*, Partition*);
static int NonSingDeg(int, Candidate*, Partition*);
static int NonSingDegPlus1(Candidate*, Partition*, int, TracesVars*);
static void NonSingDegPlus2(Candidate*, Partition*, int, TracesVars*);
static void Edge_Delete(int, int, Candidate*, TracesVars*);
static boolean VerifyPart(int*, int, int, int*);
static int VerifyPerm(int*, int,int);
static boolean VerifyCand(Candidate*, int, int);
static int FirstNeighbour(int, Candidate*, Partition*, int*, int, int*, int);
static int NextNeighbour(int, Candidate*, Partition*, int*, int, int*, int);
static sparsegraph* copy_sg_structure(sparsegraph*, sparsegraph*);


static const unsigned int fuzz1[] = {037541, 061532, 005257, 026416};
static const unsigned int fuzz2[] = {006532, 070236, 035523, 062437};

#define FUZZ1(x) ((x) ^ fuzz1[(x)&3])
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3])

#define STATS_INIT stats_arg->grpsize1 = 1; \
stats_arg->grpsize2 = 0; \
stats_arg->numorbits = n; \
stats_arg->treedepth= 0; \
stats_arg->numgenerators = 0; \
stats_arg->numnodes = 1; \
stats_arg->errstatus = 0; \
stats_arg->interrupted = 0; \
stats_arg->canupdates = 0; \
stats_arg->peaknodes = 0; \

#define TIME_INIT tv->autchk = 0; \
tv->expaths = 0; \
tv->schreier1 = 0; \
tv->schreier2 = 0; \
tv->schreier3 = 0;

#define MASHCOMM(l, i) ((l) + FUZZ1(i))
#define MASHNONCOMM(l, i) (FUZZ2(l) + (i))
#define MASH(l, i) ((((l) ^ 065435) + (i)) & 077777)
#define MASH1(l, i) ((l + (i*i)) & 077777)
#define CLEANUP(l) ((int)((l) % 0x7FFF))
#define SS(n, sing, plur)  (n), ((n) == 1?(sing):(plur))

#define SETMARK(Arr, Mrk) if (Mrk > (NAUTY_INFINITY-2)) { memset(Arr, 0, n*sizeof(int)); Mrk = 0; } Mrk++;

#define COPYNODE(W, V) { \
memcpy(W->lab, V->lab, n*sizeof(int)); \
memcpy(W->invlab, V->invlab, n*sizeof(int)); \
W->code = V->code; \
W->singcode = V->singcode; \
W->do_it = V->do_it; }

#define NEXTLINE fprintf(outfile, "\n");

#define PRINTCHAR(c) fprintf(outfile, "%s", c);

#define PRINTCAND(V, Lev) PRINTCHAR(" ") for (tmp=1; tmp<=Lev; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);}

#define PRINTCANDBIG(V, Lev) { PRINTCHAR(" ") \
for (tmp=1; tmp<=5; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);} \
fprintf(outfile, "... "); \
for (tmp=Lev-4; tmp<=Lev; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);} }

#define LINE(K, c) { PRINTCHAR(c) for (tmp=1; tmp<=K; tmp++) {fprintf(outfile, c);} }

#define TRACE_CHECK(Tr, Ind, Arg, End) { TracePos = Tr+Ind; \
if (newtrace) { \
*TracePos = Arg; \
} \
else { \
if (Ind < *End) { \
if (*TracePos != Arg) { \
if (*TracePos > Arg) { \
return FALSE; \
} \
else { \
*TracePos = Arg; \
newtrace = TRUE; \
} \
} \
} \
else { \
*TracePos = Arg; \
newtrace = TRUE; \
} \
} \
Ind++; }

#define SAMETRACE_CHECK(Tr, Ind, Arg, End) { TracePos = Tr+Ind; \
if (Ind < *End) { \
if (*TracePos != Arg) { \
return FALSE; \
} \
} \
else { \
return FALSE; \
} \
Ind++; }

#define NEWPART(P) { P = malloc(sizeof(*(P))); \
if (P == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
P->cls = malloc(n*sizeof(int)); \
if (P->cls == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
P->inv = malloc(n*sizeof(int)); \
if (P->inv == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
P->code = -1; \
P->cells = 0; }

#define NEWPARTSPINE(Lev) { if (Lev > 3) { \
Spine[Lev].part = malloc(sizeof(*(Spine[Lev].part))); \
if (Spine[Lev].part == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
Spine[Lev].part->cls = Spine[Lev-3].part->cls; \
Spine[Lev].part->inv = Spine[Lev-3].part->inv; \
Spine[Lev-3].part->cls = Spine[Lev-3].part->inv = NULL; \
Spine[Lev].part->code = -1; \
Spine[Lev].part->cells = 0; \
} \
else { \
NEWPART(Spine[Lev].part) \
} }

#define FREEPART(Part) { if (Part) { \
if (Part->cls) free(Part->cls); \
if (Part->inv) free(Part->inv); \
free(Part); } \
}

#define FREECAND(Cand) { if (Cand) { \
if (Cand->lab) free(Cand->lab); \
if (Cand->invlab) free(Cand->invlab); \
free(Cand); \
} }

#define COPYPART(P, Q) { memcpy(P->cls, Q->cls, n*sizeof(int)); \
memcpy(P->inv, Q->inv, n*sizeof(int)); \
P->cells = Q->cells; \
P->code = Q->code; } \

#define ADDTONEXTLEVEL { if (SpineTL->listend) { \
(SpineTL->listend)->next = NewCandidate(n, &GarbList, TRUE); \
if ((tv->compstage < 2) && (SpineTL->listcounter <= (NAUTY_INFINITY-2))) SpineTL->listcounter++; \
SpineTL->listend = (SpineTL->listend)->next; \
CopyCand(SpineTL->listend, NextCand, n, NULL, NULL); \
} \
else { \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
if (tv->compstage < 2) SpineTL->listcounter = 1; \
SpineTL->listend = SpineTL->liststart; \
CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL); \
} }

#define ORBITSIZES { memset(OrbSize, 0, n*sizeof(int)); \
for (i=0; i<n; i++) { \
OrbSize[tv->orbits[i]]++; \
} }

#define CURRORBITSIZES { memset(CurrOrbSize, 0, n*sizeof(int)); \
for (i=SpineTL->tgtcell; i<SpineTL->tgtend; i++) { \
CurrOrbSize[tv->currorbit[CurrCand->lab[i]]]++; \
} }

#define EXITFROMSTAGE0REFINE { PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "-=="); \
CurrCand->indnum--; \
RemoveFromLevel(tv->tolevel, tv->treedepth, tv->strategy, FALSE); \
tv->compstage = 1; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->fromlevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = TRUE; \
return 0; }

#define EXITFROMSTAGE0EXPATH2 { PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "=-="); \
tv->compstage = 1; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->tolevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = FALSE; \
return 0; }

#define EXITFROMSTAGE0EXPATH1 { PRINT_RETURN PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "==-"); \
if (SpineTL->liststart) { \
AuxCand = SpineTL->liststart; \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL); \
SpineTL->liststart->next = AuxCand; \
} \
else { \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
SpineTL->listend = SpineTL->liststart; \
SpineTL->liststart->next = NULL; \
CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL); \
} \
tv->compstage = 1; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->tolevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = FALSE; \
return 0; }

#define UPDATE_LINELGTH { if (tv->options->verbosity >= 2) { \
if (tv->tolevel < 12) { \
tv->linelgth = (tv->digits+1)*tv->tolevel+16; \
} \
else { \
tv->linelgth = (tv->digits+1)*10+20; \
} \
} }

#define PRINT_LINE { if ((tv->options->verbosity >= 1) && (tv->strategy == 0)) { \
if (!ti->newtc) \
{ \
if (tv->options->verbosity >= 2) { LINE(tv->linelgth, "-"); \
NEXTLINE} \
} \
ti->newtc = FALSE; \
} }

#define PRINT_LINE_PLUS(Lev) { if ((tv->options->verbosity >= 1) && (tv->strategy == 0)) { \
if (!ti->newtc) \
{ \
if (tv->options->verbosity >= 2) { LINE(min(tv->linelgth, 80), "-"); \
fprintf(outfile, " ");} \
fprintf(outfile, "level %d:  %d cell%s; %d singleton%s; target cell: %d; %d orbit%s; %lu node%s (%lu kept); %d update%s;", \
Lev, SS(Spine[Lev].part->cells, "", "s"), SS(Spine[Lev].singend, "", "s"), Spine[Lev].tgtsize, SS(tv->stats->numorbits, "", "s"), \
SS(Spine[Lev].levelcounter, "", "s"), Spine[Lev].keptcounter, SS(Spine[Lev].updates, "", "s")); \
NEXTLINE \
} \
ti->newtc = FALSE; \
} }

#define PRINT_CANDIDATE(Cand, Lev) { \
for (tmp = Cand->name, cu = 0; tmp > 0; tmp /= 10, ++cu) {} \
for (tmp = Lev, cu1 = 0; tmp > 0; tmp /= 10, ++cu1) {} \
cu = 14-cu-cu1; \
LINE(cu, "-") \
fprintf(outfile, " %d, %d) ", Lev % 10000, Cand->name % 10000000); \
if (Lev < 12) { \
PRINTCAND(Cand, Lev) \
} \
else { \
PRINTCANDBIG(Cand, Lev) \
} \
PRINTCHAR("| ")	\
fprintf(outfile, "{%x, %x} ", Cand->code, Cand->singcode); \
}

#define PRINT_EXPPATHSTEP(Cand, Boo) { if (tv->options->verbosity >= 2) { \
if ((tv->tolevel_tl-tv->tolevel < 6) || (NextPart->cells == tv->finalnumcells)) { \
fprintf(outfile, "%d ", tv->indiv_vtx+labelorg); \
if (tv->options->verbosity >= 2) {if (Boo) fprintf(outfile, "{%x} ", Cand->code); else fprintf(outfile, "{interr.(%d)} ", NextPart->cells);} \
else { if (!Boo) fprintf(outfile, "{interr.(%d)} ", NextPart->cells);} \
if (NextPart->cells == tv->finalnumcells) { \
fprintf(outfile, "(%d) ", tv->tolevel_tl); \
} \
} \
else { \
if (tv->tolevel_tl-tv->tolevel == 6) { \
fprintf(outfile, "... "); \
} \
} \
} }

#define PRINT_RETURN { if (tv->options->verbosity >= 2) { \
fprintf(outfile, "\n"); \
} }

#define PRINT_FROM_VERB(Verb) { if (tv->options->verbosity >= Verb) { \
fprintf(outfile, "FROM: "); \
PRINTCAND(CurrCand, tv->fromlevel) \
fprintf(outfile, " do_it: %d, indnum: %d, stnode->index: %d ", CurrCand->do_it, CurrCand->indnum, CurrCand->stnode->index); \
PRINT_RETURN \
} }

#define PRINT_NOTMIN_VERB(Verb) { if (tv->options->verbosity >= Verb) { \
fprintf(outfile, " is NOT minimal in orbits (1, %d) [%d]; ", gom_level, CurrCand->lab[Spine[gom_level+1].tgtpos]+labelorg); \
fprintf(outfile, "at lev %d, orb[%d] = %d.\n", gom_level+1, CurrCand->lab[Spine[gom_level+1].tgtpos]+labelorg, tv->currorbit[CurrCand->lab[Spine[gom_level+1].tgtpos]]+labelorg); } }

#define PRINT_SKIPPED_VERB(Verb) { if (tv->options->verbosity >= Verb) \
fprintf(outfile, " skipped (0) (orbit[%d] = %d)\n", \
NextCand->vertex+labelorg, tv->currorbit[NextCand->vertex]+labelorg); }

#define PRINT_REFINE_VERB(Verb) { if (tv->options->verbosity >= Verb) \
fprintf(outfile, " REFINE (orbit[%d] = %d)\n", NextCand->vertex+labelorg, tv->currorbit[NextCand->vertex]+labelorg); }

#define PRINT_INDIV_VERB(Verb) { if (tv->options->verbosity >= Verb) { \
PRINTCAND(CurrCand, tv->fromlevel) \
fprintf(outfile, "| "); \
fprintf(outfile, tv->digstring, NextCand->vertex+labelorg); \
} }

#define SPECIALGENERATORS { if (tv->options->generators) addpermutation(ring, AUTPERM, n); \
tv->stats->numgenerators++; \
tv->specialgens++; \
if (tv->options->writeautoms) { \
fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators); \
writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n); \
} \
if (tv->options->userautomproc) { \
(*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n); \
} }

#define UPDATEMIN(A, B) if (B < A) A = B;

#define PAIRORBJOIN(A, V) { if (A != V) { \
PrmPairs[tv->permInd].arg = A; \
PrmPairs[tv->permInd].val = V; \
tv->permInd++; \
orbjoin_sp_pair(tv->orbits, OrbList, n, A, V, &tv->stats->numorbits); \
MakeTree(A, V, sg, n, tv, FALSE); \
} }

#define SETPAIRS(A, V) { if (A != V) { \
PrmPairs[tv->permInd].arg = A; \
PrmPairs[tv->permInd].val = V; \
tv->permInd++; \
} }

#define SETPAIRSAUT(A, V) { if ((A != V) && (AUTPERM[A] != V)) { \
AUTPERM[A] = V; \
PrmPairs[tv->permInd].arg = A; \
PrmPairs[tv->permInd].val = V; \
tv->permInd++; \
} }

#define SETPAIRSAUTANDTREE(arg, val) { if (tv->build_autom) { SETPAIRSAUT(arg, val) } \
if (arg != val) orbjoin_sp_pair(tv->orbits, OrbList, n, arg, val, &tv->stats->numorbits); \
MakeTree(arg, val, sg_orig, n, tv, FALSE); }

#define SETPAIRSAUTANDTREE_PREPROC(arg, val) { if (tv->build_autom) { SETPAIRSAUT(arg, val) } \
if (arg != val) orbjoin_sp_pair(tv->orbits, OrbList, n, arg, val, &tv->stats->numorbits); \
MakeTree(arg, val, sg, n, tv, FALSE); }

#define PRINTF2(A, B) if (tv->options->verbosity > 2) printf(A, B)
#define PRINTF2_2(A, B, C) if (tv->options->verbosity > 2) printf(A, B, C)
#define PRINTF2_3(A, B, C, D) if (tv->options->verbosity > 2) printf(A, B, C, D)
#define PRINTF2_4(A, B, C, D, E) if (tv->options->verbosity > 2) printf(A, B, C, D, E)

#define NEIGHCOUNT_SING HitClsInd = 0; \
labi = lab[ind0]; \
iend1int = TheGraph[labi].d; \
nghb =  TheGraph[labi].e; \
for (j1int = 0; j1int < iend1int; ++j1int) { \
    k = nghb[j1int]; \
    value = Part->inv[InvLab[k]]; \
    if (cls[value] > 1) { \
        if (Markers[value] != tv->mark) { \
            HitCls[HitClsInd++] = value; \
            Markers[value] = tv->mark; \
            ElmHitCll[value] = value; \
        } \
        HitVtx[ElmHitCll[value]++] = k; \
    } \
    else { \
        switch (variation) { \
            case 1: \
                longcode = MASHCOMM(longcode, value); \
                break; \
            default: \
                break; \
        } \
    } \
}

#define NEIGHCOUNT_SPARSE if (cls[ind0] == n) { \
    for (i = 0; i < n; i++) { \
        NghCounts[i] = TheGraph[i].d; \
    } \
    HitCls[0] = 0; \
    HitClsInd = 1; \
    ElmHitCll[0] = 0; \
    for (i = 0; i < n; i++) { \
        NghCounts[i] = TheGraph[i].d; \
        if (TheGraph[i].d) { \
            HitVtx[ElmHitCll[0]++] = i; \
        } \
    } \
} \
else { \
    HitClsInd = 0; \
    for (i = ind0; i < ind2; i++) { \
        labi = lab[i]; \
        nghb = TheGraph[labi].e; \
        j = TheGraph[labi].d; \
        while (j-- > 0) { \
            k = nghb[j]; \
            if (MarkHitVtx[k] == tv->mark) { \
                NghCounts[k]++; \
            } \
            else { \
                value = Part->inv[InvLab[k]]; \
                if (cls[value] > 1) { \
                    MarkHitVtx[k] = tv->mark; \
                    NghCounts[k] = 1; \
                    if (Markers[value] != tv->mark) { \
                        HitCls[HitClsInd++] = value; \
                        Markers[value] = tv->mark; \
                        HitVtx[value] = k; \
                        ElmHitCll[value] = 1; \
                    } \
                    else { \
                        HitVtx[value+ElmHitCll[value]++] = k; \
                    } \
                } \
                else { \
                    switch (variation) { \
                        case 1: \
                            longcode = MASHCOMM(longcode, value); \
                            break; \
                        default: \
                            break; \
                    } \
                } \
            } \
        } \
    } \
}

#define NEIGHCOUNT_DENSE_SPARSE \
HitClsInd = 0; \
for (i = ind0; i < ind2; i++) { \
    labi = lab[i]; \
    iend1int = TheGraph[labi].d; \
    nghb =  TheGraph[labi].e; \
    for (j1int = 0; j1int < iend1int; ++j1int) { \
        k = nghb[j1int]; \
        (NghCounts[k])++; \
        value = Part->inv[InvLab[k]]; \
        if (Markers[value] != tv->mark) { \
            if (cls[value] > 1) HitCls[HitClsInd++] = value; \
            Markers[value] = tv->mark; \
        } \
    } \
}

#define NEIGHCOUNT_DENSE_DENSE for (i = ind0; i < ind2; i++) { \
    labi = lab[i]; \
    iend1int = TheGraph[lab[i]].d; \
    nghb = TheGraph[labi].e; \
    for (j1int = 0; j1int < iend1int; ++j1int) { \
        k = nghb[j1int]; \
        (NghCounts[k])++; \
    } \
}

#if !MAXN
DYNALLSTAT(int, AUTPERM, AUTPERM_sz);
DYNALLSTAT(int, BreakSteps, BreakSteps_sz);
DYNALLSTAT(int, CurrOrbSize, CurrOrbSize_sz);
DYNALLSTAT(int, CurrRefCells, CurrRefCells_sz);
DYNALLSTAT(int, fix, fix_sz);
DYNALLSTAT(int, IDENTITY_PERM, IDENTITY_PERM_sz);
DYNALLSTAT(int, Markers, Markers_sz);
DYNALLSTAT(int, TreeMarkers, TreeMarkers_sz);
DYNALLSTAT(int, AutMarkers, AutMarkers_sz);
DYNALLSTAT(int, MarkHitVtx, MarkHitVtx_sz);
DYNALLSTAT(int, MultRefCells, MultRefCells_sz);
DYNALLSTAT(int, NghCounts, NghCounts_sz);
DYNALLSTAT(int, OrbSize, OrbSize_sz);
DYNALLSTAT(int, OrbList, OrbList_sz);
DYNALLSTAT(pair, PrmPairs, PrmPairs_sz);
DYNALLSTAT(int, TempOrbList, TempOrbList_sz);
DYNALLSTAT(int, RefCells, RefCells_sz);
DYNALLSTAT(searchtrie*, RefPath, RefPath_sz);
DYNALLSTAT(int, Singletons, Singletons_sz);
DYNALLSTAT(int, SplCls, SplCls_sz);
DYNALLSTAT(int, SplCnt, SplCnt_sz);
DYNALLSTAT(int, SplPos, SplPos_sz);
DYNALLSTAT(int, StackMarkers, StackMarkers_sz);
DYNALLSTAT(int, TheTrace, TheTrace_sz);
DYNALLSTAT(int, TheTraceCC, TheTraceCC_sz);
DYNALLSTAT(int, TheTraceSplNum, TheTraceSplNum_sz);
DYNALLSTAT(int, TheTraceSteps, TheTraceSteps_sz);
DYNALLSTAT(int, TEMPLAB, TEMPLAB_sz);
DYNALLSTAT(int, TEMPINVLAB, TEMPINVLAB_sz);
DYNALLSTAT(int, WorkArray, WorkArray_sz);
DYNALLSTAT(int, WorkArray0, WorkArray0_sz);
DYNALLSTAT(int, WorkArray1, WorkArray1_sz);
DYNALLSTAT(int, WorkArray2, WorkArray2_sz);
DYNALLSTAT(int, WorkArray3, WorkArray3_sz);
DYNALLSTAT(int, WorkArray4, WorkArray4_sz);
DYNALLSTAT(int, WorkArray5, WorkArray5_sz);
DYNALLSTAT(int, WorkArray6, WorkArray6_sz);
DYNALLSTAT(int, WorkArray7, WorkArray7_sz);
DYNALLSTAT(int, Neighbs1, Neighbs1_sz);
DYNALLSTAT(int, Neighbs2, Neighbs2_sz);
DYNALLSTAT(int, TreeStack, TreeStack_sz);
DYNALLSTAT(TracesSpine, Spine, Spine_sz);
DYNALLSTAT(trie*, TrieArray, TrieArray_sz);
DYNALLSTAT(grph_strct, TheGraph, TheGraph_sz);
#else
static TLS_ATTR int AUTPERM[MAXN];
static TLS_ATTR int BreakSteps[MAXN];
static TLS_ATTR int CurrOrbSize[MAXN];
static TLS_ATTR int CurrRefCells[MAXN];
static TLS_ATTR int fix[MAXN];
static TLS_ATTR int IDENTITY_PERM[MAXN];
static TLS_ATTR int Markers[MAXN];
static TLS_ATTR int TreeMarkers[MAXN];
static TLS_ATTR int AutMarkers[MAXN];
static TLS_ATTR int MarkHitVtx[MAXN];
static TLS_ATTR int MultRefCells[MAXN];
static TLS_ATTR int NghCounts[MAXN];
static TLS_ATTR int OrbSize[MAXN];
static TLS_ATTR int OrbList[MAXN];
static TLS_ATTR pair PrmPairs[MAXN];
static TLS_ATTR int TempOrbList[MAXN];
static TLS_ATTR int RefCells[MAXN];
static TLS_ATTR searchtrie* RefPath[MAXN];
static TLS_ATTR int Singletons[MAXN];
static TLS_ATTR int SplCls[MAXN];
static TLS_ATTR int SplCnt[MAXN];
static TLS_ATTR int SplPos[MAXN];
static TLS_ATTR int StackMarkers[MAXN];
static TLS_ATTR int TheTrace[MAXN];
static TLS_ATTR int TheTraceCC[MAXN];
static TLS_ATTR int TheTraceSplNum[MAXN];
static TLS_ATTR int TheTraceSteps[MAXN];
static TLS_ATTR int TEMPLAB[MAXN];
static TLS_ATTR int TEMPINVLAB[MAXN];
static TLS_ATTR int WorkArray[MAXN];
static TLS_ATTR int WorkArray0[MAXN];
static TLS_ATTR int WorkArray1[MAXN];
static TLS_ATTR int WorkArray2[MAXN];
static TLS_ATTR int WorkArray3[MAXN];
static TLS_ATTR int WorkArray4[MAXN];
static TLS_ATTR int WorkArray5[MAXN];
static TLS_ATTR int WorkArray6[MAXN];
static TLS_ATTR int WorkArray7[MAXN];
static TLS_ATTR int TreeStack[MAXN];
static TLS_ATTR TracesSpine Spine[MAXN];
static TLS_ATTR trie* TrieArray[MAXN];
static TLS_ATTR grph_strct TheGraph[MAXN];
static TLS_ATTR int Neighbs1[MAXN];
static TLS_ATTR int Neighbs2[MAXN];
#endif

static TLS_ATTR FILE *outfile;

/* Brendan's SCHREIER */
static TLS_ATTR schreier  *gpB;				/* This will point to the Schreier structure */
static TLS_ATTR permnode  *gensB;			/* This will point to the stored generators */

static TLS_ATTR Candidate *GarbList, *SpOrd, *SpCyc, *SpSwp;
static TLS_ATTR Partition *SpPart1, *SpPart2;
static TLS_ATTR TracesSpine *SpineTL, *SpineFL, *SpineTL_tl;
static TLS_ATTR trie *trieroot, *trieref;
static TLS_ATTR int *TempOrbits = NULL;
static TLS_ATTR sparsegraph redgraph;

/* Check if the permutations in the list gens are automorphisms, 
 * also set mark and refcount fields and initialise orbits. */
static int
given_gens(sparsegraph *g, permnode *gens, int *orbits, boolean digraph) {
	int i, m, n, norbs;
	permnode *pn;
	
	n = g->nv;
	for (i = 0; i < n; ++i) orbits[i] = i;
	memcpy(IDENTITY_PERM, orbits, n*sizeof(int));
	norbs = n;
	
	if (!gens) return norbs;
	
	m = SETWORDSNEEDED(n);
	pn = gens;
	do {
		if (!isautom_sg((graph*)g, pn->p, digraph, m, n)) {
			fprintf(ERRFILE, "Input permutation is not an automorphism\n");
			exit(1);
		}
	    norbs = orbjoin(orbits, pn->p, n);
	    pn->mark = 1;
	    pn->refcount = 0;
	    pn = pn->next;
	} while (pn != gens);
	
	return norbs;
}

void
Traces(sparsegraph *g_arg, int *lab, int *ptn,
	   int *orbits_arg, TracesOptions *options_arg, TracesStats *stats_arg,
	   sparsegraph *canong_arg) {
	int i, j;
	int tmp;
    
    trielist *STStart, *STAux;
    searchtrie *TrieNode;
    int retval;
    Partition *CurrPart, *NextPart;
    Candidate *CurrCand, *NextCand, *BestCand, *AuxCand;
	
	const int n = g_arg->nv;
	const int m = SETWORDSNEEDED(n);
	
    if (g_arg->nv > (NAUTY_INFINITY-2))
    {
        fprintf(ERRFILE, "Traces: need n <= %d, but n=%d\n\n",
                NAUTY_INFINITY-2, g_arg->nv);
        return;
    }
    
	struct TracesVars *tv = malloc(sizeof(struct TracesVars));
    if (tv == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
	struct TracesInfo *ti = malloc(sizeof(struct TracesInfo));
    if (ti == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
    
	trieroot = NULL;
	NextCand = GarbList = NULL;
	
    tv->canlist = 0;
    tv->compstage = 0;
    tv->expathlength = n;
    tv->finalnumcells = n;
    tv->firstpathlength = 0;
    tv->linelgth = 0;
    tv->name = 0;
    tv->maxspineorblevel = 0;
    tv->nfix = 0;
    tv->options = options_arg;
    tv->orbits = orbits_arg;
    tv->permInd = 0;
    
    if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc)
        tv->build_autom = TRUE;
    else
        tv->build_autom = FALSE;
    
    //    tv->smalldeglevel = n;
    PRINTF2("Traces 1<: finalnumcells: %d\n", tv->finalnumcells);
    
    tv->specialgens = 0;
    tv->stats = stats_arg;
    tv->treedepth = 0;
    tv->gotonode = NULL;
    tv->input_graph = tv->graph = g_arg;
    tv->cangraph = canong_arg;
    tv->mark = tv->stackmark = tv->treemark = tv->autmark = tv->markcell1 = tv->markcell2 = 1073741825;
    tv->conta0 = tv->conta1 = tv->conta2 = tv->conta3 = tv->conta4 = tv->conta5 = tv->conta6 = tv->conta7 = tv->contatc = 0;
    //    tv->auxtime1 = tv->auxtime2 = tv->auxtime3 = tv->auxtime4 = tv->auxtime5 = 0;
    
	if (tv->options->strategy == 0) {
		tv->steps = n;
		tv->strategy = 0;
	}
	else {
		tv->strategy = 1;
        tv->steps = tv->options->strategy;
        if (tv->steps > n) {
            tv->steps = n;
        }
	}
    
#if !MAXN
    DYNALLOC1(int, AUTPERM, AUTPERM_sz, n, "Traces");
    DYNALLOC1(int, BreakSteps, BreakSteps_sz, n, "Traces");
    DYNALLOC1(int, CurrOrbSize, CurrOrbSize_sz, n, "Traces");
    DYNALLOC1(int, CurrRefCells, CurrRefCells_sz, n, "Traces");
    DYNALLOC1(int, fix, fix_sz, n, "Traces");
    DYNALLOC1(int, IDENTITY_PERM, IDENTITY_PERM_sz, n, "Traces");
	DYNALLOC1(int, Markers, Markers_sz, n, "Traces");
	DYNALLOC1(int, TreeMarkers, TreeMarkers_sz, n, "Traces");
	DYNALLOC1(int, AutMarkers, AutMarkers_sz, n, "Traces");
	DYNALLOC1(int, MarkHitVtx, MarkHitVtx_sz, n, "Traces");
    DYNALLOC1(int, MultRefCells, MultRefCells_sz, n, "Traces");
	DYNALLOC1(int, NghCounts, NghCounts_sz, n, "Traces");
    DYNALLOC1(int, OrbSize, OrbSize_sz, n, "Traces");
    DYNALLOC1(int, OrbList, OrbList_sz, n, "Traces");
    DYNALLOC1(pair, PrmPairs, PrmPairs_sz, n, "Traces");
    DYNALLOC1(int, TempOrbList, TempOrbList_sz, n, "Traces");
    DYNALLOC1(int, RefCells, RefCells_sz, n, "Traces");
	DYNALLOC1(int, Singletons, Singletons_sz, n, "Traces");
	DYNALLOC1(int, SplCls, SplCls_sz, n, "Traces");
	DYNALLOC1(int, SplCnt, SplCnt_sz, n, "Traces");
	DYNALLOC1(int, SplPos, SplPos_sz, n, "Traces");
	DYNALLOC1(int, StackMarkers, StackMarkers_sz, n, "Traces");
    DYNALLOC1(int, TheTrace, TheTrace_sz, n+10, "Traces");
    DYNALLOC1(int, TheTraceCC, TheTraceCC_sz, n, "Traces");
    DYNALLOC1(int, TheTraceSplNum, TheTraceSplNum_sz, n, "Traces");
    DYNALLOC1(int, TheTraceSteps, TheTraceSteps_sz, n+10, "Traces");
    DYNALLOC1(int, TEMPLAB, TEMPLAB_sz, n, "Traces");
    DYNALLOC1(int, TEMPINVLAB, TEMPINVLAB_sz, n, "Traces");
	DYNALLOC1(int, WorkArray, WorkArray_sz, n, "Traces");
    DYNALLOC1(int, WorkArray0, WorkArray0_sz, n, "Traces");
    DYNALLOC1(int, WorkArray1, WorkArray1_sz, n, "Traces");
    DYNALLOC1(int, WorkArray2, WorkArray2_sz, n, "Traces");
    DYNALLOC1(int, WorkArray3, WorkArray3_sz, n, "Traces");
    DYNALLOC1(int, WorkArray4, WorkArray4_sz, n, "Traces");
    DYNALLOC1(int, WorkArray5, WorkArray5_sz, n, "Traces");
    DYNALLOC1(int, WorkArray6, WorkArray6_sz, n, "Traces");
    DYNALLOC1(int, WorkArray7, WorkArray7_sz, n, "Traces");
    DYNALLOC1(int, TreeStack, TreeStack_sz, n, "Traces");
    DYNALLOC1(TracesSpine, Spine, Spine_sz, n, "Traces");
    DYNALLOC1(trie*, TrieArray, TrieArray_sz, n, "Traces");
    DYNALLOC1(grph_strct, TheGraph, TheGraph_sz, n, "Traces");
#endif
    
	outfile = (tv->options->outfile == NULL ? stdout : tv->options->outfile);
    
	SpOrd = SpCyc = SpSwp = NULL;
	SpPart1 = SpPart2 = NULL;
	
    if (tv->options->verbosity >= 2) {
        for (i = n, tv->digits = 0; i > 0; i /= 10, ++tv->digits) {}
        sprintf(tv->digstring, "%s%dd ", "%", tv->digits);
    }
	
#define PERMSTACK WorkArray1
#define CYCLES WorkArray1
#define HitCls WorkArray2
#define CYCOLR WorkArray2
#define HitVtx WorkArray3
#define CYLGTH WorkArray3
#define CStack WorkArray4
#define CYMULT WorkArray4
#define HitCount WorkArray5
#define ElmHitCll WorkArray5
#define CYCPOS WorkArray5
#define CYCHIT TempOrbList
#define LGHATTR RefCells
#define CYCREP MultRefCells
#define TempOrbSize TEMPLAB
#define AutomCount TEMPINVLAB
#define CanonIndices MarkHitVtx
#define NSFCells NghCounts
#define TreeNodes AutMarkers
#define CellMarkers1 WorkArray6
#define CellMarkers2 WorkArray7
    
	/* Initialize group and statistics */
    STATS_INIT
    if (tv->options->verbosity >= 2) {
        TIME_INIT
    }
    
	if (tv->options->generators) {
		tv->stats->numorbits = given_gens(g_arg, *tv->options->generators, orbits_arg, tv->options->digraph);
		newgroup(&gpB, NULL, n);
		gensB = *tv->options->generators;
		expandschreier(gpB, &gensB, n);
		ti->thegrouphaschanged = TRUE;
		ti->identitygroup = FALSE;
        // FARE ORBLIST ???
	}
	else {
		newgroup(&gpB, &gensB, n);
		for (i = 0; i < n; ++i) {
            orbits_arg[i] = i;
        }
        
		memcpy(IDENTITY_PERM, orbits_arg, n*sizeof(int));
		memcpy(OrbList, IDENTITY_PERM, n*sizeof(int));
		tv->stats->numorbits = n;
		ti->thegrouphaschanged = FALSE;
		ti->identitygroup = TRUE;
	}
    
	tv->currorbit = gpB->orbits;
	
    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
    tv->permInd = 0;
    
	memset(fix, 0, n*sizeof(int));
	memset(TheTraceCC, 0, n*sizeof(int));
	/* ran_init(1234);  any long int as an argument */
    
	/* The graph is sparse? */
	if (g_arg->nde < n || g_arg->nde / n < n / (g_arg->nde / n)) {
		ti->thegraphisparse = TRUE;
	}
	else {
		ti->thegraphisparse = FALSE;
	}
	
    tv->preprocessed = 0;
    ti->deg_one = FALSE;
    retval = 0;
    
	/* Initialize candidate, partition, cells, orbits */
	NEWPART(Spine[0].part);
	Spine[0].part->code = -1;
	
	NEWPART(NextPart);
	CurrPart = Spine[0].part;
	
	CurrCand = NewCandidate(n, &GarbList, TRUE);
	memset(CurrPart->inv, 0, n*sizeof(int));
	
	CurrCand->singcode = 0;
    TempOrbits = NULL;
    STStart = NULL;
    
	if (tv->options->defaultptn) {
		memcpy(CurrCand->lab, IDENTITY_PERM, n*sizeof(int));
		memcpy(CurrCand->invlab, IDENTITY_PERM, n*sizeof(int));
		CurrPart->cells = 1;
		CurrPart->cls[0] = n;
		TheTrace[0] = 0;
	}
	else {
		memcpy(CurrCand->lab, lab, n*sizeof(int));
		CurrPart->cells = 0;
		j = 0;
		for (i = 0; i < n; i++) {
			if (j) {
				CurrPart->inv[i] = j;
			}
			CurrCand->invlab[CurrCand->lab[i]] = i;
			if (!ptn[i]) {
				CurrPart->cls[j] = i-j+1;
				if (CurrPart->cls[j] == 1) {
					CurrCand->singcode = MASHCOMM(CurrCand->singcode, CurrCand->lab[j]);
				}
				TheTrace[CurrPart->cells++] = j;
				j = i+1;
			}
		}
	}
    
    ti->deg_one = traces_degree_refine(g_arg, CurrCand, n, CurrPart, tv, ti, 0);
    copy_sg_structure(g_arg, &redgraph);
    
    tv->graph = &redgraph;
    
    tv->maxdeg = 0;
    for (i=0; i<n; i++) {
        TheGraph[i].d = g_arg->d[i];
        if (TheGraph[i].d > tv->maxdeg) {
            tv->maxdeg = TheGraph[i].d;
        }
        TheGraph[i].e = g_arg->e + g_arg->v[i];
        TheGraph[i].one = FALSE;
    }
    
#if !MAXN
    DYNALLOC1(int, Neighbs1, Neighbs1_sz, tv->maxdeg, "Traces");
    DYNALLOC1(int, Neighbs2, Neighbs2_sz, tv->maxdeg, "Traces");
#endif
    
    if (ti->deg_one) {
        for (i=0; i<n; i++) {
            TheGraph[i].e = tv->graph->e + g_arg->v[i];
        }
        memcpy(tv->graph->e, g_arg->e, tv->graph->elen*sizeof(int));
        tv->preprocessed = Preprocess(g_arg, &gensB, CurrCand, n, CurrPart, tv);
    }
    
    /* Initialization of Spine structure */
	SpineFL = Spine;
	SpineFL->tgtcell = SpineFL->tgtpos = 0;
	SpineFL->tgtend = n;
	SpineFL->tgtfrom = -1;
	SpineFL->trcstart = 0;
	SpineFL->trcend = CurrPart->active = CurrPart->cells;
	SpineFL->ccstart = SpineFL->ccend = 0;
	SpineFL->stpstart = SpineFL->stpend = 0;
	SpineFL->singstart = SpineFL->singend = 0;
	SpineFL->thetracexists = FALSE;
	SpineFL->liststart = SpineFL->listend = CurrCand;
    SpineFL->levelcounter = 1;
    SpineFL->keptcounter = 1;
    SpineFL->updates = 1;
	
	/* Further initializations */
	tv->maxtreelevel = 0;
	tv->tolevel = 0;
	tv->tcell = 0;
	UPDATE_LINELGTH
    
	/* First refinement */
	if (tv->preprocessed < 2)
        traces_refine(CurrCand, n, CurrPart, tv, ti, 0, FALSE);
    
    CurrCand->name = 0;
    
	if (CurrPart->cells == n) {
		tv->stats->canupdates++;
        /* CANONICAL FORM ? */
		if (tv->options->getcanon) {
			memcpy(lab, CurrCand->lab, n*sizeof(int));
			if (canong_arg) memcpy(CurrPart->inv, CurrCand->invlab, n*sizeof(int));
            if (tv->options->usercanonproc != NULL)
            {
                (*tv->options->usercanonproc)((graph*)g_arg, lab, (graph*)canong_arg, tv->stats->canupdates, CurrCand->code, m, n);
            }
		}
	}
	else {
        STStart = searchtrie_new(n, tv);
        CurrCand->stnode = tv->strielist->triearray;
		NextCand = NewCandidate(n, &GarbList, TRUE);
		SpineFL->listcounter = 1;
		tv->tcellevel = 0;
//        printf("SET(0) tv->tcellevel to %d\n", tv->tcellevel);
		ti->newtc = FALSE;
        ti->exitfromref = FALSE;
        Spine[1].levelcounter = 0;
        Spine[1].updates = 0;
        Spine[1].tgtfrom = 0;
        
        memset(WorkArray, 0, n*sizeof(int));
        
		do {
			tv->fromlevel = tv->tolevel;
			SpineFL = Spine+tv->fromlevel;
            
			if (CurrCand) {
				switch (tv->compstage) {
					case 0: retval = CompStage0(CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
						break;
					case 1:
                        if (!TempOrbits) {
                            memcpy(WorkArray1, tv->currorbit, n*sizeof(int));
                            TempOrbits = WorkArray1;
                        }
                        memset(TempOrbSize, 0, n*sizeof(int));
                        for (i=SpineTL->tgtcell; i<SpineTL->tgtend; i++) {
                            TempOrbSize[TempOrbits[CurrCand->lab[i]]]++;
                        }
                        retval = CompStage1(CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
                        break;
					case 2:
                        retval = CompStage2(CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
                        break;
					default:
						break;
				}
                if (retval == NAUTY_ABORTED)
                    tv->stats->errstatus = NAUABORTED;
                else if (retval == NAUTY_KILLED) {
                    tv->stats->errstatus = NAUKILLED;
                    return;
                }
			}
            
			/* NEXT CANDIDATE */
			if (ti->thereisnextlevel) {
				if (tv->nextlevel != tv->fromlevel) {
					UPDATE_LINELGTH
				}
				tv->tolevel = tv->nextlevel;
				CurrCand = Spine[tv->nextlevel].liststart;
				CurrPart = Spine[tv->nextlevel].part;
			}
		}
		while (ti->thereisnextlevel);
        
        if (!retval) {
            if (tv->compstage) {
                memset(CurrOrbSize, 0, n*sizeof(int));
                for (i=0; i<n; i++) {
                    CurrOrbSize[TempOrbits[i]]++;
                }
            }
            
            if (!tv->options->getcanon) {
                if (tv->compstage) {
                    tv->maxtreelevel++;
                    TrieNode = Spine[tv->maxtreelevel].liststart->stnode;
                    if (TrieNode->father) {
                        TrieNode = TrieNode->father;
                    }
                    if (!ti->exitfromref) {
                        TrieNode->index = 0;
                        for (i=1; i<AutomCount[0]; i++) {
                            TrieNode->index += CurrOrbSize[AutomCount[i]];
                        }
                    }
                }
            }
            
            if (tv->maxtreelevel) PRINT_LINE_PLUS(tv->maxtreelevel);
            
            AuxCand = Spine[tv->maxtreelevel].liststart;
            while (!AuxCand->do_it) {
                AuxCand = AuxCand->next;
            }
            
            /* CANONICAL FORM ? */
            if (tv->options->getcanon) {
                BestCand = AuxCand;
                AuxCand = AuxCand->next;
                while (AuxCand) {
                    if (AuxCand->do_it) {
                        if (comparelab_tr(g_arg, BestCand->lab, BestCand->invlab, AuxCand->lab, AuxCand->invlab, Spine[tv->maxtreelevel].part->cls, Spine[tv->maxtreelevel].part->inv) == 1) {
                            BestCand = AuxCand;
                            tv->stats->canupdates++;
                            if (tv->options->usercanonproc != NULL)
                            {
                                (*tv->options->usercanonproc)((graph*)g_arg, lab, (graph*)canong_arg, tv->stats->canupdates, CurrCand->code, m, n);
                            }
                        }
                    }
                    AuxCand = AuxCand->next;
                }
                
                grouporderplus(g_arg, BestCand, Spine[tv->maxtreelevel].part, &gensB, &(tv->stats->grpsize1), &(tv->stats->grpsize2), n, tv, ti);
                
                if (tv->options->verbosity >= 2) {
                    LINE(tv->linelgth+4, "*")
                    NEXTLINE
                    fprintf(outfile, "Canonical:");
                    PRINTCAND(BestCand, tv->maxtreelevel)
                    PRINT_RETURN
                }
                memcpy(lab, BestCand->lab, n*sizeof(int));
                if (canong_arg) memcpy(CurrPart->inv, BestCand->invlab, n*sizeof(int));
            }
            else {
                grouporderplus(g_arg, AuxCand, Spine[tv->maxtreelevel].part, &gensB, &(tv->stats->grpsize1), &(tv->stats->grpsize2), n, tv, ti);
            }
            
            if (tv->options->verbosity >= 2) {
                if (tv->linelgth < 40) {
                    tv->linelgth = 40;
                }
                LINE(tv->linelgth+4, "*")
                NEXTLINE
            }
            
        }
    }
    tv->stats->treedepth = tv->treedepth;
    if (Spine[tv->treedepth].part->code == -1) {
        tv->stats->treedepth--;
    }
    
    if (tv->options->verbosity >= 2) {
        fprintf(outfile, "group time: %.2f, %.2f, %.2f; total: %.2f; exp_paths time: %.2f; aut_check time: %.2f\n%lu refinement%s interrupted by trace comparison (%s); special cells: %d\n------", tv->schreier1, tv->schreier2, tv->schreier3, tv->schreier1+tv->schreier2+tv->schreier3,
                tv->expaths, tv->autchk, SS(tv->stats->interrupted, "", "s"), (tv->options->strategy == 0 ? "breadth-first" : "depth-first"), tv->specialgens);
        PRINT_RETURN
    }
    if (tv->options->verbosity >= 2) fprintf(outfile, "CPYCAND(0): %d, ID<-TMPORB(1): %d, LAB(2): %d, PART(3): %d, TEMP(4)->: %d, TEMP(5)<-: %d, CHKFORAUT(6): %d, ISAUT(7): %d, ContaTC: %d\n", tv->conta0, tv->conta1, tv->conta2, tv->conta3, tv->conta4, tv->conta5, tv->conta6, tv->conta7, tv->contatc);
    
    if (tv->options->getcanon && canong_arg) {
        canong_arg->nv  = g_arg->nv;
        canong_arg->nde  = g_arg->nde;
        SG_ALLOC(*canong_arg, g_arg->nv, g_arg->nde, "traces canong");
        updatecan_tr(g_arg, canong_arg, lab, CurrPart->inv, 0);
    }
    
    if (tv->options->generators) {
        deleteunmarked(&gensB);
        *tv->options->generators = gensB;
        freeschreier(&gpB, NULL);
    }
    else {
        freeschreier(&gpB, &gensB);
        schreier_freedyn();
    }
    
    while (STStart) {
        STAux = STStart;
        free(STAux->triearray);
        STStart = STStart->next;
        free(STAux);
    }
    
    tv->canlist = 0;
    for (i=0; i<=tv->treedepth; i++) {
        if (Spine[i].liststart) {
            tv->canlist += FreeList(Spine[i].liststart, TRUE);
            Spine[i].liststart = Spine[i].listend = NULL;
        }
    }
    
    if (GarbList) {
        tv->stats->peaknodes = FreeList(GarbList, FALSE);
    }
    
    FREECAND(NextCand)
    FREECAND(SpOrd)
    FREECAND(SpCyc)
    FREECAND(SpSwp)
    FREEPART(NextPart)
    FREEPART(SpPart1)
    FREEPART(SpPart2)
    
    if (!tv->options->getcanon && trieroot) {
        for (i=0; i<=tv->triepos; i++) {
            free(TrieArray[i]);
        }
    }
    
    for (i=0; i <= tv->treedepth; i++) {
        FREEPART(Spine[i].part)
    }
    
    CurrCand = GarbList = NULL;
    tv->stats->peaknodes += tv->canlist;
    
    if (tv->graph != g_arg) {
        SG_FREE(redgraph);
    }
    free(tv);
    free(ti);
    traces_freedyn();
    return;
}

int traces_degree_refine(sparsegraph *sg, 
                         Candidate *Cand, 
                         int n, 
                         Partition *Part, 
                         struct TracesVars* tv, 
                         struct TracesInfo *ti, 
                         int num_indv) {
    
    int i, j, k;
    int ind0, ind1, ind2;
    int value, sc, iend;
//    int *sge;
    int HitClsInd, SplInd, SplCntInd, newcell, TraceInd;
    int *NghCounts, *HitCls;
//    boolean deg_one = FALSE;
    
    NghCounts = sg->d;
    HitCls = TheTrace;
    //PrintVect(HitCls, 0, Part->cells, 0);
    TraceInd = Part->cells;
    
    if (tv->mark > (NAUTY_INFINITY-2)) {
        memset(Markers, 0, n*sizeof(int));
        memset(MarkHitVtx, 0, n*sizeof(int));
        tv->mark = 0;
    }
    tv->mark++;
    
    memcpy(ElmHitCll, Part->cls, n*sizeof(int));
    memcpy(HitVtx, Cand->lab, n*sizeof(int));
    HitClsInd = TraceInd;
    
    SplInd = 0;
    SplCls[0] = n;
//    printf("HitClsInd: %d\n", HitClsInd);
    for (j = 0; j < HitClsInd; j++) {
        ind1 = HitCls[j];
        if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < Part->cls[ind1])) {
            SplCls[SplInd++] = ind1;
        }
        else {
            ind2 = ind1+Part->cls[ind1];
//            printf("ind1: %d; ind2: %d\n", ind1, ind2);
            value = NghCounts[Cand->lab[ind1++]];
            for (i = ind1; i < ind2; i++)
            {
//                if (sg->d[Cand->lab[i]] == 1) {
////                    printf("deg_one true A\n");
//                    deg_one = TRUE;
//                }
                if (NghCounts[Cand->lab[i]] != value)
                {
//                    printf("break @ %d\n", i);
                    SplCls[SplInd++] = HitCls[j];
                    break;
                }
            }
        }
    }
    
    if (SplInd) {
        
        /* Sorting the cells to be split */
        switch (SplInd) {
            case 0:
            case 1:
                break;
            case 2:
                if (SplCls[0] > SplCls[1]) {
                    value = SplCls[0];
                    SplCls[0] = SplCls[1];
                    SplCls[1] = value;
                }
                break;
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
                for (k = 1; k < SplInd; ++k) {
                    value = SplCls[k];
                    i = k - 1;
                    while ((i >= 0) && (value < SplCls[i])) {
                        SplCls[i + 1] = SplCls[i];
                        --i;
                    }
                    SplCls[i + 1] = value;
                }
                break;
            default:
                quickSort(SplCls, SplInd);
                break;
        }
        
        for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
            ind0 = SplCls[sc];
            ind1 = ind0 + Part->cls[ind0];
            SplCntInd = 0;
            if (ElmHitCll[ind0] < Part->cls[ind0]) {
                SplCnt[SplCntInd++] = 0;
                SplPos[0] = Part->cls[ind0] - ElmHitCll[ind0];
            }
            
            /* According to the degree of C */
            /* compute how many vertices in C will be placed into the same new cell */
            iend = ind0 + ElmHitCll[ind0];
            for (i = ind0; i < iend; i++) {
                value = NghCounts[HitVtx[i]];
                if (Markers[value] != tv->mark) {
                    Markers[value] = tv->mark;
                    SplCnt[SplCntInd++] = value;
                    SplPos[value] = 1;
                }
                else {
                    SplPos[value]++;
                }
            }
            tv->mark++;
            
            /* Sort the values deriving from the previous step */
            switch (SplCntInd) {
                case 0:
                case 1:
                    break;
                case 2:
                    if (SplCnt[0] > SplCnt[1]) {
                        value = SplCnt[0];
                        SplCnt[0] = SplCnt[1];
                        SplCnt[1] = value;
                    }
                    break;
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    for (k = 1; k < SplCntInd; ++k) {
                        value = SplCnt[k];
                        i = k - 1;
                        while ((i >= 0) && (value < SplCnt[i])) {
                            SplCnt[i + 1] = SplCnt[i];
                            --i;
                        }
                        SplCnt[i + 1] = value;
                    }
                    break;
                default:
                    quickSort(SplCnt, SplCntInd);
                    break;
            }
            
            Part->cells += SplCntInd-1;
            
            /* Split the cell C and update the information for sizes of new cells */
            i = ind0;
            
            for (k = 0; k < SplCntInd; k++) {
                value = SplPos[SplCnt[k]];
                Part->cls[i] = value;
                SplPos[SplCnt[k]] = i;
                i += value;
                if (i < ind1) {
                    //CStack[++CStackInd] = i;  // CSTACK??
                    TheTrace[TraceInd++] = i;
                }
            }
            //PrintVect(SplPos, 0, n, 0);
            /* Permute elements of the cell C */
            iend = ind0 + ElmHitCll[ind0];
            
            for (i = ind0; i < iend; i++) {
                value = HitVtx[i];
                j = SplPos[NghCounts[value]]++;    /* where HitVtx[i] goes */
                k = Cand->invlab[value];           /* where HitVtx[i] is in lab */
                //printf("%d %d\n", j, k);
                Cand->lab[k] = Cand->lab[j];
                Cand->lab[j] = value;
                Cand->invlab[value] = j;
                Cand->invlab[Cand->lab[k]] = k;
                //NghCounts[value] = 0;
            }
            
            //printf("deg(%d): %d\n", Cand->lab[ind0]+labelorg, sg->d[Cand->lab[ind0]]);
//            if (sg->d[Cand->lab[ind0]] == 1) {
//                printf("deg_one true B\n");
//                deg_one = TRUE;
//            }
            /* Reconstruct the cell C and update the inverse partition */
            newcell = ind1 - ElmHitCll[ind0];
            i = newcell;
            ind2 = newcell+Part->cls[newcell]-1;
            do {
                Part->inv[i] = newcell;
                if (i == ind2) {
                    newcell = i+1;
                    if (newcell < n) ind2 = newcell+Part->cls[newcell]-1;
                }
            }
            while (++i < ind1);
            
            //            for (i = ind0, k = 0; k < SplCntInd; i+=Part->cls[i], k++) {
            //                if (sg->d[Cand->lab[i]] == 1) {
            //                    printf("MOVE %d (vtx: %d) into stack[%d]\n", i, Cand->lab[i]+labelorg, ind);
            //                    CStack[ind++] = i;
            //                }
            //                if (Part->cls[i] == 1) {
            //                    Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
            //                    //                            if (newtrace)
            //                    Singletons[SingInd] = i;
            //                    SingInd++;
            //                    //                                printf("(C %d %d) ", i, Cand->lab[i]+labelorg);
            //                }
            //            }
            
        }
    }
    
    //    j=Part->cls[0];
    //    for (i=0; i<n; i++) {
    //        if (i==j) {
    //            printf("| ");
    //            j+=Part->cls[j];
    //        }
    //        printf("%d ", sg->d[Cand->lab[i]]);
    //    }
    //    printf("\n");
    for (i=0; i<n; i += Part->cls[i]) {
        if (sg->d[Cand->lab[i]] == 1) {
            return TRUE;
        }
    }
//    return deg_one;
    return FALSE;
}

int traces_refine(Candidate *Cand, 
				  //int m, 
				  int n, 
				  Partition *Part, 
				  struct TracesVars* tv, 
				  struct TracesInfo *ti, 
                  int num_indv, 
                  boolean make_code) {
	int i, j, k, sc, ind0, ind1, ind2, ind3, ind4, jk, tlp1, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd, SingInd;
//	size_t j1, iend1;
	int j1int, iend1int;
	unsigned int longcode;
	int newtrace = FALSE;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, *TraceEnd, Traceccend, *Tracestpend;
	int BigCell, BigCellPos, BigCellSize;
	boolean TraceCell = FALSE;
    int *nghb;
    int conta;
    const int variation = 0;
		
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	SpineTL = Spine+tv->tolevel;
	TraceEnd = &(SpineTL->trcend);
	Traceccend = SpineTL->ccend;
	Tracestpend = &(SpineTL->stpend);
	TraceCCInd = SpineTL->ccstart;
	TraceStepsInd = SpineTL->stpstart;
        
    SingInd = SpineTL->singstart + num_indv;
    
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
    //printf("TraceStart: %d\n", SpineTL->trcstart);
//    if (ti->deg_one && (SpineTL->trcstart == 0)) {
//        CStackInd = 0;
//        for (i = 0; i < Part->active; i++) {
//            if (TheGraph[lab[TheTrace[i]]].d > 1) {
//                CStack[++CStackInd] = TheTrace[i];
//                StackMarkers[CStackInd] = tv->stackmark;
//            }
//        }
//        printf("Cstack elems: %d/%d\n", CStackInd, Part->active);
//    }
//    else {
//        memcpy(CStack+1, TheTrace+SpineTL->trcstart, (Part->active)*sizeof(int));
//        CStackInd = Part->active;
//        for (i = 1; i <= CStackInd; i++) {
//            StackMarkers[CStack[i]] = tv->stackmark;
//        }
//    }

	memcpy(CStack+1, TheTrace+SpineTL->trcstart, (Part->active)*sizeof(int));
	CStackInd = Part->active;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
    TraceInd = SpineTL->trcstart+Part->active;
//    printf("\n[%u]", longcode);
	if (!SpineTL->thetracexists) {
		newtrace = TRUE;
	}
    conta=0;
    
//    printf("CStackInd: %d\n", CStackInd);
	while (CStackInd > 0) {
        
        //PrintVect(Singletons, 0, SingInd, 0);
        
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
        //		if (Part->cells == tv->finalnumcells) break;    ?????
        if (Part->cells == n) break;
		
//        printf("CStack: ");
//        PrintVect(CStack, 1, CStackInd+1, 0);

		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
        
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		StackMarkers[ind0] = 0;
		    
//        printf("CC is %d (vtx: %d, deg: %d --> %d)\n", ind0, lab[ind0]+labelorg, sg->d[lab[ind0]], TheGraph[lab[ind0]].d);
        //        printf("%7d) CC is %d (vtx: %d, deg: %d --> %d)\n", ++conta, ind0, lab[ind0]+labelorg, sg->d[lab[ind0]], TheGraph[lab[ind0]].d);
		if (!newtrace) {
			TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
        }
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {			/* SINGLETON CURRENT CELL CASE */
//            HitClsInd = 0;
//            if (sg->d[lab[ind0]] > 0) {    // RIVEDERE OVUNQUE
//                j1 = sg->v[lab[ind0]];
//                iend1 = j1+sg->d[lab[ind0]];
//                //
//                //printf("sg->v[%d]: %lu (deg: %d), iend1: %lu\n", lab[ind0]+labelorg, sg->v[lab[ind0]], sg->d[lab[ind0]], iend1);
//                for (; j1 < iend1; ++j1) {
//                    k = sg->e[j1];
//                    value = Part->inv[InvLab[k]];
//                    if (cls[value] > 1) {
//                        if (Markers[value] != tv->mark) {
//                            HitCls[HitClsInd++] = value;
//                            Markers[value] = tv->mark;
//                            ElmHitCll[value] = value;
//                        }
//                        HitVtx[ElmHitCll[value]++] = k;
//                    }
//                }
//            }
//            tv->mark++;
            NEIGHCOUNT_SING;
            
//            HitClsInd = 0;
//            labi = lab[ind0];
//			iend1int = TheGraph[labi].d;
//            nghb =  TheGraph[labi].e;
//            
//			for (j1int = 0; j1int < iend1int; ++j1int) {
//				k = nghb[j1int];
//				value = Part->inv[InvLab[k]];
//				if (cls[value] > 1) {
//					if (Markers[value] != tv->mark) {
//						HitCls[HitClsInd++] = value;
//						Markers[value] = tv->mark;
//						ElmHitCll[value] = value;
//					}
//					HitVtx[ElmHitCll[value]++] = k;
//				}
//			}
			tv->mark++;
            
            //            j1 = sg->v[labi];
            //            iend1 = j1+sg->d[labi];
            //            printf(" ]\nVerifica; &j1: %p\n[ ", &(sg->e[j1]));
            //            for (; j1 < iend1; ++j1) {
            //                printf("%d ", sg->e[j1]);
            //            }
            //            printf(" ]\n");
            
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
            /* SINGLETON CC CASE */
			if (SplInd) {
				if (newtrace) {
                    TheTraceCC[TraceCCInd] = ind0;
				}
                if (!TraceCell) {
                    TheTraceCC[TraceCCInd] = ind0;
                    newtrace = TRUE;
                }
                
                TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
                
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					i = ind1+cls[ind1]-ElmHitCll[ind1];
                    TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
				}
				
				/* REARRANGE THE CELLS */
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					cls[ind1] = cls[ind1]-ElmHitCll[ind1];
					newcell = ind1+cls[ind1];
					cls[newcell] = ElmHitCll[ind1];
					Part->cells++;
					
					if (StackMarkers[ind1] != tv->stackmark) {
						if (cls[newcell] < cls[ind1]) {
							CStack[++CStackInd] = newcell;
//                            printf("Add %d to stack (a) ", newcell);
//                            PrintVect(CStack, 1, CStackInd+1, 0);
							StackMarkers[newcell] = tv->stackmark;
						}
						else {
							CStack[++CStackInd] = ind1;
							StackMarkers[ind1] = tv->stackmark;
//                            printf("Add %d to stack (b) ", ind1);
//                            PrintVect(CStack, 1, CStackInd+1, 0);
						}
					}
					else {
						CStack[++CStackInd] = newcell;
//                        printf("Add %d to stack (c) ", newcell);
//                        PrintVect(CStack, 1, CStackInd+1, 0);
						StackMarkers[newcell] = tv->stackmark;
					}
					
					SplitCell = HitVtx+ind1;
					ind3 = cls[newcell];
					LabCell = lab+newcell;
					
					for (jk = 0; jk < ind3; jk++) {
						k = SplitCell[jk];
						i = LabCell[jk];
						Part->inv[newcell+jk] = newcell;
						lab[InvLab[k]] = i;
						InvLab[i] = InvLab[k];
						LabCell[jk] = k;
						InvLab[k] = newcell+jk;
					}
					if (cls[ind1] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[ind1]);
                        if (newtrace) Singletons[SingInd] = ind1;
                        SingInd++;
                        //                        printf("(A %d %d) ", ind1, Cand->lab[ind1]+labelorg);
                    }
					if (cls[newcell] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[newcell]);
                        //                        printf("(B %d %d) ", newcell, Cand->lab[newcell]+labelorg);
                        if (newtrace) Singletons[SingInd] = newcell;
                        SingInd++;
                    }
				}
			}
			else {
                if ((!newtrace) && TraceCell) {
//                    printf("return a\n");
                    return 0;
				}
			}
		}
		else {
			if (ti->thegraphisparse) {
//				if (cls[ind0] == n) {
//					memcpy(NghCounts, sg->d, n*sizeof(int));
//					HitCls[0] = 0;
//					HitClsInd = 1;
//					ElmHitCll[0] = 0;
//					for (i = 0; i < n; i++) {
//                        NghCounts[i] = TheGraph[i].d;
//						//if (sg->d[i]) {
//                        if (TheGraph[i].d) {
//							HitVtx[ElmHitCll[0]++] = i;
//						}
//					}
//				}
//				else {
//					HitClsInd = 0;
//                    for (i = ind0; i < ind2; i++) {
//                        labi = lab[i];
////                        nghb = sg->e+sg->v[labi];
////                        j = sg->d[labi];
//                        nghb = TheGraph[labi].e;
//                        j = TheGraph[labi].d;
//
//                        while (j-- > 0) {
//                            k = nghb[j];
//                            if (MarkHitVtx[k] == tv->mark) {
//                                NghCounts[k]++;
//                            }
//                            else {
//                                value = Part->inv[InvLab[k]];
//                                if (cls[value] > 1) {
//                                    MarkHitVtx[k] = tv->mark;
//                                    NghCounts[k] = 1;
//                                    if (Markers[value] != tv->mark) {
//                                        HitCls[HitClsInd++] = value;
//                                        Markers[value] = tv->mark;
//                                        HitVtx[value] = k;
//                                        ElmHitCll[value] = 1;
//                                    }
//                                    else {
//                                        HitVtx[value+ElmHitCll[value]++] = k;
//                                    }
//                                }
//                            }
//						}
//					}
//                    
//                }
                
				NEIGHCOUNT_SPARSE;
				tv->mark++;				
                
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = ind1;
					}
					else {
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
				}
                
                /* SPARSE CASE */
				if (SplInd) {
					if (newtrace) {
                        TheTraceCC[TraceCCInd] = ind0;
					}
                    if (!TraceCell) {
                        TheTraceCC[TraceCCInd] = ind0;
                        newtrace = TRUE;
                    }
                    TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+n, &Traceccend)
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
						ind0 = SplCls[sc];
						ind1 = ind0 + cls[ind0];
						SplCntInd = 0;
						if (ElmHitCll[ind0] < cls[ind0]) {
							SplCnt[SplCntInd++] = 0;
							SplPos[0] = cls[ind0] - ElmHitCll[ind0];
						}
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = NghCounts[HitVtx[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						tv->mark++;
						
						if (SplCntInd) {
                            TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+n, Tracestpend)
						}
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
//                                printf("Add %d to stack (d) ", i);
//                                PrintVect(CStack, 1, CStackInd+1, 0);
								StackMarkers[i] = tv->stackmark;
                                TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
							}
						}
						
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						/* Permute elements of the cell C */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = HitVtx[i];
							j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
							k = InvLab[value];				/* where HitVtx[i] is in lab */
							lab[k] = lab[j];
							lab[j] = value;
							InvLab[value] = j;
							InvLab[lab[k]] = k;
							NghCounts[value] = 0;
						}
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind1 - ElmHitCll[ind0];
						i = newcell;
						ind2 = newcell+cls[newcell]-1;
						do {
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
						
						for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
							if (cls[i] == 1) {
								Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                if (newtrace) Singletons[SingInd] = i;
                                SingInd++;
                                //                                printf("(C %d %d) ", i, Cand->lab[i]+labelorg);
                            }
						}
						
					}
				}
				else {
					if ((!newtrace) && TraceCell) {
//                        printf("return b\n");
                        return 0;
					}
				}
				
			}
			else {
//				if (sg->d[lab[ind0]] > n/cls[ind0]) {
                if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
//						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
//						memset(NghCounts, 0, n*sizeof(int));
//						HitClsInd = 0;
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//								value = Part->inv[InvLab[k]];
//								if (Markers[value] != tv->mark) {
//									if (cls[value] > 1) HitCls[HitClsInd++] = value;
//									Markers[value] = tv->mark;
//								}
//							}
//						}
//					}
//					
//					tv->mark++;

                    memset(NghCounts, 0, n*sizeof(int));
//                    HitClsInd = 0;
//                    for (i = ind0; i < ind2; i++) {
////                        j1 = sg->v[lab[i]];
//                        labi = lab[i];
//                        iend1int = TheGraph[labi].d;
//                        nghb =  TheGraph[labi].e;
//                        for (j1int = 0; j1int < iend1int; ++j1int) {
//                            k = nghb[j1int];
//                            (NghCounts[k])++;
//                            value = Part->inv[InvLab[k]];
//                            if (Markers[value] != tv->mark) {
//                                if (cls[value] > 1) HitCls[HitClsInd++] = value;
//                                Markers[value] = tv->mark;
//                            }
//                        }
//                    }
                        NEIGHCOUNT_DENSE_SPARSE;
                }
                
                tv->mark++;

                    
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
                    /* DENSE-SPARSE CASE */
					if (SplInd) {
						if (newtrace) {
                            TheTraceCC[TraceCCInd] = ind0;
						}
                        if (!TraceCell) {
                            TheTraceCC[TraceCCInd] = ind0;
                            newtrace = TRUE;
                        }
                        TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+2*n, &Traceccend)
						
						/* Sorting the cells to be split */
						switch (SplInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCls[0] > SplCls[1]) {
									value = SplCls[0];
									SplCls[0] = SplCls[1];
									SplCls[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplInd; ++k) {
									value = SplCls[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCls[i])) {
										SplCls[i + 1] = SplCls[i];
										--i;
									}
									SplCls[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCls, SplInd);
								break;
						}
						
						for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
							ind0 = SplCls[j];
							ind1 = ind0+cls[ind0];
							SplCntInd = 0;
							
							/* According to the numbers of neighbors of C into the current cell */
							/* compute how many vertices in C will be placed into the same new cell */
							for (i = ind0; i < ind1; i++) {
								value = NghCounts[lab[i]];
								if (Markers[value] != tv->mark) {
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = 1;
								}
								else {
									SplPos[value]++;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
                                TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+2*n, Tracestpend)
							}
							
							/* Sort the values deriving from the previous step */
							switch (SplCntInd) {
								case 0:
								case 1:
									break;
								case 2:
									if (SplCnt[0] > SplCnt[1]) {
										value = SplCnt[0];
										SplCnt[0] = SplCnt[1];
										SplCnt[1] = value;
									}
									break;
								case 3:
								case 4:
								case 5:
								case 6:
								case 7:
								case 8:
									for (k = 1; k < SplCntInd; ++k) {
										value = SplCnt[k];
										i = k - 1;
										while ((i >= 0) && (value < SplCnt[i])) {
											SplCnt[i + 1] = SplCnt[i];
											--i;
										}
										SplCnt[i + 1] = value;
									}
									break;
								default:
									quickSort(SplCnt, SplCntInd);
									break;
							}
							
							Part->cells += SplCntInd-1;
							
							/* Split the cell C and update the information for sizes of new cells */
							/* Put the new cells into the stack */
							i = ind0;
							if (StackMarkers[i] != tv->stackmark) {
								BigCellSize = 0;
							}
							for (k = 0; k < SplCntInd; k++) {
								value = SplPos[SplCnt[k]];
								cls[i] = value;
								
								if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
									BigCell = i;
									BigCellPos = CStackInd;
									BigCellSize = cls[i];
								}
								SplPos[SplCnt[k]] = i;
								i += value;
								if (i < ind1) {
									CStack[++CStackInd] = i;
//                                    printf("Add %d to stack (e) ", i);
//                                    PrintVect(CStack, 1, CStackInd+1, 0);
									StackMarkers[i] = tv->stackmark;
                                    TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
								}
							}
							
							if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
								CStack[BigCellPos] = ind0;
								StackMarkers[BigCell] = 0;
								StackMarkers[ind0] = tv->stackmark;
							}
							
							/* Permute elements of the cell C */
							i = ind0;
							do {
								SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
							}
							while(++i < ind1);
							
							/* Reconstruct the cell C and update the inverse partition */
							newcell = ind0;
							i = ind0;
							ind2 = newcell+cls[newcell]-1;
							do {
								lab[i] = SplCnt[i];
								InvLab[lab[i]] = i;
								
								Part->inv[i] = newcell;
								if (i == ind2) {
									newcell = i+1;
									if (newcell < n) ind2 = newcell+cls[newcell]-1;
								}
							}
							while (++i < ind1);
							
							for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
								if (cls[i] == 1) {
									Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                    if (newtrace) Singletons[SingInd] = i;
                                    SingInd++;
                                    //                                    printf("(D %d %d) ", i, Cand->lab[i]+labelorg);
                                }
							}
							
						}
					}
					else {
						if ((!newtrace) && TraceCell) {
//                            printf("return c\n");
                            return 0;
						}
					}
					
				}
				else {
					if (cls[ind0] == n) {
//						memcpy(NghCounts, sg->d, n*sizeof(int));
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
//						
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//							}
//						}

//						for (i = ind0; i < ind2; i++) {
//                            labi = lab[i];
//							iend1int = TheGraph[lab[i]].d;
//                            nghb = TheGraph[labi].e;
//							for (j1int = 0; j1int < iend1int; ++j1int) {
//								k = nghb[j1int];
//								(NghCounts[k])++;
//							}
//						}
                        NEIGHCOUNT_DENSE_DENSE;
					}
					SplInd = 0;
					ind4 = 0;
                    while (ind4 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind4+cls[ind4];
						if (cls[ind4] > 1) {
							
                            
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind4]];
							for (i = ind4+1; i < ind1; i++) {
								if (NghCounts[lab[i]] != value)
								{
									SplInd++;
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind4;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							tv->mark++;
							
                            if (SplInd && !TraceCell) newtrace = TRUE;
                            
							if (SplCntInd) {
                                TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+3*n, Tracestpend)
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind4;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									
									if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
//                                        printf("Add %d to stack (f) ", i);
//                                        PrintVect(CStack, 1, CStackInd+1, 0);
										StackMarkers[i] = tv->stackmark;
                                        TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
									}
								}
								if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
									CStack[BigCellPos] = ind4;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind4] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind4;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind4;
								i = ind4;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
								
								for (i = ind4; i < ind1; i+=cls[i]) {
									if (cls[i] == 1) {
										Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                        if (newtrace) Singletons[SingInd] = i;
                                        SingInd++;
                                        //                                        printf("(E %d %d) ", i, Cand->lab[i]+labelorg);
                                    }
								}
							}
						}
						ind4 = ind1;
					}
                    
                    /* DENSE-DENSE CASE */
					if (SplInd) {
                        if (!TraceCell) {
                            newtrace = TRUE;
                        }
                        TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+3*n, &Traceccend)
						if (newtrace) {
                            TheTraceCC[TraceCCInd-1] = ind0;
						}
					}
					else {
						if ((!newtrace) && TraceCell) {
//                            printf("return d\n");
                            return 0;
						}
					}
					
				}
			}
		}
	}  /* end while (CStackInd > 0) */
    if (make_code) {
//        PrintVect(TheTrace, 0, TraceInd, 0);
//        PrintPartition(lab, Part->cls, n, labelorg, 2681);
//        for (i=SpineTL->trcstart; i < TraceInd; i++) {
//            ind0 = TheTrace[i];
//            longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
//            //            printf("[%u]", longcode);
//            j1 = sg->v[lab[ind0]];
//            iend1 = j1+sg->d[lab[ind0]];
//            for (; j1 < iend1; ++j1) {
//                k = sg->e[j1];
//                value = Part->inv[InvLab[k]];
//                longcode = MASHCOMM(longcode, value);
//                //                printf("[%u]", longcode);
//            }
//        }
        for (i=SpineTL->trcstart; i < TraceInd; i++) {
            ind0 = TheTrace[i];
            longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
            labi = lab[ind0];
            iend1int = TheGraph[labi].d;
            nghb = TheGraph[labi].e;
            for (j1int = 0; j1int < iend1int; ++j1int) {
                k = nghb[j1int];
                value = Part->inv[InvLab[k]];
                longcode = MASHCOMM(longcode, value);
            }
        }
    }
 //   printf("Complete refinement\n");
	Part->code = Cand->code = CLEANUP(longcode);
    tlp1 = tv->tolevel+1;
	if (newtrace) {
        if ((tlp1 < n) && (tlp1 > tv->treedepth)) {
			tv->treedepth = tlp1;
            if (tv->strategy) {
                NEWPART(Spine[tlp1].part)
			}
			else {
                NEWPARTSPINE(tlp1)
			}
			Spine[tlp1].liststart = Spine[tlp1].listend = NULL;
			Spine[tlp1].listcounter = 0;
		}
		*TraceEnd = TraceInd;
		SpineTL->ccend = TraceCCInd;
		SpineTL->singend = SingInd;
        //        printf("Set SpineTL->singend = %d @ %d\n", SpineTL->singend, tv->tolevel);
		*Tracestpend = TraceStepsInd;
		
		SpineTL->thetracexists = TRUE;
		if (tlp1 < n) {
			Spine[tlp1].ccstart = TraceCCInd;
			Spine[tlp1].singstart = Spine[tlp1].singend = SingInd;
            //            printf("Set SpineTL->singstart = %d, SpineTL->singend = %d @ %d\n", Spine[tlp1].singstart, Spine[tlp1].singend, tlp1);
			Spine[tlp1].stpstart = TraceStepsInd;
			Spine[tlp1].thetracexists = FALSE;
			Spine[tlp1].part->code = -1;
		}
		return 2;
	}
	else {
        if (TraceInd < *TraceEnd) {
            return 0;
        }
		if (Cand->code > SpineTL->part->code) {
			return 2;
		}
		else {
			if (Cand->code < SpineTL->part->code) {
				return 0;
			}
		}
		return 1;
	}
}

void traces_refine_notrace(Candidate *Cand, 
						   //int m, 
						   int n, 
						   Partition *Part, 
						   struct TracesVars* tv, 
						   struct TracesInfo *ti) {
	int i, i1, j, k, sc, ind0, ind1, ind2, ind3, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd;
//	size_t j1, iend1;
    int iend1int, j1int;
	unsigned int longcode;
	int Split = 0;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *SplitCell, *LabCell;
	int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    const int variation = 1;
    
//    printf("ref_notrace(%d); \n", Part->cells);
 	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	CStackInd = 1;
    CStack[1] = tv->tcellexpath+cls[tv->tcellexpath];
	
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	
	while (CStackInd > 0)
	{
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
        
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];		/* Current Cell */
		longcode = MASHNONCOMM(longcode, ind0);
		StackMarkers[ind0] = 0;
		
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {	/* SINGLETON CURRENT CELL CASE */
//			HitClsInd = 0;
//			j1 = sg->v[lab[ind0]];
//			iend1 = j1+sg->d[lab[ind0]];
//			
//			for (; j1 < iend1; ++j1) {
//				k = sg->e[j1];
//				value = Part->inv[InvLab[k]];
//				if (cls[value] > 1) {
//					if (Markers[value] != tv->mark) {
//						HitCls[HitClsInd++] = value;
//						Markers[value] = tv->mark;
//						ElmHitCll[value] = value;
//					}
//					HitVtx[ElmHitCll[value]++] = k;
//				}
//				else {
//					longcode = MASHCOMM(longcode, value);
//				}
//			}

            NEIGHCOUNT_SING;
			tv->mark++;
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}

			switch (SplInd) {
				case 0:
				case 1:
					break;
				case 2:
					if (SplCls[0] > SplCls[1]) {
						value = SplCls[0];
						SplCls[0] = SplCls[1];
						SplCls[1] = value;
					}
					break;
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
					for (k = 1; k < SplInd; ++k) {
						value = SplCls[k];
						i = k - 1;
						while ((i >= 0) && (value < SplCls[i])) {
							SplCls[i + 1] = SplCls[i];
							--i;
						}
						SplCls[i + 1] = value;
					}
					break;
				default:
					quickSort(SplCls, SplInd);
					break;
			}
			
			/* REARRANGE THE CELL */
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				cls[ind1] = cls[ind1]-ElmHitCll[ind1];
				newcell = ind1+cls[ind1];
				cls[newcell] = ElmHitCll[ind1];
				
				Part->cells++;
				
				if (StackMarkers[ind1] != tv->stackmark) {
					if (cls[newcell] < cls[ind1]) {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					else {
						CStack[++CStackInd] = ind1;
						StackMarkers[ind1] = tv->stackmark;
					}
				}
				else {
					CStack[++CStackInd] = newcell;
					StackMarkers[newcell] = tv->stackmark;
				}
				
				SplitCell = HitVtx+ind1;
				ind3 = cls[newcell];
				LabCell = lab+newcell;
				
				for (i1 = 0; i1 < ind3; i1++) {
					k = SplitCell[i1];
					i = LabCell[i1];
					Part->inv[newcell+i1] = newcell;
					lab[InvLab[k]] = i;
					InvLab[i] = InvLab[k];
					LabCell[i1] = k;
					InvLab[k] = newcell+i1;
				}
                if (cls[ind1] == 1) {
                    Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[ind1]);
                }
                if (cls[newcell] == 1) {
                    Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[newcell]);
                }
			}
		}
		else {
			if (ti->thegraphisparse) {
//				if (cls[ind0] == n) {
//					memcpy(NghCounts, sg->d, n*sizeof(int));
//					HitCls[0] = 0;
//					HitClsInd = 1;
//					ElmHitCll[0] = 0;
//					for (i = 0; i < n; i++) {
//						if (sg->d[i]) {
//							HitVtx[ElmHitCll[0]++] = i;
//						}
//					}
//				}
//				else {
//					HitClsInd = 0;
//                    for (i = ind0; i < ind2; i++) {
//                        labi = lab[i];
//                        nghb = sg->e+sg->v[labi];
//                        j = sg->d[labi];
//                        //for (; j--;) {
//                        while (j-- > 0) {
//                            k = nghb[j];
//                            if (MarkHitVtx[k] == tv->mark) {
//                                NghCounts[k]++;
//                            }
//                            else {
//                                value = Part->inv[InvLab[k]];
//                                if (cls[value] > 1) {
//                                    MarkHitVtx[k] = tv->mark;
//                                    NghCounts[k] = 1;
//                                    if (Markers[value] != tv->mark) {
//                                        HitCls[HitClsInd++] = value;
//                                        Markers[value] = tv->mark;
//                                        HitVtx[value] = k;
//                                        ElmHitCll[value] = 1;
//                                    }
//                                    else {
//                                        HitVtx[value+ElmHitCll[value]++] = k;
//                                    }
//                                }
//                                else {
//                                    longcode = MASHCOMM(longcode, value);
//                                }
//                            }
//						}
//					}
//                    
//                }
				
                NEIGHCOUNT_SPARSE;
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = ind1;
					}
					else {
						ind2 = ind1+cls[ind1];
						Split = FALSE;
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								Split = TRUE;
								break;
							}
						}
						if (!Split) {
							longcode = MASHCOMM(longcode, ind1);
						}
					}
				}
				/* Sorting the cells to be split */
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
					ind0 = SplCls[sc];
					ind1 = ind0 + cls[ind0];
					SplCntInd = 0;
					if (ElmHitCll[ind0] < cls[ind0]) {
						SplCnt[SplCntInd++] = 0;
						SplPos[0] = cls[ind0] - ElmHitCll[ind0];
					}
					
					/* According to the numbers of neighbors of C into the current cell */
					/* compute how many vertices in C will be placed into the same new cell */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = NghCounts[HitVtx[i]];
						if (Markers[value] != tv->mark) {
							Markers[value] = tv->mark;
							SplCnt[SplCntInd++] = value;
							SplPos[value] = 1;
						}
						else {
							SplPos[value]++;
						}
					}
					tv->mark++;
					
					/* Sort the values deriving from the previous step */
					switch (SplCntInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCnt[0] > SplCnt[1]) {
								value = SplCnt[0];
								SplCnt[0] = SplCnt[1];
								SplCnt[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplCntInd; ++k) {
								value = SplCnt[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCnt[i])) {
									SplCnt[i + 1] = SplCnt[i];
									--i;
								}
								SplCnt[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCnt, SplCntInd);
							break;
					}
					Part->cells += SplCntInd-1;
					
					/* Split the cell C and update the information for sizes of new cells */
					/* Put the new cells into the stack */
					i = ind0;
					if (StackMarkers[i] != tv->stackmark) {
						BigCellSize = 0;
					}
					for (k = 0; k < SplCntInd; k++) {
						value = SplPos[SplCnt[k]];
						cls[i] = value;
						if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
							BigCell = i;
							BigCellPos = CStackInd;
							BigCellSize = cls[i];
						}
						SplPos[SplCnt[k]] = i;
						i += value;
						if (i < ind1) {
							CStack[++CStackInd] = i;
							StackMarkers[i] = tv->stackmark;
						}
					}
					
					if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
						CStack[BigCellPos] = ind0;
						StackMarkers[BigCell] = 0;
						StackMarkers[ind0] = tv->stackmark;
					}
					/* Permute elements of the cell C */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = HitVtx[i];
						j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
						k = InvLab[value];				/* where HitVtx[i] is in lab */
						lab[k] = lab[j];
						lab[j] = value;
						InvLab[value] = j;
						InvLab[lab[k]] = k;
						NghCounts[value] = 0;
					}
					
					/* Reconstruct the cell C and update the inverse partition */
					newcell = ind1 - ElmHitCll[ind0];
					i = newcell;
					ind2 = newcell+cls[newcell]-1;
					do {
						Part->inv[i] = newcell;
						if (i == ind2) {
							newcell = i+1;
							if (newcell < n) ind2 = newcell+cls[newcell]-1;
						}
					}
					while (++i < ind1);
                    
                    for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                        if (cls[i] == 1) {
                            Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                        }
                    }
                    
				}
			}
			else {
//				if (sg->d[lab[ind0]] > n/cls[ind0]) {
                if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
//						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
//						HitClsInd = 0;
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//								value = Part->inv[InvLab[k]];
//								if (Markers[value] != tv->mark) {
//									if (cls[value] > 1) HitCls[HitClsInd++] = value;
//									Markers[value] = tv->mark;
//								}
//							}
//						}
                        NEIGHCOUNT_DENSE_SPARSE
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
						ind0 = SplCls[j];
						ind1 = ind0+cls[ind0];
						SplCntInd = 0;
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						for (i = ind0; i < ind1; i++) {
							value = NghCounts[lab[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						
						tv->mark++;
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
							}
						}
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						
						/* Permute elements of the cell C */
						i = ind0;
						do {
							SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
						}
						while(++i < ind1);
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind0;
						i = ind0;
						ind2 = newcell+cls[newcell]-1;
						do {
							lab[i] = SplCnt[i];
							InvLab[lab[i]] = i;
							Part->inv[i] = newcell;
							
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
                        
                        for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                            if (cls[i] == 1) {
                                Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                            }
                        }
                        
					}
				}
				else {
					if (cls[ind0] == n) {
//						memcpy(NghCounts, sg->d, n*sizeof(int));
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//							}
//						}
                        NEIGHCOUNT_DENSE_DENSE;
					}
					
					ind0 = 0;
					while (ind0 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind0+cls[ind0];
						if (cls[ind0] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind0]];
							for (i = ind0+1; i < ind1; i++)
							{
								if (NghCounts[lab[i]] != value)
								{
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind0;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							
							if (SplCntInd) {
								tv->mark++;
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind0;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
									}
								}
								if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
									CStack[BigCellPos] = ind0;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind0] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind0;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind0;
								i = ind0;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
                                
                                for (i = ind0; i < ind1; i+=cls[i]) {
									if (cls[i] == 1) {
										Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                                    }
								}
                                
							}
						}
						ind0 = ind1;
					}
				}
			}
		}
	}  /* end while (CStackInd > 0) */
	
	Cand->code = CLEANUP(longcode);
	return;
}

void traces_refine_maketrie(Candidate *Cand, 
							//int m, 
							int n, 
							Partition *Part, 
							struct TracesVars* tv, 
							struct TracesInfo *ti) {
	int i, i1, j, k, sc, ind0, ind1, ind2, ind3, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd;
//	size_t j1, iend1;
	int j1int, iend1int;
	unsigned int longcode;
	int Split = 0;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *SplitCell, *LabCell;
	int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    const int variation = 1;
    
	//TcSize = Spine[tv->tolevel].tgtsize;
	
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	CStack[1] = Spine[tv->tolevel_tl].tgtpos;
	CStackInd = 1;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	
	while (CStackInd > 0)
	{
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
        //		if (Part->cells == tv->finalnumcells) break;    ?????
        if (Part->cells == n) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		longcode = MASHNONCOMM(longcode, ind0);
		StackMarkers[ind0] = 0;
		
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
//			HitClsInd = 0;
//			j1 = sg->v[lab[ind0]];
//			iend1 = j1+sg->d[lab[ind0]];
//			
//			for (; j1 < iend1; ++j1) {
//				k = sg->e[j1];
//				value = Part->inv[InvLab[k]];
//				if (cls[value] > 1) {
//					if (Markers[value] != tv->mark) {
//						HitCls[HitClsInd++] = value;
//						Markers[value] = tv->mark;
//						ElmHitCll[value] = value;
//					}
//					HitVtx[ElmHitCll[value]++] = k;
//				}
//				else {
//					longcode = MASHCOMM(longcode, value);
//				}
//			}

            NEIGHCOUNT_SING;
			tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
			switch (SplInd) {
				case 0:
				case 1:
					break;
				case 2:
					if (SplCls[0] > SplCls[1]) {
						value = SplCls[0];
						SplCls[0] = SplCls[1];
						SplCls[1] = value;
					}
					break;
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
					for (k = 1; k < SplInd; ++k) {
						value = SplCls[k];
						i = k - 1;
						while ((i >= 0) && (value < SplCls[i])) {
							SplCls[i + 1] = SplCls[i];
							--i;
						}
						SplCls[i + 1] = value;
					}
					break;
				default:
					quickSort(SplCls, SplInd);
					break;
			}
			
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				i = ind1+cls[ind1]-ElmHitCll[ind1];
				trieref = trie_make(trieref, i, n, tv);
			}
			
			/* REARRANGE THE CELL */
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				cls[ind1] = cls[ind1]-ElmHitCll[ind1];
				newcell = ind1+cls[ind1];
				cls[newcell] = ElmHitCll[ind1];
				
				Part->cells++;
				
				if (StackMarkers[ind1] != tv->stackmark) {
					if (cls[newcell] < cls[ind1]) {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					else {
						CStack[++CStackInd] = ind1;
						StackMarkers[ind1] = tv->stackmark;
					}
				}
				else {
					CStack[++CStackInd] = newcell;
					StackMarkers[newcell] = tv->stackmark;
				}
				
				SplitCell = HitVtx+ind1;
				ind3 = cls[newcell];
				LabCell = lab+newcell;
				
				for (i1 = 0; i1 < ind3; i1++) {
					k = SplitCell[i1];
					i = LabCell[i1];
					Part->inv[newcell+i1] = newcell;
					lab[InvLab[k]] = i;
					InvLab[i] = InvLab[k];
					LabCell[i1] = k;
					InvLab[k] = newcell+i1;
				}
			}
		}
		else {
			if (ti->thegraphisparse) {
//				if (cls[ind0] == n) {
//					memcpy(NghCounts, sg->d, n*sizeof(int));
//					HitCls[0] = 0;
//					HitClsInd = 1;
//					ElmHitCll[0] = 0;
//					for (i = 0; i < n; i++) {
//						if (sg->d[i]) {
//							HitVtx[ElmHitCll[0]++] = i;
//						}
//					}
//				}
//				else {
//					HitClsInd = 0;
//                    for (i = ind0; i < ind2; i++) {
//                        labi = lab[i];
//                        nghb = sg->e+sg->v[labi];
//                        j = sg->d[labi];
//                        //for (; j--;) {
//                        while (j-- > 0) {
//                            k = nghb[j];
//                            if (MarkHitVtx[k] == tv->mark) {
//                                NghCounts[k]++;
//                            }
//                            else {
//                                value = Part->inv[InvLab[k]];
//                                if (cls[value] > 1) {
//                                    MarkHitVtx[k] = tv->mark;
//                                    NghCounts[k] = 1;
//                                    if (Markers[value] != tv->mark) {
//                                        HitCls[HitClsInd++] = value;
//                                        Markers[value] = tv->mark;
//                                        HitVtx[value] = k;
//                                        ElmHitCll[value] = 1;
//                                    }
//                                    else {
//                                        HitVtx[value+ElmHitCll[value]++] = k;
//                                    }
//                                }
//                                else {
//                                    longcode = MASHCOMM(longcode, value);
//                                }
//                            }
//						}
//					}
//                    
//                }
                
				NEIGHCOUNT_SPARSE;
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = ind1;
					}
					else {
						ind2 = ind1+cls[ind1];
						Split = FALSE;
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								Split = TRUE;
								break;
							}
						}
						if (!Split) {
							longcode = MASHCOMM(longcode, ind1);
						}
					}
				}
				/* Sorting the cells to be split */
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
					ind0 = SplCls[sc];
					ind1 = ind0 + cls[ind0];
					SplCntInd = 0;
					if (ElmHitCll[ind0] < cls[ind0]) {
						SplCnt[SplCntInd++] = 0;
						SplPos[0] = cls[ind0] - ElmHitCll[ind0];
					}
					
					/* According to the numbers of neighbors of C into the current cell */
					/* compute how many vertices in C will be placed into the same new cell */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = NghCounts[HitVtx[i]];
						if (Markers[value] != tv->mark) {
							Markers[value] = tv->mark;
							SplCnt[SplCntInd++] = value;
							SplPos[value] = 1;
						}
						else {
							SplPos[value]++;
						}
					}
					
					tv->mark++;
					
					/* Sort the values deriving from the previous step */
					switch (SplCntInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCnt[0] > SplCnt[1]) {
								value = SplCnt[0];
								SplCnt[0] = SplCnt[1];
								SplCnt[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplCntInd; ++k) {
								value = SplCnt[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCnt[i])) {
									SplCnt[i + 1] = SplCnt[i];
									--i;
								}
								SplCnt[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCnt, SplCntInd);
							break;
					}
					Part->cells += SplCntInd-1;
					
					/* Split the cell C and update the information for sizes of new cells */
					/* Put the new cells into the stack */
					i = ind0;
					if (StackMarkers[i] != tv->stackmark) {
						BigCellSize = 0;
					}
					for (k = 0; k < SplCntInd; k++) {
						value = SplPos[SplCnt[k]];
						cls[i] = value;
						if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
							BigCell = i;
							BigCellPos = CStackInd;
							BigCellSize = cls[i];
						}
						SplPos[SplCnt[k]] = i;
						i += value;
						if (i < ind1) {
							CStack[++CStackInd] = i;
							StackMarkers[i] = tv->stackmark;
							trieref = trie_make(trieref, i, n, tv);
						}
					}
					
					if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
						CStack[BigCellPos] = ind0;
						StackMarkers[BigCell] = 0;
						StackMarkers[ind0] = tv->stackmark;
					}
					/* Permute elements of the cell C */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = HitVtx[i];
						j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
						k = InvLab[value];				/* where HitVtx[i] is in lab */
						lab[k] = lab[j];
						lab[j] = value;
						InvLab[value] = j;
						InvLab[lab[k]] = k;
						NghCounts[value] = 0;
					}
					
					/* Reconstruct the cell C and update the inverse partition */
					newcell = ind1 - ElmHitCll[ind0];
					i = newcell;
					ind2 = newcell+cls[newcell]-1;
					do {
						Part->inv[i] = newcell;
						if (i == ind2) {
							newcell = i+1;
							if (newcell < n) ind2 = newcell+cls[newcell]-1;
						}
					}
					while (++i < ind1);
				}
			}
			else {
//				if (sg->d[lab[ind0]] > n/cls[ind0]) {
                if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
//						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
//						HitClsInd = 0;
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//								value = Part->inv[InvLab[k]];
//								if (Markers[value] != tv->mark) {
//									if (cls[value] > 1) HitCls[HitClsInd++] = value;
//									Markers[value] = tv->mark;
//								}
//							}
//						}
                        NEIGHCOUNT_DENSE_SPARSE;
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
						ind0 = SplCls[j];
						ind1 = ind0+cls[ind0];
						SplCntInd = 0;
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						for (i = ind0; i < ind1; i++) {
							value = NghCounts[lab[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						
						tv->mark++;
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
								trieref = trie_make(trieref, i, n, tv);
							}
						}
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						
						/* Permute elements of the cell C */
						i = ind0;
						do {
							SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
						}
						while(++i < ind1);
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind0;
						i = ind0;
						ind2 = newcell+cls[newcell]-1;
						do {
							lab[i] = SplCnt[i];
							InvLab[lab[i]] = i;
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
					}
				}
				else {
					if (cls[ind0] == n) {
//						memcpy(NghCounts, sg->d, n*sizeof(int));
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//							}
//						}
                        NEIGHCOUNT_DENSE_DENSE;
					}
					
					ind0 = 0;
					while (ind0 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind0+cls[ind0];
						if (cls[ind0] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind0]];
							for (i = ind0+1; i < ind1; i++)
							{
								if (NghCounts[lab[i]] != value)
								{
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind0;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							
							if (SplCntInd) {
								tv->mark++;
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind0;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
									}
								}
								if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
									CStack[BigCellPos] = ind0;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind0] = tv->stackmark;
									trieref = trie_make(trieref, i, n, tv);
								}
								
								/* Permute elements of the cell C */
								i = ind0;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind0;
								i = ind0;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
							}
						}
						ind0 = ind1;
					}
				}
			}
		}
	}
	
	Cand->code = CLEANUP(longcode);
	return;
}

int traces_refine_comptrie(Candidate *Cand, 
						   //int m, 
						   int n, 
						   Partition *Part, 
						   struct TracesVars* tv, 
						   struct TracesInfo *ti) {
	int i, i1, j, k, sc, ind0, ind1, ind2, ind3, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd;
//	size_t j1, iend1;
	int j1int, iend1int;
	unsigned int longcode;
	int Split = 0;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *SplitCell, *LabCell;
	int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    const int variation = 1;
    
	//TcSize = Spine[tv->tolevel].tgtsize;
	
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	CStack[1] = Spine[tv->tolevel_tl].tgtpos;
	CStackInd = 1;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	while (CStackInd > 0)
	{
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
        //		if (Part->cells == tv->finalnumcells) break;    ?????
        if (Part->cells == n) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];		/* Current Cell */
		longcode = MASHNONCOMM(longcode, ind0);
		StackMarkers[ind0] = 0;
		
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
//			HitClsInd = 0;
//			j1 = sg->v[lab[ind0]];
//			iend1 = j1+sg->d[lab[ind0]];
//			
//			for (; j1 < iend1; ++j1) {
//				k = sg->e[j1];
//				value = Part->inv[InvLab[k]];
//				if (cls[value] > 1) {
//					if (Markers[value] != tv->mark) {
//						HitCls[HitClsInd++] = value;
//						Markers[value] = tv->mark;
//						ElmHitCll[value] = value;
//					}
//					HitVtx[ElmHitCll[value]++] = k;
//				}
//				else {
//					longcode = MASHCOMM(longcode, value);
//				}
//			}
//			

			NEIGHCOUNT_SING;
            tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
			switch (SplInd) {
				case 0:
				case 1:
					break;
				case 2:
					if (SplCls[0] > SplCls[1]) {
						value = SplCls[0];
						SplCls[0] = SplCls[1];
						SplCls[1] = value;
					}
					break;
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
					for (k = 1; k < SplInd; ++k) {
						value = SplCls[k];
						i = k - 1;
						while ((i >= 0) && (value < SplCls[i])) {
							SplCls[i + 1] = SplCls[i];
							--i;
						}
						SplCls[i + 1] = value;
					}
					break;
				default:
					quickSort(SplCls, SplInd);
					break;
			}
			
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				i = ind1+cls[ind1]-ElmHitCll[ind1];
				trieref = trie_comp(trieref, i);
				if (trieref == NULL) return 0;
			}
			
			/* REARRANGE THE CELL */
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				cls[ind1] = cls[ind1]-ElmHitCll[ind1];
				newcell = ind1+cls[ind1];
				cls[newcell] = ElmHitCll[ind1];
				
				Part->cells++;
				if (StackMarkers[ind1] != tv->stackmark) {
					if (cls[newcell] < cls[ind1]) {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					else {
						CStack[++CStackInd] = ind1;
						StackMarkers[ind1] = tv->stackmark;
					}
				}
				else {
					CStack[++CStackInd] = newcell;
					StackMarkers[newcell] = tv->stackmark;
				}
				SplitCell = HitVtx+ind1;
				ind3 = cls[newcell];
				LabCell = lab+newcell;
				
				for (i1 = 0; i1 < ind3; i1++) {
					k = SplitCell[i1];
					i = LabCell[i1];
					Part->inv[newcell+i1] = newcell;
					lab[InvLab[k]] = i;
					InvLab[i] = InvLab[k];
					LabCell[i1] = k;
					InvLab[k] = newcell+i1;
				}
			}
		}
		else {
			if (ti->thegraphisparse) {
//				if (cls[ind0] == n) {
//					memcpy(NghCounts, sg->d, n*sizeof(int));
//					HitCls[0] = 0;
//					HitClsInd = 1;
//					ElmHitCll[0] = 0;
//					for (i = 0; i < n; i++) {
//						if (sg->d[i]) {
//							HitVtx[ElmHitCll[0]++] = i;
//						}
//					}
//				}
//				else {
//					HitClsInd = 0;
//                    for (i = ind0; i < ind2; i++) {
//                        labi = lab[i];
//                        nghb = sg->e+sg->v[labi];
//                        j = sg->d[labi];
//                        //for (; j--;) {
//                        while (j-- > 0) {
//                            k = nghb[j];
//                            if (MarkHitVtx[k] == tv->mark) {
//                                NghCounts[k]++;
//                            }
//                            else {
//                                value = Part->inv[InvLab[k]];
//                                if (cls[value] > 1) {
//                                    MarkHitVtx[k] = tv->mark;
//                                    NghCounts[k] = 1;
//                                    if (Markers[value] != tv->mark) {
//                                        HitCls[HitClsInd++] = value;
//                                        Markers[value] = tv->mark;
//                                        HitVtx[value] = k;
//                                        ElmHitCll[value] = 1;
//                                    }
//                                    else {
//                                        HitVtx[value+ElmHitCll[value]++] = k;
//                                    }
//                                }
//                                else {
//                                    longcode = MASHCOMM(longcode, value);
//                                }
//                            }
//						}
//					}
//                    
//                }
				
                NEIGHCOUNT_SPARSE;
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = ind1;
					}
					else {
						ind2 = ind1+cls[ind1];
						Split = FALSE;
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								Split = TRUE;
								break;
							}
						}
						if (!Split) {
							longcode = MASHCOMM(longcode, ind1);
						}
					}
				}
				/* Sorting the cells to be split */
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
					ind0 = SplCls[sc];
					ind1 = ind0 + cls[ind0];
					SplCntInd = 0;
					if (ElmHitCll[ind0] < cls[ind0]) {
						SplCnt[SplCntInd++] = 0;
						SplPos[0] = cls[ind0] - ElmHitCll[ind0];
					}
					
					/* According to the numbers of neighbors of C into the current cell */
					/* compute how many vertices in C will be placed into the same new cell */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = NghCounts[HitVtx[i]];
						if (Markers[value] != tv->mark) {
							Markers[value] = tv->mark;
							SplCnt[SplCntInd++] = value;
							SplPos[value] = 1;
						}
						else {
							SplPos[value]++;
						}
					}
					
					tv->mark++;
					
					/* Sort the values deriving from the previous step */
					switch (SplCntInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCnt[0] > SplCnt[1]) {
								value = SplCnt[0];
								SplCnt[0] = SplCnt[1];
								SplCnt[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplCntInd; ++k) {
								value = SplCnt[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCnt[i])) {
									SplCnt[i + 1] = SplCnt[i];
									--i;
								}
								SplCnt[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCnt, SplCntInd);
							break;
					}
					Part->cells += SplCntInd-1;
					
					/* Split the cell C and update the information for sizes of new cells */
					/* Put the new cells into the stack */
					i = ind0;
					if (StackMarkers[i] != tv->stackmark) {
						BigCellSize = 0;
					}
					for (k = 0; k < SplCntInd; k++) {
						value = SplPos[SplCnt[k]];
						cls[i] = value;
						if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
							BigCell = i;
							BigCellPos = CStackInd;
							BigCellSize = cls[i];
						}
						SplPos[SplCnt[k]] = i;
						i += value;
						if (i < ind1) {
							CStack[++CStackInd] = i;
							StackMarkers[i] = tv->stackmark;
							trieref = trie_comp(trieref, i);
							if (trieref == NULL) return 0;
						}
					}
					
					if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
						CStack[BigCellPos] = ind0;
						StackMarkers[BigCell] = 0;
						StackMarkers[ind0] = tv->stackmark;
					}
					/* Permute elements of the cell C */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = HitVtx[i];
						j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
						k = InvLab[value];				/* where HitVtx[i] is in lab */
						lab[k] = lab[j];
						lab[j] = value;
						InvLab[value] = j;
						InvLab[lab[k]] = k;
						NghCounts[value] = 0;
					}
					
					/* Reconstruct the cell C and update the inverse partition */
					newcell = ind1 - ElmHitCll[ind0];
					i = newcell;
					ind2 = newcell+cls[newcell]-1;
					do {
						Part->inv[i] = newcell;
						if (i == ind2) {
							newcell = i+1;
							if (newcell < n) ind2 = newcell+cls[newcell]-1;
						}
					}
					while (++i < ind1);
				}
			}
			else {
//				if (sg->d[lab[ind0]] > n/cls[ind0]) {
                if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
//						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
//						HitClsInd = 0;
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//								value = Part->inv[InvLab[k]];
//								if (Markers[value] != tv->mark) {
//									if (cls[value] > 1) HitCls[HitClsInd++] = value;
//									Markers[value] = tv->mark;
//								}
//							}
//						}
                        NEIGHCOUNT_DENSE_SPARSE;
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
						ind0 = SplCls[j];
						ind1 = ind0+cls[ind0];
						SplCntInd = 0;
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						for (i = ind0; i < ind1; i++) {
							value = NghCounts[lab[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						
						tv->mark++;
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
								trieref = trie_comp(trieref, i);
								if (trieref == NULL) return 0;
							}
						}
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						
						/* Permute elements of the cell C */
						i = ind0;
						do {
							SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
						}
						while(++i < ind1);
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind0;
						i = ind0;
						ind2 = newcell+cls[newcell]-1;
						do {
							lab[i] = SplCnt[i];
							InvLab[lab[i]] = i;
							Part->inv[i] = newcell;
							
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
					}
				}
				else {
					if (cls[ind0] == n) {
//						memcpy(NghCounts, sg->d, n*sizeof(int));
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//							}
//						}
                        NEIGHCOUNT_DENSE_DENSE;
					}
					
					ind0 = 0;
					while (ind0 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind0+cls[ind0];
						if (cls[ind0] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind0]];
							for (i = ind0+1; i < ind1; i++)
							{
								if (NghCounts[lab[i]] != value)
								{
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind0;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							
							if (SplCntInd) {
								tv->mark++;
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind0;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
									}
								}
								if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
									CStack[BigCellPos] = ind0;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind0] = tv->stackmark;
									trieref = trie_comp(trieref, i);
									if (trieref == NULL) return 0;
								}
								
								/* Permute elements of the cell C */
								i = ind0;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind0;
								i = ind0;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
							}
						}
						ind0 = ind1;
					}
				}
			}
		}
	}
	
	Cand->code = CLEANUP(longcode);
	return 1;
}

int traces_refine_sametrace(Candidate *Cand, 
							//int m, 
							int n, 
							Partition *Part, 
							struct TracesVars* tv, 
							struct TracesInfo *ti) {
	int i, j, k, sc, ind0, ind1, ind2, ind3, ind4, jk, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd;
//	size_t j1, iend1;
	int j1int, iend1int;
	unsigned int longcode;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, *TraceEnd, Traceccend, *Tracestpend;
	int BigCell, BigCellPos, BigCellSize;
	boolean TraceCell = FALSE;
    int *nghb;
    const int variation = 0;
    
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	SpineTL = Spine+tv->tolevel;
	TraceEnd = &(SpineTL->trcend);
	Traceccend = SpineTL->ccend;
	Tracestpend = &(SpineTL->stpend);
	TraceCCInd = SpineTL->ccstart;
	TraceStepsInd = SpineTL->stpstart;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	memcpy(CStack+1, TheTrace+SpineTL->trcstart, (Part->active)*sizeof(int));
	CStackInd = Part->active;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	TraceInd = SpineTL->trcstart+Part->active;
	
	while (CStackInd > 0)
	{
        
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
        //		if (Part->cells == tv->finalnumcells) break;    ?????
        if (Part->cells == n) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		StackMarkers[ind0] = 0;
        
        
		TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
        
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
//			HitClsInd = 0;
//			j1 = sg->v[lab[ind0]];
//			iend1 = j1+sg->d[lab[ind0]];
//			
//			for (; j1 < iend1; ++j1) {
//				k = sg->e[j1];
//				value = Part->inv[InvLab[k]];
//				if (cls[value] > 1) {
//					if (Markers[value] != tv->mark) {
//						HitCls[HitClsInd++] = value;
//						Markers[value] = tv->mark;
//						ElmHitCll[value] = value;
//					}
//					HitVtx[ElmHitCll[value]++] = k;
//				}
//			}

			NEIGHCOUNT_SING;
            tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
            
            /* SINGLETON CC CASE */
			if (SplInd) {
				SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
                
				
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					i = ind1+cls[ind1]-ElmHitCll[ind1];
                    SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
				}
				
				/* REARRANGE THE CELLS */
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					cls[ind1] = cls[ind1]-ElmHitCll[ind1];
					newcell = ind1+cls[ind1];
					cls[newcell] = ElmHitCll[ind1];
					Part->cells++;
					
					if (StackMarkers[ind1] != tv->stackmark) {
						if (cls[newcell] < cls[ind1]) {
							CStack[++CStackInd] = newcell;
							StackMarkers[newcell] = tv->stackmark;
						}
						else {
							CStack[++CStackInd] = ind1;
							StackMarkers[ind1] = tv->stackmark;
						}
					}
					else {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					
					SplitCell = HitVtx+ind1;
					ind3 = cls[newcell];
					LabCell = lab+newcell;
					
					for (jk = 0; jk < ind3; jk++) {
						k = SplitCell[jk];
						i = LabCell[jk];
						Part->inv[newcell+jk] = newcell;
						lab[InvLab[k]] = i;
						InvLab[i] = InvLab[k];
						LabCell[jk] = k;
						InvLab[k] = newcell+jk;
					}
					if (cls[ind1] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[ind1]);
                    }
					if (cls[newcell] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[newcell]);
					}
				}
			}
			else {
				if (TraceCell) {
                    return FALSE;
				}
			}
			
		}
		else {
			if (ti->thegraphisparse) {
//				if (cls[ind0] == n) {
//					memcpy(NghCounts, sg->d, n*sizeof(int));
//					HitCls[0] = 0;
//					HitClsInd = 1;
//					ElmHitCll[0] = 0;
//					for (i = 0; i < n; i++) {
//						if (sg->d[i]) {
//							HitVtx[ElmHitCll[0]++] = i;
//						}
//					}
//				}
//				else {
//					HitClsInd = 0;
//                    for (i = ind0; i < ind2; i++) {
//                        labi = lab[i];
//                        nghb = sg->e+sg->v[labi];
//                        j = sg->d[labi];
//                        //for (; j--;) {
//                        while (j-- > 0) {
//                            k = nghb[j];
//                            if (MarkHitVtx[k] == tv->mark) {
//                                NghCounts[k]++;
//                            }
//                            else {
//                                value = Part->inv[InvLab[k]];
//                                if (cls[value] > 1) {
//                                    MarkHitVtx[k] = tv->mark;
//                                    NghCounts[k] = 1;
//                                    if (Markers[value] != tv->mark) {
//                                        HitCls[HitClsInd++] = value;
//                                        Markers[value] = tv->mark;
//                                        HitVtx[value] = k;
//                                        ElmHitCll[value] = 1;
//                                    }
//                                    else {
//                                        HitVtx[value+ElmHitCll[value]++] = k;
//                                    }
//                                }
//                            }
//						}
//					}
//                    
//                }
				
                NEIGHCOUNT_SPARSE;
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = ind1;
					}
					else {
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
				}
				
                /* SPARSE CASE */
				if (SplInd) {
                    SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+n, &Traceccend)
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
						ind0 = SplCls[sc];
						ind1 = ind0 + cls[ind0];
						SplCntInd = 0;
						if (ElmHitCll[ind0] < cls[ind0]) {
							SplCnt[SplCntInd++] = 0;
							SplPos[0] = cls[ind0] - ElmHitCll[ind0];
						}
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = NghCounts[HitVtx[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						tv->mark++;
						
						if (SplCntInd) {
                            SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+n, Tracestpend)
						}
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
                                SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
							}
						}
						
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						/* Permute elements of the cell C */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = HitVtx[i];
							j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
							k = InvLab[value];				/* where HitVtx[i] is in lab */
							lab[k] = lab[j];
							lab[j] = value;
							InvLab[value] = j;
							InvLab[lab[k]] = k;
							NghCounts[value] = 0;
						}
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind1 - ElmHitCll[ind0];
						i = newcell;
						ind2 = newcell+cls[newcell]-1;
						do {
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
						
						for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
							if (cls[i] == 1) {
								Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
							}
						}
						
					}
				}
				else {
                    if (TraceCell) {
                        return FALSE;
                    }
				}
				
			}
			else {
//				if (sg->d[lab[ind0]] > n/cls[ind0]) {
                if (TheGraph[lab[ind0]].d > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
//						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
//						HitClsInd = 0;
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//								value = Part->inv[InvLab[k]];
//								if (Markers[value] != tv->mark) {
//									if (cls[value] > 1) HitCls[HitClsInd++] = value;
//									Markers[value] = tv->mark;
//								}
//							}
//						}
                        NEIGHCOUNT_DENSE_SPARSE;
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
                    /* DENSE-SPARSE CASE */
					if (SplInd) {
                        SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+2*n, &Traceccend)
						
						/* Sorting the cells to be split */
						switch (SplInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCls[0] > SplCls[1]) {
									value = SplCls[0];
									SplCls[0] = SplCls[1];
									SplCls[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplInd; ++k) {
									value = SplCls[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCls[i])) {
										SplCls[i + 1] = SplCls[i];
										--i;
									}
									SplCls[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCls, SplInd);
								break;
						}
						
						for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
							ind0 = SplCls[j];
							ind1 = ind0+cls[ind0];
							SplCntInd = 0;
							
							/* According to the numbers of neighbors of C into the current cell */
							/* compute how many vertices in C will be placed into the same new cell */
							for (i = ind0; i < ind1; i++) {
								value = NghCounts[lab[i]];
								if (Markers[value] != tv->mark) {
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = 1;
								}
								else {
									SplPos[value]++;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
                                SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+2*n, Tracestpend)
							}
							
							/* Sort the values deriving from the previous step */
							switch (SplCntInd) {
								case 0:
								case 1:
									break;
								case 2:
									if (SplCnt[0] > SplCnt[1]) {
										value = SplCnt[0];
										SplCnt[0] = SplCnt[1];
										SplCnt[1] = value;
									}
									break;
								case 3:
								case 4:
								case 5:
								case 6:
								case 7:
								case 8:
									for (k = 1; k < SplCntInd; ++k) {
										value = SplCnt[k];
										i = k - 1;
										while ((i >= 0) && (value < SplCnt[i])) {
											SplCnt[i + 1] = SplCnt[i];
											--i;
										}
										SplCnt[i + 1] = value;
									}
									break;
								default:
									quickSort(SplCnt, SplCntInd);
									break;
							}
							
							Part->cells += SplCntInd-1;
							
							/* Split the cell C and update the information for sizes of new cells */
							/* Put the new cells into the stack */
							i = ind0;
							if (StackMarkers[i] != tv->stackmark) {
								BigCellSize = 0;
							}
							for (k = 0; k < SplCntInd; k++) {
								value = SplPos[SplCnt[k]];
								cls[i] = value;
								if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
									BigCell = i;
									BigCellPos = CStackInd;
									BigCellSize = cls[i];
								}
								SplPos[SplCnt[k]] = i;
								i += value;
								if (i < ind1) {
									CStack[++CStackInd] = i;
									StackMarkers[i] = tv->stackmark;
                                    SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
								}
							}
							
							if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
								CStack[BigCellPos] = ind0;
								StackMarkers[BigCell] = 0;
								StackMarkers[ind0] = tv->stackmark;
							}
							
							/* Permute elements of the cell C */
							i = ind0;
							do {
								SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
							}
							while(++i < ind1);
							
							/* Reconstruct the cell C and update the inverse partition */
							newcell = ind0;
							i = ind0;
							ind2 = newcell+cls[newcell]-1;
							do {
								lab[i] = SplCnt[i];
								InvLab[lab[i]] = i;
								Part->inv[i] = newcell;
								if (i == ind2) {
									newcell = i+1;
									if (newcell < n) ind2 = newcell+cls[newcell]-1;
								}
							}
							while (++i < ind1);
							
							for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
								if (cls[i] == 1) {
									Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
								}
							}
							
						}
					}
					else {
                        if (TraceCell) {
                            return FALSE;
                        }
					}
					
				}
				else {
					if (cls[ind0] == n) {
//						memcpy(NghCounts, sg->d, n*sizeof(int));
                        for (i = 0; i < n; i++) {
                            NghCounts[i] = TheGraph[i].d;
                        }
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
//						for (i = ind0; i < ind2; i++) {
//							j1 = sg->v[lab[i]];
//							iend1 = j1+sg->d[lab[i]];
//							for (; j1 < iend1; ++j1) {
//								k = sg->e[j1];
//								(NghCounts[k])++;
//							}
//						}
                        
                        NEIGHCOUNT_DENSE_DENSE;
					}
					SplInd = 0;
					ind4 = 0;
					while (ind4 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind4+cls[ind4];
						if (cls[ind4] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind4]];
							for (i = ind4+1; i < ind1; i++) {
								if (NghCounts[lab[i]] != value)
								{
									SplInd++;
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind4;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
                                SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+3*n, Tracestpend)
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind4;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
                                        SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
									}
								}
								if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
									CStack[BigCellPos] = ind4;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind4] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind4;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind4;
								i = ind4;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
								
								for (i = ind4, k = 0; k < SplCntInd; i+=cls[i], k++) {
									if (cls[i] == 1) {
										Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
									}
								}
								
							}
						}
						ind4 = ind1;
					}
					
                    /* DENSE-DENSE CASE */
					if (SplInd) {
                        SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+3*n, &Traceccend)
                    }
					else {
                        if (TraceCell) {
                            return FALSE;
                        }
					}
                    
				}
			}
		}
	}  /* end while (CStackInd > 0) */
    
//    for (i=SpineTL->trcstart; i < TraceInd; i++) {
//		ind0 = TheTrace[i];
//		longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
//        j1 = sg->v[lab[ind0]];
//		iend1 = j1+sg->d[lab[ind0]];
//		for (; j1 < iend1; ++j1) {
//			value = Part->inv[InvLab[sg->e[j1]]];
//			longcode = MASHCOMM(longcode, value);
//        }
//	}
    for (i=SpineTL->trcstart; i < TraceInd; i++) {
        ind0 = TheTrace[i];
        longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
        labi = lab[ind0];
        iend1int = TheGraph[labi].d;
        nghb = TheGraph[labi].e;
        for (j1int = 0; j1int < iend1int; ++j1int) {
            k = nghb[j1int];
            value = Part->inv[InvLab[k]];
            longcode = MASHCOMM(longcode, value);
        }
    }

    Part->code = Cand->code = CLEANUP(longcode);
    if ((Cand->code != SpineTL->part->code) || (TraceInd != *TraceEnd)) return FALSE;
	return TRUE;
}

void refine_tr(sparsegraph *sg, int *lab, int *ptn, int *numcells, int *code, TracesOptions *options) {

	int i, j;
	int TraceEnd;
	
	struct TracesVars tvar, *tv;
	struct TracesInfo tinf, *ti;
	
    Partition *CurrPart;
    Candidate *CurrCand;
	
	const int n = sg->nv;
//	const int m = SETWORDSNEEDED(n);
	
    if (n > (NAUTY_INFINITY-2))
    {
        fprintf(ERRFILE, "Traces: need n <= %d, but n=%d\n\n", 
                NAUTY_INFINITY-2, n);
        return;
    }
    
#if !MAXN
	DYNALLOC1(int, Markers, Markers_sz, n, "refine_tr");
	DYNALLOC1(int, MarkHitVtx, MarkHitVtx_sz, n, "refine_tr");
	DYNALLOC1(int, NghCounts, NghCounts_sz, n, "refine_tr");
	DYNALLOC1(int, SplPos, SplPos_sz, n, "refine_tr");
	DYNALLOC1(int, SplCls, SplCls_sz, n, "refine_tr");
	DYNALLOC1(int, SplCnt, SplCnt_sz, n, "refine_tr");
	DYNALLOC1(int, StackMarkers, StackMarkers_sz, n, "refine_tr");
	DYNALLOC1(int, TheTrace, TheTrace_sz, n+10, "refine_tr");
	DYNALLOC1(int, TheTraceSteps, TheTraceSteps_sz, n+10, "refine_tr");
	DYNALLOC1(int, TheTraceCC, TheTraceCC_sz, n, "refine_tr");
	DYNALLOC1(int, TheTraceSplNum, TheTraceSplNum_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray2, WorkArray2_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray3, WorkArray3_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray4, WorkArray4_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray5, WorkArray5_sz, n, "refine_tr");
	DYNALLOC1(TracesSpine, Spine, Spine_sz, n, "refine_tr");
#endif
	
#define HitCls WorkArray2
#define HitVtx WorkArray3
#define CStack WorkArray4
#define ElmHitCll WorkArray5
	
	outfile = (options->outfile == NULL ? stdout : options->outfile);
	
	tv = &tvar;
	ti = &tinf;
	/* The graph is sparse? */
	if (sg->nde < n || sg->nde / n < n / (sg->nde / n)) {
		ti->thegraphisparse = TRUE;
	}
	else {
		ti->thegraphisparse = FALSE;
	}
    
	/* Initialize candidate, partition, cells, orbits */
	CurrPart = malloc(sizeof(Partition));
	if (CurrPart == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrPart->cls = malloc(n*sizeof(int));
	if (CurrPart->cls == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrPart->inv = malloc(n*sizeof(int));
	if (CurrPart->inv == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrPart->code = -1;
	
	CurrCand = malloc(sizeof(Candidate));
	if (CurrCand == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrCand->lab = malloc(n*sizeof(*CurrCand->lab));
	if (CurrCand->lab == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrCand->invlab = malloc(n*sizeof(*CurrCand->invlab));
	if (CurrCand->invlab == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrCand->do_it = TRUE;
	CurrCand->indnum = 0;
	CurrCand->code = 0;
	CurrCand->next = NULL;
	
	memset(CurrPart->cls, 0, n*sizeof(int));
	memset(CurrPart->inv, 0, n*sizeof(int));
	memcpy(CurrCand->lab, lab, n*sizeof(int));
	memcpy(CurrCand->invlab, lab, n*sizeof(int));
	CurrPart->cells = 0;
	
	j = 0;
	for (i = 0; i < n; i++) {
		if (j) {
			CurrPart->inv[i] = j;
		}
		CurrCand->invlab[CurrCand->lab[i]] = i;
		if (!ptn[i]) {
			CurrPart->cls[j] = i-j+1;
			TheTrace[CurrPart->cells++] = j;
			j = i+1;
		}
	}
	
	/* Further initializations */
	tv->mark = tv->stackmark = 1073741825;
	tv->maxtreelevel = 0;
	tv->tolevel = 0;
	
	TraceEnd = traces_refine_refine(sg, CurrCand, n, CurrPart, tv, ti);
	
	for (i = CurrPart->active; i < TraceEnd; i++) {
		ptn[TheTrace[i]-1] = 0;
	}
	memcpy(lab, CurrCand->lab, n*sizeof(int));
	*code = CurrCand->code;
	*numcells = CurrPart->cells;
    
	FREECAND(CurrCand)
	FREEPART(CurrPart)
    
#if !MAXN
    DYNFREE(Markers, Markers_sz);
    DYNFREE(MarkHitVtx, MarkHitVtx_sz);
    DYNFREE(NghCounts, NghCounts_sz);
    DYNFREE(SplCls, SplCls_sz);
    DYNFREE(SplCnt, SplCnt_sz);
    DYNFREE(SplPos, SplPos_sz);
    DYNFREE(StackMarkers, StackMarkers_sz);
    DYNFREE(TheTrace, TheTrace_sz);
    DYNFREE(TheTraceCC, TheTraceCC_sz);
    DYNFREE(TheTraceSplNum, TheTraceSplNum_sz);
    DYNFREE(TheTraceSteps, TheTraceSteps_sz);
    DYNFREE(WorkArray2, WorkArray2_sz);
    DYNFREE(WorkArray3, WorkArray3_sz);
    DYNFREE(WorkArray4, WorkArray4_sz);
    DYNFREE(WorkArray5, WorkArray5_sz);
    DYNFREE(Spine, Spine_sz);
#endif
}

int traces_refine_refine(sparsegraph *sg, 
						 Candidate *Cand, 
						 //int m, 
						 int n, 
						 Partition *Part, 
						 struct TracesVars* tv, 
						 struct TracesInfo *ti) {
	int i, j, k, sc, ind0, ind1, ind2, ind3, ind4, jk, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd;
	size_t j1, iend1;
	unsigned int longcode;
	int newtrace = FALSE;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, TraceEnd, Traceccend, Tracestpend;
	int BigCell, BigCellPos, BigCellSize;
	boolean TraceCell = FALSE;
    int *nghb;
    
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	TraceEnd = Part->active = Part->cells;
	Traceccend = 0;
	Tracestpend = 0;
	TraceCCInd = 0;
	TraceStepsInd = 0;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	memcpy(CStack+1, TheTrace, (Part->active)*sizeof(int));
	CStackInd = Part->active;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	TraceInd = Part->active;
	
	newtrace = TRUE;
	
	while (CStackInd > 0)
	{
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
		if (Part->cells == n) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		StackMarkers[ind0] = 0;
        if (!newtrace) {
			TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
        }
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
            HitClsInd = 0;
			j1 = sg->v[lab[ind0]];
			iend1 = j1+sg->d[lab[ind0]];
			
			for (; j1 < iend1; ++j1) {
				k = sg->e[j1];
				value = Part->inv[InvLab[k]];
				if (cls[value] > 1) {
					if (Markers[value] != tv->mark) {
						HitCls[HitClsInd++] = value;
						Markers[value] = tv->mark;
						ElmHitCll[value] = value;
					}
					HitVtx[ElmHitCll[value]++] = k;
				}
			}
			tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
			if (SplInd) {
				if (newtrace) {
					TheTraceCC[TraceCCInd] = ind0;
                }
				TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
				
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					i = ind1+cls[ind1]-ElmHitCll[ind1];
					TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
				}
				
				/* REARRANGE THE CELLS */
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					cls[ind1] = cls[ind1]-ElmHitCll[ind1];
					newcell = ind1+cls[ind1];
					cls[newcell] = ElmHitCll[ind1];
					Part->cells++;
					
					if (StackMarkers[ind1] != tv->stackmark) {
						if (cls[newcell] < cls[ind1]) {
							CStack[++CStackInd] = newcell;
							StackMarkers[newcell] = tv->stackmark;
						}
						else {
							CStack[++CStackInd] = ind1;
							StackMarkers[ind1] = tv->stackmark;
						}
					}
					else {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					
					SplitCell = HitVtx+ind1;
					ind3 = cls[newcell];
					LabCell = lab+newcell;
					
					for (jk = 0; jk < ind3; jk++) {
						k = SplitCell[jk];
						i = LabCell[jk];
						Part->inv[newcell+jk] = newcell;
						lab[InvLab[k]] = i;
						InvLab[i] = InvLab[k];
						LabCell[jk] = k;
						InvLab[k] = newcell+jk;
					}
				}
			}
			else {
				if ((!newtrace) && TraceCell) {
					return 0;
				}
			}
			
		}
		else {
			if (ti->thegraphisparse) {
                if (cls[ind0] == n) {
					memcpy(NghCounts, sg->d, n*sizeof(int));
					HitCls[0] = 0;
					HitClsInd = 1;
					ElmHitCll[0] = 0;
					for (i = 0; i < n; i++) {
						if (sg->d[i]) {
							HitVtx[ElmHitCll[0]++] = i;
						}
					}
				}
                else {
					HitClsInd = 0;
                    for (i = ind0; i < ind2; i++) {
                        labi = lab[i];
                        nghb = sg->e+sg->v[labi];
                        j = sg->d[labi];
                        //for (; j--;) {
                        while (j-- > 0) {
                            k = nghb[j];
                            if (MarkHitVtx[k] == tv->mark) {
                                NghCounts[k]++;
                            }
                            else {
                                value = Part->inv[InvLab[k]];
                                if (cls[value] > 1) {
                                    MarkHitVtx[k] = tv->mark;
                                    NghCounts[k] = 1;
                                    if (Markers[value] != tv->mark) {
                                        HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                        HitVtx[value] = k;
                                        ElmHitCll[value] = 1;
                                    }
                                    else {
                                        HitVtx[value+ElmHitCll[value]++] = k;
                                    }
                                }
                            }
						}
					}
                }
                
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = ind1;
					}
					else {
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
				}
				
				if (SplInd) {
					if (newtrace) {
                        TheTraceCC[TraceCCInd] = ind0;
                    }
					TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
						ind0 = SplCls[sc];
						ind1 = ind0 + cls[ind0];
						SplCntInd = 0;
						if (ElmHitCll[ind0] < cls[ind0]) {
							SplCnt[SplCntInd++] = 0;
							SplPos[0] = cls[ind0] - ElmHitCll[ind0];
						}
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = NghCounts[HitVtx[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						tv->mark++;
						
						if (SplCntInd) {
							TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd, &Tracestpend)
						}
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
								TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
							}
						}
						
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						/* Permute elements of the cell C */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = HitVtx[i];
							j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
							k = InvLab[value];				/* where HitVtx[i] is in lab */
							lab[k] = lab[j];
							lab[j] = value;
							InvLab[value] = j;
							InvLab[lab[k]] = k;
							NghCounts[value] = 0;
						}
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind1 - ElmHitCll[ind0];
						i = newcell;
						ind2 = newcell+cls[newcell]-1;
						do {
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
					}
				}
				else {
					if ((!newtrace) && TraceCell) {
						return 0;
					}
				}
				
			}
			else {
				if (sg->d[lab[ind0]] > n/cls[ind0]) {
					Sparse = FALSE;
                }
				else {
					Sparse = TRUE;
                }
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						HitClsInd = 0;
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
								value = Part->inv[InvLab[k]];
								if (Markers[value] != tv->mark) {
									if (cls[value] > 1) HitCls[HitClsInd++] = value;
									Markers[value] = tv->mark;
								}
							}
						}
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					if (SplInd) {
						if (newtrace) {
                            TheTraceCC[TraceCCInd] = ind0;
                        }
						TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
						
						/* Sorting the cells to be split */
						switch (SplInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCls[0] > SplCls[1]) {
									value = SplCls[0];
									SplCls[0] = SplCls[1];
									SplCls[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplInd; ++k) {
									value = SplCls[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCls[i])) {
										SplCls[i + 1] = SplCls[i];
										--i;
									}
									SplCls[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCls, SplInd);
								break;
						}
						
						for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
							ind0 = SplCls[j];
							ind1 = ind0+cls[ind0];
							SplCntInd = 0;
							
							/* According to the numbers of neighbors of C into the current cell */
							/* compute how many vertices in C will be placed into the same new cell */
							for (i = ind0; i < ind1; i++) {
								value = NghCounts[lab[i]];
								if (Markers[value] != tv->mark) {
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = 1;
								}
								else {
									SplPos[value]++;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
								TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd, &Tracestpend)
							}
							
							/* Sort the values deriving from the previous step */
							switch (SplCntInd) {
								case 0:
								case 1:
									break;
								case 2:
									if (SplCnt[0] > SplCnt[1]) {
										value = SplCnt[0];
										SplCnt[0] = SplCnt[1];
										SplCnt[1] = value;
									}
									break;
								case 3:
								case 4:
								case 5:
								case 6:
								case 7:
								case 8:
									for (k = 1; k < SplCntInd; ++k) {
										value = SplCnt[k];
										i = k - 1;
										while ((i >= 0) && (value < SplCnt[i])) {
											SplCnt[i + 1] = SplCnt[i];
											--i;
										}
										SplCnt[i + 1] = value;
									}
									break;
								default:
									quickSort(SplCnt, SplCntInd);
									break;
							}
							
							Part->cells += SplCntInd-1;
							
							/* Split the cell C and update the information for sizes of new cells */
							/* Put the new cells into the stack */
							i = ind0;
							if (StackMarkers[i] != tv->stackmark) {
								BigCellSize = 0;
							}
							for (k = 0; k < SplCntInd; k++) {
								value = SplPos[SplCnt[k]];
								cls[i] = value;
								if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
									BigCell = i;
									BigCellPos = CStackInd;
									BigCellSize = cls[i];
								}
								SplPos[SplCnt[k]] = i;
								i += value;
								if (i < ind1) {
									CStack[++CStackInd] = i;
									StackMarkers[i] = tv->stackmark;
									TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
								}
							}
							
							if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
								CStack[BigCellPos] = ind0;
								StackMarkers[BigCell] = 0;
								StackMarkers[ind0] = tv->stackmark;
							}
							
							/* Permute elements of the cell C */
							i = ind0;
							do {
								SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
							}
							while(++i < ind1);
							
							/* Reconstruct the cell C and update the inverse partition */
							newcell = ind0;
							i = ind0;
							ind2 = newcell+cls[newcell]-1;
							do {
								lab[i] = SplCnt[i];
								InvLab[lab[i]] = i;
								Part->inv[i] = newcell;
								if (i == ind2) {
									newcell = i+1;
									if (newcell < n) ind2 = newcell+cls[newcell]-1;
								}
							}
							while (++i < ind1);
						}
					}
					else {
						if ((!newtrace) && TraceCell) {
							return 0;
						}
					}
					
				}
				else {
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
							}
						}
					}
					SplInd = 0;
					ind4 = 0;
					while (ind4 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind4+cls[ind4];
						if (cls[ind4] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind4]];
							for (i = ind4+1; i < ind1; i++) {
								if (NghCounts[lab[i]] != value)
								{
									SplInd++;
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind4;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
								TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd, &Tracestpend)
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind4;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
										TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
									}
								}
								if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
									CStack[BigCellPos] = ind4;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind4] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind4;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind4;
								i = ind4;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
							}
						}
						ind4 = ind1;
					}
					
					if (SplInd) {
						if (newtrace) {
                            TheTraceCC[TraceCCInd] = ind0;
                        }
						TheTraceSplNum[TraceCCInd] = SplInd;
						TraceCCInd++;
					}
					else {
						if ((!newtrace) && TraceCell) {
							return 0;
						}
					}
					
				}
			}
		}
	} /* end while (CStackInd > 0) */
	
	for (i=0; i < TraceInd; i++) {
		ind0 = TheTrace[i];
		longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
		j1 = sg->v[lab[ind0]];
		iend1 = j1+sg->d[lab[ind0]];
		for (; j1 < iend1; ++j1) {
			value = Part->inv[InvLab[sg->e[j1]]];
			longcode = MASHCOMM(longcode, value);
		}
	}
	Part->code = Cand->code = CLEANUP(longcode);
	return TraceInd;
}

void quickSort(int *arr, int elements) {
	
#define MAX_LEVELS 300
	
	int piv, beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;
	int k, value;
	
	beg[0] = 0;
	end[0] = elements;
	while (i>= 0) {
		L = beg[i];
		R = end[i]-1;
		if (L<R-8) {
			piv = arr[(L+R)/2];
			arr[(L+R)/2] = arr[L];
			arr[L] = piv;
			while (L<R) {
				while (arr[R]>= piv && L<R) R--;
				if (L<R) arr[L++] = arr[R];
				while (arr[L]<= piv && L<R) L++;
				if (L<R) arr[R--] = arr[L];
			}
			arr[L] = piv;
			beg[i+1] = L+1;
			end[i+1] = end[i]; end[i++] = L;
			if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
				swap = beg[i];
				beg[i] = beg[i-1];
				beg[i-1] = swap;
				swap = end[i];
				end[i] = end[i-1];
				end[i-1] = swap;
			}
		}
		else {
			i--;
		}
	}
	for (k = 1; k < elements; ++k) {
		value = arr[k];
		i = k - 1;
		while ((i >= 0) && (value < arr[i])) {
			arr[i + 1] = arr[i];
			--i;
		}
		arr[i + 1] = value;
	}
}

struct Candidate *NewCandidate(int n, Candidate **GarbList, int Mrk) {
	struct Candidate *Cand;
	
	if (*GarbList) {
		Cand = *GarbList;
		*GarbList = (*GarbList)->next;
	}
	else {
		Cand = malloc(sizeof(*Cand));
		if (Cand == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
		Cand->lab = malloc(n*sizeof(*Cand->lab));
		if (Cand->lab == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
		Cand->invlab = malloc(n*sizeof(*Cand->invlab));
		if (Cand->invlab == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
	}
	Cand->do_it = Mrk;
	Cand->indnum = 0;
	Cand->code = 0;
	Cand->next = NULL;
	Cand->stnode = NULL;
	Cand->sortedlab = FALSE;
	return Cand;
}

int FreeList(Candidate *List, int cond) {
	Candidate *Temp;
	int conta = 0;
	int conta1 = 0;
	while (List) {
		if (List->do_it == cond) {
			conta1++;
		}
		conta++;
		Temp = List;
		if (List->lab) free(List->lab);
		if (List->invlab) free(List->invlab);
		List = List->next;
		free(Temp);
	}
	
	if (cond) {
		return conta1;
	}
	else {
		return conta;
	}
}

int FixBase(int *fix, struct TracesVars *tv, Candidate *Cand, int from, int to) {
	int i, j, k, go, nfix;
    nfix = j = 0;
	go = TRUE;
	for (i = from; i < to; i++) {
		k = Cand->lab[Spine[i+1].tgtpos];
		if (go && (nfix < tv->nfix) && (fix[nfix] == k)) {
			j++;
		}
		else {
			fix[nfix] = k;
			if (go) go = FALSE;
		}
		nfix++;
	}
	tv->nfix = nfix;
    return j;
}

void traces_freedyn(void) {
	/* Free the static dynamic memory used by Traces */
#if !MAXN
    DYNFREE(AUTPERM, AUTPERM_sz);
    DYNFREE(BreakSteps, BreakSteps_sz);
    DYNFREE(CurrOrbSize, CurrOrbSize_sz);
    DYNFREE(CurrRefCells, CurrRefCells_sz);
    DYNFREE(fix, fix_sz);
    DYNFREE(IDENTITY_PERM, IDENTITY_PERM_sz);
    DYNFREE(Markers, Markers_sz);
    DYNFREE(TreeMarkers, TreeMarkers_sz);
    DYNFREE(AutMarkers, AutMarkers_sz);
    DYNFREE(MarkHitVtx, MarkHitVtx_sz);
    DYNFREE(MultRefCells, MultRefCells_sz);
    DYNFREE(NghCounts, NghCounts_sz);
    DYNFREE(OrbSize, OrbSize_sz);
    DYNFREE(OrbList, OrbList_sz);
    DYNFREE(PrmPairs, PrmPairs_sz);
    DYNFREE(TempOrbList, TempOrbList_sz);
    DYNFREE(RefCells, RefCells_sz);
    DYNFREE(RefPath, RefPath_sz);
    DYNFREE(Spine, Spine_sz);
    DYNFREE(Singletons, Singletons_sz);
    DYNFREE(SplCls, SplCls_sz);
    DYNFREE(SplCnt, SplCnt_sz);
    DYNFREE(SplPos, SplPos_sz);
    DYNFREE(StackMarkers, StackMarkers_sz);
    DYNFREE(TheTrace, TheTrace_sz);
    DYNFREE(TheTraceCC, TheTraceCC_sz);
    DYNFREE(TheTraceSplNum, TheTraceSplNum_sz);
    DYNFREE(TheTraceSteps, TheTraceSteps_sz);
    DYNFREE(TreeStack, TreeStack_sz);
    DYNFREE(TrieArray, TrieArray_sz);
	DYNFREE(TEMPLAB, TEMPLAB_sz);
	DYNFREE(TEMPINVLAB, TEMPINVLAB_sz);
	DYNFREE(WorkArray, WorkArray_sz);
    DYNFREE(WorkArray0, WorkArray0_sz);
    DYNFREE(WorkArray1, WorkArray1_sz);
    DYNFREE(WorkArray2, WorkArray2_sz);
    DYNFREE(WorkArray3, WorkArray3_sz);
    DYNFREE(WorkArray4, WorkArray4_sz);
    DYNFREE(WorkArray5, WorkArray5_sz);
    DYNFREE(WorkArray6, WorkArray6_sz);
    DYNFREE(WorkArray7, WorkArray7_sz);
    DYNFREE(TheGraph, TheGraph_sz);
    DYNFREE(Neighbs1, Neighbs1_sz);
    DYNFREE(Neighbs2, Neighbs2_sz);
#endif
}

void factorial(double *size1, int *size2, int k) {
    int i;
//    writegroupsize(outfile, *size1, *size2);
//    printf(" * [%d]\n", k);
    for(i = k; i; i--) {
        MULTIPLY(*size1, *size2, i);
    }
    //    printf("grpsize: ");
    //    writegroupsize(outfile, *size1, *size2);
    //    printf("\n");
}

void factorial2(double *size1, int *size2, int k) {
    int i;
//    writegroupsize(outfile, *size1, *size2);
//    printf(" * (%d)\n", k);
    for(i = k; i > 0; i -= 2) {
        MULTIPLY(*size1, *size2, i);
    }
    //    printf("grpsize: ");
    //    writegroupsize(outfile, *size1, *size2);
    //    printf("\n");
}

int CheckForAutomorphisms(Candidate *CurrCand, Candidate *NextCand, 
                          struct TracesVars* tv, struct TracesInfo* ti, 
                          int m, int n, Partition* Part) {
	Candidate *CheckAutList;
	int i, j, k, tgt_level, numtemporbits; //, j1;
	int CheckLevel, CheckLevelEnd;
	int temp, tmp, tmp1, arg, arg1, val, val1; //, lev;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    
    CheckLevel = 0;
    CheckLevelEnd = 0;
    temp = 0;
    tv->gotonode = NULL;
    tv->conta6++;
    
	switch (tv->compstage) {
		case 0:
            if (tv->strategy) {
                CheckLevel = CheckLevelEnd = tv->maxtreelevel;
            }
            else {
                if ((Spine[tv->tolevel].part)->cells == tv->finalnumcells) {
                    CheckLevel = CheckLevelEnd = tv->maxtreelevel;
                }
                else {
                    CheckLevel = 1;
                    if ((Spine[tv->maxtreelevel].part)->cells == tv->finalnumcells) {
                        CheckLevelEnd = tv->maxtreelevel - 1;
                    }
                    else {
                        CheckLevelEnd = tv->maxtreelevel;
                    }
                }
            }
			break;
		case 1:
			CheckLevel = CheckLevelEnd = tv->tolevel;
			break;
		case 2:
			if (m || (tv->tolevel == tv->maxtreelevel+1)) {
				CheckLevel = CheckLevelEnd = tv->maxtreelevel+1;
			}
			else {
				CheckLevel = 1;
				if ((Spine[tv->maxtreelevel].part)->cells == tv->finalnumcells) {
					CheckLevelEnd = tv->maxtreelevel - 1;
				}
				else {
					CheckLevelEnd = tv->maxtreelevel;
				}
			}
			break;
		default:
			break;
	}

    while (CheckLevel <= CheckLevelEnd) {
		CheckAutList = Spine[CheckLevel].liststart;
		while (CheckAutList) {
			if (CheckAutList->do_it && lookup(CheckAutList->stnode) && (CheckAutList != NextCand) && (CheckAutList != CurrCand)) {
                if (CheckAutList->code == NextCand->code) {
                    SETMARK(Markers, tv->mark)
                    if (Part->cells == n) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        for (i = 0; i < n; i++) {
                            arg = NextCand->lab[i];
                            val = CheckAutList->lab[i];
                            SETPAIRSAUT(arg, val)
                        }
                    }
                    else {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        SETMARK(CellMarkers1, tv->markcell1)
                        SETMARK(CellMarkers2, tv->markcell2)
                        for (i=0; i<n; i+=Part->cls[i]) {
                            if (Part->cls[i] == 1) {
                                arg = NextCand->lab[i];
                                val = CheckAutList->lab[i];
                                SETPAIRSAUT(arg, val)
                                if ((TheGraph[arg].d > 1) && (tv->input_graph->d[arg] != TheGraph[arg].d))
                                    MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                            }
                            else {
                                k = i;
                                for (j=i; j<i+Part->cls[i]; j++) {
                                    arg = arg1 = NextCand->lab[j];
                                    if (CellMarkers1[arg] != tv->markcell1) {
                                        CellMarkers1[arg] = tv->markcell1;
                                        while ((CellMarkers2[CheckAutList->lab[k]] == tv->markcell2) && (k < i+Part->cls[i])) {
                                            k++;
                                        }
                                        if (k < i+Part->cls[i]) {
                                            val = val1 = CheckAutList->lab[k];
                                            CellMarkers2[val] = tv->markcell2;
                                            SETPAIRSAUT(arg, val)
                                            if ((TheGraph[arg].d > 1) && (tv->input_graph->d[arg] != TheGraph[arg].d))
                                                MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                            tmp = FirstNeighbour(arg, NextCand, Part, CellMarkers1, tv->markcell1, &arg, n);
                                            if (tmp) {
                                                CellMarkers1[arg] = tv->markcell1;
                                                tmp = FirstNeighbour(val, CheckAutList, Part, CellMarkers2, tv->markcell2, &val, n);
                                                CellMarkers2[val] = tv->markcell2;
                                                SETPAIRSAUT(arg, val)
                                                if ((TheGraph[arg].d > 1) && (tv->input_graph->d[arg] != TheGraph[arg].d))
                                                    MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                                while (tmp) {
                                                    tmp = NextNeighbour(arg, NextCand, Part, CellMarkers1, tv->markcell1, &arg, n);
                                                    if (tmp) {
                                                        CellMarkers1[arg] = tv->markcell1;
                                                        tmp = NextNeighbour(val, CheckAutList, Part, CellMarkers2, tv->markcell2, &val, n);
                                                        if (tmp) {
                                                            CellMarkers2[val] = tv->markcell2;
                                                            SETPAIRSAUT(arg, val)
                                                            if ((TheGraph[arg].d > 1) && (tv->input_graph->d[arg] != TheGraph[arg].d))
                                                                MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                                        }
                                                    }
                                                }
                                                arg = arg1;
                                                val = val1;
                                                do {
                                                    tmp = NextNeighbour(arg, NextCand, Part, CellMarkers1, tv->markcell1, &arg, n);
                                                    if (tmp) {
                                                        CellMarkers1[arg] = tv->markcell1;
                                                        tmp = NextNeighbour(val, CheckAutList, Part, CellMarkers2, tv->markcell2, &val, n);
                                                        if (tmp) {
                                                            CellMarkers2[val] = tv->markcell2;
                                                            SETPAIRSAUT(arg, val)
                                                            if ((TheGraph[arg].d > 1) && (tv->input_graph->d[arg] != TheGraph[arg].d))
                                                                MakeTree(arg, val, tv->input_graph, n, tv, TRUE);
                                                        }
                                                    }
                                                } while (tmp);
                                            }
                                        }
                                    }
                               }
                            }
                        }
                    }

                    if (isautom_sg_pair((graph*)tv->input_graph, AUTPERM, tv->options->digraph, m, n, tv)) {
                        if (!findperm(gensB, AUTPERM, n)) {
                            if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                            addgenerator(&gpB, &gensB, AUTPERM, n);
                            if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                            if (tv->options->verbosity >= 2) {
								fprintf(outfile, "[A (%d, %d)] ", CheckLevel, CheckAutList->name);
							}
							tv->stats->numgenerators++;
							orbjoin_sp_perm(tv->orbits, AUTPERM, OrbList, n, &tv->stats->numorbits);
							ti->thegrouphaschanged = TRUE;
							ti->identitygroup = FALSE;
							if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
								PRINT_RETURN
							}
							if (tv->options->writeautoms) {
								fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators);
								writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
							}
							if (tv->options->userautomproc) {
                                (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
							}
						}
						else {
							if (tv->options->verbosity >= 2) {
								fprintf(outfile, "[A* (%d, %d)] ", CheckLevel, CheckAutList->name);
							}
						}
                        TrieCandFrom = NULL;
                        TrieCheckFrom = CheckAutList->stnode;
                        if (CurrCand->stnode->level <= 1) {
                            tgt_level = CurrCand->stnode->level + 1;
                            while (TrieCheckFrom->level > tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if (tv->tolevel <= TrieCheckFrom->level) {
                                tgt_level = tv->tolevel;
                                while (TrieCheckFrom->level != tgt_level) {
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                            else {
                                TrieCandFrom = CurrCand->stnode;
                                tgt_level = TrieCheckFrom->level;
                                while (TrieCandFrom->level != tgt_level) {
                                    TrieCandFrom = TrieCandFrom->father;
                                }
                            }
                        }
                        if (TrieCandFrom) {
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                                TrieCandFrom = CurrCand->stnode;
                                TrieCheckFrom = TrieCheckFrom->father;
                                while (TrieCandFrom->father != TrieCheckFrom->father) {
                                    TrieCandFrom = TrieCandFrom->father;
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                        }
                        
                        while (TrieCheckFrom->goes_to) {
                            TrieCheckFrom = TrieCheckFrom->goes_to;
                        }
                        
						for (temp=1; temp<=tv->tolevel; temp++) { 
							if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
								break;
							}
						}

                        if (temp == tv->tolevel) {
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n,
                                                        PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            } else {
                                orbjoin(tv->currorbit, AUTPERM, n);
                            }
                        }
                        
                        switch (tv->compstage) {
                            case 0:
                                if (tv->strategy && (tv->steps == 1)) {
                                    RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                    }
                                    else {
                                        TrieCheckFrom->index++;
                                    }
                                    NextCand->do_it = FALSE;
                                }
                                else {
                                    if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                        CheckAutList->do_it = FALSE;
                                        if (TrieCandFrom) {
                                            TrieCandFrom->index += TrieCheckFrom->index;
                                            tv->newindex = 0;
                                            TrieCheckFrom->goes_to = TrieCandFrom;
                                        }
                                        else {
                                            if (CurrCand->stnode->level > 1) {
                                                tv->newgotonode = TrieCheckFrom;
                                                tv->newindex = TrieCheckFrom->index;
                                            }
                                            else {
                                                tv->newgotonode = NULL;
                                                tv->newindex = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (TrieCandFrom) {
                                            TrieCheckFrom->index += TrieCandFrom->index;
                                            TrieCandFrom->goes_to = TrieCheckFrom;
                                        }
                                        else {
                                            TrieCheckFrom->index++;
                                        }
                                        NextCand->do_it = FALSE;
                                    }
                                }
                                break;
                            case 1:
                                TrieCheckFrom->index ++;
                                tv->gotonode = TrieCheckFrom;
                                break;
                            case 2:
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                }
                                else {
                                    TrieCheckFrom->index++;
                                }
                                if (temp == tv->maxtreelevel) {
                                    tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                    if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                }
                                break;
                            default:
                                break;
                        }
						return temp;
					}
				}
			}
			CheckAutList = CheckAutList->next;
		}
		CheckLevel++;
	}
	return FALSE;
}


int CheckForSingAutomorphisms(Candidate *CurrCand, Partition *NextPart, Candidate *NextCand, 
							  struct TracesVars* tv, struct TracesInfo* ti, 
							  int m, int n) {
	int i, j, temp, tmp, tmp1, result, tgt_level, numtemporbits;
	TracesSpine *SpineTL;
	Candidate *CheckAutList;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    
	SpineTL = Spine+tv->tolevel;
	CheckAutList = SpineTL->liststart;
	tv->gotonode = NULL;
    
	result = 0;
	while (CheckAutList != NULL) {
        if (CheckAutList->do_it && (CheckAutList->stnode->father == CurrCand->stnode)) {
            if (CheckAutList->firstsingcode == NextCand->firstsingcode) {  // hash crea collisioni?
                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                if ((tv->tolevel == 1) && (Spine[0].part->cells == 1)) tmp = 2; else tmp = tv->tolevel;
                if (TreeFyTwo(tmp, CheckAutList, NextCand, NextPart, n, tv, ti)) {
                    if (isautom_sg((graph*)tv->input_graph, AUTPERM, tv->options->digraph, m, n)) {
                        if (!findperm(gensB, AUTPERM, n)) {
                            if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                            addgenerator(&gpB, &gensB, AUTPERM, n);
                            if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                            result = CheckAutList->name;
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n, 
                                                        PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            }
                            tv->stats->numgenerators++;
							orbjoin_sp_perm(tv->orbits, AUTPERM, OrbList, n, &tv->stats->numorbits);
                            
                            ti->thegrouphaschanged = TRUE;
                            ti->identitygroup = FALSE;
                            if (tv->options->verbosity >= 2) fprintf(outfile, "[a(%d)] ", CheckAutList->name);
                            if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
                                PRINT_RETURN
                            }
                            if (tv->options->writeautoms) {
                                fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators);
                                writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
                            }
                            if (tv->options->userautomproc) {
                                (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
                            }
                        }
                        else {
                            if (tv->options->verbosity >= 2) {
                                fprintf(outfile, "[a*]");
                            }
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n, 
                                                        PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            }
                            else {
                                orbjoin(tv->currorbit, AUTPERM, n);
                            }
                            result = -CheckAutList->name;
                        }
                        
                        TrieCandFrom = NULL;
                        TrieCheckFrom = CheckAutList->stnode;
                        if (CurrCand->stnode->level <= 1) {
                            tgt_level = CurrCand->stnode->level + 1;
                            while (TrieCheckFrom->level > tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if (tv->tolevel <= TrieCheckFrom->level) {
                                tgt_level = tv->tolevel;
                                while (TrieCheckFrom->level != tgt_level) {
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                            else {
                                TrieCandFrom = CurrCand->stnode;
                                tgt_level = TrieCheckFrom->level;
                                while (TrieCandFrom->level != tgt_level) {
                                    TrieCandFrom = TrieCandFrom->father;
                                }
                            }
                        }
                        if (TrieCandFrom) {
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                                TrieCandFrom = CurrCand->stnode;
                                TrieCheckFrom = TrieCheckFrom->father;
                                while (TrieCandFrom->father != TrieCheckFrom->father) {
                                    TrieCandFrom = TrieCandFrom->father;
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                        }
                        
                        while (TrieCheckFrom->goes_to) {
                            TrieCheckFrom = TrieCheckFrom->goes_to;
                        }
                        
                        for (temp=1; temp<=tv->tolevel; temp++) {
                            if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
                                break;
                            }
                        }
                        
                        switch (tv->compstage) {
                            case 0:
                                if (tv->strategy && (tv->steps == 1)) {
                                    RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                    }
                                    else {
                                        TrieCheckFrom->index++;
                                    }
                                    NextCand->do_it = FALSE;
                                }
                                else {
                                    if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                        CheckAutList->do_it = FALSE;
                                        if (TrieCandFrom) {
                                            TrieCandFrom->index += TrieCheckFrom->index;
                                            tv->newindex = 0;
                                            TrieCheckFrom->goes_to = TrieCandFrom;
                                        }
                                        else {
                                            if (CurrCand->stnode->level > 1) {
                                                tv->newgotonode = TrieCheckFrom;
                                                tv->newindex = TrieCheckFrom->index;
                                            }
                                            else {
                                                tv->newgotonode = NULL;
                                                tv->newindex = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (TrieCandFrom) {
                                            TrieCheckFrom->index += TrieCandFrom->index;
                                            TrieCandFrom->goes_to = TrieCheckFrom;
                                        }
                                        else {
                                            TrieCheckFrom->index++;
                                        }
                                        NextCand->do_it = FALSE;
                                    }
                                }
                                break;
                            case 1:
                                TrieCheckFrom->index++;
                                tv->gotonode = TrieCheckFrom;
                                break;
                            case 2:
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                }
                                else {
                                    TrieCheckFrom->index++;
                                }
                                if (temp == tv->maxtreelevel) {
                                    tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                    if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                }
                                break;
                            default:
                                break;
                        }
                        return result;
                    }
                }
            }
        }
		CheckAutList = CheckAutList->next;
	}
	return result;
}

int CheckForMatching(Candidate *CurrCand, Candidate *NextCand, Partition *Part, struct TracesVars* tv, struct TracesInfo* ti, int m, int n) {
	int i, j, vtx, vtx1, temp, tmp1, tgt_level, numtemporbits;
	TracesSpine *SpineTL;
	Candidate *CheckAutList;
    int *cls;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    boolean CodeVerify;
    
	SpineTL = Spine+tv->tolevel;
	CheckAutList = SpineTL->liststart;
	cls = Part->cls;
    numtemporbits = 0;
    tv->gotonode = NULL;
    
//    unsigned int samplecode, samplecode1;
    
	while (CheckAutList != NULL) {
        if (CheckAutList->do_it && (CheckAutList->singcode == NextCand->singcode)) {
//            PrintPartition(NextCand->lab, Part->cls, n, labelorg, 7514);
//            PrintPartition(CheckAutList->lab, Part->cls, n, labelorg, 7515);
//            printf("-------------\nSings: ");
            TrieCheckFrom = CheckAutList->stnode->father;
            TrieCandFrom = CurrCand->stnode;
            while (TrieCandFrom != TrieCheckFrom) {
                TrieCandFrom = TrieCandFrom->father;
                TrieCheckFrom = TrieCheckFrom->father;
            }
            
//            PrintVect(Singletons, Spine[TrieCheckFrom->level+1].singstart, SpineTL->singend,0);
//            samplecode = samplecode1 = 0;
            
            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
            SETMARK(Markers,tv->mark)
            CodeVerify = TRUE;
            for (i=Spine[TrieCheckFrom->level+1].singstart; i<SpineTL->singend; i++) {
                Markers[NextCand->lab[Singletons[i]]] = tv->mark;
            }
            for (i=Spine[TrieCheckFrom->level+1].singstart; i<SpineTL->singend; i++) {
                vtx1 = CheckAutList->lab[Singletons[i]];
                if (Markers[vtx1] != tv->mark) {
                    CodeVerify = FALSE;
                    break;
                }
                vtx = NextCand->lab[Singletons[i]];
//                samplecode = MASHCOMM(samplecode,vtx);
//                samplecode1 = MASHCOMM(samplecode1,vtx1);
//                printf("(%d > %u)(%d > %u); ",vtx,samplecode,vtx1,samplecode1);
                SETPAIRSAUT(vtx, vtx1)
                MakeTree(vtx, vtx1, tv->input_graph, n, tv, TRUE);
            }
//            printf("[%d %d]\n",samplecode,samplecode1);
            tv->conta7++;
            if (CodeVerify) {
                if (isautom_sg_pair((graph*)tv->input_graph, AUTPERM, tv->options->digraph, m, n, tv)) {
                    
                    if (!findperm(gensB, AUTPERM, n)) {
                        if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                        if (tv->options->generators) addpermutation(&gensB, AUTPERM, n);  // !!!!!!!!
//                        addgenerator(&gpB, &gensB, AUTPERM, n);
                        //addpermutation(&gensB, AUTPERM, n);
                        
                        if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                        
//                        printf("CheckAutList->stnode->father: %p, CurrCand->stnode: %p; ",CheckAutList->stnode->father,CurrCand->stnode);
                        if (CheckAutList->stnode->father == CurrCand->stnode) {
//                        if ((CheckAutList->stnode->father == CurrCand->stnode) && TempOrbits) {
//                            printf("TempOrbits: %p, tv->permInd: %d\n", TempOrbits, tv->permInd);
                            if (TempOrbits) {
                                if (tv->compstage == 0) {
                                    for (j=0; j<tv->permInd; j++) {
                                        orbjoin_sp_pair(TempOrbits, TempOrbList, n, PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                    }
                                }
                                else {
                                    orbjoin(TempOrbits, AUTPERM, n);
                                }
                            } else {
                                orbjoin(tv->currorbit, AUTPERM, n);
                            }
                        }

                        tv->stats->numgenerators++;
                        
                        for (j=0; j<tv->permInd; j++) {
                            orbjoin_sp_pair(tv->orbits, OrbList, n, PrmPairs[j].arg, PrmPairs[j].val, &tv->stats->numorbits);
                        }
                        
                        ti->thegrouphaschanged = TRUE;
                        if (tv->options->verbosity >= 2) fprintf(outfile, "[M(%d)] ", CheckAutList->name);
                        if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
                            PRINT_RETURN
                        }
                        if (tv->options->writeautoms) {
                            fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators);
//                            printf("\n");
//                            PrintVect(CheckAutList->lab,0,n,labelorg);
//                            PrintVect(NextCand->lab,0,n,labelorg);
//                            if (VerifyPerm(AUTPERM,n,7563)>0)
                                writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
                        }
                        if (tv->options->userautomproc) {
                            (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
                        }
                        
                    }
                    else {
                        if (tv->options->verbosity >= 2) {
                            fprintf(outfile, "[M*]");
                        }
                        if (TempOrbits) {
                            if (tv->compstage == 0) {
                                for (j=0; j<tv->permInd; j++) {
                                    orbjoin_sp_pair(TempOrbits, TempOrbList, n,
                                                    PrmPairs[j].arg, PrmPairs[j].val, &numtemporbits);
                                }
                            }
                            else {
                                orbjoin(TempOrbits, AUTPERM, n);
                            }
                        }
                        //                    result = -CheckAutList->name;
                    }
                    
                    TrieCandFrom = NULL;
                    TrieCheckFrom = CheckAutList->stnode;
                    if (CurrCand->stnode->level <= 1) {
                        tgt_level = CurrCand->stnode->level + 1;
                        while (TrieCheckFrom->level > tgt_level) {
                            TrieCheckFrom = TrieCheckFrom->father;
                        }
                    }
                    else {
                        if (tv->tolevel <= TrieCheckFrom->level) {
                            tgt_level = tv->tolevel;
                            while (TrieCheckFrom->level != tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            TrieCandFrom = CurrCand->stnode;
                            tgt_level = TrieCheckFrom->level;
                            while (TrieCandFrom->level != tgt_level) {
                                TrieCandFrom = TrieCandFrom->father;
                            }
                        }
                    }
                    if (TrieCandFrom) {
                        while (TrieCandFrom->father != TrieCheckFrom->father) {
                            TrieCandFrom = TrieCandFrom->father;
                            TrieCheckFrom = TrieCheckFrom->father;
                        }
                    }
                    else {
                        if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                            TrieCandFrom = CurrCand->stnode;
                            TrieCheckFrom = TrieCheckFrom->father;
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                    }
                    
                    while (TrieCheckFrom->goes_to) {
                        TrieCheckFrom = TrieCheckFrom->goes_to;
                    }
                    
                    for (temp=1; temp<=tv->tolevel; temp++) {
                        if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
                            break;
                        }
                    }
                    switch (tv->compstage) {
                        case 0:
                            if (tv->strategy && (tv->steps == 1)) {
                                RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                }
                                else {
                                    TrieCheckFrom->index++;
                                }
                                NextCand->do_it = FALSE;
                            }
                            else {
                                if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                    CheckAutList->do_it = FALSE;
                                    if (TrieCandFrom) {
                                        TrieCandFrom->index += TrieCheckFrom->index;
                                        tv->newindex = 0;
                                        TrieCheckFrom->goes_to = TrieCandFrom;
                                    }
                                    else {
                                        if (CurrCand->stnode->level > 1) {
                                            tv->newgotonode = TrieCheckFrom;
                                            tv->newindex = TrieCheckFrom->index;
                                        }
                                        else {
                                            tv->newgotonode = NULL;
                                            tv->newindex = 0;
                                        }
                                    }
                                }
                                else {
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                    }
                                    else {
                                        TrieCheckFrom->index++;
                                    }
                                    NextCand->do_it = FALSE;
                                }
                            }
                            break;
                        case 1:
                            TrieCheckFrom->index ++;
                            tv->gotonode = TrieCheckFrom;
                            break;
                        case 2:
                            if (TrieCandFrom) {
                                TrieCheckFrom->index += TrieCandFrom->index;
                                TrieCandFrom->goes_to = TrieCheckFrom;
                            }
                            else {
                                TrieCheckFrom->index++;
                            }
                            if (temp == tv->maxtreelevel) {
                                tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                            }
                            break;
                        default:
                            break;
                    }
                    return temp;
                }
            }
        }
        CheckAutList = CheckAutList->next;
	}
    return FALSE;
}

void Individualize(Partition *NextPart, Candidate *NextCand, int K, int Tc, int Cl, int Pos) {
	int i, j;
	
	NextCand->do_it = TRUE;
	if (NextPart->cls[Tc] > 1) {
		NextPart->cells = Cl+1;
		NextPart->active = 1;
		NextPart->cls[Tc]--;
		NextPart->cls[Pos] = 1;
    }
	NextPart->inv[Pos] = Pos;
	
	j = NextCand->lab[Pos];
	i = NextCand->invlab[K];
	NextCand->lab[Pos] = K;
	NextCand->invlab[K] = Pos;
	NextCand->lab[i] = j;
	NextCand->invlab[j] = i;
	return;
}

boolean TargetCell(Candidate *TargCand, Partition *Part, int n, struct TracesVars* tv, int Lv) {
	int TCell = -1, TCSize = 1;
	int i;
//    printf("TC; ");
	if (Lv < tv->tcellevel) {
		tv->tcell = Spine[Lv+1].tgtcell;
		return TRUE;
	}
	else {
		while (TCell < 0) {
            for (i = Spine[Lv].tgtcell; i < Spine[Lv].tgtend; i += Part->cls[i]) {
                if (Part->cls[i] > TCSize) {
                    if ((NonSingDeg(TargCand->lab[i], TargCand, Part) > 2) && (TheGraph[TargCand->lab[i]].d < n-1)) {
						TCSize = Part->cls[i];
						TCell = i;
					}
				}
			}
			Lv--;
			if ((Lv < 0) && (TCell < 0)) return FALSE;
		}
        tv->tcell = TCell;
		return TRUE;
	}
}

//boolean TargetCellFirstPath(Candidate *TargCand, Partition *Part, struct TracesVars* tv) {
//	int n, TCell, TCSize;
//	int Lv, i, Lev, vtx, vtx_d;
//    int loopstart, loopend;
//    
//    n = tv->input_graph->nv;
//    if (Part->cells == n) {
//		return 0;
//	}
//    Lev = tv->tolevel_tl;
//    Lv = tv->tolevel_tl;
//    TCell = -1;
//    TCSize = 1;
//    while (TCell < 0) {
//        
//        loopstart = Spine[Lv].tgtcell;
//        loopend = Spine[Lv].tgtend;
//        
//        i = loopstart;
//        
//        while (i < loopend) {
//            if (Part->cls[i] > TCSize) {
//                vtx = TargCand->lab[i];
//                vtx_d = TheGraph[vtx].d;
//                
//                if ((vtx_d > 2) && (vtx_d < n-1)) {
//                    if ((NonSingDeg(vtx, TargCand, Part) > 2) && (vtx_d < n-1)) {
//                        TCSize = Part->cls[i];
//                        TCell = i;
//                        if (TCSize == WorkArray[Lv]) {
//                            break;
//                        }
//                    }
//                }
//            }
//            i += Part->cls[i];
//        }
//        
//        Lv = Spine[Lv].tgtfrom;
//        if ((Lv < 0) && (TCell < 0)) {
//            tv->finalnumcells = Part->cells;
//            return FALSE;
//        }
//    }
//    tv->tcellexpath = tv->lastcell = TCell;
//    tv->tolevel_tl++;
//    Spine[tv->tolevel_tl].tgtfrom = tv->lastlev = Lv+1;
//    Spine[tv->tolevel_tl].tgtcell = tv->tcellexpath;
//    Spine[tv->tolevel_tl].tgtsize = WorkArray[Lv+1] = TCSize;
//    Spine[tv->tolevel_tl].tgtend = Spine[tv->tolevel_tl].tgtcell + TCSize;
//    Spine[tv->tolevel_tl].tgtpos = Spine[tv->tolevel_tl].tgtend - 1;
//    tv->tcellevel = tv->tolevel_tl;
//    
//    if (Lv+1 != Lev) {
//        BreakSteps[Lev] = ++tv->brkstpcount;
//        if (Spine[tv->tolevel].liststart) {
//            if (!Spine[tv->tolevel].liststart->firstsingcode) {
//                Spine[tv->tolevel].liststart->firstsingcode = Spine[tv->tolevel].liststart->pathsingcode;
//            }
//        }
//    }
//    
//    return TRUE;
//}

boolean TargetCellFirstPath(Candidate *TargCand, Partition *Part, struct TracesVars* tv) {
	int n, TCell, TCSize, TCell1, TCSize1;
	int Lv, i, Lev, vtx, vtx_d;
    int loopstart, loopend;
    boolean divided;
    
    n = tv->input_graph->nv;
    if (Part->cells == n) {
		return 0;
	}
    Lev = tv->tolevel_tl;
    Lv = tv->tolevel_tl;
    TCell = TCell1 = -1;
    TCSize = TCSize1 = 1;
    
//    printf("TCFP (lv: %d) cells: %d\n",Lv,Part->cells);
    //    printf("TCFP (lv: %d) lastlev: %d, lastcell: %d\n",Lv,tv->lastlev,tv->lastcell);
    
    while (TCell < 0) {
        loopstart = Part->inv[Spine[Lv].tgtcell];
        divided = FALSE;
        
        if (Lv == tv->lastlev) {
            loopstart = Part->inv[tv->lastcell];
            divided = TRUE;
        }
        
        i = loopstart;
        loopend = Spine[Lv].tgtend;
        
//        printf("TCFP search in [%d...%d]\n",loopstart,loopend);
        while (i < loopend) {
            //            if (Part->cls[i] >1) {
            //                printf("try %d (nonsingdeg: %d); cell: ",i,NonSingDeg(TargCand->lab[i], TargCand, Part));
            //                PrintVect(TargCand->lab+i,0,Part->cls[i],labelorg);
            //            }
            if (Part->cls[i] > TCSize) {
                vtx = TargCand->lab[i];
                vtx_d = TheGraph[vtx].d;
                if ((vtx_d > 2) && (vtx_d < n-1)) {
                    if ((NonSingDeg(vtx, TargCand, Part) > 2) && (vtx_d < n-1)) {
                        TCSize = Part->cls[i];
                        TCell = i;
                        if (TCSize == WorkArray[Lv]) {
                            break;
                        }
                    }
                }
            }
            i += Part->cls[i];
            if (divided && (i == loopend)) {
                i = loopstart = Spine[Lv].tgtcell;
                loopend = tv->lastcell;
                divided = FALSE;
                TCSize1 = TCSize;
                TCell1 = TCell;
                TCell = -1;
                TCSize = 1;
//                printf("TC(div) search in [%d...%d]\n",loopstart,loopend);
            }
        }
        
        if (TCSize1 > TCSize) {
            TCell = TCell1;
            TCSize = TCSize1;
        }
        
        if (TCell < 0) {
            if (Lv == 0) {
                tv->finalnumcells = Part->cells;
                return FALSE;
            } else {
//                printf("Lv: %d ",Lv);
                Lv = Spine[Lv].tgtfrom;
//                printf("--> %d\n",Lv);
            }
        }
        
        //        if ((Lv < 0) && (TCell < 0)) {
        //            tv->finalnumcells = Part->cells;
        //            return FALSE;
        //        }
    }
    tv->tcellexpath = tv->lastcell = TCell;
    tv->tolevel_tl++;
    //    Spine[tv->tolevel_tl].tgtfrom = tv->lastlev = Lv+1;   // XXXZZZ
    Spine[tv->tolevel_tl].tgtfrom = tv->lastlev = Lv;   // XXXZZZ
//    printf("Set Spine[%d].tgtfrom to %d\n",tv->tolevel_tl,Spine[tv->tolevel_tl].tgtfrom);
    Spine[tv->tolevel_tl].tgtcell = tv->tcellexpath;
    Spine[tv->tolevel_tl].tgtsize = WorkArray[Lv] = TCSize;
//    Spine[tv->tolevel_tl].tgtsize = WorkArray[Lv+1] = TCSize;
    Spine[tv->tolevel_tl].tgtend = Spine[tv->tolevel_tl].tgtcell + TCSize;
    Spine[tv->tolevel_tl].tgtpos = Spine[tv->tolevel_tl].tgtend - 1;
    tv->tcellevel = tv->tolevel_tl;
    //    printf("SET(a) tv->tcellevel to %d\n", tv->tcellevel);
//    printf("TCFP at %d is %d: ",Lev, tv->tcellexpath);
//    PrintVect(TargCand->lab+tv->tcellexpath,0,Part->cls[tv->tcellexpath],labelorg);
    
//    for (i=0; i<=tv->tcellevel; i++) {
//        printf("Targets[%d]. tgtfrom: %d, tgtcell: %d, tgtsize: %d, tgtend: %d, tgtpos: %d\n", i,
//               Spine[i].tgtfrom, Spine[i].tgtcell, Spine[i].tgtsize, Spine[i].tgtend, Spine[i].tgtpos);
//    }

    
    //if (Lv+1 != Lev) {
    if (Lv != Lev) {
        BreakSteps[Lev] = ++tv->brkstpcount;
        if (Spine[tv->tolevel].liststart) {
            if (!Spine[tv->tolevel].liststart->firstsingcode) {
                Spine[tv->tolevel].liststart->firstsingcode = Spine[tv->tolevel].liststart->pathsingcode;
            }
        }
    }
    return TRUE;
}

int TargetCellExpPath(Candidate *TargCand, Partition *Part, struct TracesVars* tv) {
//	int TCell = -1, TCSize = 1;
	int Lv, n; //, i;
    
//    printf("TCEP, tv->tcellevel: %d, tv->tolevel_tl: %d\n", tv->tcellevel, tv->tolevel_tl);
//    n = sg->nv;
    n = tv->input_graph->nv;
    if (Part->cells == n) {
		return 0;
	}
    
    Lv = tv->tolevel_tl+1;
    SpineTL_tl = Spine+Lv;
	if (tv->tolevel_tl < tv->tcellevel) {
        tv->tcellexpath = Part->inv[SpineTL_tl->tgtcell];
        tv->tolevel_tl++;
        if (Part->cls[tv->tcellexpath] == 1) {
            return TargetCellExpPath(TargCand, Part, tv);
        } else {
            return 1+((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
        }
    }

//    Lv = tv->tolevel_tl+1;
//    SpineTL_tl = Spine+Lv;
////    printf("a) tv->tolevel_tl: %d, tv->tcellevel: %d\n", tv->tolevel_tl, tv->tcellevel);
//	if (tv->tolevel_tl < tv->tcellevel) {
//        tv->tcellexpath = Part->inv[SpineTL_tl->tgtcell];
//        tv->tolevel_tl++;
//        return 1+((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
//    }

	else {
//        printf("TCE->TCF @ %d (%d cells)\n", tv->tolevel_tl, Part->cells);
        if (TargetCellFirstPath(TargCand, Part, tv)) {
            return 1+((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
        }
        else {
            return 0;
        }
        //printf("TCEP; ");
//        while (TCell < 0) {
//            Lv--;
//            for (i = Part->inv[Spine[Lv].tgtcell]; i < Spine[Lv].tgtend; i += Part->cls[i]) {
//                //if ((sg->d[TargCand->lab[i]] > 2) && (sg->d[TargCand->lab[i]] > 2)) {
//                //if (sg->d[TargCand->lab[i]] > 2) {
//                if (Part->cls[i] > TCSize) {
//                    if (NonSingDeg(sg, TargCand->lab[i], TargCand, Part) > 2) {
//						TCSize = Part->cls[i];
//						TCell = i;
//					}
//				}
//			}
//		}
//        
//        tv->tolevel_tl++;
//		tv->tcellevel = tv->tolevel_tl;
//        Spine[tv->tolevel_tl].tgtcell = tv->tcellexpath = TCell;
//        Spine[tv->tolevel_tl].tgtfrom = Lv+1;
//        Spine[tv->tolevel_tl].tgtsize = Part->cls[TCell];
//		Spine[tv->tolevel_tl].tgtend = TCell+Part->cls[TCell];
//		Spine[tv->tolevel_tl].tgtpos = Spine[tv->tolevel_tl].tgtend - 1;
        
//        return ((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
	}
}

boolean SelectNextLevel(int n, struct TracesVars *tv) {
//    printf("SelectNextLevel @ %d (smalldeglevel: %d)\n", tv->compstage, tv->smalldeglevel);
	switch (tv->compstage) {
		case 2:
//            PRINTF2("SelectNextLevel 1?(<n): tv->smalldeglevel: %d; ", tv->smalldeglevel);
//            PRINTF2("tv->maxtreelevel: %d\n", tv->maxtreelevel);
//			if (tv->smalldeglevel < n) {
//                tv->nextlevel = tv->smalldeglevel - 1;
//			}
//			else {
//				tv->nextlevel = tv->maxtreelevel;
//			}
            tv->nextlevel = tv->maxtreelevel;
			while (tv->nextlevel >=0) {
				if (Spine[tv->nextlevel].liststart) {
					break;
				}
				tv->nextlevel--;
			}
			if (tv->nextlevel < 0) {
				return FALSE;
			}
			break;
		default:
			switch (tv->strategy) {
				case 0:
					tv->nextlevel = tv->fromlevel;
//                    printf("set tv->nextlevel = %d\n", tv->nextlevel);
					while (!Spine[tv->nextlevel].liststart) {
						(tv->nextlevel)++;
					}
                    PRINTF2("SelectNextLevel 1?: finalnumcells: %d; ", tv->finalnumcells);
                    PRINTF2("Spine[tv->nextlevel].part->cells: %d; ", Spine[tv->nextlevel].part->cells);
                    PRINTF2("tv->maxtreelevel: %d; ", tv->maxtreelevel);
                    PRINTF2("tv->nextlevel: %d\n", tv->nextlevel);
					if ((Spine[tv->nextlevel].part->cells == tv->finalnumcells) || (tv->nextlevel > tv->maxtreelevel)) {
                        return FALSE;
					}
					break;
				case 1:
					tv->nextlevel = tv->maxtreelevel;
                    PRINTF2("SelectNextLevel 2?: finalnumcells: %d; ", tv->finalnumcells);
                    PRINTF2("Spine[tv->nextlevel].part->cells: %d; ", Spine[tv->nextlevel].part->cells);
                    if (Spine[tv->nextlevel].part->cells == tv->finalnumcells) {
						(tv->nextlevel)--;
					}
                    while (tv->nextlevel >= 0) {
						if (Spine[tv->nextlevel].liststart) {
                            break;
						}
						tv->nextlevel--;
                    }
					if (tv->nextlevel < 0) {
						return FALSE;
					}
					break;
				default:
					break;
			}
			break;
	}
//    printf("nextlevel is %d\n", tv->nextlevel);
    return TRUE;
}

boolean TreeFyTwo(int From, Candidate *Cand1, Candidate *Cand2, Partition *Part, int n, 
                  struct TracesVars* tv, struct TracesInfo *ti) {
	int i, i1, i2, j1, j2, k;
    int vtx1, vtx2, ngh1, ngh2, arg, val;
    int *tgtc1, *tgtc2, *adj1, *adj2;
	int iend;
    
//    printf("TF2; ");
    SETMARK(Markers, tv->mark)
    //SETMARK(TreeMarkers, tv->treemark)
	i2=0;
    
    //Neighbs1
    
    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
    i1 = Spine[From].tgtsize;
    tgtc1 = Cand1->lab+Spine[From].tgtcell;
    tgtc2 = Cand2->lab+Spine[From].tgtcell;
    for (i=0; i<i1; i++) {
        arg = tgtc1[i];
        val = tgtc2[i];
        if ((Markers[arg] != tv->mark) && (Markers[val] != tv->mark)) {
            SETPAIRSAUT(arg, val)
            SETPAIRSAUT(val, arg)
            Markers[arg] = tv->mark;
            Markers[val] = tv->mark;
        }
//        TreeMarkers[val] = tv->treemark;
    }
    
//    for (i=0; i<i1; i++) {
//        printf("[%d, %d]", PrmPairs[i].arg, PrmPairs[i].val);
//    }
//    printf("\n");
    
    //memcpy(TreeStack, Cand1->lab+Spine[From].tgtcell, i1*sizeof(int));
//    PrintVect(TreeStack, 0, i1, labelorg);
//    return FALSE;
//    for (i=0; i<i1; i++) {
//        printf("[[%d -- %d]]\n", Cand1->lab[TreeStack[i]]+labelorg, Cand2->lab[TreeStack[i]]+labelorg);
//        SETPAIRSAUT(Cand1->lab[TreeStack[i]], Cand2->lab[TreeStack[i]])
//        //AUTPERM[Cand1->lab[TreeStack[i]]] = Cand2->lab[TreeStack[i]];
//    }
    //while (i2 < i1) {
    while (i2 < tv->permInd) {
		vtx1 = PrmPairs[i2].arg;
		vtx2 = PrmPairs[i2++].val;
        adj1 = TheGraph[vtx1].e;
        adj2 = TheGraph[vtx2].e;
        iend = TheGraph[vtx1].d;
        j1 = j2 = 0;
        for (k=0; k < iend; k++) {
            ngh1 = adj1[k];
			if (Markers[ngh1] != tv->mark) {
                Neighbs1[j1++] = Cand1->invlab[ngh1];
            }
            ngh2 = adj2[k];
			if (Markers[ngh2] != tv->mark) {
                Neighbs2[j2++] = Cand2->invlab[ngh2];
            }
            
            //            printf("k: %d\n", k+labelorg);
            //            PrintVect(Cand1->lab+Part->inv[Cand1->invlab[k]], 0, Part->cls[Part->inv[Cand1->invlab[k]]], labelorg);
            //            PrintVect(Cand2->lab+Part->inv[Cand1->invlab[k]], 0, Part->cls[Part->inv[Cand1->invlab[k]]], labelorg);
//            if ((k != Cand2->lab[Cand1->invlab[k]]) && (AUTPERM[k] == k)) {
//                printf("(%d -- %d)\n", k+labelorg, Cand2->lab[Cand1->invlab[k]]+labelorg);
//                SETPAIRSAUT(k, Cand2->lab[Cand1->invlab[k]])
//                //AUTPERM[k] = Cand2->lab[Cand1->invlab[k]];
//                TreeStack[i1++] = Cand1->invlab[k];
//                //                printf("MOVE %d into stack\n", Cand1->invlab[k]);
//            }
            
		}
        
//        printf("%6d: [ ", vtx1);
//        for (k=0; k<j1; k++) {
//            printf("%d (%d); ", Neighbs1[k], Cand1->lab[Neighbs1[k]]+labelorg);
//        }
//        printf("]\n");
//        printf("%6d: [ ", vtx2);
//        for (k=0; k<j2; k++) {
//            printf("%d (%d); ", Neighbs2[k], Cand2->lab[Neighbs2[k]]+labelorg);
//        }
//        printf("]\n");
        
        k = tv->permInd; // solo per stampare dopo
        if (j1 == j2) {
            quickSort(Neighbs1, j1);
            quickSort(Neighbs2, j2);
            for (i=0; i<j1; i++) {
                arg = Cand1->lab[Neighbs1[i]];
                val = Cand2->lab[Neighbs2[i]];
                if ((Markers[arg] != tv->mark) && (Markers[val] != tv->mark)) {
                    SETPAIRSAUT(arg, val)
                    SETPAIRSAUT(val, arg)
                    Markers[arg] = tv->mark;
                    Markers[val] = tv->mark;
                }
            }

//            for (i=k; i<tv->permInd; i++) {
//                printf("[%d, %d]", PrmPairs[i].arg, PrmPairs[i].val);
//            }
//            printf("\n");
            
        }
    }
    return TRUE;

//    return (i1 != Spine[From].tgtsize);
}

void ExperimentalStep(Partition *NextPart, Candidate *NextCand, 
					  TracesVars *tv, TracesInfo *ti, int m, int n) {
	int i, iend, min, tmp;
    
	SpineTL_tl = Spine+tv->tolevel_tl;
	NextPart->active = 1;
	
	/* EXPERIMENTAL PATH INDIVIDUALIZATION AND REFINEMENT */
    if (tv->answ == 2) {
        min = NextCand->lab[tv->tcellexpath];
        tmp = tv->tcellexpath;
        iend = tv->tcellexpath + NextPart->cls[tv->tcellexpath];
        for (i=tv->tcellexpath + 1; i<iend ; i++) {
            if (NextCand->lab[i] < min) {
                min = NextCand->lab[i];
                tmp = i;
            }
        }
    }
    else {
        tmp = tv->tcellexpath+KRAN(NextPart->cls[tv->tcellexpath]);
    }
//    printf("Indiv %d; ", NextCand->lab[tmp]+labelorg);
    if (NextPart->cls[tv->tcellexpath] == 2) {
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tv->tcellexpath]);
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tv->tcellexpath+1]);
    }
    else {
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tmp]);
    }
    
    tv->indiv_vtx = NextCand->lab[tmp];
    Individualize(NextPart, NextCand, NextCand->lab[tmp], tv->tcellexpath, NextPart->cells, tv->tcellexpath + NextPart->cls[tv->tcellexpath]-1);
    
	tv->stats->numnodes++;
	if (tv->compstage == 0) {
//		traces_refine_notrace(g_arg, 
        traces_refine_notrace(NextCand, 
							  //m, 
							  n, 
							  NextPart, tv, ti);
        //PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 8410);
	}
	else {
		if (tv->tolevel_tl == tv->maxtreelevel+1) {
			trieref = trieroot;
//            tv->answ = traces_refine_comptrie(g_arg, 
			tv->answ = traces_refine_comptrie(NextCand, 
											  //m, 
											  n, 
											  NextPart, tv, ti);
			if (tv->answ == 0 ) {
				tv->stats->interrupted++;
			}
		}
		else {
			//traces_refine_notrace(g_arg, 
            traces_refine_notrace(NextCand, 
								  //m, 
								  n, 
								  NextPart, tv, ti);
		}
	}
}

trie* trie_new(int n, struct TracesVars* tv) {
	TrieArray[0] = malloc(n*sizeof(trie));
	if (TrieArray[0] == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	TrieArray[0][0].first_child = TrieArray[0][0].next_sibling = NULL;
	tv->triepos = 0;
	tv->trienext = 1;
	return TrieArray[0];
}

struct trie *trie_make(trie *t, int value, int n, struct TracesVars* tv) {
	trie *t1;
	t1 = t;
	if (tv->trienext == n) {
		tv->trienext = 0;
		tv->triepos++;
		TrieArray[tv->triepos] = malloc(n*sizeof(trie));
		if (TrieArray[tv->triepos] == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
	}
	if (t->first_child) {
		t = t->first_child;
		if (value < t->value) {
			t1->first_child = &TrieArray[tv->triepos][tv->trienext++];
			t1->first_child->next_sibling = t;
			t1->first_child->first_child = NULL;
			t = t1->first_child;
			t->value = value;
			return t;
		}
		while (value > t->value) {
			t1 = t;
			if (t->next_sibling) {
				t = t->next_sibling;
			}
			else break;
		}
		if (value == t->value) {
			return t;
		}
		t1->next_sibling = &TrieArray[tv->triepos][tv->trienext++];
		t1->next_sibling->first_child = t1->next_sibling->next_sibling = NULL;
		if (t != t1) {
			t1->next_sibling->next_sibling = t;
		}
		t = t1->next_sibling;
	}
	else {
		t->first_child = &TrieArray[tv->triepos][tv->trienext++];
		t = t->first_child;
		t->first_child = t->next_sibling = NULL;
	}
	t->value = value;
	return t;
}

struct trie *trie_comp(trie *t, int value) {
	if (t->first_child) {
		t = t->first_child;
		while (t) {
			if  (value != t->value) {
				t = t->next_sibling;
			}
			else {
				break;
			}
		}
		return t;
	}
	else {
		return NULL;
	}
}

void CopyCand(Candidate *W, Candidate *V,int n, int *lab, int *invlab) {
	if (lab) {
        memcpy(W->lab, lab, n*sizeof(int));
        memcpy(W->invlab, invlab, n*sizeof(int));
    }
    else {
        memcpy(W->lab, V->lab, n*sizeof(int));
        memcpy(W->invlab, V->invlab, n*sizeof(int));
    }
	W->name = V->name;
    W->vertex = V->vertex;
	W->code = V->code;
	W->singcode = V->singcode;
	W->firstsingcode = V->firstsingcode;
	W->do_it = V->do_it;
    W->sortedlab = FALSE;
}

void RemoveFromLevel(int from, int to, int strategy, boolean reinit) {
	int i;
    
	for (i=from; i<=to; i++) {
		if (Spine[i].listend) {
			(Spine[i].listend)->next = GarbList;
			GarbList = Spine[i].liststart;
			Spine[i].liststart = Spine[i].listend = NULL;
		}
        if (strategy == 0 || reinit) {
            Spine[i].listcounter = 0;
            if (i>from) {
                Spine[i].thetracexists = FALSE;
                Spine[i].part->code = -1;
            }
        }
	}
}

int CompStage0(Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand, 
               int m, int n, struct TracesVars* tv, struct TracesInfo *ti) {
	int i, j, i1, j2, k, cu, cu1, num_indv;
	int temp, tmp, auxcode, search_vtx, gom_level;
	boolean closeloop, firstsing, has_nexttcell;
	Candidate *SpTLliststart, *AuxCand;
    searchtrie *TreeNode, *TreeNode1, *TreeNode2;
	
#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt) return NAUTY_KILLED;
#else
    if (nauty_kill_request) return NAUTY_KILLED;
#endif
    
    PRINT_FROM_VERB(3)
    
	if (TargetCell(CurrCand, CurrPart, n, tv, tv->tolevel)) {
        ++tv->tolevel;
		SpineTL = Spine+tv->tolevel;
		SpineTL->tgtcell = tv->tcell;
        SpineTL->tgtsize = CurrPart->cls[tv->tcell];
		SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
		SpineTL->tgtpos = SpineTL->tgtend - 1;
	}
	else {
        tv->finalnumcells = CurrPart->cells;
		ti->thereisnextlevel = SelectNextLevel(n, tv);
		return 0;
	}
    
    tv->newgotonode = NULL;
	
	/*  CANDIDATE */
	temp = CurrCand->lab[Spine[1].tgtpos];
	k = SpineTL->tgtend;
	
    TreeNode = CurrCand->stnode;
    while (TreeNode) {
        if (TreeNode->goes_to) {
            CurrCand->do_it = FALSE;
            break;
        }
        TreeNode = TreeNode->father;
    }
    
	if (CurrCand->do_it) {
        if ((tv->orbits[temp] == temp) || tv->tolevel == 1) {
            ti->minimalinorbits = TRUE;
            if ((!ti->identitygroup) && (((Spine[tv->fromlevel].liststart != Spine[tv->fromlevel].listend) && (CurrPart->cls[tv->tcell] > 6))
                || tv->strategy || (tv->expathlength <=10))) {

                TempOrbits = NULL;
                tv->samepref = FixBase(fix, tv, CurrCand, 0, tv->fromlevel);
                if (tv->samepref != tv->nfix || ti->thegrouphaschanged) {
//                    printf("GETORBITSMIN\n");

                    if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
                    gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit, CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
                    if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
                    
                    ti->thegrouphaschanged = FALSE;
                    
                    if (gom_level < tv->nfix) {
                        PRINT_NOTMIN_VERB(3)
                        
                        TreeNode = CurrCand->stnode;
                        j2 = CurrCand->lab[Spine[gom_level+1].tgtpos];
                        i1 = tv->currorbit[j2];
                        for (j=0; j < tv->nfix - gom_level; j++) {
                            TreeNode = TreeNode->father;
                        }
                        TreeNode1 = TreeNode->first_child;
                        while (TreeNode1) {
                            if (TreeNode1->vtx == i1) {
                                break;
                            }
                            TreeNode1 = TreeNode1->next_sibling;
                        }
                        if (TreeNode1) {
                            while (TreeNode1->goes_to) {
                                TreeNode1 = TreeNode1->goes_to;
                            }
                            TreeNode2 = TreeNode1->next_sibling;
                            while (TreeNode2->vtx != j2) {
                                TreeNode2 = TreeNode2->next_sibling;
                            }
                            TreeNode1->index += TreeNode2->index;
                            TreeNode2->goes_to = TreeNode1;
                            
                            ti->minimalinorbits = FALSE;
                        }
                        else {
                            tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                        }
                    }
                }
                else {
                    tv->currorbit = findcurrorbits(gpB, tv->nfix);
                }
            }
            else {
                TempOrbits = WorkArray1;
                memcpy(TempOrbits, IDENTITY_PERM, n*sizeof(int));  // XXXYYY
                memcpy(TempOrbList, IDENTITY_PERM, n*sizeof(int));  // XXXYYY
//                printf("A memcpy(TempOrbits, IDENTITY_PERM), memcpy(TempOrbList, IDENTITY_PERM)\n");
//                for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
//                    TempOrbits[CurrCand->lab[i]] = TempOrbList[CurrCand->lab[i]] = CurrCand->lab[i];
//                }
                //                printf("TEMP: ");
                //                PrintVect(TempOrbits, 0, n, labelorg);
                //                printf("LIST: ");
                //                PrintVect(TempOrbList, 0, n, labelorg);
                
                tv->conta1++;
                tv->currorbit = TempOrbits;
            }
            
            if (ti->minimalinorbits) {
                //tv->auxtime2 -= CPUTIME;
                memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));  // XXXYYY
                memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));  // XXXYYY
                //tv->auxtime2 += CPUTIME;
//                printf("B memcpy(NextCand->lab, CurrCand->lab), memcpy(NextCand->invlab, CurrCand->invlab)\n");
                tv->conta2++;
                auxcode = CurrCand->code;
                SpineTL->trcstart = CurrPart->cells;
                TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
                if (!CurrCand->sortedlab) {
                    quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                    for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                        CurrCand->invlab[CurrCand->lab[i]] = i;
                    }
                    CurrCand->sortedlab = TRUE;
                }
                
                tv->indivstart = tv->tcell+CurrCand->indnum;
                tv->indivend = tv->indivstart+tv->steps;
                if (tv->indivend > SpineTL->tgtend) {
                    tv->indivend = SpineTL->tgtend;
                }
                
                temp = CurrCand->lab[tv->indivstart];
//                printf("%d: ",tv->indivstart);
//                PrintVect(CurrCand->lab,tv->indivstart, tv->indivend,labelorg);
                for (k = tv->indivstart; k < tv->indivend; k++) {
//                    if (TempOrbits) {
//                        PrintVect(TempOrbits,0,n,labelorg);
//                        PrintVect(TempOrbList,0,n,labelorg);
//                        printf("---\n");
//                    }
                    
//                    putorbits(outfile, tv->currorbit, 0, n);
                    
                    CurrCand->indnum++;
                    NextCand->singcode = CurrCand->singcode;
                    NextCand->vertex = CurrCand->lab[k];
                    NextCand->name = ++tv->name;
                    if (NextCand->name == (NAUTY_INFINITY-2)) {
                        NextCand->name = tv->name = 1;
                    }
                    
                    PRINT_INDIV_VERB(3)
                    if (tv->currorbit[NextCand->vertex] != NextCand->vertex) {
                        PRINT_SKIPPED_VERB(3)
                        
                        search_vtx = tv->currorbit[NextCand->vertex];
                        TreeNode = CurrCand->stnode;
                        if (TreeNode->first_child) {
                            TreeNode = TreeNode->first_child;
                            while (TreeNode) {
                                if (TreeNode->vtx == search_vtx) {
                                    break;
                                }
                                TreeNode = TreeNode->next_sibling;
                            }
                            if (TreeNode) {
                                while (TreeNode->goes_to) {
                                    TreeNode = TreeNode->goes_to;
                                }
                                TreeNode->index++;
                                continue;
                            }
                        }
                    }
                    PRINT_REFINE_VERB(3)
                    
                    //tv->auxtime3 -= CPUTIME;
//                    PrintVect(NextPart->cls,0,n,0);
//                    PrintVect(CurrPart->cls,0,n,0);
//                    PrintVect(NextPart->inv,0,n,0);
//                    PrintVect(CurrPart->inv,0,n,0);
//                    if (tv->conta3 > 0) {
//                        for (i=0; i<n; i+=CurrPart->cls[i]) {
//                            if (CurrPart->cls[i] != NextPart->cls[i]) {
//                                printf("Cell %d split into ",i);
//                                for (j=i; j<i+CurrPart->cls[i]; j+=NextPart->cls[j]) {
//                                    printf("%d ",j);
//                                }
//                                printf("\n");
//                            }
//                        }
//                    }
                    memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));  // XXXYYY
                    memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));  // XXXYYY
//                    printf("C memcpy(NextPart->cls, CurrPart->cls), memcpy(NextPart->inv, CurrPart->inv)\n");
                    tv->conta3++;
                    //                    printf("C ");
                    //tv->auxtime3 += CPUTIME;
                    
                    //                    printf("Target cell: %d ", tv->tcell);
                    //                    PrintVect(NextCand->lab+tv->tcell, 0, NextPart->cls[tv->tcell], labelorg);
                    if (NextPart->cls[tv->tcell] == 2) {
                        num_indv = 2;
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                        //                        printf("((%d)) ", CurrCand->lab[tv->tcell]+labelorg);
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                        //                        printf("((%d)) ", CurrCand->lab[tv->tcell+1]+labelorg);
                        if (SpineTL->singstart == SpineTL->singend) {
                            Singletons[SpineTL->singend++] = tv->tcell;
                            Singletons[SpineTL->singend++] = tv->tcell+1;
                        }
                    }
                    else {
                        num_indv = 1;
                        NextCand->singcode = MASHCOMM(NextCand->singcode, NextCand->vertex);
                        //                        printf("((%d)) ", NextCand->vertex+labelorg);
                        if (SpineTL->singstart == SpineTL->singend) {
                            Singletons[SpineTL->singend++] = tv->tcell + NextPart->cls[tv->tcell] - 1;
                            //Singletons[SpineTL->singend++] = tv->tcell + NextPart->cls[tv->tcell] - 1;
                        }
                    }
                    //                    printf("SpineTL->singend: %d\n", SpineTL->singend);
                    
                    Individualize(NextPart, NextCand, NextCand->vertex, tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                    tv->stats->numnodes++;
//                    printf("Ref(%d) ", NextCand->vertex+labelorg); PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 8886);

                    //tv->answ = traces_refine(g_arg, 
                    tv->answ = traces_refine(NextCand, 
                                             //m, 
                                             n, 
                                             NextPart, tv, ti, num_indv, TRUE);
//                    printf("--> "); PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 8893);

//
//                    if (!tv->strategy && !tv->options->getcanon && (NextPart->cells == tv->finalnumcells) && (tv->tolevel > 1)) {
//                        temp = 0;
//                        for (i=0; i<tv->tolevel; i++) {
//                            temp += Spine[i].listcounter;
//                        }
//                        if (temp > 5) {
//                            EXITFROMSTAGE0REFINE
//                        }
//                    }
//                    printf("ANSW: %d; tv->tolevel: %d, tv->tolevel_tl: %d, tv->tcellevel: %d\n", tv->answ, tv->tolevel, tv->tolevel_tl, tv->tcellevel);
                    switch (tv->answ) {
                        case 0:				/* Interrupted refinement: do not add to the list */
                            tv->stats->interrupted++;
                            SpineTL->levelcounter++;
                            break;
                        case 1 :			/* The same trace has been found once more : add to the list */
                            SpineTL->levelcounter++;
                            
                            NextCand->do_it = TRUE;
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel);
                            
                            tv->tolevel_tl = tv->tolevel;
                            NextCand->pathsingcode = NextCand->singcode;
                            NextCand->firstsingcode = 0;
                            
                            if (tv->steps > 1) {
                                    closeloop = CheckForMatching(CurrCand, NextCand, NextPart, tv, ti, m, n);
                                //printf("return CFM\n");
                                    if (NextCand->do_it) {
                                        firstsing = TRUE;
                                        if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                        
                                        /* EXPERIMENTAL PATH */
                                        //while (NextPart->cells < tv->finalnumcells) {
                                        while (NextPart->cells < n) {
//                                            printf("firstsing: %d, BreakSteps[tv->tolevel]: %d; ", firstsing, BreakSteps[tv->tolevel]);
                                            if (firstsing && BreakSteps[tv->tolevel]) {
                                                firstsing = FALSE;
                                                NextCand->firstsingcode = NextCand->pathsingcode;
                                                if (CheckForSingAutomorphisms(CurrCand, NextPart, NextCand, tv, ti, m, n))
                                                    if (!NextCand->do_it) {
//                                                        printf("brk1; ");
                                                        break;
                                                    }
                                            }
                                            if (!TargetCellExpPath(NextCand, NextPart, tv)) {
                                                NextCand->firstsingcode = NextCand->pathsingcode;
//                                                printf("brk2; ");
                                                break;
                                            }
//                                            printf("ES1 (%d cells); ", NextPart->cells);
                                            ExperimentalStep(NextPart, NextCand, tv, ti, m, n);
                                            PRINT_EXPPATHSTEP(NextCand, TRUE)
//                                            printf("(%d cells) ", NextPart->cells);
                                        }
                                        
                                        //if (NextPart->cells == tv->finalnumcells) {
//                                            UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                        //}
                                        if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                    }
                                    else {
                                        if (closeloop < tv->tolevel) k = SpineTL->tgtend;
                                        PRINT_RETURN
                                        break;
                                    }
                                //}
                                
                                //if (!tv->strategy && !tv->options->getcanon && (NextPart->cells == tv->finalnumcells) && (tv->tolevel_tl == tv->tolevel + 1)) {
                                if (!tv->strategy && !tv->options->getcanon && (tv->tolevel_tl == tv->tolevel + 1)) {
//                                    printf("esco aaa (tv->tolevel: %d)\n", tv->tolevel);
                                    tv->levelfromCS0 = tv->tolevel;
                                    tv->maxtreelevel = tv->tolevel_tl;
                                    if (tv->tolevel == 1) {
                                        tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                        EXITFROMSTAGE0EXPATH1
                                    }
                                    else {
                                        temp = 0;
                                        for (i=0; i<tv->tolevel; i++) {
                                            temp += Spine[i].listcounter;
                                        }
                                        if (temp > 5) {
                                            tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                            EXITFROMSTAGE0EXPATH1
                                        }
                                    }
                                }
                                
                                /* ANY AUTOMORPHISM? */
                                if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                tv->newindex = 0;
                                if (NextCand->do_it) {
//                                    printf("CFA-8972\n");
                                    closeloop = CheckForAutomorphisms(CurrCand, NextCand, tv, ti, m, n, NextPart);
                                    if (!NextCand->do_it && closeloop < tv->tolevel) k = SpineTL->tgtend;
                                }
                                if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                
                                if (NextCand->do_it) {
                                    ADDTONEXTLEVEL;
                                    SpineTL->keptcounter++;
                                    searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                }
                            }
                            else {
                                if (BreakSteps[tv->tolevel]) {
                                    NextCand->firstsingcode = NextCand->pathsingcode;
                                    if (CheckForSingAutomorphisms(CurrCand, NextPart, NextCand, tv, ti, m, n))
                                        if (!NextCand->do_it) {
                                            PRINT_RETURN
                                            break;
                                        }
                                }
                                
                                /* ANY AUTOMORPHISM? */
                                if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                tv->newindex = 0;
                                if (NextCand->do_it) {
//                                    printf("CFA-8998\n");
                                    closeloop = CheckForAutomorphisms(CurrCand, NextCand, tv, ti, m, n, NextPart);
                                    if (!NextCand->do_it && closeloop < tv->tolevel) k = SpineTL->tgtend;
                                }
                                if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                
                                if (NextCand->do_it) {
                                    ADDTONEXTLEVEL;
                                    SpineTL->keptcounter++;
                                    searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                }
                            }
                            PRINT_RETURN
                            break;
                        case 2 :	/* Delete the old list and start a new one: a better trace has been found */
//                            printf("AAA tv->tolevel: %d, tv->treedepth: %d\n",tv->tolevel, tv->treedepth);
//                            PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 8694);

//                            for (i=0; i<tv->tcellevel; i++) {
//                                printf("Targets[%d]. tgtfrom: %d, tgtcell: %d, tgtsize: %d, tgtend: %d, tgtpos: %d\n", i,
//                                       Spine[i].tgtfrom, Spine[i].tgtcell, Spine[i].tgtsize, Spine[i].tgtend, Spine[i].tgtpos);
//                            }

                            tv->tolevel_tl = tv->tolevel;
                            has_nexttcell = FALSE;
                            if (NextPart->cells == n) {   // ???? CONTROLLARE
                                tv->stats->canupdates++;
                                if (tv->options->usercanonproc != NULL)
                                {
                                    (*tv->options->usercanonproc)((graph*)tv->input_graph, NextCand->lab, (graph*)tv->cangraph, tv->stats->canupdates, NextCand->code, m, n);
                                }
                            }

                            if (tv->tolevel > tv->treedepth) {
                                tv->treedepth = tv->tolevel;
                                if (tv->strategy) {
                                    NEWPART(SpineTL->part);
                                }
                                else {
                                    NEWPARTSPINE(tv->tolevel);
                                }
                            }

                            if (!tv->strategy && (tv->tolevel > 1) && !SpineTL->liststart) {
                                /* First Candidate at current level */
//                                memcpy(NextCand->lab, TEMPLAB, n*sizeof(int));  // XXXYYY
//                                memcpy(NextCand->invlab, TEMPINVLAB, n*sizeof(int));  // XXXYYY
//                                tv->conta4++;
                                tv->maxtreelevel = tv->tolevel;
                                
                                SpineTL->liststart = NewCandidate(n, &GarbList, TRUE);
                                SpineTL->listend = SpineTL->liststart;

//                                if (NextPart->cells != tv->finalnumcells) NextCand->code = auxcode;
                                tv->conta0++;
//                                CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL);  // XXXYYY
                                CopyCand(SpineTL->liststart, NextCand, n, TEMPLAB, TEMPINVLAB);  // XXXYYY
//                                printf("a CopyCand(SpineTL->liststart, NextCand)\n");
                                if (NextPart->cells != tv->finalnumcells) SpineTL->liststart->code = auxcode;
                                COPYPART(SpineTL->part, NextPart);  // XXXYYY
//                                printf("b COPYPART(SpineTL->part, NextPart)\n");
                                tv->newindex = 0;
                                tv->newst_stage1 = searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                
                                SpineTL->listcounter = 1;
                                SpTLliststart = SpineTL->liststart;
                                
                                i = tv->tolevel;
                                if (tv->brkstpcount) {
                                    while ((i<n) && !BreakSteps[i]) {
                                        i++;
                                    }
                                    if (i<n) SpineTL->liststart->firstsingcode = Spine[i].singcode;
                                }
                                
                                SpineTL->updates = 1;
                                SpineTL->levelcounter = 1;
                                SpineTL->keptcounter = 1;
                                PRINT_LINE_PLUS(tv->fromlevel)
                                
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(SpineTL->liststart, tv->tolevel);
                                PRINT_RETURN;

                                if (!tv->strategy && !tv->options->getcanon && (tv->tolevel+1 == tv->firstpathlength)) {
//                                    printf("esco bbb (tv->tolevel: %d, tv->firstpathlength: %d)\n", tv->tolevel, tv->firstpathlength);
                                    if ((tv->tolevel == 1) && (CurrPart->cls[tv->tcell] > 5)) {
                                        //printf("TC: %d (%d), size: %d\n", tv->tcell, NextCand->lab[tv->tcell]+labelorg, CurrPart->cls[tv->tcell]);
                                        EXITFROMSTAGE0EXPATH2;
                                    }
                                    else {
                                        temp = 0;
                                        for (i=0; i<tv->tolevel; i++) {
                                            temp += Spine[i].listcounter;
                                        }
                                        if (temp > 5) {
                                            EXITFROMSTAGE0EXPATH2;
                                        }
                                    }
                                }
                            }
                            else {
                                memset(WorkArray, 0, n*sizeof(int));
                                
                                tv->lastcell = tv->lastlev = -1;
//                                Spine[tv->tolevel].tgtfrom = 0;
                                has_nexttcell = TargetCellFirstPath(NextCand, NextPart, tv);
//                                printf("has_nexttcell: %d; \n", has_nexttcell);
//                                PrintVect(NextCand->lab, Spine[tv->tolevel_tl].tgtcell, Spine[tv->tolevel_tl].tgtend, labelorg);
                                if (!has_nexttcell) {
                                    tv->stats->canupdates++;
                                    if (tv->options->usercanonproc != NULL)
                                    {
//                                        (*tv->options->usercanonproc)((graph*)g_arg, NextCand->lab, (graph*)tv->cangraph, tv->stats->canupdates, NextCand->code, m, n);
                                        (*tv->options->usercanonproc)((graph*)tv->input_graph, NextCand->lab, (graph*)tv->cangraph, tv->stats->canupdates, NextCand->code, m, n);
                                    }
                                }

//                                tv->maxtreelevel = tv->tolevel;
                                tv->tcellevel = tv->maxtreelevel = tv->tolevel;
//                                printf("SET(b) tv->tcellevel to %d\n", tv->tcellevel);
                                SpineTL->levelcounter++;
                                SpineTL->updates++;
                                SpineTL->keptcounter = 1;
                                
                                RemoveFromLevel(tv->tolevel, tv->treedepth, tv->strategy, TRUE);
                                SpineTL->liststart = NewCandidate(n, &GarbList, TRUE);
                                SpineTL->listend = SpineTL->liststart;
                                
                                tv->conta0++;
                                CopyCand(SpineTL->liststart, NextCand, n, NULL, NULL);  // XXXYYY
//                                printf("c CopyCand(SpineTL->liststart, NextCand)\n");
                                COPYPART(SpineTL->part, NextPart); // XXXYYY
//                                printf("d COPYPART(SpineTL->part, NextPart)\n");
                                tv->newindex = 0;
                                
                                tv->newst_stage1 = searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                
                                SpineTL->listcounter = 1;
                                SpTLliststart = SpineTL->liststart;
                                
                                SpTLliststart->pathsingcode = SpineTL->singcode = SpTLliststart->singcode;
                                SpTLliststart->firstsingcode = 0;
                                
                                PRINT_LINE
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(SpTLliststart, tv->tolevel);
                                
                                //tv->tolevel_tl = tv->tolevel;  // spostato sopra
                                
                                memset(BreakSteps, 0, n*sizeof(int));
                                tv->brkstpcount = 0;
                                
                                if (tv->steps > 1) {
//                                    printf("CS0, case 2e\n");

                                    if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                    
                                    /* EXPERIMENTAL PATH */
//                                    tv->lastcell = tv->lastlev = -1;
//                                    Spine[tv->tolevel].tgtfrom = 0;
//                                    printf("CS0, case tcfp ??? (NextPart->cells: %d, tv->finalnumcells: %d)\n", NextPart->cells, tv->finalnumcells);
                                    PRINTF2("CStage0 2: %d\n", tv->finalnumcells);
                                    tv->finalnumcells = n;
//                                    tv->smalldeglevel = n;
                                    PRINTF2("CStage0 2<: %d\n", tv->finalnumcells);

                                    //while (NextPart->cells < tv->finalnumcells) {
                                    while (has_nexttcell) {
//                                        printf("ExLoop; ");
//                                        //printf("CS0, case tcfp\n");
//                                        if (!TargetCellFirstPath(tv->graph, SpTLliststart, NextPart, n, tv)) {
//                                            tv->tolevel_tl--;
//                                            if (tv->tolevel_tl > 6) {
//                                                PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
//                                            }
//                                            break;
//                                        }
//                                        printf("ES2; \n");
                                        ExperimentalStep(NextPart, SpTLliststart, tv, ti, m, n);
                                        PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
                                        
                                        Spine[tv->tolevel_tl].singcode = SpTLliststart->pathsingcode;
                                        has_nexttcell = TargetCellFirstPath(SpTLliststart, NextPart, tv);
//                                        printf("Has NC: %d; ", has_nexttcell);
//                                        if (!has_nexttcell) {
//                                            //tv->tolevel_tl--;
//                                            if (tv->tolevel_tl > 6) {
//                                                PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
//                                            }
////                                            break;
//                                        }

                                    }
                                    if (NextPart->cells < n) {
                                        PRINTF2("CStage0 3: %d\n", tv->finalnumcells);
//                                        tv->smalldeglevel = tv->tolevel;
                                        tv->finalnumcells = NextPart->cells;
                                        PRINTF2("CStage0 3<: %d\n", tv->finalnumcells);
                                    }

//                                    while (TargetCellFirstPath(tv->graph, SpTLliststart, NextPart, n, tv)) {
////                                        printf("CS0, case tcfp\n");
//                                        ExperimentalStep(g_arg, NextPart, SpTLliststart, tv, ti, m, n);
//                                        PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
//                                        
//                                        Spine[tv->tolevel_tl].singcode = SpTLliststart->pathsingcode;
//                                    }
//
//                                    tv->tolevel_tl--;
//                                    if (tv->tolevel_tl > 6) {
//                                        PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
//                                    }
                                    PRINTF2("CS0 2?: finalnumcells: %d\n", tv->finalnumcells);
                                    if (NextPart->cells == tv->finalnumcells) {
                                        UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                    }
                                    
                                    if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                    //PRINTF2("CStage0 4: %d\n", tv->smalldeglevel);
//                                    if (tv->tolevel_tl > tv->smalldeglevel) {
//                                        tv->smalldeglevel = tv->tolevel_tl;
//                                    }
//                                    if ((tv->finalnumcells < n) && (tv->tolevel_tl < tv->smalldeglevel)) {
//                                        tv->smalldeglevel = tv->tolevel_tl;
//                                    }
                                    //PRINTF2("CStage0 4<: %d\n", tv->smalldeglevel);

                                    tv->firstpathlength = tv->tolevel_tl;
//                                    printf("set tv->firstpathlength = %d\n", tv->firstpathlength);
                                    PRINT_RETURN
                                    if (!tv->strategy && !tv->options->getcanon && (NextPart->cells == tv->finalnumcells) && (tv->tolevel_tl == tv->tolevel + 1)) {
//                                        printf("esco ccc (tv->tolevel: %d)\n", tv->tolevel);
                                        tv->maxtreelevel = tv->tolevel_tl;
                                        if ((tv->tolevel == 1) && (CurrPart->cls[tv->tcell] > 5)) {
                                            //printf("TC: %d (%d), size: %d\n", tv->tcell, NextCand->lab[tv->tcell]+labelorg, CurrPart->cls[tv->tcell]);
                                            //tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                            EXITFROMSTAGE0EXPATH2
                                        }
                                        else {
                                            temp = 0;
                                            for (i=0; i<tv->tolevel; i++) {
                                                temp += Spine[i].listcounter;
                                            }
                                            if (temp > 5) {
                                                //tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                                EXITFROMSTAGE0EXPATH2
                                            }
                                        }
                                    }
                                    //tv->auxtime5 -= CPUTIME;
//                                    TEMPLAB = SpTLliststart->lab;
//                                    TEMPINVLAB = SpTLliststart->invlab;
                                    memcpy(TEMPLAB, SpTLliststart->lab, n*sizeof(int));  // XXXYYY
                                    memcpy(TEMPINVLAB, SpTLliststart->invlab, n*sizeof(int));  // XXXYYY
//                                    printf("D memcpy(TEMPLAB, SpTLliststart->lab), memcpy(TEMPINVLAB, SpTLliststart->invlab)\n");
                                    //tv->auxtime5 += CPUTIME;
                                    //                                    printf("E ");
                                    tv->conta5++;
                                }
                                else {
//                                    printf("CS0, case 2f\n");

                                    PRINT_RETURN
                                }
                            }
                            
                            break;
                        default:
                            break;
                    }
                } /* end for */
            }
        }
    }
    
	/* REMOVE CURRENT CANDIDATE */
	if (SpineFL->liststart && (k >= SpineTL->tgtend)) {
		SpineFL->liststart = CurrCand->next;
		if (CurrCand->next == NULL) {
			SpineFL->listend = NULL;
		}
		SpineFL->listcounter--;
		CurrCand->next = GarbList;
		GarbList = CurrCand;
	}
    //	ti->thereisnextlevel = SelectNextLevel(g_arg, n, tv);
	ti->thereisnextlevel = SelectNextLevel(n, tv);
	return 0;
}

int CompStage1(Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand, 
               int m, int n, 
               struct TracesVars* tv, struct TracesInfo *ti) {
	int i, k, cu, cu1, tmp, gom_level, search_vtx;
    searchtrie *TreeNode;
    
    //printf("CS1 (nkr: %d); ", nauty_kill_request);
#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt) return NAUTY_KILLED;
#else
    if (nauty_kill_request) return NAUTY_KILLED;
#endif
    
    CurrCand->stnode = tv->newst_stage1;
//    printf("CS1 start; finalnumcells: %d\n", tv->finalnumcells);

    tv->tolevel++;
	SpineTL = Spine+tv->tolevel;
	tv->tcell = SpineTL->tgtcell;
    SpineTL->levelcounter = 0;
    SpineTL->keptcounter = 0;
    SpineTL->updates = 1;
    
	
	if (tv->options->verbosity >= 2) {
		LINE(tv->linelgth-3, "=");
		NEXTLINE
	}
	
	memset(RefCells, 0, n*sizeof(int));
	memset(MultRefCells, 0, n*sizeof(int));
	ti->thegrouphaschanged = TRUE;
	
	/*  CANDIDATE */
	memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
	memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
	NextCand->do_it = TRUE;
	SpineTL->trcstart = CurrPart->cells;
	
	tv->indivstart = tv->tcell;
	tv->indivend = SpineTL->tgtend;
//	if (g_arg->d[CurrCand->lab[tv->indivstart]] == 1) {
    if (TheGraph[CurrCand->lab[tv->indivstart]].d == 1) {
		tv->indivstart = SpineTL->tgtend-1;
	}
    
    FixBase(fix, tv, NextCand, 0, tv->fromlevel);
    
	if (!ti->identitygroup) {
        if (tv->options->verbosity >= 2) tv->schreier2 -= CPUTIME;
        tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
        if (tv->options->verbosity >= 2) tv->schreier2 += CPUTIME;
    }
    else {
        if (n / CurrPart->cls[tv->tcell] < 256) {
            memcpy(tv->currorbit, IDENTITY_PERM, n*sizeof(int));
        }
        else {
            for (k = tv->indivstart; k < tv->indivend; k++) {
                tv->currorbit[CurrCand->lab[k]] = CurrCand->lab[k];
            }
        }
    }
    
    if (!CurrCand->sortedlab) {
        quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
        for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
            CurrCand->invlab[CurrCand->lab[i]] = i;
        }
        CurrCand->sortedlab = TRUE;
    }
//    printf("Indiv: ");
//    PrintVect(CurrCand->lab, tv->indivstart, tv->indivend, labelorg);
	for (k = tv->indivstart; k < tv->indivend; k++) {
        NextCand->vertex = CurrCand->lab[k];
        NextCand->name = ++tv->name;
        if (NextCand->name == (NAUTY_INFINITY-2)) {
            NextCand->name = tv->name = 1;
        }
//        printf("tv->currorbit[%d]: %d\n", CurrCand->lab[k]+labelorg, tv->currorbit[CurrCand->lab[k]]+labelorg);
        if (tv->currorbit[CurrCand->lab[k]] != CurrCand->lab[k]) {
            search_vtx = tv->currorbit[NextCand->vertex];
            TreeNode = CurrCand->stnode;
            if (TreeNode->first_child) {
                TreeNode = TreeNode->first_child;
                while (TreeNode) {
                    if (TreeNode->vtx == search_vtx) {
                        break;
                    }
                    TreeNode = TreeNode->next_sibling;
                }
                if (TreeNode) {
                    while (TreeNode->goes_to) {
                        TreeNode = TreeNode->goes_to;
                    }
                    TreeNode->index++;
                }
            }
			continue;
		}
        PRINT_REFINE_VERB(3)
		memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
		memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
		
		Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
		
		tv->stats->numnodes++;
        SpineTL->levelcounter++;
        tv->tolevel_tl = tv->tolevel;
		trieref = trieroot;
        SpineTL->levelcounter++;
        
//		traces_refine_maketrie(g_arg, 
        traces_refine_maketrie(NextCand, 
							   //m, 
							   n, 
							   NextPart, tv, ti);
		
//        PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 9434);
//        PrintVect(NextPart->inv, 0, n, 0);

		RefCells[CurrCand->lab[k]] = NextPart->cells;
//        printf("cells: %d; ", NextPart->cells);
        PRINTF2("CS1 1?: finalnumcells: %d\n", tv->finalnumcells);
		if (NextPart->cells == tv->finalnumcells) {
			tv->answ=1;
			if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                
//                if (NextPart->cells < tv->finalnumcells) {
                if (NextPart->cells != tv->finalnumcells) {
                    if (!Spine[tv->tolevel].part) {
                        NEWPART(Spine[tv->tolevel].part)
                    }
                }
            
			/* ANY AUTOMORPHISM? */
			if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
            
            PRINTF2("CS1 2?: finalnumcells: %d\n", tv->finalnumcells);
            if (NextPart->cells == tv->finalnumcells) {
//                printf("CFA-9424\n");
                CheckForAutomorphisms(CurrCand, NextCand, tv, ti, m, n, NextPart);
            }
            else {
                CheckForSingAutomorphisms(CurrCand, NextPart, NextCand, tv, ti, m, n);
            }
            
			if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
			
			PRINT_RETURN
			
			/* ADD TO NEXT LEVEL */
			if (tv->answ == 1) {
                SpineTL->keptcounter++;
                if (!Spine[tv->tolevel].listend) COPYPART(Spine[tv->tolevel].part, NextPart);
				ADDTONEXTLEVEL;
                searchtrie_make(CurrCand, SpineTL->listend, n, tv);
			}
		}
	} /* end for */
    PRINTF2("CS1 3: finalnumcells: %d\n", tv->finalnumcells);
	for (k = tv->indivstart; k < tv->indivend; k++) {
        MultRefCells[RefCells[tv->currorbit[CurrCand->lab[k]]] % n]++;
        //		MultRefCells[RefCells[tv->currorbit[CurrCand->lab[k]]] % tv->finalnumcells]++;
        //		MultRefCells[RefCells[tv->currorbit[CurrCand->lab[k]]]]++;
	}
//    MultRefCells[tv->finalnumcells] = 0;
    
	if (tv->options->verbosity >= 2) {
        if (MultRefCells[0]) {
            fprintf(outfile, tv->digstring, n);
            fprintf(outfile, "cells: %d\n", MultRefCells[0]);
        }
        for (k=1; k<n; k++) {
			if (MultRefCells[k]) {
				fprintf(outfile, tv->digstring, k);
				fprintf(outfile, "cells: %d\n", MultRefCells[k]);
			}
		}
	}
	
#if !MAXN
    DYNALLOC1(searchtrie*, RefPath, RefPath_sz, tv->tolevel, "Traces-CS1");
#endif
    
    TreeNode = CurrCand->stnode;
    while (TreeNode) {
        RefPath[TreeNode->level] = TreeNode;
        TreeNode = TreeNode->father;
    }
    
	/* REMOVE CURRENT CANDIDATE */
	SpineFL->liststart = CurrCand->next;
	if (CurrCand->next == NULL) {
		SpineFL->listend = NULL;
		SpineFL->listcounter = 1;
	}
	SpineFL->listcounter--;
	CurrCand->next = GarbList;
	GarbList = CurrCand;
	
	if (tv->options->verbosity >= 2) {
		LINE(tv->linelgth, "=");
		NEXTLINE
	}
	tv->compstage = 2;
	tv->steps = n;
    
	if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
//    printf("* ");
	gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit, 
							 CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
	if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
	ORBITSIZES
    //	ti->thereisnextlevel = SelectNextLevel(g_arg, n, tv);
	ti->thereisnextlevel = SelectNextLevel(n, tv);
    PRINTF2("CS1 4: finalnumcells: %d\n", tv->finalnumcells);
	SpineTL->part->cells = tv->finalnumcells;
    
    AutomCount[0] = 2;
    AutomCount[1] = CurrCand->vertex;
    
    return 0;
}

int CompStage2(Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand, 
               int m, int n, 
               struct TracesVars* tv, struct TracesInfo *ti) {
    //	int i, j, i1, j2, k, cu, cu1, vertex, search_vtx, gom_level;
	int i, j, i1, j2, k, cu, cu1, vertex, gom_level;
	int temp, tmp, autom;
	Candidate *AuxCand;
    searchtrie *TreeNode, *TreeNode1, *TreeNode2;
    int *CuOrb;
	
//    printf("CS2 (level: %d)\n", tv->tolevel);
#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt) return NAUTY_KILLED;
#else
    if (nauty_kill_request) return NAUTY_KILLED;
#endif
    
    autom = 0;
	
    TreeNode = CurrCand->stnode;
    tv->cand_level = 0;
    //tv->tolevel = tv->levelfromCS0;
    
    while (TreeNode) {
        if (TreeNode->goes_to) {
            CurrCand->do_it = FALSE;
        }
        if (!tv->cand_level && TreeNode == RefPath[TreeNode->level]) {
            tv->cand_level = TreeNode->level;
        }
        TreeNode = TreeNode->father;
    }
    if (tv->cand_level+1 == tv->maxtreelevel) {
        ti->useTempOrbits1 = TRUE;
    }
    else {
        ti->useTempOrbits1 = FALSE;
    }
    if (tv->cand_level == tv->fromlevel) {
        ti->useTempOrbits2 = TRUE;
    }
    else {
        ti->useTempOrbits2 = FALSE;
    }
    
    PRINT_FROM_VERB(3)
    
    if (CurrCand->do_it) {
        if (tv->tolevel == 0) {
            tv->fromlevel = tv->tolevel;
            SpineFL = Spine+tv->fromlevel;
            vertex = Spine[tv->maxtreelevel+1].liststart->lab[Spine[1].tgtpos];
            k = n;
            
            //if (TargetCell(g_arg, CurrCand, CurrPart, n, tv, tv->tolevel)) {
            if (TargetCell(CurrCand, CurrPart, n, tv, tv->tolevel)) {
                ++tv->tolevel;
                SpineTL = Spine+tv->tolevel;
                SpineTL->tgtcell = tv->tcell;
                SpineTL->tgtsize = CurrPart->cls[tv->tcell];
                SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
                SpineTL->tgtpos = SpineTL->tgtend - 1;
            }
            else {
                PRINTF2("CStage2 1: %d\n", tv->finalnumcells);
//                tv->smalldeglevel = tv->tolevel;
                tv->finalnumcells = CurrPart->cells;
                PRINTF2("CStage2 1<: %d\n", tv->finalnumcells);
                return 0;
            }
            
            memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
            memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
            SpineTL->trcstart = CurrPart->cells;
            TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
            
            tv->indivstart = tv->tcell+CurrCand->indnum;
            tv->indivend = tv->indivstart+tv->steps;
            if (tv->indivend > SpineTL->tgtend) {
                tv->indivend = SpineTL->tgtend;
            }
            memset(CurrRefCells, 0, n*sizeof(int));
            ti->thegrouphaschanged = TRUE;
            
            if (!CurrCand->sortedlab) {
                quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                    CurrCand->invlab[CurrCand->lab[i]] = i;
                }
                CurrCand->sortedlab = TRUE;
            }
//            printf("Indiv: ");
//            PrintVect(CurrCand->lab, tv->indivstart, tv->indivend, labelorg);

            for (k = tv->indivstart; k < tv->indivend; k++) {
                if ((tv->orbits[CurrCand->lab[k]] == CurrCand->lab[k]) && ((tv->finalnumcells < n) || (OrbSize[tv->orbits[CurrCand->lab[k]]] >= OrbSize[tv->orbits[vertex]]))) {
                // CONTROLLARE CURRORBITSIZES !!!
//                putorbits(outfile, tv->orbits, 0, n);
//                printf("Condition OrbSize[tv->orbits[%d]] >= OrbSize[tv->orbits[%d]]: %d\n", 
//                       CurrCand->lab[k]+labelorg, vertex+labelorg, OrbSize[tv->orbits[CurrCand->lab[k]]] >= OrbSize[tv->orbits[vertex]]);
//                PrintVect(OrbSize, 0, n, 0);
                //if (tv->orbits[CurrCand->lab[k]] == CurrCand->lab[k]) {
                    
                    CurrCand->indnum++;
                    NextCand->singcode = CurrCand->singcode;
                    NextCand->vertex = CurrCand->lab[k];
                    NextCand->name = ++tv->name;
                    if (NextCand->name == (NAUTY_INFINITY-2)) {
                        NextCand->name = tv->name = 1;
                    }
                    
                    if (ti->thegrouphaschanged) {
                        if (tv->fromlevel == tv->maxtreelevel) {
                            CURRORBITSIZES
                        }
                        ti->thegrouphaschanged = FALSE;
                    }
                    
                    if (tv->currorbit[CurrCand->lab[k]] != CurrCand->lab[k]) {
                        continue;
                    }
                    
                    memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                    memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                    if (NextPart->cls[tv->tcell] == 2) {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                    }
                    else {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[k]+labelorg);
                    }
                    
                    Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                    
                    tv->stats->numnodes++;
                    Spine[tv->tolevel+1].levelcounter++;
                    if (tv->fromlevel == tv->maxtreelevel) {
                        tv->tolevel_tl = tv->tolevel;
                        trieref = trieroot;
                        
//                        tv->answ = traces_refine_comptrie(g_arg, 
                        tv->answ = traces_refine_comptrie(NextCand, 
                                                          //m, 
                                                          n, 
                                                          NextPart, tv, ti);
                        
//                        PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 9675);
//                        printf("answ: %d; interr: %lu; ", tv->answ, tv->stats->interrupted);
//                        PrintVect(NextPart->inv, 0, n, 0);
                        
                        if (tv->answ) {
                            if (NextPart->cells != tv->finalnumcells) {
                            //if (NextPart->cells != n) {
                                //CurrRefCells[NextPart->cells] += CurrOrbSize[CurrCand->lab[k]];
                                CurrRefCells[NextPart->cells % n] += CurrOrbSize[CurrCand->lab[k]];
                                //if (CurrRefCells[NextPart->cells] > MultRefCells[NextPart->cells]) {
                                if (CurrRefCells[NextPart->cells % n] > MultRefCells[NextPart->cells % n]) {
                                    k = n;
                                    break;
                                }
                                continue;
                            }
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                }
                    }
                    else {
//                        tv->answ = traces_refine_sametrace(g_arg, 
                        tv->answ = traces_refine_sametrace(NextCand, 
                                                           //m, 
                                                           n, 
                                                           NextPart, tv, ti);

//                        PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 9697);
//                        printf("answ: %d; interr: %lu; ", tv->answ, tv->stats->interrupted);
//                        PrintVect(NextPart->inv, 0, n, 0);
                        
                        if (tv->answ) {
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                if (tv->tolevel == tv->maxtreelevel) {
                                    tv->tolevel_tl = tv->tolevel;
                                    if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                    //TargetCellExpPath(g_arg, NextCand, NextPart, n, tv);
                                    //TargetCellExpPath(tv->graph, NextCand, NextPart, n, tv);
                                    TargetCellExpPath(NextCand, NextPart, tv);
//                                    printf("ES3; ");
                                    ExperimentalStep(NextPart, NextCand, tv, ti, m, n);
                                    PRINT_EXPPATHSTEP(NextCand, tv->answ)
                                    PRINTF2("CS2 1?: finalnumcells: %d\n", tv->finalnumcells);
                                    if (NextPart->cells == tv->finalnumcells) {
                                        UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                    }
                                    
                                    if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                    if (!tv->answ) {
                                        PRINT_RETURN
                                    }
                                }
                        }
                    }
                    if (tv->answ) {
                        //printf("CF2; NextPart->cells: %d, tv->finalnumcells: %d\n", NextPart->cells, tv->finalnumcells);
                        PRINTF2("CS2 2?: finalnumcells: %d\n", tv->finalnumcells);
                        if (NextPart->cells == tv->finalnumcells) {
                            if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                            temp = (tv->tolevel_tl == tv->tolevel+1);
//                            printf("CFA-9690\n");
                            autom = CheckForAutomorphisms(CurrCand, NextCand, 
                                                          tv, ti, temp, n, NextPart);
                            if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                            
                            if (ti->thegrouphaschanged) {
                                ORBITSIZES
                            }
                        }
                        PRINT_RETURN
                        
                        /* ADD TO NEXT LEVEL */
                        PRINTF2_2("CS2 3?: cells: %d, finalnumcells: %d\n", NextPart->cells, tv->finalnumcells);
//                        if ((NextPart->cells < tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                        if ((NextPart->cells != tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                            ADDTONEXTLEVEL;
                            searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                        }
                    }
                    else {
                        tv->stats->interrupted++;
                    }
                    if (tv->fromlevel == tv->maxtreelevel) {
                        k = n;
                        break;
                    }
                }
            } /* end for */
        }
        else {
//            printf("AAA (tv->cand_level: %d)\n", tv->cand_level);
            temp = CurrCand->lab[Spine[1].tgtpos];
            vertex = Spine[tv->maxtreelevel+1].liststart->lab[Spine[1].tgtpos];
            k = n;
//            putorbits(outfile, tv->orbits, 0, n);
//            printf("condition OrbSize[tv->orbits[%d]] >= OrbSize[tv->orbits[%d]]: %d\n", 
//                   CurrCand->lab[k]+labelorg, vertex+labelorg, OrbSize[tv->orbits[CurrCand->lab[Spine[1].tgtpos]]] >= OrbSize[tv->orbits[vertex]]);
//            PrintVect(OrbSize, 0, n, 0);
//            printf("tv->orbits[%d] == %d\n", temp+labelorg, tv->orbits[temp]+labelorg);
            if (tv->cand_level || ((tv->orbits[temp] == temp) && ((tv->finalnumcells < n) || (OrbSize[tv->orbits[CurrCand->lab[Spine[1].tgtpos]]] >= OrbSize[tv->orbits[vertex]])))) {
//            if (tv->cand_level || (tv->orbits[temp] == temp)) {
//                printf("BBB\n");
                tv->fromlevel = tv->tolevel;
                SpineFL = Spine+tv->fromlevel;
                
                //if (TargetCell(g_arg, CurrCand, CurrPart, n, tv, tv->tolevel)) {
                if (TargetCell(CurrCand, CurrPart, n, tv, tv->tolevel)) {
                    tv->tcellevel = ++tv->tolevel;
//                    printf("SET(c) tv->tcellevel to %d\n", tv->tcellevel);  //???
                    SpineTL = Spine+tv->tolevel;
                    SpineTL->tgtcell = tv->tcell;
                    SpineTL->tgtsize = CurrPart->cls[tv->tcell];
                    SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
                    SpineTL->tgtpos = SpineTL->tgtend - 1;
                }
                else {
                    PRINTF2("CStage2 2: %d\n", tv->finalnumcells);
//                    tv->smalldeglevel = tv->tolevel;
                    tv->finalnumcells = CurrPart->cells;
                    PRINTF2("CStage2 2<: %d\n", tv->finalnumcells);
                    return 0;
                }
                ti->minimalinorbits = TRUE;
                
                if (!ti->identitygroup) {
                    
                    if (ti->useTempOrbits1 && ti->useTempOrbits2) {
                        CuOrb = TempOrbits;
                    }
                    else {
                        FixBase(fix, tv, CurrCand, 0, tv->fromlevel);
                        if (ti->useTempOrbits1 && tv->fromlevel == tv->maxtreelevel) {
                            tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                            CuOrb = tv->currorbit;
                        }
                        else {
                            if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
                            gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit, 
                                                     CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
                            if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
                            //                            printf("@ %.2f ", tv->schreier1);
                            //                            PrintVect(fix, 0, tv->nfix, labelorg);
                            CuOrb = tv->currorbit;
                            if (gom_level < tv->nfix) {
                                PRINT_NOTMIN_VERB(3)
                                if (ti->useTempOrbits1) {
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == CuOrb[tv->currorbit[CurrCand->vertex]]) break;
                                    if (i < AutomCount[0]) {
                                        AutomCount[AutomCount[0]++] = CurrCand->vertex;
                                    }
                                    ti->minimalinorbits = FALSE;
                                }
                                else {
                                    TreeNode = CurrCand->stnode;
                                    j2 = CurrCand->lab[Spine[gom_level+1].tgtpos];
                                    i1 = tv->currorbit[j2];
                                    for (j=0; j < tv->nfix - gom_level; j++) {
                                        TreeNode = TreeNode->father;
                                    }
                                    TreeNode1 = TreeNode->first_child;
                                    while (TreeNode1) {
                                        if (TreeNode1->vtx == i1) {
                                            break;
                                        }
                                        TreeNode1 = TreeNode1->next_sibling;
                                    }
                                    if (TreeNode1) {
                                        while (TreeNode1->goes_to) {
                                            TreeNode1 = TreeNode1->goes_to;
                                        }
                                        TreeNode2 = TreeNode->first_child;
                                        while (TreeNode2->vtx != j2) {
                                            TreeNode2 = TreeNode2->next_sibling;
                                        }
                                        
                                        TreeNode1->index += TreeNode2->index;
                                        TreeNode2->goes_to = TreeNode1;
                                        
                                        ti->minimalinorbits = FALSE;
                                    }
                                    else {
                                        tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                                    }
                                }
                            }
                        }
                    }
                    ti->thegrouphaschanged = FALSE;
                }
                else {
                    CuOrb = IDENTITY_PERM;
                }
                
                if (ti->minimalinorbits) {
                    memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
                    memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
                    SpineTL->trcstart = CurrPart->cells;
                    TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
                    
                    tv->indivstart = tv->tcell+CurrCand->indnum;
                    tv->indivend = tv->indivstart+tv->steps;
                    if (tv->indivend > SpineTL->tgtend) {
                        tv->indivend = SpineTL->tgtend;
                    }
                    memset(CurrRefCells, 0, n*sizeof(int));
                    if (!ti->identitygroup) ti->thegrouphaschanged = TRUE;
                    
                    if (!CurrCand->sortedlab) {
                        quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                        for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                            CurrCand->invlab[CurrCand->lab[i]] = i;
                        }
                        CurrCand->sortedlab = TRUE;
                    }
//                    printf("INDIV: ");
//                    PrintVect(CurrCand->lab, tv->indivstart, tv->indivend, labelorg);
    
                    for (k = tv->indivstart; k < tv->indivend; k++) {
                        CurrCand->indnum++;
                        NextCand->singcode = CurrCand->singcode;
                        NextCand->vertex = CurrCand->lab[k];
                        NextCand->name = ++tv->name;
                        if (NextCand->name == (NAUTY_INFINITY-2)) {
                            NextCand->name = tv->name = 1;
                        }
                        
                        if (ti->thegrouphaschanged) {
                            if (tv->fromlevel == tv->maxtreelevel) {
                                CURRORBITSIZES
                            }
                            ti->thegrouphaschanged = FALSE;
                        }
                        
                        
                        if (CuOrb[CurrCand->lab[k]] != CurrCand->lab[k]) {
                            
                            continue;
                        }
                        
                        memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                        memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                        if (NextPart->cls[tv->tcell] == 2) {
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                        }
                        else {
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[k]);
                        }
                        
                        Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                        
                        tv->stats->numnodes++;
                        Spine[tv->tolevel+1].levelcounter++;
                        if (tv->fromlevel == tv->maxtreelevel) {
                            tv->tolevel_tl = tv->tolevel;
                            trieref = trieroot;
                            
//                            tv->answ = traces_refine_comptrie(g_arg, 
                            tv->answ = traces_refine_comptrie(NextCand, 
                                                              //m, 
                                                              n, 
                                                              NextPart, tv, ti);

//                            PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 9928);
//                            printf("answ: %d; interr: %lu; ", tv->answ, tv->stats->interrupted);
//                            PrintVect(NextPart->inv, 0, n, 0);
                            
                            if (tv->answ) {
                                PRINTF2("CS2 4?: finalnumcells: %d\n", tv->finalnumcells);
                                if (NextPart->cells != tv->finalnumcells) {
//                                    CurrRefCells[NextPart->cells] += CurrOrbSize[CurrCand->lab[k]];
//                                    if (CurrRefCells[NextPart->cells] > MultRefCells[NextPart->cells]) {
                                    CurrRefCells[NextPart->cells % n] += CurrOrbSize[CurrCand->lab[k]];
                                    if (CurrRefCells[NextPart->cells % n] > MultRefCells[NextPart->cells % n]) {
                                        k = n;
                                        break;
                                    }
                                    continue;
                                }
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel);
                            }
                        }
                        else
                        {
//                            tv->answ = traces_refine_sametrace(g_arg, 
                            tv->answ = traces_refine_sametrace(NextCand, 
                                                               //m, 
                                                               n, 
                                                               NextPart, tv, ti);

//                            PrintPartition(NextCand->lab, NextPart->cls, n, labelorg, 9951);
//                            printf("answ: %d; interr: %lu; ", tv->answ, tv->stats->interrupted);
//                            PrintVect(NextPart->inv, 0, n, 0);

                            if (tv->answ) {
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                    if (tv->tolevel == tv->maxtreelevel) {
                                        tv->tolevel_tl = tv->tolevel;
                                        if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                        //if (TargetCellExpPath(g_arg, NextCand, NextPart, n, tv)) {
                                        //if (TargetCellExpPath(tv->graph, NextCand, NextPart, n, tv)) {
                                        if (TargetCellExpPath(NextCand, NextPart, tv)) {
//                                            printf("ES4; ");
                                            ExperimentalStep(NextPart, NextCand, tv, ti, m, n);
                                            PRINT_EXPPATHSTEP(NextCand, tv->answ)
                                            PRINTF2("CS2 5?: finalnumcells: %d\n", tv->finalnumcells);
                                            if (NextPart->cells == tv->finalnumcells) {
                                                UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                            }
                                            
                                        }
                                        if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                        if (!tv->answ) {
                                            PRINT_RETURN
                                        }
                                    }
                            }
                        }
                        if (tv->answ) {
                            PRINTF2("CS2 6?: finalnumcells: %d\n", tv->finalnumcells);
                            if (NextPart->cells == tv->finalnumcells) {
                                if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                temp = (tv->tolevel_tl == tv->tolevel+1);
//                                printf("CFA-9931\n");
                                autom = CheckForAutomorphisms(CurrCand, NextCand, 
                                                              tv, ti, temp, n, NextPart);
                                if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                if (autom) {
                                    for (i=autom; i<=tv->maxtreelevel; i++) {
                                        AuxCand = Spine[i].liststart;
                                        while (AuxCand && Prefix(AuxCand, NextCand, autom)) {
                                            AuxCand->do_it = FALSE;
                                            AuxCand = AuxCand->next;
                                        }
                                    }
                                    if (autom == tv->tolevel) {
                                        autom = 0;
                                    }
                                }
                                
                                if (ti->thegrouphaschanged) {
                                    ORBITSIZES
                                }
                            }
                            PRINT_RETURN
                            
                            /* ADD TO NEXT LEVEL */
                            PRINTF2("CS2 7?: finalnumcells: %d\n", tv->finalnumcells);
//                            if ((NextPart->cells < tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                            if ((NextPart->cells != tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                                ADDTONEXTLEVEL;
                                searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                            }
                        }
                        else {
                            tv->stats->interrupted++;
                        }
                        if (autom) {
                            k = n;
                            autom = 0;
                            break;
                        }
                        if (tv->fromlevel == tv->maxtreelevel) {
                            k = n;
                            break;
                        }
                    } /* end for */
                    TreeNode = RefPath[tv->maxtreelevel];
                    //search_vtx = TreeNode->vtx;
                }
            }
            else SpineTL = &Spine[tv->tolevel+1];
        }
        
    }
    
	/* REMOVE CURRENT CANDIDATE */
	if (!CurrCand->do_it || k >= SpineTL->tgtend) {
		SpineFL->liststart = CurrCand->next;
		if (CurrCand->next == NULL) {
			SpineFL->listend = NULL;
		}
		CurrCand->next = GarbList;
		GarbList = CurrCand;
	}
    //    ti->thereisnextlevel = SelectNextLevel(g_arg, n, tv);
    ti->thereisnextlevel = SelectNextLevel(n, tv);
    return 0;
}

void grouporderplus(sparsegraph *sg_orig, Candidate *Cand, Partition *Part, permnode **ring, 
                    double *grpsize1, int *grpsize2, int n, TracesVars *tv, TracesInfo *ti) {
    
	int i, i1, j, j0, j2, k, k1, k2, w, w1, w2, c, c1, c2, n1, n2;
    int prev, step, start, counts, StInd, CyInd, cycnum; //, orbrep, cellord;
	int tmp, temp, halfsize, nghcell, numvertices;
    int arg, val;
//    int *cylab, *cycol;
    //	TracesSpine SpineSDegLev;
    searchtrie *TrieNode;
    int NSFCInd, ind;
    boolean do_ngh = FALSE;
//    printf("Cand-Part grouporderplus\n");
//    PrintVect(Cand->lab, 0, n, labelorg);
//    PrintVect(Cand->invlab, 0, n, 0);
//    PrintVect(Part->cls, 0, n, 0);
//    PrintVect(Part->inv, 0, n, 0);

    // IL CASO numvertices-1 DEVE ESSERE TRATTATO SEPARATAMENTE DAL CASO deg=0 E DEVE DIVENTARE IL PRIMO CASO
//    PrintVect(sg->d, 0, n, 0);
//    PrintVect(sg_orig->d, 0, n, 0);
    
    numvertices = n;
    //    VerifyCand(Cand, n, 9923);
    //    return;
    memcpy(CanonIndices, IDENTITY_PERM, n*sizeof(int));
    memset(TreeNodes, 0, n*sizeof(int));
    
    TrieNode = Spine[tv->maxtreelevel].liststart->stnode;
    if (tv->options->verbosity >= 2) {
        fprintf(outfile, "-->> ");
        while (TrieNode->father) {
            fprintf(outfile, "%d (%d), ", TrieNode->name, TrieNode->index);
            if (TrieNode->father->name) {
                MULTIPLY(tv->stats->grpsize1, tv->stats->grpsize2, TrieNode->index);
            }
            else {
                temp = spinelementorbsize(tv->orbits, Spine[tv->maxtreelevel].liststart->lab+Spine[1].tgtcell, Spine[1].tgtsize, TrieNode->vtx);
                MULTIPLY(*grpsize1, *grpsize2, temp);
                fprintf(outfile, "orbcount vtx %d from cell %d (%d): %d\n", 
                        TrieNode->vtx+labelorg, Spine[1].tgtcell, Spine[1].tgtsize, temp);
            }
            TrieNode = TrieNode->father;
        }
    }
    else {
        while (TrieNode->father) {
            if (TrieNode->father->name) {
                MULTIPLY(tv->stats->grpsize1, tv->stats->grpsize2, TrieNode->index);
            }
            else {
                temp = spinelementorbsize(tv->orbits, Spine[tv->maxtreelevel].liststart->lab+Spine[1].tgtcell, Spine[1].tgtsize, TrieNode->vtx);
                MULTIPLY(*grpsize1, *grpsize2, temp);
            }
            TrieNode = TrieNode->father;
        }
    }
    //    printf("grouporderplus 1?(<n): %d\n", tv->smalldeglevel);
//    printf("grouporderplus cells: %d\n", Part->cells);
    if (Part->cells < n) {
        
        if (!ti->deg_one) {
            memcpy(tv->graph->e, sg_orig->e, tv->graph->elen*sizeof(int));
            for (i=0; i<n; i++) {
                TheGraph[i].e = tv->graph->e + sg_orig->v[i];
//                PrintVect(TheGraph[i].e, 0, TheGraph[i].d, labelorg);
            }
        }
//        PrintVect(sg_orig->e, 0, tv->graph->elen, labelorg);
//        PrintVect(tv->graph->e, 0, tv->graph->elen, labelorg);
        
        //    if (tv->smalldeglevel < n) {
        //        SpineSDegLev = Spine[tv->smalldeglevel];
        NSFCInd = 0;
        
        
        //        printf("A0\n");
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 9961);
        
        //        counts=0;
        //        //PrintPartition(Cand->lab, Part->cls, n, labelorg, 9338);
        //        for (i=0; i<n; i += Part->cls[i]) {
        //            tmp = Cand->lab[i];
        //            if ((sg->d[tmp] > 2) && (Part->cls[i] > 1)) {
        //                printf("%d ", tmp+labelorg);
        //                counts++;
        //            }
        //        }
        //        printf("Numero vertici: %d (+ %d = %d)\n", counts, Part->cells, counts+Part->cells);
        
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10017);
        
        //        memcpy(sg->e, sg_orig->e, sg->nde*sizeof(int));
        //        for (i = 0; i < n; i++) {
        //            SortByCell(sg->e+sg->v[i], 0, sg_orig->d[i], Cand->invlab);
        //        }
        //        printf("A2\n");
        //        VerifyCand(Cand, n, 9980);
        //        printf("Before MakeCanTree\n");
        //        PrintVect(Cand->lab, 0, n, 0);
        //        PrintVect(Cand->invlab, 0, n, 0);
        /* Trees */
        //        if (tv->options->getcanon && tv->preprocessed) {
        //            //memcpy(CanonIndices, IDENTITY_PERM, n*sizeof(int));
        //            //PrintPartition(Cand->lab, Part->cls, n, labelorg, 9848);
        //            //            printf("PART->INV: ");
        //            //            PrintVect(Part->inv, 0, n, 0);
        //            //printf("INV: ");
        //            //PrintVect(Part->inv, 0, n, 0);
        //            //            for (i=0; i<n; i++) {
        //            //                tmp = Cand->lab[i];
        //            //                printf("%d: %d < %d\n", tmp+labelorg, sg->d[tmp], sg_orig->d[tmp]);
        //            //            }
        //            //printf("treestarts: [ ");
        //            //            for (ind = 0; ind < NSFCInd; ind++) {
        //            //                i = NSFCells[ind];
        //            for (i = 0; i < n; i += Part->cls[i]) {
        //                //tmp = Cand->lab[i];
        //                tmp = Cand->lab[i];
        //                if ((sg->d[tmp] >= 0) && (sg->d[tmp] < sg_orig->d[tmp])) {
        //                    for (j=i; j<i+Part->cls[i]; j++) {
        //                        //                        printf("MCT(%d)\n", Cand->lab[j]+labelorg);
        //                        //MakeCanTree(Cand->lab[j], sg_orig, sg, n, Cand, Part, tv);
        //                        MakeCanTree(Cand->lab[j], sg_orig, sg, n, Cand, Part, tv);
        //                    }
        //                }
        //            }
        //            //printf("]\n");
        //        }
//        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10018);
        /* Trees */
        if (tv->options->getcanon && tv->preprocessed) {
            for (i = 0; i < n; i += Part->cls[i]) {
                if (Part->cls[i] == 1) {
                    tmp = Cand->lab[i];
//                    if ((sg->d[tmp] >= 0) && (sg->d[tmp] < sg_orig->d[tmp])) {
                    if ((TheGraph[tmp].d >= 0) && (TheGraph[tmp].d < sg_orig->d[tmp])) {
                        for (j=i; j<i+Part->cls[i]; j++) {
                            MakeCanTree(Cand->lab[j], sg_orig, n, Cand, Part, tv);
                        }
                    }
                }
            }
        }
        
        memset(Singletons, 0, n*sizeof(int));
//        printf("BEFORE BEFORE\n");
//        for (i=0; i<n; i++) {
//            printf("%2d: %2d ", i+labelorg, TheGraph[i].d);
//            PrintVect(TheGraph[i].e, 0, sg_orig->d[i], labelorg);
//        }
        for (i = 0; i < n; i += Part->cls[i]) {
            if (Part->cls[i] > 1) {
                if (TheGraph[Cand->lab[i]].d > 2) {
                    for (j=i; j<i+Part->cls[i]; j++) {
                        Singletons[Cand->lab[j]] = 2;
                    }
                }
//                NSFCells[NSFCInd++] = i;
            }
            else {
                Singletons[Cand->lab[i]] = 1;
//                numvertices--;
            }
        }

//        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10031);
//        PrintVect(Singletons, 0, n, 0);
//        printf("BEFORE\n");
//        for (i=0; i<n; i++) {
//            printf("%2d: %2d ", i+labelorg, TheGraph[i].d);
//            PrintVect(TheGraph[i].e, 0, sg_orig->d[i], labelorg);
//        }

        for (i = 0; i < n; i += Part->cls[i]) {
//            if (TheGraph[Cand->lab[i]].d > 0) {
                if (Part->cls[i] > 1) {
//                printf("NsDp(%d)\n", Cand->lab[i]+labelorg);
                  if (TheGraph[Cand->lab[i]].d > 2) NonSingDegPlus1(Cand, Part, i, tv);
                    NSFCells[NSFCInd++] = i;
                }
                else {
                    NonSingDegPlus2(Cand, Part, i, tv);
                    numvertices--;
                }
//            }
        }

//        printf("AFTER\n");
//        for (i=0; i<n; i++) {
//            printf("%2d: %2d ", i+labelorg, TheGraph[i].d);
//            PrintVect(TheGraph[i].e, 0, sg_orig->d[i], labelorg);
//        }
//        return;

        //        for (ind = 0; ind < NSFCInd; ind++) {
//                    i = NSFCells[ind];
////                    PrintVect(Cand->lab+i, 0, Part->cls[i], labelorg);
//                }
        
        //        VerifyCand(Cand, n, 9993);
        /* Degree 2 and at least one nghb with deg > 2 */
        //        printf("A\n");
        SETMARK(StackMarkers, tv->stackmark)
        for (ind = 0; ind < NSFCInd; ind++) {
            i = NSFCells[ind];
            SETMARK(Markers, tv->mark)
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
//                if ((sg->d[tmp] == 2) && ((sg->d[sg->e[sg->v[tmp]]] > 2) || ((sg->d[sg->e[sg->v[tmp]+1]] > 2)))) {
                if ((TheGraph[tmp].d == 2) && ((TheGraph[TheGraph[tmp].e[0]].d > 2) || ((TheGraph[TheGraph[tmp].e[1]].d > 2)))) {
                    //                    printf("2>2: ");
//                    n1 = sg->e[sg->v[tmp]];
//                    n2 = sg->e[sg->v[tmp]+1];
                    n1 = TheGraph[tmp].e[0];
                    n2 = TheGraph[tmp].e[1];
//                    if (sg->d[n1] > 2) {
                    if (TheGraph[n1].d > 2) {
//                        if (sg->d[n2] > 2) {
                        if (TheGraph[n2].d > 2) {
                            if (Cand->invlab[n1] < Cand->invlab[n2]) {
                                start = n1;
                            }
                            else {
                                start = n2;
                            }
                        }
                        else {
                            start = n1;
                        }
                    }
                    else {
                        start = n2;
                    }
                    counts = 0;
                    StInd = 0;
                    for (j=i; j<i+Part->cls[i]; j++) {
                        step = Cand->lab[j];
                        if (Markers[step] != tv->mark) {
                            prev = start;
                            counts++;
                            do {
                                Markers[step] = tv->mark;
                                PERMSTACK[StInd++] = step;
//                                if (sg->e[sg->v[step]] != prev) {
                                if (TheGraph[step].e[0] != prev) {
                                    prev = step;
//                                    step = sg->e[sg->v[step]];
                                    step = TheGraph[step].e[0];
                                }
                                else {
                                    prev = step;
//                                    step = sg->e[sg->v[step]+1];
                                    step = TheGraph[step].e[1];
                                }
//                            } while (sg->d[step] == 2);
                            } while (TheGraph[step].d == 2);
//                            if (sg->d[step] == 1) {
                            if (TheGraph[step].d == 1) {
                                PERMSTACK[StInd++] = step;
                            }
                        }
                    }
                    
                    //                    printf("counts: %d, PERMSTACK: ", counts);
                    //                    PrintVect(PERMSTACK, 0, StInd, labelorg);
                    
                    if (counts == Part->cls[i]) {
                        factorial(grpsize1, grpsize2, Part->cls[i]);
                        //                        if (tv->build_autom) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        for (k=0; k<StInd/counts; k++) {
                            i1 = PERMSTACK[k];
                            for (j0=0; j0<counts-1; j0++) {
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                            }
                            SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                        }
                        SPECIALGENERATORS
                        if (counts > 2) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                    }
                    else {
                        factorial2(grpsize1, grpsize2, Part->cls[i]);
                        //if (tv->build_autom) {
                        for (j=0; j<counts; j++) {
                            j0 = j*(StInd/counts);
                            k1 = (j+1)*(StInd/counts);
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                SETPAIRSAUTANDTREE(PERMSTACK[k], PERMSTACK[i1])
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 1) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<counts-1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 2) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                        }
                    }
                    //                    PrintPartition(Cand->lab, Part->cls, n, labelorg, 10272);
                    for (j=0; j<StInd; j++) {
                        //printf("Place %d\n", PERMSTACK[j]+labelorg);
                        Place(PERMSTACK[j], Cand, Part);
//                        if ((sg->d[PERMSTACK[j]] >= 0) && (sg->d[PERMSTACK[j]] < sg_orig->d[PERMSTACK[j]])) {
                        if ((TheGraph[PERMSTACK[j]].d >= 0) && (TheGraph[PERMSTACK[j]].d < sg_orig->d[PERMSTACK[j]])) {
                            MakeCanTree(PERMSTACK[j], sg_orig, n, Cand, Part, tv);
                        }
                    }
                    //                    PrintPartition(Cand->lab, Part->cls, n, labelorg, 10280);
                }
            }
        }
        
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10285);
        
        //        VerifyCand(Cand, n, 10221);
        
        /* Degree 2 and at least one nghb with == 1 */
        //        printf("B\n");
        for (ind = 0; ind < NSFCInd; ind++) {
            //            printf("##### ");
            //            PrintVect(Cand->lab, Part->inv[Cand->invlab[1629595]], 
            //                      Part->inv[Cand->invlab[1629595]]+Part->cls[Part->inv[Cand->invlab[1629595]]], labelorg);
            SETMARK(Markers, tv->mark)
            i = NSFCells[ind];
            //            printf("Cell: ");
            //            PrintVect(Cand->lab, i, i+Part->cls[i], labelorg);
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
//                if ((sg->d[tmp] == 2) && ((sg->d[sg->e[sg->v[tmp]]] == 1) || ((sg->d[sg->e[sg->v[tmp]+1]] == 1)))) {
                if ((TheGraph[tmp].d == 2) && ((TheGraph[TheGraph[tmp].e[0]].d == 1) || ((TheGraph[TheGraph[tmp].e[1]].d == 1)))) {
                    //                    printf("2-1: ");
                    counts = 0;
                    StInd = 0;
                    for (j=i; j<i+Part->cls[i]; j++) {
                        step = Cand->lab[j];
                        if (Markers[step] != tv->mark) {
//                            n1 = sg->e[sg->v[step]];
//                            n2 = sg->e[sg->v[step]+1];
                            n1 = TheGraph[step].e[0];
                            n2 = TheGraph[step].e[1];
                            //                            printf("n1: %d (deg: %d); ", n1+labelorg, sg->d[n1]);
                            //                            PrintVect(sg->e+sg->v[n1], 0, sg->d[n1], labelorg);
                            //                            printf("n2: %d (deg: %d); ", n2+labelorg, sg->d[n2]);
                            //                            PrintVect(sg->e+sg->v[n2], 0, sg->d[n2], labelorg);
//                            if (sg->d[n1] == 1) {
//                                if (sg->d[n2] == 1) {
                            if (TheGraph[n1].d == 1) {
                                if (TheGraph[n2].d == 1) {
                                    if (Cand->invlab[n1] < Cand->invlab[n2]) {
                                        start = n1;
                                    }
                                    else {
                                        start = n2;
                                    }
                                }
                                else {
                                    start = n1;
                                }
                            }
                            else {
                                start = n2;
                            }
                            PERMSTACK[StInd++] = start;
                            prev = start;
                            counts++;
//                            do {
//                                Markers[step] = tv->mark;
//                                PERMSTACK[StInd++] = step;
//                                if (sg->e[sg->v[step]] != prev) {
//                                    prev = step;
//                                    step = sg->e[sg->v[step]];
//                                } else {
//                                    prev = step;
//                                    step = sg->e[sg->v[step]+1];
//                                }
//                            } while (sg->d[step] == 2);
                            do {
                                Markers[step] = tv->mark;
                                PERMSTACK[StInd++] = step;
                                if (TheGraph[step].e[0] != prev) {
                                    prev = step;
                                    step = TheGraph[step].e[0];
                                } else {
                                    prev = step;
                                    step = TheGraph[step].e[1];
                                }
                            } while (TheGraph[step].d == 2);
                            PERMSTACK[StInd++] = step;
                        }
                    }
                    
                    //                    printf("counts: %d; PERMSTACK: ", counts);
                    //                    PrintVect(PERMSTACK, 0, StInd, labelorg);
                    
                    if (counts == Part->cls[i]) {
                        if (Part->inv[Cand->invlab[PERMSTACK[0]]] != Part->inv[Cand->invlab[PERMSTACK[StInd/counts-1]]]) {
                            factorial(grpsize1, grpsize2, Part->cls[i]);
                        }
                        else {
                            factorial2(grpsize1, grpsize2, 2*Part->cls[i]);
                            for (j=0; j<counts; j++) {
                                //       if (tv->build_autom) {
                                j0 = j*(StInd/counts);
                                k1 = (j+1)*(StInd/counts);
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[k], PERMSTACK[i1])
                                }
                                SPECIALGENERATORS
                                //           }
                            }
                        }
                        //if (tv->build_autom) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        for (k=0; k<StInd/counts; k++) {
                            i1 = PERMSTACK[k];
                            for (j0=0; j0<counts-1; j0++) {
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                            }
                            SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                        }
                        SPECIALGENERATORS
                        //}
                        if (counts > 2) {
                            //if (tv->build_autom) {
                            //if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                            //}
                        }
                    }
                    else {
                        factorial2(grpsize1, grpsize2, Part->cls[i]);
                        for (j=0; j<counts; j++) {
                            //if (tv->build_autom) {
                            j0 = j*(StInd/counts);
                            k1 = (j+1)*(StInd/counts);
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                SETPAIRSAUTANDTREE(PERMSTACK[k], PERMSTACK[i1])
                            }
                            SPECIALGENERATORS
                            //}
                        }
                        if (counts > 1) {
                            // if (tv->build_autom) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<counts-1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                            //  }
                        }
                        if (counts > 2) {
                            //if (tv->build_autom) {
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<1; j0++) {
                                    SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], PERMSTACK[(j0+1)*(StInd/counts)+k])
                                }
                                SETPAIRSAUTANDTREE(PERMSTACK[j0*(StInd/counts)+k], i1)
                            }
                            SPECIALGENERATORS
                            // }
                        }
                    }
                    //                    for (j=0; j<StInd; j++) {
                    //                        Place(PERMSTACK[j], Cand, Part);
                    //                    }
                    for (j=0; j<StInd; j++) {
                        //printf("Place %d\n", PERMSTACK[j]+labelorg);
                        Place(PERMSTACK[j], Cand, Part);
//                        if ((sg->d[PERMSTACK[j]] >= 0) && (sg->d[PERMSTACK[j]] < sg_orig->d[PERMSTACK[j]])) {
                        if ((TheGraph[PERMSTACK[j]].d >= 0) && (TheGraph[PERMSTACK[j]].d < sg_orig->d[PERMSTACK[j]])) {
                            MakeCanTree(PERMSTACK[j], sg_orig, n, Cand, Part, tv);
                        }
                    }
                    
                    //                    i1 = StInd/counts;
                    //                    if (Part->inv[Cand->invlab[PERMSTACK[0]]] != Part->inv[Cand->invlab[PERMSTACK[StInd/counts-1]]]) {
                    //                        for (j=0; j<i1; j++) {
                    //                            orbrep = n;
                    //                            for (k=0; k<StInd; k+=i1) {
                    //                                tmp = PERMSTACK[k+j];
                    //                                if (tv->orbits[tmp] < orbrep) {
                    //                                    orbrep = tv->orbits[tmp];
                    //                                }
                    //                            }
                    //                            cellord = Part->inv[Cand->invlab[PERMSTACK[j]]];
                    //                            StackMarkers[orbrep] = tv->stackmark;
                    //                            for (k=0; k<StInd; k+=i1) {
                    //                                tmp = PERMSTACK[k+j];
                    //                                Cand->lab[cellord] = tmp;
                    //                                Cand->invlab[tmp] = cellord++;
                    //                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                    //                                    tv->stats->numorbits--;
                    //                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                    //                                }
                    //                                tv->orbits[tmp] = orbrep;
                    //                            }
                    //                        }
                    //                    }
                    //                    else {
                    //                        for (j=0; j<i1/2; j++) {
                    //                            orbrep = n;
                    //                            for (k=0, k1=i1-1; k<StInd; k+=i1, k1+=i1) {
                    //                                tmp = PERMSTACK[k+j];
                    //                                temp = PERMSTACK[k1-j];
                    //                                if (tv->orbits[tmp] < orbrep) {
                    //                                    orbrep = tv->orbits[tmp];
                    //                                }
                    //                                if (tv->orbits[temp] < orbrep) {
                    //                                    orbrep = tv->orbits[temp];
                    //                                }
                    //                            }
                    //                            cellord = Part->inv[Cand->invlab[PERMSTACK[j]]];
                    //                            StackMarkers[orbrep] = tv->stackmark;
                    //                            for (k=0, k1=i1-1; k<StInd; k+=i1, k1+=i1) {
                    //                                tmp = PERMSTACK[k+j];
                    //                                Cand->lab[cellord] = tmp;
                    //                                Cand->invlab[tmp] = cellord++;
                    //                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                    //                                    tv->stats->numorbits--;
                    //                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                    //                                }
                    //                                tv->orbits[tmp] = orbrep;
                    //
                    //                                temp = PERMSTACK[k1-j];
                    //                                Cand->lab[cellord] = temp;
                    //                                Cand->invlab[temp] = cellord++;
                    //                                if (StackMarkers[tv->orbits[temp]] != tv->stackmark) {
                    //                                    tv->stats->numorbits--;
                    //                                    StackMarkers[tv->orbits[temp]] = tv->stackmark;
                    //                                }
                    //                                tv->orbits[temp] = orbrep;
                    //                            }
                    //                        }
                    //                        if (i1 % 2) {
                    //                            orbrep = n;
                    //                            for (k=j; k<StInd; k+=i1) {
                    //                                tmp = PERMSTACK[k];
                    //                                if (tv->orbits[tmp] < orbrep) {
                    //                                    orbrep = tv->orbits[tmp];
                    //                                }
                    //                            }
                    //                            cellord = Part->inv[Cand->invlab[PERMSTACK[j]]];
                    //                            StackMarkers[orbrep] = tv->stackmark;
                    //                            j0 = 0;
                    //                            for (k=j; k<StInd; k+=i1) {
                    //                                tmp = PERMSTACK[k];
                    //                                Cand->lab[cellord+j0] = tmp;
                    //                                Cand->invlab[tmp] = cellord+j0++;
                    //                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                    //                                    tv->stats->numorbits--;
                    //                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                    //                                }
                    //                                tv->orbits[tmp] = orbrep;
                    //                            }
                    //                        }
                    //                    }
                    //                    Propagate(sg, Cand, Part, i, i+Part->cls[i], n);
                }
            }
        }
        
        //        VerifyCand(Cand, n, 10459);
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10526);
        
        /* Cycles */
        //        printf("C\n");
        for (ind = 0; ind < NSFCInd; ind++) {
            //            printf("##### ");
            //            PrintVect(Cand->lab, Part->inv[Cand->invlab[1629595]], 
            //                      Part->inv[Cand->invlab[1629595]]+Part->cls[Part->inv[Cand->invlab[1629595]]], labelorg);
            i = NSFCells[ind];
            //for (i = 0; i < n; i += Part->cls[i]) {
            SETMARK(Markers, tv->mark)
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
//                if (sg->d[tmp] == 2) {
                if (TheGraph[tmp].d == 2) {
//                    printf("2cy: ");
//                    PrintVect(tv->orbits, 0, n, labelorg);
                    //                    printf("%d\n", tmp+labelorg);
                    //                    PrintVect(Cand->lab, i, i+Part->cls[i], labelorg);
                    CyInd = StInd = cycnum = 0;
                    for (j=i; j<i+Part->cls[i]; j++) {
                        start = Cand->lab[j];
                        if (Markers[start] != tv->mark) {
                            counts = 1;
                            CYCLES[StInd] = start;
                            CYCOLR[StInd++] = Part->inv[Cand->invlab[start]];
                            Markers[start] = tv->mark;
//                            k = Cand->invlab[sg->e[sg->v[start]]];
//                            k1 = Cand->invlab[sg->e[sg->v[start]+1]];
//                            if (Part->inv[k] < Part->inv[k1]) {
//                                step = sg->e[sg->v[start]];
//                            }
//                            else {
//                                step = sg->e[sg->v[start]+1];
//                            }
                            k = Cand->invlab[TheGraph[start].e[0]];
                            k1 = Cand->invlab[TheGraph[start].e[1]];
                            if (Part->inv[k] < Part->inv[k1]) {
                                step = TheGraph[start].e[0];
                            }
                            else {
                                step = TheGraph[start].e[1];
                            }
                            prev = start;
                            do {
                                counts++;
                                Markers[step] = tv->mark;
                                CYCLES[StInd] = step;
                                CYCOLR[StInd++] = Part->inv[Cand->invlab[step]];
                                //                                if (sg->e[sg->v[step]] != prev) {
                                //                                    prev = step;
                                //                                    step = sg->e[sg->v[step]];
                                //                                }
                                //                                else {
                                //                                    prev = step;
                                //                                    step = sg->e[sg->v[step]+1];
                                //                                }
//                                printf("step: %d; start: %d\n", step+labelorg, start+labelorg);
                                //                                printf("D %7d: ", step+labelorg);
                                //PrintVect(TheGraph[step].e, 0, sg_orig->d[step], labelorg);
 
                                if (TheGraph[step].e[0] != prev) {
                                    prev = step;
                                    step = TheGraph[step].e[0];
//                                    if (TheGraph[step].d == 0) {
//                                        printf("DEG(%d) is 0\n", step+labelorg);
//                                    }
                                }
                                else {
                                    prev = step;
                                    step = TheGraph[step].e[1];
//                                    if (TheGraph[step].d == 0) {
//                                        printf("DEG(%d) is 0\n", step+labelorg);
//                                    }
                                }
                            } while (step != start);
                            CYLGTH[CyInd++] = counts;
                            cycnum++;
                        }
                    }
                    
                    CYCPOS[0] = 0;
                    for (j=1; j<CyInd; j++) {
                        CYCPOS[j] = CYCPOS[j-1]+CYLGTH[j-1];
                    }
                    memcpy(WorkArray, CYLGTH, CyInd*sizeof(int));
                    sort2ints(WorkArray, CYCPOS, CyInd);
                    
                    k = 0;
                    for (i1=0; i1<CyInd; i1++) {
                        k1 = CYCOLR[k];
                        k2 = CYCOLR[k+1];
                        for (j=1; j<=CYLGTH[i1]/2; j++) {
                            w1 = CYCOLR[j+k];
                            w2 = CYCOLR[j+1+k];
                            if ((w1 == k1) && (w2 == k2)) {
                                //                                if (tv->build_autom) {
                                //                                    printf("AA1 ");
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (w=0; w<CYLGTH[i1]; w++) {
                                    if (CYCOLR[w+k] == CYCOLR[((w+j) % CYLGTH[i1]) + k]) {
                                        SETPAIRSAUTANDTREE(CYCLES[w+k], CYCLES[((w+j) % CYLGTH[i1]) + k])
                                    }
                                    else {
                                        break;
                                    }
                                }
                                if (w == CYLGTH[i1]) { SPECIALGENERATORS }
                                if (w == CYLGTH[i1]) {
                                    MULTIPLY(*grpsize1, *grpsize2, CYLGTH[i1]/j);
                                    break;
                                }
                            }
                        }
                        
//                        PrintVect(CYCLES, 0, StInd, labelorg);
//                        PrintVect(CYCOLR, 0, StInd, labelorg);
//                        PrintVect(CYLGTH, 0, CyInd, 0);
                        //                        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10953);
                        
                        //                        printf("Part->cls[%d]: %d; Part->cls[%d]: %d\n", k1, Part->cls[k1], k2, Part->cls[k2]);
                        if (Part->cls[k1] >= Part->cls[k2]) {
                            for (j=CYLGTH[i1]-1; j>0; j--) {
                                w1 = CYCOLR[j % CYLGTH[i1] + k];
                                w2 = CYCOLR[(j-1) % CYLGTH[i1] + k];
                                //                                printf("AA w1: %d, w2: %d; k1: %d; k2: %d\n", w1, w2, k1, k2);
                                if ((w1 == k1) && (w2 == k2)) {
                                    //                                    if (tv->build_autom) {
                                    //                                        printf("AA2 ");
                                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                    for (w=0; w<CYLGTH[i1]; w++) {
                                        SETPAIRSAUTANDTREE(CYCLES[w+k], CYCLES[((j-w+(w>j)*CYLGTH[i1]) % CYLGTH[i1]) + k])
                                    }
                                    SPECIALGENERATORS
                                    MULTIPLY(*grpsize1, *grpsize2, 2);
                                    break;
                                }
                            }
                        }
                        else {
                            j=CYLGTH[i1]-1;
                            w2 = CYCOLR[j % CYLGTH[i1] + k];
                            //                                printf("BB w1: %d, w2: %d; k1: %d; k2: %d\n", w1, w2, k1, k2);
                            if (w2 == k2) {
                                //                                    if (tv->build_autom) {
                                //                                        printf("AA2 ");
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (w=1; w<CYLGTH[i1]; w++) {
                                    SETPAIRSAUTANDTREE(CYCLES[w+k], CYCLES[CYLGTH[i1]-w+k])
                                }
                                SPECIALGENERATORS
                                MULTIPLY(*grpsize1, *grpsize2, 2);
                                //break;
                            }
                        }
                        k += CYLGTH[i1];
                    }
                    k = 0;
                    for (i1=0; i1<CyInd; i1++) {
                        if (CYLGTH[i1] > 0) {
                            CYMULT[0] = k;
                            k1 = k;
                            counts = 1;
                            for (j0=i1+1; j0<CyInd; j0++) {
                                k1 += abs(CYLGTH[j0]);
                                if (CYLGTH[j0] == CYLGTH[i1]) {
                                    CYMULT[counts++] = k1;
                                    CYLGTH[j0] = -CYLGTH[j0];
                                }
                            }
                            if (counts > 1) {
                                //                                if (tv->build_autom) {
                                //                                    printf("BBB ");
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                for (j0=0; j0<CYLGTH[i1]; j0++) {
                                    for (j2 = 0; j2<counts-1; j2++) {
                                        SETPAIRSAUTANDTREE(CYCLES[CYMULT[j2]+j0], CYCLES[CYMULT[j2+1]+j0])
                                    }
                                    SETPAIRSAUTANDTREE(CYCLES[CYMULT[j2]+j0], CYCLES[CYMULT[0]+j0])
                                }
                                SPECIALGENERATORS
                                if (counts > 2) {
                                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                    //VerifyId(AUTPERM, n);
                                    for (j0=0; j0<CYLGTH[i1]; j0++) {
                                        SETPAIRSAUTANDTREE(CYCLES[CYMULT[1]+j0], CYCLES[CYMULT[0]+j0])
                                        if (tv->build_autom) {
                                            SETPAIRSAUT(CYCLES[CYMULT[0]+j0], CYCLES[CYMULT[1]+j0])
                                        }
                                        MakeTree(CYCLES[CYMULT[0]+j0], CYCLES[CYMULT[1]+j0], sg_orig, n, tv, FALSE);
                                    }
                                    SPECIALGENERATORS
                                }
                                factorial(grpsize1, grpsize2, counts);
                            }
                        }
                        k += abs(CYLGTH[i1]);
                        CYLGTH[i1] = -CYLGTH[i1];
                    }
                    
                    for (c1=0; c1<CyInd; c1++) {
                        c = CYCPOS[c1]+WorkArray[c1];
                        for (c2=CYCPOS[c1]; c2<c; c2++) {
                            //                            printf("%d ", CYCLES[c2]);
                            Place(CYCLES[c2], Cand, Part);
//                            if ((sg->d[CYCLES[c2]] >= 0) && (sg->d[CYCLES[c2]] < sg_orig->d[CYCLES[c2]])) {
                            if ((TheGraph[CYCLES[c2]].d >= 0) && (TheGraph[CYCLES[c2]].d < sg_orig->d[CYCLES[c2]])) {
                                MakeCanTree(CYCLES[c2], sg_orig, n, Cand, Part, tv);
                            }
                            
                        }
                    }
                    //                    PrintPartition(Cand->lab, Part->cls, n, labelorg, 10838);
                    
                    //Propagate(sg, Cand, Part, i, i+Part->cls[i], n);
                }
            }
        }
//        PrintVect(tv->orbits, 0, n, labelorg);

        //        VerifyCand(Cand, n, 10725);
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10793);
        
        /* Degree 1, and nghb too */
        //        printf("D\n");
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10797);
        SETMARK(Markers, tv->mark)
        for (ind = 0; ind < NSFCInd; ind++) {
            //            printf("##### ");
            //            PrintVect(Cand->lab, Part->inv[Cand->invlab[1629595]], 
            //                      Part->inv[Cand->invlab[1629595]]+Part->cls[Part->inv[Cand->invlab[1629595]]], labelorg);
            i = NSFCells[ind];
            //for (i = 0; i < n; i += Part->cls[i]) {
            //            PrintVect(Cand->lab+i, 0, Part->cls[i], labelorg);
			if (Part->cls[i] > 1) {
				tmp = Cand->lab[i];
//                if (Part->cls[i] > 1) {
//                tmp = Cand->lab[i];
//                if ((sg->d[tmp] == 1) && (sg->d[sg->e[sg->v[tmp]]] == 1) && (i == Part->inv[Cand->invlab[sg->e[sg->v[tmp]]]])) {
                if ((TheGraph[tmp].d == 1) && (TheGraph[TheGraph[tmp].e[0]].d == 1) && (i == Part->inv[Cand->invlab[TheGraph[tmp].e[0]]])) {
                    //PrintPartition(Cand->lab, Part->cls, n, labelorg, 10898);
                    //                        printf("1-1: ");
                    //                    printf("\n");
                    factorial2(grpsize1, grpsize2, Part->cls[i]);
                    /* the cell has size two */
                    if (Part->cls[i] == 2) {
                        val = Cand->lab[i+1];
                        //                            if (tv->build_autom) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        arg = tmp;
                        SETPAIRSAUTANDTREE(arg, val)
                        SETPAIRSAUTANDTREE(val, arg)
                        //                                if (tv->build_autom) {
                        //                                    MakeTree(val, tmp, sg_orig, sg, n, tv, FALSE);
                        //                                }
                        SPECIALGENERATORS
                    }
                    else {
                        /* the cell has size greater than two */
                        //                            if (tv->build_autom) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        SETMARK(Markers, tv->mark)
                        halfsize = Part->cls[i]/2;
                        i1 = 0;
                        for (j=i; j<i+Part->cls[i]; j++) {
                            if (Markers[Cand->lab[j]] != tv->mark) {
//                                Markers[sg->e[sg->v[Cand->lab[j]]]] = tv->mark;
                                Markers[TheGraph[Cand->lab[j]].e[0]] = tv->mark;
                                PERMSTACK[i1] = Cand->lab[j];
//                                PERMSTACK[i1+halfsize] = sg->e[sg->v[Cand->lab[j]]];
                                PERMSTACK[i1+halfsize] = TheGraph[Cand->lab[j]].e[0];
                                i1++;
                            }
                        }
//                        PrintVect(PERMSTACK, 0, 2*i1, labelorg);
                        temp = PERMSTACK[0];
                        for (j=0; j<Part->cls[i]-1; j++) {
                            SETPAIRSAUTANDTREE(PERMSTACK[j], PERMSTACK[j+1])
                        }
                        SETPAIRSAUTANDTREE(PERMSTACK[j], temp)
                        SPECIALGENERATORS
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        //                        for (j=halfsize; j<i+Part->cls[i]-1; j++) {
                        //                            PERMSTACK[j] = PERMSTACK[j+1];
                        //                        }
//                        printf("Halfsize: %d\n", halfsize);
//                        PrintVect(PERMSTACK, 0, 2*i1, labelorg);
                        memmove(PERMSTACK+halfsize, PERMSTACK+halfsize+1, (halfsize-1)*sizeof(int));  // was memcpy
//                        PrintVect(PERMSTACK, 0, 2*i1, labelorg);
                        temp = PERMSTACK[1];
                        for (j=1; j<Part->cls[i]-2; j++) {
                            SETPAIRSAUTANDTREE(PERMSTACK[j], PERMSTACK[j+1])
                        }
                        SETPAIRSAUTANDTREE(PERMSTACK[j], temp)
//                        printf("AAA ");
                        SPECIALGENERATORS
                    }
                    //                        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10934);
                    
                    //                        for (j=0; j<StInd; j++) {
                    //                            printf("Place %d\n", Cand->lab[j]+labelorg);
                    //                            Place(Cand->lab[j], Cand, Part);
                    //                            if ((sg->d[Cand->lab[j]] >= 0) && (sg->d[Cand->lab[j]] < sg_orig->d[Cand->lab[j]])) {
                    //                                MakeCanTree(Cand->lab[j], sg_orig, sg, n, Cand, Part, tv);
                    //                            }
                    //                        }
                    //                        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10943);
                    
                    SETMARK(Markers, tv->mark)
                    for (j=i; j<i+Part->cls[i]; j++) {
                        temp = Cand->lab[j];
                        if (Markers[temp] != tv->mark) {
//                            if ((sg->d[temp] >= 0) && (sg->d[temp] < sg_orig->d[temp])) {
                            if ((TheGraph[temp].d >= 0) && (TheGraph[temp].d < sg_orig->d[temp])) {
                                //printf("a) MCT @ %d\n", temp+labelorg);
                                MakeCanTree(temp, sg_orig, n, Cand, Part, tv);
                            }
                            tmp = Cand->lab[j+1];
//                            Markers[sg->e[sg->v[temp]]] = tv->mark;
                            Markers[TheGraph[temp].e[0]] = tv->mark;
//                            i1 = Cand->invlab[sg->e[sg->v[temp]]];
                            i1 = Cand->invlab[TheGraph[temp].e[0]];
//                            Cand->lab[j+1] = sg->e[sg->v[temp]];
                            Cand->lab[j+1] = TheGraph[temp].e[0];
//                            if ((sg->d[sg->e[sg->v[temp]]] >= 0) && (sg->d[sg->e[sg->v[temp]]] < sg_orig->d[sg->e[sg->v[temp]]])) {
                            if ((TheGraph[TheGraph[temp].e[0]].d >= 0) && (TheGraph[TheGraph[temp].e[0]].d < sg_orig->d[TheGraph[temp].e[0]])) {
                                //printf("b) MCT @ %d\n", sg->e[sg->v[temp]]+labelorg);
//                                MakeCanTree(sg->e[sg->v[temp]], sg_orig, sg, n, Cand, Part, tv);
                                MakeCanTree(TheGraph[temp].e[0], sg_orig, n, Cand, Part, tv);
                            }
//                            Cand->invlab[sg->e[sg->v[temp]]] = j+1;
                            Cand->invlab[TheGraph[temp].e[0]] = j+1;
                            Cand->lab[i1] = tmp;
                            Cand->invlab[tmp] = i1;
                        }
                    }
                    
                    //Propagate(sg, Cand, Part, i, i+Part->cls[i], n);
                }
                //                }
			}
		}
        
        //        VerifyCand(Cand, n, 10870);
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 10966);
        
        /* Degree 0 */
        //        printf("E\n");
        for (ind = 0; ind < NSFCInd; ind++) {
            //            printf("##### ");
            //            PrintVect(Cand->lab, Part->inv[Cand->invlab[1629595]], 
            //                      Part->inv[Cand->invlab[1629595]]+Part->cls[Part->inv[Cand->invlab[1629595]]], labelorg);
            i = NSFCells[ind];
            //for (i = 0; i < n; i += Part->cls[i]) {
            if (Part->cls[i] > 1) {
                tmp = Cand->lab[i];
//                printf("tmp: %d, sg->d[tmp]: %d, numvertices: %d\n", tmp+labelorg, sg->d[tmp], numvertices);
//                if (sg->e != NULL) {
//                    PrintVect(Part->inv, 0, n, 0);
//                    PrintVect(Cand->lab, 0, n, labelorg);
//                    PrintVect(Cand->invlab, 0, n, 0);
//                    printf("tmp: %d, sg->d[tmp]: %d\n", tmp, sg->d[tmp]);
//                    printf("tmp: %d, sg->e: %p, sg->e[sg->v[tmp]]: %d, exp: %d\n", 
//                           tmp+labelorg, sg->e, sg->e[sg->v[tmp]]+labelorg, Part->inv[Cand->invlab[sg->e[sg->v[tmp]]]]);
//                }
//                if ((sg->e != NULL) && (sg->d[tmp] != 0) && (sg->d[tmp] != numvertices-1))
                if ((TheGraph[0].e != NULL) && (TheGraph[tmp].d != 0) && (TheGraph[tmp].d != numvertices-1))
//                    nghcell = Part->inv[Cand->invlab[sg->e[sg->v[tmp]]]]; else nghcell = i;
                    nghcell = Part->inv[Cand->invlab[TheGraph[tmp].e[0]]]; else nghcell = i;
//                if ((sg->d[tmp] == 0) ||
//                    ((sg->d[tmp] == numvertices-1) && (sg->d[tmp] > 2)) ||
//                    ((sg->d[tmp] == 1) && (sg->d[sg->e[sg->v[tmp]]] == 1) && (i < nghcell))) {
                if ((TheGraph[tmp].d == 0) ||
                    ((TheGraph[tmp].d == numvertices-1) && (TheGraph[tmp].d > 2)) ||
                    ((TheGraph[tmp].d == 1) && (TheGraph[TheGraph[tmp].e[0]].d == 1) && (i < nghcell))) {
//                    PrintVect(Cand->lab+i, 0, Part->cls[i], labelorg);
                    //                    printf(" 0 : ");
                    do_ngh = FALSE;
//                    if ((sg->d[tmp] == 1) && (sg->d[sg->e[sg->v[tmp]]] == 1) && (i != nghcell)) {
                    if ((TheGraph[tmp].d == 1) && (TheGraph[TheGraph[tmp].e[0]].d == 1) && (i != nghcell)) {
                        //                        printf("XXX ");
                        do_ngh = TRUE;
                    }
                    //                    if (tv->build_autom) {
                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                    for (j=i; j<i+Part->cls[i]-1; j++) {
                        arg = Cand->lab[j];
                        val = Cand->lab[j+1];
                        SETPAIRSAUTANDTREE(arg, val)
                        if (do_ngh) {
//                            SETPAIRSAUTANDTREE(sg->e[sg->v[arg]], sg->e[sg->v[val]])
                            SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])

                        }
                    }
                    arg = Cand->lab[j];
                    val = tmp;
                    SETPAIRSAUTANDTREE(arg, val)
                    if (do_ngh) {
//                        SETPAIRSAUTANDTREE(sg->e[sg->v[arg]], sg->e[sg->v[val]])
                        SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])
                    }
                    SPECIALGENERATORS
                    if (Part->cls[i] > 2) {
                        //                        if (tv->build_autom) {
                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                        arg = tmp;
                        val = Cand->lab[i+1];
                        SETPAIRSAUTANDTREE(arg, val)
                        if (do_ngh) {
//                            SETPAIRSAUTANDTREE(sg->e[sg->v[arg]], sg->e[sg->v[val]])
                            SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])
                        }
                        arg = Cand->lab[i+1];
                        val = tmp;
                        SETPAIRSAUTANDTREE(arg, val)
                        if (do_ngh) {
//                            SETPAIRSAUTANDTREE(sg->e[sg->v[arg]], sg->e[sg->v[val]])
                            SETPAIRSAUTANDTREE(TheGraph[arg].e[0], TheGraph[val].e[0])
                        }
                        SPECIALGENERATORS
                    }
                    factorial(grpsize1, grpsize2, Part->cls[i]);
                    if (do_ngh) {
                        //                        PrintPartition(Cand->lab, Part->cls, n, labelorg, 11299);
                        for (j=i; j<i+Part->cls[i]; j++) {
//                            temp = sg->e[sg->v[Cand->lab[j]]];
                            temp = TheGraph[Cand->lab[j]].e[0];
                            Cand->lab[nghcell] = temp;
                            Cand->invlab[temp] = nghcell;
                            nghcell++;
                        }
                        //                        PrintPartition(Cand->lab, Part->cls, n, labelorg, 11306);
                    }
                    
                    k = i+Part->cls[i];
                    for (j=i; j<k; j++) {
                        //                        printf("Place %d\n", Cand->lab[j]+labelorg);
                        Place(Cand->lab[j], Cand, Part);
//                        if ((sg->d[Cand->lab[j]] >= 0) && (sg->d[Cand->lab[j]] < sg_orig->d[Cand->lab[j]])) {
                        if ((TheGraph[Cand->lab[j]].d >= 0) && (TheGraph[Cand->lab[j]].d < sg_orig->d[Cand->lab[j]])) {

                            MakeCanTree(Cand->lab[j], sg_orig, n, Cand, Part, tv);
                            if (do_ngh) {
                                MakeCanTree(TheGraph[Cand->lab[j]].e[0], sg_orig, n, Cand, Part, tv);
                            }
                        }
                    }
                    //                    PrintPartition(Cand->lab, Part->cls, n, labelorg, 10280);
                }
                //Propagate(sg, Cand, Part, i, i+Part->cls[i], n);
            }
        }
        
        //        PrintPartition(Cand->lab, Part->cls, n, labelorg, 11012);
        //PrintPartition(Cand->lab, Part->cls, n, labelorg, 10165);
        //        VerifyCand(Cand, n, 10971);
        //        printf("G\n");
        
        /* Trees */
        //        if (tv->options->getcanon && tv->preprocessed) {
        //            for (i = 0; i < n; i += Part->cls[i]) {
        //                tmp = Cand->lab[i];
        //                if ((sg->d[tmp] >= 0) && (sg->d[tmp] < sg_orig->d[tmp])) {
        //                    for (j=i; j<i+Part->cls[i]; j++) {
        //                        MakeCanTree(Cand->lab[j], sg_orig, sg, n, Cand, Part, tv);
        //                    }
        //                }
        //            }
        //        }
        
        
	}
    
    /* Orbit Count */
    //printf("Orbits count: %d -> ", tv->stats->numorbits);
    SETMARK(Markers, tv->mark)
    i1=0;
    for (c1=0; c1<n; c1++) {
        if (Markers[tv->orbits[c1]] != tv->mark) {
            i1++;
            Markers[tv->orbits[c1]] = tv->mark;
        }
    }
    tv->stats->numorbits = i1;
    //printf("%d\n", tv->stats->numorbits);
    return;
}

boolean Prefix(Candidate *Cand1, Candidate *Cand2, int k) {
	int i;
	for (i=1; i<=k; i++) {
		if (Cand1->lab[Spine[k].tgtpos] != Cand2->lab[Spine[k].tgtpos]) {
			break;
		}
	}
	return (i>k);
}

boolean findperm(permnode *pn, int *p, int n) {
    permnode *rn;
    if (!pn) {
        return FALSE;
    }
    rn = pn;
    do {
        if (!memcmp(rn->p, p, n*sizeof(int))) {
            return TRUE;
        }
        rn = rn->next;
    } while (rn != pn);
    return FALSE;
}

int spinelementorbsize(int *orbits, int *lab, int size, int elem) {
    int i, j, val;
    j = 0;
    val = orbits[elem];
    for (i = 0; i < size; ++i) {
        if (orbits[lab[i]] == val) ++j;
    }
    return j;
}

//void Propagate(sparsegraph *sg, Candidate *Cand, Partition *Part, int start, int end, int n) {
//    int i, vtx, n1, n2;
//    
//    for (i=start; i<end; i++) {
//        Part->cls[i] = 1;
//        Part->inv[i] = i;
//    }
//    Part->cells = Part->cells + end - start - 1;
//    vtx = Cand->lab[start];
//    if (sg->d[vtx] == 2) {
//        n1 = sg->e[sg->v[vtx]];
//        n2 = sg->e[sg->v[vtx]+1];
//        if ((Part->cls[Part->inv[Cand->invlab[n1]]] > 1) && (sg->d[n1] < n-1))
//            Propagate(sg, Cand, Part, Part->inv[Cand->invlab[n1]], Part->inv[Cand->invlab[n1]]+Part->cls[Part->inv[Cand->invlab[n1]]], n);
//        if ((Part->cls[Part->inv[Cand->invlab[n2]]] > 1) && (sg->d[n2] < n-1))
//            Propagate(sg, Cand, Part, Part->inv[Cand->invlab[n2]], Part->inv[Cand->invlab[n2]]+Part->cls[Part->inv[Cand->invlab[n2]]], n);
//        return;
//    }
//}

trielist *searchtrie_new(int n, struct TracesVars *tv) {
	tv->strielist = malloc(sizeof(struct trielist));
	if (tv->strielist == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
    tv->strielist->prev = tv->strielist->next = NULL;
    tv->strielist->triearray = malloc(n*sizeof(searchtrie));
	if (tv->strielist->triearray == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	tv->strielist->triearray[0].father = tv->strielist->triearray[0].first_child = NULL;
	tv->strielist->triearray[0].next_sibling = tv->strielist->triearray[0].last_child = NULL;
    tv->strielist->triearray[0].goes_to = NULL;
    tv->strielist->triearray[0].index = 1;
    tv->strielist->triearray[0].name = tv->strielist->triearray[0].level = 0;
    tv->strielist->triearray[0].vtx = n;
    
	tv->strienext = 1;
    return tv->strielist;
}

searchtrie *searchtrie_make(Candidate *CurrCand, Candidate *NextCand, int n, struct TracesVars *tv) {
    searchtrie *st;
	if (tv->strienext == n) {
		tv->strienext = 0;
		tv->strielist->next = malloc(sizeof(struct trielist));
		if (tv->strielist->next == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
        tv->strielist->next->prev = tv->strielist;
        tv->strielist = tv->strielist->next;
        tv->strielist->next = NULL;
        tv->strielist->triearray = malloc(n*sizeof(searchtrie));
        if (tv->strielist->triearray == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
    }
    st = &(tv->strielist->triearray[tv->strienext]);
	st->father = CurrCand->stnode;
	st->name = NextCand->name;
	st->index = tv->newindex+1;
	st->vtx = NextCand->vertex;
	st->level = tv->tolevel;
    st->first_child = st->next_sibling = st->last_child = st->goes_to = NULL;
    if (st->father) {
        if (st->father->first_child) {
            st->father->last_child->next_sibling = st;
            st->father->last_child = st;
        }
        else {
            st->father->first_child = st->father->last_child = st;
        }
    }
    NextCand->stnode = st;
    if (tv->newgotonode) {
        tv->newgotonode->goes_to = st;
    }
    if (tv->gotonode) {
        st->goes_to = tv->gotonode;
        tv->gotonode = NULL;
    }
    tv->strienext++;
    return st;
}

boolean lookup(searchtrie *t) {
    searchtrie *TreeNode;
    TreeNode = t;
    while (TreeNode->level >= 1) {
        if (TreeNode->goes_to) {
            return FALSE;
        }
        TreeNode = TreeNode->father;
    }
    return TRUE;
}

int *findcurrorbits(schreier *gp, int k) {
    int i;
    schreier *sh;
    
    sh = gp;
    for (i = 0; i < k; i++) {
        sh = sh->next;
    }
    return sh->orbits;
}

int Preprocess(sparsegraph *sg, 
                    permnode **ring, 
                    Candidate *Cand, 
                    int n, 
                    Partition *Part, 
                    struct TracesVars* tv) {
    int i, curr_cell, ind, ind0, ind1, ind2, j, j0, k; //, cell, cell1, vtx, ngh, cellsize, cell1size, cell1end, temp, sons, d_ngh, d_vtx, pos1;
    int *sge;
    int HitClsInd, labi, nghb, value, SplInd, SplCntInd, sc, iend, CStackInd, SingInd, newcell, TraceInd;
    
    
    CStackInd = 0;
    for (i = 0; i < n; i += Part->cls[i]) {
        //if (sg->d[Cand->lab[i]] == 1) {
        if (TheGraph[Cand->lab[i]].d == 1) {
            CStack[CStackInd++] = i;
        }
    }
    
    TraceInd = Part->cells;
//    printf("CStackInd: %d\n", CStackInd);
    if (CStackInd > 0) {
        ind = 0;
        while (ind < CStackInd) {
//            printf("PP-STACK: ");
//            PrintVect(CStack, ind, CStackInd, 0);
//            PrintPartition(Cand->lab, Part->cls, n, labelorg, 11142);
            if (tv->mark > (NAUTY_INFINITY-2)) {
                memset(Markers, 0, n*sizeof(int));
                memset(MarkHitVtx, 0, n*sizeof(int));
                tv->mark = 0;
            }
            tv->mark++;
            
            curr_cell = CStack[ind++];
            ind2 = curr_cell+Part->cls[curr_cell];
            HitClsInd = 0;
            for (i = curr_cell; i < ind2; i++) {
                labi = Cand->lab[i];
                nghb = *(TheGraph[labi].e);
                
                if (TheGraph[nghb].d != 1) {
                    //sg1->d[labi] = -1;
                    TheGraph[labi].d = -1;
                    TheGraph[labi].one = TRUE;
                }
                
                if (MarkHitVtx[nghb] == tv->mark) {
                    NghCounts[nghb]++;
                }
                else {
                    value = Part->inv[Cand->invlab[nghb]];
                    MarkHitVtx[nghb] = tv->mark;
                    NghCounts[nghb] = 1;
                    if (Markers[value] != tv->mark) {
                        HitCls[HitClsInd++] = value;
                        Markers[value] = tv->mark;
                        HitVtx[value] = nghb;
                        ElmHitCll[value] = 1;
                    }
                    else {
                        HitVtx[value+ElmHitCll[value]++] = nghb;
                    }
                }
            }
            
            tv->mark++;
            
            switch (HitClsInd) {
                case 0:
                case 1:
                    break;
                case 2:
                    if (HitCls[0] > HitCls[1]) {
                        value = HitCls[0];
                        HitCls[0] = HitCls[1];
                        HitCls[1] = value;
                    }
                    break;
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    for (k = 1; k < HitClsInd; ++k) {
                        value = HitCls[k];
                        i = k - 1;
                        while ((i >= 0) && (value < HitCls[i])) {
                            HitCls[i + 1] = HitCls[i];
                            --i;
                        }
                        HitCls[i + 1] = value;
                    }
                    break;
                default:
                    quickSort(HitCls, HitClsInd);
                    break;
            }

            SplInd = 0;
            SplCls[0] = n;
            for (j = 0; j < HitClsInd; j++) {
                ind1 = HitCls[j];
                if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < Part->cls[ind1])) {
                    SplCls[SplInd++] = ind1;
                }
                else {
                    ind2 = ind1+Part->cls[ind1];
                    value = NghCounts[Cand->lab[ind1++]];
                    for (i = ind1; i < ind2; i++) {
                        if (NghCounts[Cand->lab[i]] != value) {
                            SplCls[SplInd++] = HitCls[j];
                            break;
                        }
                    }
                    if (i == ind2) {
                        ind1 = HitCls[j];
                        if (TheGraph[Cand->lab[ind1]].d != 1) {   // !!!!
                            for (i = ind1; i < ind2; i++) {
                                value = Cand->lab[i];
                                Edge_Delete(value, NghCounts[value], Cand, tv);
//                                if (NghCounts[value]>1) {
//                                    printf("%d!\n", NghCounts[value]);
//                                    printf("A %d\n", value+labelorg);
//                                    writegroupsize(outfile, tv->stats->grpsize1, tv->stats->grpsize2);
//                                    printf("\n");
//                                }
                                
                                
                                sge = TheGraph[value].e+TheGraph[value].d;
                                if (NghCounts[value]>1) {
//                                    printf("A %d: ", value+labelorg);
//                                    PrintVect(sge, 0, NghCounts[value], labelorg);
//                                    
                                    factorial(&(tv->stats->grpsize1), &(tv->stats->grpsize2), NghCounts[value]);
                                    if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                    for (j0=0; j0<NghCounts[value]-1; j0++) {
                                        SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[j0+1])
//                                        printf("[%d %d]", sge[j0], sge[j0+1]);
                                    }
                                    SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[0])
//                                    printf("[%d %d]\n", sge[j0], sge[0]);
                                    SPECIALGENERATORS
                                    if (NghCounts[value] > 2) {
                                        if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                        SETPAIRSAUTANDTREE_PREPROC(sge[0], sge[1])
//                                        printf("[%d %d]", sge[0], sge[1]);
                                        if (tv->build_autom) {
                                            SETPAIRSAUT(sge[1], sge[0])
//                                            printf("[%d %d]\n", sge[1], sge[0]);
                                        }
                                        MakeTree(sge[1], sge[0], sg, n, tv, FALSE);
                                        SPECIALGENERATORS
                                    }
                                }
                            }
                            if (TheGraph[Cand->lab[ind1]].d == 1) {
                                CStack[CStackInd++] = ind1;
//                                printf("A %d\n", CStack[CStackInd-1]);
                            }
                        }
                    }
                }
            }
            
            if (SplInd) {
//                printf("SplCls:");
//                PrintVect(SplCls, 0, SplInd, 0);
                /* Sorting the cells to be split */
//                switch (SplInd) {
//                    case 0:
//                    case 1:
//                        break;
//                    case 2:
//                        if (SplCls[0] > SplCls[1]) {
//                            value = SplCls[0];
//                            SplCls[0] = SplCls[1];
//                            SplCls[1] = value;
//                        }
//                        break;
//                    case 3:
//                    case 4:
//                    case 5:
//                    case 6:
//                    case 7:
//                    case 8:
//                        for (k = 1; k < SplInd; ++k) {
//                            value = SplCls[k];
//                            i = k - 1;
//                            while ((i >= 0) && (value < SplCls[i])) {
//                                SplCls[i + 1] = SplCls[i];
//                                --i;
//                            }
//                            SplCls[i + 1] = value;
//                        }
//                        break;
//                    default:
//                        quickSort(SplCls, SplInd);
//                        break;
//                }
                
                for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
                    ind0 = SplCls[sc];
                    ind1 = ind0 + Part->cls[ind0];
                    SplCntInd = 0;
                    if (ElmHitCll[ind0] < Part->cls[ind0]) {
                        SplCnt[SplCntInd++] = 0;
                        SplPos[0] = Part->cls[ind0] - ElmHitCll[ind0];
                    }
                    
                    /* According to the numbers of neighbors of C into the current cell */
                    /* compute how many vertices in C will be placed into the same new cell */
                    iend = ind0 + ElmHitCll[ind0];
                    for (i = ind0; i < iend; i++) {
                        value = NghCounts[HitVtx[i]];
                        if (Markers[value] != tv->mark) {
                            Markers[value] = tv->mark;
                            SplCnt[SplCntInd++] = value;
                            SplPos[value] = 1;
                        }
                        else {
                            SplPos[value]++;
                        }
                    }
                    tv->mark++;
                    
                    /* Sort the values deriving from the previous step */
                    switch (SplCntInd) {
                        case 0:
                        case 1:
                            break;
                        case 2:
                            if (SplCnt[0] > SplCnt[1]) {
                                value = SplCnt[0];
                                SplCnt[0] = SplCnt[1];
                                SplCnt[1] = value;
                            }
                            break;
                        case 3:
                        case 4:
                        case 5:
                        case 6:
                        case 7:
                        case 8:
                            for (k = 1; k < SplCntInd; ++k) {
                                value = SplCnt[k];
                                i = k - 1;
                                while ((i >= 0) && (value < SplCnt[i])) {
                                    SplCnt[i + 1] = SplCnt[i];
                                    --i;
                                }
                                SplCnt[i + 1] = value;
                            }
                            break;
                        default:
                            quickSort(SplCnt, SplCntInd);
                            break;
                    }
                    
                    Part->cells += SplCntInd-1;
                    
                    /* Split the cell C and update the information for sizes of new cells */
                    /* Put the new cells into the stack */
                    i = ind0;
                    for (k = 0; k < SplCntInd; k++) {
                        value = SplPos[SplCnt[k]];
                        Part->cls[i] = value;
                        SplPos[SplCnt[k]] = i;
                        i += value;
                        if (i < ind1) {
                            TheTrace[TraceInd++] = i;
//                            printf("a %d\n", TheTrace[TraceInd-1]);
                        }
                    }
                    
                    /* Permute elements of the cell C */
                    iend = ind0 + ElmHitCll[ind0];
                    
                    for (i = ind0; i < iend; i++) {
                        value = HitVtx[i];
//                        sge = TheGraph[value].e;
//                        if (NghCounts[value]>1) {
//                            printf("%d!\n", NghCounts[value]);
//                            writegroupsize(outfile, tv->stats->grpsize1, tv->stats->grpsize2);
//                            printf("\n");
//                            printf("b %d: ", value+labelorg);
//                            PrintVect(sge, 0, TheGraph[value].d, labelorg);
//                            factorial(&(tv->stats->grpsize1), &(tv->stats->grpsize2), NghCounts[value]);
//                        }
                        Edge_Delete(value, NghCounts[value], Cand, tv);
                        sge = TheGraph[value].e+TheGraph[value].d;
                        if (NghCounts[value] > 1) {
                            factorial(&(tv->stats->grpsize1), &(tv->stats->grpsize2), NghCounts[value]);
//                            printf("B %d: ", value+labelorg);
//                            PrintVect(sge, 0, NghCounts[value], labelorg);
                            
                            if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                            for (j0=0; j0<NghCounts[value]-1; j0++) {
                                SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[j0+1])
//                                printf("[%d %d]a ", sge[j0], sge[j0+1]);
                            }
                            SETPAIRSAUTANDTREE_PREPROC(sge[j0], sge[0])
//                            printf("[%d %d]b \n", sge[j0], sge[0]);
                            SPECIALGENERATORS
                            if (NghCounts[value] > 2) {
                                if (tv->permInd) ResetAutom(tv->permInd, n, tv);
                                SETPAIRSAUTANDTREE_PREPROC(sge[0], sge[1])
//                                printf("[%d %d]c ", sge[0], sge[1]);
                                if (tv->build_autom) {
                                    SETPAIRSAUT(sge[1], sge[0])
//                                    printf("[%d %d]d \n", sge[1], sge[0]);
                                }
                                MakeTree(sge[1], sge[0], sg, n, tv, FALSE);
                                SPECIALGENERATORS
                            }
                        }
                                                
                        j = SplPos[NghCounts[value]]++;         /* where HitVtx[i] goes */
                        k = Cand->invlab[value];				/* where HitVtx[i] is in lab */
                        Cand->lab[k] = Cand->lab[j];
                        Cand->lab[j] = value;
                        Cand->invlab[value] = j;
                        Cand->invlab[Cand->lab[k]] = k;
                        NghCounts[value] = 0;
                    }
                    
                    /* Reconstruct the cell C and update the inverse partition */
                    newcell = ind1 - ElmHitCll[ind0];
                    i = newcell;
                    ind2 = newcell+Part->cls[newcell]-1;
                    do {
                        Part->inv[i] = newcell;
                        if (i == ind2) {
                            newcell = i+1;
                            if (newcell < n) ind2 = newcell+Part->cls[newcell]-1;
                        }
                    }
                    while (++i < ind1);
                    
                    for (i = ind0, k = 0; k < SplCntInd; i+=Part->cls[i], k++) {
                        //if (Part->cls[i] < 100) PrintVect(Cand->lab, i, i+Part->cls[i], labelorg);
                        if ((k > 0) || (SplCnt[0] > 0)) {
                            if (TheGraph[Cand->lab[i]].d == 1) {
                                CStack[CStackInd++] = i;
//                                printf("B %d\n", CStack[CStackInd-1]);
                            }
                        }
                        if (Part->cls[i] == 1) {
                            Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                            //Singletons[SingInd] = i;   ??????
                            SingInd++;
                        }
                    }
                    
                }
            }
        }
        // VA RIVISTA!!!
//        for (i=0; i<n; i++) {
//            if (TheGraph[i].d > 2) {
//                return 1;
//            }
//        }
//        return 2;
        return 1;
    }
    else {
        return 0;
    }
}

void PrintVect(int *v, int z, int n, int l) {
	int i;
	printf("[");
	for (i = z; i<n; i++)
		printf(" %2d", v[i]+l);
	printf(" ]\n");
	return;
}

void MakeTree(int v1, int v2, sparsegraph *sg, int n, struct TracesVars* tv, boolean forceautom) {
    int ind, vtx1, vtx2, ngh1, ngh2, trind, deg0, deg1;
    size_t j1; //, v_vtx1, v_vtx2;
    int *sge1, *sge2;
    boolean build_autom;
    
//    printf("MT (%d %d)\n", v1+labelorg, v2+labelorg);
    if (v1 == v2) return;
    build_autom = tv->build_autom || forceautom;
    trind = 2;
    ind = 0;
    TreeStack[0] = v1;
    TreeStack[1] = v2;
    SETMARK(TreeMarkers, tv->treemark);
    
    while (ind < trind) {
        vtx1 = TreeStack[ind++];
        vtx2 = TreeStack[ind++];
//        printf("vtx1 %d, vtx2 %d (dg: %d -> %d, %d -> %d)\n", vtx1+labelorg, vtx2+labelorg, sg->d[vtx1], sg1->d[vtx1], sg->d[vtx2], sg1->d[vtx2]);
        TreeMarkers[vtx1] = tv->treemark;
        TreeMarkers[vtx2] = tv->treemark;
//        v_vtx1 = sg->v[vtx1];
//        v_vtx2 = sg->v[vtx2];

        
//        deg0 = max(sg1->d[vtx1], 0);
//        deg1 = sg->d[vtx1];
//        sge1 = sg1->e+v_vtx1;
//        sge2 = sg1->e+v_vtx2;

        deg0 = max(TheGraph[vtx1].d, 0);
        deg1 = sg->d[vtx1];
        sge1 = TheGraph[vtx1].e;
        sge2 = TheGraph[vtx2].e;

        //        printf("ngh %d: ", vtx1+labelorg);
        //        PrintVect(sge1, deg0, deg1, labelorg);
        //        printf("ngh %d: ", vtx2+labelorg);
        //        PrintVect(sge2, deg0, deg1, labelorg);
        //        switch (which) {
        //            case 1:
        //                TreeMarkers[sg1->e[vtx1]] = tv->treemark;
        //                TreeMarkers[sg1->e[vtx2]] = tv->treemark;
        //                break;
        //
        //            default:
        //                break;
        //        }
        //        printf("-----\n");
        //        PrintVect(sg1->e+v_vtx1, 0, deg1, labelorg);
        //        PrintVect(sg1->e+v_vtx2, 0, deg1, labelorg);
        //        PrintVect(sg1->e+v_vtx1, deg0, deg1, labelorg);
        //        PrintVect(sg1->e+v_vtx2, deg0, deg1, labelorg);
        //        printf("=====\n");
//        printf("deg0: %d; deg1: %d\n", deg0, deg1);
        for (j1 = deg0; j1 < deg1; j1++) {
            ngh1 = sge1[j1];
            ngh2 = sge2[j1];
//            printf("mt?: (%d %d); TreeMarkers[ngh1]: %d, TreeMarkers[ngh2]: %d, tv->treemark: %d\n", ngh1+labelorg, ngh2+labelorg, TreeMarkers[ngh1], TreeMarkers[ngh2], tv->treemark);
            if ((TreeMarkers[ngh1] != tv->treemark) && (ngh1 != ngh2)) {
                TreeStack[trind++] = ngh1;
                TreeStack[trind++] = ngh2;
//                printf("MOVE %d and %d INTO THE STACK\n", ngh1+labelorg, ngh2+labelorg);
                if (ngh1 != ngh2) {
                    //printf("build_autom: %d\n", build_autom);
                    if (build_autom) {
                        AUTPERM[ngh1] = ngh2;
                        PrmPairs[tv->permInd].arg = ngh1;
                        PrmPairs[tv->permInd].val = ngh2;
                        //                        printf("pp(%d): (%d %d)\n", tv->permInd, ngh1+labelorg, ngh2+labelorg);
                        tv->permInd++;
                    }
                    orbjoin_sp_pair(tv->orbits, OrbList, n, 
                                    ngh1, ngh2, &tv->stats->numorbits);
                }
                //printf("mt: (%d %d)\n", ngh1+labelorg, ngh2+labelorg);
            }
        }
        //        printf("\n");
    }
//    printf("\n");
    return;
}

void MakeCanTree(int v1, sparsegraph *sg_orig, int n, Candidate *Cand, Partition *Part, struct TracesVars* tv) {
    int ind, vtx, ngh, trind, deg0, deg1; //, vtxto, vtxpos;
    size_t j1; //, v_vtx;
    int *sge1; //, *sge2;
    
//    VerifyCand(Cand, n, 11632);
//    if (TreeNodes[v1]) {
//        return;
//    }
//    TreeNodes[v1] = TRUE;
    trind = 1;
    ind = 0;
    TreeStack[0] = v1;
    SETMARK(TreeMarkers, tv->treemark);
//    printf("MakeCanTree @ %d\n", v1+labelorg);

    while (ind < trind) {
        vtx = TreeStack[ind++];
        if (TreeNodes[vtx]) {
            return;
        }

        //if (sg1->d[vtx] == -1) {
        if (TheGraph[vtx].d == -1) {
//            printf("Place %d\n", vtx+labelorg);
            Place(vtx, Cand, Part);
            TreeNodes[vtx] = TRUE;
        }

        //        vtxpos = Cand->invlab[vtx];
        //        printf("vtxpos(%d) = %d; ", vtx+labelorg, vtxpos);
        //        vtxto = CanonIndices[Part->inv[vtxpos]]++;
        //        printf("goesto: %d (lab: %d)\n", vtxto, Cand->lab[vtxto]+labelorg);
        //        //printf("vtx %d (cell: %d, dg: %d -> %d)   %d <-> %d\n", vtx+labelorg, vtxto, sg_orig->d[vtx], sg1->d[vtx], Cand->lab[vtxpos]+labelorg, Cand->lab[vtxto]+labelorg);
        //        if (Cand->lab[vtxpos] != Cand->lab[vtxto]) {
        //            Cand->lab[vtxpos] = Cand->lab[vtxto];
        //            Cand->lab[vtxto] = vtx;
        //            Cand->invlab[Cand->lab[vtxpos]] = vtxpos;
        //            Cand->invlab[Cand->lab[vtxto]] = vtxto;
        //        }
        //        if (Part->cls[vtxto] > 1) {
        //            Part->cls[vtxto+1] = Part->cls[vtxto]-1;
        //            Part->cls[vtxto] = 1;
        //        }
        
        //        PrintVect(Cand->lab, 0, n, labelorg);
        
        TreeMarkers[vtx] = tv->treemark;
//        v_vtx = sg_orig->v[vtx];
//        deg0 = max(sg1->d[vtx], 0);
        deg0 = max(TheGraph[vtx].d, 0);
        deg1 = sg_orig->d[vtx];
//        sge1 = sg1->e+v_vtx;
        sge1 = TheGraph[vtx].e;
        
        for (j1 = deg0; j1 < deg1; j1++) {
            ngh = sge1[j1];
//            if ((sg1->d[ngh] == -1) && (TreeMarkers[ngh] != tv->treemark)) {
            if ((TheGraph[ngh].d == -1) && (TreeMarkers[ngh] != tv->treemark)) {
                TreeStack[trind++] = ngh;
//                printf("ngh %d (dg: %d -> %d)\n", ngh+labelorg, sg_orig->d[ngh], sg1->d[ngh]);
            }
        }
    }
//    printf("return MakeCanTree\n");
//    VerifyCand(Cand, n, 11632);
    return;
}

int max(int u, int v) {
    if (u > v) {
        return u;
    }
    else {
        return v;
    }
}

int min(int u, int v) {
    if (u < v) {
        return u;
    }
    else {
        return v;
    }
}

void putgraphplus_sg(FILE *f, sparsegraph *sg, sparsegraph *sg1, int linelength) {
    int i, n, curlen, slen;
    int *d, *e;
    size_t *v, j;
    char s[60];
    
    n = sg->nv;
    SG_VDE(sg, v, d, e);
//    e = sg1->e;
    for (i = 0; i < n; ++i)
    {
        fprintf(f, "%3d : ", i+labelorg);
        curlen = 7;
        
        //for (j = v[i]; j < v[i]+d[i]; ++j)
        for (j = 0; j < d[i]; ++j)
        {
//            slen = itos(e[j]+labelorg, s);
            slen = itos(TheGraph[v[i]].e[j] + labelorg, s);
            if (linelength > 0 && curlen + slen + 1 > linelength)
            {
                putstring(f, "\n  ");
                curlen = 2;
            }
            PUTC(' ', f);
            putstring(f, s);
            curlen += slen + 1;
//            if (sg1->d[i] == j-v[i]) {
            if (TheGraph[i].d == j-v[i]) {
                printf(" *");
            }
        }
        putstring(f, ";\n");
    }
//    printf(" SG->d: ");
//    PrintVect(sg->d, 0, n, 0);
//    printf("SG1->d: ");
//    PrintVect(sg1->d, 0, n, 0);
    
}

void orbjoin_sp_perm(int *orbits, int *map, int *list, int n, int *numorbs) {
    int i, j1, j2, k1, k2;
//    printf("OJ_PERM\n");
// controllare se numorbs contiene esattamente ilnumero corrente di orbite
    for (i = 0; i < n; ++i)
        if (map[i] != i)
        {
            //            printf("i: %d, map: %d; -- > ", i+labelorg, map[i]+labelorg);
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[map[i]];
            while (orbits[j2] != j2) j2 = orbits[j2];
            //
//            printf("j1: %d, j2: %d\n", j1+labelorg, j2+labelorg);
//            putorbits(outfile, orbits, 0, n);
//            printf("TEMP2: ");
//            PrintVect(TempOrbits, 0, n, labelorg);
//            printf("LIST2: ");
//            PrintVect(TempOrbList, 0, n, labelorg);
            
            k1 = j1;
            k2 = j2;
            if (k1 < k2) {
                (*numorbs)--;
                while (OrbList[j2] != k2) {
                    //printf("%d %d, ", j2+labelorg, OrbList[j2]+labelorg);
                    orbits[j2] = k1;
                    j2 = OrbList[j2];
                }
                //printf("\n");
                orbits[j2] = k1;
                k1 = OrbList[k1];
                //k2 = OrbList[k2];
                //printf("OrbList[%d] = %d\n", j2+labelorg, k1+labelorg);
                OrbList[j2] = k1;
                //printf("OrbList[%d] = %d\n", j1+labelorg, k2+labelorg);
                OrbList[j1] = k2;
            }
            else if (k1 > k2) {
                (*numorbs)--;
                while (OrbList[j1] != k1) {
                    //printf("%d %d, ", j1+labelorg, OrbList[j1]+labelorg);
                    orbits[j1] = k2;
                    j1 = OrbList[j1];
                }
                //printf("\n");
                orbits[j1] = k2;
                //k1 = OrbList[k1];
                k2 = OrbList[k2];
                //printf("OrbList[%d] = %d\n", j1+labelorg, k2+labelorg);
                OrbList[j1] = k2;
                //printf("OrbList[%d] = %d\n", j2+labelorg, k1+labelorg);
                OrbList[j2] = k1;
            }
        }
    
    //    j1 = 0;
    //    for (i = 0; i < n; ++i)
    //        if (orbits[i] == i) ++j1;
    //    PrintVect(IDENTITY_PERM, 0, n, labelorg);
    //    PrintVect(orbits, 0, n, labelorg);
    //    PrintVect(OrbList, 0, n, labelorg);
    //    return j1;
}

void orbjoin_sp_pair(int *orbits, int *list, int n, int u, int v, int *numorbs) {
    int j1, j2, k1, k2;

//    printf("orbs: ");
//    PrintVect(orbits, 0, n, labelorg);
//    printf("list: ");
//    PrintVect(list, 0, n, labelorg);
//    printf("arg: %d, val: %d; -- > ", u+labelorg, v+labelorg);
    
    j1 = orbits[u];
    while (orbits[j1] != j1) j1 = orbits[j1];
    j2 = orbits[v];
    while (orbits[j2] != j2) j2 = orbits[j2];
 
//    printf("arg: %d, val: %d\n", j1+labelorg, j2+labelorg);
    
    k1 = j1;
    k2 = j2;
    if (k1 < k2) {
        (*numorbs)--;
        while (list[j2] != k2) {
            orbits[j2] = k1;
            j2 = list[j2];
        }
        orbits[j2] = k1;
        k1 = list[k1];
        list[j2] = k1;
        list[j1] = k2;
    }
    else if (k1 > k2) {
        (*numorbs)--;
        while (list[j1] != k1) {
            orbits[j1] = k2;
            j1 = list[j1];
        }
        orbits[j1] = k2;
        k2 = list[k2];
        list[j1] = k2;
        list[j2] = k1;
    }

//    PrintVect(orbits, 0, n, labelorg);
//    PrintVect(list, 0, n, labelorg);
    //    return j1;
    //    for (j=tv->tcell; j<tv->tcell+Spine[tv->fromlevel].part->cls[tv->tcell]; j++) {
    //        printf("((%d %d)); ", CurrCand->lab[j]+labelorg, orbits[CurrCand->lab[j]]+labelorg);
    //    }
    //    printf("\n");
    //printf("numorbs: %d\n", *numorbs);
    return;
}

boolean isautom_sg_pair(graph *g, int *p, boolean digraph, int m, int n, struct TracesVars *tv) {
    int *d, *e;
    size_t *v;
    int i, k, pi, di;
    size_t vi, vpi, j;
    
    SG_VDE(g, v, d, e);
    
    //    printf("ISAUT by ");
    //    for (k = 0; k < tv->permInd; ++k) {
    //        printf("[%d %d]", PrmPairs[k].arg, PrmPairs[k].val);
    //    }
    //    printf("\n");
    for (k = 0; k < tv->permInd; ++k)
        //if (p[i] != i || digraph)
    {
        i = PrmPairs[k].arg;
        pi = p[i];
        di = d[i];
        if (d[pi] != di) return FALSE;
        
        vi = v[i];
        vpi = v[pi];
        //            printf("[%d %d]\n", i+labelorg, pi+labelorg);
        SETMARK(AutMarkers, tv->autmark)
        for (j = 0; j < di; ++j) AutMarkers[p[e[vi+j]]] = tv->autmark;
        for (j = 0; j < di; ++j) if (AutMarkers[e[vpi+j]] != tv->autmark) {
            //printf("wrong @ %d (from %d -> %d)\n", e[vpi+j]+labelorg, i+labelorg, pi+labelorg);
            return FALSE;
        }
    }
    
    return TRUE;
}

void SetAutom(int q, int n, struct TracesVars *tv) {
    int i;
    //    if (q) {
    //        ResetAutom(q, n, tv);
    //    }
    for (i=0; i<q; i++) {
        AUTPERM[PrmPairs[i].arg] = PrmPairs[i].val;
        //printf("%d -> %d\n", PrmPairs[i].arg+labelorg, PrmPairs[i].val+labelorg);
    }
    return;
}

void ResetAutom(int q, int n, struct TracesVars *tv) {
    int i;
    
    if (n/q < 256) {
        memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
    }
    else {
        for (i=0; i<q; i++) {
            AUTPERM[PrmPairs[i].arg] = PrmPairs[i].arg;
        }
    }
    //VerifyId(AUTPERM, n);
    tv->permInd = 0;
    return;
}

boolean VerifyId(int *p, int n) {
    int i, r;
    r = TRUE;
    for (i=0; i<n; i++) {
        if (p[i] != i) {
            printf("p[%d] = %d\n", i, p[i]);
            r = FALSE;
        }
    }
    return r;
}

void PrintPartition(int *v, int *cls, int n, int l, int line) {
	int i, j;
	fprintf(outfile, "[ ");
	for (i=0; i<n; i+=cls[i]) {
        if ((cls[i]<=0) || i>=n) {
            printf("WRONG");
            break;
        }
		for (j=i; j<i+cls[i]; j++) {
			fprintf(outfile, "%d ", v[j]+l);
		}
		if ((i+cls[i])<n) fprintf(outfile, "| ");
	}
	fprintf(outfile, "] at line %d\n", line);
	return;
}

void Place(int vtx, Candidate *Cand, Partition *Part) {
    int vtxto, vtxpos;
    
    vtxpos = Cand->invlab[vtx];
//    printf("vtxpos(%d) = %d; ", vtx+labelorg, vtxpos);
    vtxto = CanonIndices[Part->inv[vtxpos]]++;
//    printf("goesto: %d (lab: %d)\n", vtxto, Cand->lab[vtxto]+labelorg);
    //printf("vtx %d (cell: %d, dg: %d -> %d)   %d <-> %d\n", vtx+labelorg, vtxto, sg->d[vtx], sg1->d[vtx], Cand->lab[vtxpos]+labelorg, Cand->lab[vtxto]+labelorg);
    if (Cand->lab[vtxpos] != Cand->lab[vtxto]) {
        Cand->lab[vtxpos] = Cand->lab[vtxto];
        Cand->lab[vtxto] = vtx;
        Cand->invlab[Cand->lab[vtxpos]] = vtxpos;
        Cand->invlab[Cand->lab[vtxto]] = vtxto;
    }
    if (Part->cls[vtxto] > 1) {
        //printf("Part->cls[%d]: %d ", vtxto, Part->cls[vtxto]);
        Part->cls[vtxto+1] = Part->cls[vtxto]-1;
        Part->cls[vtxto] = 1;
        //printf("--> %d %d\n", Part->cls[vtxto], Part->cls[vtxto+1]);
    }
}

int NonSingDeg(int vtx, Candidate *Cand, Partition *Part) {
    int *e_vtx;
    int i, deg, retdeg;
    
    retdeg = TheGraph[vtx].d;
    deg = retdeg;
    e_vtx = TheGraph[vtx].e;
    for (i=0; i<deg; i++) {
        if (Part->cls[Part->inv[Cand->invlab[e_vtx[i]]]] == 1) {
            retdeg--;
        }
    }
    return retdeg;
}

int NonSingDegPlus1(Candidate *Cand, Partition *Part, int cell, TracesVars *tv) {
    
    int *e_vtx;
    int vtx, sing;
    int i, j, deg, retdeg, n, singcount;
    
    n = tv->input_graph->nv;
    singcount = 0;
    
    SETMARK(StackMarkers, tv->stackmark)
    
    //    printf("Cell: ");
    //    PrintVect(Cand->lab, cell, cell+Part->cls[cell], labelorg);
    for (j=cell; j<cell+Part->cls[cell]; j++) {
        vtx = Cand->lab[j];
        deg = TheGraph[vtx].d;
        retdeg = 0;
        e_vtx = TheGraph[vtx].e;
        
        for (i=0; i<deg; i++) {
            //            printf("ngh: %d (%d)\n", e_vtx[i]+labelorg, Part->cls[Part->inv[Cand->invlab[e_vtx[i]]]]);
            //            if (Part->cls[Part->inv[Cand->invlab[e_vtx[i]]]] != 1) {
            if (Singletons[e_vtx[i]] != 1) {
                e_vtx[retdeg++] = e_vtx[i];
            }
            else {
                if (StackMarkers[e_vtx[i]] != tv->stackmark) {
                    sing = e_vtx[i];
                    WorkArray2[singcount] = Part->inv[Cand->invlab[sing]];
                    WorkArray[singcount++] = sing;
                    //                    printf("%d in cell %d (#%d)\n", sing, Part->inv[Cand->invlab[sing]], singcount);
                    StackMarkers[e_vtx[i]] = tv->stackmark;
                    //                    singdeg = 0;
                    //                    deg1 = TheGraph[sing].d;
                    //                    e_sing = TheGraph[sing].e;
                    //                    //                    printf("%d: %d; Neigh: ", sing+labelorg, deg1);
                    //                    //                    PrintVect(e_sing, 0, deg1, labelorg);
                    //                    for (k=0; k<deg1; k++) {
                    //                        //                        printf("*");
                    //                        if (Part->inv[Cand->invlab[e_sing[k]]] != cell) {
                    //                            //                            printf("%d;", e_sing[k]+labelorg);
                    //                            e_sing[singdeg++] = e_sing[k];
                    //                        }
                    //                        //                        else printf("[%d];", e_sing[k]+labelorg);
                    //                    }
                    //                    //                    printf("\n");
                }
                //                TheGraph[sing].d = singdeg;
                //                printf("%d: %d; Neigh: ", sing+labelorg, deg1);
                //                PrintVect(e_sing, 0, singdeg, labelorg);
            }
        }
        if (j == cell) {
            sort2ints(WorkArray2, WorkArray, singcount);
        }
        
//        printf("Nghb of %d Sings and their cells\n", vtx+labelorg);
//        PrintVect(WorkArray, 0, singcount, labelorg);
//        PrintVect(WorkArray2, 0, singcount, 0);
        
        //memcpy(e_vtx, WorkArray, retdeg*sizeof(int));
        if (deg != retdeg) {
            //            printf("Copy [0..%d) ", singcount);
            //            PrintVect(WorkArray, 0, singcount, labelorg);
            memcpy(e_vtx+retdeg, WorkArray, singcount*sizeof(int));
            //            printf("reduce nghbhood of %d to %d\n", vtx+labelorg, retdeg);
            //            PrintVect(e_vtx, 0, TheGraph[vtx].d, labelorg);
            //            TheGraph[vtx].d = sg->d[vtx] = retdeg;  // XXXYYY
            TheGraph[vtx].d = retdeg;
        }
    }
    return retdeg;
}

void NonSingDegPlus2(Candidate *Cand, Partition *Part, int cell, TracesVars *tv) {
    
    int *e_sing;
    int sing;
    int k, deg1, singdeg, singcount;
    
//    n = tv->input_graph->nv;
    singcount = 0;
    
    
    //    printf("Cell: ");
    //    PrintVect(Cand->lab, cell, cell+Part->cls[cell], labelorg);
    sing = Cand->lab[cell];
    singdeg = 0;
    deg1 = TheGraph[sing].d;
    e_sing = TheGraph[sing].e;
    //                    printf("%d: %d; Neigh: ", sing+labelorg, deg1);
    //                    PrintVect(e_sing, 0, deg1, labelorg);
    for (k=0; k<deg1; k++) {
        //                        printf("*");
        if (Singletons[e_sing[k]] != 2) {
            e_sing[singdeg++] = e_sing[k];
        }
    }
    TheGraph[sing].d = singdeg;
}

void Edge_Delete(int vertex, int sons, Candidate *Cand, TracesVars *tv) {
    int d_vtx, j1, temp;
    int *sge;

    if (TheGraph[vertex].d <= 1) {
//        printf("edge delete return\n");
        return;
    }

    d_vtx = TheGraph[vertex].d = TheGraph[vertex].d - sons;
//    printf("edge delete from %d; newdeg: %d (%d sons)\n", vertex+labelorg, d_vtx, sons);
    sge = TheGraph[vertex].e;
//    PrintVect(sge, 0, sg->d[vertex], labelorg);

    for (j1=0; j1<d_vtx; j1++) {
        if (TheGraph[sge[j1]].one) {
            while (TheGraph[sge[TheGraph[vertex].d]].d == -1) {
                (TheGraph[vertex].d)++;
            }
            temp = sge[j1];
            sge[j1] = sge[TheGraph[vertex].d];
            sge[TheGraph[vertex].d] = temp;
        }
    }
//    if (sons>1) {
//        printf("%d: ", vertex+labelorg);
//        PrintVect(TheGraph[vertex].e+d_vtx, 0, sons, labelorg);
//    }
//    sg1->d[vertex] = TheGraph[vertex].d = d_vtx;
    TheGraph[vertex].d = d_vtx;  // XXXYYY
//    printf("--> ");
//    PrintVect(sge, 0, sg->d[vertex], labelorg);
}

boolean VerifyPart(int *Part, int start, int end, int *lab) {
    int i;
    for (i=start; i<end; i+=Part[i]) {
        if (Part[i] == 0 || i>=end) {
            return FALSE;
        }
    }
    return TRUE;
}

int VerifyPerm(int *perm, int n,int where) {
    int i;
    memset(Markers, 0, n*sizeof(int));
    
    for (i=0; i<n; i++) {
        if ((perm[i] >= n) || (Markers[perm[i]])) {
            fprintf(stderr,"wrong permutation @ %d\n",where);
            PrintVect(perm,0,i+1,labelorg);
        }
        Markers[perm[i]] = TRUE;
    }
    return TRUE;
}

boolean VerifyCand(Candidate *Cand, int n, int line) {
    int i, k;
    for (i=0; i<n; i++) {
        k=Cand->lab[i];
        if (Cand->invlab[k] != i) {
            printf("Cand->invlab wrong at %d (vtx: %d), line %d\n", i, k, line);
            PrintVect(Cand->lab, 0, n, 0);
            PrintVect(Cand->invlab, 0, n, 0);
            return FALSE;
        }
    }
    return TRUE;
}

//boolean Same_Cell(int *cell1, int *cell2, int k, TracesVars *tv, int n) {
//    int i;
//    
//    SETMARK(CellMarkers, tv->autmark)
//    for (i=0; i<k; i++) CellMarkers[cell1[i]] = tv->autmark;
//    for (i=0; i<k; i++) if (CellMarkers[cell2[i]] != tv->autmark) return FALSE;
//    return TRUE;
//}

int FirstNeighbour(int vtx, Candidate *Cand, Partition *Part, int* Markers, int mark, int *ngh, int n) {
//    int *d, *e, *e_vtx;
    int *e_vtx;
//    size_t *v;
    int i, k, deg;
    int ngh1, ngh2, cell1, cell2;
    
    k = 0;
//    SG_VDE(sg, v, d, e);  //????
    
    deg = TheGraph[vtx].d;
    e_vtx = TheGraph[vtx].e;
    
//    if (deg == sg->nv-1) {
    if (deg == n-1) {
        return 0;
    }
//    printf("Neighbors of %d: ", vtx+labelorg);
//    PrintVect(e_vtx, 0, deg, labelorg);
    
    for (i=0; i<deg; i++) {
        //printf("Markers[%d]: %d (mark: %d)\n", e_vtx[i], Markers[e_vtx[i]], mark);
        if (Markers[e_vtx[i]] != mark) {
            cell1 = Part->inv[Cand->invlab[e_vtx[i]]];
            if (Part->cls[cell1] > 1) {
                ngh1 = e_vtx[i++];
                k++;
                break;
            }
        }
    }
    for (; i<deg; i++) {
        //printf("MARKERS[%d]: %d (mark: %d)\n", e_vtx[i], Markers[e_vtx[i]], mark);
        if (Markers[e_vtx[i]] != mark) {
            cell2 = Part->inv[Cand->invlab[e_vtx[i]]];
            if (Part->cls[cell2] > 1) {
                ngh2 = e_vtx[i];
                k++;
                break;
            }
        }
    }
    //    printf("K is %d\n", k);
    switch (k) {
        case 0:
            break;
            
        case 1:
            *ngh = ngh1;
            break;
            
        case 2:
            //            printf("(%d in cell1: %d, %d in cell2: %d)\n", ngh1+labelorg, cell1, ngh2+labelorg, cell2);
            if (cell1 < cell2) {
                *ngh = ngh1;
            }
            else {
                *ngh = ngh2;
            }
            break;
            
        default:
            break;
    }
    return k;
}

int NextNeighbour(int vtx, Candidate *Cand, Partition *Part, int* Markers, int mark, int *ngh, int n) {
//    int *d, *e, *e_vtx;
    int *e_vtx;
//    size_t *v;
    int i, deg; //, col;
//    int ngh1, ngh2, cell1, cell2;
    int cell1;
    
//    SG_VDE(sg, v, d, e);  //????
    
    deg = TheGraph[vtx].d;
    e_vtx = TheGraph[vtx].e;
    
//    if (deg == sg->nv-1) {
    if (deg == n-1) {
        return 0;
    }

//    printf("Vicini di %d: ", vtx+labelorg);
//    PrintVect(e_vtx, 0, deg, labelorg);
    for (i=0; i<deg; i++) {
        //printf("MARKERS[%d]: %d (mark: %d)\n", e_vtx[i]+labelorg, Markers[e_vtx[i]], mark);
        if (Markers[e_vtx[i]] != mark) {
            cell1 = Part->inv[Cand->invlab[e_vtx[i]]];
            if (Part->cls[cell1] > 1) {
                *ngh = e_vtx[i];
                break;
            }
        }
    }
    if (i<deg) return 1; else return 0;
}

sparsegraph* copy_sg_structure(sparsegraph *sg1, sparsegraph *sg2)
{
    int *d1, *e1, *d2, *e2;
    int i, n;
    size_t *v1, *v2, k;
    
    if (!sg2)
    {
        if ((sg2 = (sparsegraph*)ALLOCS(1, sizeof(sparsegraph))) == NULL)
        {
            fprintf(ERRFILE, "copy_sg: malloc failed\n");
            exit(1);
        }
        SG_INIT(*sg2);
    }
    
    SG_VDE(sg1, v1, d1, e1);
    
    n = sg1->nv;
    
    k = 0;
    for (i = 0; i < n; ++i)
        if (v1[i]+d1[i]>k) k = v1[i] + d1[i];
    
    SG_ALLOC(*sg2, n, k, "copy_sg malloc");
    SG_VDE(sg2, v2, d2, e2);
    
    sg2->nv = n;
    sg2->nde = sg1->nde;
//    memcpy(v2, v1, n*sizeof(size_t));
//    memcpy(d2, d1, n*sizeof(int));
//    memcpy(e2, e1, k*sizeof(int));
    
    return sg2;
}
