/* genspecialg.c  version 1.0; B D McKay, July 6, 2015 */

#define USAGE "genspecialg \n\
[-s|-g|-z] [-p#|-c#|-e#|-k#|-b#,#|-q#|-P#,#|C#,#...|G#,#...] [outfile]"

#define HELPTEXT \
" Generate one particular graph.\n\
     #  : size parameter called n in the descriptions\n\
\n\
    -s : Write in sparse6 format (default)\n\
    -g : Write in graph6 format\n\
    -z : Make digraph versions and write in digraph6 format\n\
\n\
    If defined, the digraph version is shown in parentheses:\n\
    -p#   : path (directed path) on n vertices.\n\
    -c#   : cycle (directed cycle) on n vertices.\n\
    -e#   : empty graph (digraph with loops only) on n vertices.\n\
    -k#   : complete graph (with loops) on n vertices\n\
    -b#,# : complete bipartite graph (directed l->r) on n vertices\n\
    -P#,# : generalized Petersen graph; usual one is -P5,2\n\
    -q#   : hypercube on 2^n vertices and degree n.\n\
    -C#,#... : circulant graph (digraph).\n\
    -G#,#... : (directed) grid, use negative values for open directions\n"

/* Ideas: multipartite, knesser, full trees */

#include "gtools.h"

#define MAXARGS 1000  /* Maximum argument list for multi-argument parameters */

static long args[MAXARGS];

/**************************************************************************/

static int
vnumber(long *dimen, int *index, int ndimen)
{
    int i,v;
    
    v = 0;
    for (i = 0; i < ndimen; ++i)
        v = v*dimen[i] + index[i];

    return v;
}

/**************************************************************************/

static void
makepath(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i;
    size_t *v,k;

    if (n < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -p\n");

    if (digraph) SG_ALLOC(*sg,n,n-1,"genspecialg");
    else         SG_ALLOC(*sg,n,2UL*n-2,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph || n == 1)
    {
	sg->nv = n;
	sg->nde = n-1;

	for (i = 0; i < n-1; ++i)
        {
	    d[i] = 1;
	    v[i] = i;
	    e[i] = i+1;
        }
	d[n-1] = 0;
	v[n-1] = 0;
    }
    else
    {
	sg->nv = n;
	sg->nde = 2*n-2;

	d[0] = 1;
	v[0] = 0;
	e[0] = 1;
	for (i = 1, k = 1; i < n-1; ++i, k += 2)
	{
	    d[i] = 2;
	    v[i] = k;
	    e[k++] = i-1;
	    e[k++] = i+2;
	}
	d[n-1] = 1;
	v[n-1] = k;
	e[k] = n-2;
    }
}

/**************************************************************************/

static void
makecycle(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i;
    size_t *v,k;

    if (!digraph && (n < 1 || n == 2 || n > NAUTY_INFINITY-2))
        gt_abort(">E genspecialg: bad argument for -c\n");
    if (digraph && (n < 1 || n > NAUTY_INFINITY-2))
        gt_abort(">E genspecialg: bad argument for -zc\n");

    if (digraph) SG_ALLOC(*sg,n,n,"genspecialg");
    else         SG_ALLOC(*sg,n,2UL*n,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph || n == 1)
    {
	sg->nv = n;
	sg->nde = n;

	for (i = 0; i < n-1; ++i)
        {
	    d[i] = 1;
	    v[i] = i;
	    e[i] = i+1;
        }
	d[n-1] = 1;
	v[n-1] = n-1;
	e[n-1] = 0;
    }
    else
    {
	sg->nv = n;
	sg->nde = 2UL*n;

	d[0] = 2;
	v[0] = 0;
	e[0] = 1;
	e[1] = n-1;

	for (i = 1; i < n-1; ++i)
        {
	    d[i] = 2;
	    v[i] = 2UL*i;
	    e[2UL*i] = i-1;
	    e[2UL*i+1] = i+1;
        }
	d[n-1] = 2;
	v[n-1] = 2UL*n-2;
	e[2UL*n-2] = 0;
	e[2UL*n-1] = n-2;
    }
}

/**************************************************************************/

static void
makecomplete(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j;
    size_t *v,k;

    if (n < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -k\n");

    if (digraph) SG_ALLOC(*sg,n,n*(size_t)n,"genspecialg");
    else         SG_ALLOC(*sg,n,n*(size_t)(n-1),"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph)
    {
        sg->nv = n;
	sg->nde = n*(size_t)n;

	for (i = 0, k = 0; i < n; ++i, k += n)
	{
	    d[i] = n;
	    v[i] = k;
	    for (j = 0; j < n; ++j) e[k+j] = j;
        }
    }
    else
    {
        sg->nv = n;
	sg->nde = n*(size_t)(n-1);

	for (i = 0, k = 0; i < n; ++i)
	{
	    d[i] = n-1;
	    v[i] = k;
	    for (j = 0; j < n; ++j)
               if (j != i) e[k++] = j;
        }
    }
}

/**************************************************************************/

static void
makeempty(long n, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i;
    size_t *v;

    if (n < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -e\n");

    if (digraph) SG_ALLOC(*sg,n,n,"genspecialg");
    else         SG_ALLOC(*sg,n,0,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph)
    {
        sg->nv = n;
	sg->nde = n;

	for (i = 0; i < n; ++i)
	{
	    d[i] = 1;
	    v[i] = i;
	    e[i] = i;
        }
    }
    else
    {
        sg->nv = n;
	sg->nde = 0;

	for (i = 0; i < n; ++i)
	{
	    d[i] = 0;
	    v[i] = 0;
        }
    }
}

/**************************************************************************/

static void
makehypercube(long deg, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j;
    size_t *v,k,nv;

    if (deg < 0 || deg > 30)
        gt_abort(">E genspecialg: bad argument for -q\n");
    if (digraph)
        gt_abort(">E genspecialg: -zq is not implemented\n");

    nv = 1UL << deg;
    SG_ALLOC(*sg,nv,deg*nv,"genspecialg");

    SG_VDE(sg,v,d,e);

    sg->nv = nv;
    sg->nde = deg*nv;

    for (i = 0, k = 0; i < nv; ++i, k += deg)
    {
	d[i] = deg;
	v[i] = k;
	for (j = 0; j < deg; ++j) e[k+j] = i ^ (1<<j);
    }
}

/**************************************************************************/

static void
makegrid(long *dim, int ndim, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,deg,n,oldn;
    size_t *v,k;
    boolean closed[30];
    int index[30];

    n = 1;
    deg = 0;
    for (i = 0; i < ndim; ++i)
    {
        if (dim[i] >= -1 && dim[i] <= 1)
            gt_abort(">E genspecialg: -G dimensions must be at least 2\n");
	if (dim[i] == 2 && !digraph)
            gt_abort(">E genspecialg: -G dimen 2 is only ok for digraphs\n");

	closed[i] = (dim[i] > 0);
	if (dim[i] < 0) dim[i] = -dim[i];

	oldn = n;
        n *= dim[i];
	if (n < 0 || n / dim[i] != oldn)
	    gt_abort(">E genspecialg: -G size is too big\n");

	if (digraph || dim[i] == 2) ++deg;
        else                        deg += 2;

        index[i] = 0;
    }

    if (n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: -G size is too big\n");

    SG_ALLOC(*sg,n,deg*(size_t)n,"genspecialg");

    SG_VDE(sg,v,d,e);

    sg->nv = n;
    sg->nde = deg*(size_t)n;

    k = 0;
    for (i = 0; i < n; ++i)
    {
	v[i] = k;
	for (j = 0; j < ndim; ++j)
	{
	    if (index[j] < dim[j]-1)
	    {
		++index[j];
		e[k++] = vnumber(dim,index,ndim);
		--index[j];
	    }
	    if (!digraph && index[j] > 0)
	    {
		--index[j];
		e[k++] = vnumber(dim,index,ndim);
		++index[j];
	    }
	    if (closed[j] && index[j] == dim[j]-1)
	    {
		index[j] = 0;
		e[k++] = vnumber(dim,index,ndim);
		index[j] = dim[j]-1;
	    }
	    if (closed[j] && !digraph && index[j] == 0)
	    {
		index[j] = dim[j]-1;
		e[k++] = vnumber(dim,index,ndim);
		index[j] = 0;
	    }
	}

        d[i] = k - v[i];

	for (j = ndim; --j >= 0;)
	{
	    if (index[j] != dim[j]-1)
	    {
		++index[j];
		break;
	    }
	    else
		index[j] = 0;
	}
    }
}

/**************************************************************************/

static void
makecirculant(long n, long *conn, int nconn, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,deg;
    size_t *v,k;

    if (nconn > 0 && conn[0] <= 0)
        gt_abort(">E genspecialg: -C connections must be nonzero\n");

    for (i = 1; i < nconn; ++i)
	if (conn[i] <= conn[i-1])
	    gt_abort(">E genspecialg: -C connections must be increasing\n");

    if (nconn == 0)
	deg = 0;
    else
    {
        if (digraph)
	{
	    if (conn[nconn-1] >= n) gt_abort(
                 ">E genspecialg: -C connections must be 1..n-1\n");
	    deg = nconn;
	}
	else
        {
            if (conn[nconn-1] > n/2) gt_abort(
                 ">E genspecialg: -C connections must be 1..n/2\n");
	    deg = 2*nconn - (2*conn[nconn-1]==n);
	}
    }

    SG_ALLOC(*sg,n,deg*n,"genspecialg");

    SG_VDE(sg,v,d,e);
    sg->nv = n;
    sg->nde = deg*n;

    for (i = 0; i < n; ++i)
    {
        d[i] = deg;
	v[i] = deg*(size_t)i;
    }
 
    for (i = 0; i < n; ++i)
    {
	k = v[i];
	for (j = 0; j < nconn; ++j)
	{
	    e[k++] = (i + conn[j]) % n;
	    if (!digraph && 2*conn[j] != n)
		e[k++] = (i - conn[j] + n) % n;
	}
    }
}

/**************************************************************************/

static void
makegenpetersen(long n1, long n2, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,n;
    size_t *v,k;

    if (digraph) gt_abort(">E no digraph version of -P is implemented\n");

    n = 2*n1;
    if (n < 1 || n1 > NAUTY_INFINITY/2-1 || n2 < 1 || 2*n2 >= n1)
	gt_abort(">E -Pm,k needs m>0,0<k<m/2; or m too large\n");

    SG_ALLOC(*sg,n,3UL*n,"genspecialg");

    SG_VDE(sg,v,d,e);
    sg->nv = n;
    sg->nde = 3UL*n;

    for (i = 0; i < n; ++i)
    {
        d[i] = 3;
	v[i] = 3UL*i;
    }

    for (i = 0; i < n1; ++i)
    {
	k = v[i];
	e[k] = (i + 1) % n1;
	e[k+1] = (i + n1 - 1) % n1;
	e[k+2] = i + n1;
    }
    
    for (i = 0; i < n1; ++i)
    {
	k = v[n1+i];
	e[k] = n1 + (i + n2) % n1;
        e[k+1] = n1 + (i - n2 + n1) % n1;
	e[k+2] = i;
    }
} 

/**************************************************************************/

static void
makecompletebipartite(long n1, long n2, boolean digraph, sparsegraph *sg)
{
    int *d,*e,i,j,n;
    size_t *v,k;

    n = n1 + n2;

    if (n1 < 1 || n2 < 1 || n > NAUTY_INFINITY-2)
        gt_abort(">E genspecialg: bad argument for -b\n");

    if (digraph) SG_ALLOC(*sg,n,n1*n2,"genspecialg");
    else         SG_ALLOC(*sg,n,2*n1*n2,"genspecialg");

    SG_VDE(sg,v,d,e);

    if (digraph)
    {
	sg->nv = n;
	sg->nde = n1*n2;

	for (i = 0, k = 0; i < n1; ++i)
	{
	    d[i] = n2;
	    v[i] = k;
	    for (j = n1; j < n; ++j) e[k++] = j;
	}
	for (i = n1; i < n; ++i)
	{
	    d[i] = 0;
	    v[i] = 0;
        }
    }
    else
    {
	sg->nv = n;
	sg->nde = 2*n1*n2;

	for (i = 0, k = 0; i < n1; ++i)
	{
	    d[i] = n2;
	    v[i] = k;
	    for (j = n1; j < n; ++j) e[k++] = j;
	}
	for (i = n1; i < n; ++i)
	{
	    d[i] = n1;
	    v[i] = k;
	    for (j = 0; j < n1; ++j) e[k++] = j;
        }
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    int n,codetype;
    int argnum,i,j;
    char *arg,sw;
    boolean badargs;
    boolean Cswitch,Pswitch,gswitch,sswitch,zswitch;
    boolean pswitch,cswitch,eswitch,kswitch,bswitch,qswitch,Gswitch;
    long size;
    static FILE *outfile;
    char *outfilename;
    sparsegraph sg;
    boolean usesparse,digraph;
    long Pargs[2],bargs[2];
    int nPargs,nbargs,nCargs,nGargs;

    HELP;

    gswitch = sswitch = zswitch = Pswitch = FALSE;
    pswitch = cswitch = eswitch = kswitch = FALSE;
    Gswitch = Cswitch = bswitch = qswitch = FALSE;

    outfilename = NULL;

    argnum = 0;
    badargs = FALSE;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('z',zswitch)
                else SWLONG('p',pswitch,size,"genspecialg -p")
                else SWLONG('c',cswitch,size,"genspecialg -c")
                else SWLONG('e',eswitch,size,"genspecialg -e")
                else SWLONG('k',kswitch,size,"genspecialg -k")
                else SWLONG('q',qswitch,size,"genspecialg -q")
                else SWSEQUENCE('b',",",bswitch,bargs,2,
                                nbargs,"genspecialg -b")
                else SWSEQUENCE('P',",",Pswitch,Pargs,2,
                                nPargs,"genspecialg -P")
                else SWSEQUENCE('C',",",Cswitch,args,MAXARGS,
                                nCargs,"genspecialg -C")
                else SWSEQUENCE('G',",",Gswitch,args,30,
                                nGargs,"genspecialg -G")
                else badargs = TRUE;
            }
        }
        else
        {
            ++argnum;
            if (argnum == 1) outfilename = arg;
            else             badargs = TRUE;
        }
    }

    if ((gswitch!=0) + (sswitch!=0) + (zswitch!=0) > 1)
        gt_abort(">E genspecialg: -gsz are incompatible\n");
 
    if ((pswitch!=0) + (cswitch!=0) + (eswitch!=0) + (kswitch!=0)
           + (bswitch!=0) + (qswitch!=0) + (Pswitch!=0) 
           + (Cswitch!= 0) + (Gswitch!=0) > 1)
        gt_abort(">E genspecialg: -pckbqPCG are incompatible\n");

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (gswitch)      codetype = GRAPH6;
    else if (zswitch) codetype = DIGRAPH6;
    else              codetype = SPARSE6;

    if (!outfilename || outfilename[0] == '-')
    {
        outfilename = "stdout";
        outfile = stdout;
    }
    else if ((outfile = fopen(outfilename,"w")) == NULL)
    {
        fprintf(stderr,"Can't open output file %s\n",outfilename);
        gt_abort(NULL);
    }

    SG_INIT(sg);

    if (pswitch)
        makepath(size,zswitch,&sg);
    else if (cswitch)
        makecycle(size,zswitch,&sg);
    else if (kswitch)
        makecomplete(size,zswitch,&sg);
    else if (eswitch)
        makeempty(size,zswitch,&sg);
    else if (qswitch)
        makehypercube(size,zswitch,&sg);
    else if (bswitch)
    {
	if (nbargs != 2) gt_abort(">E genspecialg: -b needs two arguments\n");
        makecompletebipartite(bargs[0],bargs[1],zswitch,&sg);
    }
    else if (Pswitch)
    {
	if (nPargs != 2) gt_abort(">E genspecialg: -P needs two arguments\n");
        makegenpetersen(Pargs[0],Pargs[1],zswitch,&sg);
    }
    else if (Cswitch)
        makecirculant(args[0],args+1,nCargs-1,zswitch,&sg);
    else if (Gswitch)
    {
	if (nGargs < 2)
            gt_abort(">E genspecialg: -G needs at least two arguments\n");
        makegrid(args,nGargs,zswitch,&sg);
    }

    if (codetype == GRAPH6)        writeg6_sg(outfile,&sg);
    else if (codetype == DIGRAPH6) writed6_sg(outfile,&sg);
    else                           writes6_sg(outfile,&sg);
    
    exit(0);
}
