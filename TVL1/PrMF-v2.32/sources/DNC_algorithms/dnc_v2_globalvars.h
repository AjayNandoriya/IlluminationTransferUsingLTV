/* global variables */

	/* Matlab input and output */
	#define mtxImage_IN		prhs[0]
	#define dbLambda_IN		prhs[1]
	#define mtxLambda_IN	prhs[1]
	#define type_IN			prhs[2]
	#define mtxImage_OUT	plhs[0]

	double lambda;					/* parameter lambda if a scalar lambda is given */
	double *Lambda;					/* parameter Lambda if a matrix Lambda is given */
	double twolambda;				/* = 2*lambda */
	double *neighborType;				/* pointing to the uint32 matrix giving neighbor graph topology */
	edge_t *out_edges;					/* absolute positions of neighbor edges */
	int n_edges;					/* # of edges in neighborhood relationship */
	int edge_size;					/* =4 for 2D images, =5 for 3D images */
	int mm;		/* # of rows of the input image */
	int nn;		/* # of columns of the input image */
	int ll;		/* # of 3rd dimension depth of the input image */

	/* graph initialization part */
	cType vtcap;

	/* max flow part */
	long   _n;                    /* number of nodes */
	long   _m;                    /* number of arcs */
	long &n=_n, &m=_m;	/* copy of n and m */
	node   *nodes;               /* array of nodes */
	arc    *arcs;                /* array of arcs */
	bucket *buckets;             /* array of buckets */
	bucket *stateBuckets;		 /* array of buckets to store nodes of different states */
	node   *source;              /* source node pointer */
	node   *sink;                /* sink node pointer */
	long   dMax;                 /* maximum distance label */
	long   aMax;                 /* maximum active node label */
	long   aMin;                 /* minimum active node label */
	double flow;                 /* flow value */
	#define EPSCONST 0.00001
	double eps;					 /* small epsilon */

	node	*sentinelNode;        /* end of the node list marker */
	arc		*stopA;                /* used in forAllArcs */
	float	globUpdtFreq;          /* global update frequency */

	/* counters part */
	unsigned long pushCnt;          /* number of pushes */
	unsigned long relabelCnt;       /* number of relabels */
	unsigned long updateCnt;        /* number of updates */
	unsigned long gapCnt;           /* number of gaps */
	unsigned long gNodeCnt;         /* number of nodes after gap */

	/* time part */
	clock_t start_time, times[5];
	#ifdef _DEBUG_TIMING
	clock_t updTime, mxfTime, tmpTime;
	#endif

	/* parametric update part */
	int			depth;			/* depth of the binary tree "tree" */
	int			*pTree;			/* pointers to tree node */
	//unsigned int current_lvl,last_lvl;	/* current and last working levels */
	unsigned int lvlCnt;		/* number of distinct levels */
	unsigned int lowest_lvl, highest_lvl, scnd_lowest_lvl, totl_lvl;	/* lowest and highest levels */
	//unsigned int actUpdCnt;	/* number of the updates that do generate some active nodes */ 

	/* graph topology/construction part */
	unsigned char *nArcMtx;				/* matrix: number of arcs that each node has */
	unsigned char *stMtx;				/* matrix: bit1=1 if having s arc; bit2=1 if having t arc */
	int	nTypes;							/* # of cells in typeCell input */
	#define MAXTYPES 10
	double *typeCells[MAXTYPES];		/* each points to a double matrix for one node type */
	int typeCellsSize[MAXTYPES];		/* each tells the number of edges of one cell type */
	unsigned char *typeMtx;				/* matrix: type of node in typeCells */
		
/* general macros */

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i,a) for (a = i->first.pointer, stopA = (i+1)->first.pointer; a != stopA; a++)

#define nNode( i ) ( (i) - nodes )
#define nArc( a )  ( ( a == NULL )? -1 : (a) - arcs )

void free_up();

#if _MSC_VER

#define Calloc_n_Check(var,type,units) \
{\
	if (var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
		throw std::bad_alloc(buffer);\
	}\
	var = (type*)calloc(units,sizeof(type));\
	if (!var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
		throw std::bad_alloc(buffer);\
	}\
}

#define Malloc_n_Check(var,type,units) \
{\
	if (var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
		throw std::bad_alloc(buffer);\
	}\
	var = (type*)malloc((units)*sizeof(type));\
	if (!var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
		throw std::bad_alloc(buffer);}\
}

#else

#define Calloc_n_Check(var,type,units) \
{\
	if (var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
		throw std::runtime_error(buffer);\
	}\
	var = (type*)calloc(units,sizeof(type));\
	if (!var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
		throw std::runtime_error(buffer);\
	}\
}

#define Malloc_n_Check(var,type,units) \
{\
	if (var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
		throw std::runtime_error(buffer);\
	}\
	var = (type*)malloc((units)*sizeof(type));\
	if (!var)\
	{\
		char buffer[100];\
		sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
		throw std::runtime_error(buffer);}\
}

#endif

/* 
   bucket macros:
   bucket's active node list is singly-linked
     operations aAdd, aRemove (from the front)
   bucket's inactive list is doubly-linked
     operations iAdd, iDelete (from arbitrary position)
*/

long i_dist;
node *i_next, *i_prev;


/* aAdd(l,i) from front */
#define aAdd(l,i)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  i_dist = i->d.dist;\
  if (i_dist < aMin)\
    aMin = i_dist;\
  if (i_dist > aMax)\
    aMax = i_dist;\
  if (dMax < aMax)\
    dMax = aMax;\
}


/* i must be the first element */
#define aRemove(l,i)\
{\
assert(l->firstActive != sentinelNode);\
l->firstActive = i->bNext;\
}

#define iAdd(l,i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l,i)\
{\
  i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = i_next;\
    i_next->bPrev = sentinelNode;\
  }\
  else {\
    i_prev = i->bPrev;\
    i_prev->bNext = i_next;\
    i_next->bPrev = i_prev;\
  }\
}

/* Bucket operations for nodes of different states */

#define aStateAdd(l,i)\
{\
i->bNext = l->firstActive;\
l->firstActive = i;\
}

#define iStateAdd(l,i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

/* Pixel bucket macro */
#define pxAdd(b,p,i)\
{\
	(p)->pNode = (i);\
	(p)->bNext = (b)->first;\
	(b)->first = (p);\
}

