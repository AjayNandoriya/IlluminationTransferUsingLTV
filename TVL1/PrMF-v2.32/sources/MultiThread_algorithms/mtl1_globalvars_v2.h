/* global variables */

/* timing part */
clock_t up_start_time, up_times[3];
clock_t dn_start_time, dn_times[3];
#ifdef _EXTRA_TIMING	// global def
clock_t up_updTime, up_mxfTime, up_tmpTime;
clock_t dn_updTime, dn_mxfTime, dn_tmpTime;
clock_t up_testTime1, up_testTime2;
clock_t dn_testTime1, dn_testTime2;
#endif
/* counters part */
unsigned long up_pushCnt,dn_pushCnt;			/* number of pushes */
unsigned long up_relabelCnt,dn_relabelCnt;		/* number of relabels */
unsigned long up_updateCnt,dn_updateCnt;		/* number of updates */
unsigned long up_gapCnt,dn_gapCnt;				/* number of gaps */
unsigned long up_gNodeCnt,dn_gNodeCnt;			/* number of nodes after gap */
unsigned long up_dNodeRelabelToN,dn_dNodeRelabelToN;		/* number of node-times checked for update, but node with label=n */
#ifdef LABEL_STATS	// def global stats
unsigned long up_minFinalLabel,dn_minFinalLabel;	/* minimum final distance label */
unsigned long up_maxFinalLabel,dn_maxFinalLabel;	/* maximum final distance label */
double up_avgFinalLabel,dn_avgFinalLabel;	/* sum of final distance labels */
#endif


#define EPSCONST 0.00001
double eps;					 /* small epsilon */

class threadStates
{
private:
	enum
	{
		NotScheduled,
		NotStarted,		
		Started,
		ParaMaxFlowFinished,
		Error
	} state;
	bool multithreaded;
	std::string message;

public:
	threadStates():state(NotStarted),multithreaded(false) { }

	void SetStarted(){state = Started;}
	void SetFinished(){state = ParaMaxFlowFinished;}
	void FoundError(){state = Error;}
	void SetMThreaded(){multithreaded = true;}
	void SetError(const char *msg){state=Error;message.append(msg);}

	void Clear() {state=NotStarted; multithreaded=false; message.erase(message.begin(),message.end());/*message.clear();*/}

	bool QFinished() const {return (state==ParaMaxFlowFinished);}
	bool QStarted() const {return (state!=NotStarted);}
	bool QMThreaded() const {return multithreaded;}
	bool QErrorFree() const {return (state!=Error);}
	const char *QMessage() const {return message.c_str();}

};


class CFlow_Graph
{
private:
	enum
	{
		Up,
		Down
	}	direction;			    /* up or down direction of \mu */

	node   *nodes;               /* array of nodes */
	arc    *arcs;                /* array of arcs */
	bucket *buckets;             /* array of buckets */
	node   *source;              /* source node pointer */
	node   *sink;                /* sink node pointer */

	node	*sentinelNode;        /* end of the node list marker */

#ifdef _EXTRA_TIMING	// in-class var definition
	clock_t testTime1, testTime2;	/* unspecified test timing variables */
#endif

	/* algorithm part */
	long   dMax;                 /* maximum distance label */
	long   aMax;                 /* maximum active node label */
	long   aMin;                 /* minimum active node label */
	double flow;                 /* flow value */
	/* counters part */
	unsigned long pushCnt;			   /* number of pushes */
	unsigned long relabelCnt;       /* number of relabels */
	unsigned long updateCnt;		   /* number of updates */
	unsigned long gapCnt;		       /* number of gaps */
	unsigned long gNodeCnt;           /* number of nodes after gap */
#ifdef LABEL_STATS	// in-class var definition
	/* final distance label of a node = the distance label of the node immediately before it's lifted to "n" */
	unsigned long minFinalLabel;	/* current minimum final distance label */
	unsigned long maxFinalLabel;	/* current maximum final distance label */
	unsigned long sumFinalLabel;	/* current sum of final distance labels */
	void finalLabelUpdate(unsigned long lab)
	{	if (lab<minFinalLabel)  minFinalLabel=lab;
		if (lab>maxFinalLabel)	maxFinalLabel=lab;
		sumFinalLabel += lab;
	}
#endif
#ifdef SPECIAL_CNT
	unsigned long pushCntSp;          /* number of pushes */
	unsigned long relabelCntSp;       /* number of relabels */
	unsigned long gapCntSp;           /* number of gaps */
	unsigned long gNodeCntSp;         /* number of nodes after gap */
#endif
	/* level part */
	unsigned int last_lvl;	/* last working levels */

	/* working variables */
#ifdef GLOBAL_UPDATE	// in-class definition
	unsigned int workSinceUpdate;
	unsigned long nm;
#endif
	long i_dist;
	node *i_next, *i_prev;
	arc	 *stopA;                /* used in forAllArcs */

public:
	long   n;                    /* number of nodes */
	long   m;                    /* number of arcs */
	unsigned current_lvl;		 /* current working level */

	CFlow_Graph();
	~CFlow_Graph();
	void generalInit();
	void initialize_graph_up();     /* initialize graph for up direction */
	void initialize_graph_dn();   /* initialize graph for down direction */
	void post_initialization();

	void para_max_flow_up();
	void para_max_flow_dn();

	void checkFlow();

	//void printgraph();
	void PrintStats();
	void TransStats();
	void printIntermediateSolutions ();

private:
	/* initialization part */
	void initialize_graph_TV();
	void initialize_fidelity();
	inline cType initialize_sv_vt_cap(node*);
	//inline cType __cdecl initialize_sv_vt_cap_dn(node*);
	/* network flow part */
	void liftNodeToN(node*);
	void globalUpdate();
	void update_network(int);
	void stageOne();
	void discharge(node *);
	int relabel(node *);
	int gap(node *);
	void checkMincut();

	//friend void CFlow_Graph::para_max_flow_up();
	//friend void CFlow_Graph::para_max_flow_down();
};

/* global variables of input data */

unsigned int G_current_lvl, G_rev_current_lvl;
threadStates UpThreadState, DnThreadState;
#ifdef _MSC_VER
	HANDLE  DnThreadRunMutex;	/* "Keep Running" mutex of the Down thread*/
	HANDLE  InitGraphMutex;		/* Initializing graph mutex */
	HANDLE  UpFinishInitEvent, DnFinishInitEvent;
#else
#endif
node *gpnodes;
arc *gparcs;
bucket* gpbuckets;

//unsigned int lvlCnt;		/* number of distinct levels */
unsigned int lowest_lvl, highest_lvl, scnd_lowest_lvl;	/* lowest and highest levels */

int _n;		/* total # of nodes, including s and t */
int _m;		/* total # of edges */
int mm;		/* # of rows of the input image */
int nn;		/* # of columns of the input image */
int ll;		/* # of 3rd dimension depth of the input image */

/* parametric update part */
pxBucket	*pxBuckets;		/* for organizing pixels by intensity */
pixel		*pixels;		/* elements of pxBuckets */

/* Matlab input and output */
#define mtxImage_IN		prhs[0]
#define dbLambda_IN		prhs[1]
#define mtxLambda_IN	prhs[1]
#define type_IN			prhs[2]
#define thread_IN		prhs[3]
#define mtxImage_OUT	plhs[0]

imgType *inputdata;				/* pointer of actual input image */
double lambda;					/* parameter lambda */
double twolambda;				/* = 2*lambda */
double *Lambda;					/* parameter Lambda if a matrix Lambda is given */
double *neighborType;				/* pointing to the uint32 matrix giving neighbor graph topology */
edge_t *out_edges;					/* absolute positions of neighbor edges */
int n_edges;					/* # of edges in neighborhood relationship */
int edge_size;					/* =4 for 2D images, =5 for 3D images */
int threadSwitch;
imgType *outputdata;			/* pointer of actual output image */


/* general macros */
using namespace std;

#ifdef _MSC_VER
	#define Calloc_n_Check(var,type,units) \
	{\
		if (var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
			throw bad_alloc(buffer);\
		}\
		var = (type*)calloc(units,sizeof(type));\
		if (!var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
			throw bad_alloc(buffer);\
		}\
	}

	#define Malloc_n_Check(var,type,units) \
	{\
		if (var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
			throw bad_alloc(buffer);\
		}\
		var = (type*)malloc((units)*sizeof(type));\
		if (!var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
			throw bad_alloc(buffer);}\
	}
#else
	#define Calloc_n_Check(var,type,units) \
	{\
		if (var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
			throw runtime_error(buffer);\
		}\
		var = (type*)calloc(units,sizeof(type));\
		if (!var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
			throw runtime_error(buffer);\
		}\
	}

	#define Malloc_n_Check(var,type,units) \
	{\
		if (var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Trying to allocate existing pointer "#var".\n");\
			throw runtime_error(buffer);\
		}\
		var = (type*)malloc((units)*sizeof(type));\
		if (!var)\
		{\
			char buffer[100];\
			sprintf(buffer,"Cannot allocate "#var" (%d units * %d bytes) = %d bytes.\n",sizeof(type),units,sizeof(type)*units);\
			throw runtime_error(buffer);}\
	}
#endif

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i,a) for (a = i->first.pointer, stopA = (i+1)->first.pointer; a != stopA; a++)

#define nNode( i ) ( (i) - nodes )
#define nArc( a )  ( ( a == NULL )? -1 : (a) - arcs )


/* 
   bucket macros:
   bucket's active node list is singly-linked
     operations aAdd, aRemove (from the front)
   bucket's inactive list is doubly-linked
     operations iAdd, iDelete (from arbitrary position)
*/


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
if (l->firstActive != sentinelNode)\
	l->firstActive->bPrev = sentinelNode;\
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

#define pushArc(i,j,a,delta,jD)\
{\
	pushCnt ++;\
	if (i->excess - a->resCap >= eps)\
		delta = a->resCap;\
	else\
		delta = i->excess;\
	a->resCap -= delta;\
	a->rev.pointer->resCap += delta;\
	if (j != sink){\
		if (j->excess < eps) \
		{\
			lj = buckets + jD;\
			iDelete(lj,j);\
			aAdd(lj,j);\
		}\
	}\
j -> excess += delta;\
i -> excess -= delta;\
if (i->excess < eps) break;\
}
