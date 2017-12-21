#pragma once

/* ====== type definition part ====== */

#ifdef DOUBLEFLOW
typedef	/* excess flow at nodes */
	double excessType;	/* change to double if needed */
typedef /* capacity type */
	double cType;	/* change to double if needed */
#else
typedef	/* excess flow at nodes */
	unsigned long excessType;	/* change to double if needed */

typedef /* capacity type */
	unsigned long cType;	/* change to double if needed */
#endif

/* the two unions below are used for saving graph topology to file */
/* before saving, pointer should be converted to offset */
/* after loading, offset should be converted to pointer */
union arc_p
{
	struct arcSt		*pointer;				/* pointer to arc */
	#if defined(__PPC64__) || defined(_WIN64)
		long long		offset;					/* relative offset to arcs */
		//#error 64-bit compiling
	#else
		long			offset;					/* relative offset to arcs */
	#endif
};

union node_p
{
	struct nodeSt	*pointer;				/* pointer to node */
	#if defined(__PPC64__) || defined(_WIN64)
		long long		offset;					/* relative offset to nodes */
		//#error 64-bit compiling
	#else
		long			offset;					/* relative offset to nodes */
	#endif
};

typedef  /* arc */
	struct arcSt
{
	cType           resCap;          /* residual capasity */
	union node_p    head;           /* arc head */
	union arc_p     rev;            /* reverse arc */
} arc;

typedef  /* node */
	struct nodeSt
{
	union arc_p     first;			  /* first outgoing arc */
	arc             *current;         /* current outgoing arc */
	excessType      excess;           /* excess at the node if >0, <0 residual cap to t */
	union dist_u
	{
		long		dist;               /* distance label */
		long		neighbor;			/* neighbor switches */
	}d;
	struct nodeSt   *bNext;           /* next node in bucket */
	struct nodeSt   *bPrev;           /* previous node in bucket */
    #ifdef DNC_TV_L2
	   unsigned short	state;
    #endif
} node;

typedef /* bucket */
	struct bucketSt
{
	node             *firstActive;      /* first node with positive excess */
	node             *firstInactive;    /* first node with zero excess */
} bucket;

#ifdef MT_TV_L1
	typedef /* pixel */
		struct pixelSt
	{
		int idxNode;
		struct pixelSt	*bNext;
	} pixel;

	typedef /* bucket of pixels */
		struct pxBucketSt
	{
		pixel	*first;		/* pointing to the first pixel in bucket */
	} pxBucket;
#endif

typedef /* edge */
	struct edgeSt
{
	double weight;
	int 	   rel_pos;
} edge_t;

