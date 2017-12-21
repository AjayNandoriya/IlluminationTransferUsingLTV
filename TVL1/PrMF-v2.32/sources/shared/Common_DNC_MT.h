#pragma once

#define VERSION 2.32

/* ===== functions part ===== */

/* This file contains functions that are similar in both DNC and MT codes */
/* Function list:

1. procInput: process and check MATLAB mexfile input and output arguments, assign dimesion, #nodes/arcs, lambda, eps, etc.
2. printHeader.

*/

void procInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], int &_nDim, const int *&_dim)
{
	/* input argument meanings:
	nlhs: # of MATLAB funtion output arguments
	plhs: MATLAB output argument array pointer
	nlhs: # of MATLAB funtion input arguments
	plhs: MATLAB input argument array pointer

	_nDim: number of dimensions. 2 or 3, for 2D or 3D, respectively.
	_dim: dimension array, storing the data size in each dimension
	*/
	
	int ii;

	/* input/output argument # check */
	if (nrhs != N_INPUT_ARG) {mexPrintf("This function requires %d input arguements.", N_INPUT_ARG); mexErrMsgTxt("The number of input arguements is incorrect.");}
	if (nlhs != 1) {mexErrMsgTxt("This function requires 1 output arguement.");}

	/* input 1 */
	if (!mxIsClass(mtxImage_IN,strMATLABCLASS)) {mexPrintf("Input shoud be of type: %s\n",strMATLABCLASS); mexErrMsgTxt("The 1st input arguement has wrong type.");}
	if (mxIsSparse(mtxImage_IN)) {mexErrMsgTxt("The input matrix cannot be sparse.");}
	inputdata = (imgType *)mxGetData(mtxImage_IN);
	{
		if (( _nDim=mxGetNumberOfDimensions(mtxImage_IN)) > 3) {mexErrMsgTxt("The 1st input arguement (input image) must be either two or three dimensional.");}
		_dim = (int *)mxGetDimensions(mtxImage_IN);
		mm = _dim[0]; nn = _dim[1];
		if (_nDim==3) ll = _dim[2]; else ll = 1;
		if (mm<=1 || nn<=1) {mexErrMsgTxt("The 1st and 2nd dimensions of input image must be at least 2.");}
		_n = mm*nn*ll + 2;	/* # of nodes = # of pixels + |{s,t}| */
		_m = 0;				/* # of arcs is unknown */
	}
	
	/* input 2 */
	#ifdef MATRIX_Lambda
		if (!mxIsDouble(mtxLambda_IN)) {mexErrMsgTxt("The 2nd input arguement (Lambda) must be a double matrix.");}
		/* check dimension */
		if (mxGetNumberOfDimensions(mtxLambda_IN) != _nDim) {mexErrMsgTxt("The 1st and 2nd input arguements must have the same dimension.");}
		{
			const int * _LambdaDim;
			_LambdaDim = (int *)mxGetDimensions(mtxLambda_IN);
			if (_nDim==3)
			{
				if (_LambdaDim[0] != _dim[0] || _LambdaDim[1] != _dim[1] || _LambdaDim[2] != _dim[2]) {mexErrMsgTxt("The 1st and 2nd input arguements must have the same dimension.");}
			}
			else
			{
				if (_LambdaDim[0] != _dim[0] || _LambdaDim[1] != _dim[1]) {mexErrMsgTxt("The 1st and 2nd input arguements must have the same dimension.");}
			}
		}
		Lambda = (double *)mxGetData(mtxLambda_IN); lambda = 0; twolambda = 0;
		eps = EPSCONST;
	#else
		if (!mxIsDouble(dbLambda_IN)) {mexErrMsgTxt("The 2nd input arguement (Lambda) must be a double scalar.");}
		lambda = mxGetScalar(dbLambda_IN);
		twolambda = 2*lambda;
		if (lambda <= 0) {mexErrMsgTxt("The 2nd input arguement (Lambda) must be strictly positive.");}
		eps = (EPSCONST < lambda*0.01)?EPSCONST:(lambda*0.01);
	#endif

	/* input 3 */
	if (!mxIsDouble(type_IN)) {mexErrMsgTxt("The 3rd input arguement (neighborhood description) must be a double maxtrix.");}
	if ( mxGetNumberOfDimensions(type_IN) != 2) {mexErrMsgTxt("The 3rd input arguement (neighborhood description) must be two dimensional.");}
	n_edges = mxGetN(type_IN); edge_size = mxGetM(type_IN);
	if ( ll==1 && edge_size != 4) {mexErrMsgTxt("The 3rd input arguement (neighborhood description) must have 4 rows for 2D images.");}
	if ( ll>1 && edge_size != 5) {mexErrMsgTxt("The 3rd input arguement (neighborhood description) must have 5 rows for 3D images.");}
	if ( n_edges > 32) {mexErrMsgTxt("At most 32 out-going arcs are allowed for each pixel (can be extended later)!");}
	#ifndef NO_DISP
	mexPrintf("Each node is connected to %d neighbors.\n",n_edges);
	#endif
	neighborType = (double *)mxGetData(type_IN);
	{	/* store the absolute offset information in out_edges */
		double *edge;
		Malloc_n_Check(out_edges,edge_t,n_edges);
		if (ll==1)
		{
			for (edge=neighborType, ii=0; ii<n_edges; ii++, edge+=edge_size)
			{
				if (edge[0] != 1) {free_up();mexErrMsgTxt("Only type 1 (out-going) arcs are allowed in the neighborhood specification.");}
				if (((int)edge[2])!=edge[2] || ((int)edge[3])!=edge[3]) {free_up();mexErrMsgTxt("Rows 3 and 4 of the 3rd input arguement must use integer values.");}
				out_edges[ii].rel_pos = ((int)edge[2])+((int)edge[3])*mm;
				out_edges[ii].weight = edge[1];
			}
		}
		else
		{
			for (edge=neighborType, ii=0; ii<n_edges; ii++, edge+=edge_size)
			{
				if (edge[0] != 1) {free_up();mexErrMsgTxt("Only type 1 (out-going) arcs are allowed in the neighborhood specification.");}
				if (((int)edge[2])!=edge[2] || ((int)edge[3])!=edge[3] || ((int)edge[4])!=edge[4]) {free_up();mexErrMsgTxt("Rows 3, 4, and 5 of the 3rd input arguement must use integer values.");}
				out_edges[ii].rel_pos = ((int)edge[2])+(((int)edge[3])+((int)edge[4])*nn)*mm;
				out_edges[ii].weight = edge[1];
			}
		}
	}

	#ifdef MT_TV_L1
	/* input 4 */
		if (!mxIsDouble(thread_IN)) {mexErrMsgTxt("The 4th input arguement (thread switches) must be a double scalar.");}
		threadSwitch = (int)mxGetScalar(thread_IN);
		if (threadSwitch < 0 || threadSwitch > 3) {mexErrMsgTxt("The 4th input arguement must be (binary: UP/DN thread switch) 00, 01, 10, or 11.");}
	#endif

	return;
}

void printHeader()
{
	mexPrintf("\nParametric maximum flow algorithm (PrMF) for the %s model.\nWotao Yin, CAAM, Rice University, and Donald Goldfarb, IEOR, Columbia University, 2007.\n",PROBLEM_TYPE);
	mexPrintf("Version %4.2f\n\n", VERSION);

if (0){
	mexPrintf("Enabled options: no_source_sink_arcs");
#ifdef INIT_UPDATE
	mexPrintf(" / initial backward breadth-first search");
#endif
#ifdef GLOBAL_UPDATE
	mexPrintf(" / global backward breadth-first search with frequency");
#endif
#ifdef OFFSET_CURRENT
	mexPrintf(" / use \"short\" type of \"current\"");
#endif
#ifdef RESET_NODE_LABELS
	mexPrintf(" / reset all labels of T nodes to 1 after going to next level");
#endif
#ifdef SEP_STATES
	mexPrintf(" / separate the solutions of different states");
#endif

	mexPrintf(".\n");
}

	return;
}

/* definition of Macro annArc for inserting arcs with head, resCap, and cap */
/* add (ntail,nhead) forward arc, but leave "rev" empty */
#define annArc(ntail,nhead,fcap) \
	{\
		parc=(ntail)->current;\
		{	/* add new arc, but leave "rev" empty */\
			parc->head.offset=(nhead); parc->resCap=fcap;\
			(ntail)->current++;\
		}\
	}

#if defined(MT_TV_L1)
void CFlow_Graph::initialize_graph_TV()
#elif defined(DNC_TV_L2)
void initialize_graph_TV()
#else
#error no matched macro
#endif
{

	int i,j,k,ii,jj;
	int xx,yy,zz;
	int n_arcs;
	double *edge;
	node *pnode;
	arc *parc, *end_arc;
	bucket *l1;

/* allocate node space */
	Calloc_n_Check(nodes,node,n+1);
	source = nodes+(n-2);
	sink = nodes+(n-1);
	sentinelNode = nodes+n;

/* allocate buckets */
	Malloc_n_Check(buckets,bucket,n+1);
	for (l1 = buckets; l1 < buckets + (n+1); l1++)
	{	 l1 -> firstActive   = sentinelNode; l1 -> firstInactive  = sentinelNode;	}

/* determine # of out-going arcs for every node */
	pnode=nodes;
	pnode->first.offset=0;
	n_arcs=0;

	if (ll==1)
	{
		for (j=0;j<nn;j++)
		for (i=0;i<mm;i++)
		{
			for (edge=neighborType, ii=0; ii<n_edges; ii++,edge+=edge_size)
			{
				/* assuming there's only case 1 */
				pnode->d.neighbor <<= 1;		// left shift: add a bit for this edge
				xx = (int)edge[2]+i; yy = (int)edge[3]+j;
				if (xx>=0 && xx<mm && yy>=0 && yy<nn)
				{
					n_arcs++;
					pnode->d.neighbor++;	// set last bit to 1
				}
			}
			pnode++;
			pnode->first.offset=n_arcs;
			assert(pnode<=nodes+n);
		}	
	}
	else
	{
		for (k=0;k<ll;k++)
		for (j=0;j<nn;j++)
		for (i=0;i<mm;i++)
		{
			for (edge=neighborType, ii=0; ii<n_edges; ii++,edge+=edge_size)
			{
				/* assuming there's only case 1 */
				pnode->d.neighbor <<= 1;		// left shift: add a bit for this edge
				xx = (int)edge[2]+i; yy = (int)edge[3]+j; zz = (int)edge[4]+k;
				if (xx>=0 && xx<mm && yy>=0 && yy<nn && zz>=0 && zz<ll)
				{
					n_arcs++;
					pnode->d.neighbor++;	// set last bit to 1
				}
			}
			pnode++;
			pnode->first.offset=n_arcs;
			assert(pnode<=nodes+n);
		}
	}	/* now, n_arcs = # of (u,v) arcs */
	/* the last pnode was "source" */
	/* sink */
	pnode++;
	pnode->first.offset=n_arcs;
	/* the sentinelNode */
	pnode++;
	pnode->first.offset=(m = n_arcs);
	assert(pnode==nodes+n);


/* allocate arcs space */
	Calloc_n_Check(arcs,arc,m+1);

/* insert arcs */
/* for each pixel node v, no (s,v) or (v,t) arc is used */
	{
		/* add arcs for v nodes */
		/* leave "rev" empty, which is assigned later */
		for (pnode=nodes,jj=0; pnode!=source; pnode++,jj++)
		{
			pnode->current = arcs+pnode->first.offset;
			for(ii=n_edges-1; ii>=0 ; ii--)
			{
				if (pnode->d.neighbor & 1)
					annArc(pnode,(jj+out_edges[ii].rel_pos),out_edges[ii].weight);	// this step uses "current"
				pnode->d.neighbor >>= 1;	// right shift by 1
			}
			//pnode->current = arcs+pnode->first.offset;	// reset current to the first arc
		}
		
	}

/* set "rev" for all arcs of pixel nodes */
	{
		arc *parc_search;

		parc = arcs;
		for (ii=0,pnode=nodes; ii<n-2; ii++,pnode++)
		{
			// we have parc = arc+pnode->first.offset
			end_arc = arcs+(pnode+1)->first.offset;
			for (; parc<end_arc; parc++)
			{
				if (parc->rev.offset == 0)
				{
					/* this node's rev hasn't been assigned */
					parc_search = arcs+(nodes+parc->head.offset)->first.offset;
					while (parc_search->head.offset != ii)
					{
						parc_search++;
						assert(parc_search != arcs+(nodes+parc->head.offset+1)->first.offset);
					}
					parc->rev.offset = parc_search-arcs;
					parc_search->rev.offset = parc-arcs;
					
					assert(parc<=arcs+n_arcs);
				}
			}
		}
	}

/* save nodes and arcs to file */
/* not implemented, but recall that ->current is a pointer, not offset */	
	return;
}


/*	1. graph initialization for the fidelity terms
	2. converts offsets to pointers 
	3. calls initialize_sv_vt_cap, different for DNC_TV_L2 and MT_TV_L1, to get the residual cap/excess for each non-terminal node
*/
#if defined(MT_TV_L1)
void CFlow_Graph::initialize_fidelity()
#elif defined(DNC_TV_L2)
void initialize_fidelity()
#else
#error no matched macro
#endif
{
	arc *end_arc;
	arc *parc;
	bucket *l1;
	node *end_node;
	node *pnode;
	cType sv_minus_vt;

/* convert offsets to pointers */
/* contains initialization to the fidelity terms */
	{
	/* arcs part */
		end_arc = arcs+m;
		for (parc = arcs; parc < end_arc; parc++)
		{
			parc->head.pointer = nodes+parc->head.offset;
			parc->rev.pointer = arcs+parc->rev.offset;
		}
		assert(parc == arcs+m);

	/* nodes part */
		l1 = buckets+1;
		/* consider pixel nodes only */
		end_node = nodes+n-2;
		for (pnode = nodes; pnode < end_node; pnode++)
		{
			pnode->d.dist = 1;
			pnode->first.pointer = arcs + pnode->first.offset;
			pnode->current = pnode->first.pointer;
			sv_minus_vt = initialize_sv_vt_cap(pnode);
			if (sv_minus_vt > eps)
			{
				pnode->excess = sv_minus_vt;
				aAdd(l1,pnode);
				/* here, we are ignoring the excess at sink */
			}
			else if (sv_minus_vt <= eps)
			{
				pnode->excess = sv_minus_vt;
				iAdd(l1,pnode);
				/* here, we are ignoring the excess at sink */
			}
		}
		/* node s */
		assert(pnode == source);
		pnode->d.dist = n;
		pnode->current = pnode->first.pointer = arcs+pnode->first.offset;
		pnode->excess = 0;
		pnode++;
		/* node t */
		assert(pnode == sink);
		pnode->d.dist = 0;
		pnode->current = pnode->first.pointer = arcs+pnode->first.offset;
		pnode++;
		/* sentinel node */
		assert(pnode == sentinelNode);
		pnode->current = pnode->first.pointer = arcs+pnode->first.offset;
	
	}

	return;
}

#if defined(MT_TV_L1)
void CFlow_Graph::generalInit()
#elif defined(DNC_TV_L2)
void generalInit()
#else
#error no matched macro
#endif
{
	pushCnt  = 0;           /* number of pushes */
	relabelCnt   = 0;       /* number of relabels */
	updateCnt    = 0;       /* number of updates */
	gapCnt   = 0;           /* number of gaps */
	gNodeCnt = 0;           /* number of nodes after gap */
#ifdef SPECIAL_CNT
	pushCntSp  = 0;           /* number of pushes */
	relabelCntSp   = 0;       /* number of relabels */
	gapCntSp   = 0;           /* number of gaps */
	gNodeCntSp = 0;           /* number of nodes after gap */
#endif

	aMax = 0; aMin = n; dMax = 0;

	return;
}

/* Check whether the forward arcs over the cut are all saturated */
#if defined(MT_TV_L1)
void CFlow_Graph::checkMincut()
#elif defined(DNC_TV_L2)
void checkMincut()
#else
#error no matched macro
#endif
{
	node *i;
	arc *a;

	for(i=nodes;i!=source;i++)
	{
		if (i->d.dist >= n)
		{
			if (i->excess < -eps) {mexPrintf("ERROR: There exists an unsaturated node-t cut arc\n");return;}
			forAllArcs(i,a)
			{
				if(a->head.pointer->d.dist < n && a->resCap >= eps)
				{
					mexPrintf("ERROR: There exists an unsaturated cut arc: res_c(%d,%d)=%4.2f\n",i-nodes,a->head.pointer-nodes,a->resCap);return;
				}
			}
		}
	}
}

/* differ to relabel by adding "inputdata[i-nodes] = (imgType)last_lvl" below*/
#if defined(MT_TV_L1)
inline int CFlow_Graph::relabel(node *i)
#elif defined(DNC_TV_L2)
inline long relabel(node *i)
#else
#error no matched macro
#endif

{

	node  *j;
	long  minD;     /* minimum d of a node reachable from i */
	arc   *minA;    /* an arc which leads to the node with minimal d */
	arc   *a;
	bucket *ll = buckets+n;

	assert(i->excess > eps);

	relabelCnt++;

	minD = n;
	minA = NULL;

	/* find the minimum */
	forAllArcs(i,a)
	{
		//if (a -> resCap > 0) {
		if (a -> resCap > eps)
		{
			j = a -> head.pointer;
			if (j->d.dist < minD)
			{
				minD = j->d.dist;
				minA = a;
			}
		}
	}

	minD++;
	  
	if (minD < n) 
	{

		i->d.dist = minD;
		i->current = minA;

		if (dMax < minD) dMax = minD;
	}
	else
	{
		liftNodeToN(i);
	} /* end of minD < n */

#ifdef _DEBUG
	if (minD == n)
	{
		forAllArcs(i,a)
		{
			  assert(i->excess>-eps);
			  assert(a->head.pointer->d.dist == n || a->resCap < eps);
		}
	}
#endif

	return ( minD );

} /* end of relabel */

#if defined(MT_TV_L1)
inline int CFlow_Graph::gap(node *pnode)
#elif defined(DNC_TV_L2)
inline int gap(node *pnode)
#else
#error no matched macro
#endif
{
	bucket *l;
	bucket *ll=buckets+n;
	bucket *emptyB=buckets+pnode->d.dist;
	node  *i; 
	long  r;           /* index of the bucket before l  */
	int   cc;          /* cc = 1 if no nodes with positive excess before
		      		the gap */

	gapCnt ++;
	r = ((long)( emptyB - buckets )) - 1;

	/* set labels of nodes beyond the gap to "infinity" */
	for ( l = emptyB + 1; l <= buckets + dMax; l ++ ) {
	for ( i = l -> firstInactive; i != sentinelNode; i = i -> bNext ) 
	{
		liftNodeToN(i);
		i -> d.dist = n;
		gNodeCnt ++;
	}

	assert(l->firstActive == sentinelNode);

	l -> firstInactive = sentinelNode;
	}

	cc = ( aMin > r ) ? 1 : 0;

	dMax = r;
	aMax = r;

	liftNodeToN(pnode);

	return ( cc );
}

#if defined(MT_TV_L1)
inline void CFlow_Graph::discharge (node *i)
#elif defined(DNC_TV_L2)
inline void discharge (node *i)
#else
#error no matched macro
#endif
{

  node  *j;                 /* successor of i */
  long  jD;                 /* d of the next bucket */
  arc   *a;                 /* current arc (i,j) */
  cType  delta;
  arc *stopA;
  bucket *l;                /* i's bucket */
  bucket *lj;               /* j's bucket */
  bucket *ll = buckets+n;  /* Bucket n */


  assert(i->excess > eps);
  assert(i != sink);

  do {
	l = buckets + i->d.dist;		/* the bucket where the current node was in */

	/* scanning node-t arc: do nothing */
	assert(i->excess > -eps);

	/* scanning other arcs outgoing from  i  */
	jD = i->d.dist - 1;
	for (a = i->current, stopA = (i+1)->first.pointer; a != stopA; a++)
	{
		if (a -> resCap > eps)
		{
			j = a -> head.pointer;
			assert(j!=sink && j!=source);
			#ifdef DNC_TV_L2
				assert(i->state == j->state);
			#endif

			if (j->d.dist == jD)
			{
				pushCnt ++;
				if (i->excess >= a->resCap)
					delta = a->resCap;
				else
					delta = i->excess;

				a->resCap -= delta;
				a->rev.pointer->resCap += delta;

				if (j->excess < eps && j->excess + delta >= eps)
				{	lj = buckets + j->d.dist; iDelete(lj,j); aAdd(lj,j);	}
				j -> excess += delta;
				i -> excess -= delta;

				if (i->excess < eps) break;
			} /* a is admissible */
		} /* a is not saturated */
	} /* end of scanning arcs from i */

	if (a == stopA)
	{
		/* i must be relabeled */
		if ((l -> firstActive == sentinelNode) && (l -> firstInactive == sentinelNode))
		{   gap(i); }
		else
		{   relabel(i);}

		if (i->d.dist == n) break;
	}
	else
	{
		/* i no longer active */
		i->current = a;

		/* put i on inactive list */
		iAdd(l,i);

		break;
	}
  } while (1);
}


#if defined(MT_TV_L1)
void CFlow_Graph::stageOne ()
#elif defined(DNC_TV_L2)
void stageOne()
#else
#error no matched macro
#endif
{

	node   *i;
	bucket  *l;             /* current bucket */

	#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
  		globalUpdate ();
	#endif

	/* main loop */
	while ( aMax >= aMin )
	{
		l = buckets + aMax;
		i = l->firstActive;

		if (i == sentinelNode)
			aMax--;
		else
		{
			aRemove(l,i);
 			assert(i->excess > eps && i->d.dist < n);

			#ifdef _EXTRA_TIMING // post outgoing arc scan
				clock_t start = clock();
			#endif

			discharge(i);
			
			#ifdef _EXTRA_TIMING // post outgoing arc scan
				testTime2 += (clock() - start);
			#endif
		}

		#ifdef GLOBAL_UPDATE
			/* is it time for global update? */
			if (workSinceUpdate * GLOB_UPDT_FREQ > nm)
			{
				globalUpdate ();
				workSinceUpdate = 0;
			}
		#endif
	} /* end of the main loop */
}

#if defined(MT_TV_L1)
void CFlow_Graph::printIntermediateSolutions ()
#elif defined(DNC_TV_L2)
void printIntermediateSolutions()
#else
#error no matched macro
#endif
{
	if (n<10 && m<30)
	{	/* print node-arc structure */
		int kk;
		node *pnode;
		arc *parc;
		for(pnode=nodes;pnode!=source;pnode++)
		{
			kk=(int)(pnode-nodes);
			mexPrintf("Node #%d:tr[CE%4.2f]",kk, pnode->excess);
			if (pnode->d.dist < n)	mexPrintf("T\t");	else mexPrintf("S\t");
			forAllArcs(pnode,parc)
			{
				mexPrintf("#%d:(%d,%d)[C%4.2f]->#%d | ",parc-arcs,kk,parc->head.pointer-nodes,parc->resCap,parc->rev.pointer-arcs);
			}
			mexPrintf("\n");
		}
		/* print Bucket structure */
		bucket *ll;
		for (ll = buckets; ll <= buckets+n; ll++)
		{
			mexPrintf("Bucket %d Active: ", ll-buckets);
			for(pnode = ll->firstActive; pnode != sentinelNode; pnode=pnode->bNext)
				mexPrintf("%d ", (int)(pnode-nodes));
			mexPrintf("; Inactive: ", ll-buckets);
			for(pnode = ll->firstInactive; pnode != sentinelNode; pnode=pnode->bNext)
				mexPrintf("%d ", (int)(pnode-nodes));
			mexPrintf("\n");
		}
	}
}

/* performs a series of checks on the results of maxflow */
#if defined(MT_TV_L1)
void CFlow_Graph::checkFlow ()
#elif defined(DNC_TV_L2)
void checkFlow()
#else
#error no matched macro
#endif
{
	node *i;
	arc *a;
	bucket *l;

	/* check if you have a flow (pseudoflow) */
	/* check arc flows */
	forAllNodes(i)
	{
		forAllArcs(i,a)
		{
			if ((a->resCap <= (-eps))
			|| (a->rev.pointer->resCap <= (-eps))) {
				mexPrintf("ERROR: bad arc flow: #%d, %f %f\n",a-arcs,a->resCap,a->rev.pointer->resCap);
				mexPrintf("       node: #%d, excess: %f\n",i-nodes,i->excess);
				return;
			}
		}
	}

	/* check if there's any active node with d<n */
	for (l = buckets; l < buckets + n; l++) {
		if (l->firstActive != sentinelNode) {mexPrintf("ERROR: There exists an active node with d<n\n");return;}
	}

	/* check if mincut is saturated */
#ifdef PARAFLOW
	forAllNodes(i)
	{
		if (i->d.dist >= n)
		{
			if (i->excess < -eps) {mexPrintf("ERROR: There exists an unsaturated node-t cut arc\n");return;}
			forAllArcs(i,a)
			{
				if(a->head.pointer->d.dist < n && a->resCap >= eps) 
				{
					mexPrintf("ERROR: There exists an unsaturated cut arc: res_c(%d(d%d),%d(d%d))=%4.2f\n",i-nodes,i->d.dist,a->head.pointer-nodes,a->head.pointer->d.dist,a->resCap);return;
				}
			}
		}
	}
#else
	cType sum;
	aMax = dMax = 0;
	for (l = buckets; l < buckets + n; l++) {
		l->firstActive = sentinelNode;
		l->firstInactive = sentinelNode;
	}
	if (source->d.dist < n) {
		mexPrintf("ERROR: the solution is not optimal\n");
		return;
	}

	/* caculate total flow to t */
	sum=0;
	forAllArcs(sink,a)
	{
		sum += a->resCap;
	}
#endif
	
	mexPrintf("Solution checks (feasible and optimal)\n");
	return;
}
