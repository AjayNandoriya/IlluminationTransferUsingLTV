/* This file contains the core functions for computing a parametric maximum flow */

inline void CFlow_Graph::liftNodeToN(node* pnode)
{
	pnode->d.dist = n;
	outputdata[pnode-nodes]=(imgType)last_lvl;
}

void CFlow_Graph::PrintStats()
{
	mexPrintf ("Statistics:\n");
	mexPrintf ("c pushes:      %10ld\n", pushCnt);
	mexPrintf ("c relabels:    %10ld\n", relabelCnt);
	mexPrintf ("c gaps:        %10ld\n", gapCnt);
	mexPrintf ("c gap nodes:   %10ld\n", gNodeCnt);
	mexPrintf ("c updates:     %10ld\n", updateCnt);
	mexPrintf ("\n");

	return;
}

/* global update via backward breadth first search from the sink */
void CFlow_Graph::globalUpdate ()
{

  node  *i, *j;       /* node pointers */
  arc   *a;           /* current arc pointers  */
  bucket *l, *jL;          /* bucket */
  long curDist, jD;
  long state;

  /* initialization */

  forAllNodes(i)
    i -> d.dist = n;
  sink -> d.dist = 0;

  // clear buckets
  for (l = buckets; l <= buckets + dMax; l++) {
    l -> firstActive   = sentinelNode;
    l -> firstInactive  = sentinelNode;
  }

  dMax = aMax = 0;
  aMin = n;

  /* breadth first search */

  // add sink to bucket zero

  iAdd(buckets, sink);

  // main loop
  for (curDist = 0; 1; curDist++) {

    state = 0;
    l = buckets + curDist;
    jD = curDist + 1;	// tmp distance for j
    jL = l + 1;			// tmp bucket for j

	// if current bucket is empty before any processing, global update terminates
    if ((l->firstActive == sentinelNode) && 
	(l->firstInactive == sentinelNode))
      break;

    while (1) {

      switch (state) {
      case 0: 
	i = l->firstInactive;	// take out first inactive
	state = 1;
	break;
      case 1:
	i = i->bNext;			// take out next inactive
	break;
      case 2:
	i = l->firstActive;		// take out first active
	state = 3;
	break;
      case 3:
	i = i->bNext;			// take out next active
	break;
      default: 
	assert(0);
	break;
      }
      
    // after finishing all inactive nodes, go to active bucket
    if (i == sentinelNode) 
	{
		if (state == 1) {
			state = 2;
			continue;
		}
		else {
			assert(state == 3);
			break;
		}
    }

      /* scanning arcs incident to node i */
      forAllArcs(i,a) {
	if (a->rev.pointer->resCap > 0 ) {
	  j = a->head.pointer;
	  if (j->d.dist == n) {
	    j->d.dist = jD;
	    j->current = j->first.pointer;
	    if (jD > dMax) dMax = jD;
	    
	    if (j->excess > 0) {
	      /* put into active list */
	      aAdd(jL,j);
	    }
	    else {
	      /* put into inactive list */
	      iAdd(jL,j);
	    }
	  }
	}
      } /* node i is scanned */ 
    }
  }

} /* end of global update */

/* Update graph for the input FF's in pxBuckets[last_lvl] now switched from 1 to 0 */
void CFlow_Graph::update_network(int level)	/* this level is the LAST level */
{
	node *pNode;
	pixel *pPixel;
	arc *pArc;
	bucket *buk;

	updateCnt++;

	assert(pxBuckets[level].first != NULL);
	for (pPixel=pxBuckets[level].first; pPixel != NULL; pPixel=pPixel->bNext)
	{
		pNode = nodes+pPixel->idxNode;
		if (pNode->d.dist < n)	/* then this node is not with s (labeled as "n") yet */
		{
			#ifdef MATRIX_Lambda
				twolambda = 2*Lambda[pNode-nodes];
			#endif
		
			/* decrase c(v,t) by 2*lambda */
			pArc = pNode->first.pointer;
			if (pNode->excess < eps && pNode->excess+twolambda > eps)
			{
				buk = buckets + pNode->d.dist;
				iDelete(buk,pNode);
				aAdd(buk,pNode);
			}
			pNode->excess += twolambda;

		}

	}
}

/* main parametric max-flow iteraion */
void CFlow_Graph::para_max_flow_up()
{
	node *i;
	bucket *l;

	assert(direction == Up);
	UpThreadState.SetStarted();

	last_lvl = lowest_lvl;
	for (current_lvl=lowest_lvl+1; current_lvl <= highest_lvl; current_lvl++)
	{
		if (pxBuckets[current_lvl].first != NULL)
		{
			if (current_lvl >= G_rev_current_lvl)
			{
				//mexPrintf ("! Up   para-max-flow stopped at levels (current=%d, last=%d)\n", current_lvl, last_lvl);
				break;
			}

			//mexPrintf("Up thread: process level %d\n", current_lvl);
			G_current_lvl = current_lvl;

			#ifdef _EXTRA_TIMING
			up_tmpTime = clock();
			#endif

			/* parametric network update */
			update_network(last_lvl);

			#ifdef _EXTRA_TIMING
			up_updTime += clock()-up_tmpTime;
			up_tmpTime = clock();
			#endif

			/* preflow algorithm */
			stageOne();

			#ifdef _EXTRA_TIMING
			up_mxfTime += clock()-up_tmpTime;
			#endif


			#ifdef _DEBUG
			/* check min-cut */
			checkMincut();
			#endif

			last_lvl = current_lvl;
		}
	}

	UpThreadState.SetFinished();

	if(DnThreadState.QStarted() && DnThreadState.QErrorFree())
	{	/* If the other thread has started, need to wait for it to finish before signing value
		   to the nodes that are not yet assigned values by either threads */
#ifdef _MSC_VER
		WaitForSingleObject( DnThreadRunMutex, INFINITE );
#else
#endif
	
		for (l=buckets+1; l<=buckets+dMax;l++)	/* for all nodes with dist labels = 1 to dMax */
		{
			for (i=l->firstInactive; i!= sentinelNode; i=i->bNext)
				if (outputdata[i-nodes] < (imgType)last_lvl)
					outputdata[i-nodes] = (imgType)last_lvl;
		}

#ifdef _MSC_VER
		ReleaseMutex( DnThreadRunMutex );
#else
#endif
	}
	else
	{	/* However, if the other thread has not started, we do not need to wait */
		/* Assign intensity values to pixels still in T in the output image */
		assert(current_lvl == highest_lvl+1);
		for (l=buckets+1; l<=buckets+dMax;l++)	/* for all nodes with dist labels = 1 to dMax */
		{
			for (i=l->firstInactive; i!= sentinelNode; i=i->bNext)
				outputdata[i-nodes] = (imgType)last_lvl;
		}
	}
	return;
}

/* main parametric max-flow iteraion */
void CFlow_Graph::para_max_flow_dn()
{
	assert(direction == Down);

#ifdef _MSC_VER
	DnThreadRunMutex = CreateMutex( NULL, TRUE, NULL );
#else
#endif
	DnThreadState.SetStarted();

	last_lvl = highest_lvl+1;
	for (current_lvl=highest_lvl; current_lvl > lowest_lvl; current_lvl--)
	{
		if (pxBuckets[current_lvl].first != NULL)
		{
			if (current_lvl <= G_current_lvl)
			{
				//mexPrintf ("! Down para-max-flow stopped at levels (current=%d, last=%d)\n", current_lvl, last_lvl);
				break;
			}

			//mexPrintf("Down thread: process level %d\n", current_lvl);
			G_rev_current_lvl = current_lvl;

			#ifdef _EXTRA_TIMING
			dn_tmpTime = clock();
			#endif

			/* IMPORTANT: last_lvl is used in the max-flow functions to set output data */
			last_lvl = current_lvl;
			/* parametric network update */
			update_network(last_lvl);

			#ifdef _DEBUG
			mexPrintf("---End of Network Update level %d---\n",G_rev_current_lvl);
			printIntermediateSolutions();
			#endif

			#ifdef _EXTRA_TIMING
			dn_updTime += clock()-dn_tmpTime;
			dn_tmpTime = clock();
			#endif

			/* preflow algorithm */
			stageOne();

			#ifdef _DEBUG
			mexPrintf("---End of Max-Flow level %d---\n",G_rev_current_lvl);
			printIntermediateSolutions();
			#endif

			#ifdef _EXTRA_TIMING
			dn_mxfTime += clock()-dn_tmpTime;
			#endif


			#ifdef _DEBUG
			/* check min-cut */
			checkMincut();
			#endif
		}

	}

	DnThreadState.SetFinished();

	if( ! (UpThreadState.QStarted() && UpThreadState.QErrorFree()) )	/* if the UP thread has started and errorfree, leave the final job to UP thread */
	{
		/* assign intensity values to pixels in the output image */
		assert(current_lvl == lowest_lvl);
		bucket *l;
		node *i;

		for (l=buckets+1; l<=buckets+dMax;l++)	/* for all nodes with dist labels = 1 to dMax */
		{
			for (i=l->firstInactive; i!= sentinelNode; i=i->bNext)
				outputdata[i-nodes] = (imgType)lowest_lvl;
		}
	}

#ifdef _MSC_VER
	ReleaseMutex( DnThreadRunMutex );
#else
#endif
	return;
}

/* initialize_sv_vt_cap and initialize_sv_vt_cap_db are the same */
/* specific for the UP thread of MT_TV_L1 */
inline cType CFlow_Graph::initialize_sv_vt_cap(node* pnode)
{
	/*	c(s,v) = lambda * (1 - 1{f_ij >= \mu} )
		c(v,t) = lambda * 1{f_ij >= \mu}

		for the lowest \mu, therefore, c(s,v)=0 and c(v,t) = lambda
	*/
	
	#ifdef MATRIX_Lambda
		lambda = Lambda[pnode-nodes];
		if (lambda < 0) {throw std::runtime_error("Found negative Lambda value");}
		return(-lambda); 
	#else
		return(-lambda);
	#endif
}

void CFlow_Graph::initialize_graph_up()
{
	direction = Up;
	
	/* finish the TV part of graph */
	/* pointers are stored as offsets */
	initialize_graph_TV();

	/* pass the graph topology to DN thread through global variables */
#ifdef _MSC_VER
	if (threadSwitch == 3)
	{
		_m = m;
		gpnodes = nodes;
		gparcs = arcs;
		gpbuckets = buckets;
		
		SetEvent(UpFinishInitEvent);
		WaitForSingleObject(DnFinishInitEvent, INFINITE);
	}
#endif

	/* convert offsets to pointers */
	/* contains initialization */
	initialize_fidelity();

	post_initialization();

	return;
}

void CFlow_Graph::initialize_graph_dn()
{
	direction = Down;

//#ifdef _MSC_VER
//	if (threadSwitch == 3)
//	{
//		/* wait for Up thread to finish initialization first */
//		WaitForSingleObject(UpFinishInitEvent, INFINITE);
//
//		m=_m;
//		assert(m>0 && gpnodes!=NULL && gparcs!=NULL);
//
//		/* allocate space */
//		Malloc_n_Check(nodes,node,n+1);
//		memcpy(nodes,gpnodes,(n+1)*sizeof(node));
//		source = nodes+(n-2);
//		sink = nodes+(n-1);
//		sentinelNode = nodes+n;
//
//		Malloc_n_Check(arcs,arc,m+1);
//		memcpy(arcs,gparcs,(m+1)*sizeof(arc));
//
//		/* allocate buckets */
//		Malloc_n_Check(buckets,bucket,n+1);
//		for (bucket *l1 = buckets; l1 < buckets + (n+1); l1++)
//		{	 l1 -> firstActive   = sentinelNode; l1 -> firstInactive  = sentinelNode;	}
//
//		/* finished coping data */
//		SetEvent(DnFinishInitEvent);
//	}
//	else
//#endif
	{
		/* finish the TV part of graph */
		/* pointers are stored as offsets */
		initialize_graph_TV();
	}

	/* convert offsets to pointers */
	/* contains initialization */
	initialize_fidelity();

	return;
}

void CFlow_Graph::TransStats()
{
	if (direction == Up)
	{
		up_pushCnt=pushCnt;         
		up_relabelCnt=relabelCnt;   
		up_updateCnt=updateCnt;     
		up_gapCnt=gapCnt;           
		up_gNodeCnt=gNodeCnt;       

#ifdef _EXTRA_TIMING // transfer UP stats
		up_testTime1 = testTime1;
		up_testTime2 = testTime2;
#endif
#ifdef CYCLE_DSCVR	// transfer UP stats
		up_cycleCnt = cycleCnt;
		up_abandonCnt = abandonCnt;
#endif
#ifdef ISOLATED_DSCVR // transfer UP stats
		up_isoCallCnt = isoCallCnt;
		up_isoFailCnt = isoFailCnt;
		up_isoFailNodeCnt = isoFailNodeCnt;
		up_isoFailNodeMinCnt = isoFailNodeMinCnt;
		up_isoFailNodeMaxCnt = isoFailNodeMaxCnt;
		up_isoSuccCnt = isoSuccCnt;
		up_isoSuccNodeCnt = isoSuccNodeCnt;
		up_isoSuccNodeMaxCnt = isoSuccNodeMaxCnt;
	#ifdef _EXTRA_TIMING
		up_searchTotalTime = searchTotalTime;
		up_searchSuccTime = searchSuccTime;
		up_searchFailTime = searchFailTime;
	#endif
#endif
#ifdef LABEL_STATS	// transfer UP stats
		up_minFinalLabel = minFinalLabel;
		up_maxFinalLabel = maxFinalLabel;
		up_avgFinalLabel = ((double)sumFinalLabel)/(n-2);
#endif
	}
	else
	{
		dn_pushCnt=pushCnt;         
		dn_relabelCnt=relabelCnt;   
		dn_updateCnt=updateCnt;     
		dn_gapCnt=gapCnt;           
		dn_gNodeCnt=gNodeCnt;       
#ifdef _EXTRA_TIMING // transfer DN stats
		dn_testTime1 = testTime1;
		dn_testTime2 = testTime2;
#endif
#ifdef CYCLE_DSCVR	// transfer DN stats
		dn_cycleCnt = cycleCnt;
		dn_abandonCnt = abandonCnt;
#endif
#ifdef ISOLATED_DSCVR // transfer DN stats
		dn_isoCallCnt = isoCallCnt;
		dn_isoFailCnt = isoFailCnt;
		dn_isoFailNodeCnt = isoFailNodeCnt;
		dn_isoFailNodeMinCnt = isoFailNodeMinCnt;
		dn_isoFailNodeMaxCnt = isoFailNodeMaxCnt;
		dn_isoSuccCnt = isoSuccCnt;
		dn_isoSuccNodeCnt = isoSuccNodeCnt;
		dn_isoSuccNodeMaxCnt = isoSuccNodeMaxCnt;
	#ifdef _EXTRA_TIMING
		dn_searchTotalTime = searchTotalTime;
		dn_searchSuccTime = searchSuccTime;
		dn_searchFailTime = searchFailTime;
	#endif
#endif
#ifdef LABEL_STATS	// transfer DN stats
		dn_minFinalLabel = minFinalLabel;
		dn_maxFinalLabel = maxFinalLabel;
		dn_avgFinalLabel = ((double)sumFinalLabel)/(n-2);
#endif
	}
	return;
}

CFlow_Graph::CFlow_Graph()
{
	nodes=NULL;
	arcs=NULL;
	buckets=NULL;

	n=0;
	m=0;

#ifdef _EXTRA_TIMING
	testTime1 = 0;
	testTime2 = 0;
#endif

	flow = 0;

	return;
}

CFlow_Graph::~CFlow_Graph()
{
	if (nodes) free(nodes);
	if (arcs) free(arcs);
	if (buckets) free(buckets);

	nodes=NULL;
	arcs=NULL;
	buckets=NULL;

	return;
}

/* free memory */
void free_up()
{
	if (pixels) free(pixels);
	if (pxBuckets) free(pxBuckets);
	if (out_edges) free(out_edges);

	pixels=NULL; pxBuckets=NULL; out_edges=NULL;

	UpThreadState.Clear();
	DnThreadState.Clear();
}

#ifdef _MSC_VER
unsigned __stdcall UpThread( void* pArguments )
#else
void *UpThread( void* pArguments )
#endif
{
	try
	{
		CFlow_Graph G;			/* the main graph */
		try
		{
			G.n = _n; G.m = _m;

			up_start_time = clock();	/* timing */

			if (DnThreadState.QFinished()) goto EndUpThread;	// if the other thread is almost finished, stop this thread

			G.generalInit();

			up_times[0] = clock();	/* timing */

			if (DnThreadState.QFinished()) goto EndUpThread;	// if the other thread is almost finished, stop this thread

			G.initialize_graph_up();			/* build main graph */
			G.post_initialization();

			up_times[1] = clock();	/* timing */

		#ifdef _DEBUG	
			G.printIntermediateSolutions();
		#endif

			if (DnThreadState.QFinished()) goto EndUpThread;	// if the other thread is almost finished, stop this thread

			G.para_max_flow_up();

			up_times[2] = clock();	/* timing */

		#ifdef _DEBUG
			G.checkFlow();
		#endif

		}
		catch(const std::bad_alloc &ba)
		{
			up_times[0] = up_times[1] = up_times[2] = clock();

			UpThreadState.SetError(ba.what());

			char buffer[100];
			sprintf(buffer,"Graph has %d nodes and %d arcs, and needs %d bytes of memory.\n", G.n, G.m,sizeof(node)*(G.n+1)+sizeof(bucket)*(G.n+1)+sizeof(arc)*(G.m+1)+sizeof(cType)*G.m);
			UpThreadState.SetError(buffer);
		}
		catch(const std::exception &ba)
		{
			up_times[0] = up_times[1] = up_times[2] = clock();

			UpThreadState.SetError(ba.what());

			char buffer[100];
			sprintf(buffer,"Graph has %d nodes and %d arcs, and needs %d bytes of memory.\n", G.n, G.m,sizeof(node)*(G.n+1)+sizeof(bucket)*(G.n+1)+sizeof(arc)*(G.m+1)+sizeof(cType)*G.m);
			UpThreadState.SetError(buffer);
		}
		G.TransStats();
	}
	catch(const std::exception &ba)
	{
		up_times[0] = up_times[1] = up_times[2] = clock();

		UpThreadState.SetError(ba.what());
	}



EndUpThread:
	if (UpThreadState.QMThreaded())
	{	// if this is a standalone thread 
		#ifdef _MSC_VER
			_endthreadex( 0 );
			return 0;
		#else
			return NULL;
		#endif
	}
	else
		return NULL;
}

#ifdef _MSC_VER
unsigned __stdcall DnThread( void* pArguments )
#else
void *DnThread( void* pArguments )
#endif
{
	try
	{
		CFlow_Graph G_rev;	/* the main graph with reversed edges and terminal nodes*/
		try
		{
			G_rev.n = _n; G_rev.m = _m;

			dn_start_time = clock();	/* timing */

			if (UpThreadState.QFinished()) goto EndDnThread;	// if the other thread is almost finished, stop this thread

			G_rev.generalInit();

			dn_times[0] = clock();	/* timing */

			if (UpThreadState.QFinished()) goto EndDnThread;	// if the other thread is almost finished, stop this thread

			G_rev.initialize_graph_dn();	/* build the reversed graph */
			G_rev.post_initialization();

			dn_times[1] = clock();	/* timing */

		#ifdef _DEBUG	
			G_rev.printIntermediateSolutions();
		#endif

			if (UpThreadState.QFinished()) goto EndDnThread;	// if the other thread is almost finished, stop this thread

			G_rev.para_max_flow_dn();

			dn_times[2] = clock();	/* timing */

		#ifdef _DEBUG
			G_rev.checkFlow();
		#endif
		}
		catch(const std::bad_alloc &ba)
		{
			dn_times[0] = dn_times[1] = dn_times[2] = clock();

			DnThreadState.SetError(ba.what());

			char buffer[100];
			sprintf(buffer,"Graph has %d nodes and %d arcs, and needs %d bytes of memory.\n", G_rev.n, G_rev.m,sizeof(node)*(G_rev.n+1)+sizeof(bucket)*(G_rev.n+1)+sizeof(arc)*(G_rev.m+1)+sizeof(cType)*G_rev.m);
			DnThreadState.SetError(buffer);
		}
		catch(const std::exception &ba)
		{
			dn_times[0] = dn_times[1] = dn_times[2] = clock();

			DnThreadState.SetError(ba.what());

			char buffer[100];
			sprintf(buffer,"Graph has %d nodes and %d arcs, and needs %d bytes of memory.\n", G_rev.n, G_rev.m,sizeof(node)*(G_rev.n+1)+sizeof(bucket)*(G_rev.n+1)+sizeof(arc)*(G_rev.m+1)+sizeof(cType)*G_rev.m);
			DnThreadState.SetError(buffer);
		}

		G_rev.TransStats();
	}
	catch(const std::exception &ba)
	{
		dn_times[0] = dn_times[1] = dn_times[2] = clock();

		DnThreadState.SetError(ba.what());
	}

EndDnThread:
	if (DnThreadState.QMThreaded())
	{	// if this is a standalone thread 
		#ifdef _MSC_VER
			_endthreadex( 0 );
			return 0;
		#else
			return NULL;
		#endif
	}
	else
		return NULL;
}

void printStatistics(clock_t &mex_start)
{
	double ticks = 1/(double)CLOCKS_PER_SEC;

	mexPrintf("Main loop total time  = %f sec\n\n",((double)(clock()-mex_start))*ticks);

	if (UpThreadState.QStarted())
	{
		if (!UpThreadState.QErrorFree())
			mexPrintf("UP Thread Preflow algorithm ERROR:\n%s\n",UpThreadState.QMessage());

		mexPrintf ("UP Thread Preflow algorithm statistics:\n");
		mexPrintf (" updates:     %10ld\n", up_updateCnt);
		mexPrintf (" pushes:      %10ld\n", up_pushCnt);
		mexPrintf (" relabels:    %10ld\n", up_relabelCnt);
		mexPrintf (" gaps:        %10ld\n", up_gapCnt);
		mexPrintf (" gap nodes:   %10ld\n", up_gNodeCnt);
		mexPrintf (" rlb nodes:   %10ld\n", up_dNodeRelabelToN);
#ifdef CYCLE_DSCVR	// display UP stats
		mexPrintf (" cycles found:%10ld\n", up_cycleCnt);
		mexPrintf (" abandoned:   %10ld\n", up_abandonCnt);
#endif
#ifdef ISOLATED_DSCVR // display UP stats
		mexPrintf (" iso calls:   %10ld\n", up_isoCallCnt);
		mexPrintf (" iso fails:   %10ld\n", up_isoFailCnt);
		mexPrintf (" iso f_nodes: %10ld\n", up_isoFailNodeCnt);
		mexPrintf (" iso f_min:   %10ld\n", up_isoFailNodeMinCnt);
		mexPrintf (" iso f_max:   %10ld\n", up_isoFailNodeMaxCnt);
		mexPrintf (" iso succs:   %10ld\n", up_isoSuccCnt);
		mexPrintf (" iso s_nodes: %10ld\n", up_isoSuccNodeCnt);
		mexPrintf (" iso s_max:   %10ld\n", up_isoSuccNodeMaxCnt);
#endif
#ifdef LABEL_STATS	// display UP stats 
		mexPrintf (" min fn_dist: %10ld\n", up_minFinalLabel);
		mexPrintf (" max fn_dist: %10ld\n", up_maxFinalLabel);
		mexPrintf (" avg fn_dist: %10f\n", up_avgFinalLabel);
#endif
		mexPrintf (" graph init time      = %f sec\n", ((double)(up_times[0]-up_start_time))*ticks);
		mexPrintf (" preflow init time    = %f sec\n", ((double)(up_times[1]-up_times[0]))*ticks);
		mexPrintf (" total para-max-flow  = %f sec\n", ((double)(up_times[2]-up_times[1]))*ticks);
		#ifdef _EXTRA_TIMING
		mexPrintf ("  this includes:\n");
		mexPrintf ("  %5ld preflows time = %f sec\n", up_updateCnt, ((double)(up_mxfTime))*ticks);
		mexPrintf ("  %5ld updates time  = %f sec\n", up_updateCnt, ((double)(up_updTime))*ticks);
		mexPrintf ("  test timing 1      = %f sec\n", ((double)up_testTime1)*ticks);
		mexPrintf ("  test timing 2      = %f sec\n", ((double)up_testTime2)*ticks);
		mexPrintf ("  iso totl tm: %f sec\n",  ((double)up_searchTotalTime)*ticks);
		mexPrintf ("  iso succ tm: %f sec\n",  ((double)up_searchSuccTime)*ticks);
		mexPrintf ("  iso fail tm: %f sec\n",  ((double)up_searchFailTime)*ticks);
		#endif
		mexPrintf ("\n");
	}

	if (DnThreadState.QStarted())
	{
		if (!DnThreadState.QErrorFree())
			mexPrintf("DN Thread Preflow algorithm ERROR:\n%s\n",DnThreadState.QMessage());

		mexPrintf ("DN Thread Preflow algorithm statistics:\n");
		mexPrintf (" updates:     %10ld\n", dn_updateCnt);
		mexPrintf (" pushes:      %10ld\n", dn_pushCnt);
		mexPrintf (" relabels:    %10ld\n", dn_relabelCnt);
		mexPrintf (" gaps:        %10ld\n", dn_gapCnt);
		mexPrintf (" gap nodes:   %10ld\n", dn_gNodeCnt);
		mexPrintf (" rlb nodes:   %10ld\n", dn_dNodeRelabelToN);
#ifdef CYCLE_DSCVR	// display DN stats
		mexPrintf (" cycles found:%10ld\n", dn_cycleCnt);
		mexPrintf (" abandoned:   %10ld\n", dn_abandonCnt);
#endif
#ifdef ISOLATED_DSCVR // display DN stats
		mexPrintf (" iso calls:   %10ld\n", dn_isoCallCnt);
		mexPrintf (" iso fails:   %10ld\n", dn_isoFailCnt);
		mexPrintf (" iso f_nodes: %10ld\n", dn_isoFailNodeCnt);
		mexPrintf (" iso f_min:   %10ld\n", dn_isoFailNodeMinCnt);
		mexPrintf (" iso f_max:   %10ld\n", dn_isoFailNodeMaxCnt);
		mexPrintf (" iso succs:   %10ld\n", dn_isoSuccCnt);
		mexPrintf (" iso s_nodes: %10ld\n", dn_isoSuccNodeCnt);
		mexPrintf (" iso s_max:   %10ld\n", dn_isoSuccNodeMaxCnt);
#endif
#ifdef LABEL_STATS	// display DN stats
		mexPrintf (" min fn_dist: %10ld\n", dn_minFinalLabel);
		mexPrintf (" max fn_dist: %10ld\n", dn_maxFinalLabel);
		mexPrintf (" avg fn_dist: %10f\n", dn_avgFinalLabel);
#endif
		mexPrintf (" graph init time      = %f sec\n", ((double)(dn_times[0]-dn_start_time))*ticks);
		mexPrintf (" preflow init time    = %f sec\n", ((double)(dn_times[1]-dn_times[0]))*ticks);
		mexPrintf (" total para-max-flow  = %f sec\n", ((double)(dn_times[2]-dn_times[1]))*ticks);
		#ifdef _EXTRA_TIMING
		mexPrintf ("  this includes:\n");
		mexPrintf ("  %5ld preflows time = %f sec\n", dn_updateCnt, ((double)(dn_mxfTime))*ticks);
		mexPrintf ("  %5ld updates time  = %f sec\n", dn_updateCnt, ((double)(dn_updTime))*ticks);
		mexPrintf ("  test timing 1      = %f sec\n", ((double)dn_testTime1)*ticks);
		mexPrintf ("  test timing 2      = %f sec\n", ((double)dn_testTime2)*ticks);
		mexPrintf ("  iso totl tm: %f sec\n",  ((double)dn_searchTotalTime)*ticks);
		mexPrintf ("  iso succ tm: %f sec\n",  ((double)dn_searchSuccTime)*ticks);
		mexPrintf ("  iso fail tm: %f sec\n",  ((double)dn_searchFailTime)*ticks);
		#endif
		mexPrintf ("\n");
	}

	return;
}

/* add pixels to buckets according to their input intensity */
void initLevelBuckets()
{
	int ii;

	Calloc_n_Check(pxBuckets,pxBucket,(GRAYLEVEL));
	Calloc_n_Check(pixels,pixel,(mm*nn*ll));

/* Pixel bucket macro */
#define pxAdd(b,p,i)\
{\
	(p)->idxNode = (i);\
	(p)->bNext = (b)->first;\
	(b)->first = (p);\
}

	lowest_lvl=GRAYLEVEL-1;highest_lvl=0;

	for (ii=0;ii<mm*nn*ll;ii++)
	{
		if (lowest_lvl>inputdata[ii]) lowest_lvl = inputdata[ii];
		if (highest_lvl<inputdata[ii]) highest_lvl = inputdata[ii];

		pxAdd(pxBuckets+inputdata[ii],pixels+ii,ii);
	}

	return;
}

void CFlow_Graph::post_initialization()
{
#ifdef GLOBAL_UPDATE	// initialization
	workSinceUpdate = 0;
	assert(n>0 && m>0);
	nm = G_ALPHA*n+m;
#endif
#ifdef LABEL_STATS	// initialization
	minFinalLabel = n;
	maxFinalLabel = 0;
	sumFinalLabel = 0;
#endif

	return;
}
