/* This file contains the functions that are specifically written for solving the following model:

min TV(u) + lambda * || u - f ||_L2

where

(1) Handles both:
	TV_l2 is discretized at each pixel as sqrt(|u_i+1,j - u_i,j|^2 + |u_i,j+1 - u_i,j|^2)

(2) f and u are both 8/16/p-bit unsigned integer matrix taking values

*/

/* global variable */
imgType *inputdata;				/* pointer of actual input image */
imgType *outputdata;			/* pointer of actual output image */
imgType	*intensities;			/* store the existence of intensity values */
double	*tree;					/* binary tree to store leveled intensity values */

void printIntermediateSolutions();
void stageOne ();
void checkMincut();
void initialize_graph_TV();
void initialize_bTree();
void initialize_fidelity();
inline cType initialize_sv_vt_cap(node*);

inline void liftNodeToN(node* pnode)
{
	pnode->d.dist = n;
}

/* free memory */
void free_up()
{
	if (out_edges)	free(out_edges);
	if (nodes)		free(nodes);
	if (arcs)		free(arcs);
	if (buckets)	free(buckets);
	if (stateBuckets) free(stateBuckets);
	if (tree)		free(tree);
}


#ifdef  SEP_STATES
	inline void update_network_ver1_para()
	{
		node *pNode;
		bucket *buk;
		bucket *ll = buckets+1;

	#ifdef SPECIAL_CNT_NODES
		int nNodesSp = 0;
	#endif

		updateCnt++;

		int oldState,newState;
		cType fChange;
		arc *a;

		for(pNode = nodes; pNode != source; pNode++)
		{	
	#ifdef MATRIX_Lambda
			twolambda = 2*Lambda[pNode-nodes];
	#endif
			if (pNode->d.dist >= _n)
			{	/* this is a node on the S side of the cut, i.e., u_ij >= l_k,
				we are going to increase level. */

				oldState = pNode->state;
				pNode->state = newState = 2*oldState+2;	// take the right child, increase level
				assert(newState <= GRAYLEVEL-1);

				if (tree[newState] != tree[oldState])
				{	/* increase node-t capacity by 2*lambda*( tree[newState] - tree[oldState] ) */

					/* bring the node back to network */
					assert(pNode->d.dist == _n);
					pNode->d.dist = 1;
					pNode->current = pNode->first.pointer;
					outputdata[pNode-nodes] = oldState+1;	// this pixel takes a value at least as large as the one give by "tree[oldState]+0.5"

					fChange = twolambda*( tree[newState] - tree[oldState] );
					assert(fChange>0);

					assert(pNode->excess > -eps);	// must have excess but no node-t capacity
					/* make the increase in node-t capacity or decrease of flow */
					pNode->excess -= fChange;

					buk = stateBuckets+newState;
					if (pNode->excess >= eps)
					{	aStateAdd(buk,pNode);}
					else
					{	iStateAdd(buk,pNode);}

				}

			}
			else
			{	/* this is a node on the T side of the cut, u_ij < l_ij, will decrease level */

				oldState = pNode->state;
				pNode->state = newState = 2*oldState+1;	// take the left child, decrease level
				assert(newState <= GRAYLEVEL-1);

				if (tree[newState] != tree[oldState])
				{	/* decrease node-t capacity by 2*lambda*( tree[oldState] - tree[newState] ) */

					fChange = twolambda*( tree[oldState] - tree[newState] );
					assert(fChange>0);
					assert(pNode->excess < eps);
					pNode->d.dist = 1;
					pNode->current = pNode->first.pointer;	// need to reset current arc

					/* make the decrease */
					pNode->excess += fChange;

					buk = stateBuckets+newState;
					if (pNode->excess >= eps)
					{	
						aStateAdd(buk,pNode);
					}
					else
					{	iStateAdd(buk,pNode);}

				}

			}
		}

		/*---------- done above ----------------*/
		/* block cross-state arcs */
		for(a=arcs; a<arcs+_m; a++)
		{
			if (a->resCap > eps && a->head.pointer->state != a->rev.pointer->head.pointer->state)
			{
				assert(a->rev.pointer->resCap < eps);
				a->resCap = 0;
			}
		}
	}
#else
	inline void update_network_ver1_para()
	{
		node *pNode;
		bucket *buk;
		bucket *ll = buckets+1;

	#ifdef SPECIAL_CNT_NODES
		int nNodesSp = 0;
	#endif

		updateCnt++;

		int oldState,newState;
		cType fChange;
		arc *a;

		/* decrease node-t capacities by lambda */
		for(pNode = nodes; pNode != source; pNode++)
		{	
	#ifdef MATRIX_Lambda
			twolambda = 2*Lambda[pNode-nodes];
	#endif
			if (pNode->d.dist >= _n)
			{	/* this is a node on the S side of the cut, i.e., u_ij >= l_k,
				we are going to increase level. */

				oldState = pNode->state;
				pNode->state = newState = 2*oldState+2;	// take the right child, increase level

				if (tree[newState] != tree[oldState])
				{	/* increase node-t capacity by 2*lambda*( tree[newState] - tree[oldState] ) */

					/* bring the node back to network */
					assert(pNode->d.dist == _n);
					pNode->d.dist = 1;
					pNode->current = pNode->first.pointer;
					outputdata[pNode-nodes] = oldState+1;	// this pixel takes a value at least as large as the one give by "tree[oldState]+0.5"

					fChange = twolambda*( tree[newState] - tree[oldState] );
					assert(fChange>0);

					assert(pNode->excess > -eps);	// must have excess but no node-t capacity
					/* make the increase in node-t capacity or decrease of flow */
					pNode->excess -= fChange;
					if (pNode->excess >= eps)
					{	aAdd(ll,pNode);}
					else
					{	iAdd(ll,pNode);}

				}

			}
			else
			{	/* this is a node on the T side of the cut, u_ij < l_ij, will decrease level */

				oldState = pNode->state;
				pNode->state = newState = 2*oldState+1;	// take the left child, decrease level

				if (tree[newState] != tree[oldState])
				{	/* decrease node-t capacity by 2*lambda*( tree[oldState] - tree[newState] ) */

					fChange = twolambda*( tree[oldState] - tree[newState] );
					assert(fChange>0);
					assert(pNode->excess < eps);
					pNode->current = pNode->first.pointer;	// need to reset current arc

					/* make the decrease */
					pNode->excess += fChange;
					if (pNode->excess >= eps)
					{
						buk = buckets+pNode->d.dist; iDelete(buk,pNode); aAdd(buk,pNode);
					}

					#ifdef RESET_NODE_LABELS
						/* reset the labels of nodes to 1 */
						pNode->current = pNode->first.pointer;
						buk = buckets + pNode->d.dist;
						pNode->d.dist = 1;
						if (pNode->excess >= eps)
						{	aRemove(buk,pNode);	aAdd(ll,pNode);	}
						else
						{	iDelete(buk,pNode);	iAdd(ll,pNode);	}
					#endif
				}

			}
		}

	#ifdef RESET_NODE_LABELS
		aMax = aMin = dMax = 1;
	#endif

		/*---------- done above ----------------*/
		/* block cross-state arcs */
		for(a=arcs; a<arcs+_m; a++)
		{
			if (a->resCap > eps && a->head.pointer->state != a->rev.pointer->head.pointer->state)
			{
				assert(a->rev.pointer->resCap < eps);
				a->resCap = 0;
			}
		}
	}
#endif

inline void update_network_last()
{
	node *pNode;

	/* decrease node-t capacities by lambda */
	for(pNode=nodes;pNode!=source;pNode++)
	{	
		if (pNode->d.dist >= _n)
			outputdata[pNode-nodes] = pNode->state+1;	// this pixel actually takes a value "tree[pNode->state]+0.5"
	}
}

/* add pixels to buckets according to their input intensity */
void initLevelBuckets(imgType *F, int mm, int nn)
{
	return;
}

/* assign value to pixels */
void finish_network_ver1_para()
{
	int ii;
	node *pnode;
	double *tmp_tree = tree-1;
	
	pnode = nodes;

	for (ii=0;ii<_n-2;ii++)
		outputdata[ii] = (imgType)((outputdata[ii]>0)?(tmp_tree[outputdata[ii]]+0.5):lowest_lvl);

	return;
}

/* main parametric max-flow iteraion */
void para_max_flow()
{
	int current_dep;
	#ifdef SEP_STATES
		unsigned int state = 0;
		unsigned int max_state = 0;
		bucket* ll;
	#endif
	lvlCnt = 1;
	for (current_dep=0; current_dep < depth; current_dep++)
	{
		#ifdef _DEBUG_TIMING
		tmpTime = clock();
		#endif

#ifdef SEP_STATES
		max_state = (max_state<<1) + 1;
		if (current_dep == 0)
		{
			stageOne();
			/* clear buckets */
			for (ll = buckets+1; ll<=buckets+dMax; ll++)
			{
				assert(ll->firstActive == sentinelNode);
				ll->firstInactive = sentinelNode;
			}
			aMax = aMin = dMax = 1;
			state=1;
		}
		else
		{
			for (;state<max_state;state++)
			{
				if(stateBuckets[state].firstActive != sentinelNode)
				{
					//if ((state % 2) == 0)
					//{	/* S nodes */
						/* transfer bucket */
						buckets[1].firstActive = stateBuckets[state].firstActive;
						buckets[1].firstInactive = stateBuckets[state].firstInactive;
					//}
					//else
					//{	/* T nodes */
					//	/* transfer bucket with labels */
					//	for (pNode = stateBuckets[state].firstActive; pNode != sentinelNode; pNode = stateBuckets[state].firstActive)
					//	{	ll=stateBuckets+state; aStateRemove(ll,pNode); 
					//		ll=buckets+pNode->d.dist; aAdd(ll,pNode);}
					//	for (pNode = stateBuckets[state].firstInactive; pNode != sentinelNode; pNode = stateBuckets[state].firstInactive)
					//	{	ll=stateBuckets+state; iStateDelete(ll,pNode); 
					//	    ll=buckets+pNode->d.dist; iAdd(ll,pNode); if (pNode->d.dist>dMax) dMax=pNode->d.dist;}
					//}

					/* call algorithm */
					stageOne();
					/* clear buckets */
					for (ll = buckets+1; ll<=buckets+dMax; ll++)
					{
						assert(ll->firstActive == sentinelNode);
						ll->firstInactive = sentinelNode;
					}
					aMax = aMin = dMax = 1;
				}
			}
		}
#else
		/* preflow algorithm */
		stageOne();
#endif

		#ifdef _DEBUG
			mexPrintf("---End of Max-flow level %d---\n",current_dep);
			printIntermediateSolutions();
		#endif
		#ifdef _DEBUG_TIMING
		{
			mexPrintf ("-------------------------------------\n");
			mexPrintf ("c Preflow algorithm statistics:\n");
			mexPrintf ("c pushes:      %10ld\n", pushCnt);
			mexPrintf ("c relabels:    %10ld\n", relabelCnt);
			mexPrintf ("c gaps:        %10ld\n", gapCnt);
			mexPrintf ("c gap nodes:   %10ld\n", gNodeCnt);
			pushCnt=0; relabelCnt=0; gapCnt=0; gNodeCnt=0;
		}
		mexPrintf("Itr #%d Preflow time = %f sec\n",current_dep,(((double)(clock()-tmpTime))*0.001));
		mxfTime += clock()-tmpTime;
		#endif

		#ifdef _DEBUG
		/* check min-cut */
		checkMincut();
		#endif

		/* parametric network update */
		if (current_dep != depth-1)
			update_network_ver1_para();
		else
			update_network_last();

		#ifdef _DEBUG
			mexPrintf("---After network update level %d---\n",current_dep);
			printIntermediateSolutions();
		#endif

		#ifdef _DEBUG_TIMING
		updTime += clock()-tmpTime;
		#endif

		//last_lvl = current_lvl;
		lvlCnt++;
	}

	/* assign intensity values to pixels in the output image */
	finish_network_ver1_para();
}

void gen_tree_nodes(double low, double high, int dep)
{
	double mid;

	if (dep == depth-1)
		tree[pTree[dep]++] = (low+high)*0.5;
	else
	{
		mid = low+ceil((high-low+1)*0.5)-1.0;
		tree[pTree[dep]++] = mid;
		if (low<mid-eps)
			gen_tree_nodes(low,mid-1.0,dep+1);
		else
			gen_tree_nodes(low,mid,dep+1);
		if (mid<high-eps)
			gen_tree_nodes(mid+1.0,high,dep+1);
		else
			gen_tree_nodes(mid,high,dep+1);
	}
}

void initialize_graph()
{
	/* generate the binary tree for storing intensity values \mu needed to be solved for */
    #ifdef _DEBUG
         mexPrintf("initialize_bTree...\n"); mexEvalString("drawnow;");
    #endif
	initialize_bTree();
	if (lowest_lvl == highest_lvl)
	{	return;
	}


	/* finish the TV part of graph */
	/* pointers are stored as offsets */
    #ifdef _DEBUG
         mexPrintf("initialize_graph_TV...\n"); mexEvalString("drawnow;");
    #endif
	initialize_graph_TV();

    #ifdef SEP_STATES
    {
        #ifdef _DEBUG
             mexPrintf("stateBuckets...\n"); mexEvalString("drawnow;");
        #endif
        int pw = 1;
        int ii;
        bucket *l1;
        for (ii=0; ii<depth; ii++) { pw *= 2; }  /* compute 2^depth */
        /* generate state buckets */
		Malloc_n_Check(stateBuckets, bucket, pw-1);
		for (l1 = stateBuckets; l1 < stateBuckets + pw-1; l1++)
		{	 l1 -> firstActive   = sentinelNode; l1 -> firstInactive  = sentinelNode;	}
    }
    #endif

    #ifdef _DEBUG
         mexPrintf("initialize_fidelity...\n"); mexEvalString("drawnow;");
    #endif        
    /* convert offsets to pointers */
	/* contains initialization */
	initialize_fidelity();

	return;
} 

void initialize_bTree()
{
	int ii;
	bucket *l1;
	imgType itensity;
	
	/* calculate lowest_lvl, highest_lvl */
	{
		lowest_lvl=GRAYLEVEL;highest_lvl=0;

		for (ii=0;ii<mm*nn*ll;ii++)
		{
			itensity = inputdata[ii];
			if (lowest_lvl>itensity) {lowest_lvl = itensity;}
			if (highest_lvl<itensity) {highest_lvl = itensity;}
		}
		if (lowest_lvl == highest_lvl)
		{	
			return;
		}
	}

	/* generate a binary tree for storing leveled intensity */
	{
		int nLvl = highest_lvl-lowest_lvl;
		int pw = 1;
		
		depth = (int)ceil(log((double)(nLvl+1))/log((double)2.0));

		/* let pTree[i] points to the first empty slot of ith level tree node */
		Malloc_n_Check(pTree, int, depth);

		for (ii=0; ii<depth; ii++)
		{
			pTree[ii] = pw-1;
			pw *= 2;	// pw = 2^depth
		}

		Malloc_n_Check(tree, double, pw-1)

		/* generate tree nodes */
		gen_tree_nodes((double)(lowest_lvl)+0.5,(double)(highest_lvl)-0.5,0);

		free(pTree); pTree = NULL;
	}

	#ifndef MATRIX_Lambda
		vtcap = (cType)twolambda*((cType)tree[0]);
	#endif
}

/* specific for DNC_TV_L2 */
inline cType initialize_sv_vt_cap(node* pnode)
{
	#ifdef MATRIX_Lambda
		twolambda = 2*Lambda[pnode-nodes];
		if (twolambda < 0) {throw std::runtime_error("Found negative Lambda value");}
		vtcap = (cType)twolambda*((cType)tree[0]);
	#else
		assert(vtcap!=0);
	#endif
	
	return( twolambda*(inputdata[pnode-nodes])-vtcap); 
}

