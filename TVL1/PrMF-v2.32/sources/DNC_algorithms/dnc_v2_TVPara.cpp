/* This file compiles and works correctly for 8-bit and 16-bit ROF model
   using a divide-_n-conquer strategy */


/*! Matlab libraries */
#include <mex.h>
#include <matrix.h>
/*! C/C++ running-time libraries */
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <new>
#include <stdexcept>
#include <vector>
#include <stack>

//#define _DEBUG_TIMING
//#define RESET_NODE_LABELS	// during update reset the labels of all T nodes to 1 (active or inactive)
#define SEP_STATES			// separate the solutions of each state

/*! Type and global headers */
#define DOUBLEFLOW		/* because lambda is of the double type */

/* problem specific inputdata/function configuration */
/* pick one and ONLY one problem type below */
/* Parametric max flow options */
//#define DNC_TV_L2_8bit
//#define DNC_TV_L2_16bit
//#define MATRIX_Lambda		// on-off switch for matrix Lambda input

/*!!!!----- one cannot use more than 16bit with this code because "unsigned short state;" ----!!!!*/

#if !(defined(DNC_TV_L2_1bit) || defined(DNC_TV_L2_8bit) || defined(DNC_TV_L2_16bit))
    #error You must define one of the preprocess macros: DNC_TV_L2_8bit, DNC_TV_L2_16bit. In addition, you can define MATRIX_Lambda if needed.
#endif

#if defined(DNC_TV_L2_1bit)
	/* Definitions for 1-bit logical images */
	#define DNC_TV_L2
	#define N_INPUT_ARG 3
	#define PROBLEM_TYPE "1-bit TV_L2"
	#define GRAYLEVEL 2
	#define MATLABCLASS mxLOGICAL_CLASS
	#define strMATLABCLASS "logical"
	typedef mxLogical imgType;

#elif defined(DNC_TV_L2_8bit)
	/* Definitions for 8-bit greyscale images */
	#define DNC_TV_L2
	#define N_INPUT_ARG 3
	#define PROBLEM_TYPE "8-bit TV_L2"
	#define GRAYLEVEL 256
	#define MATLABCLASS mxUINT8_CLASS
	#define strMATLABCLASS "uint8"
	typedef unsigned char imgType;

#elif defined(DNC_TV_L2_16bit)
	/* Definitions for 16-bit greyscale images */
	#define DNC_TV_L2
	#define N_INPUT_ARG 3
	#define PROBLEM_TYPE "16-bit TV_L2"
	#define GRAYLEVEL 65536
	#define MATLABCLASS mxUINT16_CLASS
	#define strMATLABCLASS "uint16"
	typedef unsigned short imgType;

#else
	#error You must define one of the preprocess macros: DNC_TV_L2_8bit, DNC_TV_L2_16bit. In addition, you can define MATRIX_Lambda if needed.
#endif


/* parametric maximum flow algorithms */
#include"../shared/Common_DNC_MT_types.h"		/* type and structure definitions */
#include "dnc_v2_globalvars.h"	/* global variable and macro definitions */
#include "dnc_v2_TV_L2_core.h"
#include "../shared/Common_DNC_MT.h"

//! Matlab entry function
/*!
	\param nlhs number of left-hand-side output arguments
	\param plhs mxArray of output arguments
	\param nrhs number of right-hand-side input arguments
	\param prhs mxArray of input arguments


*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int _nDim;		// # of dimensions of input image
	const int *_dim;// dimension array of input image

	/* clear global pointers */
	{ out_edges=NULL; neighborType=NULL; nodes=NULL; arcs=NULL; buckets=NULL; stateBuckets=NULL; tree=NULL; pTree=NULL; Lambda=NULL;
	}

	try
	{
	/* print head */
	#ifndef NO_DISP
	printHeader();
	#endif

    #ifdef _DEBUG
         mexPrintf("procInput...\n"); mexEvalString("drawnow;");
    #endif        
	/* input inputdata handling */
	procInput(nlhs,plhs,nrhs,prhs,_nDim,_dim);

	start_time = clock();	/* timing */
	
    #ifdef _DEBUG
         mexPrintf("initLevelBuckets...\n"); mexEvalString("drawnow;");
    #endif        
	/* initialize level buckets */
	initLevelBuckets(inputdata,mm,nn);	/* empty for the TV/L2 problem */

    #ifdef _DEBUG
         mexPrintf("generalInit...\n"); mexEvalString("drawnow;");
    #endif        
	/* initialization counters and etc. */
	generalInit();	/* same for both L1 and L2 problems */

    #ifdef _DEBUG
         mexPrintf("initialize_graph...\n"); mexEvalString("drawnow;");
    #endif        
	/* initialize graph */
	initialize_graph();

	times[0] = clock();	/* timing */

	times[1] = clock();	/* timing */

	#ifdef _DEBUG
		if (lowest_lvl < highest_lvl)
		{	/* non-trivial case */
			printIntermediateSolutions();
		}
	#endif

	/* generate output inputdata structure */
	mtxImage_OUT = mxCreateNumericArray(_nDim,_dim,MATLABCLASS,mxREAL);
	outputdata = (imgType *)mxGetData(mtxImage_OUT);

	if (lowest_lvl == highest_lvl)
	{
        #ifdef _DEBUG
             mexPrintf("trivial...\n"); mexEvalString("drawnow;");
        #endif        
		/* trivial case */
		for (int ii=0; ii<mm*nn*ll; ii++) outputdata[ii] = (imgType)lowest_lvl;
		#ifndef NO_DISP
			mexPrintf(" Output is equal to input");
		#endif
	}
	else
	{
		times[2] = clock();	/* timing */

		#ifdef _DEBUG_TIMING
		updTime = 0; mxfTime = 0;
		#endif

        #ifdef _DEBUG
             mexPrintf("para_max_flow...\n"); mexEvalString("drawnow;");
        #endif        
		/* calculate parametric max-flow */
		para_max_flow();
        #ifdef _DEBUG
             mexPrintf("done...\n"); mexEvalString("drawnow;");
        #endif        

		times[3] = clock();	/* timing */


		#ifdef _DEBUG
		/* check final max-flow and print counters*/
		checkFlow();
		#endif

		/* Statistics */

		#ifndef NO_DISP
		double ticks = 1/(double)CLOCKS_PER_SEC;

	//	mexPrintf(" 1st section time = %f \n 2nd section time = %f \n 3rd section time = %f \n", ((double)times[0])*ticks, ((double)times[1])*ticks, ((double)times[2])*ticks);
	//	mexPrintf(" 1st section time = %f \n 2nd section time = %f \n 3rd section time = %f \n 4th section time = %f \n", ((double)times[0])*ticks, ((double)times[1])*ticks, ((double)times[2])*ticks, ((double)times[3])*ticks);
		mexPrintf(" Graph_init time = %f \n Pixel_bckt time = %f \n Genrl_init section time = %f \n Param_flow time = %f (with %d updates) \n Total time = %f \n", ((double)(times[0]-start_time))*ticks, ((double)(times[1]-times[0]))*ticks, ((double)(times[2]-times[1]))*ticks, ((double)(times[3]-times[2]))*ticks, updateCnt, ((double)(times[3]-start_time))*ticks);
		mexPrintf (" Preflow algorithm statistics:\n");
		mexPrintf (" pushes:      %10ld\n", pushCnt);
		mexPrintf (" relabels:    %10ld\n", relabelCnt);
		mexPrintf (" gaps:        %10ld\n", gapCnt);
		mexPrintf (" gap nodes:   %10ld\n", gNodeCnt);
		mexPrintf ("\n");

		#ifdef _DEBUG_TIMING
		mexPrintf(" Main loop time (the sum of two times below) = %f sec\n",((double)(updTime))*ticks);
		mexPrintf(" Total Preflow time = %f sec\n",((double)(mxfTime))*ticks);
		mexPrintf(" %d updates' time = %f sec\n", updateCnt, ((double)(updTime-mxfTime))*ticks);
		//mexPrintf(" # of the updates generating active nodes: %d\n", actUpdCnt);
		//mexPrintf(" # of total node-update: %d\n", updateCnt*mm*nn);
		#endif
		#endif // #ifndef NO_DISP
	}
	}
	catch(const std::bad_alloc &ba)
	{
		mexPrintf("Run-time ERROR: ");
		mexPrintf(ba.what());
		mexPrintf("Graph has %d nodes and %d arcs, and needs %d bytes of memory.\n", _n, _m,sizeof(node)*(_n+1)+sizeof(bucket)*(_n+1)+sizeof(arc)*(_m+1)+sizeof(cType)*_m);
	}
	catch(const std::exception &ba)
	{
		mexPrintf("Run-time ERROR: ");
		mexPrintf(ba.what());
	}


    #ifdef _DEBUG
         mexPrintf("free_up...\n"); mexEvalString("drawnow;");
    #endif        
	free_up();

	return;
}
