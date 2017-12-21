/* This file contains the functions that are specifically written for solving the following model:

min TV(u) + lambda * || u - f ||_L1

where

(1) Handles both:
	TV_l1, which is discretized at each pixel as |u_i+1,j-u_i,j|+|u_i,j+1 - u_i,j|

(2) f and u are both 8-bit unsigned integer matrix taking values from 0,1,...,255

*/

/*! Matlab libraries */
#include <mex.h>
#include <matrix.h>
//#include <mclcppclass.h>
/*! C++ time library */
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <string>
#include <vector>
#include <new>
#include <stdexcept>
#include <queue>
#include <stack>

#ifdef _MSC_VER
 #include <windows.h>
 #include <process.h>
#else
// #include <pthread.h>
// use "-lpthread" with gcc
#endif

//#define	INIT_UPDATE /* do initial global update after initial pushes */
//#define GLOBAL_UPDATE	/* global backward breadth-first search w/ frequency */
#ifdef GLOBAL_UPDATE	/* parameter settings */
	#define G_ALPHA 6
	#define G_BETA 12
	#define GLOB_UPDT_FREQ 1
#endif
//#define _EXTRA_TIMING

#define PARAFLOW		/* only applies to a part of code in flow check */
#define DOUBLEFLOW		/* because lambda is of the double type */
//#define LABEL_STATS	/* report distance label statistics */

#if defined(LABEL_STATS) && defined(GLOBAL_UPDATE)
	#error DO NOT enable LABEL_STATS and GLOBAL_UPDATE at the same time
#endif

/* problem specific data/function configuration */
/* pick one and ONLY one problem type below */
/* Parametric max flow options */
//#define MT_TV_L1_8bit
//#define MT_TV_L1_16bit

//#define MATRIX_Lambda

#if defined(MT_TV_L1_1bit)
	#define MT_TV_L1
	#define N_INPUT_ARG 4
	#define PROBLEM_TYPE "1-bit (binary) TV_L1"
	#define GRAYLEVEL 2
	#define MATLABCLASS mxLOGICAL_CLASS
	#define strMATLABCLASS "logical"
	typedef mxLogical imgType;

#elif defined(MT_TV_L1_8bit)
	#define MT_TV_L1
	#define N_INPUT_ARG 4
	#define PROBLEM_TYPE "8-bit TV_L1"
	#define GRAYLEVEL 256
	#define MATLABCLASS mxUINT8_CLASS
	#define strMATLABCLASS "uint8"
	typedef unsigned char imgType;

#elif defined(MT_TV_L1_16bit)
	#define MT_TV_L1
	#define N_INPUT_ARG 4
	#define PROBLEM_TYPE "16-bit TV_L1"
	#define GRAYLEVEL 65536
	#define MATLABCLASS mxUINT16_CLASS
	#define strMATLABCLASS "uint16"
	typedef unsigned short imgType;

#endif

/* type and structure definitions */
#include "../shared/Common_DNC_MT_types.h"
/* global variables and macros */
#include "mtl1_globalvars_v2.h"	/* global variable and macro definitions */
/* parametric maximum flow algorithms */
#include "mtl1_Core_2D_v2.h"
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
	int ii;			// index
	int _nDim;			// # of dimensions of input image
	const int *_dim;	// dimension array of input image

	try
	{
	clock_t mex_start = clock();

	/* clear global pointers */
	{
		/* these will be allocated, so also need to be freed up */
		out_edges=NULL; pxBuckets=NULL; pixels=NULL;
		/* these will be pointed to input data, so do not need to be freed up */
		neighborType=NULL; inputdata=NULL; 
	}

	/* print header */
	#ifndef NO_DISP
	printHeader();
	#endif

	/* check and read input */
	procInput(nlhs, plhs, nrhs, prhs, _nDim, _dim);

	/* initialize level buckets, needed to TV/L1 only */
	try{
		initLevelBuckets();
	}catch(std::bad_alloc &ba)
	{
		mexPrintf("ERROR: %s\n",ba.what());
		free_up();
		return;
	}

	/* generate output data structure */
	mtxImage_OUT = mxCreateNumericArray(_nDim,_dim,MATLABCLASS,mxREAL);
	outputdata = (imgType *)mxGetData(mtxImage_OUT);

	if (lowest_lvl == highest_lvl)
	{	/* trivial case */
		for (ii=0; ii<mm*nn*ll; ii++) outputdata[ii] = (imgType)lowest_lvl;
	}
	else
	{	/* nontrivial case */

		#ifdef _EXTRA_TIMING
		up_updTime = 0; up_mxfTime = 0;
		dn_updTime = 0; dn_mxfTime = 0;
		#endif

		/* calculate parametric max-flow */
		G_current_lvl = lowest_lvl;
		G_rev_current_lvl = highest_lvl+1;

		switch (threadSwitch)
		{
		case 2:	// binary 10
			UpThread(NULL);
			break;
		case 1: // binary 01
			DnThread(NULL);
			break;
		case 3:	// binary 11, 2-threaded

	#ifdef _MSC_VER
			HANDLE hThreads[2];
			unsigned threadID1;
			unsigned threadID2;
			DWORD iObj;

			DnThreadRunMutex = NULL;

			UpThreadState.SetMThreaded();
			DnThreadState.SetMThreaded();

			//InitGraphMutex = CreateMutex( NULL, TRUE, NULL );
			UpFinishInitEvent = CreateEvent( NULL, TRUE, FALSE, NULL );
			DnFinishInitEvent = CreateEvent( NULL, TRUE, FALSE, NULL );
			
			/* Create the second thread.*/
			hThreads[0] = (HANDLE)_beginthreadex( NULL, 0, &UpThread, NULL, 0, &threadID1 );
			/* Create the third thread.*/
			hThreads[1] = (HANDLE)_beginthreadex( NULL, 0, &DnThread, NULL, 0, &threadID2 );

			iObj = WaitForMultipleObjects(2, hThreads, TRUE, INFINITE);	// wait for both threads to terminate
			CloseHandle( hThreads[0] );
			CloseHandle( hThreads[1] );
			//CloseHandle( InitGraphMutex );
			CloseHandle( UpFinishInitEvent );
			CloseHandle( DnFinishInitEvent );
	#else
			throw std::runtime_error("2-threaded version is incomplete under xNIX system.");

			/* below is ok except that I haven't had a similar "event" system work under xNIX system */
			/* therefore, we disallow multithreading under xNIX system */

			//pthread_t threadID1;
			//pthread_t threadID2;

			//pthread_create(&threadID1, NULL, UpThread, NULL);
			//pthread_create(&threadID2, NULL, DnThread, NULL);

			//pthread_join(threadID2, NULL);
			//pthread_join(threadID1, NULL);
	#endif
			break;
		default:
			0;
		}
	}

	/* print final Statistics */
	#ifndef NO_DISP
	printStatistics(mex_start);
	#endif
	}
	catch(const std::exception &ba)
	{
		mexPrintf("Run-time ERROR: ");
		mexPrintf(ba.what());
	}


	free_up();

	return;
}
