/*
 * das3stepSfunction.c: Basic 'C' template for a level 2 S-function.
 *
 * Copyright 1990-2018 The MathWorks, Inc.
 */


/*
 * You must specify the S_FUNCTION_NAME as the name of your S-function
 * (i.e. replace sfuntmpl_basic with the name of your S-function).
 */

#define S_FUNCTION_NAME  das3stepSfunction
#define S_FUNCTION_LEVEL 2
#define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"
#include "mex.h"
#include "das3mex.h"



/* Error handling
 * --------------
 *
 * You should use the following technique to report errors encountered within
 * an S-function:
 *
 *       ssSetErrorStatus(S,"Error encountered due to ...");
 *       return;
 *
 * Note that the 2nd argument to ssSetErrorStatus must be persistent memory.
 * It cannot be a local variable. For example the following will cause
 * unpredictable errors:
 *
 *      mdlOutputs()
 *      {
 *         char msg[256];         {ILLEGAL: to fix use "static char msg[256];"}
 *         sprintf(msg,"Error due to %s", string);
 *         ssSetErrorStatus(S,msg);
 *         return;
 *      }
 *
 */

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);  /* 0 parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);   /* 0 continuous states */
    ssSetNumDiscStates(S, NSTATES); /* 298 discrete states */

    if (!ssSetNumInputPorts(S, 1)) return; /* one input port */
    ssSetInputPortWidth(S, 0, NMUS); /* 138 inputs are the muscle excitations */
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
    ssSetInputPortComplexSignal(S,0,COMPLEX_NO);
    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     */
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 2)) return; /* two output ports */
    ssSetOutputPortWidth(S, 0, NDOF); /* 11 joint angle outputs */
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S,0,COMPLEX_NO);
    ssSetOutputPortWidth(S, 1, NMUS); /* 138 muscle force outputs */
    ssSetOutputPortDataType(S, 1, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S,1,COMPLEX_NO);

    // initialize sample times
    ssSetNumSampleTimes(S, PORT_BASED_SAMPLE_TIMES); /* port-based sample times */
    //ssSetNumSampleTimes(S, 1); /* block-based sample times */
    ssSetInputPortSampleTime(S, 0, 0);
    ssSetInputPortOffsetTime(S, 0, 0);
    ssSetOutputPortSampleTime(S, 0, 0);
    ssSetOutputPortOffsetTime(S, 0, 0);
    ssSetOutputPortSampleTime(S, 1, 0.003);
    ssSetOutputPortOffsetTime(S, 1, 0);
    
    // initialize work vectors
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    
    // initialize D-work vectors
    ssSetNumDWork(S, 1);
    ssSetDWorkName(S, 0, "x");
    ssSetDWorkWidth(S, 0, NSTATES);
    ssSetDWorkDataType(S, 0, SS_DOUBLE);
    ssSetDWorkComplexSignal(S, 0, COMPLEX_NO);
    ssSetDWorkUsageType(S, 0, SS_DWORK_USED_AS_DSTATE);

    /* Specify the operating point save/restore compliance to be same as a 
     * built-in block */
    ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);

    ssSetOptions(S, 0);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, 0.0);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}



#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
  /* Function: mdlInitializeConditions ========================================
   * Abstract:
   *    In this function, you should initialize the continuous and discrete
   *    states for your S-function block.  The initial states are placed
   *    in the state vector, ssGetContStates(S) or ssGetRealDiscStates(S).
   *    You can also perform any other initialization activities that your
   *    S-function may require. Note, this routine will be called at the
   *    start of simulation and if it is present in an enabled subsystem
   *    configured to reset states, it will be call when the enabled subsystem
   *    restarts execution to reset the states.
   */
  static void mdlInitializeConditions(SimStruct *S)
  {
    real_T *x = (real_T*) ssGetDWork(S,0);  
    real_T *u = ssGetInputPortRealSignal(S,0);
    real_T exF[2], F_SCAP[2][3], F_GH[3], f[NSTATES];
    int_T  i, j, status;
    FILE *fid;
	
	// =============== initialization ====================

	// put zeros in exF (arm support force input)
	exF[0] = 0.0;
	exF[1] = 0.0;

	// load equilibrium.txt to get the initial condition
	// this file is generated by >> das3test('equilibrium')
	mexPrintf("das3steptest: reading initial condition from equilibrium.txt...\n");
	fid = fopen("equilibrium.txt", "r");
	if (fid == NULL) {
		mexPrintf("Error opening equilibrium.txt\n");
		return;
	}
	for (i=0; i<NSTATES; i++) fscanf(fid, "%lf", &x[i]);
	fclose(fid);

    // initialize the das3step function
	mexPrintf("das3steptest: initializing das3step...\n");
	status = das3step(x, u, 0.0, exF, F_GH, F_SCAP, f);
	if (status != 0) {
		mexPrintf("das3steptest: error in initialization of das3step, status=%d\n", status);
		return;
	}
  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#undef MDL_START  /* Change to #define to include function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {

  }
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    int i;
    real_T       *y1 = ssGetOutputPortSignal(S,0); /* set up a pointer for the first output port */
    real_T       *y2 = ssGetOutputPortSignal(S,1); /* set up a pointer for the second output port */
    real_T *x = (real_T*) ssGetDWork(S,0);     /* set up a pointer for the state vector */
    real_T       *u  = ssGetInputPortRealSignal(S,0); /* set up a pointer for the input vector */
    for (i=0; i<NDOF; i++) y1[i] = x[i]; /* joint angles */ 
    for (i=0; i<NMUS; i++) y2[i] = u[i]; /* this is not yet the muscle force */   
}   



#define MDL_UPDATE  /* Change to #undef to remove function */
//#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
        real_T *x = (real_T*) ssGetDWork(S,0);  
        real_T *u = ssGetInputPortRealSignal(S,0);
        real_T exF[2], F_SCAP[2][3], F_GH[3], f[NSTATES];
        int_T  i, j, status;
        //InputRealPtrsType uPtrs    = ssGetInputPortRealSignalPtrs(S,0);
        
        UNUSED_ARG(tid); /* not used in single tasking mode */
        
        //for(i=0; i<NSTATES; i++) x[i] = 3;
        
        status = das3step(x, u, 0.003, exF, F_GH, F_SCAP, f);
        
  }
//#endif /* MDL_UPDATE */



#undef MDL_DERIVATIVES  /* Change to #define to include function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  static void mdlDerivatives(SimStruct *S)
  {
  }
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
}


/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
