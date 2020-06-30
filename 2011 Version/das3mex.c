/*=================================================================
 *
 * das3mex.c
 *
 * Implicit differential equation for DAS3 model : f(x, dx/dt, u) = 0
 
 * This is the source code for the MEX function das3mex.mex32
 * The musculoskeletal model is documented in das3implicit.doc
 * The function specification and usage instructions are in das3mex.m
 *
 * Copyright 2010-2011 Case Western Reserve University
 *
 *=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "das3mex.h"

// comment out the next line to output full Jacobians rather than sparse
// #define SPARSEJACOBIANS

// M_PI is known in gcc but not in Visual C++ 2008
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

// The following is only used in the vector normalization function
#define VERY_SMALL 1e-12

// size of the model (some other size constants such as NDOF are in das3mex.h)
#define NMUS 138					// number of muscle elements 
#define NSTATES 2*NDOF+2*NMUS		// number of system state variables, two for each DOF, two for each muscle
#define MAXMUSDOF 	10              // max number of DOF between a muscle's origin and insertion 
#define MAXPOLTERMS 35              // max number of polynomial terms in muscle geometry 
#define NAMLEN  	40              // length of names 

// define struct that holds muscle properties
// A rough estimate of memory use: MAXPOLTERMS*(MAXMUSDOF+1)*4(bytes)*4(length and 3 GH direction cosines) = 35*11*16 = 7 kB per muscle
// Total: 7kb x 138 = about 1 MB.  Not too bad at all.
typedef struct {
    char   Name[NAMLEN];    				// Name of muscle 
	double Lceopt; 							// Optimal length of CE (m)
    double Pennopt;         				// Pennation angle at optimal length of CE, in radians
	double Width;							// Width of CE force-length relationship relative to Lceopt
	double Fmax;							// Maximal isometric force of CE (N)
	double Vmax;							// Max. contraction velocity in Lceopt/s
	double Tact, Tdeact;					// Activation and deactivation time constants
	double gmax;							// Maximum eccentric force
	double SEEslack;						// Slack length of the SEE (m)
	double PEEslack;						// Slack length of the PEE, relative to Lceopt
	double umax;							// Strain of SEE at Fmax load
	double krel;                            // Stiffness of PEE, force/Fmax at elongation of Width*Lceopt
	int    nmusdof;              			// Number of DOFs between origin and insertion 
    int    musdof[MAXMUSDOF];    			// List of DOFs between origin and insertion 
    int    npolterms;            			// Number of terms in polynomial 
    double  polcoeff[MAXPOLTERMS];        	// Polynomial coefficients 
    int    expon[MAXPOLTERMS][MAXMUSDOF];   // Polynomial exponents 
    int    crossGH;            				// 1 if the muscle cross GH, 0 otherwise 
	
    // the following parameters are only relevant if the muscle crosses GH
    int    nFpolterms[3];            		// number of terms in force vector polynomial 
    double  xFpolcoeff[MAXPOLTERMS];        // polynomial coefficients for the x direction of force 
    double  yFpolcoeff[MAXPOLTERMS];        // polynomial coefficients for the y direction of force 
    double  zFpolcoeff[MAXPOLTERMS];        // polynomial coefficients for the z direction of force 
    int    xFexpon[MAXPOLTERMS][MAXMUSDOF]; // polynomial exponents for the x direction of force 
    int    yFexpon[MAXPOLTERMS][MAXMUSDOF]; // polynomial exponents for the y direction of force 
    int    zFexpon[MAXPOLTERMS][MAXMUSDOF]; // polynomial exponents for the z direction of force 
	
	// the following parameters are derived from other parameters during initialization
	double c1,c2;							// Activation-deactivation rate constants
	double c3;								// Continuity parameter for eccentric force-velocity relationship
	double kPEE;							// Stiffness parameter of PEE, relative to Fmax/Lceopt^2
	double kSEE;							// Stiffness parameter of SEE, in Fmax/m^2	
} muscleprop;

typedef struct {
   char Name[NAMLEN];      // Name of joint (= name of DOF)
   double Lowerlim;        // Lower joint limit
   double Upperlim;        // Upper joint limit 
   double Midpoint;			// halfway between limits
}jointprop;

// some muscle properties that are always the same for each muscle
static double HillA = 0.25;			// Normalized Hill parameter a/Fmax for F-v relationship (usually 0.25)
static double Gmax = 1.5;			// Maximum eccentric muscle force, relative to Fmax

// global variables
static muscleprop musdef;
static muscleprop muscles[NMUS];		// contains all muscle properties
static jointprop joints[NDOF];			// contains all joint properties
static param_struct param;				// contains parameters that must be passed to Autolev code
static int initialized = 0;
static int nonzeromomentarms;			// number of nonzero moment arms (nonzeros in dL/dq), for memory allocation
   
// Passive joint properties
static double JointK1 = 1.0;			// passive joint stiffness within range of motion (linear, Nm/rad) nominal value: 10
static double JointK2 = 5000.0;			// passive joint stiffness outside range of motion (quadratic, Nm/rad^2), nominal value: 5000
static double JointB = 1.0;				// passive joint damping (Nm s/rad)

// Zeros
static double zeros[NSTATES];

// function declarations
int ReadMuscles(char error_msg[]);
double MusclePath(muscleprop *m, double ang[NDOF], double dLm_dang[NDOF], double dLm_dangdang[NDOF][NDOF]);
int find_joint(char *name);
int read_muscle(FILE *biofile, muscleprop *muscle, char error_msg[]);
int read_biofile(char error_msg[]);
void get_force_vector(muscleprop *m, double q[], double nforce_vec[]);
double normalize_vector(double vector[], double norm_vector[]);
double MuscleDynamics(muscleprop *m,										// muscle parameters (input)
	double a, double Lce, double Lcedot, double Lm,							// the input variables
	double *df_da, double *df_dLce, double *df_dLcedot, double *df_dLm, 	// the gradients (output)													
	double *force, double *dforce_dLm, double *dforce_dLce);

// get pointer to a MEX input vector
double* getvec(const mxArray *prhs[], int argnum, int len) {
	mwIndex nrows, ncols;
	nrows = mxGetM(prhs[argnum-1]);
	ncols = mxGetN(prhs[argnum-1]);
	if (!mxIsDouble(prhs[argnum-1]) || mxIsComplex(prhs[argnum-1]) ) {
		printf("Input %d must be double\n", argnum);
		mexErrMsgTxt("das3mex: Incorrect type for input.");
	}
	if ((nrows != len) || (ncols != 1) ) {
		printf("Input %d must be a %d x 1 column vector\n", argnum, len);
		mexErrMsgTxt("das3mex: Incorrect size for input.");
	}
	return mxGetPr(prhs[argnum-1]);
}

// ===================================================================================
// ReadMuscles: read muscle parameters from BIO file
// ===================================================================================
int ReadMuscles(char error_msg[]) {
	
	int i;
	// This function should verify that the number of muscles on the file is equal to NMUS

	// This function will count the total number of nonzero moment arms, over all muscles.
	// This is needed for memory allocation of the sparse Jacobians
	nonzeromomentarms = 0;
	
	if (read_biofile(error_msg))
		return 1;
				
// Preprocessing and error checking
	for (i=0; i<NMUS; i++) {
		muscles[i].kSEE = 1.0/(muscles[i].umax*muscles[i].umax*muscles[i].SEEslack*muscles[i].SEEslack);
		muscles[i].kPEE = muscles[i].krel / (muscles[i].Width * muscles[i].Width);
		muscles[i].c2 = 1.0/muscles[i].Tdeact;
		muscles[i].c1 = 1.0/muscles[i].Tact - muscles[i].c2;
		if (muscles[i].c2 < 0 || muscles[i].c1 < 0) {
			sprintf("Error in muscle %s: time constants must be positive, and deactivation > activation.\n", muscles[i].Name);
			return 1;
		}
		muscles[i].gmax = Gmax;		// for now, gmax is same for all muscles
		muscles[i].c3 = muscles[i].Vmax * HillA * (muscles[i].gmax - 1.0) / (HillA + 1.0);
		
		// calculate how many non-zero moment arms we have in the entire model
		nonzeromomentarms = nonzeromomentarms + muscles[i].nmusdof;
	}

	return 0;	
}

// ===================================================================================
// MusclePath: calculate muscle-tendon length and its derivatives w.r.t. joint angles
// ===================================================================================
double MusclePath(muscleprop *mus, double q[NDOF], double dL_dq[NDOF], double dL_dqdq[NDOF][NDOF]) {

	// returns the length Lm and its gradient dLm/dq for muscle m
	
	int i,j,k,m,kdof,mdof;
	double L, term, dterm;

	// initialize length and all derivatives to 0 
	L = 0.0;
	for (i=0; i<NDOF; i++) {
		dL_dq[i] = 0.0;    
		for (j=0; j<NDOF; j++)		
			dL_dqdq[i][j] = 0.0;       
	}		
	
	// add contributions from each polynomial term
    for(i=0; i < mus->npolterms; i++) {
	
		// add this term's contribution to the muscle length 
		term = mus->polcoeff[i];
		for(j=0; j < mus->nmusdof; j++) {
			mdof = mus->musdof[j];
			for(k=0; k < mus->expon[i][j]; k++)
				term = term * q[mdof];			// this creates polcoeff[i] * product of all q_j to the power expon[i][j] 
		}
		L = L + term;	

		// first derivatives of L with respect to all q's
		for (k=0; k < mus->nmusdof; k++) {
			kdof = mus->musdof[k];
			if ((mus->expon[i][k] > 0) && (q[kdof] != 0.0)) {		// derivative with respect to q_k is zero unless exponent is 1 or higher and q is not zero
				dterm = mus->expon[i][k]*term/q[kdof];
				dL_dq[kdof] = dL_dq[kdof] + dterm;
			
				// second derivatives dL_dqk_dqm for m not equal to k
				for (m=0; m < k; m++) {									// no need for m > k because the dL_dqdq matrix is symmetric
					mdof = mus->musdof[m];
					if ((mus->expon[i][m] > 0) && (q[mdof] != 0.0))
						dL_dqdq[kdof][mdof] = dL_dqdq[kdof][mdof] + mus->expon[i][m]*dterm/q[mdof];
						dL_dqdq[mdof][kdof] = dL_dqdq[kdof][mdof] ;		// symmetry
				}
				
				// second derivatives dL_dqk_dqk
				if (mus->expon[i][k] > 1)
					dL_dqdq[kdof][kdof] = dL_dqdq[kdof][kdof] + (mus->expon[i][k]-1)*dterm/q[kdof];
					
			}
		}
    }
	
	return L;
}

// ===================================================================================
// find_joint: given a joint name, determine the dof number
// returns -1 if name was not found
// and -2 if the name is "end"
// ===================================================================================
int find_joint(char *name) {

    int i;
	char *end = {"end"};
	
	if (strcmp(name, end) == 0) {
	  return -2;
	}
    for (i=0; i<NDOF; i++)
        if (strcmp(name, joints[i].Name) == 0) 
            return i;
    return -1;              
}

// ===================================================================================
// read_muscle: read one muscle from BIO file
// ===================================================================================
int read_muscle(FILE *biofile, muscleprop *muscle, char error_msg[]) {

	char s[25], s2[25], geopart[25];         /* to hold strings read from file */
	char *end = {"end"};
	char *vmrel = {"vmrel"};
	char *krel = {"krel"};
	char *fmax = {"fmax"};
	char *umax = {"umax"};
	char *lceopt = {"lceopt"};
	char *pennopt = {"pennopt"};
	char *lslack = {"lslack"};
	char *peeslack = "PEEslack";
	char *width = {"width"};
	char *dofs = {"dofs"};
	char *parameters = {"parameters"};
	char *time_constants = {"time_constants"};
	char *GHfvectors = {"GH_force_vectors"};
	int i, j, k, dofcount, dof;
	double penn_degrees;

	muscle->crossGH = 0; // this is the default: the muscle does not cross GH. If it does, crossGH will be 1.
 
	do {
		// read muscle property keyword
		fscanf (biofile, "%s", s);
		
		// test for known keywords
		if (strcmp(s,vmrel) == 0)      fscanf(biofile, "%lf", &muscle->Vmax);
		else if (strcmp(s,krel) == 0)       fscanf(biofile,"%lf", &muscle->krel);
		else if (strcmp(s,umax) == 0)       fscanf(biofile,"%lf", &muscle->umax);
		else if (strcmp(s,time_constants) == 0) {
            fscanf(biofile, "%lf", &muscle->Tact);
            fscanf(biofile, "%lf", &muscle->Tdeact);
		}
		else if (strcmp(s,fmax) == 0)		fscanf(biofile, "%lf", &muscle->Fmax);
		else if (strcmp(s,lceopt) == 0)     fscanf(biofile, "%lf", &muscle->Lceopt);
		else if (strcmp(s,lslack) == 0)     fscanf(biofile, "%lf", &muscle->SEEslack);
		else if (strcmp(s,pennopt) == 0) {                 
			fscanf(biofile, "%lf", &penn_degrees); // pennation angle at optimal fiber length is in degrees in the bio file
			muscle->Pennopt = penn_degrees * M_PI / 180;
		}
		else if (strcmp(s,peeslack) == 0)    fscanf(biofile, "%lf", &muscle->PEEslack);
		else if (strcmp(s,width) == 0)       fscanf(biofile, "%lf", &muscle->Width);
		else if (strcmp(s,dofs) == 0) {
            dofcount = -1;
            do {
                fscanf(biofile, "%s", geopart);
				if (strcmp(geopart,"GH_y")==0) muscle->crossGH = 1;
                dof = find_joint(geopart);
				if (dof == -1) {
					sprintf(error_msg,"Bio file error: Unknown joint name: %s",geopart);
					return 1;
				}
				if (dof >= 0) {
					dofcount++;
					muscle->musdof[dofcount] = dof;
				}
			} while (strcmp(geopart,end) != 0);
			
			muscle->nmusdof = dofcount + 1;
			fscanf(biofile, "%s", s2);
			if (strcmp(s2,parameters) != 0) {
				sprintf(error_msg,"Bio file error: Keyword 'parameters' missing in polynomial muscle description.");
				return 1;
			}
			fscanf(biofile, "%d", &muscle->npolterms);
			
			for (i=0; i<muscle->npolterms; i++) {
				for (j=0; j<muscle->nmusdof; j++) {
					// read exponents 
			        fscanf(biofile, "%d", &muscle->expon[i][j]);
				}
				fscanf(biofile, "%lf", &muscle->polcoeff[i]);
			}			
		}
 
		else if (strcmp(s,GHfvectors) == 0) {
		// this muscle crosses GH so we need to read the force vector polynomial terms 
			for (i=0;i<3;i++) { 
				// three polynomials: for x, y, and z components of the force
				fscanf(biofile, "%s", s2);
				if (strcmp(s2,parameters) != 0) {
					sprintf(error_msg,"Bio file error: Keyword 'parameters' missing in polynomial muscle force description.");
					return 1;
				}
				fscanf(biofile, "%d", &muscle->nFpolterms[i]);
				for (k=0; k<muscle->nFpolterms[i]; k++) {
					for (j=0; j<muscle->nmusdof; j++) /* read exponents */ {
						if (i==0)
			           		fscanf(biofile, "%d", &muscle->xFexpon[k][j]);
						if (i==1)
			           		fscanf(biofile, "%d", &muscle->yFexpon[k][j]);
						if (i==2)
			           		fscanf(biofile, "%d", &muscle->zFexpon[k][j]);
					}
					if (i==0)
						fscanf(biofile, "%lf", &muscle->xFpolcoeff[k]);
					if (i==1)
						fscanf(biofile, "%lf", &muscle->yFpolcoeff[k]);
					if (i==2)
						fscanf(biofile, "%lf", &muscle->zFpolcoeff[k]);
				}
			}
		}	  
	} while (strcmp(s,end) != 0);  // stop muscle loop at "end"
 
	return 0;
}

// ===================================================================================
// read_biofile: get joint and muscle information from BIO file
// ===================================================================================

int read_biofile(char error_msg[]) {
	FILE *biofile;
	int njoints, nmuscle, nmuscles;
	double tempvar;
	char keyword[25], jointopt[25], musclename[NAMLEN];

	char *joint = {"joint"};
	char *limits = {"limits"};
	char *muscle = {"muscle"};
	char *def = {"default"};
	char *end = {"end"};
	char *distal_body = {"distal_body"};
	
	// Here we hardcode those model parameters that must be shared with the Autolev code
	// We may instead read them from .bio file in a future version of this code
	
	// scapula-thorax contact force model
	param.epscontact 	= 0.01;		// epsilon for scapula-thorax contact model ~ 0.01 is recommended
	param.kcontact   	= 2e2;		// stiffness k (N/m) in scapula-thorax contact model (2e4 is recommended)
	
	// conoid force model
	// the stiffness of 105N/mm is from Harris et al. 2000 (Am J Sports Med vol. 28, pp. 103-8)
	// and the length of 19.4mm from Harris et al. 2001 (J Shoulder Elbow Surgery, vol. 10, pp. 585-8)
	// was adjusted to 15.9mm so that the force-displacement curve matches the one in Harris et al. 2000
	param.epsconoid		= 0.001;	// epsilon for conoid force-length model, 2 mm is recommended
	param.conK			= 80000;	// stiffness of conoid (105000 N/m is recommended)
	param.conL			= 0.0174;	// length (m) of conoid ligament
	param.conOx =  0.1365; param.conOy =   0.0206; param.conOz = 0.0136;	// coordinates of conoid origin in clavicle frame
	param.conIx = -0.0536; param.conIy =  -0.0009; param.conIz = -0.0266;	// coordinates of conoid insertion in scapula frame

	// morphological parameters
	param.SCx =  0.0014; param.SCy =   0.1334; param.SCz = -0.0563;	// coordinates of SC joint centre in thorax frame
	param.ACx =  0.1909; param.ACy =   0.0077; param.ACz = -0.0026;	// coordinates of AC joint centre in clavicle frame
	param.GHx = -0.0073; param.GHy =  -0.0316; param.GHz = -0.0126;	// coordinates of GH joint centre in scapula frame
	param.ULx =  0.0058; param.ULy =  -0.2926; param.ULz =  0.0049;	// coordinates of humero-ulnar joint centre in humerus frame
	param.RDx =  0.0139; param.RDy =  -0.0196; param.RDz = -0.0124;	// coordinates of ulnar-radial joint centre in ulna frame
	param.HDx =  0.0317; param.HDy =  -0.2344; param.HDz =  0.0103;	// coordinates of radio-carpal joint centre in radius frame
	param.Ax  =  0.1470; param.Ay  =   0.2079; param.Az  =  0.0944;	// radii of the thorax ellipsoid
//	param.TSx = -0.1297; param.TSy =  -0.0086; param.TSz =  0.0313; // coordinates of TS (palpated) in scapula frame
//	param.AIx = -0.1245; param.AIy =  -0.1260; param.AIz =  0.0313;	// coordinates of AI (palpated) in scapula frame
	// actually we want the projections of TS and AI onto the thorax ellipsoid, not the actual TS and AI
	param.TSx = -0.1305; param.TSy =  -0.0181; param.TSz = -0.0046;	// coordinates of TS projection in scapula frame
	param.AIx = -0.1267; param.AIy =  -0.1246; param.AIz =  0.0052;	// coordinates of AI projection in scapula frame
	param.FlexAxisX =  0.999766; param.FlexAxisY =  0.008446; param.FlexAxisZ = -0.019909;	// Elbow flexion axis orientation, in humerus frame
	param.PronAxisX = -0.113295; param.PronAxisY =  0.974059; param.PronAxisZ =  0.195892;	// Forearm pronation axis orientation, in ulna frame
	
	// mass of the six bodies
	param.Mass_Clavicle 	= 0.1560;	
	param.Mass_Scapula 		= 0.7054;	
	param.Mass_Humerus 		= 2.0519;
	param.Mass_Ulna 		= 0.5464;	
	param.Mass_Radius 		= 0.5464;	
	param.Mass_Hand 		= 0.5250;	

	// center of mass of the six bodies
	param.ClavCMx =  0.0983; param.ClavCMy =  0.0167; param.ClavCMz =  0.0042;		
	param.ScapCMx = -0.0530; param.ScapCMy = -0.0229; param.ScapCMz = -0.0247;
	param.HumCMx  =  0.0064; param.HumCMy  = -0.0776; param.HumCMz  = -0.0036;
	param.UlnaCMx = -0.0003; param.UlnaCMy = -0.0634; param.UlnaCMz = -0.0039;
	param.RadCMx  =  0.0106; param.RadCMy  = -0.1179; param.RadCMz  = -0.0129;
	param.HandCMx =  0.0006; param.HandCMy = -0.0905; param.HandCMz = -0.0365;

	// moments of inertia of the six bodies
	param.ClavIxx=0.0000064; param.ClavIyy=0.0000263; param.ClavIzz=0.0000243; param.ClavIxy=-0.0000067; param.ClavIyz=-0.0000095, param.ClavIxz=0.0000029;		
	param.ScapIxx=0.0010000; param.ScapIyy=0.0010000; param.ScapIzz=0.0010000; param.ScapIxy= 0.0; param.ScapIyz= 0.0; param.ScapIxz=0.0;
	param.HumIxx =0.0132000; param.HumIyy =0.0019880; param.HumIzz =0.0132000; param.HumIxy = 0.0; param.HumIyz = 0.0; param.HumIxz =0.0;
	param.UlnaIxx=0.0030585; param.UlnaIyy=0.0004533; param.UlnaIzz=0.0030585; param.UlnaIxy= 0.0; param.UlnaIyz= 0.0; param.UlnaIxz=0.0;
	param.RadIxx =0.0030585; param.RadIyy =0.0004533; param.RadIzz =0.0030585; param.RadIxy = 0.0; param.RadIyz = 0.0; param.RadIxz =0.0;
	param.HandIxx=0.0006387; param.HandIyy=0.0001904; param.HandIzz=0.0006387; param.HandIxy= 0.0; param.HandIyz= 0.0; param.HandIxz=0.0;

	// Now process the .bio file
	if ((biofile = fopen("das3.bio","r")) == NULL)			
	{
		sprintf(error_msg,"Could not open bio file: das3.bio"); 
		return 1;
	}

	njoints = -1;
	nmuscle = -1;

    while (!feof(biofile)) {
	    fscanf(biofile, "%s", keyword);
								
		// read a joint definition
	    if (strcmp(keyword,joint) == 0) {

		     njoints = njoints + 1;
		     fscanf(biofile, "%s", joints[njoints].Name);  // read in name of joint
		     do {
			       fscanf(biofile, "%s", jointopt);
				   if      (strcmp(jointopt,distal_body) == 0) 
						 fscanf(biofile, "%s", jointopt); // distal body, ignore
		           else if (strcmp(jointopt, limits) == 0) {
				         fscanf(biofile, "%lf", &tempvar);
				         joints[njoints].Lowerlim = tempvar * M_PI /180;
				         fscanf(biofile, "%lf", &tempvar);
				         joints[njoints].Upperlim = tempvar * M_PI / 180;
						 joints[njoints].Midpoint = (joints[njoints].Lowerlim + joints[njoints].Upperlim)/2;
		           }
		           else if (strcmp(jointopt,end) == 0); 
		           else {
					    sprintf(error_msg,"Bio file error: Unexpected keyword in joint %s: %s",joints[njoints].Name,jointopt);
						fclose(biofile);
						return 1;
			       }
		     } while (strcmp(jointopt,end) != 0);
	    }
		
		// read a muscle element
		else if (strcmp(keyword,muscle) == 0) {
		    fscanf(biofile, "%s", musclename);  // get muscle name 
		    if (strcmp(musclename,def) == 0) {
				// store this in the "default" muscle
				sprintf(musdef.Name,musclename);
			    if (read_muscle(biofile, &musdef, error_msg))  {        
					fclose(biofile);
					return 1;
				}	
			}
		    else {
			    nmuscle = nmuscle + 1;
			    muscles[nmuscle] = musdef;         // set all properties to default initially 
			    sprintf(muscles[nmuscle].Name, musclename);   // store muscle name
			    if (read_muscle(biofile, &muscles[nmuscle], error_msg)) {
					fclose(biofile);
					return 1;	
				}					
		    }
		}
		
		else 
			fscanf(biofile, "\n");
	}
	 
	fclose(biofile);
	
	nmuscles = nmuscle + 1;
	// make sure nmuscles = NMUS	 
	if (nmuscles!=NMUS) {
		sprintf(error_msg,"Error in bio file: Number of muscles = %i, while NMUS = %i",nmuscles,NMUS);
		return 1;
	}

	return 0;
}

// ==========================================================================================
// get_force_vector: calculates the muscle's GH force vector at the current state
// 						this is how much 1 N of muscle force would contribute to the GH force
// ==========================================================================================
void get_force_vector(muscleprop *m, double q[], double nforce_vec[]) {

	int i,j,k;
	double term, qnonzero;
	double force_vec[3];

	for (i=0; i<3 ; i++) force_vec[i] = 0.0;            // initialize the force vector to 0

	// x component
	for(i=0; i < m->nFpolterms[0]; i++) {
		// add this term's contribution to the force component
		term = m->xFpolcoeff[i];
		for(j=0; j < m->nmusdof; j++) {
			qnonzero = q[m->musdof[j]];
			for(k=0; k < m->xFexpon[i][j]; k++)
				term = term * qnonzero;	// this creates q to the power expon[][]
		}
		force_vec[0] = force_vec[0] + term;
	}
	
	// y component
	for(i=0; i < m->nFpolterms[1]; i++) {
		// add this term's contribution to the force component
		term = m->yFpolcoeff[i];
		for(j=0; j < m->nmusdof; j++) {
			qnonzero = q[m->musdof[j]];
			for(k=0; k < m->yFexpon[i][j]; k++)
				term = term * qnonzero;	// this creates q to the power expon[][]
		}
		force_vec[1] = force_vec[1] + term;
	}

	// z component	
	for(i=0; i < m->nFpolterms[2]; i++) {
		// add this term's contribution to the force component
		term = m->zFpolcoeff[i];
		for(j=0; j < m->nmusdof; j++) {
			qnonzero = q[m->musdof[j]];
			for(k=0; k < m->zFexpon[i][j]; k++)
				term = term * qnonzero;	// this creates q to the power expon[][]
		}
		force_vec[2] = force_vec[2] + term;
	}
	normalize_vector(force_vec, nforce_vec);
	return;
}

// ===================================================================================
//	normalize_vector: normalizes a 1x3 vector
//	(from Dynamics Pipeline)
// ===================================================================================

double normalize_vector(double vector[], double norm_vector[]) {
   double magnitude;
   magnitude = sqrt((vector[0]*vector[0])+(vector[1]*vector[1])+(vector[2]*vector[2]));

   if (magnitude < VERY_SMALL)
   {
      norm_vector[0] = vector[0];
      norm_vector[1] = vector[1];
      norm_vector[2] = vector[2];
   }
   else
   {
      norm_vector[0] = vector[0]/magnitude;
      norm_vector[1] = vector[1]/magnitude;
      norm_vector[2] = vector[2]/magnitude;
   }

   return (magnitude);
}

// ===================================================================================================
// MuscleDynamics: the implicit muscle dynamics f()=0, returns f(a,s,sdot,Lm) and its derivatives
// For pennated muscle model, the state variable s is Lce*cos(p)
// ===================================================================================================
double MuscleDynamics(muscleprop *m,									// muscle parameters (input)
	double a, double s, double sdot, double Lm,							// the input variables
	double *df_da, double *df_ds, double *df_dsdot, double *df_dLm, 	// the gradients (output)													
	double *force, double *dforce_dLm, double *dforce_ds) {				// muscle force (output) and derivatives

	double Lce, cosp, dLce_ds, dcosp_ds, b;
	double Lcedot, dLcedot_dsdot, dLcedot_ds;
	double F1, dF1_dLce, dF1_ds;
	double f,x,k1;
	double F2, dF2_dLcedot, dF2_dsdot, dF2_ds;
	double F3, dF3_dLce, dF3_ds;
	double F4, dF4_ds, dF4_dLm;
	double F5, dF5_dsdot;
		
	// Compute Lce and cos(p) from state s using the constant volume constraint: Lce*sin(p) = Lceopt*sin(popt)
	b = sin(m->Pennopt);				// corrected from Lceopt * sin(popt)
	Lce = sqrt(s*s + b*b);
	cosp = s/Lce;
	dLce_ds = cosp;
	dcosp_ds = b*b/Lce/Lce/Lce;
	
	// Compute Lcedot and its derivatives with respect to sdot and s
	Lcedot = sdot*cosp;
	dLcedot_dsdot = cosp;
	dLcedot_ds = sdot*dcosp_ds;
		
	// F1 is the normalized isometric force-length relationship at max activation
	x = (Lce - 1.0)/m->Width;
	F1 = exp(-x*x);
	dF1_dLce = -2.0*x*F1 / m->Width;
	dF1_ds = dF1_dLce * dLce_ds;
		
	// F2 is the normalized force-velocity relationship
	if (Lcedot < 0) {
		// concentric contraction
		x = m->Vmax - Lcedot/HillA;
		F2 = (m->Vmax + Lcedot)/x;
		dF2_dLcedot = (1.0 + F2/HillA)/x;
	}
	else {
		// eccentric contraction
		x = Lcedot + m->c3;
		F2 = (m->gmax*Lcedot + m->c3) / x;
		dF2_dLcedot = (m->gmax - F2)/x;
	}
	dF2_dsdot =  dF2_dLcedot * dLcedot_dsdot;
	dF2_ds = dF2_dLcedot * dLcedot_ds;
	
	// F3 is the PEE force-length relationship
	k1 = 10.0/m->Fmax*m->Lceopt;	// stiffness of the linear term is 10 N/m, convert to Fmax/Lceopt units	
	x = (Lce - m->PEEslack);		// elongation of PEE, relative to Lceopt
	F3 = k1*x;						// low stiffness linear term
	dF3_dLce = k1;
	if (x>0) {						// add quadratic term for positive elongation						
		F3 = F3 + m->kPEE*x*x;
		dF3_dLce = dF3_dLce + 2*m->kPEE*x;
	}
	dF3_ds = dF3_dLce * dLce_ds;
	
	//  F4 is the SEE force-length relationship
	k1 = 10.0/m->Fmax;			// stiffness of the linear term is 10 N/m, convert to Fmax/meter	
	x = Lm - s * m->Lceopt - m->SEEslack;			// elongation of SEE, in meters
	F4 = k1*x;										// low stiffness linear term
	dF4_ds = -k1*m->Lceopt;
	dF4_dLm = k1;
	if (x>0) {										// add quadratic term for positive deformation
		F4 = F4 + m->kSEE*x*x;
		dF4_ds = dF4_ds - 2 * m->kSEE * m->Lceopt * x;
		dF4_dLm = dF4_dLm + 2 * m->kSEE * x;
	}

	// F5 is viscous damping in the projected CE (0.001 of Fmax at Lceopt/s) to ensure df/dLcedot is never zero
	// this is only really needed if we want to solve Lcedot explicitly for all possible muscle states (including a=0)
	F5 = .001*sdot;
	dF5_dsdot = .001;
		
	// Compute f, the force imbalance in the muscle contraction, and its derivatives
	f = F4 - (a*F1*F2 + F3)*cosp - F5;
	*df_da = -F1*F2*cosp;
	*df_ds = dF4_ds - (a*(dF1_ds*F2 + F1*dF2_ds) + dF3_ds)*cosp - (a*F1*F2 + F3)*dcosp_ds;
	*df_dsdot = -a*F1*dF2_dsdot*cosp - dF5_dsdot;
	*df_dLm = dF4_dLm;
	
	// Muscle force is the force in SEE
	*force = m->Fmax*F4;
	*dforce_dLm = m->Fmax * dF4_dLm;
	*dforce_ds = m->Fmax * dF4_ds;
	
	// Return the imbalance
	return f;
}

// =========================================================================
// mexFunction: this is the actual MEX function interface
// =========================================================================
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// working variables
	int i, j, k;
	double ang, angvel, d;
	char *command;
	int do_dynamics = 0;
	int do_derivatives = 0;			// do we need to output derivatives?
	int do_ghforce = 0;				// do we need to output GH force?
	int do_LCEopt = 0;
	int do_SEEslack = 0;
	int do_stick = 0;
	int do_scapulacontact = 0;
	int do_moments = 0;
	int do_fscapula = 0;
	int do_GHF_only = 0;
	int do_externalangles = 0;
	int do_muscleforces = 0;
	int do_momentarms = 0;
	int do_lengths = 0;
	int do_name = 0;
	int do_crossGH = 0;
	double GHforce_vec[3];
	static char error_msg[100];	

	// internal muscle variables
	double Lm[NMUS];						// muscle+tendon length, based on skeleton kinematic state
	double dLm_dang[NMUS][NDOF];			// derivatives of Lm with respect to joint angles
	double dLm_dangdang[NMUS][NDOF][NDOF]; 	// second derivative of Lm w.r.t. joint angles
	double force[NMUS];						// muscle forces
	double g[NMUS];							// muscle force imbalance
	double dg_da[NMUS], dg_dLce[NMUS], dg_dLcedot[NMUS], dg_dLm[NMUS]; 	// derivatives of muscle imbalance
	double h[NMUS];							// activation dynamics imbalance
	double dh_da[NMUS], dh_dadot[NMUS], dh_du[NMUS];
	double dforce_dLm[NMUS], dforce_dLce[NMUS];
		
	// internal multibody dynamics variables, used to communicate with Autolev code
	double *q, *qd, *qdd;
	double mom[NDOF];				// joint moments
	double zero[NDOF];				// dynamics imbalance
	double dz_dq[NDOF][NDOF];		// and its derivatives
	double dz_dqd[NDOF][NDOF];
	double dz_dqdd[NDOF][NDOF];
	double dz_dmom[NDOF][NDOF];
	double dmom_dang[NDOF][NDOF];
	double dmom_dangvel[NDOF];		// is a diagonal matrix (no coupling between joints), so no need for full storage
	double dmom_dLce[NDOF][NMUS];
	double F_GH[3];
	double F_SCAP[2][3];
	double Stick[NSTICK][3];
	double qTH[3], mTH[3], exF[2];

	// MEX function pointers to inputs and outputs
	mxDouble *x, *xdot, *u, *M, *forF;
    mwSize nrows = NSTATES;		// if we just use NSTATES, multiplication by NSTATES does not work for full jacobian indexing!
	mxDouble *f, *FGHtotal, *stickdata;
	mxDouble *df_dx, *df_dxdot, *df_du, *MA;
	mwIndex *df_dx_irs, *df_dx_jcs;
	mwIndex *df_dxdot_irs,*df_dxdot_jcs;
	mwIndex *df_du_irs, *df_du_jcs;
	mwIndex *MA_irs, *MA_jcs;
	mwIndex Ndf_dx, Ndf_dxdot, Ndf_du, NMA;
	mxDouble *mforces;				// muscle forces
	mxDouble *moments;				// joint moments
	mxDouble *LCEopt;
	mxDouble *SEEslack;
	
	// If there is no input argument, we should initialize the model
	if (nrhs == 0) {
		if (nlhs != 0)
			mexErrMsgTxt("das3mex: No output argument allowed for initialization.");
		printf("*****************************************************\n");
		printf("*                      DAS3MEX                      *\n");  
		printf("*  (c) 2010-2012 Case Western Reserve University    *\n");
		printf("*****************************************************\n");
		printf("Initializing model...\n");
		if (ReadMuscles(error_msg)) {
			mexErrMsgTxt(error_msg);
			return;
		}
		for (i=0; i<NSTATES; i++) zeros[i] = 0.0;
		initialized = 1959;
		return;
	}
	
	// Give error message if model was not initialized
	if (initialized != 1959)
		mexErrMsgTxt("das3mex: initialization must be done first.");
		
	// If first argument is a string, see what the user wants.
	// If it is not a string, user wants dynamics
    if (mxIsChar(prhs[0])) {
		do_dynamics = 0;
		command = mxArrayToString(prhs[0]);
		if (strcmp(command,"LCEopt") == 0)
			do_LCEopt = 1;
		else if (strcmp(command,"SEEslack") == 0)
			do_SEEslack = 1;
		else if (strcmp(command,"Stick") == 0)
			do_stick = 1;
		else if (strcmp(command,"Scapulacontact") == 0)
			do_scapulacontact = 1;
		else if (strcmp(command,"Jointmoments") == 0)
			do_moments = 1;
		else if (strcmp(command,"Muscleforces") == 0)
			do_muscleforces = 1;
		else if (strcmp(command,"Momentarms") == 0)
			do_momentarms = 1;
		else if (strcmp(command,"Musclelengths") == 0)
			do_lengths = 1;
		else if (strcmp(command,"Musclename") == 0)
			do_name = 1;
		else if (strcmp(command,"crossGH") == 0)
			do_crossGH = 1;
		else {
			printf("Command given was: %s\n", command);
			mexErrMsgTxt("das3mex: unknown command.");
		}
	}
	else {
		do_dynamics = 1;
	}

	// Get the inputs we need, and some commands can be completed already
	
	if (do_LCEopt) {
		if (nrhs != 1) {
			mexErrMsgTxt("No inputs must be given for das3mex command LCEopt.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("das3mex command LCEopt only has one output.");
		}
		// Create NMUS x 1 matrix for the LCEopt output
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		LCEopt = mxGetPr(plhs[0]);
		// Copy LCEopt from the muscle structs
		for (i=0; i<NMUS; i++) *LCEopt++ = muscles[i].Lceopt;
		return;									// and we are done!
	}
	else if (do_SEEslack) {
		if (nrhs != 1) {
			mexErrMsgTxt("No inputs must be given for das3mex command SEEslack.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("das3mex command SEEslack only has one output.");
		}
		// Create NMUS x 1 matrix for the SEEslack output
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		SEEslack = mxGetPr(plhs[0]);
		// Copy SEEslack from the muscle structs
		for (i=0; i<NMUS; i++) *SEEslack++ = muscles[i].SEEslack;
		return;									// and we are done!
	}
	else if (do_name) {
		if (nrhs != 2) {
			mexErrMsgTxt("One input (muscle number) must be given for das3mex command 'Musclename'.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("das3mex command 'Musclename' only has one output.");
		}
		i = (int) *mxGetPr(prhs[1]);
		if (i<1 || i>NMUS) {
			mexErrMsgTxt("Muscle number out of range in das3mex command 'Musclename'.");
		}				
		plhs[0] = mxCreateString(muscles[i-1].Name);
		return;
	}
	else if (do_crossGH) {
		if (nrhs != 2) {
			mexErrMsgTxt("One input (muscle number) must be given for das3mex command 'crossGH'.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("das3mex command 'crossGH' only has one output.");
		}
		i = (int) *mxGetPr(prhs[1]);
		if (i<1 || i>NMUS) {
			mexErrMsgTxt("Muscle number out of range in das3mex command 'crossGH'.");
		}				
		plhs[0] = mxCreateDoubleScalar(muscles[i-1].crossGH);
		return;
	}
	else if (do_dynamics) {
		if (nrhs < 3) {
			mexErrMsgTxt("Three inputs (x,xdot,u) are needed for das3mex dynamics.");
		}
		if (nrhs > 5) {
			mexErrMsgTxt("No more than five inputs (x,xdot,u,M,exF) are allowed for das3mex dynamics.");
		}
		// determine which optional outputs were requested with dynamics
		do_derivatives = 0;
		do_ghforce = 0;
		do_fscapula = 0;
		do_externalangles = 0;
		if (nlhs == 4) {
			do_derivatives = 1;
		}
		else if (nlhs == 5) {
			do_derivatives = 1;
			do_ghforce = 1;
		}
		else if (nlhs == 6) {
			do_derivatives = 1;
			do_ghforce = 1;
			do_fscapula = 1;
		}
		else if (nlhs == 7) {
			do_derivatives = 1;
			do_ghforce = 1;
			do_fscapula = 1;
			do_externalangles = 1;
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One, four, five, six, or seven outputs must be given for das3mex dynamics.");
		}
		x 		= getvec(prhs, 1, NSTATES);
		xdot 	= getvec(prhs, 2, NSTATES);
		u 		= getvec(prhs, 3, NMUS);
	}
	else {					// all other commands have one input (state x) and one output
		if (nrhs != 2) {
			printf("Command %s needs exactly one input (x).\n", command);
			mexErrMsgTxt("Incorrect input for das3mex.");
		}
		if (nlhs > 1) {
			printf("Command %s has only one output.\n", command);
			mexErrMsgTxt("Incorrect output request for das3mex.");
		}
		x 		= getvec(prhs, 2,NSTATES);
		xdot	= zeros;
		u 		= zeros;
	}
		
	// Compute the muscle dynamics, and get muscle forces
	for(i=0; i<NMUS; i++) {	
		// Calculate imbalance of activation dynamics
		double rate = muscles[i].c1*u[i] + muscles[i].c2;
		h[i] = xdot[2*NDOF+NMUS+i] - rate*(u[i] - x[2*NDOF+NMUS+i]);
		if (do_derivatives) {
			dh_da[i] = rate;
			dh_dadot[i] = 1.0;
			dh_du[i] = -rate - muscles[i].c1*(u[i] - x[2*NDOF+NMUS+i]);
		}
			
		// Calculate muscle length Lm and derivatives dLm/dq from generalized coordinates in x
		Lm[i] = MusclePath(&muscles[i], x, dLm_dang[i], dLm_dangdang[i]);

		// Calculate muscle force imbalance, normalized to Fmax
		g[i] = MuscleDynamics(&muscles[i],
			x[2*NDOF+NMUS+i],		// active state of muscle i
			x[2*NDOF+i],			// Lce of muscle i
			xdot[2*NDOF+i],			// Lcedot
			Lm[i],					// muscle length
			&dg_da[i], &dg_dLce[i], &dg_dLcedot[i], &dg_dLm[i],
			&force[i], &dforce_dLm[i], &dforce_dLce[i]);				
	}

	// Does user only want muscle forces?
	if (do_muscleforces) {
		// Create NMUS x 1 matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);
		// Copy muscle forces from the force[] array that we already have
		for (i=0; i<NMUS; i++) {
			*mforces++ = force[i];
		}
		return;				// done!
	}

	// Does user only want muscle lengths?
	if (do_lengths) {
		// Create NMUS x 1 matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);
		// Copy lengths from the Lm[] array that we already have
		for (i=0; i<NMUS; i++) {
			*mforces++ = Lm[i];
		}
		return;				// done!
	}

	// Does user only want moment arms?
	if (do_momentarms) {
		// Create NMUS x NDOF sparse matrix for result
		plhs[0] = mxCreateSparse(NMUS, NDOF, nonzeromomentarms, 0);
		MA = mxGetPr(plhs[0]);
		MA_irs = mxGetIr(plhs[0]);
		MA_jcs = mxGetJc(plhs[0]);
		NMA = 0;
		
		// copy moment arms from dLm_dang and store them in sparse NMUS x NDOF matrix
		// must be done column-wise
		for (j=0; j<NDOF; j++) {			
			MA_jcs[j] = NMA;					// store element number where this column starts
			for (i=0; i<NMUS; i++) {
				if (dLm_dang[i][j] != 0) {
					MA_irs[NMA] = i;				// store row number of this matrix element
					MA[NMA] = -dLm_dang[i][j];		// store the value of this matrix element
					NMA++;		
				}
			}
		}
		MA_jcs[NDOF] = NMA;		// store final element number for the matrix
		return;					// done!
	}
	
	// Compute the joint moments
	for (i=0; i<NDOF; i++) {
		// initialize derivatives to zero
		if (do_derivatives) {
			for (j=0; j<NDOF; j++) dmom_dang[i][j] = 0.0;
			for (j=0; j<NMUS; j++) dmom_dLce[i][j] = 0.0;
			dmom_dangvel[i] = 0.0;
		}
		
		// start with passive joint moment
		ang = x[i];									// joint angle i is state variable i
		angvel = x[NDOF+i];							// and the corresponding angular velocity
		mom[i] = -JointK1*(ang-joints[i].Midpoint) - JointB*angvel;		// a small linear elastic moment at all angles, plus damping
		if (do_derivatives) {
			dmom_dang[i][i] = -JointK1;				// and its derivatives
			dmom_dangvel[i] = -JointB;
		}
		d = joints[i].Lowerlim - ang;				// are we below min angle?
		if (d > 0) {								// yes, add quadratic term
			mom[i] = mom[i] + JointK2*d*d;	
			if (do_derivatives) {						// and its derivative
				dmom_dang[i][i] = dmom_dang[i][i] - 2*JointK2*d;
			}
		}
		d = ang - joints[i].Upperlim;
		if (d > 0) {								// are we above max angle?
			mom[i] = mom[i] - JointK2*d*d;			// yes, add quadratic term
			if (do_derivatives) {
				dmom_dang[i][i] = dmom_dang[i][i] - 2*JointK2*d;
			}
		}
		
		// add the muscle moments
		for (j=0; j<NMUS; j++) {
			mom[i] = mom[i] - dLm_dang[j][i]*force[j];		// moment arm is -dLm/dang
			if (do_derivatives) {
				for (k=0; k<NDOF; k++) {
					dmom_dang[i][k] = dmom_dang[i][k] - dLm_dang[j][i]*dforce_dLm[j]*dLm_dang[j][k] - dLm_dangdang[j][i][k]*force[j];
				}
				dmom_dLce[i][j] = dmom_dLce[i][j] - dLm_dang[j][i]*dforce_dLce[j];
			}
		}
	}
	
	// Does user only want joint moments?
	if (do_moments) {
		// Create NDOF x 1 matrix for the function output
		plhs[0] = mxCreateDoubleMatrix(NDOF, 1, mxREAL);
		moments = mxGetPr(plhs[0]);
		// Copy moments data from the mom[] array that we already have
		for (i=0; i<NDOF; i++) {
			*moments++ = mom[i];
		}
		return;				// done!
	}
	
	// did the user supply the 5 extra moments as 4th input of MEX function?
	if (nrhs > 3) {
		M = getvec(prhs, 4, 5);
		// put the first three moments into mTH:
		mTH[0] = M[0];
		mTH[1] = M[1];
		mTH[2] = M[2];
		// and add the last two moments to the last two elements of mom (flexion and pronation moments):
		mom[NDOF-2] = mom[NDOF-2] + M[3];
		mom[NDOF-1] = mom[NDOF-1] + M[4];
	}
	else {
		mTH[0] = 0;
		mTH[1] = 0;
		mTH[2] = 0;	
	}

	// did the user supply the mobile arm support force as 5th input of MEX function?
	if (nrhs > 4) {
		forF = getvec(prhs, 5, 2);
		exF[0] = forF[0];
		exF[1] = forF[1];
	}
	else {
		exF[0] = 0;
		exF[1] = 0;
	}
	
	// Call the C function that was generated by Autolev and cleaned by autolevclean.exe
	q = &x[0];
	qd = &x[NDOF];	
	qdd = &xdot[NDOF];
	das3_al(&param, q, qd, qdd, mom, mTH, exF, zero, dz_dq, dz_dqd, dz_dqdd, dz_dmom, F_GH, F_SCAP, Stick, qTH);
	
	// for (i=0; i<NDOF; i++) {
		// for (j=0; j<NDOF; j++) printf("%8.3f ",dz_dmom[i][j]);
		// printf("\n");
	// }
	
	// Does user only want stick figure data?
	if (do_stick) {
		// Create 3 x NSTICK matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(3, NSTICK, mxREAL);
		stickdata = mxGetPr(plhs[0]);

		// Copy stickdata from the matrix that came out of Autolev
		for (i=0; i<NSTICK; i++) {
			// fill column i
			for (j=0; j<3; j++) {
				*stickdata++ = Stick[i][j];
			}
		}
		return;			// done!
	}
	
	// Does user only want scapula contact?
	if (do_scapulacontact) {
		// Create 2 x 1 matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);
		// AI and TS positions are in Stick calculated by Autolev
		// TS is the 4th bony landmark
		// AI is the 5th
		for (i=3; i<5; i++) {
			*mforces++ = (Stick[i][0]/param.Ax)*(Stick[i][0]/param.Ax)+(Stick[i][1]/param.Ay)*(Stick[i][1]/param.Ay)+(Stick[i][2]/param.Az)*(Stick[i][2]/param.Az)-1;
		}
		return;				// done!
	}
	
	// Assemble the MEX function output for GH force
	if (do_ghforce) {
		// Create 3x1 matrix for the 5th return argument, the total GH reaction force
		plhs[4] = mxCreateDoubleMatrix(3, 1, mxREAL);
		FGHtotal = mxGetPr(plhs[4]);

		// first get the net force from Autolev
		for (j=0; j<3; j++) {
			FGHtotal[j] = F_GH[j];
		}
		
		// add muscle contributions
		for (i=0; i<NMUS; i++) {
			if (muscles[i].crossGH) {  
				// only include muscles that cross GH
				get_force_vector(&muscles[i], q, GHforce_vec);
				for (j=0; j<3; j++) {
					FGHtotal[j] += force[i] * GHforce_vec[j];
				}
			}
		}
	}
	
	// Does user want scapula-thorax contact forces?
	if (do_fscapula) {
		// Create 3x2 matrix for the 6th return argument, the scapula-thorax contact forces
		plhs[5] = mxCreateDoubleMatrix(3, 2, mxREAL);
		f = mxGetPr(plhs[5]);
		for (i=0; i<2; i++) {
			for (j=0; j<3; j++) {
				*f++ = F_SCAP[i][j];
			}
		}
	}
	
	// Does user want thorax-humerus angles?
	if (do_externalangles) {
		// Create 3x1 matrix for the 7th return argument, the thorax-humerus angles
		plhs[6] = mxCreateDoubleMatrix(3, 1, mxREAL);
		f = mxGetPr(plhs[6]);
		for (i=0; i<3; i++) {
			*f++ = qTH[i];
		}
	}

	// Create matrix for the 1st return argument (f)
	plhs[0] = mxCreateDoubleMatrix(NSTATES, 1, mxREAL);
	f = mxGetPr(plhs[0]);
		
	// the first NDOF rows of the implicit differential equation (IDE) are: qdot-dq/dt = 0
	for (i=0; i<NDOF; i++) f[i] = x[NDOF+i] - xdot[i];
	// the next NDOF rows of the IDE are the equations of motion from Autolev
	for (i=0; i<NDOF; i++) f[NDOF+i] = zero[i];
	// the next NMUS rows of the IDE are the muscle contraction dynamics
	for (i=0; i<NMUS; i++) f[2*NDOF+i] = g[i];
	// the final NMUS rows of the IDE are the muscle activation dynamics: da/dt - (u-a)(c1*u + c2) = 0
	for (i=0; i<NMUS; i++) f[2*NDOF+NMUS+i] =  h[i];
	
	
	// Create matrix for the Jacobian results
	if (do_derivatives) {	
        #ifdef SPARSEJACOBIANS
            // The sparse Jacobians have to be filled in column order, using Matlab sparse data structure
            // --------Second output of MEX function: Jacobian df/dx
            plhs[1] = mxCreateSparse(NSTATES,NSTATES, 
                NDOF*(1 + 2*NDOF + NMUS) + 3*NMUS + nonzeromomentarms, 0);
            // TODO: Ton, is it possible to do the memory allocations (mxCreate...) only during initialization,
            // store the pointer as a global variable and keep using the same pointer.  This should speed things up.
            df_dx = mxGetPr(plhs[1]);
            df_dx_irs = mxGetIr(plhs[1]);
            df_dx_jcs = mxGetJc(plhs[1]);
            Ndf_dx = 0;

            // derivatives with respect to q, columns 1..NDOF of x
            for (i=0; i<NDOF; i++) {			
                df_dx_jcs[i] = Ndf_dx;				// store element number where this column starts

                // derivatives of ZERO with respect to q are in rows NDOF+1 to 2*NDOF
                for (j=0; j<NDOF; j++) {
                    df_dx_irs[Ndf_dx] = NDOF+j;		// store row number of this matrix element
                    df_dx[Ndf_dx] = dz_dq[j][i];	// store the value of this matrix element
                    // add the contributions dz/dmom * dmom/dq
                    for (k=0; k<NDOF; k++) {
                        df_dx[Ndf_dx] = df_dx[Ndf_dx] + dz_dmom[j][k]*dmom_dang[k][i];
                    }
                    Ndf_dx++;
                }

                // derivatives of muscle imbalance with respect to q are in rows 2*NDOF+1 to 2*NDOF+NMUS
                for (j=0; j<NMUS; j++) {
                    if (dLm_dang[j][i] != 0) {						// this only exists if muscle j crosses joint i
                        df_dx_irs[Ndf_dx] = 2*NDOF+j;				// store row number of this matrix element
                        df_dx[Ndf_dx] = dg_dLm[j]*dLm_dang[j][i];	// store the value of this matrix element
                        Ndf_dx++;
                    }
                }
            }

            // derivatives with respect to qdot, columns NDOF+1 to 2*NDOF of x
            for (i=0; i<NDOF; i++) {
                df_dx_jcs[NDOF+i] = Ndf_dx;		// store element number where this column starts

                // derivatives of (qdot-dq/dt) with respect to qdot are diagonal, in rows 1 to NDOF
                df_dx_irs[Ndf_dx] = i;		// store row number of this matrix element
                df_dx[Ndf_dx] = 1.0;		// store the value of this matrix element
                Ndf_dx++;

                // derivatives of ZERO with respect to qdot are in rows NDOF+1 to 2*NDOF
                for (j=0; j<NDOF; j++) {
                    df_dx_irs[Ndf_dx] = NDOF+j;			// store row number of this matrix element
                    df_dx[Ndf_dx] = dz_dqd[j][i];		// store the value of this matrix element
                    // add the contributions dz/dmom * dmom/dqdot (remember dmom_dangvel is diagonal)
                    df_dx[Ndf_dx] = df_dx[Ndf_dx] + dz_dmom[j][i]*dmom_dangvel[i];
                    Ndf_dx++;
                }

            }

            // derivatives with respect to Lce, columns 2*NDOF+1 to 2*NDOF+NMUS
            for (i=0; i<NMUS; i++) {			
                df_dx_jcs[2*NDOF+i] = Ndf_dx;		// store element number where this column starts

                // derivatives of ZERO with respect to Lce are in rows NDOF+1 to 2*NDOF
                for (j=0; j<NDOF; j++) {
                    df_dx_irs[Ndf_dx] = NDOF+j;		// store row number of this matrix element
                    df_dx[Ndf_dx] = 0.0;			// store the value of this matrix element
                    for (k=0; k<NDOF; k++) {
                        df_dx[Ndf_dx] = df_dx[Ndf_dx] + dz_dmom[j][k]*dmom_dLce[k][i];
                    }
                    Ndf_dx++;
                }

                // derivatives of muscle force balance with respect to Lce are diagonal, rows 2*NDOF+1 to 2*NDOF+NMUS
                df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
                df_dx[Ndf_dx] = dg_dLce[i];			// store the value of this matrix element
                Ndf_dx++;				
            }

            // derivatives with respect to Act, columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
            for (i=0; i<NMUS; i++) {			
                df_dx_jcs[2*NDOF+NMUS+i] = Ndf_dx;		// store element number where this column starts

                // derivatives of muscle force balance with respect to Act are diagonal, rows 2*NDOF+1 to 2*NDOF+NMUS
                df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
                df_dx[Ndf_dx] = dg_da[i];		// store the value of this matrix element
                Ndf_dx++;				

                // derivatives of activation dynamics with respect to Act are diagonal, rows 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
                df_dx_irs[Ndf_dx] = 2*NDOF+NMUS+i;		// store row number of this matrix element
                df_dx[Ndf_dx] = dh_da[i];		// store the value of this matrix element
                Ndf_dx++;				
            }
            df_dx_jcs[NSTATES] = Ndf_dx;		// store final element number

            // --------Third output of MEX function: Jacobian df/dxdot
            plhs[2] = mxCreateSparse(NSTATES, NSTATES, NDOF*(NDOF+1) + 2*NMUS, 0);
            df_dxdot = mxGetPr(plhs[2]);
            df_dxdot_irs = mxGetIr(plhs[2]);
            df_dxdot_jcs = mxGetJc(plhs[2]);
            Ndf_dxdot = 0;

            // derivatives with respect to dq/dt, columns 1..NDOF of xdot
            for (i=0; i<NDOF; i++) {			
                df_dxdot_jcs[i] = Ndf_dxdot;			// store element number where this column starts

                // derivatives of (qdot-dq/dt) with respect to dq/dt are diagonal, in rows 1 to NDOF
                df_dxdot_irs[Ndf_dxdot] = i;		// store row number of this matrix element
                df_dxdot[Ndf_dxdot] = -1.0;			// store the value of this matrix element
                Ndf_dxdot++;
            }

            // derivatives with respect to dqdot/dt, columns NDOF+1 to 2*NDOF of xdot
            for (i=0; i<NDOF; i++) {
                df_dxdot_jcs[NDOF+i] = Ndf_dxdot;		// store element number where this column starts

                // derivatives of ZERO with respect to qdd are in rows NDOF to 2*NDOF
                for (j=0; j<NDOF; j++) {
                    df_dxdot_irs[Ndf_dxdot] = NDOF+j;		// store row number of this matrix element
                    df_dxdot[Ndf_dxdot] = dz_dqdd[j][i];	// store the value of this matrix element
                    Ndf_dxdot++;
                }
            }

            // derivatives with respect to Lcedot, columns 2*NDOF+1 to 2*NDOF+NMUS
            for (i=0; i<NMUS; i++) {			
                df_dxdot_jcs[2*NDOF+i] = Ndf_dxdot;		// store element number where this column starts

                // derivatives of muscle force balance with respect to Lcedot are diagonal, rows 2*NDOF to 2*NDOF+NMUS
                df_dxdot_irs[Ndf_dxdot] = 2*NDOF+i;		// store row number of this matrix element
                df_dxdot[Ndf_dxdot] = dg_dLcedot[i];	// store the value of this matrix element
                Ndf_dxdot++;
            }

            // derivatives with respect to Actdot, columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
            for (i=0; i<NMUS; i++) {			
                df_dxdot_jcs[2*NDOF+NMUS+i] = Ndf_dxdot;		// store element number where this column starts

                // derivatives of activation dynamics with respect to Actdot are diagonal, rows 2*NDOF+NMUS to 2*NDOF+2*NMUS
                df_dxdot_irs[Ndf_dxdot] = 2*NDOF+NMUS+i;		// store row number of this matrix element
                df_dxdot[Ndf_dxdot] = dh_dadot[i];				// store the value of this matrix element
                Ndf_dxdot++;				
            }
            df_dxdot_jcs[NSTATES] = Ndf_dxdot;		// store final element number

            // --------Fourth output of MEX function: Jacobian df/du
            plhs[3] = mxCreateSparse(NSTATES, NMUS, NMUS, 0);
            df_du = mxGetPr(plhs[3]);
            df_du_irs = mxGetIr(plhs[3]);
            df_du_jcs = mxGetJc(plhs[3]);
            Ndf_du = 0;

            // derivatives with respect to u of each muscle
            for (i=0; i<NMUS; i++) {			
                df_du_jcs[i] = Ndf_du;		// store element number where this column starts

                // derivatives of activation dynamics with respect to u are diagonal, rows 2*NDOF+NMUS to 2*NDOF+2*NMUS
                df_du_irs[Ndf_du] = 2*NDOF+NMUS+i;		// store row number of this matrix element
                df_du[Ndf_du] = dh_du[i];				// store the value of this matrix element
                Ndf_du++;				
            }
            df_du_jcs[NMUS] = Ndf_du;		// store final element number
        #else
            // this section of the code generates non-sparse Jacobians, outputs 2,3,4 of the MEX function
            plhs[1] = mxCreateDoubleMatrix(NSTATES, NSTATES, mxREAL);   // df/dx
            plhs[2] = mxCreateDoubleMatrix(NSTATES, NSTATES, mxREAL);	// df/dxdot
            plhs[3] = mxCreateDoubleMatrix(NSTATES, NMUS,    mxREAL);   // df/du
			df_dx    = mxGetPr(plhs[1]);
            df_dxdot = mxGetPr(plhs[2]);	
            df_du    = mxGetPr(plhs[3]);

            // the first NDOF rows of the implicit differential equation (IDE) are: qdot-dq/dt = 0
            for (i=0; i<NDOF; i++) {
                // f[i] = x[NDOF+i] - xdot[i];   // this is the f[] output for these rows
                df_dx[i + nrows*(NDOF+i)] = 1.0;  // df/dx row i, column NDOF+1
                df_dxdot[i + nrows*i]   = -1.0; // df/dxdot row i, column i
            }
			
			// the next NDOF rows of the IDE are the equations of motion from Autolev, the "zero" expressions
            for (i=0; i<NDOF; i++) {
                // f[NDOF+i] = zero[i];         // this is the f[] output for these rows
                
                // put dzero/dq in df/dx columns 1..NDOF
                for (j=0; j<NDOF; j++) {
                    mwIndex ii = NDOF+i + nrows*j;       // index of row NDOF+i, column j
                    df_dx[ii] = dz_dq[i][j];                // df/dx row NDOF+i, column j
                    // add the contributions dz/dmom * dmom/dq to this element of df/dx
                    for (k=0; k<NDOF; k++) {
                        df_dx[ii] = df_dx[ii] + dz_dmom[i][k]*dmom_dang[k][j];
                    }
                }
                
                // put dzero/dqot in df/dx columns NDOF+1 to 2*NDOF
                for (j=0; j<NDOF; j++) {
                    mwIndex ii = NDOF+i + nrows*(NDOF+j);
                    df_dx[ii] = dz_dqd[i][j];  // df/dxdot, row NDOF+i, column NDOF+j
                    // add the contributions dz/dmom * dmom/dqdot (knowing that dmom_dangvel is diagonal)
                    df_dx[ii] = df_dx[ii] + dz_dmom[i][j]*dmom_dangvel[j];
                }
                
                // put dzero/dLce in df/dx columns 2*NDOF+1 to 2*NDOF+NMUS
                for (j=0; j<NMUS; j++) {
                    mwIndex ii = NDOF+i + nrows*(2*NDOF+j);   // index of row NDOF+1, column 2*NDOF+j
                    df_dx[ii] = 0.0;                            // initialize this matrix element
                    for (k=0; k<NDOF; k++) {                    // contributions from joint moments
                        df_dx[ii] = df_dx[ii] + dz_dmom[i][k]*dmom_dLce[k][j];
                    }
                }
                
                // put dzero/dqdd in df/dxdot columns NDOF+1 to 2*NDOF
                for (j=0; j<NDOF; j++) {
                    mwIndex ii = NDOF+i + nrows*(NDOF+j);   // index of row NDOF+i, column NDOF+j
                    df_dxdot[ii] = dz_dqdd[i][j];	
                }
            }

            // the next NMUS rows of the IDE are the muscle contraction dynamics
            for (i=0; i<NMUS; i++) {
                // f[i] = g[i]; 
                
                // put dg/dq in df/dx columns 1 to NDOF
                for (j=0; j<NDOF; j++) {
                    mwIndex ii = 2*NDOF+i + nrows*j;     // index of row 2*NDOF+i, column j
                    df_dx[ii] = dg_dLm[i]*dLm_dang[i][j];  // store the value of this matrix element
                }
                
                // put dg/dLce in df/dx columns 2*NDOF+1 to 2*NDOF+NMUS
                // there is only one element on each row, from Lce of the same muscle
                df_dx[2*NDOF+i + nrows*(2*NDOF+i)] = dg_dLce[i];

                // put dg/dAct in df/dx columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
                // there is only one element on each row, from active state of the same muscle
                df_dx[2*NDOF+i + nrows*(2*NDOF+NMUS+i)] = dg_da[i];		
       
                // put dg/dLcedot in df/dxdot columns 2*NDOF+1 to 2*NDOF+NMUS
                // there is only one element on each row, from Lcedot of the same muscle
                df_dxdot[2*NDOF+i + nrows*(2*NDOF+i)] = dg_dLcedot[i];
            }
                
            // the final NMUS rows of the IDE are the muscle activation dynamics: h = da/dt - (u-a)(c1*u + c2) = 0
            for (i=0; i<NMUS; i++) {
                // f[2*NDOF+NMUS+i] =  h[i];
                
                // put dh/da in df/dx columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
                // only one element on each row (dh/da from same muscle)
                df_dx[2*NDOF+NMUS+i + nrows*(2*NDOF+NMUS+i)] = dh_da[i];		
                
                // put dh/dadot in df/dxdot columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
                // only one element on each row (dh/dadot from same muscle)   
                df_dxdot[2*NDOF+NMUS+i + nrows*(2*NDOF+NMUS+i)] = dh_dadot[i];	
                
                // put dh/du in df/du columns 1 to NMUS
                df_du[2*NDOF+NMUS+i + nrows*i] = dh_du[i];	

			}
		#endif
	}
}
