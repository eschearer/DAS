/*
	autolevclean.c

	Program to extract clean C source code from C file generated by Autolev program das3.al
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "das3mex.h"

#define STRLEN 100 				// source lines are never longer than 100 characters

int main() {

	FILE *fid1, *fid2;
	char line[STRLEN];
	int i, copying;
	char *ptr;

// open the raw source file that came from Autolev
	if ((fid1 = fopen("das3_al_raw.c","r")) == NULL) {
		printf("Could not open das3_al_raw.c\n");
		exit(1);
	}
	
// write the clean C file
	if ((fid2 = fopen("das3_al.c","w")) == NULL) {
		printf("Could not write das3_al.c\n");
		exit(1);
	}
	
	// write the function header
	fprintf(fid2,"// This file was generated by autolevclean.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include \"das3mex.h\"\n");
	fprintf(fid2,"#include <math.h>\n\n");
	fprintf(fid2,"void das3_al(param_struct* par, double q[NDOF], double qd[NDOF], double qdd[NDOF],\n");
	fprintf(fid2,"   double mTH[3], double exF[2], double Zero[NDOF], double dz_dq[NDOF][NDOF],\n");
	fprintf(fid2,"   double dz_dqd[NDOF][NDOF], double dz_dqdd[NDOF][NDOF],\n");
	fprintf(fid2,"   double F_GH[3], double F_SCAP[2][3],\n");
	fprintf(fid2,"   double Stick[NSTICK][3], double qTH[3]) {\n");
	
	// make a macro for the "sign" function
	fprintf(fid2,"#define sign(x) ((x > 0) - (x < 0))\n");
	
	// generate C code to copy q, qd, qdd into scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);		
		fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);	
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);	
	}
	
	// declare some internal variables
	fprintf(fid2,"\tdouble MHx, MHy, MHz;\n");
	
	// generate C code to copy the 3 thorax-humerus into scalar variables
	fprintf(fid2,"\tdouble MTHy = mTH[0];\n");
	fprintf(fid2,"\tdouble MTHz = mTH[1];\n");
	fprintf(fid2,"\tdouble MTHyy = mTH[2];\n");

	// generate C code to copy forearm force description (arm support) into scalar variables
	fprintf(fid2,"\tdouble distF = exF[0];\n");
	fprintf(fid2,"\tdouble ampF = exF[1];\n");

	// generate C code to declare the temporary variables used in contact force model and conoid force model
	fprintf(fid2,"\tdouble PxTS,PyTS,PzTS,FxTS,FyTS,FzTS,FTS,FminusTS;\n");				
	fprintf(fid2,"\tdouble PxAI,PyAI,PzAI,FxAI,FyAI,FzAI,FAI,FminusAI;\n");				
	fprintf(fid2,"\tdouble LOI, Stretch, StretchPositive;\n");			

	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// make sure there is enough room for all Zs
		}
		else if (strcmp(line, "/* Evaluate constants */\n") == 0) 				// Z[] code starts here
			copying = 1;
		else if (strcmp(line, "/* Evaluate output quantities */\n") == 0)	// Z[] code ends here
			copying = 0;
		else if (strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0)		// kinematics code starts here
			copying = 1;
		else if (strcmp(line, "  Encode[0] = 0.0;\n") == 0)	// and stop here
			copying = 0;
		else if (copying) {
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
	}
	
// close the input file
	fclose(fid1);
	
// close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);
		
}			// END OF MAIN PROGRAM
