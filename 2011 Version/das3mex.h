// das3mex.h

// This file must be included in both das3mex.c and das3_al.c to ensure that declarations
// and constants remain consistent between them.

#define NDOF 11			// number of kinematic degrees of freedom
#define NSTICK 9		// number of stick figure points
#define NMUS 138		// number of muscles
#define NSTATES 298		// number of states (2*NDOF + 2*NMUS)

// struct with model parameters that must be shared with Autolev code
typedef struct {

	// parameters for scapula-thorax contact model
	double epscontact;	// epsilon for scapula-thorax contact model
	double kcontact;	// stiffness (N/m)

	// parameters for conoid model
	double epsconoid;			// epsilon for conoid force-length model, 2 mm is recommended
	double conK, conL;			// Stiffness (N/m) and length (cm) of conoid ligament
	double conOx,conOy,conOz;   // coordinates of conoid origin in clavicle frame
	double conIx,conIy,conIz;   // coordinates of conoid insertion in scapula frame

	// morphological parameters
	double SCx, SCy, SCz;	// coordinates of SC in thorax frame
	double ACx, ACy, ACz;	// coordinates of AC in clavicle frame
	double GHx, GHy, GHz;	// coordinates of GH in scapula frame
	double ULx, ULy, ULz;	// coordinates of UL in humerus frame
	double RDx, RDy, RDz;	// coordinates of RD in ulna frame
	double HDx, HDy, HDz;	// coordinates of HD in radius frame
	double Ax, Ay, Az;		// radii of the thorax ellipsoid
	double TSx, TSy, TSz;	// coordinates of TS in scapula frame
	double AIx, AIy, AIz;	// coordinates of AI in scapula frame
	double FlexAxisX, FlexAxisY, FlexAxisZ;		// Elbow flexion axis orientation, in humerus frame
	double PronAxisX, PronAxisY, PronAxisZ;		// Forearm pronation axis orientation, in ulna frame
	
	// mass of the six bodies
	double Mass_Clavicle;	
	double Mass_Scapula;	
	double Mass_Humerus;
	double Mass_Ulna;	
	double Mass_Radius;	
	double Mass_Hand;	

	// center of mass of the six bodies
	double ClavCMx, ClavCMy, ClavCMz;		
	double ScapCMx, ScapCMy, ScapCMz;
	double HumCMx,  HumCMy,  HumCMz;
	double UlnaCMx, UlnaCMy, UlnaCMz;
	double RadCMx,  RadCMy,  RadCMz;
	double HandCMx, HandCMy, HandCMz;

	// moments of inertia of the six bodies
	double ClavIxx, ClavIyy, ClavIzz, ClavIxy, ClavIyz, ClavIxz;		
	double ScapIxx, ScapIyy, ScapIzz, ScapIxy, ScapIyz, ScapIxz;
	double HumIxx,  HumIyy,  HumIzz , HumIxy , HumIyz , HumIxz;
	double UlnaIxx, UlnaIyy, UlnaIzz, UlnaIxy, UlnaIyz, UlnaIxz;
	double RadIxx,  RadIyy,  RadIzz , RadIxy , RadIyz , RadIxz;
	double HandIxx, HandIyy, HandIzz, HandIxy, HandIyz, HandIxz;
	
} param_struct;


// prototype for the Autolev C function
void das3_al(
	param_struct* param,				// input: pointer to struct containing parameter
	double q[NDOF], 					// input: joint angles
	double qd[NDOF], 					// input: angular velocities
	double qdd[NDOF],					// input: angular accelerations
	double mTH[3],						// input: joint moments in "thorax-humerus" joint
	double exF[2],						// input: simulated arm support: [distance from elbow, force amplitude]
	double Zero[NDOF],					// output: dynamics imbalance
	double dz_dq[NDOF][NDOF], 			// output: Jacobian of Zero with respect to q
	double dz_dqd[NDOF][NDOF], 			// output: Jacobian of Zero with respect to qdot
	double dz_dqdd[NDOF][NDOF], 		// output: Jacobian of Zero with respect to qdotdot
	double F_GH[3],				    	// output: net reaction force in GH joint
	double F_SCAP[2][3],				// output: scapula contact force at TS and AI
	double Stick[NSTICK][3],			// output: 9-point stick figure data
	double qTH[3]);						// output: thorax humerus angles using YZY sequence
	
// prototype for the das3step function in das3mex.c
int das3step(
	double x[NSTATES],		// input/output: state before and after the time step
	double uin[NMUS], 		// input: muscle controls
	double hh,				// input: time step size (must be small enough to ensure stability and accuracy
	double exF[2],			// input: vertical force of amplitude exF(2) applied to the ulna at a distance of exF(1) from the elbow (to simulate a mobile arm support)
	double F_GH[3],			// output: reaction force in the GH joint (for detecting dislocation)
	double F_SCAP[2][3],	// output: forces at the two scapula-thorax contact points
	double f[NSTATES]);		// output: dynamics residuals (when they get large, the solver is not working well)





